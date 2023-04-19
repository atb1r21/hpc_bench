!=============================================================================!
!                            S O L V A T I O N                                !
!=============================================================================!
! This module defines the implicit solvation model outlined in the following  !
! papers and references therein:                                              !
!                                                                             !
! [1] DA Scherlis, J-L Fattebert, F Gygi, M Cococcioni, and N Marzari         !
!     Journal of Chemical Physics 124, 074103 (2006)                          !
! [2] O Andreussi, I Dabo, N Marzari,                                         !
!     Journal of Chemical Physics 136, 064102 (2012).                         !
!                                                                             !
!-----------------------------------------------------------------------------!
! Code for CASTEP written by Hatem H Helal, starting in 8/2007.               !
!                                                                             !
! ONETEP implementation by Jacek Dziedzic, 04-05/2010.                        !
! Made parallel-ready by Jacek Dziedzic, 04-05/2010.                          !
!                                                                             !
! Please report bugs to jd12g09@soton.ac.uk                                   !
!-----------------------------------------------------------------------------!

module is_solvation

  use constants, only: DP
  use is_poisson, only: POISSON_PROBLEM, SOLVATION_ENERGY_TERMS
  use utils, only: utils_trace_in, utils_trace_out
#ifdef HAVE_DL_MG
  use dl_mg  !! External dependency
#endif

  implicit none

  ! Everything is private ...
  private

  !... unless exposed here.
  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !                              a n d   t y p e s                            !
  !---------------------------------------------------------------------------!

  public :: implicit_solvent_exit
  public :: implicit_solvent_hartree
  public :: initialize_solvation_problem

  public :: have_rho_ion     ! ndmh: to allow resetting of rho_ion to zero
  public :: have_initial_eps ! jd: to allow resetting from energy_and_force
  public :: implicit_solvent_energy_terms   ! jd: Ugly state needed to print out
                                            ! final decomposition of solvation
                                            ! energies in energy_and_force

  logical, save       :: have_initial_eps = .false.
  logical, save       :: have_rho_ion     = .false.

  type(SOLVATION_ENERGY_TERMS), save :: implicit_solvent_energy_terms

  !---------------------------------------------------------------------------!
  !                     P r i v a t e   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!

  type(POISSON_PROBLEM),  save               :: initial_problem

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_hartree(phi_tot, rho_elec, grid, cell, &
       E_Hartree, solvation_terms, elements, no_dump)

    !=========================================================================!
    ! This is the main driver for the implicit solvation method which will    !
    ! use the parameters 'solvation_method' and 'dielectric_model' to setup   !
    ! the solvation problem with the correct dielectric cavity and select     !
    ! either the direct or corrective approaches to solving the Poisson eq in !
    ! dielectric.  The dielectric is either self-consistently updated along   !
    ! along with the potential or else it is held fixed.  For the latter it   !
    ! is essential that the starting charge density is something sensible     !
    ! since this is used to initialize the dielectric.                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi_tot,          intent=out, total potential: hartree + dielectric   !
    !   rho_elec,         intent=in,  the input electronic density            !
    !   E_Hartree, (opt)  intent=out, calculated Hartree energy               !
    !   solvation_terms, (opt, inout), if passed, its components will be      !
    !                                  populated with a breakdown of energy   !
    !                                  terms relevant in solvation.           !
    !   no_dump, (opt)     intent=in,  logical parameter to avoid 3D dumps    !
    ! Returns:                                                                !
    !   The (molecular) Hartree energy in the presence of dielectric.         !
    !   Note that in the case of self-consistently changing dielectric, this  !
    !   energy is different than the one obtained by integrating the returned !
    !   phi_tot with the molecular density. This is because the returned      !
    !   phi_tot contains additional contributions (V_eps, V_apolar) from the  !
    !   derivative of the permittivity wrt changing density.                  !
    ! NB: This subroutine is spin-agnostic.                                   !
    ! Q. Why are E_apolar_cavitation and E_apolar_disrep (parts of the        !
    !    optional solvation_terms) calculated here?                           !
    ! A. Because current_problem lives only for the duration of this routine. !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009.                              !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010.                       !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010.                      !
    ! Updated to separate the apolar term by Jacek Dziedzic, 08/2016.         !
    ! Adjusted for Boltzmann solvation by Jacek Dziedzic, 11/2019.            !
    !=========================================================================!

    use bibliography, only: bibliography_cite
    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use is_poisson, only: deallocate_poisson_problem
    use rundat, only: pub_is_check_solv_energy_grad, pub_is_dielectric_model, &
         pub_is_solvation_method, pub_is_include_apolar, &
         pub_is_dielectric_function
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort
    use ion, only: ELEMENT

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(ELEMENT), intent(in), optional  :: elements(:)
    real(kind=DP), intent(out)  :: phi_tot(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    real(kind=DP), intent(in)   :: rho_elec(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    real(kind=DP), intent(out), optional :: E_Hartree
    type(SOLVATION_ENERGY_TERMS), intent(inout), optional :: solvation_terms
    logical, intent(in), optional :: no_dump

    ! Internal variables
    type(POISSON_PROBLEM) :: current_problem
    real(kind=DP) :: E_elec_fixed_vac_bak

    !------------------------------------------------------------------------

    call utils_trace_in('implicit_solvent_hartree')
    call bibliography_cite('IS')

    if (pub_is_dielectric_function=='SOFT_SPHERE') then
       call bibliography_cite('Fisicaro_Soft_Sphere')
       call bibliography_cite('Alvarez_vdW_Radii')
    end if

    ! jd: Allocate arrays and evaluate dielectric.
    call initialize_solvation_problem(rho_elec, current_problem, grid, cell, &
         elements)

    select case(pub_is_solvation_method)
    case ('DIRECT')
       ! Directly solve the Poisson equation in dielectric to find the
       ! solvation potential
       call direct_solvation_potential(current_problem,grid,cell,elements,no_dump)

    case ('CORRECTIVE')
       ! Use the corrective potential method to get the solvation potential
       call utils_abort('Corrective approach not (yet) implemented in ONETEP.')

    case default
       call utils_abort('Unrecognized is_solvation_method in &
            &implicit_solvent_hartree.')
    end select

    if(present(solvation_terms)) then
       solvation_terms = current_problem%solvation_terms
    end if

    ! Now we have the solvation potential we need to calculate the
    ! dielectric potential and energy correction
    select case(pub_is_dielectric_model)
    case ('FIX_INITIAL') ! The simple case...
       if(present(solvation_terms)) then
          solvation_terms%E_apolar_cavitation = &
               initial_problem%solvation_terms%E_apolar_cavitation
          solvation_terms%E_apolar_disrep = &
               initial_problem%solvation_terms%E_apolar_disrep
          solvation_terms%E_apolar_sasa = &
               initial_problem%solvation_terms%E_apolar_sasa
          solvation_terms%E_apolar_sav = &
               initial_problem%solvation_terms%E_apolar_sav
       end if

    case ('SELF_CONSISTENT') ! The more involved case...
       ! jd: Include V_eps correction (eq. (6) in [1]) in phi_mol.
       !     Energy uses uncorrected phi and has been calculated earlier, so OK.
       call calc_diel_correction(current_problem,grid)

       current_problem%phi_mol = current_problem%phi_mol + current_problem%V_eps

       ! jd: Include V_apolar correction (eq. (14) in [1]) in phi_mol.
       !     Energy uses uncorrected phi and has been calculated earlier, so OK.
       if(pub_is_include_apolar) then
          current_problem%phi_mol = &
               current_problem%phi_mol + current_problem%V_apolar
       end if

    case default
       call utils_abort('Unrecognized is_dielectric_model in &
            &implicit_solvent_hartree.')
    end select

    ! Finally we set up pot variable for returning
    phi_tot = current_problem%phi_mol

    ! Check electrostatic energy gradient if necessary
    if(pub_is_check_solv_energy_grad) then
       call check_elec_energy_gradient(current_problem,grid,cell,elements)
    end if

    ! Write out data files for later analysis
    call write_solvation_data(current_problem,grid,cell)

    ! All done so deallocate solvation problem that is no longer needed
    call deallocate_poisson_problem(current_problem,pub_is_solvation_method)

    if(present(E_Hartree)) E_Hartree = current_problem%E_Hartree

    ! jd: Store the most recent cavitation, dispersion-repulsion and Boltzmann
    !     energies in ugly module-wide state, so that they are accessible at
    !     the end of energy_and_force_calculate() for pretty-printing solvation
    !     energy components in auto solvation. Do not overwrite the in-vacuum
    !     term.
    if(present(solvation_terms)) then
       E_elec_fixed_vac_bak = implicit_solvent_energy_terms%E_elec_fixed_vac
       implicit_solvent_energy_terms = solvation_terms
       implicit_solvent_energy_terms%E_elec_fixed_vac = E_elec_fixed_vac_bak
       implicit_solvent_energy_terms%cavity_surface_area=&
            current_problem%cavity_surface_area
       implicit_solvent_energy_terms%cavity_volume=&
            current_problem%cavity_volume
    end if

    call utils_trace_out('implicit_solvent_hartree')

  end subroutine implicit_solvent_hartree
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_exit
    !=========================================================================!
    ! Cleans up after implicit solvent.                                       !
    ! Currently a no-op. Boltzmann solvation is cleaned up separately, and    !
    ! after every force evaluation (as steric pot needs to be recalculated).  !
    !=========================================================================!

    implicit none

    !------------------------------------------------------------------------

  end subroutine implicit_solvent_exit
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine initialize_solvation_problem(den,problem,grid,cell,elements)
    !=========================================================================!
    ! This subroutine will allocate the necessary globally available arrays   !
    ! and initialize the dielectric depending on the chosen method.           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   den,              intent(in),    the input electron density.          !
    !   problem,          intent(inout), the POISSON_PROBLEM to initialize.   !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    ! Updated to support periodic BCs by James C. Womack, 05/2017             !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE
    use is_poisson, only: POISSON_PROBLEM, allocate_poisson_problem
    use is_smeared_ions, only: rho_ion
    use rundat, only: pub_is_dielectric_model, pub_is_solvation_method, &
         pub_output_detail, pub_multigrid_in_use, pub_multigrid_bc_is_periodic
    use utils, only: utils_abort, utils_assert
    use simulation_cell, only: CELL_INFO
    use ion, only: ELEMENT

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(ELEMENT), intent(in)   :: elements(:)
    real(kind=DP), intent(in)   :: den(grid%ld1, grid%ld2, grid%max_slabs12)
    type(POISSON_PROBLEM), intent(inout) :: problem

    !------------------------------------------------------------------------

    call utils_trace_in('initialize_solvation_problem')

    ! JCW: Check that multigrid solver is in use
    ! JCW: Important, since pub_multigrid_bc_is_periodic is undefined when
    ! JCW: pub_multigird_in_use is .false.
    call utils_assert(pub_multigrid_in_use,"Error in initialize_solvation_problem: &
         &Multigrid solver must be active for solvation calculation to proceed &
         &(pub_multigrid_in_use must be .true.")

    ! Start by allocating solvation problem arrays
    call allocate_poisson_problem(problem,pub_is_solvation_method,grid)

    ! Check if initial problem has been allocated and do so if needed
    if(.not. initial_problem%is_allocated) then
       call allocate_poisson_problem(initial_problem,'INITIAL',grid)
    end if

    ! jd: Set up the rho_elec array
    problem%rho_elec = den

    ! Retrieve ion density if its already been stored in initial data
    if (have_rho_ion) then
       problem%rho_ion = initial_problem%rho_ion
    ! Otherwise this is our first time through so calculate and store it
    else
       problem%rho_ion = rho_ion
       initial_problem%rho_ion = problem%rho_ion
       have_rho_ion = .true.
    end if

    ! Find out if we are using a previously stored dielectric
    select case(pub_is_dielectric_model)
    case('FIX_INITIAL')

      if (have_initial_eps) then
        ! The dielectric is fixed and we already have the initial data
        ! so we simply retrieve the saved dielectric and cavitation data
        problem%eps_full = initial_problem%eps_full
        problem%eps_half = initial_problem%eps_half

        problem%cavity_surface_area = initial_problem%cavity_surface_area
        problem%cavity_volume = initial_problem%cavity_volume
        problem%solvation_terms = initial_problem%solvation_terms

        if(pub_on_root .and. pub_output_detail >= VERBOSE) then
           write(stdout,'(/a)') 'Reusing previously generated solute cavity.'
        end if

      else

        ! We save initial electronic charge density
        initial_problem%rho_elec = problem%rho_elec

        ! Calculate dielectric on full grid
        call calc_dielectric_medium(initial_problem%rho_elec, &
             problem%eps_full, grid, cell, elements)
        initial_problem%eps_full = problem%eps_full

        ! Calculate dielectric on half grids
        call interpolate_dielectric_medium(initial_problem%rho_elec, &
             problem%eps_half, grid, cell, elements)
        initial_problem%eps_half = problem%eps_half

        ! ...and also calculate the surface area of the current cavity and
        !    the corresponding apolar terms.
        call calc_apolar_term(initial_problem, grid, pub_multigrid_bc_is_periodic)
        ! JCW: (use the periodicity for the multigrid solver to determine how the
        ! JCW: density gradients required in the apolar term will be computed)

        problem%cavity_surface_area = initial_problem%cavity_surface_area
        problem%cavity_volume = initial_problem%cavity_volume
        problem%solvation_terms = initial_problem%solvation_terms

        ! Remember we've done that
        have_initial_eps = .true.

      end if

    case('SELF_CONSISTENT')

      ! We can go ahead and use current charge density to update dielectric
      call calc_dielectric_medium(problem%rho_elec,problem%eps_full,grid, &
           cell,elements)
      call interpolate_dielectric_medium(problem%rho_elec,problem%eps_half,grid, &
           cell,elements)

      ! Calculate the surface area of the current cavity and the
      ! corresponding apolar terms.
      ! jd: In the SCF case, also calculate V_apolar (eq. (14) in [1]).
      call calc_apolar_term(problem, grid, pub_multigrid_bc_is_periodic)
      ! JCW: (use the periodicity for the multigrid solver to determine how the
      ! JCW: density gradients required in the apolar term will be computed)

    case default
      call utils_abort('Unrecognized pub_is_dielectric_model.')
    end select

    call utils_trace_out('initialize_solvation_problem')

  end subroutine initialize_solvation_problem
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine corrective_solvation_potential(solv_problem)
    !=========================================================================!
    ! A stub for future implementation of the corrective solvation potential. !
    !=========================================================================!

    use is_poisson, only: POISSON_PROBLEM
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(POISSON_PROBLEM), intent(inout) :: solv_problem

    ! -----------------------------------------------------------------------

    call utils_abort('Corrective solvation potential not implemented in onetep')
    solv_problem = solv_problem ! jd: Kill 'solv_problem usused' compiler warn.

  end subroutine corrective_solvation_potential
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine direct_solvation_potential(solv_problem,grid,cell,elements,no_dump)
    !=========================================================================!
    ! This subroutine uses the solvation type data to calculate the           !
    ! Hartree energy and potential in the presence of implicit solvent by     !
    ! directly solving the Poisson equation in dielectric.                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    !   no_dump, (opt)     intent=in,  logical parameter to avoid 3D dumps    !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use multigrid_methods, only: multigrid_calculate_hartree
    use rundat, only: pub_is_bulk_permittivity
    use simulation_cell, only: CELL_INFO
    use ion, only: ELEMENT

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(POISSON_PROBLEM), intent(inout) :: solv_problem
    type(ELEMENT), intent(in), optional :: elements(:)
    logical, intent(in), optional :: no_dump

    !------------------------------------------------------------------------

    call utils_trace_in('direct_solvation_potential')

    solv_problem%rho_elec = solv_problem%rho_elec + solv_problem%rho_ion
    solv_problem%E_Hartree = multigrid_calculate_hartree( &
         solv_problem%phi_mol, &                            ! output
         solv_problem%rho_elec, &  ! NB, rho_ion added earlier
         grid, cell, &
         uniform_eps = pub_is_bulk_permittivity, &          ! input
         eps_full = solv_problem%eps_full, &                ! input
         eps_half = solv_problem%eps_half, &                ! input
         solvation_terms = solv_problem%solvation_terms, &  ! output
         elements = elements, no_dump=no_dump)

    solv_problem%rho_elec = solv_problem%rho_elec - solv_problem%rho_ion

    call utils_trace_out('direct_solvation_potential')

  end subroutine direct_solvation_potential
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine calc_dielectric_medium(den,dielectric,grid,cell,elements)
    !=========================================================================!
    ! For a given electronic density (den) this subroutine will calculate the !
    ! dielectric functional (dielectric) for later use in solving the Poisson !
    ! equation in the presence of a dielectric.                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   den,              intent=in,  the input electron density              !
    !   dielectric        intent=out, the dielectric we seek                  !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use simulation_cell, only: CELL_INFO
    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE
    use rundat, only: pub_is_dielectric_function, pub_output_detail, &
         pub_multigrid_bc_is_periodic
    use utils, only: utils_abort
    use ion, only: ELEMENT

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(ELEMENT), intent(in)   :: elements(:)
    real(kind=DP), dimension(grid%ld1,grid%ld2, &
         grid%max_slabs12), intent(in)  :: den
    real(kind=DP), dimension(grid%ld1,grid%ld2, &
         grid%max_slabs12), intent(out) :: dielectric
    type(CELL_INFO), intent(in) :: cell

    !------------------------------------------------------------------------

    call utils_trace_in('calc_dielectric_medium')

    ! jd: Fill the dielectric array, using elemental functions
    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
       write(stdout,'(/a)',advance='no') &
            'Generating solute cavity for implicit solvation ...'
    end if

    select case(pub_is_dielectric_function)
    case('FGF')
       dielectric = eps_function_fattebert_gygi(den)
    case('GAUSSIAN')
       dielectric = eps_function_gaussian(den)
    case('ANDREUSSI')
       dielectric = eps_function_andreussi(den)
    case('SOFT_SPHERE')
       call calculate_eps_soft_sphere(dielectric, elements, &
            grid, cell, pub_multigrid_bc_is_periodic, 0)
    case default
       call utils_abort('Unrecognized dielectric fn in calc_dielectric_medium')
    end select

    ! jd: Set the dielectric to unity inside the cores, in case the electronic
    !     density drops too low because of the pseudos, which would create
    !     a spurious cavity.
    call exclude_cavity(dielectric,1.0_DP,grid)

    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
       write(stdout,'(a)') ' done'
    end if

    call utils_trace_out('calc_dielectric_medium')

  end subroutine calc_dielectric_medium
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine interpolate_dielectric_medium(den,eps_half,grid,cell,elements)
    !=========================================================================!
    ! For a given electronic density (den) this subroutine will use FFTs to   !
    ! interpolate the density onto 'half grids' and then calculates the       !
    ! dielectric functional (dielectric) for later use in solving the Poisson !
    ! equation in the presence of a dielectric.                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   den,              intent=in,  the input electron density              !
    !   eps_half          intent=out, the dielectric functional we seek       !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use simulation_cell, only: CELL_INFO
    use comms, only: pub_my_proc_id
    use constants, only: DP
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use geometry, only: magnitude
    use rundat, only: pub_is_dielectric_function, pub_is_bulk_permittivity, &
         pub_multigrid_bc_is_periodic
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    use ion, only: ELEMENT

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(ELEMENT), intent(in)   :: elements(:)
    real(kind=DP), dimension(grid%ld1, &
         grid%ld2, grid%max_slabs12),  intent(in)  :: den
    real(kind=DP), dimension(grid%ld1, &
         grid%ld2, grid%max_slabs12,3),intent(out) :: eps_half

    ! Internal variables
    integer                                           :: ig1, ig2, ig3, dir
    complex(kind=DP), allocatable, dimension(:,:,:)   :: rho_recip
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: phases
    real(kind=DP), allocatable, dimension(:,:,:,:)    :: rho_interp
    real(kind=DP), dimension(3)                       :: half_grids
    integer                                           :: ierr ! jd: Error flag
    real(kind=DP)                                     :: d1f, d2f, d3f
    real(kind=DP)                                     :: gvec(3)
    character(len=*), parameter :: myself = 'interpolate_dielectric_medium'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    ! Allocate memory
    allocate(rho_recip(grid%ld3, grid%ld2, &
         grid%max_slabs23),stat=ierr)
    call utils_alloc_check(myself,'rho_recip',ierr)
    allocate(phases(grid%ld3, grid%ld2, &
         grid%max_slabs23,3),stat=ierr)
    call utils_alloc_check(myself,'phases',ierr)
    allocate(rho_interp(grid%ld1, grid%ld2, &
         grid%max_slabs12,3),stat=ierr)
    call utils_alloc_check(myself,'rho_interp',ierr)

    ! Store half-grid distances
    d1f = magnitude(grid%da1)
    d2f = magnitude(grid%da2)
    d3f = magnitude(grid%da3)
    half_grids(1) = d1f * 0.5_DP
    half_grids(2) = d2f * 0.5_DP
    half_grids(3) = d3f * 0.5_DP

    ! Calculate the needed phase shifts to get the density on the half grids
    ! NB: As long as pub_recip_grid orders indices differently from eps_half,
    !     there is no way this can be done cache-efficiently.
    do ig1 = 1, grid%num_slabs23
       do ig2 = 1, grid%ld2
          do ig3 = 1, grid%ld3
             call cell_grid_recip_pt(gvec,ig1 + &
                  grid%first_slab23(pub_my_proc_id) - 1,ig2,ig3, &
                  grid)
             do dir = 1, 3
                phases(ig3,ig2,ig1,dir) = exp(cmplx(0.0_DP, &
                     gvec(dir) * half_grids(dir),kind=DP))
             end do
          end do
       end do
    end do

    ! Now transform density to reciprocal space
    call fourier_apply_cell_forward(den,rho_recip,grid)

    ! Compute product of rho(G)*phase shift, transform to real space to get
    ! density in rho_interp
    do dir=1, 3
       phases(:,:,:,dir) = phases(:,:,:,dir) * rho_recip(:,:,:)
       call fourier_apply_cell_backward(rho_interp(:,:,:,dir), &
            phases(:,:,:,dir),grid)
    end do

    ! Initialize dielectric to bulk @removeme
    eps_half(:,:,:,:) = pub_is_bulk_permittivity

    ! Fill in eps_half using the interpolated density for Fattebert-Gygi or
    ! the distance based dielectric function for soft sphere.
    select case(pub_is_dielectric_function)
    case('FGF')
       eps_half = eps_function_fattebert_gygi(rho_interp)
    case('GAUSSIAN')
       eps_half = eps_function_gaussian(rho_interp)
    case('ANDREUSSI')
       eps_half = eps_function_andreussi(rho_interp)
    case('SOFT_SPHERE')
       do dir=1,3
          call calculate_eps_soft_sphere(eps_half(:,:,:,dir), elements, &
               grid, cell, pub_multigrid_bc_is_periodic, dir)
       end do
    case default
       call utils_abort('Unrecognized dielectric fn in calc_dielectric_medium')
    end select

    ! jd: Set the dielectric to unity inside the cores, in case the electronic
    !     density drops too low because of the pseudos, which would create
    !     a spurious cavity.
    do dir=1, 3
       call exclude_cavity(eps_half(:,:,:,dir),1.0_DP,grid)
    end do

    ! jd: Clean up
    deallocate(rho_interp,stat=ierr)
    call utils_dealloc_check(myself,'rho_interp',ierr)
    deallocate(phases,stat=ierr)
    call utils_dealloc_check(myself,'phases',ierr)
    deallocate(rho_recip,stat=ierr)
    call utils_dealloc_check(myself,'rho_recip',ierr)

    call utils_trace_out(myself)

  end subroutine interpolate_dielectric_medium
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine calc_apolar_term(problem, grid, is_periodic)
    !=========================================================================!
    ! This subroutine will take a given electron density stored in 'den' to   !
    ! calculate the surface area of the associated dielectric cavity and the  !
    ! cavitation energy.  We base this calculation upon the relationship      !
    ! between the density and the dielectric functional.  The calculation is  !
    ! done by computing the difference in volume between two nearby density   !
    ! isovalues (giving the volume of a 'shell' of charge) and dividing       !
    ! by the local thickness of this charge film or shell - giving the surface!
    ! area of the associated dielectric cavity.                               !
    ! ----------------------------------------------------------------------- !
    ! Note on computing the gradient of the density (JCW, 19/07/17):          !
    ! The calculation of the apolar term requires the gradient of the         !
    ! density. How this is calculated is determined by the BCs (as determined !
    ! by the is_periodic array). In fully open BCs the gradient is calculated !
    ! using finite differences, while in fully periodic BCs the gradient is   !
    ! calculated in reciprocal space. This is necessary because ONETEP's      !
    ! finite difference routines do not support periodic BCs. The use of      !
    ! finite differences for open BCs ensures consistency with the earlier    !
    ! version of the solvation model where only open BCs were available.      !
    !                                                                         !
    ! In future, the method used to compute the gradient of the density may   !
    ! be made consistent across all BC types, though further testing is       !
    ! necessary to determine the impact of this (i.e. does using the gradient !
    ! computed in reciprocal space significantly affect the results?).        !
    !                                                                         !
    ! The devel_code MG:APOLAR_USE_RECIP_GRAD=T:MG forces the density         !
    ! gradient to be computed using the reciprocal space method, regardless   !
    ! of the BCs and also allows the second derivative of the density (for    !
    ! self-consistent cavity mode) to be computed using this method, too ---  !
    ! this is prevented by default.                                           !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    !   is_periodic (input): Logical array indicating whether the grid is     !
    !                        periodic (T) or not (F) along each lattice       !
    !                        vector. This is used to determine how the        !
    !                        gradient of the density is computed.             !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    ! Added a volume calculation, Jacek Dziedzic 06/2010                      !
    ! Separated cavitation from dispersion-repulstion, Jacek Dziedzic 08/2016.!
    ! Updated to support periodic BCs by James C. Womack, 05/2017             !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP
    use finite_differences, only: finite_difference_gradient
    use multigrid_methods, only: mg
    use rundat, only: pub_is_solvent_surf_tension, pub_finite_difference_order, &
         pub_is_dielectric_model, pub_is_include_apolar, &
         pub_is_apolar_scaling_factor, pub_is_dielectric_function, &
         pub_is_bulk_permittivity, pub_is_apolar_method, pub_is_solvent_pressure, &
         pub_is_apolar_sasa_definition
    use services, only: services_gradient_on_grid
    use utils, only: utils_assert, &
         utils_alloc_check, utils_dealloc_check, utils_sanity_check, utils_abort

    ! JCW: TODO Remove once computing density gradient in reciprocal space is
    ! JCW: well tested
    use rundat, only: pub_devel_code
    use utils, only: utils_devel_code

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(POISSON_PROBLEM), intent(inout) :: problem
    logical, intent(in) :: is_periodic(3)

    ! Internal variables
    real(kind=DP), dimension(:,:,:,:), allocatable :: grad_den
    real(kind=DP), dimension(:,:,:,:), allocatable :: grad_eps
    real(kind=DP), dimension(:,:,:,:,:), allocatable :: second_der
    real(kind=DP) :: lcl_mod_grad_den, lcl_mod_grad_eps
    real(kind=DP) :: lcl_diel
    real(kind=DP) :: bulk_eps
    real(kind=DP) :: smooth_a,smooth_b,smooth_0    !smoothed density isosurfaces
    integer       :: i1, i2, i3
    integer       :: ierr
    real(kind=DP) :: cavity_surface_area, cavity_volume
    logical       :: want_scf
    real(kind=DP) :: V_apolar_here
    integer       :: i, j ! jd: x, y or z
    character(len=*), parameter :: myself = 'calc_apolar_term'

    ! JCW:  TODO Remove once computing density gradient in reciprocal space is
    ! JCW:  well tested
    logical       :: use_recip_grad

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    want_scf = (pub_is_dielectric_model == 'SELF_CONSISTENT')

    ! JCW: TODO Remove once computing density gradient in reciprocal space is
    ! JCW: well tested
    ! JCW: devel_code MG:APOLAR_USE_RECIP_GRAD=T:MG activates the computation
    ! JCW: of the density gradient for use in calculating apolar energy
    ! JCW: contributions for all BCs -- if not present, finite differences are
    ! JCW: used for open BCs and the reciprocal space method for periodic/mixed
    ! JCW: BCs
    ! JCW: This also supresses the error that prevents calculation of second
    ! JCW: derivatives using the reciprocal space approach.
    use_recip_grad = utils_devel_code(.false.,"MG","APOLAR_USE_RECIP_GRAD",&
        pub_devel_code)

    ! jd: If the apolar term is not needed, don't bother with all this
    if(.not. pub_is_include_apolar) then
       problem%solvation_terms%E_apolar_cavitation = 0D0
       problem%solvation_terms%E_apolar_disrep = 0D0
       problem%solvation_terms%E_apolar_sasa = 0D0
       problem%solvation_terms%E_apolar_sav = 0D0
       problem%solvation_terms%cavity_surface_area = 0D0
       problem%solvation_terms%cavity_volume = 0D0
       problem%cavity_surface_area = 0D0
       problem%cavity_volume = 0D0
       if(want_scf) problem%V_apolar = 0D0
       call utils_trace_out(myself)
       return
    end if

    ! Allocate memory for density gradient calculation
    if (pub_is_dielectric_function=='SOFT_SPHERE') then
       allocate(grad_eps(grid%ld1, grid%ld2, grid%max_slabs12,3),stat=ierr)
       call utils_alloc_check(myself,'grad_eps',ierr)
    else
       allocate(grad_den(grid%ld1, grid%ld2, grid%max_slabs12,3),stat=ierr)
       call utils_alloc_check(myself,'grad_den',ierr)
    end if

    ! JCW: finite_difference gradient kept in open BCs for maximum backward
    ! JCW: compatibility
    ! JCW: TODO Consider completely replacing FD gradient with reciprocal space
    ! JCW: gradient in open BCs
    if ((.not.any(is_periodic)).and.(.not.use_recip_grad)) then
       ! JCW: Fully open BCs
       ! Use finite differences to calculate density gradient
       if (pub_is_dielectric_function=='SOFT_SPHERE') then
          ! gab: Soft sphere cavity apolar term is calculated through the
          ! gradient of the dielectric as opposed to density.
          call finite_difference_gradient(grad_eps, &                   ! out
               problem%eps_full, pub_finite_difference_order, mg, grid) ! in
       else
          call finite_difference_gradient(grad_den, &                   ! out
               problem%rho_elec, pub_finite_difference_order, mg, grid) ! in
       end if
    else
       ! JCW: Periodic BCs along at least one lattice vector
       ! JCW: (or explicitly requested by devel_code)
       ! JCW: ONETEP's finite_difference module does not support periodic BCs
       ! JCW: At the boundaries of the simulation cell, forward and backward
       ! JCW: differences are used for periodic BCs when central differences could
       ! JCW: be used.
       ! JCW: --> Obtain the gradient of the density in reciprocal space instead
       if (pub_is_dielectric_function=='SOFT_SPHERE') then
          ! gab: Soft sphere cavity apolar term is calculated through the
          ! gradient of the dielectric as opposed to density.
          call services_gradient_on_grid(grad_eps,problem%eps_full,grid)
       else
          call services_gradient_on_grid(grad_den,problem%rho_elec,grid)
       end if
    end if

    ! jd: Allocate memory for second derivatives, needed only in the
    !     self-consistent version
    if(want_scf) then
       allocate(second_der(grid%ld1, grid%ld2, grid%max_slabs12,3,3),stat=ierr)
       call utils_alloc_check(myself,'second_der',ierr)
       if ((.not.any(is_periodic)).and.(.not.use_recip_grad)) then
          ! JCW: Fully open BCs
          do i=1, 3
             call finite_difference_gradient(second_der(:,:,:,:,i), &   ! out
                  grad_den(:,:,:,i), pub_finite_difference_order, &     ! in
                  mg, grid)                                             ! in
          end do
       else
          ! JCW: Periodic BCs along at least one lattice vector.
          ! JCW: (or explicitly requested by devel_code)
          ! JCW: --> Obtain the second derivative of the density in reciprocal space
          call utils_assert(use_recip_grad,"Error in "//myself//": &
               &Attempting to calculate second derivatives of the density &
               &by applying derivative operator in reciprocal space, but this &
               &is not yet supported. To use experimental code for calculating &
               &the second derivative in this way, use devel_code &
               &'MG:APOLAR_USE_RECIP_GRAD=T:MG'.")
          ! TODO This needs careful testing to check whether taking the second
          !      derivative in this way is accurate
          ! TODO Implement second derivative so that the second derivative
          !      operator is applied in reciprocal space (only 1 forward / backward
          !      FFT, rather than two!).
          do i = 1, 3
             call services_gradient_on_grid(second_der(:,:,:,:,i),&
                  grad_den(:,:,:,i),grid)
          end do
       end if
    end if

    ! jd: Calculate the cavity surface area
    ! @improveme: not cache friendly wrt grad_den
    ! @improveme: not cache friendly wrt second_der
    ! @improveme: pdh knows of a more elegant way of calculating the surf area
    cavity_surface_area = 0.0_DP
    cavity_volume = 0.0_DP
    if(want_scf) problem%V_apolar = 0.0_DP ! jd: Take care of ld-pt padding

    select case(pub_is_dielectric_function)
    case('SOFT_SPHERE')
       ! @ improveme
       ! gab: Calculates the soft sphere apolar term, which is taken as the gradient
       ! of the dielectric term. Using conditional statement here avoids placing if
       ! within the working loop itself and uncessary computations, but repeats large
       ! elements of the code which looks rather untidy. May need more elegant solution.

       bulk_eps = pub_is_bulk_permittivity

       do i3=1, grid%num_my_slabs12
          do i2=1, grid%n2
             do i1=1, grid%n1

                lcl_diel = problem%eps_full(i1,i2,i3)

                smooth_0 = ( bulk_eps - lcl_diel ) / &
                     ( bulk_eps - 1.0_DP )

                ! Calculate local value of |grad den|
                lcl_mod_grad_eps = sqrt(grad_eps(i1,i2,i3,1)**2 + &
                     grad_eps(i1,i2,i3,2)**2 + grad_eps(i1,i2,i3,3)**2)

                cavity_surface_area = cavity_surface_area + lcl_mod_grad_eps

                cavity_volume = cavity_volume + smooth_0

             end do
          end do
       end do

       cavity_surface_area = cavity_surface_area / (bulk_eps - 1.0_DP)

    case('FGF')

       select case(pub_is_apolar_sasa_definition)
       case('DENSITY')
          call calc_sasa_sav_density(grid, problem, grad_den, second_der, &
               & cavity_surface_area, cavity_volume, want_scf)
       case('ISODENSITY')
          call calc_sasa_sav_isodensity(grid, problem, grad_den, second_der, &
               & cavity_surface_area, cavity_volume, want_scf)
       case default
          call utils_abort('Unrecognized SASA calculation method')
       end select

    case('ANDREUSSI')

       select case(pub_is_apolar_sasa_definition)
       case('DENSITY')
          call calc_sasa_sav_andreussi(grid, problem, grad_den, second_der, &
               & cavity_surface_area, cavity_volume, want_scf)
       case('ISODENSITY')
          call utils_abort("The SCCS model (is_dielectric_function ANDREUSSI) &
               &mandates using is_apolar_sasa_definition DENSITY.")
       case default
          call utils_abort('Unrecognized SASA calculation method')
       end select

    case('GAUSSIAN')
       call utils_abort('GAUSSIAN dielectric function not implemented in &
            &calc_apolar_term')

    case default
       call utils_abort('Unrecognized dielectric function in calc_apolar_term')

    end select

    call comms_reduce('SUM',cavity_surface_area)
    call comms_reduce('SUM',cavity_volume)

    ! Now normalize the surface area and the volume
    cavity_surface_area = cavity_surface_area * grid%weight
    cavity_volume = cavity_volume * grid%weight

    ! Update the problem
    select case(pub_is_apolar_method)
    case('SASA')

       problem%solvation_terms%E_apolar_cavitation = &
            cavity_surface_area * pub_is_solvent_surf_tension
       problem%solvation_terms%E_apolar_disrep = &
            cavity_surface_area * &
            pub_is_solvent_surf_tension * (-(1D0-pub_is_apolar_scaling_factor))

       problem%cavity_surface_area = cavity_surface_area
       problem%cavity_volume = cavity_volume

    case('SAV')

       problem%solvation_terms%E_apolar_cavitation = &
            (cavity_surface_area * pub_is_solvent_surf_tension) + &
            (cavity_volume * pub_is_solvent_pressure)
       problem%solvation_terms%E_apolar_disrep = &
            cavity_surface_area * &
            pub_is_solvent_surf_tension * (-(1D0-pub_is_apolar_scaling_factor))
       problem%solvation_terms%E_apolar_sasa = &
            (cavity_surface_area * pub_is_solvent_surf_tension)
       problem%solvation_terms%E_apolar_sav = &
            (cavity_volume * pub_is_solvent_pressure)

       problem%cavity_surface_area = cavity_surface_area
       problem%cavity_volume = cavity_volume

    case default
       call utils_abort('Unrecognized apolar method')
    end select

    ! jd: Clean up
    if(want_scf) then
       call utils_sanity_check(problem%V_apolar(:,:,:),'V_apolar in solvation')
       deallocate(second_der,stat=ierr)
       call utils_dealloc_check(myself,'second_der',ierr)
    end if

    if (pub_is_dielectric_function=='SOFT_SPHERE') then
       deallocate(grad_eps,stat=ierr)
       call utils_dealloc_check(myself,'grad_eps',ierr)
    else
       deallocate(grad_den,stat=ierr)
       call utils_dealloc_check(myself,'grad_den',ierr)
    end if

    call utils_trace_out(myself)

  end subroutine calc_apolar_term
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine calc_sasa_sav_isodensity(grid, problem, grad_den, second_der, &
       & cavity_surface_area, cavity_volume, want_scf)
    !=========================================================================!
    ! Calculates the surface area by taking the volume between two nearby     !
    ! isodensity values (ie. volume has a 'shell' of charge) and dividing     !
    ! by the local thickness of this charge film or shell - giving the surface!
    ! area of the associated dielectric cavity. The volume is calculated as a !
    ! simple sum of the dielectric cavity function along the isodensity.      !
    !-------------------------------------------------------------------------!
    ! Moved from calc_apolar and placed here to aid readability               !
    ! Added by Gabriel Bramley, 05/21                                         !
    !=========================================================================!
    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use rundat, only: pub_is_density_threshold, pub_is_surface_thickness, &
         pub_is_solvent_surf_tension, pub_is_apolar_scaling_factor

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(POISSON_PROBLEM), intent(inout) :: problem
    real(kind=DP), intent(in),dimension(:,:,:,:), allocatable :: grad_den
    real(kind=DP), intent(in),dimension(:,:,:,:,:), allocatable :: second_der
    real(kind=DP), intent(out) :: cavity_surface_area
    real(kind=DP), intent(out) :: cavity_volume
    logical, intent(in) :: want_scf
    ! Internal variables
    real(kind=DP) :: lcl_mod_grad_den
    real(kind=DP) :: lcl_den, lcl_den_a, lcl_den_b
    ! Density thresholds which are used to define isosurfaces
    real(kind=DP) :: dthresh_a,dthresh_b
    real(kind=DP) :: smooth_a,smooth_b,smooth_0    !smoothed density isosurfaces
    integer       :: i1, i2, i3
    real(kind=DP) :: V_apolar_here
    integer       :: i, j ! jd: x, y or z

    ! Use density_threshold and surface_thickness parameters to define two
    ! isosurfaces for finite difference calculation

    dthresh_a = pub_is_density_threshold - 0.5_DP * pub_is_surface_thickness
    dthresh_b = pub_is_density_threshold + 0.5_DP * pub_is_surface_thickness

    cavity_surface_area = 0.0_DP
    cavity_volume = 0.0_DP

    do i3=1, grid%num_my_slabs12
       do i2=1, grid%n2
          do i1=1, grid%n1

             ! Calculate local value of |grad den|
             lcl_mod_grad_den = sqrt(grad_den(i1,i2,i3,1)**2 + &
                  grad_den(i1,i2,i3,2)**2 + grad_den(i1,i2,i3,3)**2)

             ! And local value of density
             lcl_den = problem%rho_elec(i1,i2,i3)

             smooth_a = smoothing_function_fattebert_gygi(lcl_den, dthresh_a)
             smooth_b = smoothing_function_fattebert_gygi(lcl_den, dthresh_b)
             smooth_0 = smoothing_function_fattebert_gygi(lcl_den, &
                  pub_is_density_threshold)

             cavity_surface_area = cavity_surface_area + (smooth_a-smooth_b) * &
                  lcl_mod_grad_den / pub_is_surface_thickness

             if(want_scf) then
                V_apolar_here = 0D0

                ! jd: Take into account only points where the magnitude of the
                !     gradient is not infinitesimally small. This avoids trouble
                !     in the corner case where a centre of a symmetrical mole-
                !     cule lies exactly on a grid point, yielding something like
                !     1E-12 in the gradient, which then makes V_apolar_here
                !     explode. The term should then be zero.
                if(abs(lcl_mod_grad_den) > 1D-8) then
                   do i=1, 3
                      do j=1, 3
                         V_apolar_here = V_apolar_here + &
                              (grad_den(i1,i2,i3,j) * grad_den(i1,i2,i3,i) * &
                                 second_der(i1,i2,i3,j,i)) / &
                                 (lcl_mod_grad_den * lcl_mod_grad_den)
                      end do
                      V_apolar_here = V_apolar_here - second_der(i1,i2,i3,i,i)
                   end do

                   V_apolar_here = V_apolar_here / lcl_mod_grad_den * &
                        (smooth_a-smooth_b) * &
                        pub_is_solvent_surf_tension * &
                        pub_is_apolar_scaling_factor / pub_is_surface_thickness
                else
                      ! jd: No-op -- if lcl_mod_grad_den is zero, we avoid /0
                      !     by letting V_apolar_here remain zero.
                end if

                problem%V_apolar(i1,i2,i3) = V_apolar_here

             end if

             cavity_volume = cavity_volume + smooth_0

          end do
       end do
    end do

  end subroutine calc_sasa_sav_isodensity

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine calc_sasa_sav_density(grid, problem, grad_den, second_der, &
       & cavity_surface_area, cavity_volume, want_scf)
    !=========================================================================!
    ! Calculates the surface area by taking the volume between two nearby     !
    ! density values (giving the volume of a 'shell' of charge) and dividing  !
    ! by the local thickness of this charge film or shell - giving the surface!
    ! area of the associated dielectric cavity. The volume is calculated as a !
    ! simple sum of the dielectric cavity function along the isodensity.      !
    !-------------------------------------------------------------------------!
    ! Moved from calc_apolar and placed here to aid readability               !
    ! Added by Gabriel Bramley, 05/21                                         !
    !=========================================================================!
    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use rundat, only: pub_is_density_threshold, pub_is_surface_thickness, &
         pub_is_solvent_surf_tension, pub_is_apolar_scaling_factor

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(POISSON_PROBLEM), intent(inout) :: problem
    real(kind=DP), intent(out) :: cavity_surface_area
    real(kind=DP), intent(out) :: cavity_volume
    real(kind=DP), intent(in),dimension(:,:,:,:), allocatable :: grad_den
    real(kind=DP), intent(in),dimension(:,:,:,:,:), allocatable :: second_der
    logical, intent(in) :: want_scf
    ! Internal variables
    real(kind=DP) :: lcl_mod_grad_den
    real(kind=DP) :: lcl_den, lcl_den_a, lcl_den_b
    ! Density thresholds which are used to define isosurfaces
    real(kind=DP) :: dthresh_a,dthresh_b
    real(kind=DP) :: smooth_a,smooth_b,smooth_0    !smoothed density isosurfaces
    integer       :: i1, i2, i3
    real(kind=DP) :: V_apolar_here
    integer       :: i, j ! jd: x, y or z

    ! Use lcl_density and surface_thickness parameters to define two
    ! isosurfaces for finite difference calculation

    cavity_surface_area = 0.0_DP
    cavity_volume = 0.0_DP

    do i3=1, grid%num_my_slabs12
       do i2=1, grid%n2
          do i1=1, grid%n1

             ! Calculate local value of |grad den|
             lcl_mod_grad_den = sqrt(grad_den(i1,i2,i3,1)**2 + &
                  grad_den(i1,i2,i3,2)**2 + grad_den(i1,i2,i3,3)**2)

             ! And local value of density
             lcl_den = problem%rho_elec(i1,i2,i3)
             lcl_den_a = problem%rho_elec(i1,i2,i3) + (0.5_DP * pub_is_surface_thickness)
             lcl_den_b = problem%rho_elec(i1,i2,i3) - (0.5_DP * pub_is_surface_thickness)

             smooth_a = smoothing_function_fattebert_gygi(lcl_den_a, &
                  pub_is_density_threshold)
             smooth_b = smoothing_function_fattebert_gygi(lcl_den_b, &
                  pub_is_density_threshold)
             smooth_0 = smoothing_function_fattebert_gygi(lcl_den, &
                  pub_is_density_threshold)

             cavity_surface_area = cavity_surface_area + (smooth_a-smooth_b) * &
                  lcl_mod_grad_den / pub_is_surface_thickness

             if(want_scf) then
                V_apolar_here = 0D0

                ! jd: Take into account only points where the magnitude of the
                !     gradient is not infinitesimally small. This avoids trouble
                !     in the corner case where a centre of a symmetrical mole-
                !     cule lies exactly on a grid point, yielding something like
                !     1E-12 in the gradient, which then makes V_apolar_here
                !     explode. The term should then be zero.
                if(abs(lcl_mod_grad_den) > 1D-8) then
                   do i=1, 3
                      do j=1, 3
                         V_apolar_here = V_apolar_here + &
                              (grad_den(i1,i2,i3,j) * grad_den(i1,i2,i3,i) * &
                                 second_der(i1,i2,i3,j,i)) / &
                                 (lcl_mod_grad_den * lcl_mod_grad_den)
                      end do
                      V_apolar_here = V_apolar_here - second_der(i1,i2,i3,i,i)
                   end do

                   V_apolar_here = V_apolar_here / lcl_mod_grad_den * &
                        (smooth_a-smooth_b) * &
                        pub_is_solvent_surf_tension * &
                        pub_is_apolar_scaling_factor / pub_is_surface_thickness
                else
                      ! jd: No-op -- if lcl_mod_grad_den is zero, we avoid /0
                      !     by letting V_apolar_here remain zero.
                end if

                problem%V_apolar(i1,i2,i3) = V_apolar_here

             end if

             cavity_volume = cavity_volume + smooth_0

          end do
       end do
    end do

  end subroutine calc_sasa_sav_density

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine calc_sasa_sav_andreussi(grid, problem, grad_den, second_der, &
       & cavity_surface_area, cavity_volume, want_scf)
    !=========================================================================!
    ! Calculates the surface area by taking the volume between two nearby     !
    ! density values (giving the volume of a 'shell' of charge) and dividing  !
    ! by the local thickness of this charge film or shell - giving the surface!
    ! area of the associated dielectric cavity. The volume is calculated as a !
    ! simple sum of the dielectric cavity function along the isodensity.      !
    ! This uses the Andreussi approach, which works like the ONETEP's         !
    ! 'density' approach, only it uses Andreussi's smoothing function.        !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2021.                                      !
    !=========================================================================!
    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use rundat, only: pub_is_density_threshold, pub_is_surface_thickness, &
         pub_is_solvent_surf_tension, pub_is_apolar_scaling_factor

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(POISSON_PROBLEM), intent(inout) :: problem
    real(kind=DP), intent(out) :: cavity_surface_area
    real(kind=DP), intent(out) :: cavity_volume
    real(kind=DP), intent(in),dimension(:,:,:,:), allocatable :: grad_den
    real(kind=DP), intent(in),dimension(:,:,:,:,:), allocatable :: second_der
    logical, intent(in) :: want_scf
    ! Internal variables
    real(kind=DP) :: lcl_mod_grad_den
    real(kind=DP) :: lcl_den, lcl_den_a, lcl_den_b
    ! Density thresholds which are used to define isosurfaces
    real(kind=DP) :: dthresh_a,dthresh_b
    real(kind=DP) :: smooth_a,smooth_b,smooth_0    !smoothed density isosurfaces
    integer       :: i1, i2, i3
    real(kind=DP) :: V_apolar_here
    integer       :: i, j ! jd: x, y or z

    ! Use lcl_density and surface_thickness parameters to define two
    ! isosurfaces for finite difference calculation

    cavity_surface_area = 0.0_DP
    cavity_volume = 0.0_DP

    do i3=1, grid%num_my_slabs12
       do i2=1, grid%n2
          do i1=1, grid%n1

             ! Calculate local value of |grad den|
             lcl_mod_grad_den = sqrt(grad_den(i1,i2,i3,1)**2 + &
                  grad_den(i1,i2,i3,2)**2 + grad_den(i1,i2,i3,3)**2)

             ! And local value of density
             lcl_den = problem%rho_elec(i1,i2,i3)

             ! jd: Note that Andreussi et al use an improved formula compared
             !     to Fattebert and Gygi.
             ! FG: sm(rho_0-Delta/2,rho) - sm(rho_0+Delta/2,rho) * grad_term
             ! A:  sm(params,rho-Delta/2) - sm(params,rho+Delta/2) * grad_term

             smooth_a = smoothing_function_andreussi(&
                  lcl_den + 0.5_DP * pub_is_surface_thickness)
             smooth_b = smoothing_function_andreussi(&
                  lcl_den - 0.5_DP * pub_is_surface_thickness)
             smooth_0 = smoothing_function_andreussi(lcl_den)

             cavity_surface_area = cavity_surface_area + (smooth_a-smooth_b) * &
                  lcl_mod_grad_den / pub_is_surface_thickness

             if(want_scf) then
                V_apolar_here = 0D0

                ! jd: Take into account only points where the magnitude of the
                !     gradient is not infinitesimally small. This avoids trouble
                !     in the corner case where a centre of a symmetrical mole-
                !     cule lies exactly on a grid point, yielding something like
                !     1E-12 in the gradient, which then makes V_apolar_here
                !     explode. The term should then be zero.
                if(abs(lcl_mod_grad_den) > 1D-8) then
                   do i=1, 3
                      do j=1, 3
                         V_apolar_here = V_apolar_here + &
                              (grad_den(i1,i2,i3,j) * grad_den(i1,i2,i3,i) * &
                                 second_der(i1,i2,i3,j,i)) / &
                                 (lcl_mod_grad_den * lcl_mod_grad_den)
                      end do
                      V_apolar_here = V_apolar_here - second_der(i1,i2,i3,i,i)
                   end do

                   V_apolar_here = V_apolar_here / lcl_mod_grad_den * &
                        (smooth_a-smooth_b) * &
                        pub_is_solvent_surf_tension * &
                        pub_is_apolar_scaling_factor / pub_is_surface_thickness
                else
                      ! jd: No-op -- if lcl_mod_grad_den is zero, we avoid /0
                      !     by letting V_apolar_here remain zero.
                end if

                problem%V_apolar(i1,i2,i3) = V_apolar_here

             end if

             cavity_volume = cavity_volume + smooth_0

          end do
       end do
    end do

  end subroutine calc_sasa_sav_andreussi

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine calc_diel_correction(solv_problem,grid)
    !=========================================================================!
    ! Given the electrostatic potential computed in the presence of a         !
    ! dielectric we compute the contribution to the Kohn-Sham potential       !
    ! arising from the dielectric depending on the charge density. This is    !
    ! necessary to ensure that the charge density and dielectric can respond  !
    ! to each other within the SCF procedure such that we achieve a tunable   !
    ! dielectric cavity for our solvated electronic structure calculation     !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    ! ----------------------------------------------------------------------- !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP, PI
    use is_poisson, only: POISSON_PROBLEM
    use rundat, only: pub_is_dielectric_model
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(POISSON_PROBLEM), intent(inout) :: solv_problem

    ! Internal variables
    real(kind=DP), dimension(:,:,:), allocatable :: deps_drho
    real(kind=DP), parameter :: inv_eight_pi = 1.0_DP/(8.0_DP * PI)
    integer :: ierr

    !------------------------------------------------------------------------

    call utils_trace_in('calc_diel_pot_and_corr')

    call utils_assert(pub_is_dielectric_model == 'SELF_CONSISTENT', &
         'Internal error in is_solvation_mod, calc_diel_correction')

    ! Allocate and get the dielectric medium derivative
    allocate(deps_drho(grid%ld1, grid%ld2, &
         grid%max_slabs12),stat=ierr)
    call utils_alloc_check('calc_diel_pot_and_corr','deps_drho',ierr)

    deps_drho = dielectric_functional_deriv(solv_problem%rho_elec)

    ! jd: Set the derivative to zero inside the cores, in case the electronic
    !     density drops too low because of the pseudos, which would create
    !     a spurious cavity. This is to make it consistent with the epsilon
    !     which was set to zero in these regions.
    call exclude_cavity(deps_drho,0.0_DP,grid)

    solv_problem%V_eps = -inv_eight_pi * deps_drho

    deallocate(deps_drho,stat=ierr)
    call utils_dealloc_check('calc_diel_pot_and_corr','deps_drho',ierr)

    ! Need to make sure the global grad_phi_sqd is set up
    call get_grad_phi_sqd(solv_problem,grid)

    solv_problem%V_eps = solv_problem%V_eps * solv_problem%grad_phi_sqd

    call utils_trace_out('calc_diel_pot_and_corr')

  end subroutine calc_diel_correction
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine get_grad_phi_sqd(solv_problem,grid)
    !=========================================================================!
    ! Compute the gradient of the solvation potential and store in global     !
    ! array which is needed in both the energy evaluation and the dielectric  !
    ! potential. Uses global flag have_grad_phi_sqd                           !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    ! ----------------------------------------------------------------------- !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use finite_differences, only: finite_difference_mod_grad_sqd
    use multigrid_methods, only: mg
    use is_poisson, only: POISSON_PROBLEM
    use rundat, only: pub_finite_difference_order

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(POISSON_PROBLEM), intent(inout) :: solv_problem

    !------------------------------------------------------------------------

    call utils_trace_in('get_grad_phi_sqd')

    if (.not. solv_problem%have_grad_phi_sqd) then
       ! jd: Get the gradient
       call finite_difference_mod_grad_sqd(solv_problem%grad_phi_sqd, &
            solv_problem%phi_mol, pub_finite_difference_order, &
            mg, grid)

       ! Lastly, set the flag so we don't have to repeat this process again
       solv_problem%have_grad_phi_sqd = .true.

    end if

    call utils_trace_out('get_grad_phi_sqd')

  end subroutine get_grad_phi_sqd
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine check_elec_energy_gradient(initial,grid,cell,elements)
    !=========================================================================!
    ! Takes the given charge density and compares the numerical electrostatic !
    ! energy gradient to the analytic.                                        !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   initial, intent(inout): the initial POISSON_PROBLEM.                  !
    ! ----------------------------------------------------------------------- !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: DP, VERBOSE, stdout
    use is_poisson, only: POISSON_PROBLEM, allocate_poisson_problem, &
         deallocate_poisson_problem
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_output_detail, pub_is_dielectric_model, &
         pub_is_solvation_method, pub_is_include_apolar
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    use ion, only: ELEMENT

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(POISSON_PROBLEM), intent(inout) :: initial
    type(ELEMENT), intent(in) :: elements(:)

    ! Internal variables
    type(POISSON_PROBLEM) :: varied
    real(kind=DP), allocatable, dimension(:,:,:) :: gradient
    real(kind=DP)                   :: actual,predicted,ratio,Q_varied
    integer :: ierr

    real(kind=DP), parameter :: alpha = 1.0D-5 ! jd: step length

    !------------------------------------------------------------------------

    call utils_trace_in('check_elec_energy_gradient')

    ! Allocate memory
    call allocate_poisson_problem(varied,pub_is_solvation_method,grid)
    allocate(gradient(grid%ld1, grid%ld2, &
         grid%max_slabs12),stat=ierr)
    call utils_alloc_check('check_elec_energy_gradient','gradient',ierr)

    ! Now copy functional derivative into gradient (V_eps = 0 for fixed eps)
    gradient = initial%phi_mol

    if(pub_is_dielectric_model == 'SELF_CONSISTENT') then
       gradient = gradient + initial%V_eps
       if(pub_is_include_apolar) then
          gradient = gradient + initial%V_apolar
       end if
    end if

    ! Now use the gradient to define a new density
    varied%rho_elec = initial%rho_elec + alpha * gradient
    varied%rho_ion = initial%rho_ion ! direct_solvation_potential relies on this

    ! Correctly evaluate the dielectric
    select case (pub_is_dielectric_model)
    case('FIX_INITIAL')
       ! hhh: Use the same eps as initial problem
       varied%eps_full = initial%eps_full
       varied%eps_half = initial%eps_half
    case('SELF_CONSISTENT') ! Must evaluate the new eps for the varied problem
       call calc_dielectric_medium(varied%rho_elec,varied%eps_full,grid,&
            cell,elements)
       call interpolate_dielectric_medium(varied%rho_elec,varied%eps_half,grid,&
            cell, elements)
    case default
       call utils_abort('Unrecognized pub_is_dielectric_model')
    end select

    ! And now we need to solve the varied solvation problem
    select case(pub_is_solvation_method)
    case ('DIRECT')
       ! Directly solve the Poisson equation in dielectric to find the solv. pot
       call direct_solvation_potential(varied,grid,cell)
    case ('CORRECTIVE')
       ! Use the corrective potential method to get the solvation potential
       call corrective_solvation_potential(varied)
    case default
      call utils_abort('Unrecognized pub_is_solvation_method')
    end select

    ! Compute the numerical energy gradient
    actual     = (varied%E_Hartree - initial%E_Hartree)/alpha

    ! Now we calculate the analytic gradient which is just
    ! int gradient * gradient dr

    ! We also calculate the total charge in rho_varied while we are at it
    predicted = integrals_product_on_grid(grid, &
         gradient, gradient)
    Q_varied = integrals_product_on_grid(grid, &
         varied%rho_elec)

    ! Calculate the ratio actual/predicted ('badness')
    ratio = actual/predicted

    ! Output
    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
      write(stdout,*) 'Electrostatic energy gradient check:'
      write(stdout,*) '  Q_varied      ', Q_varied
      write(stdout,*) '  Step length   ', alpha
      write(stdout,*) '  E_varied      ', varied%E_Hartree
      write(stdout,*) '  E_initial     ', initial%E_Hartree
      write(stdout,*) '  predicted     ', predicted
      write(stdout,*) '  actual        ', actual
      write(stdout,*) '  badness       ', ratio
    end if

    ! Store energy gradient in initial solvation problem
    initial%energy_gradient_ratio_rho = ratio
    initial%energy_gradient_ratio_eps = 0.0_DP ! jd: Removed functionality

    ! Remember
    initial%have_energy_gradient=.true.

    ! Clean up
    call deallocate_poisson_problem(varied,pub_is_solvation_method)
    deallocate(gradient,stat=ierr)
    call utils_dealloc_check('check_elec_energy_gradient','gradient',ierr)

    call utils_trace_out('check_elec_energy_gradient')

  end subroutine check_elec_energy_gradient
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine write_solvation_data(solv_problem,grid,cell)
    !=========================================================================!
    ! Writes out grid data for the solvation problem.                         !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    ! ----------------------------------------------------------------------- !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use is_poisson, only: POISSON_PROBLEM, write_poisson_problem
    use rundat, only: pub_is_dielectric_model, pub_is_solvation_method, &
         pub_is_solvation_output_detail
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(POISSON_PROBLEM), intent(in)  :: solv_problem

    !------------------------------------------------------------------------

    call utils_trace_in('write_solvation_data')

    ! Write out all potentials associated with solvation current problem
    if(pub_is_solvation_output_detail /= 'NONE') then
       call write_poisson_problem(solv_problem,'current', &
            pub_is_solvation_method, pub_is_dielectric_model,grid,cell)

       ! Also write out initial data in initial_problem if we have
       ! a fixed eps or gaussian ions
       if(have_initial_eps) then
          call write_poisson_problem(initial_problem,'initial','INITIAL', &
               'NONE',grid,cell)
       end if
    end if

    call utils_trace_out('write_solvation_data')

  end subroutine write_solvation_data
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



  !===========================================================================!
  !                           Private functions                               !
  !===========================================================================!

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elemental function eps_function_fattebert_gygi(rho_r)
    !=========================================================================!
    ! Smooth dielectric function, proposed by Jean-Luc Fattebert and          !
    ! Francois Gygi. See J Comp Chem 23 pgs 662-666.                          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rho_r,            intent=in,  the electronic charge density at r      !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal      07/09/2007                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_bulk_permittivity, pub_is_density_threshold, &
         pub_is_solvation_beta

    implicit none

    ! Return type
    real(kind=DP) :: eps_function_fattebert_gygi

    ! Arguments
    real(kind=DP), intent(in) :: rho_r

    ! Internal variables
    real(kind=DP) :: powbeta, rho

    !------------------------------------------------------------------------

    rho = rho_r
    if(rho < 0.0_DP) rho = 0.0_DP ! jd: Filter out ringing and noise

    powbeta = (rho/pub_is_density_threshold)**(2.0_DP*pub_is_solvation_beta)

    eps_function_fattebert_gygi = &
         1.0_DP + ((pub_is_bulk_permittivity-1.0_DP)/2.0_DP) * &
         ( 1.0_DP + (1.0_DP-powbeta)/(1.0_DP+powbeta) )

  end function eps_function_fattebert_gygi
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elemental function eps_function_gaussian(rho_r)
    !=========================================================================!
    ! Smooth, Gaussian-shaped dielectric function depending on only one param.!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rho_r,            intent=in,  the electronic charge density at r      !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal      21/09/2009                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_density_threshold, pub_is_bulk_permittivity

    implicit none

    ! Return type
    real(kind=DP) :: eps_function_gaussian

    ! Arguments
    real(kind=DP), intent(in) :: rho_r

    ! Internal variables
    real(kind=DP) :: rho, rho_frac

    !------------------------------------------------------------------------

    rho = rho_r
    if(rho < 0.0_DP) rho = 0.0_DP ! jd: Filter out ringing and noise

    rho_frac = rho/pub_is_density_threshold

    eps_function_gaussian = 1.0_DP + &
         (pub_is_bulk_permittivity-1.0_DP) * exp(-rho_frac**2.0_DP)

  end function eps_function_gaussian
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elemental function eps_function_andreussi(rho_r)
    !=========================================================================!
    ! Smooth permittivity function proposed by Andreussi et al [2].           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rho_r,            intent=in,  the electronic charge density at r      !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2021.08 by adapting                        !
    ! eps_function_fattebert_gygi().                                          !
    !=========================================================================!

    use constants, only: DP, PI
    use rundat, only: pub_is_bulk_permittivity, pub_is_density_min_threshold, &
         pub_is_density_max_threshold

    implicit none

    ! Return type
    real(kind=DP) :: eps_function_andreussi

    ! Arguments
    real(kind=DP), intent(in) :: rho_r

    ! Internal variables
    real(kind=DP) :: k, y, z, t_x

    !--------------------------------------------------------------------------

    ! jd: Do the corner cases first to avoid log(0) or log(-ve) from ringing.
    if(rho_r > pub_is_density_max_threshold) then
       eps_function_andreussi = 1.0_DP
    else if(rho_r < pub_is_density_min_threshold) then
       eps_function_andreussi = pub_is_bulk_permittivity
    else
       k = 2.0_DP * PI * ((log(pub_is_density_max_threshold) - log(rho_r)) / &
            (log(pub_is_density_max_threshold) - &
            log(pub_is_density_min_threshold)))
       y = sin(2.0_DP * PI * ((log(pub_is_density_max_threshold) - log(rho_r)) / &
            (log(pub_is_density_max_threshold) - &
            log(pub_is_density_min_threshold))))
       z = (log(pub_is_bulk_permittivity))/(2.0_DP * PI)
       t_x = z*(k-y)
       eps_function_andreussi = exp(t_x)
    end if


  end function eps_function_andreussi
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elemental function dielectric_functional_deriv(rho_r)
    !=========================================================================!
    ! Computes the functional derivative of the dielectric                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rho_r, intent=in,  the electronic charge density at r                 !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal 2/11/2007                                      !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_density_threshold, pub_is_bulk_permittivity, &
         pub_is_solvation_beta

    implicit none

    ! Return type
    real(kind=DP) :: dielectric_functional_deriv

    ! Arguments
    real(kind=DP), intent(in) :: rho_r

    ! Internal variables
    real(kind=DP)             :: num,denom,rho,rho_frac,twobeta

    !------------------------------------------------------------------------

    rho = rho_r
    if(rho < 0.0_DP) rho = 0.0_DP ! jd: Filter out ringing and noise

    rho_frac = rho/pub_is_density_threshold

    twobeta  = 2.0_DP*pub_is_solvation_beta

    num = (1.0_DP-pub_is_bulk_permittivity)*(twobeta*rho_frac**(twobeta-1.0_DP))

    denom = pub_is_density_threshold*(1.0_DP+rho_frac**twobeta)**2.0_DP

    dielectric_functional_deriv = num/denom

  end function dielectric_functional_deriv
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  function smoothing_function_fattebert_gygi(rho_r,thresh)
    !=========================================================================!
    ! Smooth threshold function used for computing cavity volume and surface  !
    ! in the Fattebert-Gygi model, cf. [1].                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rho_r,            intent=in,  the electronic charge density at r      !
    !   thresh,           intent=in,  the threshold parameter                 !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal      24/08/2007                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_solvation_beta

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: rho_r
    real(kind=DP), intent(in) :: thresh

    ! Internal variables
    real(kind=DP) :: powbeta, rho

    ! Return type
    real(kind=DP) :: smoothing_function_fattebert_gygi

    !------------------------------------------------------------------------

    rho = rho_r
    if(rho < 0.0_DP) rho = 0.0_DP ! jd: Filter out ringing and noise

    powbeta = (rho/thresh)**(2.0_DP*pub_is_solvation_beta)

    smoothing_function_fattebert_gygi = &
         0.5_DP*( (powbeta-1.0_DP) / (powbeta+1.0_DP) + 1.0_DP )

  end function smoothing_function_fattebert_gygi
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  function smoothing_function_andreussi(rho_r)
    !=========================================================================!
    ! Smooth threshold function used for computing cavity volume and surface  !
    ! in the Andreussi model, cf. [2], eq. (64).                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rho_r,            intent=in,  the electronic charge density at r.     !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2019.                             !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_bulk_permittivity, pub_is_density_min_threshold, &
         pub_is_density_max_threshold
    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: rho_r

    ! Return type
    real(kind=DP) :: smoothing_function_andreussi

    ! Internal variables
    real(kind=DP) :: rho
    character(len=*), parameter :: myself = 'smoothing_function_andreussi'

    !------------------------------------------------------------------------

    rho = rho_r
    if(rho < 0.0_DP) rho = 0.0_DP ! jd: Filter out ringing and noise

    if(pub_is_bulk_permittivity <= 1.0_DP) then
       call utils_abort(myself//': Illegal is_bulk_permittivity. This model &
            &requires eps_0 to be > 1, cf. (64) in [2].')
    end if

    smoothing_function_andreussi = &
         (pub_is_bulk_permittivity - eps_function_andreussi(rho)) / &
         (pub_is_bulk_permittivity - 1.0_DP)

  end function smoothing_function_andreussi
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine exclude_cavity(array, exclusion_value, grid)
    !==========================================================================!
    ! Allows the implicit solvent dielectric to be fully or partially excluded !
    ! from certain gridpoints as defined by dielectric_mask. This allows the   !
    ! user to define exclusion regions with the is_dielectric_exclusions       !
    ! block. These regions can be softened with a fermi smearing to prevent    !
    ! discontinuities in the dielectric by setting a small (~several times the !
    ! multigrid spacing) is_dielectric_exclusions_smearing in the input.       !
    !                                                                          !
    ! Also ensures that the dielectric is strictly 1.0 inside the cores. Since !
    ! we use pseudopotentials, the electronic density occasionally drops so    !
    ! low in the cores that a spurious cavity is created. We thus artificially !
    ! set the permittivity to 1 inside the core and the derivative to 0.       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   array, intent=inout:          the array to fix.                        !
    !   exclusion_value, intent=in:   the value to set                         !
    !                                 inside the excluded region.              !
    !--------------------------------------------------------------------------!
    ! Note:                                                                    !
    !   smeared_ion_generate_density() must have been called prior to that     !
    !   subroutine, as it fills in dielectric_mask.                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 04/2012, as "ensure_no_cavity_inside_cores".  !
    ! Subroutine renamed to "exclude_cavity" in March 2018, after              !
    ! Kevin K. B. Duff and Edward B. Linscott changed the dielectric mask type !
    ! from logical to real.                                                    !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use rundat, only: pub_is_core_width
    use is_smeared_ions, only: dielectric_mask

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)  :: grid
    real(kind=DP), intent(inout) :: array(grid%ld1, grid%ld2, grid%max_slabs12)
    real(kind=DP)                :: exclusion_value

    ! Local variables
    integer :: i1, i2, islab12

    !------------------------------------------------------------------------

    if (pub_is_core_width > 0.0_DP) then

       do islab12 = 1, grid%num_my_slabs12
          do i2 = 1, grid%n2
             do i1 = 1, grid%n1

                ! kkbd: Apply the dielectric mask appropriately
                ! now that it's a real
                array(i1,i2,islab12) = &
                     array(i1,i2,islab12) &
                     * (1.0_DP - dielectric_mask(i1,i2,islab12)) &
                     + dielectric_mask(i1,i2,islab12) * exclusion_value

             end do ! over i1
          end do ! over i2
       end do ! over 12-slabs

    end if
  end subroutine exclude_cavity

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine calculate_eps_soft_sphere(eps, elements, grid, cell, &
       is_periodic, shift_dim)

    !==========================================================================!
    ! Computes the dielectric cavity in terms of the soft-sphere cavity model  !
    ! as proposed by Fisicaro et al. [1]. This module provides an alternative  !
    ! to the isoelectronic dielectric cavity model of Fattebert and Gygi       !
    ! and can tune the cavity size for individual elements.                    !
    !                                                                          !
    ! \eps(r_i,R_i) = (\eps_0 - 1)*(\prod_i (h_i(||r-R_i||))) + 1              !
    !                                                                          !
    ! The dielectric cavity is defined by the product of distance functions,   !
    ! h_i for a particular point with respect to each atom, i, in the system.  !
    !                                                                          !
    ! This module is largely based on the is_smeared_density_create subroutine !
    ! as implemented in the is_smeared_ions_mod module by Jacek Dziedzic.      !
    !                                                                          !
    ! 1. Fisicaro, G. et al. J. Chem. Theory Comput. 13, 3829-3845 (2017).     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! elements(par%nat), intent(in): Provides the number of atoms to be looped !
    ! over.                                                                    !
    ! grid, intent(in): Points of the grid.                                    !
    ! cell, intent(in): Information for periodic distance.                     !
    ! is_periodic(3), intent(in): Determines whether periodic boundary         !
    ! need be applied to the distance function.                                !
    ! shift_dim, intent(in): Crude way of constructing the half grid required  !
    ! for the multigrid solver. Values 1,2 and 3 shift the grid by half a grid !
    ! space in x, y and z respectively, and 0 keeps the default grid.          !
    ! eps, intent(out): Returns array containing the dielectric values.        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic as smeared_ion_generate_density, 05/2010       !
    ! Adapted for soft spheres by Gabriel Bramley, 01/2019                     !
    !==========================================================================!

    use geometry, only              : POINT, geometry_distance, magnitude
    use ion, only                   : ELEMENT
    use cell_grid, only             : GRID_INFO
    use simulation_cell, only       : CELL_INFO, minimum_image_distance
    use parallel_strategy, only     : par=>pub_par
    use comms, only                 : pub_my_proc_id
    use constants, only             : DP
    use utils, only                 : utils_abort
    use rundat, only                : pub_is_dielectric_model, pub_threads_max

    implicit none

    ! gab: Arguments

    type(ELEMENT), intent(in)        :: elements(par%nat)
    type(GRID_INFO), intent(in)      :: grid
    type(CELL_INFO), intent(in)      :: cell
    logical, intent(in)              :: is_periodic(3)
    integer, intent(in)              :: shift_dim
    real(kind=DP), dimension(grid%ld1, &
         grid%ld2, grid%max_slabs12),intent(out) :: eps


    ! gab: Local variables and parameters

    real(kind=DP)                    :: h_i, sum_hi, dist, X, Y, Z, r_ss, &
         magx, magy, magz
    integer                          :: i1, i2, i3, islab12, ipt, p, atom_p
    type(POINT)                      :: R_I, R

    ! gab: Makes sure the dielectric model is set to fix_initial for the soft
    ! sphere model. Doesn't make sense to change cavity self-consistently when
    ! cavity is based on geometric factors alone.

    if (pub_is_dielectric_model == 'SELF_CONSISTENT') then
       call utils_abort('Soft-sphere cavity method is incompatible with &
            &self-consistent cavity generation.')
    end if

    ! gab: Sloppy way of shifting the grid onto half grids.

    select case (shift_dim)
    case(1)
       magx=0.5_DP
       magy=0.0_DP
       magz=0.0_DP
    case(2)
       magx=0.0_DP
       magy=0.5_DP
       magz=0.0_DP
    case(3)
       magx=0.0_DP
       magy=0.0_DP
       magz=0.5_DP
    case(0)
       magx=0.0_DP
       magy=0.0_DP
       magz=0.0_DP
    end select

    ! gab: INITIAL TESTING FOR VDW NOT PRESENT IN THE DATABASE
    do p = 1, par%nat
       r_ss = elements(p)%soft_sphere_radius
       if (r_ss.gt.100.0) then
          call utils_abort('Element: '//elements(p)%symbol//' does not &
               &yet have a soft sphere radius available. Sorry!')
       end if
    end do

    ! gab: Implemented for both sets of boundary conditions

    if (.not.any(is_periodic)) then

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,i3,R,R_I,dist,X,Y,Z,r_ss,    &
!$OMP       atom_p, h_i,p,sum_hi) &
!$OMP SHARED(cell,elements,grid,is_periodic,magx,magy,magz,  &
!$OMP     par,pub_my_proc_id,pub_threads_max,eps)

       do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / &
               (grid%n1*grid%n2) + 1
          i3 = grid%first_slab12(pub_my_proc_id) + islab12 - 1

          R%X = (magnitude(grid%da1)*magx)+((i1-1)*magnitude(grid%da1))
          R%Y = (magnitude(grid%da2)*magy)+((i2-1)*magnitude(grid%da2))
          R%Z = (magnitude(grid%da3)*magz)+((i3-1)*magnitude(grid%da3))

          ! gab: Declare value of sum_hi function.

          sum_hi = 1.0_DP

          do p = 1, par%nat

             R_I  = elements(p)%centre

             dist = geometry_distance(R_I, R)

             atom_p = elements(p)%atomic_number
             r_ss = elements(p)%soft_sphere_radius

             h_i = h_distance_function(dist,r_ss)

             sum_hi  = sum_hi * h_i

          end do

          eps(i1,i2,islab12) = eps_function_soft_sphere(sum_hi)

       end do

!$OMP END PARALLEL DO

    else if (all(is_periodic)) then

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,i3,R,R_I,dist,X,Y,Z,r_ss,   &
!$OMP       atom_p, h_i,p,sum_hi) &
!$OMP SHARED(cell,elements,grid,is_periodic,magx,magy,magz,  &
!$OMP     par,pub_my_proc_id,pub_threads_max,eps)

       do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / &
               (grid%n1*grid%n2) + 1
          i3 = grid%first_slab12(pub_my_proc_id) + islab12 - 1

          R%X = (magnitude(grid%da1)*magx)+((i1-1)*magnitude(grid%da1))
          R%Y = (magnitude(grid%da2)*magy)+((i2-1)*magnitude(grid%da2))
          R%Z = (magnitude(grid%da3)*magz)+((i3-1)*magnitude(grid%da3))

          ! Declare initial value of product

          sum_hi = 1.0_DP

          ! Calculate h function for each atom at the grid point.

          do p = 1, par%nat

             R_I  = elements(p)%centre
             dist = minimum_image_distance(R_I, R, cell)

             atom_p = elements(p)%atomic_number
             r_ss = elements(p)%soft_sphere_radius

             h_i     = h_distance_function(dist,r_ss)

             sum_hi  = sum_hi * h_i

          end do

          eps(i1,i2,islab12) = eps_function_soft_sphere(sum_hi)

       end do

!$OMP END PARALLEL DO

    end if

  end subroutine calculate_eps_soft_sphere

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  real(kind=DP) function eps_function_soft_sphere(h_i)

    !=========================================================================!
    ! Calculates the dielectric permittivity at a point as a function of      !
    ! distance as a smooth transition from 1 to the bulk value                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   h_i, intent(in), the distance function                                !
    !-------------------------------------------------------------------------!
    ! Written by Gabriel Bramley, 01/2019                                     !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_bulk_permittivity

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: h_i

    ! Internal variables
    real(kind=DP) :: h, delta

    h = h_i

    eps_function_soft_sphere = &
         ( ( pub_is_bulk_permittivity - 1.0_DP ) * ( h ) ) + 1.0_DP

  end function eps_function_soft_sphere

  real(kind=DP) function h_distance_function(dist, soft_sphere_radius)

    !=========================================================================!
    ! Calculates a value, used to represent the smooth variation of the       !
    ! dielectric function from the atom centre to the bulk value.             !
    !                                                                         !
    !     h(r_i, delta) = 0.5_DP * (1 + erf(||r - R_i|| - r_ss/delta))        !
    !                                                                         !
    ! This is achieved through an error function, taking the distance         !
    ! between the point on the grid and the atomic centre, while the radius   !
    ! is tuned for specific elements by the input Van Der Waals radius. The   !
    ! smoothness of this function is tuned by value of delta, with higher     !
    ! values decreasing the steepness of the transition.                      !
    !                                                                         !
    ! This function returns a value 0 =< h_i =< 1                             !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   dist, intent(in): distance from an atomic centre to the grid point.   !
    !   vdw_radius, intent(in): the distance vdw radius of the current atom   !
    !   used to define the cavity size.                                       !
    !-------------------------------------------------------------------------!
    ! Written by Gabriel Bramley, 01/2019                                     !
    !=========================================================================!

    use constants, only    : DP
    use rundat, only       : pub_is_bulk_permittivity
    use utils, only        : utils_erf
    use rundat, only       : pub_is_soft_sphere_delta

    implicit none

    ! Arguments
    real(kind=DP), intent(in)  :: dist
    real(kind=DP), intent(in)  :: soft_sphere_radius

    ! Internal variables
    real(kind=DP)              :: r_i             ! || r - R_i ||
    real(kind=DP)              :: bulk_eps, delta

    ! gab: Set values used in the function.

    r_i = dist

    bulk_eps = pub_is_bulk_permittivity
    delta = pub_is_soft_sphere_delta

    h_distance_function =  0.5_DP*(1.0_DP+utils_erf((r_i-soft_sphere_radius)/delta))

  end function h_distance_function

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module is_solvation
