!=============================================================================!
!                              P O I S S O N                                  !
!=============================================================================!
! This module defines a general data class for collecting and organizing      !
! arrays related to solving the generalized Poisson equation.  Initially      !
! developed for use with implicit solvent but can be used in other modules    !
! for straightforward data output.                                            !
!-----------------------------------------------------------------------------!
! Written by Hatem H Helal, starting in 11/2009                               !
! please report bugs to hhh23(at)cam.ac.uk                                    !
!                                                                             !
! Adapted for onetep by Jacek Dziedzic, 04-05/2010                            !
! Made parallel-ready by Jacek Dziedzic, 04-05/2010                           !
!-----------------------------------------------------------------------------!
module is_poisson

  use constants, only: DP, ANGSTROM
  use utils, only: utils_trace_in, utils_trace_out

  implicit none

  ! Everything is private ...
  private

  !... unless exposed here.
  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!

  public :: allocate_poisson_problem
  public :: deallocate_poisson_problem
  public :: write_poisson_problem
  public :: zero_solvation_terms

  !---------------------------------------------------------------------------!
  !                       P u b l i c   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!

  ! jd: Any changes to this structure must be reflected in zero_solvation_terms().
  type, public :: SOLVATION_ENERGY_TERMS
     ! jd: 1/2 \rho(r) \phi(r) aka 'molecular Hartree'
     real(kind=DP) :: E_elec_fixed

     ! jd: 1/2 \rho(r) \phi(r) aka 'molecular Hartree', in vacuum
     real(kind=DP) :: E_elec_fixed_vac

     ! jd: Apolar solvation energy terms
     real(kind=DP) :: E_apolar_cavitation
     real(kind=DP) :: E_apolar_disrep
     real(kind=DP) :: E_apolar_sasa
     real(kind=DP) :: E_apolar_sav

     ! gab: Stores the surface area and volume to pass to other modules
     real(kind=DP) :: cavity_surface_area
     real(kind=DP) :: cavity_volume

     ! jd: Energy terms present only in Boltzmann solvation. See
     !     implicit_solvent_boltzmann_free_energy_terms() and
     !     implicit_solvent_boltzmann_free_energy_terms_check() for description.
     real(kind=DP) :: E_elec_mob
     real(kind=DP) :: E_osmo
     real(kind=DP) :: E_acc
     real(kind=DP) :: E_atmo
     real(kind=DP) :: E_chempot
     real(kind=DP) :: E_pure
     real(kind=DP) :: E_elec_minus_mob
     real(kind=DP) :: E_minus_ktv_gamma_c_i_bulk
     real(kind=DP) :: E_minus_kt_int_c_i
  end type

  !---------------------------------------------------------------------------!
  ! jd: Holds all data relevant to the solution of the Poisson-Boltzmann      !
  !     equation. All arrays except eps_half represent quantities on the fine !
  !     grid and are stored in usual onetep 12slab-distributed fashion.       !
  !     eps_half holds three such arrays, with the last index indexing the    !
  !     direction in which the grid is shifted by half a spacing.             !
  !---------------------------------------------------------------------------!
  ! jd: Some of the fields in this data structure, denoted with (*), were     !
  !     retained for compatibility with the CASTEP version. They are not used !
  !     currently, but were left over in case the CORRECTIVE and/or FFT       !
  !     approach is implemented, or if write_poisson_problem is ever used.    !
  !---------------------------------------------------------------------------!
  type, public :: POISSON_PROBLEM

     ! jd: Electronic density
     real(kind=DP), allocatable, dimension(:,:,:)       :: rho_elec
     ! jd: (Smeared) ion density
     real(kind=DP), allocatable, dimension(:,:,:)       :: rho_ion
     ! jd: Dielectric permittivity
     real(kind=DP), allocatable, dimension(:,:,:)       :: eps_full
     ! jd: Dielectric permittivity, shifted by half the fine grid spacing
     !     last index picks the direction of the shift
     real(kind=DP), allocatable, dimension(:,:,:,:)     :: eps_half
     ! jd: Potential in PBCs (*)
     real(kind=DP), allocatable, dimension(:,:,:)       :: phi_pbc
     ! jd: Potential due to rho_elec+rho_ion
     real(kind=DP), allocatable, dimension(:,:,:)       :: phi_mol
     ! jd: Corrective potential (*)
     real(kind=DP), allocatable, dimension(:,:,:)       :: phi_corr
     ! jd: Correction to the energy gradient due to cavity changing shape
     !     - the electrostatic term
     real(kind=DP), allocatable, dimension(:,:,:)       :: V_eps
     ! jd: - the apolar term
     real(kind=DP), allocatable, dimension(:,:,:)       :: V_apolar
     ! jd: Temporary for |grad phi|^2
     real(kind=DP), allocatable, dimension(:,:,:)       :: grad_phi_sqd

     ! jd: Hartree energy for the above problem
     real(kind=DP)                                  :: E_Hartree
     ! jd: Apolar terms, and additional energy terms present in Boltzmann solv.
     type(SOLVATION_ENERGY_TERMS)                   :: solvation_terms

     ! jd: Gradient correctness ratio
     real(kind=DP)                                  :: energy_gradient_ratio_rho
     ! jd: (*)
     real(kind=DP)                                  :: energy_gradient_ratio_eps

     ! jd: Surface area of the dielectric cavity
     real(kind=DP)                                  :: cavity_surface_area
     ! jd: Volume of the cavity
     real(kind=DP)                                  :: cavity_volume

     ! jd: .true. if datastructures allocated
     logical                                    :: is_allocated = .false.
     ! jd: .true. if grad_phi_sqd up to date
     logical                                    :: have_grad_phi_sqd = .false.
     ! jd: .true. if energy_gradient_ratio_rho up to date
     logical                                    :: have_energy_gradient= .false.

  end type POISSON_PROBLEM

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine allocate_poisson_problem(problem,method,grid)
    !=========================================================================!
    ! This subroutine will allocate the arrays associated with the            !
    ! poisson_problem data type.  Method argument is used to allocate extra   !
    ! arrays if the corrective method is chosen.                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   problem,        intent=inout,  poisson_problem data type              !
    !   method,         intent=in,     method for solving poisson problem     !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use rundat, only: pub_is_dielectric_model
    use utils, only: utils_alloc_check

    implicit none

    ! jd: Arguments
    type(POISSON_PROBLEM), intent(inout) :: problem
    character(len=*), intent(in) :: method
    type(GRID_INFO), intent(in) :: grid

    ! jd: Internal variables
    integer :: ierr ! jd: Error flag
    integer :: ld1f, ld2f, max_slabs12
    character(len=*), parameter :: myself = 'allocate_poisson_problem'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    ld1f = grid%ld1
    ld2f = grid%ld2
    max_slabs12 = grid%max_slabs12

    allocate(problem%rho_elec(ld1f,ld2f,max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'problem%rho_elec',ierr)

    ! Allow for standard Hartree evaluation
    if(method=='FFT') then
       allocate(problem%phi_pbc(ld1f,ld2f,max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'problem%phi_pbc',ierr)
    else
       allocate(problem%rho_ion(ld1f,ld2f,max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'problem%rho_ion',ierr)
       allocate(problem%eps_full(ld1f,ld2f,max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'problem%eps_full',ierr)
       allocate(problem%eps_half(ld1f,ld2f,max_slabs12,3),stat=ierr)
       call utils_alloc_check(myself,'problem%eps_half',ierr)
       if(method /='INITIAL') then
          allocate(problem%phi_mol(ld1f,ld2f,max_slabs12),stat=ierr)
          call utils_alloc_check(myself,'problem%phi_mol',ierr)
          if(pub_is_dielectric_model == 'SELF_CONSISTENT') then
             allocate(problem%V_eps(ld1f,ld2f,max_slabs12),stat=ierr)
             call utils_alloc_check(myself,'problem%V_eps',ierr)
             allocate(problem%V_apolar(ld1f,ld2f,max_slabs12),stat=ierr)
             call utils_alloc_check(myself,'problem%V_apolar',ierr)
          end if
          allocate(problem%grad_phi_sqd(ld1f,ld2f,max_slabs12),stat=ierr)
          call utils_alloc_check(myself,'problem%grad_phi_sqd',ierr)
          ! Need more arrays if we are using the corrective potential method
          if(method=='CORRECTIVE') then
             allocate(problem%phi_pbc(ld1f,ld2f,max_slabs12),stat=ierr)
             call utils_alloc_check(myself,'problem%phi_pbc',ierr)
             allocate(problem%phi_corr(ld1f,ld2f,max_slabs12),stat=ierr)
             call utils_alloc_check(myself,'problem%phi_corr',ierr)
          end if
       end if
    end if

    problem%is_allocated=.true.

    call utils_trace_out(myself)

  end subroutine allocate_poisson_problem
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine deallocate_poisson_problem(problem,method)
    !=========================================================================!
    ! This subroutine will deallocate the arrays associated with the          !
    ! poisson_problem data type.  Method argument is used to deallocate       !
    ! extra arrays if the corrective method is chosen.                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   problem,        intent=inout,  poisson_problem data type              !
    !   method,         intent=in,     method for solving poisson problem     !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use rundat, only: pub_is_dielectric_model
    use utils, only: utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(POISSON_PROBLEM), intent(inout) :: problem
    character(len=*), intent(in) :: method

    ! jd: Local variables
    integer :: ierr ! jd: Error flag
    character(len=*), parameter :: myself = 'deallocate_poisson_problem'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    !@improveme: should make an attempt to deallocate in reverse order
    if(method=='FFT') then
       deallocate(problem%phi_pbc,stat=ierr)
       call utils_dealloc_check(myself,'problem%phi_pbc',ierr)
    else
       deallocate(problem%rho_ion,stat=ierr)
       call utils_dealloc_check(myself,'problem%rho_ion',ierr)
       deallocate(problem%eps_full,stat=ierr)
       call utils_dealloc_check(myself,'problem%eps_full',ierr)
       deallocate(problem%eps_half,stat=ierr)
       call utils_dealloc_check(myself,'problem%eps_half',ierr)
       if(method /='INITIAL') then
          deallocate(problem%phi_mol,stat=ierr)
          call utils_dealloc_check(myself,'problem%phi_mol',ierr)
          if(pub_is_dielectric_model == 'SELF_CONSISTENT') then
             deallocate(problem%V_apolar,stat=ierr)
             call utils_dealloc_check(myself,'problem%V_apolar',ierr)
             deallocate(problem%V_eps,stat=ierr)
             call utils_dealloc_check(myself,'problem%V_eps',ierr)
          end if
          deallocate(problem%grad_phi_sqd,stat=ierr)
          call utils_dealloc_check(myself,'problem%grad_phi_sqd',&
               ierr)
          ! Need more arrays if we are using the corrective potential method
          if(method=='CORRECTIVE') then
             deallocate(problem%phi_pbc,stat=ierr)
             call utils_dealloc_check(myself,'problem%phi_pbc',&
                  ierr)
             deallocate(problem%phi_corr,stat=ierr)
             call utils_dealloc_check(myself,'problem%phi_corr',&
                  ierr)
          end if
       end if
    end if

    deallocate(problem%rho_elec,stat=ierr)
    call utils_dealloc_check(myself,'problem%rho_elec',ierr)

    problem%is_allocated=.false.

    call utils_trace_out(myself)

  end subroutine deallocate_poisson_problem
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine write_poisson_problem(problem,root,method,model,grid,cell)
    !=========================================================================!
    ! This subroutine will write out a bunch of data files associated with    !
    ! the given poisson problem.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   problem,        intent=in,     poisson_problem data type              !
    !   root,           intent=in,     'initial', 'current', etc.             !
    !   method,         intent=in,     method for solving poisson problem     !
    !   model,          intent=in,     the dielectric model being used        !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use simulation_cell, only: CELL_INFO
    use visual, only: visual_scalarfield

    implicit none

    ! jd: Arguments
    type(POISSON_PROBLEM), intent(in) :: problem
    character(len=*), intent(in) :: method
    character(len=*), intent(in) :: root
    character(len=*), intent(in) :: model
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell

    ! jd: Local variables
    character(len=200) :: suffix
    integer, save :: counter = 0

    !------------------------------------------------------------------------

    call utils_trace_in('write_poisson_problem')

    counter = counter+1
    write(suffix,'(i20)') counter
    suffix = adjustl(suffix)

    if(method=='FFT') then
       ! Write out the charge density
       call visual_scalarfield(problem%rho_elec, &
            grid, cell, 'Charge density','_'//root//'_rho_'//suffix, &
            conversion_factor = ANGSTROM**3)
       ! Write out the periodic potential
       call visual_scalarfield(problem%phi_pbc, grid, cell, &
            'Periodic potential found by FFTs','_'//root//'_phi_pbc_'//suffix)
    end if

    if(method=='INITIAL') then
       ! Only want to write out the charge density used to initialize the
       ! dielectric, this could either be the initial electronic density
       ! or else a sum of Gaussians on ions
       call visual_scalarfield(problem%rho_elec, &
            grid, cell, 'Initial charge density',root//'_rho_'//suffix, &
            conversion_factor = ANGSTROM**3)
    end if

    if(method=='DIRECT' .or. method=='CORRECTIVE') then
       ! Write out the solvation potential
       call visual_scalarfield(problem%phi_mol, grid, cell, &
            'Solvation potential',root//'_phi_mol_'//suffix)

       ! Write out the dielectric on the full grid
       call visual_scalarfield(problem%eps_full, grid, cell, &
            'Dielectric on full grid',root//'_eps_full_'//suffix)

       ! Write out the charge density - scaled with 1/V to get true density
       call visual_scalarfield(problem%rho_elec, grid, cell, &
            'Electronic charge density',root//'_rho_'//suffix, &
            conversion_factor = ANGSTROM**3)

       ! Write out the ion charge density - scaled with 1/V to get true density
       call visual_scalarfield(problem%rho_ion, grid, cell, &
            'Smeared ion charge density', root//'_rho_ion_'//suffix, &
            conversion_factor = ANGSTROM**3)

       ! Write out the molecular charge density - scaled with 1/V as above
       call visual_scalarfield(problem%rho_elec+problem%rho_ion, &
            grid, cell, 'Molecular charge density = elec + ion', &
            root//'_rho_mol_'//suffix, conversion_factor = ANGSTROM**3)

       if(model=='SELF_CONSISTENT') then
          ! Write out the dielectric potential
          call visual_scalarfield(problem%V_eps, grid, cell, &
               'Dielectric potential term V_eps',root//'_V_eps_'//suffix)
          call visual_scalarfield(problem%V_apolar, grid, cell, &
               'Dielectric potential term V_apolar',root//'_V_apolar_'//suffix)
       end if

       if(method=='CORRECTIVE') then
          ! Write out the periodic potential
          call visual_scalarfield(problem%phi_pbc, grid, cell, &
               'Periodic potential found by FFTs',root//'_phi_pbc_'//suffix, &
               conversion_factor = ANGSTROM**3)

          ! Write out the corrective potential
          call visual_scalarfield(problem%phi_corr, grid, cell, &
               'Corrective potential',root//'_phi_corr_'//suffix)
       end if
    end if

    call utils_trace_out('write_poisson_problem')

  end subroutine write_poisson_problem
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine zero_solvation_terms(solvation_terms, zero_apolar, zero_vacuum)
    !=========================================================================!
    ! Zeroes the datastructure 'solvation_terms'. This helps hide the imple-  !
    ! mentation details of the structure from calling modules.                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   solvation_terms (inout): The datastructure to zero.                   !
    !   zero_apolar (in): If .true., the apolar terms are zeroed too.         !
    !   zero_vacuum (in): If .true., the vacuum elec fixed term is zeroed too.!
    !-------------------------------------------------------------------------!
    ! NB: Sometimes we only want to zero certain elements (apolar, vacuum),   !
    !     hence solvation_terms must be intent(inout).                        !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2019.                             !
    !=========================================================================!

    implicit none

    ! jd: Arguments
    type(SOLVATION_ENERGY_TERMS), intent(inout) :: solvation_terms
    logical, intent(in)                         :: zero_apolar
    logical, intent(in)                         :: zero_vacuum

    ! -------------------------------------------------------------------------

    solvation_terms%E_elec_fixed = 0.0_DP
    solvation_terms%E_elec_mob = 0.0_DP
    solvation_terms%E_osmo = 0.0_DP
    solvation_terms%E_acc = 0.0_DP
    solvation_terms%E_atmo = 0.0_DP
    solvation_terms%E_chempot = 0.0_DP
    solvation_terms%E_pure = 0.0_DP
    solvation_terms%E_elec_minus_mob = 0.0_DP
    solvation_terms%E_minus_ktv_gamma_c_i_bulk = 0.0_DP
    solvation_terms%E_minus_kt_int_c_i = 0.0_DP

    if(zero_apolar) then
       solvation_terms%E_apolar_cavitation = 0.0_DP
       solvation_terms%E_apolar_disrep = 0.0_DP
       solvation_terms%E_apolar_sasa = 0.0_DP
       solvation_terms%E_apolar_sav = 0.0_DP
       solvation_terms%cavity_surface_area = 0.0_DP
       solvation_terms%cavity_volume = 0.0_DP
    end if

    if(zero_vacuum) then
       solvation_terms%E_elec_fixed_vac = 0.0_DP
    end if

  end subroutine zero_solvation_terms
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


end module is_poisson
