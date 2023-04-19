!================================================================!
!                    ____  ____  _____ ____                      !
!                   |  _ \|  _ \| ____/ ___|                     !
!                   | | | | | | |  _|| |                         !
!                   | |_| | |_| | |___ |___                      !
!                   |____/|____/|_____\____|                     !
!                                                                !
!                 DDEC charge analysis module                    !
!                      Version 2.0 BETA                          !
!                         10 Dec 2013                            !
!                                                                !
!                    Written by Louis Lee                        !
!                                                                !
!----------------------------------------------------------------!
!                                                                !
! This module performs the DDEC Hirshfeld AIM charge analysis    !
! Called from properties_calculate in properties_mod             !
!                                                                !
! Notes:                                                         !
! 1. DDEC cutoff radius r_cut must be smaller than the size of   !
!    the simulation cell, as in the NGWF localization radius.    !
! 2. Compensation radii for compensation densities should        !
!    ideally be specified manually, as the automated routine is  !
!    really dumb and potentially dodgy.                          !
! 3. Currently when set to read refdens, only option #3 i.e.     !
!    copy values directly from file is supported. If rcut or     !
!    npts do not match program will halt. To change this         !
!    behavior modify the line at the comment 'MOD#02'            !
! 4. 'eff_nuc_charge' in ref_dens must be integral (MOD#03)      !
! 5. Aborts entire calculation if atomic solver was used to      !
!    generate reference densities and 'reshape_dens' or          !
!    'use_coredens' are enabled since the current atomic solver  !
!    is tested for pspots only and are incompatible with both    !
!    options (MOD#04)                                            !
! 6. MOD#05 - Created new array to store expected core charges   !
!    which would be the same as the corrected charges. This is   !
!    to ensure correct atomic charge can be computed from simply !
!    'Z - atom_popn'. Also, 'core_diff' might not actually be    !
!    needed and could be removed later on                        !
! 7. MOD#07 - When reading densities, allow failure to read in   !
!    reference densities. Such refdens will be marked as         !
!    uninitialized without crashing the program. Attempt to      !
!    use these uninitialized refdens will still cause an abort   !
! 8. MOD#08 - moment order calculations only work on the valence !
!    densities (same as multipole)                               !
! 9. MOD#09 - Cleaned up stuff that allows user-specified atomic !
!    configuration to be used in solver                          !
! 10.MOD#12 - Corrected printing of 'rad_dens_atom' to print ISA !
!    reference density from w_A(r) instead of w_A(r) itself      !
! 11.MOD#13 - Allow user-specified electronic configuration to   !
!    be used when determining rcomp automatically                !
!----------------------------------------------------------------!
! 12. COMMENT#00 - local_box <= grid_size                        !
! 13. COMMENT#01 - Removed 'SOLVE' option                        !
! 14. COMMENT#02 - Various message printouts                     !
! 15. COMMENT#03 - Parsing                                       !
! 16. COMMENT#04 - Solver                                        !
!================================================================!

! 1. read_dens integration tolerance?
! 2. Why is ih_ratio 1.0 if below threshold when we are taking an average at
!    the end?
! 3. What is the point of core density partitoning if in the end the total
!    charge is not used (only profile matters?)
! 4. Bug in core correction ( exit loop did not use abs(nA_correction(j) )?
! 5. Difference between IH and Yavg results
! 6. Why initialize ref_core_density as neutral density on the 1st core
!    iterator?
! 7. Why 1e-16 for all densities?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        1         2         3         4         5         6         7         8
!2345678901234567890123456789012345678901234567890123456789012345678901234567890
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ddec

  use constants, only: DP, PI, stdout
  use geometry, only: POINT
  use ion, only: ELEMENT
  use rundat, only: pub_debug, pub_debug_on_root

  implicit none

  private

  real(kind=DP), parameter :: FOUR_PI = 4.0_DP*PI

!------------------------------------------------------------------------------!
! lpl: Object that stores the numerous radial densities for each atom          !
!------------------------------------------------------------------------------!

  type DDEC_RAD_DENS

    ! lpl: Number of radial grid points & positions
    integer :: npts       ! Number of radial grid points
    real(kind=DP) :: rcut ! Cutoff radius
    real(kind=DP) :: dr   ! rcut/npts

    ! lpl: Radial grid coordinates for convenience
    real(kind=DP), allocatable :: rad(:)
    integer, allocatable :: rad_npts(:)

    ! lpl: Spherically-averaged density and weight factors
    real(kind=DP), allocatable :: total_part_dens(:) ! Total weighting factor
    real(kind=DP), allocatable :: ref_dens(:)        ! IH reference density
    ! lpl: Core density
    real(kind=DP), allocatable :: core_dens(:)

    ! lpl: Buffer to store new weighting factors
    real(kind=DP), allocatable :: new_part_dens(:)
    real(kind=DP), allocatable :: old_part_dens(:)
    real(kind=DP), allocatable :: Yavg_dens(:)
    real(kind=DP), allocatable :: Yavg_sqrt(:)
    real(kind=DP), allocatable :: Yavg_ratio(:)
    real(kind=DP), allocatable :: ih_ratio(:)

    ! lpl: Copy of atomic number for convenience
    integer :: iat

    ! lpl: Utility pointer for misc stuff. May be used in various
    !      'templatized' subroutines that operate on the DDEC_RAD_DENS
    !      object. Must be made to point to an appropriate radial density
    !      before any usage.
    real(kind=DP), pointer :: target_dens(:)

    ! lpl: Misc variables for density reshaping etc.
    real(kind=DP) :: renorm_slope ! 'u_A' in Eq. 73 of JCTC 2012
    real(kind=DP) :: sum_part_dens, sum_part_dens_sqrt
    real(kind=DP), allocatable :: exp_factor(:) ! Exponential decay factor

  end type DDEC_RAD_DENS

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! lpl: Dummy objects to act as containters, to allow allocation of arrays      !
!      and array-of-arrays inside subroutines.                                 !
!------------------------------------------------------------------------------!

  type DENSITY_PROFILE
     real(kind=DP), allocatable :: dens(:)
  end type DENSITY_PROFILE

  type RADIAL_GRID
    integer :: npts
    type(POINT), allocatable :: coord(:)
  end type RADIAL_GRID

  ! lpl: Object to store record of allocated reference densities
  type REF_ALLOC_SEQ
    integer :: num
    integer, allocatable :: alloc(:,:)
  end type REF_ALLOC_SEQ

  ! lpl: Object to store key list for 'ddec_read_atom'
  type KEY_PARAM_PAIR
     character(len=512) :: key, param
     character(len=16) :: delimiter
     logical :: initialized, essential
  end type

  type KEY_PARAM_PAIR_PTR
     type(KEY_PARAM_PAIR), pointer :: ptr
  end type

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! lpl: Object that stores IH reference densities for each species              !
!      This should be initialized once and remain unmodified thereafter        !
!------------------------------------------------------------------------------!

  type DDEC_REF_DENS

    type(ELEMENT) :: element          ! Copy of element object for convenience
    integer :: max_anion, max_cation  ! Hard limit on number of ionic states
    ! lpl: Effective nuclear charge i.e. (nucleus charge) - (number of
    !      pseudized core charges). Must be integral as the numbers stored
    !      here are used for indexing arrays. This is checked elsewhere in
    !      this module. Moreover, the c2 core densities are generated
    !      only for complete shells and integral number of electrons.
    integer :: eff_nuc_charge
    ! lpl: This depends on whether coredens is included. If not
    !      then this should be equal to 'ref_dens%element%ion_charge',
    !      otherwise it should be equal to 'ref_dens%element%atomic_number'
    real(kind=DP) :: ddec_nuc_charge

    integer :: npts       ! Number of radial grid points
    real(kind=DP) :: rcut ! Cutoff radius
    real(kind=DP) :: dr   ! rcut/npts

    ! lpl: Radial grid coordinates for convenience
    real(kind=DP), allocatable :: rad(:)

    ! lpl: Reference density profiles and c3 flags
    type(DENSITY_PROFILE), allocatable :: ion_dens(:)
    logical, allocatable :: init_status(:)
    logical, allocatable :: c3flag(:)
    ! lpl: Pointer to one of the cationic density which is the core density
    real(kind=DP), pointer :: core_dens(:)

    ! lpl: Compensation radius data
    real(kind=DP), allocatable :: rcomp(:,:)
    real(kind=DP), allocatable :: neutral_orb_info(:,:)
    integer :: neutral_norbs
    character(len=256), allocatable :: refconf(:)

    ! lpl: List of coorindates on the surface of each spherical shell for
    !      density interpolation
    type(RADIAL_GRID), allocatable :: rad_grid(:)
    ! lpl: Average radius for all points lying within each shell
    real(kind=DP), allocatable :: avg_rad(:)
    ! lpl: Index of the largest spherical shell requiring interpolation
    integer :: max_interp_shell

    ! lpl: Radial shell volume
    real(kind=DP), allocatable :: rad_shell_vol(:)

    ! lpl: Initialization flag
    logical :: initialized

  end type DDEC_REF_DENS

!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! lpl: Miscellaneous Cartesian grid related data objects                       !
!------------------------------------------------------------------------------!

  ! lpl: Integer version of POINT
  type INDEX_POINT
    integer :: ni1, ni2, ni3
  end type INDEX_POINT

  ! lpl: Severely dumbed-down version of GRID_INFO
  !      Local box is assumed to have the same GRID_INFO as the supercell
  !      throughout the code, so this info is not stored here
  type DDEC_LOCAL_BOX

    ! lpl: local_box dimensions
    integer :: n1, n2, n3
    ! lpl: Supercell Cartesian coordinates of the origin of the local_box
    !      (index (1,1,1)) for all atoms
    type(POINT), allocatable :: box_coord_origin(:)
    ! lpl: Supercell array index where the origin (index (1,1,1)) of the
    !      local_box is located for all atoms
    type(INDEX_POINT), allocatable :: box_index_origin(:)

    ! lpl: Grid data
    real(kind=DP), allocatable :: total_val(:,:,:)       ! Density
    real(kind=DP), allocatable :: total_pseudoval(:,:,:) ! Pseudodensity

    real(kind=DP), allocatable :: ref_pseudoval(:,:,:)   ! Reference (IH)
    real(kind=DP), allocatable :: Yavg_pseudoval(:,:,:)  ! Yavg Pseudodensity

    real(kind=DP), allocatable :: core_val(:,:,:)        ! Core
    real(kind=DP), allocatable :: core_pseudoval(:,:,:)  ! Core Pseudodensity

    real(kind=DP), allocatable :: aim_val(:,:,:)

    ! lpl: Array for temp data
    real(kind=DP), allocatable :: temp_val(:,:,:)

  end type DDEC_LOCAL_BOX

!------------------------------------------------------------------------------!

  public :: ddec_main
! lpl: c.f. COMMENT#06
!  public :: ddec_rmse

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------!
! lpl: Main DDEC subroutine                                                    !
!------------------------------------------------------------------------------!

  subroutine ddec_main(ddec_charge, density_fine, pseudodensity_fine, &
       ref_pseudodensity_fine, Yavg_pseudodensity_fine, core_density_fine, &
       core_pseudodensity_fine, grid, mdl, total_integrated_dens, &
       include_refdens, include_coredens, coredens_provided)

  !======================================================================!
  ! This subroutine performs the DDEC charge partitioning                !
  ! Called from properties_mod                                           !
  !----------------------------------------------------------------------!
  ! Written by Louis Lee 29-Aug-2013                                     !
  !======================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_extract_box, &
         cell_grid_deposit_box
    use comms, only: comms_barrier, comms_bcast, comms_reduce, comms_alltoall, &
         pub_my_proc_id, pub_root_proc_id, pub_on_root, pub_total_num_procs
    use constants, only: DP, PI, stdout
    use geometry, only: POINT, unit_vector, geometry_magnitude, &
         operator(.DOT.), operator(.CROSS.), operator(*), operator(+), &
         operator(-), local_displacement
    use model_type, only: MODEL
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_print_qc, pub_ddec_write, pub_ddec_rad_rcut, &
         pub_ddec_rad_npts, pub_ddec_maxit, pub_ddec_core_maxit, &
         pub_ddec_conv_thresh, pub_ddec_IH_frac, pub_ddec_hirshfeld, &
         pub_ddec_zero_thresh, pub_ddec_refdens_init, pub_ddec_multipole, &
         pub_ddec_interp_rad_dens, pub_ddec_moment, & !pub_ddec_moment_order, &
         pub_ddec_ionic_range, pub_ddec_rad_shell_mode, &
         pub_ddec_reshape_dens, pub_ddec_core_correction, &
         pub_ddec_c3_refdens, pub_ddec_core_corr_maxit, &
         pub_ddec_eff_decay_exp, pub_ddec_eff_decay_rmin, pub_ddec_aniso
    use services, only: services_linear_interpolation
    use utils, only: utils_abort, utils_alloc_check, utils_assert, &
         utils_close_unit_check, utils_dealloc_check, utils_qc_print, &
         utils_open_unit_check, utils_unit
    use timer, only: timer_clock
    use function_basis, only: FUNC_BASIS

    implicit none

    ! lpl: Inputs
    real(kind=DP), intent(out) :: ddec_charge(par%nat) ! Array of atom dQ
    type(GRID_INFO), intent(in) :: grid
    type(MODEL), intent(in) :: mdl

    real(kind=DP), target, intent(inout) :: density_fine(:,:,:)
    real(kind=DP), target, intent(inout) :: pseudodensity_fine(:,:,:)
    real(kind=DP), target, intent(inout) :: ref_pseudodensity_fine(:,:,:)
    real(kind=DP), target, intent(inout) :: Yavg_pseudodensity_fine(:,:,:)
    real(kind=DP), target, intent(inout) :: core_density_fine(:,:,:)
    real(kind=DP), target, intent(inout) :: core_pseudodensity_fine(:,:,:)

    real(kind=DP), intent(in) :: total_integrated_dens
    logical, intent(in) :: include_refdens, include_coredens, coredens_provided

    ! lpl: Local variables
    ! lpl: Pointers to input densities (i.e. input data). These pointers may
    !      be nullified in certain DDEC routines to prevent access and the
    !      potential of overwriting data that are meant to be left alone.
    !      Someone might want to check that this is unnecessary, and if so,
    !      remove them.
    real(kind=DP), pointer,  dimension(:,:,:) :: total_dens_fine, &
         total_pseudodens_fine
    real(kind=DP), pointer, dimension(:,:,:) :: ref_pseudodens_fine, &
         Yavg_pseudodens_fine
    real(kind=DP), pointer, dimension(:,:,:) :: core_dens_fine, &
         core_pseudodens_fine

    ! lpl: Indexes etc.
    integer :: iat, orig_iat, lc_iat, it, ddec_it, isp, iion

    integer :: npts_max
    real(kind=DP) :: rcut_max

    real(kind=DP) :: dq, q_total, q_renorm, q_spill, q_core_total, ddens_max
    integer :: dq_cons, ddens_max_cons

    ! lpl: Number of consecutive iterations below threshold for convergence
    integer, parameter :: dq_cons_num = 4

    ! lpl: More indexes, bounds, etc.
    type(POINT) :: grid_a1, grid_a2, grid_a3
    type(POINT) :: grid_ua1, grid_ua2, grid_ua3
    real(kind=DP) :: grid_la1, grid_la2, grid_la3
    real(kind=DP) :: grid_lda1, grid_lda2, grid_lda3
    real(kind=DP) :: vox_vol, vox_maxlen, total_grid_dens, total_core_dens

    ! lpl: Radial density objects for each atom or species
    type(DDEC_RAD_DENS), pointer :: rad_dens(:)
    type(DDEC_REF_DENS), pointer :: ref_dens(:)

    ! lpl: IH reference densities allocation record
    !      Should be removed if this is pointless
    type(REF_ALLOC_SEQ) :: ddec_ref_alloc_list

    ! lpl: Local density box
    type(DDEC_LOCAL_BOX), target :: local_box_object ! Allows pointer to obj
    type(DDEC_LOCAL_BOX), pointer :: local_box_fine
    logical :: i_need_box, i_have_box
    real(kind=DP), allocatable :: deposit_extract_buffer(:,:,:)

    ! lpl: Atomic population arrays
    real(kind=DP), allocatable :: atom_popn(:), atom_popn_prev(:)
    real(kind=DP), allocatable :: core_popn(:), core_diff(:)
    ! lpl: MOD#05
    real(kind=DP), allocatable :: val_popn(:), ddec_core_popn(:)
    ! lpl: ISA fraction
    real(kind=DP) :: ddec_ISA_frac

    ! lpl: File name prefix for DDEC-related outputs e.g. radial densities
    character(len=128) :: output_filename

    ! lpl: DDEC effective decay exponents
    real(kind=DP), allocatable :: ddec_eff_decay_exp(:,:)
    ! lpl: DDEC multipole calculation
    real(kind=DP), allocatable :: ddec_dipoles(:,:), ddec_quadrupoles(:,:)
    ! lpl: AIM moment
    real(kind=DP), allocatable :: ddec_moments(:,:)

    ! lpl: Constants - refer to JCTC 2012 for meanings
    real(kind=DP) :: sqrt_pub_ddec_zero_thresh, sqr_pub_ddec_zero_thresh
    real(kind=DP) :: overlap_thresh
    real(kind=DP), parameter :: valence_exp_factor = 1.75_DP ! 1/bohr
    real(kind=DP), parameter :: core_exp_factor = 2.0_DP     ! 1/bohr
    real(kind=DP), parameter :: max_renorm = 10.0_DP

    integer :: debug_unit
    integer :: ierr

    ! aeaa: variables used for finding off center point charges
    real(kind=DP), allocatable :: optimum_posit(:,:,:), optimum_charges(:,:)
    real(kind=DP), allocatable :: error_aniso(:), error_aniso_before(:)
    real(kind=DP), allocatable :: vdw_radius(:), aim_volume(:), &
         overall_dipole_noextra(:), overall_dipole_withpc(:)
    real(kind=DP), allocatable :: dist_atoms(:,:)

    integer, allocatable :: no_charges(:)
    integer :: i, j, ii, free_io, output_unit

    type(POINT), allocatable :: atom_ctr_point(:), atom_ctr_point_array(:)

    character(len=2), allocatable :: all_element_types(:)

    ! lpl: Start timer
    call timer_clock('ddec_main',1)

    if(pub_debug) then
       if(pub_on_root) write(stdout,'(a)') ' DEBUG: Entering ddec_main'

       debug_unit = utils_unit()
       write(output_filename,'(a,I6.6)') 'DDEC_DEBUG_', pub_my_proc_id
       open(unit=debug_unit,form="formatted", &
            file=trim(output_filename),action="write",iostat=ierr)
       call utils_open_unit_check('ddec_main','DDEC_DEBUG',ierr)
    end if

    ierr=0

  !----------------------------------------------------------------------------!
  ! lpl: Check and initialize variables                                        !
  !----------------------------------------------------------------------------!

    call utils_assert( pub_ddec_IH_frac >= 0.0_DP .and. &
         pub_ddec_IH_frac <= 1.0_DP, ' ERROR: ddec_IH_fraction must be &
         &between [0,1]')
    ddec_ISA_frac = 1.0_DP - pub_ddec_IH_frac

    call utils_assert(pub_ddec_maxit > 0,' ERROR: ddec_maxit < 0')
    call utils_assert(pub_ddec_core_maxit > 0,' ERROR: ddec_core_maxit < 0')
    call utils_assert(pub_ddec_conv_thresh > 0,' ERROR: ddec_conv_threshold&
         & must be > 0')

    sqrt_pub_ddec_zero_thresh = sqrt(pub_ddec_zero_thresh)
    sqr_pub_ddec_zero_thresh  = pub_ddec_zero_thresh**2

    ! lpl: Check and set Cartesian density pointers
    call ddec_check_3D_array(density_fine, &
         grid%ld1, grid%ld2, grid%max_slabs12)
    call ddec_check_3D_array(pseudodensity_fine, &
         grid%ld1, grid%ld2, grid%max_slabs12)
    total_dens_fine => density_fine             ! Input total/valence density
    total_pseudodens_fine => pseudodensity_fine ! Total pseudodensity

    ! lpl: If IH reference densities are required, then initialize
    !      these too. Otherwise, set the corresponding pointers to NULL
    if ( include_refdens ) then
       call ddec_check_3D_array(ref_pseudodensity_fine, &
            grid%ld1, grid%ld2, grid%max_slabs12)
       call ddec_check_3D_array(Yavg_pseudodensity_fine, &
            grid%ld1, grid%ld2, grid%max_slabs12)
       ref_pseudodens_fine => ref_pseudodensity_fine
       Yavg_pseudodens_fine => Yavg_pseudodensity_fine
    else
       nullify(ref_pseudodens_fine)
       nullify(Yavg_pseudodens_fine)
    end if

    ! lpl: Same as above, but for core densities.
   if ( include_coredens ) then
       call ddec_check_3D_array(core_density_fine, &
            grid%ld1, grid%ld2, grid%max_slabs12)
       call ddec_check_3D_array(core_pseudodensity_fine, &
            grid%ld1, grid%ld2, grid%max_slabs12)
       core_dens_fine => core_density_fine
       core_pseudodens_fine => core_pseudodensity_fine
    else
       nullify(core_dens_fine)
       nullify(core_pseudodens_fine)
    end if

    local_box_fine => local_box_object

    overlap_thresh = pub_ddec_conv_thresh/(10.0_DP*max_renorm)

  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lpl: Determine geometry of local density box                               !
  !----------------------------------------------------------------------------!
    ! lpl: Max length of voxel
    grid_a1%x = abs(grid%da1%x) + abs(grid%da2%x) + abs(grid%da3%x)
    grid_a1%y = abs(grid%da1%y) + abs(grid%da2%y) + abs(grid%da3%y)
    grid_a1%z = abs(grid%da1%z) + abs(grid%da2%z) + abs(grid%da3%z)
    vox_maxlen = geometry_magnitude(grid_a1)

    ! lpl: Lattice vectors
    grid_a1 = grid%n1*grid%da1
    grid_a2 = grid%n2*grid%da2
    grid_a3 = grid%n3*grid%da3

    ! lpl: Lattice vector magnitudes
    grid_la1 = geometry_magnitude(grid_a1)
    grid_la2 = geometry_magnitude(grid_a2)
    grid_la3 = geometry_magnitude(grid_a3)

    ! lpl: Grid unit vectors
    grid_ua1 = unit_vector(grid_a1)
    grid_ua2 = unit_vector(grid_a2)
    grid_ua3 = unit_vector(grid_a3)

    ! lpl: Voxel lattice vector
    grid_lda1 = geometry_magnitude(grid%da1)
    grid_lda2 = geometry_magnitude(grid%da2)
    grid_lda3 = geometry_magnitude(grid%da3)

    ! lpl: Volume of each voxel
    vox_vol = grid%da1.DOT.(grid%da2.CROSS.grid%da3)
    total_grid_dens = total_integrated_dens*vox_vol

  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lpl: Initialize radial density arrays for each atom or species             !
  !----------------------------------------------------------------------------!

    if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: ddec_main: Entering initialization'

    ! lpl: Print DDEC header
    if(pub_on_root) then
       write(stdout,'(a)') '================= DDEC Charges ================='
       write(stdout,'(a)') '              Version 01-Sep-2014'
       write(stdout,'(a)') '------------------------------------------------'
       write(stdout,'(a,1x,f16.8,a)') &
            ' Total integrated density = ',total_integrated_dens*vox_vol,' e'
    end if

    ! lpl: Atomic charge arrays
    allocate(atom_popn(par%nat),stat=ierr)
    call utils_alloc_check('ddec_main','atom_popn',ierr)
    allocate(core_popn(par%nat),stat=ierr)
    call utils_alloc_check('ddec_main','core_popn',ierr)
    allocate(core_diff(par%nat),stat=ierr)
    call utils_alloc_check('ddec_main','core_diff',ierr)
    allocate(val_popn(par%nat),stat=ierr)
    call utils_alloc_check('ddec_main','val_popn',ierr)
    ! lpl: MOD#05
    allocate(ddec_core_popn(par%nat),stat=ierr)
    call utils_alloc_check('ddec_main','ddec_core_popn',ierr)

    allocate(atom_popn_prev(par%nat),stat=ierr)
    call utils_alloc_check('ddec_main','atom_popn_prev',ierr)

    !aeaa - off_center part
    allocate( optimum_posit(3,5,par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','optimum_posit', ierr)
    allocate( optimum_charges(5,par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','optimum_charges', ierr)
    allocate( error_aniso(par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','error_aniso', ierr)
    allocate( error_aniso_before(par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','error_aniso_before', ierr)
    allocate( vdw_radius(par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','vdw_radius', ierr)
    allocate( aim_volume(par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','aim_volume', ierr)
    allocate( no_charges(par%nat),  stat=ierr )
    call utils_alloc_check('ddec_main','no_charges', ierr)
    allocate( atom_ctr_point(par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','atom_ctr_point', ierr)
    allocate( atom_ctr_point_array(par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','atom_ctr_point_array', ierr)
    allocate( overall_dipole_withpc(3), stat=ierr )
    call utils_alloc_check('ddec_main',' overall_dipole_withpc', ierr)
    allocate( overall_dipole_noextra(3), stat=ierr )
    call utils_alloc_check('ddec_main','overall_dipole_noextra', ierr)
    allocate( dist_atoms(par%nat, par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','dist_atoms', ierr)
    allocate( all_element_types(par%nat), stat=ierr )
    call utils_alloc_check('ddec_main','all_element_types', ierr)


    ! aeaa: Variables declared, sort out
    optimum_posit =  0.0_DP
    optimum_charges = 0.0_DP
    error_aniso = 0.0_DP
    error_aniso_before = 0.0_DP
    aim_volume = 0.0_DP
    vdw_radius = 0.0_DP
    no_charges = 0
    overall_dipole_noextra = 0.0_DP
    overall_dipole_withpc = 0.0_DP
    dist_atoms = 0.0_DP

    atom_popn = 0.0_DP
    core_popn = 0.0_DP
    core_diff = 0.0_DP
    val_popn  = 0.0_DP
    ddec_core_popn = 0.0_DP
    atom_popn_prev = 0.0_DP


    nullify(rad_dens)
    nullify(ref_dens)

    ! lpl: Initialize IH reference densities (replicated for all species on
    !      all procs)
!djc    if (pub_on_root) then
!       write(stdout,'(a)') '------------------------------------------------'
!       write(stdout,'(a)') '             IH Reference Densities             '
!       write(stdout,'(a)') '------------------------------------------------'
!djc    end if
    call internal_ref_dens_init(ref_dens)
    call comms_barrier

    ! lpl: QC test
    if ( pub_on_root ) then
       if ( pub_print_qc ) then
          do isp=1,min(3,par%num_species)
             ! lpl: Number of intervals (reuse variable 'dq_cons')
             dq_cons = max(1,ref_dens(isp)%npts/5)
             do it=(1+min(1,ref_dens(isp)%npts/4)),min(4,ref_dens(isp)%npts)
                output_filename = ''
                write(output_filename,'(3(a,I0))') 'sp_',isp,'_rad_',it, &
                     '_qn_',pub_ddec_ionic_range
                call utils_qc_print( trim(adjustl(output_filename))//&
                     '_rad',ref_dens(isp)%dr*(dq_cons*it - 0.5_DP) )
                call utils_qc_print( trim(adjustl(output_filename))//&
                     '_dens',ref_dens(isp)%&
                     ion_dens(-pub_ddec_ionic_range)%dens(it*dq_cons) )
             end do
          end do
       end if
       output_filename = ''
       dq_cons = 0
    end if
    call comms_barrier

    if(pub_on_root) &
       write(stdout,'(a)') '------------------------------------------------'

    ! lpl: Initialize weighting factors (or spherically-averaged 'partial
    !      densities'. Replicated for all atoms on all procs)
    call internal_rad_dens_init(rad_dens,ref_dens)

    ! lpl: Max number of points and DDEC cutoff radius for initializing buffer
    npts_max = maxval(rad_dens(:)%npts)
    rcut_max = maxval(rad_dens(:)%rcut)

    ! lpl: Print info on current DDEC parameters
    if(pub_on_root) then
       write(stdout,'(a,f14.7,/a,I5,/a,f14.7,a)') &
            ' Parameters: IH fraction = ', pub_ddec_IH_frac, &
            '             nshells_max = ', npts_max, &
            '             rcut_max    = ', rcut_max,' bohr'
       write(stdout,'(a,L1,3(/a,L1))') &
            '             c3_refdens       = ', pub_ddec_c3_refdens, &
            '             reshape_dens     = ', pub_ddec_reshape_dens, &
            '             include_coredens = ', include_coredens, &
            '             core_correction  = ', pub_ddec_core_correction
    end if

    ! lpl: Initialize local density box
    call internal_local_box_init(local_box_fine,grid,rcut_max, &
         include_refdens,include_coredens)

    ! lpl: Allocate buffer for cell_grid_extract_box after figuring out
    !      '12-dimensions' of local density box
    !      This buffer is also used as a 'local_box' buffer with
    !      expected dimensions of the local box
    allocate(deposit_extract_buffer(local_box_fine%n1,local_box_fine%n2, &
         grid%max_slabs12),stat=ierr)
    call utils_alloc_check('ddec_main','deposit_extract_buffer',ierr)

    ! lpl: Optionally call this subroutine, which requires local_box, to
    !      obtain the average radius of all Cartesian points lying within
    !      each shell, analogous to 'ref_dens%avg_rad'
    if ( trim(adjustl(pub_ddec_rad_shell_mode)) == 'AVERAGE' ) &
         call internal_rad_dens_avg_rad(rad_dens,ref_dens,local_box_fine)

    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: ddec_main: Leaving initialization'

  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lpl: Core DDEC iterative loop                                              !
  !----------------------------------------------------------------------------!

    if ( include_coredens ) then

    if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: ddec_main: Entering core iterator'

       ! lpl: Nullify density pointers that should not be used here
       nullify(total_dens_fine)
       nullify(total_pseudodens_fine)
       nullify(ref_pseudodens_fine)
       nullify(Yavg_pseudodens_fine)

       ! lpl: Accumulate the expected core population. Due to the finesse of
       !      the density grid, the integrated sum from the grid will not be
       !      accurate. This does not matter as the valence population
       !      will be used to describe the DDEC charges. Furthermore, the
       !      core density can be 'corrected' later on such that it integrates
       !      to the correct charge for the current grid spacing. See PAPER
       !      for details.
       total_core_dens = 0.0_DP
       do iat=1,par%nat
          total_core_dens = total_core_dens + &
               mdl%elements(par%orig_atom(iat))%atomic_number - &
               ref_dens(mdl%elements(par%orig_atom(iat))%species_number)&
               %eff_nuc_charge
       end do

       ! lpl: Initialize core density. Deposits all radial core reference
       !      densities into 'core_dens_fine' if a suitable input Cartesian
       !      core density is not provided.
       if ( .not. coredens_provided ) then
          call internal_deposit_box_all(core_dens_fine, rad_dens, 'core_dens', &
               local_box_fine, 'aim_val', deposit_extract_buffer)
       end if

       if (pub_debug_on_root)  then
          write(debug_unit,'(a,1x,f14.7)') ' SUM CORE DENS =', total_core_dens

          do orig_iat=1,par%nat
             iat = par%distr_atom(orig_iat)
             write(debug_unit,'(a,1x,e20.10)') 'MAX CORE DENSITY VALUE: '//&
                  mdl%elements(orig_iat)%species_id, &
                  maxval(rad_dens(iat)%core_dens)

             ! Print out core radial density profiles
             write(output_filename,'(a,I5.5)') 'core_atom_',orig_iat
             call ddec_write_dens( rad_dens(iat)%core_dens, &
                  rad_dens(iat)%rad, rad_dens(iat)%npts, &
                  trim(adjustl(output_filename)) )
          end do
       end if

       ! lpl: In CHARGEMOL, the radial core density of each atom is initialized
       !      as the NEUTRAL density, not core. As 'rad_dens%core_dens' currently
       !      store the core density profiles, we need to replace them.
       do iat=1,par%nat
          rad_dens(iat)%core_dens = ref_dens(mdl%elements(par%orig_atom(iat))&
               %species_number)%ion_dens(0)%dens
       end do

       ddens_max = 0.0_DP
       core_popn = 0.0_DP
       dq = pub_ddec_conv_thresh
       dq_cons = 0 ! Consecutive number of iterations with dQ < threshold
       ddens_max_cons = 0

       if (pub_debug_on_root) then
       ! Prints out the density profile of 'rad_dens%core_dens' at the start
       ! of the core iterator
       do orig_iat=1,par%nat
          iat = par%distr_atom(orig_iat)
          write(output_filename,'(a,I5.5)') 'init_core_atom_',orig_iat
          call ddec_write_dens( rad_dens(iat)%core_dens, &
               rad_dens(iat)%rad, rad_dens(iat)%npts, &
               trim(adjustl(output_filename)) )
       end do
       end if

       ! lpl: Iterative loop to partition the core density according to
       !      the ISA scheme
       do ddec_it=1,pub_ddec_core_maxit

          ! lpl: Reset DDEC charges
          atom_popn_prev = core_popn
          core_popn = 0.0_DP

          ! lpl: Reset and update promolecular densities using new
          !      radial densities
          core_pseudodens_fine = 0.0_DP
          call internal_deposit_box_all(core_pseudodens_fine, rad_dens, &
               'core_dens', local_box_fine, 'aim_val', deposit_extract_buffer)

          ! lpl: Calculate new weighting factors for each atom on this proc
          do lc_iat=1,par%max_atoms_on_proc

             ! lpl: Determine whether we need to extract density box for this atom
             if (lc_iat <= par%num_atoms_on_proc(pub_my_proc_id)) then
                i_need_box = .true.
                iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
             else
                i_need_box = .false.
                iat = par%first_atom_on_proc(pub_my_proc_id) ! Prevents segfault
             end if

             ! lpl: Extract Cartesian core density and promolecular density
             local_box_fine%core_val = 0.0_DP
             call internal_extract_box(local_box_fine%core_val, &
                  local_box_fine, iat, core_dens_fine, i_need_box, &
                  deposit_extract_buffer)
             local_box_fine%core_pseudoval = 0.0_DP
             call internal_extract_box(local_box_fine%core_pseudoval, &
                  local_box_fine, iat, core_pseudodens_fine, &
                  i_need_box, deposit_extract_buffer)

             ! lpl: Calculate weighting factors for atoms that exist on
             !      this proc
             if (i_need_box) then

                ! lpl: Don't bother calculating if there is no density
                !      to be partitioned
                if ( ref_dens(mdl%elements(par%orig_atom(iat))%species_number)&
                     %eff_nuc_charge == mdl%elements(par%orig_atom(iat))&
                     %atomic_number ) cycle

                call internal_core_iterator(rad_dens(iat), local_box_fine, &
                     core_popn(iat))

             end if ! END if (i_need_box)

          end do ! END calculate weighting factors for each atom on this proc

          ! lpl: Finished with all atoms here
          core_popn = vox_vol*core_popn
          call comms_barrier

          ! lpl: Update total core population for each atom
          call comms_reduce('SUM',core_popn)
          dq = maxval(abs(core_popn - atom_popn_prev))

          ddens_max = internal_ddens_max(rad_dens)

          if (pub_debug_on_root) then
          ! lpl: Write out atomic charges for every iteration
             write(debug_unit,'(a)',advance='no') ' DDEC_core_dQ'
             do iat=1,par%nat
                write(debug_unit,'(1x,f11.7)',advance='no') &
                     core_popn(par%orig_atom(iat))
             end do
             write(debug_unit,'(a)') ''
          end if

          ! lpl: Print DDEC loop header on the first iteration
          if(ddec_it == 1) then
             if(pub_on_root) then
                write(stdout,'(/a)') &
                     '-------------------------------------------&
                     &----------------------'
                write(stdout,'(a)') &
                     '                           DDEC Core'
                write(stdout,'(a)') &
                     '-------------------------------------------&
                     &----------------------'
                write(stdout,'(a)') &
                     '    Iteration   Max Change: Population (e)   Density (e)'
                write(stdout,'(a)') &
                     '-------------------------------------------&
                     &----------------------'
             end if
          end if

          ! lpl: Print out info of the current iteration
          if(pub_on_root) &
               write(stdout,'(6x,I6,5x,2(5x,f18.10))') ddec_it, dq, ddens_max

          ! lpl: Check convergence of core charges
          if (dq < pub_ddec_conv_thresh) then
             dq_cons = dq_cons + 1
          else
             dq_cons = 0 ! Reset if dQ > threshold
          end if

          if (ddens_max < pub_ddec_conv_thresh) then
             ddens_max_cons = ddens_max_cons + 1
          else
             ddens_max_cons = 0
          end if

          ! lpl: Terminate loop upon charge convergence
          if (dq_cons >= dq_cons_num .and. ddens_max_cons >= dq_cons_num) then
             if(pub_on_root) then
                write(stdout,'(a)') '-------------------------------&
                     &----------------------------------'
                write(stdout,'(a)') ' DDEC core density converged'
             end if

             exit
          end if

          ! lpl: Update radial core density profiles across all procs. Weighting
          !      factors 'part_dens' are replaced by the new weighting factors
          !      in 'new_part_dens' here, and no reference density is employed
          !      as we are using the ISA scheme (last flag of
          !      'internal_update_rad_dens' set to .false.)
          if ( ddec_it >= pub_ddec_core_maxit ) exit
          call internal_update_rad_dens_core(rad_dens,pub_ddec_reshape_dens)

          if (pub_debug_on_root) then
          ! Write core density profiles after each iteration
          do orig_iat=1,par%nat
             iat = par%distr_atom(orig_iat)
             write(output_filename,'(2(a,I5.5))') &
                  'core_dens_atom_',orig_iat,'_iter_',ddec_it
             call ddec_write_dens(rad_dens(iat)%core_dens, &
                  rad_dens(iat)%rad, rad_dens(iat)%npts, &
                  trim(adjustl(output_filename)))
          end do
          end if

          call comms_barrier

       end do ! END DDEC iterative loop

       q_core_total = sum(core_popn)

       if (pub_debug_on_root) then
       ! lpl: Write out core population for every iteration
          write(debug_unit,'(/a)') ' Final DDEC Core dQ:'
          do iat=1,par%nat
             write(debug_unit,'(1x,a4,1x,I5,2x,f11.7)') &
                  mdl%elements(par%orig_atom(iat))%species_id, &
                  par%orig_atom(iat), core_popn(par%orig_atom(iat))
          end do
         write(debug_unit,'(a)') '--------------------'
         write(debug_unit,'(a,1x,f11.7)') ' Total:', q_core_total
         write(debug_unit,'(a)') '--------------------'
       end if

       call comms_barrier

       if (pub_debug_on_root)  write(stdout,'(a)') &
         ' DEBUG: ddec_main: Leaving core iterator'

    end if ! END if ( include_coredens )

  !----------------------------------------------------------------------------!
  ! lpl: Core occupancy correction iterative loop                              !
  !----------------------------------------------------------------------------!

    if ( include_coredens .and. pub_ddec_core_correction ) then

    if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: ddec_main: Entering core correction'

    ! lpl: Store uncorrected core population temporarily in 'ddec_core_popn'
    !      for reporting purposes
    ddec_core_popn = core_popn

    ! lpl: Store the difference between the computed Cartesian charge and
    !      expected ionic charge into 'core_diff' (this array is called
    !     'nA_correction' in CHARGEMOL)
    core_diff = 0.0_DP
    do iat=1,par%nat

       ! lpl: 'atom_popn' temporarily stores core charge (Z - ion_charge)
       atom_popn(iat) = mdl%elements(par%orig_atom(iat))%atomic_number - &
            ref_dens(mdl%elements(par%orig_atom(iat))%species_number)%eff_nuc_charge
       core_diff(iat) = atom_popn(iat) - core_popn(iat)

       ! Print the uncorrected core populations
       if (pub_debug_on_root) write(debug_unit,'(a,1x,I4,1x,f14.7)') &
            'UNCORRECTED CORE POPN '//&
            mdl%elements(par%orig_atom(iat))%species_id, &
            par%orig_atom(iat), core_popn(iat)

    end do
    dq = maxval(abs(core_diff))

    ! lpl: Use the Cartesian density array 'pseudodensity_fine'  as a buffer
    !      array to store the new, corrected core densities. This is pointed
    !      to by the pointer 'total_pseudodens_fine'.
    total_pseudodens_fine => pseudodensity_fine

    ! lpl: Perform iterative correction
    do ddec_it=1,pub_ddec_core_corr_maxit

       ! lpl: Stores the 'corrected core density' promolecular density
       total_pseudodens_fine = 0.0_DP

       ! lpl: Compute a 'corrected core density' promolecular density
       do lc_iat=1,par%max_atoms_on_proc

          ! lpl: Whether we need to extract density box for this atom
          if (lc_iat <= par%num_atoms_on_proc(pub_my_proc_id)) then
             i_need_box = .true.
             i_have_box = .true.
             iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
          else
             i_need_box = .false.
             i_have_box = .false.
             iat = par%first_atom_on_proc(pub_my_proc_id) ! Prevents segfault
          end if

          ! lpl: Skip atoms with no core charges
          if ( atom_popn(iat) < pub_ddec_zero_thresh ) then
             i_need_box = .false.
             i_have_box = .false.
          end if

          ! lpl: Extract global core density and global promolecular density
          local_box_fine%core_val = 0.0_DP
          call internal_extract_box(local_box_fine%core_val, local_box_fine, &
               iat, core_dens_fine, i_need_box, deposit_extract_buffer)
          local_box_fine%core_pseudoval = 0.0_DP
          call internal_extract_box(local_box_fine%core_pseudoval, &
               local_box_fine, iat, core_pseudodens_fine, &
               i_need_box, deposit_extract_buffer)

          if ( i_need_box ) call internal_core_corr_iterator(local_box_fine, &
               rad_dens(iat), core_diff(iat))

          ! lpl: Deposit the corrected core densities into
          !      'total_pseudodens_fine'
          call internal_deposit_box(total_pseudodens_fine, &
               local_box_fine%aim_val, local_box_fine, iat, &
               i_have_box, deposit_extract_buffer)

       end do ! END do lc_iat=1,par%max_atoms_on_proc
       call comms_barrier

       ! lpl: Update 'core_dens_fine' with the corrected core density
       ! lpl: Swap pointers if this is too slow
       core_dens_fine = total_pseudodens_fine

       ! lpl: Calculate new weighting factors for each atom on this proc
       !      based on the 'corrected core densities' promolecular density
       !      that resides in 'core_dens_fine', update 'core_popn' and
       !      the core promolecular density
       core_popn = 0.0_DP
       do lc_iat=1,par%max_atoms_on_proc

          ! lpl: Whether we need to extract density box for this atom
          if (lc_iat <= par%num_atoms_on_proc(pub_my_proc_id)) then
             i_need_box = .true.
             iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
          else
             i_need_box = .false.
             iat = par%first_atom_on_proc(pub_my_proc_id) ! Prevents segfault
          end if

          ! lpl: Don't bother calculating if there is no core density
          if ( atom_popn(iat) < pub_ddec_zero_thresh ) i_need_box = .false.

          ! lpl: Extract core pseudodensity and corrected core density into
          !      'local_box_fine'
          local_box_fine%core_val = 0.0_DP
          call internal_extract_box(local_box_fine%core_val, local_box_fine, &
               iat, core_dens_fine, i_need_box, deposit_extract_buffer)
          local_box_fine%core_pseudoval = 0.0_DP
          call internal_extract_box(local_box_fine%core_pseudoval, &
               local_box_fine, iat, core_pseudodens_fine, &
               i_need_box, deposit_extract_buffer)

          ! lpl: Calculate new weighting factors based on the corrected
          !      core density
          if (i_need_box) call internal_core_iterator(rad_dens(iat), &
               local_box_fine, core_popn(iat))

       end do ! END calculate weighting factors for each atom on this proc
       call comms_barrier
       call comms_reduce('SUM',core_popn)

       ! lpl: Also print DDEC iterative loop header on 1st iteration
       if(ddec_it == 1) then
          if(pub_on_root) then
             ! lpl: Also print DDEC iterative loop header on 1st iteration
             write(stdout,'(/a)') '-------------------------------------------&
                  &----------------------'
             write(stdout,'(a)') '                     DDEC Core Correction'
             write(stdout,'(a)') '-------------------------------------------&
                  &----------------------'
             write(stdout,'(a)') &
                  '    Iteration     Max Population Difference (e)'
             write(stdout,'(a)') '-------------------------------------------&
                  &----------------------'
             write(stdout,'(6x,I6,5x,5x,f18.10)') ddec_it-1, dq
          end if
       end if

       ! lpl: Compute new correction factors and store in 'core_diff'
       core_diff = atom_popn - vox_vol*core_popn
       dq = maxval(abs(core_diff))

       ! lpl: Print out info on the current DDEC iteration
       if(pub_on_root) &
            write(stdout,'(6x,I6,5x,5x,f18.10)') ddec_it, dq

       if ( dq < pub_ddec_conv_thresh ) then
          if(pub_on_root) then
             write(stdout,'(a)') '-------------------------------&
                  &----------------------------------'
             write(stdout,'(a)') ' DDEC core density correction converged'
          end if
          exit
       end if

       ! lpl: Update partial core densities (in essence copying 'new_part_dens'
       !      into 'part_dens', reshaping density, and bcast to all procs
       call internal_update_rad_dens_core(rad_dens, pub_ddec_reshape_dens)

       ! lpl: Update core promolecular density (only 'core_pseudodens_fine'
       !      will be updated here) for ALL atoms
       core_pseudodens_fine = 0.0_DP
       call internal_deposit_box_all(core_pseudodens_fine, rad_dens, &
            'core_dens', local_box_fine, 'aim_val', deposit_extract_buffer)

    end do ! END do ddec_it=1,ddec_max_corr
    call comms_barrier

    q_core_total = vox_vol*sum(core_popn)

    ! lpl: Nullify pointers (not really required)
    nullify(total_pseudodens_fine)
    call comms_barrier

       if (pub_debug_on_root)  then
          write(debug_unit,'(a)') ' CORRECTED CORE POPN:'
          do orig_iat=1,par%nat
             iat = par%distr_atom(orig_iat)
             write(debug_unit,'(a,2x,I5,2x,f14.7)') &
                  mdl%elements(orig_iat)%species_id, orig_iat, &
                  vox_vol*core_popn(iat)
          end do

          write(stdout,'(a)') ' DEBUG: ddec_main: Leaving core iterator'
       end if

    else
       ! lpl: Rescale for later (if core correction wasn't called this
       !      would not have been done)
       core_popn = core_popn/vox_vol
    end if ! END if ( pub_ddec_core_correction )

  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lpl: Write out core charges                                                !
  !----------------------------------------------------------------------------!

    if ( include_coredens ) then

       ! lpl: Write charges
       if(pub_on_root) then
          write(stdout,'(a)') '------------------------------------------------'
          write(stdout,'(a)') '                  Core Charges'
          write(stdout,'(a)') '------------------------------------------------'
          write(stdout,'(a)') '   Atom       Uncorrected (e)     Corrected (e) '
          write(stdout,'(a)') '------------------------------------------------'
          do orig_iat=1,par%nat
             iat = par%distr_atom(orig_iat)
             isp = mdl%elements(orig_iat)%species_number
             if ( pub_ddec_core_correction ) then
                write(stdout,'(1x,a4,1x,I5,2x,3(5x,f11.7))') &
                     mdl%elements(orig_iat)%species_id, orig_iat, &
                     ddec_core_popn(iat), vox_vol*core_popn(iat)
             else
                write(stdout,'(1x,a4,1x,I5,2x,3(5x,f11.7))') &
                     mdl%elements(orig_iat)%species_id, orig_iat, &
                     vox_vol*core_popn(iat), 0.0_DP
             end if
          end do
          write(stdout,'(a)') '------------------------------------------------'
          write(stdout,'(a,f14.7)') ' Total uncorrected (e) : ', &
               sum(ddec_core_popn)
          write(stdout,'(a,f14.7)') ' Total corrected (e)   : ', q_core_total
          write(stdout,'(a)') '------------------------------------------------'
       end if
       call comms_bcast(pub_root_proc_id,ddec_charge)

    end if

  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lpl: Main DDEC iterative loop                                              !
  !----------------------------------------------------------------------------!

    if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: ddec_main: Entering valence iterator'

    ! lpl: Re-establish all densities as some pointers may have been
    !      used to point to irrelevant bufferes in previous calculations
    total_dens_fine => density_fine
    total_pseudodens_fine => pseudodensity_fine
    if ( include_refdens ) then
       ref_pseudodens_fine => ref_pseudodensity_fine
       Yavg_pseudodens_fine => Yavg_pseudodensity_fine
    end if

    ! lpl: Add core density to valence density
    if ( include_coredens ) then
       total_dens_fine = total_dens_fine + core_dens_fine
       call comms_barrier
    else
       if ( any(core_popn /= 0.0_DP) ) call utils_abort(" ERROR in ddec_main:&
            & any(core_popn /= 0.0_DP) when 'include_coredens = F'")
    end if

    ! lpl: Initialize weighting factors as neutral IH reference densities
    core_diff = 0.0_DP
    ddec_core_popn = 0.0_DP
    do iat=1,par%nat
       rad_dens(iat)%ref_dens = &
            ref_dens(mdl%elements(par%orig_atom(iat))&
            %species_number)%ion_dens(0)%dens

       ! lpl: Optionally initialize part_dens as 1.0_DP
       if ( pub_ddec_refdens_init ) then
          rad_dens(iat)%total_part_dens = rad_dens(iat)%ref_dens
       else
          rad_dens(iat)%total_part_dens = 1.0_DP ! Optionally initialize as 1.0
       end if

       ! lpl: Initialize Yavg density if required
       if ( include_refdens ) &
            rad_dens(iat)%Yavg_dens = rad_dens(iat)%ref_dens

       ! lpl: Initialize the exp_factor values correctly
       rad_dens(iat)%exp_factor = &
            rad_dens(iat)%exp_factor*valence_exp_factor

       ! lpl: MOD#05 - Get expected core population
       !      (Z - ref_dens%eff_nuc_charge) too
       isp =  mdl%elements(par%orig_atom(iat))%species_number
       orig_iat = par%orig_atom(iat)
       ddec_core_popn(iat) = ref_dens(isp)%ddec_nuc_charge - &
            ref_dens(isp)%eff_nuc_charge
    end do

    ddens_max = 0.0_DP
    atom_popn = 0.0_DP
    val_popn  = 0.0_DP
    dq = pub_ddec_conv_thresh
    dq_cons = 0 ! Consecutive number of iterations with dQ < threshold
    ddens_max_cons = 0

    ! lpl: DDEC iterative loop
    do ddec_it=1,pub_ddec_maxit

    ! Write radial density weighting for every iteration
    if(pub_debug_on_root) then
       do iat=1,par%nat
          output_filename = ''
          write(output_filename,'(a,I5.5,a,I5.5)') 'rad_w_atom_', &
               par%orig_atom(iat),'_iter_',ddec_it-1
          call ddec_write_dens(rad_dens(iat)%total_part_dens, &
               rad_dens(iat)%rad, rad_dens(iat)%npts, &
               trim(adjustl(output_filename)))
       end do
    end if

       ! lpl: Reset DDEC charges
       atom_popn_prev = val_popn
       atom_popn = 0.0_DP
       if ( include_coredens ) core_diff = 0.0_DP

       ! lpl: Reset and update promolecular densities
       total_pseudodens_fine = 0.0_DP
       call internal_deposit_box_all(total_pseudodens_fine, rad_dens, &
            'total_part_dens', local_box_fine, 'aim_val', &
            deposit_extract_buffer)
       if ( include_refdens ) then
          ref_pseudodens_fine   = 0.0_DP
          call internal_deposit_box_all(ref_pseudodens_fine, rad_dens, &
               'ref_dens', local_box_fine, 'aim_val', deposit_extract_buffer)
          Yavg_pseudodens_fine  = 0.0_DP
          call internal_deposit_box_all(Yavg_pseudodens_fine, rad_dens, &
               'Yavg_dens', local_box_fine, 'aim_val', deposit_extract_buffer)
       end if

       ! lpl: Calculate new weighting factors for each atom on this proc
       do lc_iat=1,par%max_atoms_on_proc

          ! lpl: Whether we need to extract density box for this atom
          if (lc_iat <= par%num_atoms_on_proc(pub_my_proc_id)) then
             i_need_box = .true.
             iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
          else
             i_need_box = .false.
             iat = par%first_atom_on_proc(pub_my_proc_id) ! Prevents segfault
          end if

          ! lpl: Extract all relevant global densities into 'local_box_fine'
          local_box_fine%total_val = 0.0_DP
          call internal_extract_box(local_box_fine%total_val, local_box_fine, &
               iat, total_dens_fine, i_need_box, deposit_extract_buffer)
          local_box_fine%total_pseudoval = 0.0_DP
          call internal_extract_box(local_box_fine%total_pseudoval, &
               local_box_fine, iat, total_pseudodens_fine, &
               i_need_box, deposit_extract_buffer)

          if ( include_refdens ) then
             local_box_fine%ref_pseudoval = 0.0_DP
             call internal_extract_box(local_box_fine%ref_pseudoval, &
                  local_box_fine, iat, ref_pseudodens_fine, &
                  i_need_box, deposit_extract_buffer)
             local_box_fine%Yavg_pseudoval = 0.0_DP
             call internal_extract_box(local_box_fine%Yavg_pseudoval, &
                  local_box_fine, iat, Yavg_pseudodens_fine, &
                  i_need_box, deposit_extract_buffer)
          end if

          if ( include_coredens ) then
             local_box_fine%core_val = 0.0_DP
             call internal_extract_box(local_box_fine%core_val, &
                  local_box_fine, iat, core_dens_fine, &
                  i_need_box, deposit_extract_buffer)
             local_box_fine%core_pseudoval = 0.0_DP
             call internal_extract_box(local_box_fine%core_pseudoval, &
                  local_box_fine, iat, core_pseudodens_fine, &
                  i_need_box, deposit_extract_buffer)
          end if

          ! lpl: Calculate weighting factors for atoms that exist on this proc
          if (i_need_box) call internal_valence_iterator(rad_dens(iat), &
               local_box_fine, atom_popn(iat), core_diff(iat))

       end do ! END calculate weighting factors for each atom on this proc

       ! lpl: Finished all atoms here
       call comms_barrier

       ! lpl: Update new total charges
       call comms_reduce('SUM',atom_popn)
       call comms_reduce('SUM',core_diff)

       ! lpl: Compute valence population. If no coredens is used
       !      then there is no need to subtract core density (and
       !      'core_popn' should be zero anyway
       if ( include_coredens ) then
          val_popn = vox_vol*(atom_popn - (core_popn - core_diff))
       else
          val_popn = vox_vol*atom_popn
       end if

       ! lpl: Renormalize charges to N (due to the DDEC cutoff radius, total
       !      partitioned charge will not be exact) and calculate max dQ change
       q_total = sum(val_popn)
       q_renorm = total_grid_dens/q_total  ! Renormalize only valence popn
       val_popn = q_renorm*val_popn        ! Scale only valence popn
       q_spill = total_grid_dens - q_total ! Spilling for valence popn
       q_total = q_renorm*q_total
       dq = maxval(abs(val_popn - atom_popn_prev))
       ddens_max = internal_ddens_max(rad_dens)

       ! lpl: Add back expected core charges
       atom_popn = val_popn + ddec_core_popn

       if (pub_debug_on_root) then
          write(debug_unit,'(a,/a,1x,I5,/a)') '----------', &
               'iter', ddec_it, '----------'
          write(debug_unit,'(1x,a5,5(1x,a15))') &
               'IAT','ATOM_POPN', 'CORE_POPN', &
               'CORE_DIFF', 'VAL_POPN', 'DDEC_CORE_POPN'
          do orig_iat=1,par%nat
             iat = par%distr_atom(orig_iat)
             write(debug_unit,'(1x,I5,5(1x,f14.7))') &
                  orig_iat, atom_popn(iat), &
                  vox_vol*core_popn(iat), vox_vol*core_diff(iat), &
                  val_popn(iat), ddec_core_popn(iat)
          end do
          write(debug_unit,'(a)') '----------'
       end if

       ! lpl: Print out classical Hirshfeld charges on the 1st DDEC iteration
       if(ddec_it == 1) then
          if(pub_on_root) then
             if(pub_ddec_hirshfeld .and. pub_ddec_refdens_init) then
                write(stdout,'(a)') &
                     '------------------------------------------------'
                write(stdout,'(a)') &
                     '           Classical Hirshfeld Charges          '
                write(stdout,'(a)') &
                     '------------------------------------------------'
                write(stdout,'(a)') &
                     '   Atom        Population (e)      Charge (e)   '
                write(stdout,'(a)') &
                     '------------------------------------------------'

                do orig_iat=1,par%nat
                   iat = par%distr_atom(orig_iat)
                   write(stdout,'(1x,a4,1x,I5,7x,f11.7,5x,f11.7)') &
                        mdl%elements(orig_iat)%species_id, orig_iat, &
                        atom_popn(iat), &
                        ref_dens(mdl%elements(orig_iat)%species_number)%&
                        ddec_nuc_charge - atom_popn(iat)
                end do
                write(stdout,'(a)') &
                     '------------------------------------------------'
                write(stdout,'(a,f14.7)') ' Total charge (e): ',q_total
                write(stdout,'(a)') &
                     '------------------------------------------------'
             else if(pub_ddec_hirshfeld .and. .not. pub_ddec_refdens_init) then
                write(stdout,'(a,/a)') ' WARNING: pub_ddec_refdens_init set to&
                     & false.', 'Classical Hirshfeld charges will not be&
                     & available.'
             end if

             ! lpl: Also print DDEC iterative loop header on 1st iteration
             write(stdout,'(/a)') '-------------------------------------------&
                  &----------------------'
             write(stdout,'(a)') '                              DDEC'
             write(stdout,'(a)') '-------------------------------------------&
                  &----------------------'
             write(stdout,'(a)') &
                  '    Iteration   Max Change: Population (e)   Density (e)'

             write(stdout,'(a)') '-------------------------------------------&
                  &----------------------'
          end if ! END if(pub_on_root)
       end if ! END if(ddec_it == 1) then

       call comms_barrier

       ! lpl: Print out info on the current DDEC iteration
       if(pub_on_root) &
            write(stdout,'(6x,I6,5x,3(5x,f18.10))') ddec_it, dq, ddens_max

       ! lpl: Check consecutive convergence
       if (dq < pub_ddec_conv_thresh) then
          dq_cons = dq_cons + 1
       else
          dq_cons = 0 ! Reset if dQ > threshold
       end if

       if (ddens_max < pub_ddec_conv_thresh) then
          ddens_max_cons = ddens_max_cons + 1
       else
          ddens_max_cons = 0
       end if

       if (pub_debug) then
       ! Write all radial density profiles
       do lc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
          iat  = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
          output_filename = ''
          write(output_filename,'(a,I5.5,a,I5.5)') 'new_part_dens_atom_', &
               par%orig_atom(iat),'_iter_',ddec_it
          call ddec_write_dens(rad_dens(iat)%new_part_dens, &
               rad_dens(iat)%rad, rad_dens(iat)%npts, &
               trim(adjustl(output_filename)),pub_my_proc_id)

          write(output_filename,'(a,I5.5,a,I5.5)') 'rad_npts_atom_', &
               par%orig_atom(iat),'_iter_',ddec_it
          call ddec_write_dens(real(rad_dens(iat)%rad_npts,kind=DP), &
               rad_dens(iat)%rad, rad_dens(iat)%npts, &
               trim(adjustl(output_filename)),pub_my_proc_id)

          if ( include_refdens ) then
             output_filename = ''
             write(output_filename,'(a,I5.5,a,I5.5)') 'Yavg_dens_atom_', &
                  par%orig_atom(iat),'_iter_',ddec_it-1
             call ddec_write_dens(rad_dens(iat)%Yavg_dens, &
                  rad_dens(iat)%rad, rad_dens(iat)%npts, &
                  trim(adjustl(output_filename)),pub_my_proc_id)

             output_filename = ''
             write(output_filename,'(a,I5.5,a,I5.5)') 'Yavg_sqrt_atom_', &
                  par%orig_atom(iat),'_iter_',ddec_it
             call ddec_write_dens(rad_dens(iat)%Yavg_sqrt, &
                  rad_dens(iat)%rad, rad_dens(iat)%npts, &
                  trim(adjustl(output_filename)),pub_my_proc_id)

             output_filename = ''
             write(output_filename,'(a,I5.5,a,I5.5)') 'Yavg_ratio_atom_', &
                  par%orig_atom(iat),'_iter_',ddec_it
             call ddec_write_dens(rad_dens(iat)%Yavg_ratio, &
                  rad_dens(iat)%rad, rad_dens(iat)%npts, &
                  trim(adjustl(output_filename)),pub_my_proc_id)

             output_filename = ''
             write(output_filename,'(a,I5.5,a,I5.5)') 'ih_ratio_atom_', &
                  par%orig_atom(iat),'_iter_',ddec_it
             call ddec_write_dens(rad_dens(iat)%ih_ratio, &
                  rad_dens(iat)%rad, rad_dens(iat)%npts, &
                  trim(adjustl(output_filename)),pub_my_proc_id)

             output_filename = ''
             write(output_filename,'(a,I5.5,a,I5.5)') 'ref_dens_atom_', &
                  par%orig_atom(iat),'_iter_',ddec_it
             call ddec_write_dens(rad_dens(iat)%ref_dens, &
                  rad_dens(iat)%rad, rad_dens(iat)%npts, &
                  trim(adjustl(output_filename)),pub_my_proc_id)
          end if
       end do
       end if

       ! lpl: Terminate loop upon charge convergence
       if ( dq_cons >= dq_cons_num .and. ddens_max_cons >= dq_cons_num &
            .and. ddec_it > 9 ) then
          if(pub_on_root) then
             write(stdout,'(a)') '-------------------------------&
                  &----------------------------------'
             write(stdout,'(a)') ' DDEC density converged'
          end if

          exit
       end if

       ! lpl: Update radial densities across all procs. Weighting
       !      factors 'part_dens' are replaced by the new weighting factors
       !      in 'new_part_dens' here
       if ( ddec_it >= pub_ddec_maxit ) then
          exit ! Halt updating if 'ddec_it' has reached max
       else if ( ddec_it <= 3 ) then
          call internal_update_rad_dens(rad_dens,ref_dens,atom_popn, &
               include_refdens,.false.)
       else
          call internal_update_rad_dens(rad_dens,ref_dens,atom_popn, &
               include_refdens,pub_ddec_reshape_dens)
       end if

    end do ! END DDEC iterative loop

    ! lpl: Check convergence after DDEC loop has terminated
    if (ddec_it >= pub_ddec_maxit .and. dq >= pub_ddec_conv_thresh &
         .and. dq_cons < dq_cons_num .and. ddens_max_cons < dq_cons_num) then
       if(pub_on_root) &
          write(stdout,'(a,I5,a)') ' WARNING: DDEC charges not converged &
               &after ', pub_ddec_maxit, ' iterations'
    end if

    if ( pub_on_root ) write(stdout,'(a)') '-------------------------------&
         &----------------------------------'

    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: ddec_main: Leaving valence iterator'

  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lpl: Write out converged DDEC charges                                      !
  !----------------------------------------------------------------------------!

    ! lpl: Write charges
    if(pub_on_root) then
       write(stdout,'(a)') '------------------------------------------------'
       write(stdout,'(a,f4.2,a)') '             DDEC Charges (X=', &
            pub_ddec_IH_frac, ')              '
       write(stdout,'(a)') '------------------------------------------------'
       write(stdout,'(a)') '   Atom        Population (e)      Charge (e)   '
       write(stdout,'(a)') '------------------------------------------------'
       do orig_iat=1,par%nat
          iat = par%distr_atom(orig_iat)
          isp = mdl%elements(orig_iat)%species_number
          ddec_charge(iat) = ref_dens(isp)%ddec_nuc_charge - atom_popn(iat)
          write(stdout,'(1x,a4,1x,I5,2x,2(5x,f11.7))') &
               mdl%elements(orig_iat)%species_id, orig_iat, atom_popn(iat), &
                    ddec_charge(iat)
       end do
       write(stdout,'(a)') '------------------------------------------------'
       write(stdout,'(a,f14.7)') ' Total charge (e) : ',q_total
       write(stdout,'(a,f14.7)') ' Spilling (e)     : ',q_spill
       write(stdout,'(a)') '------------------------------------------------'
    end if
    call comms_bcast(pub_root_proc_id,ddec_charge)

    ! lpl: QC test
    if ( pub_on_root ) then
       if ( pub_print_qc ) then
          do orig_iat=1,min(3,par%nat)
             iat = par%distr_atom(orig_iat)
             output_filename=''
             write(output_filename,'(a,I0,a)') 'atom_',orig_iat,'_ddec_popn'
             call utils_qc_print(trim(adjustl(output_filename)), &
                  ddec_charge(iat))
             ! lpl: Number of intervals (reuse variable 'dq_cons')
             dq_cons = max(1,rad_dens(iat)%npts/5)
             do it=(1+min(1,rad_dens(iat)%npts/4)),min(4,rad_dens(iat)%npts)
                output_filename = ''
                write(output_filename,'(2(a,I0))') 'atom_',orig_iat,'_rad_',it
                call utils_qc_print( trim(adjustl(output_filename))//&
                     '_rad',rad_dens(iat)%dr*(dq_cons*it - 0.5_DP) )
                call utils_qc_print( trim(adjustl(output_filename))//&
                     '_dens',rad_dens(iat)%total_part_dens(it*dq_cons) )
             end do
          end do
       end if
       output_filename = ''
    end if
    call comms_barrier

  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lpl: Calculate DDEC AIM effective decay exponents                          !
  !----------------------------------------------------------------------------!

    if ( pub_ddec_eff_decay_exp ) then

       if( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: ddec_main: Entering properties #1'

       call utils_assert(pub_ddec_eff_decay_rmin >= 0.0_DP, &
            'pub_ddec_eff_decay_rmin mist be >= 0.0_DP')

       allocate(ddec_eff_decay_exp(par%nat,3),stat=ierr)
       call utils_alloc_check('ddec_main','ddec_eff_decay_exp',ierr)

       ddec_eff_decay_exp = 0.0_DP

       do lc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
          iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

          call internal_eff_decay_exp(ddec_eff_decay_exp(iat,1), &
               ddec_eff_decay_exp(iat,2), ddec_eff_decay_exp(iat,3), &
               rad_dens(iat)%new_part_dens, rad_dens(iat), &
               pub_ddec_eff_decay_rmin)
       end do
       call comms_reduce('SUM',ddec_eff_decay_exp)

       if ( pub_on_root ) then
          write(stdout,'(a)') &
               '---------------------------------------------------------'
          write(stdout,'(a)') '              AIM Effective Decay Exponents'
          write(stdout,'(a,f10.5,a)') &
               '  Least-squares fit to exp(a-br) for r > ', &
               pub_ddec_eff_decay_rmin,' a0'
          write(stdout,'(a)') &
               '---------------------------------------------------------'
          write(stdout,'(a)') ' isp  iat   a              b              R^2'
          write(stdout,'(a)') &
               '---------------------------------------------------------'

          do orig_iat=1,par%nat
             iat = par%distr_atom(orig_iat)
             write(stdout,'(1x,a4,1x,I5,3(1x,f14.7))') &
                  mdl%elements(orig_iat)%species_id, orig_iat, &
                  ddec_eff_decay_exp(iat,1), &
                  -1.0_DP*ddec_eff_decay_exp(iat,2), ddec_eff_decay_exp(iat,3)
          end do

          write(stdout,'(a)') &
               '---------------------------------------------------------'
       end if

       deallocate(ddec_eff_decay_exp,stat=ierr)
       call utils_dealloc_check('ddec_main','ddec_eff_decay_exp',ierr)

       if( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: ddec_main: Leaving properties #1'

    end if

  !----------------------------------------------------------------------------!
  ! lpl: Calculate DDEC AIM multipoles and moments                             !
  !----------------------------------------------------------------------------!

! lpl: c.f. COMMENT#07
    if ( pub_ddec_multipole .or. pub_ddec_moment > 0 ) then

       if( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: ddec_main: Entering properties #2'

       ! lpl: Allocate dipole and quadrupole vector arrays
       if ( pub_ddec_multipole ) then
          allocate(ddec_dipoles(par%nat,3),stat=ierr)
          call utils_alloc_check('ddec_main','ddec_dipoles',ierr)
          allocate(ddec_quadrupoles(par%nat,5),stat=ierr)
          call utils_alloc_check('ddec_main','ddec_quadrupoles',ierr)
          ddec_dipoles = 0.0_DP
          ddec_quadrupoles = 0.0_DP

          call comms_barrier
       end if

       ! lpl: c.f. COMMENT#07
       if ( pub_ddec_moment >= 0 ) then
          allocate(ddec_moments(par%nat,3),stat=ierr)
          call utils_alloc_check('ddec_main','ddec_moments',ierr)
          ddec_moments = 0.0_DP

          call comms_barrier
       end if

       call comms_barrier


       do lc_iat=1,par%max_atoms_on_proc

          ! lpl: Whether we need to extract density box for this atom
          if (lc_iat <= par%num_atoms_on_proc(pub_my_proc_id)) then
             i_need_box = .true.
             iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
          else
             i_need_box = .false.
             iat = par%first_atom_on_proc(pub_my_proc_id) ! Prevents segfault
          end if

          ! lpl: Extract global density and global promolecular density
          local_box_fine%total_val = 0.0_DP
          call internal_extract_box(local_box_fine%total_val, &
               local_box_fine, iat, total_dens_fine, i_need_box, &
               deposit_extract_buffer)
          local_box_fine%total_pseudoval = 0.0_DP
          call internal_extract_box(local_box_fine%total_pseudoval, &
               local_box_fine, iat, total_pseudodens_fine, &
               i_need_box, deposit_extract_buffer)

          ! lpl: Extract core density and promolecular density if required
          ! lpl: MOD#08 - Both multipole and moment calculations work only
          !      with the valence AIM
          if ( include_coredens ) then
             local_box_fine%core_val = 0.0_DP
             call internal_extract_box(local_box_fine%core_val, &
                  local_box_fine, iat, core_dens_fine, i_need_box, &
                  deposit_extract_buffer)
             local_box_fine%core_pseudoval = 0.0_DP
             call internal_extract_box(local_box_fine%core_pseudoval, &
                  local_box_fine, iat, core_pseudodens_fine, &
                  i_need_box, deposit_extract_buffer)
          end if

          ! lpl: Calculate weighting factors for atoms that exist on
          !      this proc and store in 'local_box_fine%aim_val'
          if (i_need_box) then
             ! lpl: Total AIM density
             call internal_aim_density(local_box_fine%aim_val, &
                  local_box_fine%total_val, local_box_fine%total_pseudoval, &
                  local_box_fine, rad_dens(iat)%total_part_dens, rad_dens(iat))

             ! lpl: Subtract core density contribution for multipole calc
             ! lpl: MOD#08 - Both multipole and moment calculations work only
             !      with the valence AIM
             if ( include_coredens ) then
                ! lpl: Core AIM
                call internal_aim_density(local_box_fine%temp_val, &
                     local_box_fine%core_val, &
                     local_box_fine%core_pseudoval, &
                     local_box_fine, rad_dens(iat)%core_dens, rad_dens(iat))

                ! lpl: Remove explicit non-negative enforcement?
                local_box_fine%aim_val = &
                     local_box_fine%aim_val - local_box_fine%temp_val
                call ddec_zero_threshold(local_box_fine%aim_val, 0.0_DP, &
                     'Zero', .true.)
                !djc
             end if


             ! lpl: COMMENT#07 - 'pub_ddec_moment' and 'pub_ddec_moment_order'
             !      have been merged following Danny's advice. To disable
             !      moment calculation, 'pub_ddec_moment' must be set to a
             !      value <0. 0th moment is kept for the sake of testing
             !      the integrand, if a need for it ever arises.
             ! lpl: Calculate moment
             if ( pub_ddec_moment >= 0 ) then
                ! lpl: MOD#11 (Valence) AIM moment from 'local_box_fine%aim_val'
                call internal_moment(ddec_moments(iat,1), &
                     local_box_fine%aim_val, local_box_fine, rad_dens(iat), &
                     pub_ddec_moment)
! lpl: c.f. COMMENT#07
                call internal_moment( aim_volume(iat), local_box_fine%aim_val, &
                     local_box_fine, rad_dens(iat), 3 )


!                     pub_ddec_moment_order)

                ! lpl: MOD#11 (Core) AIM moment from 'local_box_fine%temp_val'
                if ( include_coredens ) then
                    call internal_moment(ddec_moments(iat,3), &
                         local_box_fine%temp_val, local_box_fine, &
                         rad_dens(iat), pub_ddec_moment)
! lpl: c.f. COMMENT#07
!                         rad_dens(iat), pub_ddec_moment_order)

                end if

                ! lpl: Total reference density into 'local_box_fine%temp_val'
                local_box_fine%temp_val = 0.0_DP
                call internal_rad_dens_to_box(local_box_fine%temp_val, &
                     local_box_fine, ref_dens(mdl%elements(par%orig_atom(iat))&
                     %species_number)%ion_dens(0)%dens, rad_dens(iat))

                ! lpl: Total reference density moment
                call internal_moment(ddec_moments(iat,2), &
                     local_box_fine%temp_val, &
                     local_box_fine, rad_dens(iat), pub_ddec_moment)
                ! lpl: c.f. COMMENT#07
                !                     local_box_fine, rad_dens(iat), pub_ddec_moment_order)

             end if

             ! lpl: Calculate multipoles
             if ( pub_ddec_multipole ) then
                call internal_multipole(ddec_dipoles(iat,:), &
                     ddec_quadrupoles(iat,:), &
                     local_box_fine%aim_val, local_box_fine, rad_dens(iat))


                ! aeaa: Section for finding off center point charges
                if(pub_ddec_aniso) then
                   ! aeaa: Arrays of atom centers and element type
                   do ii=1,par%nat
                      atom_ctr_point_array(ii) = mdl%elements(ii)%centre
                      all_element_types(ii) = mdl%elements(ii)%species_id
                   end do

                   ! aeaa: Atom center for atom in loop
                   atom_ctr_point(iat) = mdl%elements(par%orig_atom(iat))%centre&
                        - local_box_fine%box_coord_origin(iat)

                   call ddec_off_site_charges( optimum_posit(:,:,iat), &
                        optimum_charges(:,iat), error_aniso_before(iat), &
                        error_aniso(iat), atom_ctr_point(iat), &
                        atom_ctr_point_array, ddec_charge(iat),&
                        aim_volume(iat) , ierr, local_box_fine%aim_val, &
                        local_box_fine, &
                        mdl%elements(par%orig_atom(iat))%species_id,&
                        all_element_types, par%orig_atom(iat), grid%da1%X, &
                        grid%da2%Y, grid%da3%Z, no_charges(iat),&
                        vdw_radius(iat) )

                   !aeaa: Ensures hydrogen atoms never have extra charges added
                   if(trim(mdl%elements(par%orig_atom(iat))%species_id) == 'H') &
                        then
                      no_charges(iat) = 0
                   end if

                   ! aeaa: Atom centers added to off center positions
                   optimum_posit(1,:,iat) = optimum_posit(1,:,iat) + &
                        ( mdl%elements(par%orig_atom(iat))%centre%X * 0.529177 )

                   optimum_posit(2,:,iat) = optimum_posit(2,:,iat) + &
                        ( mdl%elements(par%orig_atom(iat))%centre%Y * 0.529177 )

                   optimum_posit(3,:,iat) = optimum_posit(3,:,iat) + &
                        ( mdl%elements(par%orig_atom(iat))%centre%Z  * 0.529177 )

                end if
             end if

          end if ! END if (i_need_box)

          call comms_barrier

       end do ! END  do lc_iat=1,par%max_atoms_on_proc

       if ( pub_ddec_multipole ) then

          call comms_reduce('SUM', ddec_dipoles)
          call comms_reduce('SUM', ddec_quadrupoles)
          ddec_dipoles = q_renorm*vox_vol*ddec_dipoles
          ddec_quadrupoles = q_renorm*vox_vol*ddec_quadrupoles

          call comms_reduce('SUM', optimum_posit )
          call comms_reduce('SUM', optimum_charges )
          call comms_reduce('SUM', error_aniso )
          call comms_reduce('SUM', error_aniso_before )
          call comms_reduce('SUM', no_charges )
          call comms_reduce('SUM', vdw_radius )
          call comms_reduce('SUM', aim_volume )
       end if

       call comms_barrier

       !lpl: c.f. COMMENT#07

       if ( pub_ddec_moment > 0 ) then
          call comms_reduce('SUM', ddec_moments)
          ddec_moments = q_renorm*vox_vol*ddec_moments
       end if

       ! lpl: QC test
       if ( pub_on_root ) then
          if ( pub_print_qc ) then
             do orig_iat=1,min(3,par%nat)
                iat = par%distr_atom(orig_iat)
                if ( pub_ddec_multipole ) then
                   output_filename = ''
                   write(output_filename,'(2(a,I0))') &
                        'atom_',orig_iat,'_dipole_',orig_iat
                   call utils_qc_print( trim(adjustl(output_filename)), &
                        ddec_dipoles(iat,orig_iat) )
                   output_filename = ''
                   write(output_filename,'(2(a,I0))') &
                        'atom_',orig_iat,'_quadrupole_',orig_iat
                   call utils_qc_print( trim(adjustl(output_filename)), &
                        ddec_quadrupoles(iat,orig_iat) )
                   output_filename = ''
                end if

! lpl: c.f. COMMENT#07
                if ( pub_ddec_moment > 0 ) then
                   write(output_filename,'(2(a,I0))') &
                        'atom_',orig_iat,'_moment_',pub_ddec_moment
! lpl: c.f. COMMENT#07
!                        'atom_',orig_iat,'_moment_',pub_ddec_moment_order
                   call utils_qc_print( trim(adjustl(output_filename)), &
                        ddec_moments(iat,2))
                end if
             end do
          end if
          output_filename = ''
       end if

       ! lpl: Print multipoles and/or moments (in a.u.)
       if(pub_on_root) then
          write(stdout,'(a)') '------------------------------------------&
               &-------------------------------------------------'
          ! lpl: Print multipoles
          if ( pub_ddec_multipole ) then
             write(stdout,'(a,f4.2,a)') '                                &
                  &DDEC Multipoles (X=',pub_ddec_IH_frac, &
                  ')                                   '
             write(stdout,'(a)') '------------------------------------------&
                  &-------------------------------------------------'
             write(stdout,'(a)') '  Atom        Dx        Dy        Dz      &
                  &  Qxy       Qxz       Qyz     Qx2-y2    Q3z2-r2  '
             write(stdout,'(a)') '------------------------------------------&
                  &-------------------------------------------------'

             do orig_iat=1,par%nat
                iat = par%distr_atom(orig_iat)
                write(stdout,'(1x,a4,I5,8(1x,f9.6))') &
                     mdl%elements(orig_iat)%species_id,orig_iat, &
                     ddec_dipoles(iat,1),ddec_dipoles(iat,2),&
                     ddec_dipoles(iat,3), &
                     ddec_quadrupoles(iat,1),ddec_quadrupoles(iat,2), &
                     ddec_quadrupoles(iat,3),ddec_quadrupoles(iat,4), &
                     ddec_quadrupoles(iat,5)
             end do
             write(stdout,'(a)') '------------------------------------------&
                  &-------------------------------------------------'
             write(stdout,'(a)') ''
          end if


          ! lpl: c.f. COMMENT#07
          ! lpl: Print moments
          if ( pub_ddec_moment > 0 ) then
             write(stdout,'(a)') '------------------------------------------&
                  &-------------------'
             write(stdout,'(a,I2,a,f4.2,a)') &
                  '           DDEC Radial Moment Order ', &
                                ! lpl: c.f. COMMENT#07
                  pub_ddec_moment, ' (X=',pub_ddec_IH_frac,')'
             !                  pub_ddec_moment_order,' (X=',pub_ddec_IH_frac,')'
             write(stdout,'(a)') '------------------------------------------&
                  &-------------------'
             write(stdout,'(a)') '  Atom        AIM-Valence     AIM-Core&
                  &        Reference'
             write(stdout,'(a)') '------------------------------------------&
                  &-------------------'

             do orig_iat=1,par%nat
                iat = par%distr_atom(orig_iat)
                write(stdout,'(1x,a4,I5,2x,3(2x,f14.7))') &
                     mdl%elements(orig_iat)%species_id,orig_iat, &
                     ddec_moments(iat,1), ddec_moments(iat,3), &
                     ddec_moments(iat,2)
             end do

             write(stdout,'(a)') '------------------------------------------&
                  &-------------------'
             write(stdout,'(a)') ''

          end if

          write(stdout,'(a)') '=================================&
               &==================================='


          write(stdout,'(a)') '------------------------------------------&
               &-------------------'
          write(stdout,'(a)') ''

          write(stdout,'(a)')  '                   Extra Point Charges '

          write(stdout,'(a)') '------------------------------------------&
               &-------------------'

          ! aeaa: Write number of charges to xyz file
          output_unit = utils_unit()

          open(output_unit, file='xyz_with_extra_point_charges.xyz')
          write(output_unit, '(I5)') ( par%nat + sum(abs(no_charges)) )
          write(output_unit, * ) ' '

          ! aeaa: Write anisotropy results to onetep file
          do  orig_iat=1,par%nat
             iat = par%distr_atom(orig_iat)

             write(stdout, '(a)' ) '           Atom Type and Number'
             write(stdout,'(16x,a4,I5)' ) mdl%elements(orig_iat)%species_id, &
                  orig_iat

             write(stdout, '(a)') '            Error in Potential '
             write(stdout, '(13x, f14.7)') error_aniso_before(iat)

             write(stdout, '(a)') '             VdW Radius Used '
             write(stdout, '(13x, f14.7)') vdw_radius(iat)

             ! aeaa: Calculate the overall molecules dipole
             overall_dipole_noextra(1) = overall_dipole_noextra(1) + &
                  ( mdl%elements(orig_iat)%centre%X) &
                  * ddec_charge(iat)
             overall_dipole_noextra(2) = overall_dipole_noextra(2) + &
                  ( mdl%elements(orig_iat)%centre%Y) &
                  * ddec_charge(iat)
             overall_dipole_noextra(3) = overall_dipole_noextra(3) + &
                  ( mdl%elements(orig_iat)%centre%Z) &
                  * ddec_charge(iat)

             ! aeaa: Write charges to file
             if ( no_charges(iat) /= 0 ) then

                write(stdout, '(a)') '               Positions'

                do i = 1, ( abs(no_charges(iat)) + 1 )
                   write(stdout, '(3(f14.7))') optimum_posit(1,i,iat), &
                        optimum_posit(2,i,iat),  optimum_posit(3,i,iat)
                end do

                write(stdout, '(a)') '    Off-Center Point Charge Charges'

                do i = 1,( abs(no_charges(iat)) + 1 )
                   write(stdout, '(13x, f14.7)') optimum_charges(i,iat)
                   ! aeaa: Molecule dipole computed below in e.bohr
                   overall_dipole_withpc(1) = overall_dipole_withpc(1) + &
                        ( optimum_posit(1,i,iat) / 0.529177_DP) &
                        * optimum_charges(i,iat)
                   overall_dipole_withpc(2) = overall_dipole_withpc(2) + &
                        ( optimum_posit(2,i,iat) / 0.529177_DP) &
                        * optimum_charges(i,iat)
                   overall_dipole_withpc(3) = overall_dipole_withpc(3) + &
                        ( optimum_posit(3,i,iat) / 0.529177_DP) &
                        * optimum_charges(i,iat)
                end do

                write(stdout,  '(a)') '       Error in Potential with Charges'
                write(stdout, '(13x, f14.7)') error_aniso(iat)

             else
                overall_dipole_withpc(1) = overall_dipole_withpc(1) + &
                     ( mdl%elements(orig_iat)%centre%X) &
                     * ddec_charge(iat)
                overall_dipole_withpc(2) = overall_dipole_withpc(2) + &
                     ( mdl%elements(orig_iat)%centre%Y) &
                     * ddec_charge(iat)
                overall_dipole_withpc(3) = overall_dipole_withpc(3) + &
                     ( mdl%elements(orig_iat)%centre%Z) &
                     * ddec_charge(iat)
             end if


             write(stdout, '(a)') '---------------------------------------&
                  &---------'

             ! aeaa: Write charges to new xyz file
             if ( no_charges(iat) == 0 ) then
                write(output_unit, '(A2, 4(f14.7))' ) &
                     mdl%elements(orig_iat)%species_id, &
                     ( mdl%elements(orig_iat)%centre%X * 0.529177_DP ), &
                     ( mdl%elements(orig_iat)%centre%Y * 0.529177_DP ), &
                     ( mdl%elements(orig_iat)%centre%Z * 0.529177_DP ), &
                     ddec_charge(iat)
             end if

             if (no_charges(iat) /= 0) then
                ! aeaa: Write center charges then other chagres
                write(output_unit, '(A2, 4(f14.7))' ) &
                     mdl%elements(orig_iat)%species_id, &
                     ( mdl%elements(orig_iat)%centre%X * 0.529177_DP ), &
                     ( mdl%elements(orig_iat)%centre%Y * 0.529177_DP ), &
                     ( mdl%elements(orig_iat)%centre%Z * 0.529177_DP ), &
                     optimum_charges(1, iat)
                do i = 2, ( abs(no_charges(iat)) + 1 )
                   write(output_unit, '(A1, 4(f14.7))' ) 'X', &
                        optimum_posit(1,i,iat), &
                        optimum_posit(2,i,iat),  optimum_posit(3,i,iat) , &
                        optimum_charges(i, iat)
                end do
             end if

          end do

          close(output_unit)

          ! aeaa: Write dipole moment of the molecule to file

          write(stdout,'(a)') ''

          write(stdout, '(a)') '        AIM Dipole Before (e.bohr)'
          write(stdout, '(a, f24.18)') '     dx  =    ', &
               overall_dipole_noextra(1)
          write(stdout, '(a, f24.18)') '     dy  =    ', &
               overall_dipole_noextra(2)
          write(stdout, '(a, f24.18)') '     dz  =    ', &
               overall_dipole_noextra(3)
          write(stdout, '(a, f24.18)') '    magnitude = ', &
               norm2(overall_dipole_noextra)

          write(stdout,'(a)') ''

          write(stdout, '(a)') '        AIM Dipole After (e.bohr)'
          write(stdout, '(a, f24.18)') '     dx  =    ', &
               overall_dipole_withpc(1)
          write(stdout, '(a, f24.18)') '     dy  =    ', &
               overall_dipole_withpc(2)
          write(stdout, '(a, f24.18)') '     dz  =    ', &
               overall_dipole_withpc(3)
          write(stdout, '(a, f24.18)') '    magnitude = ', &
               norm2(overall_dipole_withpc)


           write(stdout,'(a)') ''

       end if

       if ( allocated(ddec_moments) ) then
          deallocate(ddec_moments,stat=ierr)
          call utils_dealloc_check('ddec_main','ddec_moments',ierr)
       end if

       ! lpl: Deallocate multiple vector arrays
       if ( allocated(ddec_quadrupoles) ) then
          deallocate(ddec_quadrupoles,stat=ierr)
          call utils_dealloc_check('ddec_main','ddec_quadrupoles',ierr)
       end if
       if ( allocated(ddec_dipoles) ) then
          deallocate(ddec_dipoles,stat=ierr)
          call utils_dealloc_check('ddec_main','ddec_dipoles',ierr)
       end if

       if( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: ddec_main: Leaving properties #2'

       !aeaa - Off center charges array deallocation
       deallocate( optimum_posit, stat=ierr )
       call utils_dealloc_check('ddec_main','optimum_posit', ierr)
       deallocate( optimum_charges, stat=ierr )
       call utils_dealloc_check('ddec_main','optimum_charges', ierr)
       deallocate( error_aniso, stat=ierr )
       call utils_dealloc_check('ddec_main','error_aniso', ierr)
       deallocate( error_aniso_before, stat=ierr )
       call utils_dealloc_check('ddec_main','error_aniso_before', ierr)
       deallocate( no_charges,  stat=ierr )
       call utils_dealloc_check('ddec_main','no_charges', ierr)
       deallocate( atom_ctr_point, stat=ierr )
       call utils_dealloc_check('ddec_main','atom_ctr_point', ierr)
       deallocate( atom_ctr_point_array, stat=ierr )
       call utils_dealloc_check('ddec_main','atom_ctr_point_array', ierr)
       deallocate( vdw_radius, stat=ierr )
       call utils_dealloc_check('ddec_main','vdw_radius', ierr)
       deallocate( aim_volume, stat=ierr )
       call utils_dealloc_check('ddec_main','aim_volume', ierr)
       deallocate( overall_dipole_withpc, stat=ierr )
       call utils_dealloc_check('ddec_main',' overall_dipole_withpc', ierr)
       deallocate( overall_dipole_noextra, stat=ierr )
       call utils_dealloc_check('ddec_main','overall_dipole_noextra', ierr)
       deallocate( dist_atoms, stat=ierr )
       call utils_dealloc_check('ddec_main','dist_atoms', ierr)
       deallocate( all_element_types, stat=ierr )
       call utils_dealloc_check('ddec_main','all_element_types', ierr)

    end if ! END if ( pub_ddec_multipole .or. pub_ddec_moment )
    call comms_barrier

  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lpl: Interpolate, calculate <r>, and print radial densities                !
  !----------------------------------------------------------------------------!
  ! RAD_DENS DATA WILL BE DESTROYED HERE SO MAKE SURE THEY AREN'T NEEDED       !
  !----------------------------------------------------------------------------!

    if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: ddec_main: Entering output etc.'

    ! lpl: If interpolation is required
    if ( pub_ddec_interp_rad_dens ) then
       ! lpl: Deposit radial densities to local_box_fine, interpolate, and
       !      transfer back to radial grid

       do lc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)

          iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

          ! lpl: Deposit only the radial density into local_box_fine (local
          !      box is zeroed in 'internal_rad_dens_to_box') i.e. do not
          !      deposit the reference or core densities. We only need the
          !      density in 'local_box_fine%total_pseudoval'
          call internal_rad_dens_to_box(local_box_fine%temp_val, local_box_fine, &
               rad_dens(iat)%total_part_dens, rad_dens(iat))

          ! lpl: Interpolate densities and store into old_part_dens. Original
          !      converged density will remain in part_dens
          call internal_interp_rad_dens(rad_dens(iat)%old_part_dens, &
               rad_dens(iat)%total_part_dens, &
               rad_dens(iat), ref_dens, local_box_fine%temp_val, local_box_fine)

          ! lpl: If this doesn't work for some systems then pass to any
          !      'rad_dens%old_part_dens' on root proc
          if ( pub_ddec_write ) then
             output_filename = ''
             write(output_filename,'(a,I5.5)') &
                  'interp_dens_atom',par%orig_atom(iat)
             call ddec_write_dens(rad_dens(iat)%old_part_dens, &
                  rad_dens(iat)%rad, rad_dens(iat)%npts, &
                  trim(adjustl(output_filename)),pub_my_proc_id)
          end if

       end do
       call comms_barrier

    end if

    ! lpl: Print radial densities (each proc writes its own file)
    ! lpl: If printing from non-root is forbidden then copy new_part_dens
    !      to total_part_dens or some temp buffer before writing
    if(pub_ddec_write) then
       ! lpl: Print converged DDEC AIM and core densities
       do lc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
          iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

          output_filename = ''
          write(output_filename,'(a,I5.5)') 'rad_dens_atom',par%orig_atom(iat)
          ! lpl: MOD#12 - Print n^{ISA}_A (new_part_dens), not w_A (total_part_dens)
          call ddec_write_dens(rad_dens(iat)%new_part_dens, &
               rad_dens(iat)%rad, rad_dens(iat)%npts, &
               trim(adjustl(output_filename)),pub_my_proc_id)

          if ( include_coredens ) then
             output_filename = ''
             write(output_filename,'(a,I5.5)') &
                  'core_dens_atom_',par%orig_atom(iat)
             call ddec_write_dens(rad_dens(iat)%core_dens, &
                  rad_dens(iat)%rad, rad_dens(iat)%npts, &
                  trim(adjustl(output_filename)))
          end if
       end do

       ! lpl: Print reference IH densities used in this calculation
       do isp=1,par%num_species
          do iion=ref_dens(isp)%max_anion,ref_dens(isp)%max_cation
             ! lpl: MOD#07 - changed print condition to 'ref_dens%init_status'
             if ( ref_dens(isp)%init_status(iion) ) then
                output_filename = ''
                write(output_filename,'(A,A,A1,sp,I4.3)') 'ref_dens_', &
                      trim(adjustl(ref_dens(isp)%element%species_id)),'_',iion
                call ddec_write_dens(ref_dens(isp)%&
                     ion_dens(iion)%dens,ref_dens(isp)%rad, &
                     ref_dens(isp)%npts, trim(adjustl(output_filename)))
             end if
          end do

          ! lpl: Write out core density if available
          if ( associated(ref_dens(isp)%core_dens) ) then
             output_filename = ''
             write(output_filename,'(A,A)') 'ref_core_', &
                   trim(adjustl(ref_dens(isp)%element%species_id))
             call ddec_write_dens(ref_dens(isp)%core_dens, &
                  ref_dens(isp)%rad,ref_dens(isp)%npts, &
                  trim(adjustl(output_filename)))
          end if
       end do
    end if
    call comms_barrier

    if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: ddec_main: Leaving output etc.'

  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lpl: Finalize subroutine and deallocate arrays                             !
  !----------------------------------------------------------------------------!

    if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: ddec_main: Entering finalization'

    ! lpl: Subtracat added core density to recover original (valence) density
    if ( include_coredens ) then
       total_dens_fine = total_dens_fine - core_dens_fine
       call comms_barrier
    end if

    deallocate(deposit_extract_buffer,stat=ierr)
    call utils_dealloc_check('ddec_main','deposit_extract_buffer',ierr)

    call internal_rad_dens_destroy(rad_dens)
    call internal_ref_dens_destroy(ref_dens)
    call internal_local_box_destroy(local_box_fine)

    deallocate(atom_popn_prev,stat=ierr)
    call utils_dealloc_check('ddec_main','atom_popn_prev',ierr)

    ! lpl: MOD#05
    deallocate(ddec_core_popn,stat=ierr)
    call utils_dealloc_check('ddec_main','ddec_core_popn',ierr)

    deallocate(val_popn,stat=ierr)
    call utils_dealloc_check('ddec_main','val_popn',ierr)
    deallocate(core_popn,stat=ierr)
    call utils_dealloc_check('ddec_main','core_popn',ierr)
    deallocate(core_diff,stat=ierr)
    call utils_dealloc_check('ddec_main','core_diff',ierr)

    deallocate(atom_popn,stat=ierr)
    call utils_dealloc_check('ddec_main','atom_popn',ierr)

    nullify(local_box_fine)

    nullify(ref_dens)
    nullify(rad_dens)

    if(pub_on_root) &
       write(stdout,'(a)') '================================================'

    if (pub_debug) then
       close(unit=debug_unit,iostat=ierr)
       call utils_close_unit_check('ddec_main','DDEC_DEBUG',ierr)

       if( pub_on_root ) then
          write(stdout,'(a)') ' DEBUG: ddec_main: Leaving finalization'
          write(stdout,'(a)') ' DEBUG: Exiting ddec_main'
       end if
    end if

    ! lpl: Stop timer
    call timer_clock('ddec_main',2)

    !----------------------------------------------------------------------------!

   contains

   !===========================================================================!

     subroutine internal_aim_density(aim_val, total_val, &
          pseudo_val, local_box, part_dens, rad_dens)

       !-----------------------------------------------------------------------!
       ! lpl: Calculates AIM density                                           !
       !    : Used sparingly in several parts of the code to regenerate the    !
       !      AIM density from the stored spherically-averaged AIM             !
       !    : Probably not used in the valence iterator due to the complexity  !
       !      of having to compute many other quantities                       !
       !-----------------------------------------------------------------------!
       ! Output : 'aim_val' containing the new AIM density n_A^{new}(r)        !
       ! Input  : 'total_val' containing the total (i.e. real) density n(r)    !
       !        : 'pseudo_val' containing the pseudodensity \sum_B n_B(r)      !
       !        : 'part_dens' containing the old spherically-averaged AIM      !
       !           density n_A^{old}(r)                                        !
       ! Info   : 'local_box', 'rad_dens' - info. Required by subroutine       !
       !          called in this function                                      !
       !-----------------------------------------------------------------------!

       implicit none

       ! Output
       real(kind=DP), intent(out) :: aim_val(:,:,:)
       ! Input
       real(kind=DP), intent(in) :: total_val(:,:,:), pseudo_val(:,:,:)
       real(kind=DP), intent(in) :: part_dens(:)
       ! Info Input
       type(DDEC_LOCAL_BOX), intent(in) :: local_box
       type(DDEC_RAD_DENS) ,intent(in) :: rad_dens

       ! Local variables
       type(POINT) :: local_atom_ctr, local_grid_coord
       integer :: iat, it1, it2, it3, ir

       if ( size(aim_val,1) /= local_box%n1 .or. &
            size(aim_val,2) /= local_box%n2 .or. &
            size(aim_val,3) /= local_box%n3 ) &
            call utils_abort(' ERROR in internal_aim_density:&
                 & size(aim_val) /= local_box')

       if ( size(part_dens) /= rad_dens%npts ) &
            call utils_abort(' ERROR in internal_aim_density:&
                 & size(part_dens) /= rad_dens%npts')

       iat = rad_dens%iat

       local_atom_ctr = mdl%elements(par%orig_atom(iat))%centre &
            - local_box%box_coord_origin(iat)

       do it1=1,local_box%n1
          do it2=1,local_box%n2
             do it3=1,local_box%n3

                ! lpl: Coord w.r.t. voxel 1,1,1 centre
                local_grid_coord &
                     = real(it1 - 1,kind=DP)*grid%da1 &
                     + real(it2 - 1,kind=DP)*grid%da2 &
                     + real(it3 - 1,kind=DP)*grid%da3

                ! lpl: Atom-centred radius index
                ir = aint(geometry_magnitude(local_atom_ctr &
                     - local_grid_coord)/rad_dens%dr) + 1

                if ( ir < rad_dens%npts .and. pseudo_val(it1,it2,it3) &
                     > pub_ddec_zero_thresh ) then
                   if (pub_debug) then
                   if ( pseudo_val(it1,it2,it3) < part_dens(ir) ) &
                        call utils_abort(' ERROR in internal_aim_density:&
                             & pseudo_val(it1,it2,it3) < part_dens(ir)')
                   end if
                   ! lpl: Compute AIM density
                   ! n_A^{new}(r) = n_A^{old}(r)*[n(r)/{sum_B n_B^{old}(r)}]
                   aim_val(it1,it2,it3) = part_dens(ir)*&
                        (total_val(it1,it2,it3)/pseudo_val(it1,it2,it3))
                else
                   aim_val(it1,it2,it3) = 0.0_DP
                end if

             end do
          end do
       end do

     end subroutine internal_aim_density

   !===========================================================================!

     subroutine internal_core_iterator(rad_dens, local_box, core_popn)

       !-----------------------------------------------------------------------!
       ! lpl: ISA iterator                                                     !
       !    : Procedure to genreate AIM density for each ISA iteration         !
       !    : Primarily used in core iterator                                  !
       !-----------------------------------------------------------------------!
       ! InOut  : 'local_box_fine', 'rad_dens' - makes a mess of arrays stored !
       !          in both objects, so beware. Final output is stored in        !
       ! Output : 'core_popn' - AIM atomic population                          !
       !-----------------------------------------------------------------------!

       type(DDEC_LOCAL_BOX), intent(inout) :: local_box
       type(DDEC_RAD_DENS), intent(inout) :: rad_dens

       real(kind=DP), intent(out) :: core_popn
       type(POINT) :: local_atom_ctr, local_grid_coord
       integer :: iat, it1, it2, it3, ir

       iat = rad_dens%iat

       ! lpl: Reset radial arrays
       rad_dens%new_part_dens = 0.0_DP
       rad_dens%rad_npts = 0

       core_popn = 0.0_DP

       local_atom_ctr = mdl%elements(par%orig_atom(iat))%centre &
            - local_box%box_coord_origin(iat)

       do it1=1,local_box%n1
          do it2=1,local_box%n2
             do it3=1,local_box%n3

                ! lpl: Coord w.r.t. voxel 1,1,1 centre
                local_grid_coord &
                     = real(it1 - 1,kind=DP)*grid%da1 &
                     + real(it2 - 1,kind=DP)*grid%da2 &
                     + real(it3 - 1,kind=DP)*grid%da3

                ! lpl: Atom-centred radius index
                ir = aint(geometry_magnitude(local_atom_ctr &
                     - local_grid_coord)/rad_dens%dr) + 1

                if ( ir > rad_dens%npts .or. &
                     local_box%core_pseudoval(it1,it2,it3) &
                     < pub_ddec_zero_thresh ) cycle

                if (pub_debug) then
                if ( local_box%core_pseudoval(it1,it2,it3) &
                     < rad_dens%core_dens(ir) ) &
                     call utils_abort(' ERROR in internal_core_iterator:&
                          & local_box%core_pseudoval(it1,it2,it3) <&
                          & rad_dens%core_dens(ir)')
                end if

                ! lpl: Calculate new AIM weighting
                !      n_A^{new}(r) = n_A^{old}(r)*[n(r)/{sum_B n_B^{old}(r)}]
                rad_dens%new_part_dens(ir) = rad_dens%new_part_dens(ir) &
                     + rad_dens%core_dens(ir)*(local_box%core_val(it1,it2,it3)/&
                     local_box%core_pseudoval(it1,it2,it3))

                rad_dens%rad_npts(ir) = rad_dens%rad_npts(ir) + 1

             end do
          end do
       end do

       ! lpl: Spherically average AIM density from accumulated sum
       do ir=1,rad_dens%npts
          if ( rad_dens%rad_npts(ir) > 0 ) then
             core_popn = core_popn + rad_dens%new_part_dens(ir)
             rad_dens%new_part_dens(ir) = &
                  rad_dens%new_part_dens(ir)/rad_dens%rad_npts(ir)
          else
             rad_dens%new_part_dens(ir) = 0.0_DP
          end if
       end do

     end subroutine internal_core_iterator

     subroutine internal_core_corr_iterator(local_box, rad_dens, core_diff)

       ! lpl: Computes 'corrected_core_density' in CHARGEMOL and replaces
       !      data in 'local_val'

       ! InOut
       type(DDEC_LOCAL_BOX), intent(inout) :: local_box

       ! Input
      type(DDEC_RAD_DENS), intent(in) :: rad_dens
       real(kind=DP), intent(in) :: core_diff


       ! Local variables
       type(POINT) :: local_atom_ctr, local_grid_coord
       real(kind=DP) :: corr_fac, sum_core_aim_cubed, max_core_aim
       integer :: iat, it1, it2, it3, ir

       ! lpl: Compute AIM density for each atom
       call internal_aim_density(local_box_fine%aim_val, &
            local_box_fine%core_val, local_box_fine%core_pseudoval, &
            local_box_fine, rad_dens%core_dens, rad_dens)

       ! lpl: Compute correction factor 'corr_fac' (called 'K_factor' in
       !      CHARGEMOL)
       sum_core_aim_cubed = sum(local_box_fine%aim_val**3)
       max_core_aim = maxval(local_box_fine%aim_val)

       if ( sum_core_aim_cubed > pub_ddec_zero_thresh**3 .and. &
            max_core_aim > pub_ddec_zero_thresh ) &
            corr_fac = min(core_diff/(sum_core_aim_cubed*vox_vol), &
                 0.25_DP/(max_core_aim**2))

       iat = rad_dens%iat

       local_atom_ctr = mdl%elements(par%orig_atom(iat))%centre &
            - local_box%box_coord_origin(iat)

       do it1=1,local_box%n1
          do it2=1,local_box%n2
             do it3=1,local_box%n3

                ! lpl: Coord w.r.t. voxel 1,1,1 centre
                local_grid_coord &
                     = real(it1 - 1,kind=DP)*grid%da1 &
                     + real(it2 - 1,kind=DP)*grid%da2 &
                     + real(it3 - 1,kind=DP)*grid%da3

                ! lpl: Atom-centred radius index
                ir = aint(geometry_magnitude(local_atom_ctr &
                     - local_grid_coord)/rad_dens%dr) + 1

                if ( ir < rad_dens%npts .or. local_box%aim_val(it1,it2,it3) &
                     > sqr_pub_ddec_zero_thresh ) then
                   local_box%aim_val(it1,it2,it3) = &
                        local_box%aim_val(it1,it2,it3)/sqrt(1.0_DP - &
                        2.0_DP*corr_fac*local_box%aim_val(it1,it2,it3)**2)
                else
                   local_box%aim_val(it1,it2,it3) = 0.0_DP
                end if

             end do
          end do
       end do

     end subroutine internal_core_corr_iterator

   !===========================================================================!

     subroutine internal_valence_iterator(rad_dens, local_box, &
          atom_popn, core_diff)

       !-----------------------------------------------------------------------!
       ! lpl: Valence DDEC iterator                                            !
       !-----------------------------------------------------------------------!

       type(DDEC_LOCAL_BOX), intent(inout) :: local_box
       type(DDEC_RAD_DENS), intent(inout) :: rad_dens

       real(kind=DP), intent(out) :: atom_popn
       real(kind=DP), optional, intent(out) :: core_diff

       type(POINT) :: local_atom_ctr, local_grid_coord
       real(kind=DP) :: temp_val
       integer :: iat, it1, it2, it3, ir

       ! lpl: The entire procedure for valence AIM partitioning per iteration
       !      per atom per proc

       iat = rad_dens%iat

       ! lpl: Reset radial arrays
       rad_dens%new_part_dens = 0.0_DP
       rad_dens%rad_npts = 0
       rad_dens%renorm_slope = 0.0_DP

       if ( include_refdens ) then
          rad_dens%ih_ratio  = 0.0_DP
          rad_dens%Yavg_sqrt = 0.0_DP
          rad_dens%Yavg_ratio = 0.0_DP
       end if

       atom_popn = 0.0_DP
       if ( present(core_diff) ) core_diff = 0.0_DP

       local_atom_ctr = mdl%elements(par%orig_atom(iat))%centre &
            - local_box%box_coord_origin(iat)

       do it1=1,local_box%n1
          do it2=1,local_box%n2
             do it3=1,local_box%n3

                ! lpl: Coord w.r.t. voxel 1,1,1 centre
                local_grid_coord &
                     = real(it1 - 1,kind=DP)*grid%da1 &
                     + real(it2 - 1,kind=DP)*grid%da2 &
                     + real(it3 - 1,kind=DP)*grid%da3

                ! lpl: Atom-centred radius index
                ir = aint(geometry_magnitude(local_atom_ctr &
                     - local_grid_coord)/rad_dens%dr) + 1

                if ( ir > rad_dens%npts ) cycle

                if ( local_box%total_pseudoval(it1,it2,it3) &
                     < pub_ddec_zero_thresh ) then

                if (pub_debug) then
                if ( local_box%total_pseudoval(it1,it2,it3) &
                     < rad_dens%total_part_dens(ir) ) &
                     call utils_abort(' ERROR in internal_valence_iterator:&
                          & local_box%total_pseudoval(it1,it2,it3) <&
                          & rad_dens%total_part_dens(ir)')
                end if

                   ! lpl: Accumulate core density that might be suppressed
                   !      due to condition above
                   if ( include_coredens ) then
                      if ( .not. present(core_diff) ) call utils_abort( &
                           ' ERROR in internal_valence_iterator: core_diff&
                           & not present')
                      if ( local_box%core_pseudoval(it1,it2,it3) &
                           > pub_ddec_zero_thresh ) then
                         core_diff = core_diff + rad_dens%core_dens(ir)*&
                              (local_box%core_val(it1,it2,it3)/&
                              local_box%core_pseudoval(it1,it2,it3))
                      end if
                   end if

                   cycle

                end if

       !--------------------------------------------------------------------!
       ! lpl: ISA weighting                                                 !
       !--------------------------------------------------------------------!

                ! lpl: Calculate new AIM weighting
                !      n_A^{new}(r) = n_A^{old}(r)*[n(r)/{sum_B n_B^{old}(r)}]
                temp_val = rad_dens%total_part_dens(ir)*&
                     (local_box%total_val(it1,it2,it3)/&
                     local_box%total_pseudoval(it1,it2,it3))
                rad_dens%new_part_dens(ir) = rad_dens%new_part_dens(ir) &
                     + temp_val

                ! lpl: Calculate renorm_slope (Eq. 73 in JCTC 2012)
                !      n_A^{new}(r)*( 1 - n_A^{old}(r)/{sum_B n_B^{old}(r)} )
                rad_dens%renorm_slope = rad_dens%renorm_slope + &
                     temp_val*(1.0_DP - rad_dens%total_part_dens(ir)/&
                     local_box%total_pseudoval(it1,it2,it3))

                rad_dens%rad_npts(ir) = rad_dens%rad_npts(ir) + 1

       !--------------------------------------------------------------------!

       !--------------------------------------------------------------------!
       ! lpl: IH weighting                                                  !
       !--------------------------------------------------------------------!

                if ( include_refdens ) then

                   ! lpl: IH ratio (PtoWref in CHARGEMOL)
                   if ( local_box%ref_pseudoval(it1,it2,it3) &
                        > sqr_pub_ddec_zero_thresh ) then
                      rad_dens%ih_ratio(ir) = rad_dens%ih_ratio(ir) + &
                           local_box%total_val(it1,it2,it3)/&
                           local_box%ref_pseudoval(it1,it2,it3)
                   else
                      rad_dens%ih_ratio(ir) = rad_dens%ih_ratio(ir) + 1.0_DP
                   end if

                   ! lpl: Yavg_sqrt (sqrt_Wref in CHARGEMOL)
                   temp_val = sqrt(local_box%Yavg_pseudoval(it1,it2,it3))
                   rad_dens%Yavg_sqrt(ir) = rad_dens%Yavg_sqrt(ir) &
                        + temp_val

                   ! lpl: Yavg_ratio (wA_over_sqrt_Wref in CHARGEMOL)
                   if ( local_box_fine%Yavg_pseudoval(it1,it2,it3) &
                        > sqr_pub_ddec_zero_thresh ) then
                      rad_dens%Yavg_ratio(ir) = rad_dens%Yavg_ratio(ir) + &
                           rad_dens%Yavg_dens(ir)/temp_val
                   end if

                end if

       !--------------------------------------------------------------------!

             end do
          end do
       end do

       rad_dens%renorm_slope = vox_vol*rad_dens%renorm_slope

       ! lpl: Spherical averaging
       do ir=1,rad_dens%npts
          if ( rad_dens%rad_npts(ir) > 0 ) then
             atom_popn = atom_popn + rad_dens%new_part_dens(ir)
             rad_dens%new_part_dens(ir) = &
                  rad_dens%new_part_dens(ir)/rad_dens%rad_npts(ir)
          else
             rad_dens%new_part_dens(ir) = 0.0_DP
          end if
       end do

       if ( include_refdens ) then
          do ir=1,rad_dens%npts
             if ( rad_dens%rad_npts(ir) > 0 ) then
                rad_dens%ih_ratio(ir) = &
                     rad_dens%ih_ratio(ir)/rad_dens%rad_npts(ir)
                rad_dens%Yavg_sqrt(ir) = &
                     rad_dens%Yavg_sqrt(ir)/rad_dens%rad_npts(ir)
                rad_dens%Yavg_ratio(ir) = &
                     rad_dens%Yavg_ratio(ir)/rad_dens%rad_npts(ir)
             else
                rad_dens%ih_ratio(ir) = 1.0_DP
                rad_dens%Yavg_sqrt(ir) = sqrt_pub_ddec_zero_thresh
                rad_dens%Yavg_ratio(ir) = sqrt_pub_ddec_zero_thresh
             end if
          end do
       end if

     end subroutine internal_valence_iterator

   !===========================================================================!

     subroutine internal_deposit_box_all(global_dens_fine, rad_dens, &
          rad_dens_flag, local_box, local_box_flag, deposit_extract_buffer)

       !-----------------------------------------------------------------------!
       ! lpl: Deposits radial densities in 'rad_dens' into the global density  !
       !      array 'global_dens_fine'                                         !
       !    : Wrapper to convert all radial densities to Cartesian values      !
       !      followed by deposition into 'global_dens_fine'                   !
       !-----------------------------------------------------------------------!
       ! InOut  : 'global_dens_fine' - Global density array                    !
       ! Input  : 'rad_dens', 'rad_dens_flag' - 'rad_dens' array of all atoms. !
       !          'rad_dens_flag' specifies which radial density in 'rad_dens' !
       !          should be deposited into 'global_dens_fine'.                 !
       !        : 'local_box_fine' - required for info on dimensions.          !
       ! Buffer : 'deposit_extract_buffer' - used as a buffer in a subroutine. !
       !-----------------------------------------------------------------------!

       implicit none

       ! Output
       real(kind=DP), intent(inout) :: global_dens_fine(grid%ld1, &
         grid%ld2,grid%max_slabs12)

       ! InOut
       type(DDEC_LOCAL_BOX), intent(inout) :: local_box
       character(len=*), intent(in) :: local_box_flag

       ! Input
       type(DDEC_RAD_DENS), intent(in) :: rad_dens(par%nat)
       character(len=*), intent(in) :: rad_dens_flag

       ! Buffer
       real(kind=DP), intent(out) :: deposit_extract_buffer(:,:,:)

       ! Local variables
       integer :: iat, lc_iat
       logical :: i_have_box
       real(kind=DP), pointer :: rad_dens_ptr(:)
       real(kind=DP), pointer :: local_box_ptr(:,:,:)

       if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Entering internal_deposit_box_all'
       if (pub_debug) call timer_clock('internal_deposit_boxes_all',1)

       ! lpl: Deposit radial densities specified by 'rad_dens' and
       !      'rad_dens_flag' from all atoms into 'global_dens_fine'
       do lc_iat=1,par%max_atoms_on_proc

          ! lpl: Whether we need to extract density box for this atom
          if (lc_iat <= par%num_atoms_on_proc(pub_my_proc_id)) then
             i_have_box = .true.
             iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

             ! lpl: Set 'rad_dens_ptr' to point to the correct radial density
             !      that needs to be deposited into the global density
             call internal_set_rad_dens_ptr(rad_dens_ptr, rad_dens(iat), &
                  rad_dens_flag)
             ! lpl: Set 'local_box_ptr' to point to a buffer density into
             !      which 'rad_dens_ptr' can be deposited
             call internal_set_local_box_ptr(local_box_ptr, local_box, &
                  local_box_flag)

             ! lpl: Get radial densities from 'rad_dens_ptr' into
             !      'local_box_ptr'
             call internal_rad_dens_to_box(local_box_ptr, local_box, &
                  rad_dens_ptr, rad_dens(iat))
          else
             i_have_box = .false.
             iat = par%first_atom_on_proc(pub_my_proc_id) ! Prevents segfault
             call internal_set_local_box_ptr(local_box_ptr, local_box, &
                  local_box_flag) ! Prevents segfault
          end if

          ! lpl: Deposit radial density now in 'local_box_ptr' into
          !      global density 'global_dens_fine'
          call internal_deposit_box(global_dens_fine, local_box_ptr, &
               local_box, iat, i_have_box, deposit_extract_buffer)
       end do

       call comms_barrier


       if (pub_debug) call timer_clock('internal_deposit_boxes_all',2)
       if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_deposit_boxes_all'

     end subroutine internal_deposit_box_all

   !===========================================================================!

     subroutine internal_deposit_box(global_dens_fine, local_box_val, &
          local_box, iat, i_have_box, deposit_extract_buffer)

       !-----------------------------------------------------------------------!
       ! lpl: Deposits density from 'local_box_val' into 'global_dens_fine'.   !
       !      Dumb wrapper for 'cell_grid_deposit_box'.                        !
       !-----------------------------------------------------------------------!
       ! Output : 'global_dens_fine' - global density array                    !
       ! Input  : 'local_box_val', 'local_box' - local_box density to be       !
       !          deposited. 'local_box' contains grid info.                   !
       !        : 'iat' - the atom to which 'local_box_val' belongs.           !
       !        : 'i_have_box' - Used in cell_grid_deposit_box                 !
       ! Buffer : 'deposit_extract_buffer' - Buffer for cell_grid_deposit_box  !
       !-----------------------------------------------------------------------!

       implicit none

       ! Info Input
       type(DDEC_LOCAL_BOX), intent(in) :: local_box

       ! Output
       real(kind=DP), intent(inout) :: global_dens_fine(grid%ld1, &
            grid%ld2,grid%max_slabs12)

       ! Input
       real(kind=DP), intent(in) :: local_box_val(:,:,:)
       logical, intent(in) :: i_have_box
       integer, intent(in) :: iat

       ! Buffer
       real(kind=DP), intent(inout) :: deposit_extract_buffer(:,:,:)

       if ( size(local_box_val,1) /= local_box%n1 .or. &
            size(local_box_val,2) /= local_box%n2 .or. &
            size(local_box_val,3) /= local_box%n3 ) &
            call utils_abort(' ERROR in internal_deposit_box:&
                 & size(local_box_val) /= local_box')

       if ( size(deposit_extract_buffer,1) < local_box%n1 .or. &
            size(deposit_extract_buffer,2) < local_box%n2 .or. &
            size(deposit_extract_buffer,3) < grid%max_slabs12 ) &
            call utils_abort(' ERROR in internal_deposit_box:&
                 & size(deposit_extract_buffer) < local_box')

       call cell_grid_deposit_box(global_dens_fine, &
            local_box_val, deposit_extract_buffer, &
            grid, local_box%n1, local_box%n2, local_box%n3, &
            local_box%n1, local_box%n2, &
            local_box%box_index_origin(iat)%ni1, &
            local_box%box_index_origin(iat)%ni2, &
            local_box%box_index_origin(iat)%ni3, i_have_box, .false.)

     end subroutine internal_deposit_box

   !===========================================================================!

     subroutine internal_extract_box(local_box_val, local_box, iat, &
          global_dens_fine, i_need_box, deposit_extract_buffer)

       !-----------------------------------------------------------------------!
       ! lpl: Extracts density from 'global_dens_fine' into 'local_box_val'.   !
       !      Dumb wrapper for 'cell_grid_extract_box'.                        !
       !-----------------------------------------------------------------------!
       ! Output : 'local_box_val' - local_box density to be extracted to.      !
       ! Input  : 'local_box' - Contains grid info relevant to            !
       !          'local_box_val'.                                             !
       !        : 'iat' - the atom to which 'local_box_val' belongs.           !
       !        : 'global_dens_fine' - global density array                    !
       !        : 'i_need_box' - Used in cell_grid_extract_box                 !
       ! Buffer : 'deposit_extract_buffer' - Buffer for cell_grid_extract_box  !
       !-----------------------------------------------------------------------!

       implicit none

       ! Info Input
       type(DDEC_LOCAL_BOX), intent(in) :: local_box

       ! Output
       real(kind=DP), intent(out) :: local_box_val(:,:,:)

       ! Input
       real(kind=DP), intent(in) :: global_dens_fine(grid%ld1, &
            grid%ld2,grid%max_slabs12)
       logical, intent(in) :: i_need_box
       integer, intent(in) :: iat

       ! Buffer
       real(kind=DP), intent(inout) :: deposit_extract_buffer(:,:,:)

       if ( size(local_box_val,1) /= local_box%n1 .or. &
            size(local_box_val,2) /= local_box%n2 .or. &
            size(local_box_val,3) /= local_box%n3 ) &
            call utils_abort(' ERROR in internal_extract_box:&
                 & size(local_box_val) /= local_box')

       if ( size(deposit_extract_buffer,1) < local_box%n1 .or. &
            size(deposit_extract_buffer,2) < local_box%n2 .or. &
            size(deposit_extract_buffer,3) < grid%max_slabs12 ) &
            call utils_abort(' ERROR in internal_extract_box:&
                 & size(deposit_extract_buffer) < local_box')

       call cell_grid_extract_box(local_box_val, &
            deposit_extract_buffer, global_dens_fine, grid, &
            local_box%n1,local_box%n2,local_box%n3, &
            local_box%n1,local_box%n2, &
            local_box%box_index_origin(iat)%ni1, &
            local_box%box_index_origin(iat)%ni2, &
            local_box%box_index_origin(iat)%ni3,i_need_box,.false.)

     end subroutine internal_extract_box

   !===========================================================================!

     subroutine internal_set_local_box_ptr(local_box_ptr, local_box_fine, &
          local_box_flag)

       !-----------------------------------------------------------------------!
       ! lpl: Trivial subroutine to set a pointer 'local_box_ptr' to point to  !
       !      an array from the object 'local_box_fine' specified by           !
       !      'local_box_flag'                                                 !
       !-----------------------------------------------------------------------!
       ! Output : 'local_box_ptr' - pointer to array in 'local_box_fine'       !
       ! Input  : 'local_box_fine', 'local_box_flag' - Flag to specify array   !
       !          'local_box_fine'                                             !
       !-----------------------------------------------------------------------!

       implicit none

       ! Output
       real(kind=DP), pointer, intent(out) :: local_box_ptr(:,:,:)
       ! Input
       type(DDEC_LOCAL_BOX), target, intent(in) :: local_box_fine
       character(len=*), intent(in) :: local_box_flag

       select case (local_box_flag)
          case ('total_val')
             local_box_ptr => local_box_fine%total_val
          case ('total_pseudoval')
             local_box_ptr => local_box_fine%total_pseudoval
          case ('ref_pseudoval')
             local_box_ptr => local_box_fine%ref_pseudoval
          case ('Yavg_pseudoval')
             local_box_ptr => local_box_fine%Yavg_pseudoval
          case ('core_val')
             local_box_ptr => local_box_fine%core_val
          case ('core_pseudoval')
             local_box_ptr => local_box_fine%core_pseudoval
          case ('aim_val')
             local_box_ptr => local_box_fine%aim_val
          case default
             call utils_abort(' ERROR in internal_set_local_box_ptr: &
                  &invalid flag "'//local_box_flag//'"')
       end select

     end subroutine internal_set_local_box_ptr

   !===========================================================================!

     subroutine internal_set_rad_dens_ptr(rad_dens_ptr, rad_dens, rad_dens_flag)

       !-----------------------------------------------------------------------!
       ! lpl: Trivial subroutine to set a pointer 'rad_dens_ptr' to point to   !
       !      an array from the object 'rad_dens' specified by 'rad_dens_flag' !
       !-----------------------------------------------------------------------!
       ! Output : 'rad_dens_ptr' - pointer to array in 'rad_dens'              !
       ! Input  : 'rad_dens', 'rad_dens_flag' - Flag to specify array in       !
       !          'rad_dens'                                                   !
       !-----------------------------------------------------------------------!

       implicit none

       ! Output
       real(kind=DP), pointer, intent(out) :: rad_dens_ptr(:)

       ! Input
       type(DDEC_RAD_DENS), target, intent(in) :: rad_dens
       character(len=*), intent(in) :: rad_dens_flag

       select case (rad_dens_flag)
          case ('new_part_dens')
             rad_dens_ptr => rad_dens%new_part_dens
          case ('old_part_dens')
             rad_dens_ptr => rad_dens%old_part_dens
          case ('total_part_dens')
             rad_dens_ptr => rad_dens%total_part_dens
          case ('Yavg_dens')
             rad_dens_ptr => rad_dens%Yavg_dens
          case ('Yavg_sqrt')
             rad_dens_ptr => rad_dens%Yavg_sqrt
          case ('Yavg_ratio')
             rad_dens_ptr => rad_dens%Yavg_ratio
          case ('ih_ratio')
             rad_dens_ptr => rad_dens%ih_ratio
          case ('core_dens')
             rad_dens_ptr => rad_dens%core_dens
          case ('ref_dens')
             rad_dens_ptr => rad_dens%ref_dens
          case default
             call utils_abort(' ERROR in internal_set_rad_dens_ptr: &
                  &invalid flag "'//rad_dens_flag//'"')
       end select

     end subroutine internal_set_rad_dens_ptr

   !===========================================================================!

  !----------------------------------------------------------------------------!
  ! Radial density 'DDEC_RAD_DENS' initialization and destruction subroutines  !
  !   1. internal_rad_dens_init                                                !
  !   2. internal_rad_dens_destroy                                             !
  !----------------------------------------------------------------------------!

     subroutine internal_rad_dens_init(rad_dens,ref_dens)

       !-----------------------------------------------------------------------!
       ! lpl: Initializes radial densities DDEC_RAD_DENS for all atoms         !
       !-----------------------------------------------------------------------!
       ! InOut : rad_dens - Array of radial densities for ALL atoms            !
       ! Input : ref_dens - Array of reference densities for ALL species       !
       !-----------------------------------------------------------------------!

       implicit none

       ! lpl: Input
       type(DDEC_RAD_DENS), pointer, intent(inout) :: rad_dens(:)
       type(DDEC_REF_DENS), intent(in), target :: ref_dens(par%num_species)

       ! lpl: Local variables
       integer :: iat, it, isp

      if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Entering internal_rad_dens_init'

       ! lpl: Check if radial densities are initialized and do so if not
       if(associated(rad_dens)) call utils_abort(' ERROR: Attempting to &
            &initialize associated rad_dens in internal_rad_dens_init')
       allocate(rad_dens(par%nat),stat=ierr)
       call utils_alloc_check('internal_rad_dens_init','=>rad_dens',ierr)

       if(.not. associated(rad_dens)) call utils_abort(' ERROR: Failed to &
            &initialize rad_dens in internal_rad_dens_init')

       ! lpl: Assign radial density attributes and
       do iat=1,par%nat

          isp = mdl%elements(par%orig_atom(iat))%species_number

          ! lpl: Make sure ref_dens has been 'properly' initialized
          if ( .not. ref_dens(isp)%initialized ) call utils_abort(' ERROR: &
               &Attempted to access DDEC_REF_DENS which has not been &
               &properly initialized')

          ! lpl: Copy atomic number
          rad_dens(iat)%iat = iat

          ! lpl: Radial coordinates
          rad_dens(iat)%npts = ref_dens(isp)%npts
          rad_dens(iat)%rcut = ref_dens(isp)%rcut
          rad_dens(iat)%dr   = rad_dens(iat)%rcut/(1.0_DP*rad_dens(iat)%npts)

          ! lpl: Radial grid (used only to define where in a bin of width
          !      'dr' should the radius be taken - not used when e.g. binning
          !      Cartesian points into radial shells, so do not alter this to
          !      be something irregular (i.e. r_i - r_i-1 < 2*dr always)
          allocate(rad_dens(iat)%rad(rad_dens(iat)%npts),stat=ierr)
          call utils_alloc_check('internal_rad_dens_init', &
               'rad_dens(iat)%rad',ierr)
          rad_dens(iat)%rad = ref_dens(isp)%rad

          ! lpl: Number of points per spherical shell for each atom (stored
          !      to be used in routines related to density reshaping)
          allocate(rad_dens(iat)%rad_npts(rad_dens(iat)%npts),stat=ierr)
          call utils_alloc_check('internal_rad_dens_init', &
               'rad_dens(iat)%rad_npts',ierr)
          rad_dens(iat)%rad_npts = 0

          ! lpl: Exponential decay factor for each shell
          !      Defined as '-dr' for the complete exp decay
          !      term 'exp(-alpha*eta*dr)' (called 'exp_const' in CHARGEMOL)
          !      where 'alpha' is the 'density_decaying_exponent' and 'eta'
          !      is the 'constrait_term' in CHARGEMOL
          allocate(rad_dens(iat)%exp_factor(rad_dens(iat)%npts),stat=ierr)
          call utils_alloc_check('internal_rad_dens_init', &
               'rad_dens(iat)%exp_factor',ierr)
          rad_dens(iat)%exp_factor(1) = 1.0_DP ! 1st shell is irrelevant
          do it=2,rad_dens(iat)%npts
             rad_dens(iat)%exp_factor(it) = &
                  (rad_dens(iat)%rad(it-1) - rad_dens(iat)%rad(it))
          end do

          ! lpl: Total (core + valence) partial density
          allocate(rad_dens(iat)%total_part_dens(rad_dens(iat)%npts),stat=ierr)
          call utils_alloc_check('internal_rad_dens_init', &
               'rad_dens(iat)%total_part_dens',ierr)
          rad_dens(iat)%total_part_dens(:) = 0.0_DP

          ! lpl: Core density
          if ( include_coredens ) then
             allocate(rad_dens(iat)%core_dens(rad_dens(iat)%npts),stat=ierr)
             call utils_alloc_check('internal_rad_dens_init', &
                  'rad_dens(iat)%core_dens',ierr)
             rad_dens(iat)%core_dens = ref_dens(isp)%core_dens
          end if

          ! lpl: These arrays are not braodcasted as they are
          !      not needed universally, so allocate only on
          !      required procs where the atom resides
          if ( par%proc_of_atom(par%orig_atom(iat)) == pub_my_proc_id ) then
             ! lpl: Buffer for new partial density
             allocate(rad_dens(iat)%new_part_dens(rad_dens(iat)%npts),stat=ierr)
             call utils_alloc_check('internal_rad_dens_init', &
                  'rad_dens(iat)%new_part_dens',ierr)
             rad_dens(iat)%new_part_dens(:) = 0.0_DP

             ! lpl: Old ISA density
             allocate(rad_dens(iat)%old_part_dens(rad_dens(iat)%npts),stat=ierr)
             call utils_alloc_check('internal_rad_dens_init', &
                  'rad_dens(iat)%old_part_dens',ierr)
             rad_dens(iat)%old_part_dens(:) = 0.0_DP
          end if

          ! lpl: Partitioned IH (Yavg) density
          if ( include_refdens ) then
             allocate(rad_dens(iat)%Yavg_dens(rad_dens(iat)%npts),stat=ierr)
             call utils_alloc_check('internal_rad_dens_init', &
                  'rad_dens(iat)%Yavg_dens',ierr)
             rad_dens(iat)%Yavg_dens(:) = 0.0_DP
             ! lpl: These arrays are not braodcasted as they are
             !      not needed universally, so allocate only on
             !      required procs where the atom resides
             if ( par%proc_of_atom(par%orig_atom(iat)) == pub_my_proc_id ) then
                allocate(rad_dens(iat)%Yavg_sqrt(rad_dens(iat)%npts),stat=ierr)
                call utils_alloc_check('internal_rad_dens_init', &
                    'rad_dens(iat)%Yavg_sqrt',ierr)
                rad_dens(iat)%Yavg_sqrt(:) = 0.0_DP
                allocate(rad_dens(iat)%Yavg_ratio(rad_dens(iat)%npts),stat=ierr)
                call utils_alloc_check('internal_rad_dens_init', &
                     'rad_dens(iat)%Yavg_ratio',ierr)
                rad_dens(iat)%Yavg_ratio(:) = 0.0_DP
                allocate(rad_dens(iat)%ih_ratio(rad_dens(iat)%npts),stat=ierr)
                call utils_alloc_check('internal_rad_dens_init', &
                     'rad_dens(iat)%ih_ratio',ierr)
                rad_dens(iat)%ih_ratio(:) = 0.0_DP
             end if
          end if

          ! lpl: Reference IH (interpolated) density
          allocate(rad_dens(iat)%ref_dens(rad_dens(iat)%npts),stat=ierr)
          call utils_alloc_check('internal_rad_dens_init', &
               'rad_dens(iat)%ref_dens',ierr)
          rad_dens(iat)%ref_dens(:) = 0.0_DP

          nullify(rad_dens(iat)%target_dens)

          rad_dens(iat)%renorm_slope = 0.0_DP
          rad_dens(iat)%sum_part_dens = 0.0_DP
          rad_dens(iat)%sum_part_dens_sqrt = 0.0_DP

       end do

      if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_rad_dens_init'

     end subroutine internal_rad_dens_init

   !===========================================================================!

     subroutine internal_rad_dens_destroy(rad_dens)

       !-----------------------------------------------------------------------!
       ! lpl: Destroys radial densities DDEC_RAD_DENS for all atoms            !
       !-----------------------------------------------------------------------!
       ! InOut : rad_dens - Array of radial densities for ALL atoms            !
       !-----------------------------------------------------------------------!

       implicit none

       ! lpl: Input/Output
       type(DDEC_RAD_DENS), pointer, intent(inout) :: rad_dens(:)

       ! lpl: Local variables
       integer :: iat

      if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Entering internal_rad_dens_destroy'

       do iat=par%nat,1,-1

          deallocate(rad_dens(iat)%ref_dens,stat=ierr)
          call utils_dealloc_check('internal_rad_dens_destroy', &
               'rad_dens(iat)%ref_dens',ierr)

          if ( include_refdens ) then
             if ( par%proc_of_atom(par%orig_atom(iat)) == &
                  pub_my_proc_id ) then
                deallocate(rad_dens(iat)%ih_ratio,stat=ierr)
                call utils_dealloc_check('internal_rad_dens_destroy', &
                     'rad_dens(iat)%ih_ratio',ierr)
                deallocate(rad_dens(iat)%Yavg_ratio,stat=ierr)
                call utils_dealloc_check('internal_rad_dens_destroy', &
                     'rad_dens(iat)%Yavg_ratio',ierr)
                deallocate(rad_dens(iat)%Yavg_sqrt,stat=ierr)
                call utils_dealloc_check('internal_rad_dens_destroy', &
                     'rad_dens(iat)%Yavg_sqrt',ierr)
             end if
             deallocate(rad_dens(iat)%Yavg_dens,stat=ierr)
             call utils_dealloc_check('internal_rad_dens_destroy', &
                  'rad_dens(iat)%Yavg_dens',ierr)
          end if

          if ( par%proc_of_atom(par%orig_atom(iat)) == pub_my_proc_id ) then
             deallocate(rad_dens(iat)%old_part_dens,stat=ierr)
             call utils_dealloc_check('internal_rad_dens_destroy', &
                  'rad_dens(iat)%old_part_dens',ierr)
             deallocate(rad_dens(iat)%new_part_dens,stat=ierr)
             call utils_dealloc_check('internal_rad_dens_destroy', &
                  'rad_dens(iat)%new_part_dens',ierr)
          end if

          if ( include_coredens ) then
             deallocate(rad_dens(iat)%core_dens,stat=ierr)
             call utils_dealloc_check('internal_rad_dens_destroy', &
                  'rad_dens(iat)%core_dens',ierr)
          end if

          deallocate(rad_dens(iat)%total_part_dens,stat=ierr)
          call utils_dealloc_check('internal_rad_dens_destroy', &
               'rad_dens(iat)%total_part_dens',ierr)
          deallocate(rad_dens(iat)%exp_factor,stat=ierr)
          call utils_dealloc_check('internal_rad_dens_destroy', &
               'rad_dens(iat)%exp_factor',ierr)
          deallocate(rad_dens(iat)%rad_npts,stat=ierr)
          call utils_dealloc_check('internal_rad_dens_destroy', &
               'rad_dens(iat)%rad_npts',ierr)
          deallocate(rad_dens(iat)%rad,stat=ierr)
          call utils_dealloc_check('internal_rad_dens_destroy', &
               'rad_dens(iat)%rad',ierr)

          nullify(rad_dens(iat)%target_dens)

       end do

       deallocate(rad_dens,stat=ierr)
       call utils_dealloc_check('internal_rad_dens_destroy', &
            '=>rad_dens',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_rad_dens_destroy'

     end subroutine internal_rad_dens_destroy

   !===========================================================================!

     subroutine internal_rad_dens_avg_rad(rad_dens,ref_dens,local_box_fine)

       !-----------------------------------------------------------------------!
       ! lpl: Obtains the average radius of all Cartesian points lying within  !
       !      each radial shell                                                !
       !-----------------------------------------------------------------------!
       !                                                                       !
       !-----------------------------------------------------------------------!

       implicit none

       ! lpl: InOut
       type(DDEC_RAD_DENS), intent(inout) :: rad_dens(par%nat)

       ! lpl: Input
       type(DDEC_REF_DENS), intent(in) :: ref_dens(par%num_species)
       type(DDEC_LOCAL_BOX), intent(in) :: local_box_fine

       ! lpl: Local variables
       integer :: iat, isp, it1, it2, it3, ir
       integer :: max_interp_shell
       real(kind=DP) :: rad
       type(POINT) :: local_atom_ctr, local_grid_coord

       if ( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: Entering internal_rad_dens_avg_rad'

       do iat=1,par%nat

          rad_dens(iat)%rad_npts = 0
          rad_dens(iat)%rad = 0.0_DP

          local_atom_ctr = mdl%elements(par%orig_atom(iat))%centre &
               - local_box_fine%box_coord_origin(iat)

          do it1=1,local_box_fine%n1
             do it2=1,local_box_fine%n2
                do it3=1,local_box_fine%n3

                   ! lpl: Coord w.r.t. voxel 1,1,1 centre
                   local_grid_coord &
                        = real(it1 - 1,kind=DP)*grid%da1 &
                        + real(it2 - 1,kind=DP)*grid%da2 &
                        + real(it3 - 1,kind=DP)*grid%da3

                   ! lpl: Atom-centred radius (correct w.r.t. centre of voxel)
                   rad = geometry_magnitude(local_atom_ctr - local_grid_coord)
                   ir = aint(rad) + 1
                   rad_dens(iat)%rad(ir) = rad_dens(iat)%rad(ir) + rad
                   rad_dens(iat)%rad_npts(ir) = rad_dens(iat)%rad_npts(ir) + 1
                end do
             end do
          end do

          isp = mdl%elements(par%orig_atom(iat) )%species_number
          max_interp_shell = ref_dens(isp)%max_interp_shell
          rad_dens(iat)%rad(1:max_interp_shell) = &
               ref_dens(isp)%avg_rad(1:max_interp_shell)
          do ir=max_interp_shell+1,rad_dens(iat)%npts
             if ( rad_dens(iat)%rad_npts(ir) > 0 ) &
                rad_dens(iat)%rad(ir) = &
                     rad_dens(iat)%rad(ir)/(1.0_DP*rad_dens(iat)%rad_npts(ir))
          end do

          ! lpl: Reset just in case
          rad_dens(iat)%rad_npts = 0

       end do ! END do iat=1,par%nat

       call comms_barrier

       if ( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_rad_dens_avg_rad'

     end subroutine internal_rad_dens_avg_rad

   !===========================================================================!

   !---------------------------------------------------------------------------!
   ! lpl: Local density box 'DDEC_LOCAL_BOX' initialization and destruction    !
   !   1. internal_local_box_init                                              !
   !   2.
   !---------------------------------------------------------------------------!

     subroutine internal_local_box_init(local_box,grid,rcut_max, &
          include_refdens,include_coredens)

       !-----------------------------------------------------------------------!
       ! lpl: Initialize DDEC_LOCAL_BOX                                        !
       !-----------------------------------------------------------------------!
       ! InOut : local_box                                                     !
       ! Input : grid - GRID_INFO object containing lattice grid information   !
       !       : rcut_max - Maximum radius of DDEC shells. Used to determine   !
       !            size of the universal 'local_box'                          !
       !-----------------------------------------------------------------------!

       implicit none

       ! lpl: InOut
       type(DDEC_LOCAL_BOX), intent(inout) :: local_box

       ! lpl: Input
       type(GRID_INFO), intent(in) :: grid
       real(kind=DP), intent(in) :: rcut_max
       logical, intent(in) :: include_refdens, include_coredens

       ! lpl: Local variables
       integer :: iat
       real(kind=DP) :: sin_12, sin_23, sin_13
       integer :: n_offset_1, n_offset_2, n_offset_3
       integer :: n1_half, n2_half, n3_half

       integer :: debug1, debug2, debug3

       type(POINT) :: iat_ctr_coord
       type(POINT) :: half_box_diag, offset_coord


       if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Entering internal_local_box_init'

       ! lpl: Check if any DDEC_LOCAL_BOX info arrays have been allocated
       if ( allocated(local_box%total_val) .or. &
            allocated(local_box%total_pseudoval) .or. &
            allocated(local_box%ref_pseudoval) .or. &
            allocated(local_box%Yavg_pseudoval) .or. &
            allocated(local_box%core_val) .or. &
            allocated(local_box%core_pseudoval) .or. &
            allocated(local_box%box_coord_origin) .or. &
            allocated(local_box%box_index_origin) ) &
            call utils_abort(' ERROR: Attempting to initialize allocated &
                 &DDEC_LOCAL_BOX')

       ! lpl: Allocate DDEC_LOCAL_BOX location info w.r.t. global grid
       !      for each atom
       allocate(local_box%box_coord_origin(par%nat),stat=ierr)
       call utils_alloc_check('internal_local_box_init', &
            'local_box%box_coord_origin',ierr)
       allocate(local_box%box_index_origin(par%nat),stat=ierr)
       call utils_alloc_check('internal_local_box_init', &
            'local_box%box_index_origin',ierr)

       ! lpl: Calculate minimum box dimension required based on the largest
       !      shell of the DDEC radial densities
       sin_12 = sin(acos(grid_ua1.DOT.grid_ua2))
       sin_23 = sin(acos(grid_ua2.DOT.grid_ua3))
       sin_13 = sin(acos(grid_ua1.DOT.grid_ua3))

       ! lpl: Set local_box size (n_half = 1/2 number of points in each lattice
       !      vector direction, and add one point just in case number was
       !      rounded down
       n1_half = &
            int(ceiling(max(rcut_max/sin_12,rcut_max/sin_13)/grid_lda1)) + 1
       n2_half = &
            int(ceiling(max(rcut_max/sin_12,rcut_max/sin_23)/grid_lda2)) + 1
       n3_half = &
            int(ceiling(max(rcut_max/sin_13,rcut_max/sin_23)/grid_lda3)) + 1

       half_box_diag = real(n1_half,kind=DP)*grid%da1 &
            +  real(n2_half,kind=DP)*grid%da2 &
            +  real(n3_half,kind=DP)*grid%da3

       if(pub_debug_on_root) &
            write(debug_unit,'(a,3(1x,I5))') ' HALFBOX:',n1_half,n2_half,n3_half


       ! lpl: Get origin of local_box index (1,1,1) for each atom iat
       !      w.r.t. real simcell in terms of periodic starting
       !      indices and non-periodic Cartesian coordinates
       do iat=1,par%nat

          if(pub_debug_on_root) then
             write(debug_unit,'(a,I5)') ' IAT REAL LOCAL BOX INDEX ORIGIN: ',iat
             iat_ctr_coord = mdl%elements(par%orig_atom(iat))%centre
             debug1 = floor((iat_ctr_coord.DOT.grid_ua1)/grid_lda1)
             debug2 = floor((iat_ctr_coord.DOT.grid_ua2)/grid_lda2)
             debug3 = floor((iat_ctr_coord.DOT.grid_ua3)/grid_lda3)
             write(debug_unit,'(a,3(1x,I5))') ' IAT ORIGIN INDEX:', &
                  debug1, debug2, debug3
          end if


          ! lpl: Vector from atom center to box origin
          offset_coord = mdl%elements(par%orig_atom(iat))%centre - half_box_diag

          ! lpl: Number of grid points offset from atom center. +1 point so that
          !      atom ctr will be at n_half index
          n_offset_1 = floor((offset_coord.DOT.grid_ua1)/grid_lda1) + 1
          n_offset_2 = floor((offset_coord.DOT.grid_ua2)/grid_lda2) + 1
          n_offset_3 = floor((offset_coord.DOT.grid_ua3)/grid_lda3) + 1

          if(pub_debug_on_root) write(debug_unit,'(a,3(1x,I5))') &
               ' IAT LOCAL BOX INDEX OFFSET:',n_offset_1,n_offset_2,n_offset_3

          ! lpl: Non-periodic local_box origin in Cartesian coordinates
          !      FFTbox points refer to the middle of the voxel
          local_box%box_coord_origin(iat) &
               = (real(n_offset_1,kind=DP))*grid%da1 &
               + (real(n_offset_2,kind=DP))*grid%da2 &
               + (real(n_offset_3,kind=DP))*grid%da3

          if(pub_debug_on_root) write(debug_unit,'(a,3(1x,f14.7))') &
               ' IAT LOCAL BOX COORD ORIGIN:', &
               local_box%box_coord_origin(iat)%x, &
               local_box%box_coord_origin(iat)%y, &
               local_box%box_coord_origin(iat)%z

          ! lpl: Periodic local_box origin for use in fftbox routine
          if (n_offset_1 < 0) n_offset_1 = grid%n1 + mod(n_offset_1,grid%n1)
          if (n_offset_2 < 0) n_offset_2 = grid%n2 + mod(n_offset_2,grid%n2)
          if (n_offset_3 < 0) n_offset_3 = grid%n3 + mod(n_offset_3,grid%n3)

          if(pub_debug_on_root) then
             write(debug_unit,'(a,3(1x,I5))') &
                  ' IAT PERIODIC LOCAL BOX INDEX ORIGIN:', &
                  n_offset_1, n_offset_2, n_offset_3
             iat_ctr_coord = real(n1_half-1,kind=DP)*grid%da1 &
                  + real(n2_half-1,kind=DP)*grid%da2 &
                  + real(n3_half-1,kind=DP)*grid%da3 &
                  + local_box%box_coord_origin(iat) - iat_ctr_coord
             write(debug_unit,'(a,3(1x,f14.7))') &
                  ' IAT DISTANCE FROM LOCAL BOX EDGE: ', &
                  geometry_magnitude(iat_ctr_coord)
             write(debug_unit,'(a)') ''
          end if

          ! lpl: Global index origin of each local_box accounting for global
          !      cell periodicity. Index starts from 1.
          local_box%box_index_origin(iat)%ni1 = n_offset_1 + 1
          local_box%box_index_origin(iat)%ni2 = n_offset_2 + 1
          local_box%box_index_origin(iat)%ni3 = n_offset_3 + 1

       end do

       ! lpl: local_box size
       local_box%n1 = 2*n1_half
       local_box%n2 = 2*n2_half
       local_box%n3 = 2*n3_half

       ! lpl: Allocate local_box density grid arrays
       allocate(local_box%total_val(local_box%n1,local_box%n2, &
            local_box%n3), stat=ierr)
       call utils_alloc_check('internal_local_boxes_init', &
            'local_box%total_val',ierr)
       allocate(local_box%total_pseudoval(local_box%n1,local_box%n2, &
            local_box%n3), stat=ierr)
       call utils_alloc_check('internal_local_boxes_init', &
            'local_box%total_pseudoval',ierr)

       allocate(local_box%aim_val(local_box%n1,local_box%n2, &
            local_box%n3), stat=ierr)
       call utils_alloc_check('internal_local_boxes_init', &
            'local_box%aim_val',ierr)

       if ( include_refdens ) then
          allocate(local_box%ref_pseudoval(local_box%n1,local_box%n2, &
               local_box%n3), stat=ierr)
          call utils_alloc_check('internal_local_boxes_init', &
               'local_box%ref_pseudoval',ierr)
          allocate(local_box%Yavg_pseudoval(local_box%n1,local_box%n2, &
               local_box%n3), stat=ierr)
          call utils_alloc_check('internal_local_boxes_init', &
               'local_box%Yavg_pseudoval',ierr)
       end if

       if ( include_coredens ) then
          allocate(local_box%core_val(local_box%n1,local_box%n2, &
               local_box%n3), stat=ierr)
          call utils_alloc_check('internal_local_boxes_init', &
               'local_box%core_val',ierr)
          allocate(local_box%core_pseudoval(local_box%n1,local_box%n2, &
               local_box%n3), stat=ierr)
          call utils_alloc_check('internal_local_boxes_init', &
               'local_box%core_pseudoval',ierr)
       end if

       ! lpl: Initailize temp (buffer) array
       allocate(local_box%temp_val(local_box%n1,local_box%n2, &
            local_box%n3), stat=ierr)
       call utils_alloc_check('internal_local_boxes_init', &
            'local_box%temp_val',ierr)

       ! lpl: cell_grid_deposit_box is unable to deposit local box that is
       !      larger than the the global cell, so exit. Possible to resize
       !      local box by transferring content periodically into a smaller
       !      local_box, then deposit it into the global density grid, but
       !      this IS NOT CURRENTLY IMPLEMENTED
       ! lpl: COMMENT#00
       if ( ( local_box%n1 > grid%n1 .or. &
              local_box%n2 > grid%n2 .or. &
              local_box%n3 > grid%n3 ) ) call utils_abort(' ERROR: The &
            &local_box is larger than the supercell in &
            &internal_local_box_init')

       call comms_barrier

       if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_local_box_init'

     end subroutine internal_local_box_init

   !===========================================================================!

     subroutine internal_local_box_destroy(local_box)

       !-----------------------------------------------------------------------!
       ! lpl: Destroy DDEC_LOCAL_BOX                                           !
       !-----------------------------------------------------------------------!
       ! InOut : local_box                                                     !
       !-----------------------------------------------------------------------!

       implicit none

       type(DDEC_LOCAL_BOX), intent(inout) :: local_box

       if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Entering internal_local_box_destroy'

       deallocate(local_box%temp_val, stat=ierr)
       call utils_dealloc_check('internal_local_boxes_destroy', &
            'local_box%temp_val',ierr)

       if ( allocated(local_box%core_pseudoval) ) then
          deallocate(local_box%core_pseudoval,stat=ierr)
          call utils_dealloc_check('internal_local_box_destroy', &
               'local_box%core_pseudoval',ierr)
       end if
       if ( allocated(local_box%core_val) ) then
          deallocate(local_box%core_val,stat=ierr)
          call utils_dealloc_check('internal_local_box_destroy', &
               'local_box%core_val',ierr)
       end if

       if ( allocated(local_box%Yavg_pseudoval) ) then
          deallocate(local_box%Yavg_pseudoval,stat=ierr)
          call utils_dealloc_check('internal_local_box_destroy', &
               'local_box%Yavg_pseudoval',ierr)
       end if
       if ( allocated(local_box%ref_pseudoval) ) then
          deallocate(local_box%ref_pseudoval,stat=ierr)
          call utils_dealloc_check('internal_local_box_destroy', &
               'local_box%ref_pseudoval',ierr)
       end if

       deallocate(local_box%aim_val,stat=ierr)
       call utils_dealloc_check('internal_local_boxes_destroy', &
            'local_box%aim_val',ierr)

       deallocate(local_box%total_pseudoval,stat=ierr)
       call utils_dealloc_check('internal_local_box_destroy', &
            'local_box%total_pseudoval',ierr)
       deallocate(local_box%total_val,stat=ierr)
       call utils_dealloc_check('internal_local_box_destroy', &
            'local_box%total_val',ierr)

       deallocate(local_box%box_index_origin,stat=ierr)
       call utils_dealloc_check('internal_local_box_destroy', &
            'local_box%box_index_origin',ierr)
       deallocate(local_box%box_coord_origin,stat=ierr)
       call utils_dealloc_check('internal_local_box_destroy', &
            'local_box%box_coord_origin',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_local_box_destroy'

     end subroutine internal_local_box_destroy

   !===========================================================================!

  !----------------------------------------------------------------------------!
  ! lpl: Subrouties related to DDEC partitioning                               !
  !----------------------------------------------------------------------------!

     subroutine internal_rad_dens_to_box(local_box_val, local_box, &
          rad_dens_val, rad_dens)

       !-----------------------------------------------------------------------!
       ! lpl: Sums radial density contributions from all atoms into local_box  !
       !      grid of current atom 'current_iat'. input rad_dens and output    !
       !      local boxes from their respective pointers                       !
       !-----------------------------------------------------------------------!
       ! InOut : local_box - Will contain radial density in terms of the       !
       !            Cartesian grid                                             !
       ! Input : rad_dens - Array of radial densities for ALL atoms            !
       !-----------------------------------------------------------------------!

       implicit none

       ! lpl: Input Info
       type(DDEC_LOCAL_BOX), intent(in) :: local_box

       ! lpl: Output
       real(kind=DP), intent(out) :: local_box_val(:,:,:)
       ! lpl: Input
       type(DDEC_RAD_DENS), intent(in) :: rad_dens
       real(kind=DP), intent(in) :: rad_dens_val(:)

       ! lpl: Local variables
       integer :: iat, r_index
       integer :: it1, it2, it3

       real(kind=DP) :: r_dist
       type(POINT) :: local_atom_ctr, local_grid_coord

       if (pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: Entering internal_rad_dens_to_box'
       if (pub_debug) call timer_clock('ddec_internal_rad_dens_to_box',1)

       if ( size(local_box_val,1) /= local_box%n1 .or. &
            size(local_box_val,2) /= local_box%n2 .or. &
            size(local_box_val,3) /= local_box%n3 ) &
            call utils_abort(' ERROR in internal_rad_dens_to_box:&
                 & size(local_box_val) /= local_box')

       if ( size(rad_dens_val) /= rad_dens%npts ) &
            call utils_abort(' ERROR in internal_rad_dens_to_box:&
                 & size(rad_dens_val) /= rad_dens%npts')

       iat = rad_dens%iat

       local_atom_ctr = mdl%elements(par%orig_atom(iat))%centre &
            - local_box%box_coord_origin(iat)

       local_box_val = 0.0_DP

       ! lpl: Loop over all local_box grid points in 'local_box' to obtain
       !      total pseudodensity
       do it1=1,local_box%n1
          do it2=1,local_box%n2
             do it3=1,local_box%n3

                ! lpl: Coord w.r.t. to the origin of local_box
                local_grid_coord &
                     = real(it1 - 1,kind=DP)*grid%da1 &
                     + real(it2 - 1,kind=DP)*grid%da2 &
                     + real(it3 - 1,kind=DP)*grid%da3

                r_dist = geometry_magnitude(local_grid_coord - local_atom_ctr)

                ! lpl: Skip everything if rel_coord w.r.t. this grid > rcut
                r_index = aint(r_dist/rad_dens%dr) + 1
                if (r_index > rad_dens%npts) cycle

                ! lpl: Bin the radial densities into their
                !      appropriate local_box
                local_box_val(it1,it2,it3) = local_box_val(it1,it2,it3) &
                     + rad_dens_val(r_index)

             end do
          end do
       end do

       if (pub_debug) call timer_clock('ddec_internal_rad_dens_to_box',2)
       if (pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_rad_dens_to_box'

     end subroutine internal_rad_dens_to_box

   !===========================================================================!

     subroutine internal_update_rad_dens(rad_dens, ref_dens, atom_popn, &
          include_refdens, reshape_dens)

       !-----------------------------------------------------------------------!
       ! lpl: Update radial densities based on newly-calculated partial        !
       !      densities according to the DDEC partitioning scheme and          !
       !      broadcast to all procs                                           !
       !-----------------------------------------------------------------------!
       ! InOut : rad_dens, ref_dens - Arrays containing updated radial and IH  !
       !            reference densities for ALL atoms                          !
       ! Input : atom_popn - Array of atomic population for ALL atoms needed   !
       !            to update the IH reference densities                       !
       ! Flag  : cor_val - Whether to update core or total densities           !
       !-----------------------------------------------------------------------!

       use comms, only: pub_total_num_procs
       use utils, only: utils_unit

       implicit none

       ! lpl: InOut
       type(DDEC_RAD_DENS), intent(inout) :: rad_dens(:)
       type(DDEC_REF_DENS), intent(inout) :: ref_dens(:)

       ! lpl: Input
       real(kind=DP), intent(in) :: atom_popn(par%nat)
       logical, intent(in) :: include_refdens, reshape_dens

       ! lpl: Local variables
       integer :: it, reshape_it, lc_iat, iat, nproc
       integer, parameter :: max_reshape_it = 50
       logical :: break
       real(kind=DP) :: dens_scale, sum_part_dens_new, eta, exp_term

       integer :: internal_output_unit

       if ( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: Entering internal_update_rad_dens'
       if (pub_debug) call timer_clock('ddec_internal_update_rad_dens',1)

       internal_output_unit = utils_unit()

       ! lpl: Mix IH and ISA reference densities
       do lc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)

          iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

          ! lpl: Mix IH & ISA depending on IH fraction
          if ( include_refdens ) then
             ! lpl: Replace with new Yavg_dens from ref_dens*<ih_ratio>
             rad_dens(iat)%Yavg_dens = &
                  rad_dens(iat)%ref_dens*rad_dens(iat)%ih_ratio

             ! lpl: Update IH reference density based on new atom_popn
             call internal_update_ih_dens(rad_dens(iat)%ref_dens, &
                  ref_dens(mdl%elements(par%orig_atom(iat))%&
                  species_number), atom_popn(iat))

             ! lpl: Mix IH (Y^avg) and ISA densities, and update
             !      partial density (replaces 'part_dens')
             if(pub_ddec_IH_frac < 1.0_DP) then
                rad_dens(iat)%total_part_dens = &
                   (rad_dens(iat)%new_part_dens**ddec_ISA_frac) &
                   *(rad_dens(iat)%Yavg_dens**pub_ddec_IH_frac)
             else
                rad_dens(iat)%total_part_dens = rad_dens(iat)%Yavg_dens
             end if
          else
             rad_dens(iat)%total_part_dens = rad_dens(iat)%new_part_dens
          end if

          if (pub_debug_on_root) then
          if ( include_refdens ) then
             write(output_filename,'(2(a,I5.5))') &
                  'internal_sigma_atom_',par%orig_atom(iat),'_iter_',ddec_it
             open(unit=internal_output_unit,form="formatted", &
                  file=trim(adjustl(output_filename)),action="write")
             do it=1,rad_dens(iat)%npts
                write(internal_output_unit,'(4(2x,f14.7))') &
                    rad_dens(iat)%rad(it), rad_dens(iat)%new_part_dens(it), &
                    rad_dens(iat)%Yavg_dens(it), &
                    rad_dens(iat)%total_part_dens(it)
             end do
             close(unit=internal_output_unit)
          end if
          end if

          ! lpl: Compute necessary coefficient for reshaping
          if ( reshape_dens .and. include_refdens ) then

                ! lpl: Called 'first_integration_sum' in CHARGEMOL
                rad_dens(iat)%sum_part_dens = &
                     sum( rad_dens(iat)%total_part_dens*rad_dens(iat)%rad_npts )
                ! lpl: Called 'weighted_points' in CHARGEMOL
                rad_dens(iat)%sum_part_dens_sqrt = &
                     sum( sqrt(rad_dens(iat)%total_part_dens)*&
                     rad_dens(iat)%rad_npts )
                ! lpl: Store unmodified updated partial density in new_part_dens
                !      Called 'sigma' in CHARGEMOL
                rad_dens(iat)%new_part_dens = rad_dens(iat)%total_part_dens

          end if

       end do ! END do lc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)

       ! lpl: Reshape partial density (exp decay etc.)
       if ( reshape_dens .and. include_refdens ) then
          do reshape_it=1,max_reshape_it

             break = .true.

             do lc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
                iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

                if (pub_debug) then
                write(output_filename,'(3(a,I5.5))') 'exp_dens_attr_atom_', &
                     par%orig_atom(iat),'_iter_',ddec_it,'_trial_',reshape_it
                open(unit=internal_output_unit,form="formatted", &
                     file=trim(adjustl(output_filename)),action="write")
                end if

                ! lpl: Make sure density for 1st point is >= 0
                rad_dens(iat)%total_part_dens(1) = &
                     max(rad_dens(iat)%total_part_dens(1),0.0_DP)

                if (pub_debug) then
                write(internal_output_unit,'(1x,f14.7,4(2x,f14.7))') &
                     rad_dens(iat)%rad(1), &
                     rad_dens(iat)%total_part_dens(1), 1.0, 1.0, 1.0
                end if

                ! lpl: Exponential decay
                ! lpl: Density n(r+dr) depends on n(r) whereby
                !      n(r+dr) = min[ n(r+dr), n(r)*exp(-beta) ] and n(r+dr)
                !      is floored at 0. If n(r) = 0,  n(r+dr) = n(r+dr).
                !      n(dr) = n(dr) always
                do it=2,rad_dens(iat)%npts
                   ! lpl: Make sure density for every point is >= 0
                   rad_dens(iat)%total_part_dens(it) = &
                        max(rad_dens(iat)%total_part_dens(it),0.0_DP)
                   ! lpl: eta in Manz 2012, 'contraint_term' in CHARGEMOL
                   !      decay term is alredy pre-multiplied into 'exp_factor'
                   eta = 1.0_DP - (rad_dens(iat)%Yavg_ratio(it)/&
                        rad_dens(iat)%Yavg_sqrt(it))**2
                   exp_term = exp(eta*rad_dens(iat)%exp_factor(it))
                   if ( rad_dens(iat)%rad_npts(it) > 0 .and. &
                        rad_dens(iat)%total_part_dens(it-1) > 0.0_DP ) then
                      rad_dens(iat)%total_part_dens(it) = &
                           min(rad_dens(iat)%total_part_dens(it), &
                           rad_dens(iat)%total_part_dens(it-1)*exp_term)
                   end if

                   if (pub_debug) then
                   write(internal_output_unit,'(1x,f14.7,4(2x,f14.7))') &
                        rad_dens(iat)%rad(it), &
                        rad_dens(iat)%total_part_dens(it), &
                        eta, exp_term, rad_dens(iat)%exp_factor(it)
                   end if
                end do

                if (pub_debug) then
                close(unit=internal_output_unit)

                output_filename = ''
                write(output_filename,'(3(a,I5.5))') 'exp_dens_atom_', &
                     par%orig_atom(iat),'_iter_',ddec_it,'_trial_',reshape_it
                call ddec_write_dens(rad_dens(iat)%total_part_dens, &
                     rad_dens(iat)%rad, rad_dens(iat)%npts, &
                     trim(adjustl(output_filename)))
                end if

                ! lpl: Update sum of reshaped density
                !      Called 'second_integration_sum' in CHARGEMOL
                sum_part_dens_new = &
                     sum(rad_dens(iat)%total_part_dens*rad_dens(iat)%rad_npts)

                ! lpl: Compute 'Delta' (Eq. 76 in JCTC 2012)
                !      'dens_scale' is actually 'Delta' in JCTC 2012
                !      This looks like computing a normalization constraint
                !      weighted by 1/sqrt(sigma) (c.f. JCTC 2012 Eq. 65, 76)
                !      for which 'dens_scale' a.k.a. 'Delta' should converge
                !      to zero in the end
                if ( rad_dens(iat)%sum_part_dens_sqrt &
                     > pub_ddec_zero_thresh ) then
                   dens_scale = ( rad_dens(iat)%sum_part_dens - &
                        sum_part_dens_new)/rad_dens(iat)%sum_part_dens_sqrt
                else
                   dens_scale = 0.0_DP
                end if

                if (pub_debug) then
                write(debug_unit,'(3(a,I5.5))') 'dens_scale_atom_', &
                     par%orig_atom(iat),'_iter_',ddec_it,'_trial_',reshape_it
                write(debug_unit,'(5(2x,a,2x,f14.7))') 'sum_1  =', &
                     rad_dens(iat)%sum_part_dens, 'sum_2 =', &
                     sum_part_dens_new, &
                     'sum_sqrt_1 =', rad_dens(iat)%sum_part_dens_sqrt, &
                     'dens_scale =', dens_scale, 'renorm_slope =', &
                     rad_dens(iat)%renorm_slope
                end if

                ! lpl: Check convergence of reshaping of 'G_A(r_A)'
                !      (c.f. JCTC 2012) which occurs when 'dens_scale'
                !      (called 'Delta in Eq. 76 of JCTC 2012) is below
                !      a particular threshold
                ! lpl: If 'dens_scale' a.k.a. 'Delta' is above threshold
                !      then compute Eq. 77 in JCTC 2012 (i.e. add
                !      'Delta' x 'sqrt(sigma)' to 'G_A(r_A)')
                if ( abs(dens_scale) > pub_ddec_zero_thresh ) then
                   rad_dens(iat)%total_part_dens = &
                        rad_dens(iat)%total_part_dens + &
                        dens_scale*sqrt(rad_dens(iat)%new_part_dens)
                   break = .false.
                end if

             end do ! END do lc_iat=1,par%num_atoms_on_proc(nproc)

             if ( break ) exit

          end do ! END do iit=1,max_iit

          ! lpl: Guarantees that valence density (total - core) will be > 0
          !      Eq. 74 in JCTC 2012 - 'renorm_slope' is 'u_A'
          !      The 'dens_scale' variable is now reused as 'lambda_A'
          do lc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
             iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

             if ( rad_dens(iat)%renorm_slope > overlap_thresh ) then
                dens_scale = max( rad_dens(iat)%renorm_slope - &
                     atom_popn(iat)/rad_dens(iat)%renorm_slope,1.0_DP )
             else
                dens_scale = 1.0_DP
             end if

             if (pub_debug) then
             write(debug_unit,'(2(a,I5.5))') 'renorm_atom_', &
                  par%orig_atom(iat), '_iter_', ddec_it
             write(debug_unit,'(1(2x,a,2x,f14.7))') &
                  'renorm_scale =', dens_scale
             end if

             rad_dens(iat)%total_part_dens = &
                  dens_scale*rad_dens(iat)%total_part_dens
          end do
       end if ! END if ( reshape_dens )

       ! lpl: Synchronize densities over all procs
       do nproc=0,pub_total_num_procs-1
          do lc_iat=1,par%num_atoms_on_proc(nproc)
             iat = lc_iat + par%first_atom_on_proc(nproc) - 1

             ! lpl: Broadcast updated radial densities
             call comms_bcast(nproc,rad_dens(iat)%total_part_dens)
             if ( include_refdens ) then
                call comms_bcast(nproc,rad_dens(iat)%ref_dens)
                call comms_bcast(nproc,rad_dens(iat)%Yavg_dens)
             end if

          end do
       end do
       call comms_barrier

       if (pub_debug) call timer_clock('ddec_internal_update_rad_dens',2)
       if ( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_update_rad_dens'

     end subroutine internal_update_rad_dens

   !===========================================================================!

     subroutine internal_update_rad_dens_core(rad_dens,reshape_dens)

       ! lpl: InOut
       type(DDEC_RAD_DENS), intent(inout) :: rad_dens(:)

       ! lpl: Input
       logical, intent(in) :: reshape_dens

       ! lpl: Local variables
       integer :: it, lc_iat, iat, nproc
       real(kind=DP) :: exp_term

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' Entering internal_update_rad_dens_core'

       do nproc=0,pub_total_num_procs-1
          do lc_iat=1,par%num_atoms_on_proc(nproc)
             iat = lc_iat + par%first_atom_on_proc(nproc) - 1

             if ( nproc == pub_my_proc_id ) then
                if ( reshape_dens ) then
                   ! lpl: Core exp decay
                   rad_dens(iat)%core_dens(1) = &
                        max(rad_dens(iat)%new_part_dens(1),0.0_DP)
                   do it=2,rad_dens(iat)%npts
                      ! lpl: 'contraint_term' in CHARGEMOL or in this
                      !      case it is 'eta' that is always equal to 'b'
                      !      (Eq. 55 in JCTC 2012, and b='core_exp_factor')
                      exp_term = &
                           exp(core_exp_factor*rad_dens(iat)%exp_factor(it))
                      if ( rad_dens(iat)%rad_npts(it) > 0 .and. &
                           rad_dens(iat)%new_part_dens(it-1) > 0.0_DP ) then
                         rad_dens(iat)%core_dens(it) = &
                              max( min(rad_dens(iat)%new_part_dens(it), &
                              rad_dens(iat)%core_dens(it-1)*exp_term),0.0_DP )
                      else
                         rad_dens(iat)%core_dens(it) = &
                              max(rad_dens(iat)%new_part_dens(it),0.0_DP)
                      end if
                   end do
                else
                   rad_dens(iat)%core_dens = rad_dens(iat)%new_part_dens
                end if
             end if

             ! lpl: Broadcast updated radial densities
             call comms_bcast(nproc,rad_dens(iat)%core_dens)
          end do
       end do
       call comms_barrier

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' Leaving internal_update_rad_dens_core'

     end subroutine internal_update_rad_dens_core

   !===========================================================================!

     real(kind=DP) function internal_ddens_max(rad_dens)

       implicit none

       type(DDEC_RAD_DENS), intent(inout) :: rad_dens(:)

       integer :: iat, lc_iat

       if ( pub_debug_on_root )  write(stdout,'(a)') ' Entering internal_ddens_max'

       internal_ddens_max = 0.0_DP
       do lc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
          iat = lc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1

          internal_ddens_max = max( maxval(abs(rad_dens(iat)%new_part_dens &
               - rad_dens(iat)%old_part_dens)),internal_ddens_max )
          rad_dens(iat)%old_part_dens = rad_dens(iat)%new_part_dens
       end do
       call comms_reduce('MAX',internal_ddens_max)

       if ( pub_debug_on_root )  write(stdout,'(a)') ' Leaving internal_ddens_max'

     end function internal_ddens_max

  !===========================================================================!

     subroutine internal_debug_write_cube(local_val,local_box,iat,filename,valname)

       !-----------------------------------------------------------------------!
       ! lpl: Writes cube files of local_box for debugging which is being      !
       !      pointed to by the auxiliary pointer 'val_ptr' of local_box       !
       !-----------------------------------------------------------------------!

       use utils, only: utils_unit, utils_open_unit_check, &
            utils_close_unit_check

       implicit none

       ! lpl: Inputs
       type(DDEC_LOCAL_BOX), intent(in) :: local_box
       real(kind=DP), intent(in) :: local_val(:,:,:)
       integer, intent(in) :: iat
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: valname

       ! lpl: Local variables
       integer :: output, ierr
       integer :: it, it1, it2, it3

       if ( size(local_val,1) /= local_box%n1 .or. &
            size(local_val,2) /= local_box%n2 .or. &
            size(local_val,3) /= local_box%n3 ) &
            call utils_abort(' ERROR in internal_debug_write_cube:&
                 & size(local_val) /= local_box')

       output = utils_unit()
       open(unit=output,form="formatted", &
            file=(trim(filename)//trim(valname)//'.cube'), &
            action="write",iostat=ierr)
       call utils_open_unit_check('internal_debug_write_cube', &
            trim(filename)//trim(valname)//'.cube',ierr)

       write(output,'(a,/a)') 'DDEC DEBUG',trim(adjustl(filename))
       write(output,'(i4,3(f13.5))') 1,local_box%box_coord_origin(iat)%x, &
            local_box%box_coord_origin(iat)%y,local_box%box_coord_origin(iat)%z
       write(output,'(i4,3(f13.5))') local_box%n1,grid%da1%x,grid%da1%y,grid%da1%z
       write(output,'(i4,3(f13.5))') local_box%n2,grid%da2%x,grid%da2%y,grid%da2%z
       write(output,'(i4,3(f13.5))') local_box%n3,grid%da3%x,grid%da3%y,grid%da3%z
       write(output,'(i4,4(f13.5))') &
            mdl%elements(par%orig_atom(iat))%atomic_number, &
            real(mdl%elements(par%orig_atom(iat))%ion_charge,kind=DP), &
            mdl%elements(par%orig_atom(iat))%centre%x , &
            mdl%elements(par%orig_atom(iat))%centre%y , &
            mdl%elements(par%orig_atom(iat))%centre%z

       it = 1
       do it1=1,local_box%n1
          do it2=1,local_box%n2
             do it3=1,local_box%n3
                if(mod(it,6) == 0) write(output,'(a)') ''
                write(output,'(E13.5)',advance='no') &
                     local_val(it1,it2,it3)
                it = it + 1
             end do
          end do
       end do

       close(unit=output,iostat=ierr)
       call utils_close_unit_check('internal_debug_write_cube', &
            trim(filename)//trim(valname)//'.cube',ierr)

     end subroutine internal_debug_write_cube

  !===========================================================================!

  !----------------------------------------------------------------------------!
  ! lpl: Reference density 'DDEC_REF_DENS' initialization and destruction      !
  !      subroutines                                                           !
  !   1. internal_ref_dens_init                                                !
  !   2. internal_ref_dens_destroy                                             !
  !----------------------------------------------------------------------------!

     subroutine internal_ref_dens_init(ref_dens)

       !-----------------------------------------------------------------------!
       ! lpl: Initialize DDEC_REF_DENS for all species on all procs            !
       !-----------------------------------------------------------------------!
       ! InOut : ref_dens - Array of reference densities for ALL atoms         !
       !-----------------------------------------------------------------------!

! lpl: c.f. COMMENT#01
!       use atom, only: atom_gs_occ
       use comms, only: pub_root_proc_id
       use rundat, only: pub_ddec_ionic_range, &
! lpl: c.f. COMMENT#01
!            pub_ddec_rcomp, pub_ddec_target_radius, &
            pub_ddec_min_shell_dens, pub_ddec_c3_refdens, pub_ddec_refdens_path

       implicit none

       ! lpl: Input
       type(DDEC_REF_DENS), pointer, intent(inout) :: ref_dens(:)

       ! lpl: Local variables
       integer :: isp, iion, it, num_ions, iions(300,2)
       integer, parameter :: max_anion = -18
       integer, parameter :: max_ddec_ref_alloc = 1000
       integer, allocatable :: min_ion_state(:), max_ion_state(:)
       real(kind=DP) :: r, dr
       real(kind=DP) :: c1, c2
       character(len=80) :: error_string

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Entering internal_ref_dens_init'

       ! lpl: Number of allocated and initialized reference densities
       ddec_ref_alloc_list%num = 0

       ! lpl: Check if allocated (should not) then initialize
       if(allocated(ddec_ref_alloc_list%alloc)) call utils_abort(' ERROR: &
            &ddec_ref_alloc_list already allocated in internal_ref_dens_init')
       allocate(ddec_ref_alloc_list%alloc(max_ddec_ref_alloc,2),stat=ierr)
       call utils_alloc_check('internal_ref_dens_init', &
            'ddec_ref_alloc_list%alloc',ierr)

       ! lpl: Allocate ref_dens for all species
       if(associated(ref_dens)) call utils_abort(' ERROR: Attempting to &
            &initialize associated ref_dens in internal_ref_dens_init')
       allocate(ref_dens(par%num_species),stat=ierr)
       call utils_alloc_check('internal_ref_dens_init','=>ref_dens',ierr)

       if(.not. associated(ref_dens)) call utils_abort(' ERROR: Failed to &
            &initialize ref_dens in internal_ref_dens_init')

       allocate(min_ion_state(par%num_species),stat=ierr)
       call utils_alloc_check('internal_ref_dens_init','min_ion_state',ierr)
       allocate(max_ion_state(par%num_species),stat=ierr)
       call utils_alloc_check('internal_ref_dens_init','max_ion_state',ierr)

       do isp=1,par%num_species
          ref_dens(isp)%initialized = .false.
       end do
       call comms_barrier

       ! lpl: Initialize neutral reference densities, rcut and npts for
       !      all species
       do isp=1,par%num_species

          if ( ref_dens(isp)%initialized ) then
             error_string = ''
             write(error_string,'(a,I4,a)') ' ERROR: ref_dens(',isp, &
                  ') already initialized.'
             call utils_abort( trim(adjustl(error_string)) )
          end if

          it = -1 ! lpl: Error trapping

          ! lpl: Find example of current species
          do it=1,par%nat
             if(mdl%elements(it)%species_number == isp) exit
          end do

!djc          if (pub_on_root) then
!             write(stdout,'(a)') '=======================================&
!                  &=========='
!             write(stdout,'(a,a4)') ' Reference Density initialization&
!                  & for ',mdl%elements(it)%symbol
!             write(stdout,'(a)') '---------------------------------------&
!                  &----------'
!djc          end if
          ! lpl: Copy species to ref_dens(isp)%element
          ref_dens(isp)%element = mdl%elements(it)

          ! lpl: Initialize rcut and npts (same for all atoms for now)
          ref_dens(isp)%rcut = pub_ddec_rad_rcut
          ref_dens(isp)%npts = pub_ddec_rad_npts
          ref_dens(isp)%dr   = ref_dens(isp)%rcut/(1.0_DP*ref_dens(isp)%npts)

          ! lpl: Radial grid
          allocate(ref_dens(isp)%rad(ref_dens(isp)%npts),stat=ierr)
          call utils_alloc_check('internal_ref_dens_init', &
               'ref_dens(isp)%rad',ierr)
          ref_dens(isp)%rad = 0.0_DP

          ! lpl: Compute the hypothetical 'average radius' for each shell
          allocate(ref_dens(isp)%avg_rad(ref_dens(isp)%npts),stat=ierr)
          call utils_alloc_check('internal_ref_dens_init', &
               'ref_dens(isp)%avg_rad',ierr)

          do it=1,ref_dens(isp)%npts
             r = ((it-1)/real(ref_dens(isp)%npts,kind=DP))*ref_dens(isp)%rcut
             dr = ref_dens(isp)%dr
             ! lpl: <r> = integral[ dr, 4*PI*r^2*r ]/integral [dr, 4*PI*r^2 ]
             !      from r to r+dr i.e. <r> weighted by surface area
             ref_dens(isp)%avg_rad(it) = 0.75_DP*(2.0_DP*r + dr)*&
                  ((2.0_DP*r*(r+dr) + dr**2)/(3.0_DP*r*(r+dr) + dr**2))
          end do

          ! lpl: Determine the radial shell index 'max_interp_shell' where
          !      number of grid points is less than the minimum allowed
          !      surface point density
          dr = ref_dens(isp)%dr
          c1 = (3.0_DP/FOUR_PI)*pub_ddec_min_shell_dens*vox_vol
          if ( c1 < dr**3 ) call utils_abort(' ERROR: Unable to determine &
               &max_interp_shell in internal_ref_dens_init')
          c2 = (4.0_DP*c1 - dr**3)/(3.0_DP*dr)
          c1 = 0.5_DP*( sqrt(c2) - dr )

          if ( c1 < 0.0_DP ) call utils_abort(' ERROR: Unable to determine &
               &max_interp_shell in internal_ref_dens_init')
          ref_dens(isp)%max_interp_shell = floor(c1/dr)

          ! lpl: If 'max_interp_shell' index is larger than the number
          !      of radial shells 'npts' i.e. all shells do not meet
          !      minimum grid point density, set (max_interp_shell = npts)
          if ( ref_dens(isp)%max_interp_shell > ref_dens(isp)%npts ) &
               ref_dens(isp)%max_interp_shell = ref_dens(isp)%npts

          if( ref_dens(isp)%max_interp_shell > 0 ) then

             ! lpl: Allocate list of coordinates on the surface of
             !      the shells requiring density interpolation
             allocate(ref_dens(isp)%rad_grid(&
                  ref_dens(isp)%max_interp_shell),stat=ierr)
             call utils_alloc_check('internal_ref_dens_init', &
                  'ref_dens(isp)%rad_grid',ierr)

             ! lpl: Generate and store list of coordinates for each
             !      radial shell
             do it=1,ref_dens(isp)%max_interp_shell

                ! lpl: Average radius of shell
                r = ref_dens(isp)%avg_rad(it)

                ! lpl: Number of points required for this shell
                ref_dens(isp)%rad_grid(it)%npts = pub_ddec_min_shell_dens*&
                     (r/ref_dens(isp)%avg_rad(ref_dens(isp)%max_interp_shell))

                ! lpl: Allocate necessay number of points to be stored
                allocate( ref_dens(isp)%rad_grid(it)% &
                     coord(ref_dens(isp)%rad_grid(it)%npts), stat=ierr)
                call utils_alloc_check('internal_ref_dens_init', &
                     'ref_dens(isp)%rad_grid(it)%coord',ierr)

                ! lpl: Generate and store the list of coordinates on
                !      the surface of each radial shell
                call internal_radial_grid(ref_dens(isp)%rad_grid(it),r)

                if ( pub_debug_on_root ) write(stdout,'(1x,a,I3,a,f10.5,a,I5)') &
                     'shell=',it,'  avg_rad=',r,'  npts=', &
                     ref_dens(isp)%rad_grid(it)%npts

             end do

             ! lpl: Print info
!djc             if ( pub_on_root ) then
!                write(stdout,'(a,f10.5,a)') &
!                  ' Max interpolation radius = ', &
!                  ref_dens(isp)%max_interp_shell*ref_dens(isp)%dr,' bohr'
!                write(stdout,'(a)') '------------------------------------&
!                     &-------------'
!djc             end if
          end if

          ! lpl: Determine how radial shell should be defined (lower, mid,
          !      upper, or avg_rad)
          select case ( trim(adjustl(pub_ddec_rad_shell_mode)) )
             case ( 'LOWER' )
                do it=1,ref_dens(isp)%npts
                   ref_dens(isp)%rad(it) = &
                        ((it-1)/real(ref_dens(isp)%npts,kind=DP))&
                        *ref_dens(isp)%rcut
                end do
             case ( 'MIDDLE' )
                do it=1,ref_dens(isp)%npts
                   ref_dens(isp)%rad(it) = &
                        ((it-0.5_DP)/real(ref_dens(isp)%npts,kind=DP))&
                        *ref_dens(isp)%rcut
                end do
             case ( 'UPPER' )
                do it=1,ref_dens(isp)%npts
                   ref_dens(isp)%rad(it) = &
                        (it/real(ref_dens(isp)%npts,kind=DP))&
                        *ref_dens(isp)%rcut
                end do
             case ( 'AVERAGE' )
                ref_dens(isp)%rad = ref_dens(isp)%avg_rad
             case default
                call utils_abort(' ERROR in internal_ref_dens_init: &
                     &Unknown pub_ddec_rad_shell_mode')
          end select

          ! lpl: Radial shell volumes (follows CHARGEMOL, for calculating
          !      density from radial shells)
          allocate(ref_dens(isp)%rad_shell_vol(ref_dens(isp)%npts),stat=ierr)
          call utils_alloc_check('internal_ref_dens_init', &
               'ref_dens(isp)%rad_shell_vol',ierr)

          do it=1,ref_dens(isp)%npts
             r = it
             ref_dens(isp)%rad_shell_vol(it) = &
                  (FOUR_PI/3.0_DP)*(ref_dens(isp)%dr)**3*(3*it*(it-1)+1)
          end do

          ! lpl: Obtain effective nuclear charge (i.e. effective charge of
          !      the bare ion)
          if ( ref_dens(isp)%element%ion_charge <= &
               ref_dens(isp)%element%atomic_number ) then

             ! lpl: MOD#03 Make sure 'eff_nuc_charge' is integral
             if ( ref_dens(isp)%element%ion_charge /= &
                  nint(ref_dens(isp)%element%ion_charge) ) &
                  call utils_abort(' ERROR in internal_ref_dens_init:&
                       & ref_dens%element%ion_charge is not integral')

             ref_dens(isp)%eff_nuc_charge = &
                  nint(ref_dens(isp)%element%ion_charge)

          else
             call utils_abort(' ERROR in internal_ref_dens_init: Somehow, &
                  & ref_dens(isp)%ion_charge > atomic number')
          end if

          ! lpl: Obtain 'ddec_nuc_charge' that will be subtracted from
          !      'atom_popn' to obtain correct atomic charges
          if ( include_coredens ) then
             ref_dens(isp)%ddec_nuc_charge = ref_dens(isp)%element%atomic_number
          else
             ref_dens(isp)%ddec_nuc_charge = ref_dens(isp)%eff_nuc_charge
          end if

          ! lpl: Determine largest anion that will be allowed, determined
          !      arbitrarily by 'max_anion'
          ref_dens(isp)%max_anion = max_anion
          ! lpl: Largest cation limited by 'eff_nuc_charge'
          ref_dens(isp)%max_cation = ref_dens(isp)%eff_nuc_charge

          ! lpl: Allocate pointers to radial density profiles
          !      Initially this was supposed to allow new reference densities
          !      to be generated as needed in the middle of the DDEC iteration,
          !      but as ionic state generation turns out to be more complicated
          !      than anticipated, the current code won't do it. This is left
          !      open for future expansion. Currently pointers of index ranging
          !      from (-) 'max_anion' to (+) 'max_cation' will be created,
          !      corresponding to the charge of the ionic state. Allocation will
          !      possibly be made for the ionic state needed, and all ionic
          !      states between those currently available and the one needed.
          allocate(ref_dens(isp)%ion_dens( ref_dens(isp)%max_anion: &
               ref_dens(isp)%max_cation), stat=ierr)
          call utils_alloc_check('internal_ref_dens_init', &
               'ref_dens(isp)%ion_dens',ierr)

          ! lpl: MOD#07 - Additional flag 'ref_dens%init_status' to mark init
          !      status of reference densities for each ion
          allocate(ref_dens(isp)%init_status(ref_dens(isp)%max_anion: &
               ref_dens(isp)%max_cation), stat=ierr)
          call utils_alloc_check('internal_ref_dens_init', &
               'ref_dens(isp)%init_status',ierr)
          ref_dens(isp)%init_status = .false.

          allocate(ref_dens(isp)%c3flag(ref_dens(isp)%max_anion: &
               ref_dens(isp)%max_cation), stat=ierr)
          call utils_alloc_check('internal_ref_dens_init', &
               'ref_dens(isp)%c3flag',ierr)
          ref_dens(isp)%c3flag = .false.
          ! lpl: Neutral species is always left untouched, so vanilla version
          !      is effectively c3 (flag '2' indicates 'initialized as c3')
          ref_dens(isp)%c3flag(0) = .true.

          ! lpl: Set initial range of ionic states to be generated
          if(pub_ddec_IH_frac > 0.0_DP) then
             min_ion_state(isp) = &
                  max(-1*abs(pub_ddec_ionic_range),ref_dens(isp)%max_anion)
             max_ion_state(isp) = &
                  min(abs(pub_ddec_ionic_range),ref_dens(isp)%max_cation)
          else
             min_ion_state(isp) = 0
             max_ion_state(isp) = 0
          end if

          ! lpl: Allocate array for rcomp info
          allocate(ref_dens(isp)%rcomp(ref_dens(isp)%max_anion: &
               ref_dens(isp)%max_cation,6),stat=ierr)
          call utils_alloc_check('internal_ref_dens_init', &
               'ref_dens(isp)%rcomp',ierr)
          ! lpl: '-1' is flag for uninitalized variable
          ref_dens(isp)%rcomp = -1.0_DP
          ! lpl: Allocate array for reference density electronic config info
          !      This string is fed to the atomic solver 'conf=' string and
          !      uses the same nomenclature
          allocate(ref_dens(isp)%refconf(ref_dens(isp)%max_anion: &
               ref_dens(isp)%max_cation),stat=ierr)
          call utils_alloc_check('internal_ref_dens_init', &
               'ref_dens(isp)%refconf',ierr)
          ref_dens(isp)%refconf = ''

          call comms_barrier

          ! lpl: Obtain rcomp generation info
          if ( pub_on_root ) call internal_parse_rcomp(ref_dens(isp), &
               min_ion_state(isp),max_ion_state(isp))
          call comms_bcast( pub_root_proc_id, ref_dens(isp)%rcomp(:,:) )
          call comms_bcast( pub_root_proc_id, ref_dens(isp)%refconf )

       end do ! END do isp=1,par%num_species

       ! lpl: Check configurations and generate or read densities
       do isp=1,par%num_species

          ! lpl: Generate range of ionic states
          num_ions = 0
          iions = -999
          do iion=min_ion_state(isp),max_ion_state(isp)
             num_ions = num_ions + 1
             iions(num_ions,1) = iion ! Stores ionic number
             iions(num_ions,2) = 0    ! Stores flag to indicate read mode
          end do

          ! lpl: Mark ion to be read from separate core density
          !      file if available, only if 'pub_ddec_use_coredens'
          !      is set to .true.
          if ( include_coredens ) then
             if ( iions(num_ions,1) == ref_dens(isp)%eff_nuc_charge ) then
                iions(num_ions,2) = 1
             else
                num_ions = num_ions + 1
                iions(num_ions,1) = ref_dens(isp)%eff_nuc_charge
                iions(num_ions,2) = 1
             end if
          end if

          if (pub_debug) then
          do it=1,num_ions
             iion = iions(it,1)
             write(stdout,'(a,2(1x,I2),1x,a,1x,I2)') &
                  ' READ: ', iion, nint(ref_dens(isp)%rcomp(iion,6)), &
                  trim(adjustl(ref_dens(isp)%refconf(iion))), iions(it,2)
          end do
          end if

          nullify(ref_dens(isp)%core_dens)

! lpl: c.f. COMMENT#01
!          ! lpl: Allocate array to store neutral orbital info if
!          !      'solve' is used to generate neutral refdens (needed
!          !      to generate cationic densities using 'solve')
!          if ( nint(ref_dens(isp)%rcomp(0,6)) == 1 ) then
!             ! lpl: Allocate array for neutral orbital info
!             allocate(ref_dens(isp)%neutral_orb_info( &
!                atom_gs_occ(2,0,'Size'),3),stat=ierr)
!             call utils_alloc_check('internal_ref_dens_init', &
!                  'ref_dens(isp)%neutral_orb_info',ierr)
!             ref_dens(isp)%neutral_orb_info = -1.0_DP
!          end if

          ! lpl: Verify that all ionic configs exist
          if ( pub_on_root ) then
             write(stdout,'(a)') '=================================&
                  &==================================='
             write(stdout,'(a,a)') ' rcomp parameters for ', &
                  ref_dens(isp)%element%species_id
!djc             write(stdout,'(2x,a5,5(2x,a10))') &
!                  'ion','rcomp','rcomp_min','rcomp_max','rcomp_incr', &
!djc                  'solver_rad'
             write(stdout,'(a)') '---------------------------------&
                  &-----------------------------------'

             ! lpl: Verify all configurations
             do it=1,num_ions
                iion = iions(it,1)

                ! lpl: MOD#04 - Abort if 'reshape_dens' or 'include_coredens'
                !      are enabled and any one of the reference densities
                !      have been generated via 'SOLVE' since they are
                !      pseudized densities incompatible with DDEC/c3
                !      '1' and '3' means densities generated using the solver
                !      for valence (1) and core (3) densities
                if ( nint(ref_dens(isp)%rcomp(iion,6)) == 1 .or. &
                     nint(ref_dens(isp)%rcomp(iion,6)) == 3 ) then
                   if ( pub_ddec_reshape_dens ) call utils_abort(' ERROR in&
                        & internal_ref_dens_init: pub_ddec_reshape_dens&
                        & = .true. is incompatible with autogenerated&
                        & atomic solver reference densities')
                   if ( include_coredens ) call utils_abort(' ERROR in&
                        & internal_ref_dens_init: pub_ddec_reshape_dens&
                        & = .true. is incompatible with autogenerated&
                        & atomic solver reference densities')
                end if

                if ( nint(ref_dens(isp)%rcomp(iion,6)) <= 0 ) then
                   ! Densities with missing info
                   write(error_string,'(a,1x,I3)') &
                        ref_dens(isp)%element%species_id, iion
                   call utils_abort(' ERROR in internal_ref_dens_init:&
                        & A reference density configuration has been&
                        & left unspecified ('//trim(error_string)//')')
                end if

                if ( iions(it,2) == 0 ) then
                   ! 'iions(it,2) == 0' means valence density
                   if ( nint(ref_dens(isp)%rcomp(iion,6)) /= 1 .and. &
                        nint(ref_dens(isp)%rcomp(iion,6)) /= 2 ) then
                      write(error_string,'(a,1x,I3,1x,I2)') &
                           ref_dens(isp)%element%species_id, iion, &
                           nint(ref_dens(isp)%rcomp(iion,6))
                      call utils_abort(' ERROR in internal_ref_dens_init:&
                           & Somehow, a valence density configuation was&
                           & not declared as a valence density ('// &
                           trim(error_string)//')')
                   end if
                   write(stdout,'(a)',advance='no') '  '
                else if ( iions(it,2) == 1 ) then
                   ! 'iions(it,2) == 1' means core density
                   if ( .not. include_coredens ) then
                      write(error_string,'(a,1x,I3,1x,I2)') &
                           ref_dens(isp)%element%species_id, iion, &
                           nint(ref_dens(isp)%rcomp(iion,6))
                      call utils_abort(' ERROR in internal_ref_dens_init:&
                           & Somehow, a core density configuration has been&
                           & marked despite pub_ddec_use_coredens&
                           & = .false. ('//trim(error_string)//')')
                   end if
                   write(stdout,'(a)',advance='no') ' C'
                else
                   call utils_abort(' ERROR in internal_ref_dens_init:&
                        & Invalid flag detected in iions(it,2)')
                end if

                ! lpl: Print rcomp info read from the block and verify
                !      no core/valence mislabelling has occured
                if ( iions(it,2) == 0 ) then
                   if ( nint(ref_dens(isp)%rcomp(iion,6)) == 1 ) then
                      write(stdout,'(2x,I5,5(2x,f10.5))') iion, &
                              ref_dens(isp)%rcomp(iion,1), &
                              ref_dens(isp)%rcomp(iion,2), &
                              ref_dens(isp)%rcomp(iion,3), &
                              ref_dens(isp)%rcomp(iion,4), &
                              ref_dens(isp)%rcomp(iion,5)
                   else if ( nint(ref_dens(isp)%rcomp(iion,6)) == 2 ) then
                      write(stdout,'(I5,2x,a)') iion, &
                           trim(adjustl( ref_dens(isp)%refconf(iion) ))
                   else
                      write(error_string,'(a,1x,I3,1x,I2)') &
                           ref_dens(isp)%element%species_id, iion, &
                           nint(ref_dens(isp)%rcomp(iion,6))
                      call utils_abort(' ERROR in internal_ref_dens_init: &
                           & A valence density configuration has been &
                           & labelled as core ('//trim(error_string)//')')
                   end if
                else if ( iions(it,2) == 1 ) then
                   ! lpl: Not necessary to specify core density explicitly in
                   !      the rcomp block
                   if ( nint(ref_dens(isp)%rcomp(iion,6)) == 1 .or. &
                        nint(ref_dens(isp)%rcomp(iion,6)) == 3 ) then
                      write(stdout,'(2x,I5,5(2x,f10.5))') iion, &
                              ref_dens(isp)%rcomp(iion,1), &
                              ref_dens(isp)%rcomp(iion,2), &
                              ref_dens(isp)%rcomp(iion,3), &
                              ref_dens(isp)%rcomp(iion,4), &
                              ref_dens(isp)%rcomp(iion,5)
                   else if ( nint(ref_dens(isp)%rcomp(iion,6)) == 2 .or. &
                        nint(ref_dens(isp)%rcomp(iion,6)) ==4 ) then
                      write(stdout,'(I5,2x,a)') iion, &
                           trim(adjustl( ref_dens(isp)%refconf(iion) ))
                   end if
                end if

                ! lpl: Check that all densities are initialized in the same
                !      manner i.e. either all are read or solved. Otherwise,
                !      we can't apply the c3 reshaping as each resphaing
                !      depends on every other (valence) ion and they rely on
                !      having a consistent set of ionic densities
                if ( pub_ddec_c3_refdens ) then
                   ! If the first density was read from input (i.e. c2) then
                   ! make sure the current one is too (denoted by flag values
                   ! '2' for total (core + valence) and '4' for core densities.
                   ! Otherwise, spit out error as the current solver doees not
                   ! support c3
                   if ( nint(ref_dens(isp)%rcomp(iions(1,1),6)) == 2 ) then
                      ! Make sure current density is also being read from file
                      if ( nint(ref_dens(isp)%rcomp(iion,6)) /= 2 .and. &
                           nint(ref_dens(isp)%rcomp(iion,6)) /= 4 ) then
                           write(error_string,'(a,1x,I3,1x,I2)') &
                                   ref_dens(isp)%element%species_id, iion, &
                                   nint(ref_dens(isp)%rcomp(iion,6))
                           call utils_abort(' ERROR in internal_ref_dens_init:&
                                & Mixed density inputs detected ('//&
                                trim(error_string)//')')
                      end if
                   else
                      ! Abort if 'pub_ddec_c3_refdens = T' and solver density
                      ! was used anywhere
                      write(error_string,'(a,1x,I3,1x,I2)') &
                           ref_dens(isp)%element%species_id, iion, &
                           nint(ref_dens(isp)%rcomp(iion,6))
                      call utils_abort(' ERROR in internal_ref_dens_init:&
                           & Solver density is not compatible with c3&
                           & reshaping ('//trim(error_string)//')')
                   end if
                end if

             end do ! END do it=1,num_ions
             write(stdout,'(a)') '----------------------------&
                  &----------------------------------------'
          end if ! END if ( pub_on_root )
          call comms_bcast(pub_root_proc_id, pub_ddec_c3_refdens)

          ! lpl: Read or solve densities
          do it=1,num_ions

             iion = iions(it,1)

             ! lpl: Allocate reference density profiles on all procs
             if(.not. allocated(ref_dens(isp)%ion_dens(iion)%dens)) then
                ! lpl: Store allocation sequence
                ddec_ref_alloc_list%num = ddec_ref_alloc_list%num + 1
                if ( ddec_ref_alloc_list%num > &
                     size(ddec_ref_alloc_list%alloc,1)) &
                     call utils_abort(' ERROR: Failure in&
                          & internal_ref_dens_init')
                ddec_ref_alloc_list%alloc(ddec_ref_alloc_list%num,1) = isp
                ddec_ref_alloc_list%alloc(ddec_ref_alloc_list%num,2) = iion

                allocate(ref_dens(isp)%ion_dens(iion)%&
                     dens(ref_dens(isp)%npts), stat=ierr)
                call utils_alloc_check('internal_ref_dens_init', &
                     'ref_dens(isp)%ion_dens(iion)%dens',ierr)
             end if

             if(pub_debug_on_root) write(stdout,'(a,a4,1x,I4,1x,a,1x,e20.10)') &
                  ' INFO: ', ref_dens(isp)%element%symbol, iion, &
                  'ref_dens%dr =', ref_dens(isp)%dr

             ! lpl: COMMENT#01
             !      Options '1' and '3' are meant to call the internal solver
             !      to generate reference densities. These are no longer
             !      supported as we are using CHARGEMOL's reference densities
             !      read in from external files. The original code is still
             !      kept here for posterity, or if someone intends to use it for
             !      other purposes.

             ! lpl: COMMENT#01 - Refer to COMMENT#01a for index-keyword mapping
             ! lpl: Read or solve for valence densities
             if ( nint(ref_dens(isp)%rcomp(iion,6)) == 1 .or. &
                  nint(ref_dens(isp)%rcomp(iion,6)) == 3 ) then

                if ( pub_on_root ) write(stdout,'(a,/a,/a)') &
                     "ERROR - Solver-generated reference densities are&
                     & no longer supported. Please', 'supply reference&
                     & densities from external file (i.e. use option 'READ')"

                call utils_abort("ERROR - Solver-generated reference densities&
                     & no longer supported")

! lpl: ORIGINAL CODE BLOCK (COMMENT#01)----------------------------------------!
!
!                ! lpl: Solve using rcut = min(ngwf_cutoff,ref_dens(isp)%rcut)
!                !      as we want the reference density to decay to 0 at rcut
!                ! lpl: Solver works on pub_root_proc_id, bcasts to all procs
!                ! lpl: target radius is adjustable by user parameter
!                ! lpl: Additionally, if species-specific target radius is
!                !      specified then that will override the variable
!                !      'target_radius' so long as 'pub_ddec_target_radius' > 0
!                ! lpl: MOD#07 - allow incomplete ionic range by passing a
!                !      'ref_dens(isp)%init_status(iion)' flag that will be
!                !      checked later on
!                if(pub_ddec_target_radius <= 0.0_DP) then
!                   call ddec_solve_atom(ref_dens(isp),iion, &
!                        ref_dens(isp)%rcut, &
!                        ref_dens(isp)%element%nfunctions, pub_root_proc_id)
!                else
!                   call ddec_solve_atom(ref_dens(isp),iion, &
!                        pub_ddec_target_radius, &
!                        ref_dens(isp)%element%nfunctions, pub_root_proc_id)
!                end if
!
!                ! lpl: MOD#07 - 'ddec_solve_atom' will crash if solver fails
!                !      to initialize density (this is invoked in 'atom_solve'
!                !      and not 'ddec_solve_atom'), so
!                !      'ref_dens(isp)%init_status' will always be '.true.'
!                ref_dens(isp)%init_status(iion) = .true.
!
! lpl: ORIGINAL CODE BLOCK ----------------------------------------------------!

             else if ( nint(ref_dens(isp)%rcomp(iion,6)) == 2 .or. &
                  nint(ref_dens(isp)%rcomp(iion,6)) == 4) then

                  ! lpl: filename stored in 'ref_dens(isp)%refconf'
                  ! lpl: MOD#07 - allow incomplete ionic range
                  ! ndmh: added ddec_refdens_path option
                  if ( pub_on_root ) call ddec_read_atom(ref_dens(isp),iion, &
                       trim(adjustl(pub_ddec_refdens_path))// &
                       ref_dens(isp)%refconf(iion))
                  ! lpl: Unlike 'ddec_solve_atom' the subroutines
                  !      'ddec_read_atom' knows nothing about procs so bcast
                  !      whatever it returns to all procs
                  call comms_bcast(pub_root_proc_id, &
                       ref_dens(isp)%init_status(iion))
                  call comms_bcast(pub_root_proc_id, &
                       ref_dens(isp)%ion_dens(iion)%dens)

             else

                write(ref_dens(isp)%refconf(iion),'(a4,1x,I4)') &
                     ref_dens(isp)%element%species_id, iion
                call utils_abort(' ERROR in ddec_refdens_init: rcomp info&
                     & for '//trim(adjustl(ref_dens(isp)%refconf(iion)))//&
                     ' was not properly initialized')

             end if
             call comms_barrier
          end do ! END do it=1,num_ions

          ! lpl: MOD#07 - Sort out which densities have failed to initialize
          ! lpl: Check that neutral ion is at least initilaized
          if ( .not. ref_dens(isp)%init_status(0) ) &
               call utils_abort(' ERROR in internal_ref_dens_init: Neutral&
                    & refconf for species '//&
                    ref_dens(isp)%element%species_id//&
                    ' was not initialized')

          ! lpl: Check all ions
          do iion=min_ion_state(isp),max_ion_state(isp)
             if ( .not. ref_dens(isp)%init_status(iion) ) then
                if ( pub_on_root ) write(stdout,'(a,A4,sp,I4.3,a)') &
                       ' WARNING: refdens for ', &
                       ref_dens(isp)%element%species_id, iion, &
                       ' was not initialized'
             end if
          end do

          if ( include_coredens ) then
             iion = ref_dens(isp)%eff_nuc_charge
             if ( .not. ref_dens(isp)%init_status(iion) ) &
                  call utils_abort(' ERROR in internal_refdens_init: Core&
                       & density for '//ref_dens(isp)%element%species_id//&
                       ' was not initialized')
          end if

          call comms_barrier

          ! lpl: Modify density to 'c3' for all ions except the core density
          !      if required only if all species have been read not solved
          if ( pub_ddec_c3_refdens ) then
             ! lpl: Exclude 'max_ion_state' if it has the charge of the core
             !      i.e. same charge as 'eff_nuc_charge'
             if ( include_coredens .and. max_ion_state(isp) == &
                  ref_dens(isp)%eff_nuc_charge ) then
                call ddec_c3_dens(ref_dens(isp),min_ion_state(isp), &
                     max_ion_state(isp)-1)

                if ( ref_dens(isp)%c3flag(ref_dens(isp)%eff_nuc_charge) ) &
                     call utils_abort(' ERROR in internal_ref_dens_init: &
                     &c3 reshaping was applied to core density')
             else
                call ddec_c3_dens(ref_dens(isp),min_ion_state(isp), &
                     max_ion_state(isp))
             end if
          end if

          ! lpl: Point core_dens pointer to the core density (please make sure
          !      'internal_ref_dens_init' is called on ALL procs)
          if ( include_coredens ) ref_dens(isp)%core_dens => &
               ref_dens(isp)%ion_dens(ref_dens(isp)%eff_nuc_charge)%dens

          ! lpl: Reset static variables in 'ddec_read_atom'
          call ddec_read_atom(ref_dens(isp),0,'',.true.)

          ref_dens(isp)%initialized = .true.

          if ( pub_on_root ) write(stdout,'(a)') '==========================&
               &======================='

          call comms_barrier
       end do ! END do isp=1,par%num_species

       ! lpl: COMMENT#02 - These are meant to show whether the messy routine
       !      above initialized the reference densities correctly. Could be
       !      removed or flagged as a debug-only routine.
       ! lpl: Print out summary of initialized ions
       if ( pub_on_root ) then
          write(stdout,'(a)') '----------------------------------'
          write(stdout,'(a)') ' Reference density init summary'
          write(stdout,'(a)') '----------------------------------'
          write(stdout,'(6(1x,a))') &
               'Sp    ','Ion ','Init','Source','Type  ','c3'

          do isp=1,par%num_species

             ! lpl: Generate range of ionic states again
             num_ions = 0
             iions = -999
             do iion=min_ion_state(isp),max_ion_state(isp)
                num_ions = num_ions + 1
                iions(num_ions,1) = iion ! Stores ionic number
             end do

             ! lpl: Mark ion to be read from separate core density
             !      file if available, only if 'pub_ddec_use_coredens'
             !      is set to .true.
             if ( include_coredens ) then
                if ( iions(num_ions,1) /= ref_dens(isp)%eff_nuc_charge ) then
                   num_ions = num_ions + 1
                   iions(num_ions,1) = ref_dens(isp)%eff_nuc_charge
                end if
             end if

             ! lpl: Print info
             do it=1,num_ions

                iion = iions(it,1)

                if ( iion == min_ion_state(isp) ) then
                   write(stdout,'(1x,a4,1x,a)',advance='no') &
                        ref_dens(isp)%element%species_id,':'
                else
                   write(stdout,'(1x,6x)',advance='no')
                end if

                ! Ion charge
                write(stdout,'(1x,sp,I4.3)',advance='no') iion
                ! Init status
                write(stdout,'(3x,L1)',advance='no') &
                     ref_dens(isp)%init_status(iion)
                ! COMMENT#01a Init method
                select case ( nint(ref_dens(isp)%rcomp(iion,6)) )
                   case (1); write(stdout,'(1x,a6)',advance='no') 'SOLVE '
                   case (3); write(stdout,'(1x,a6)',advance='no') 'SOLVE '
                   case (2); write(stdout,'(1x,a6)',advance='no') 'READ  '
                   case (4); write(stdout,'(1x,a6)',advance='no') 'READ  '
                   case default
                      write(stdout,'(1x,a6)',advance='no') 'UNK   '
                end select
                ! Core flag
                select case ( nint(ref_dens(isp)%rcomp(iion,6)) )
                   case (1); write(stdout,'(1x,a6)',advance='no') 'NORMAL'
                   case (3); write(stdout,'(1x,a6)',advance='no') 'CORE  '
                   case (2); write(stdout,'(1x,a6)',advance='no') 'NORMAL'
                   case (4); write(stdout,'(1x,a6)',advance='no') 'CORE  '
                   case default
                      write(stdout,'(1x,a6)',advance='no') 'UNK   '
                end select
                ! c3 flag
                write(stdout,'(2x,L1)',advance='no') ref_dens(isp)%c3flag(iion)
                write(stdout,'(a)') ''
             end do
          end do ! END do isp=1,par%num_species

          write(stdout,'(a)') '----------------------------------'

       end if ! END if ( pub_on_root )
       call comms_barrier

       deallocate(max_ion_state,stat=ierr)
       call utils_dealloc_check('internal_ref_dens_init','max_ion_state',ierr)
       deallocate(min_ion_state,stat=ierr)
       call utils_dealloc_check('internal_ref_dens_init','min_ion_state',ierr)

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Leaving internal_ref_dens_init'

     end subroutine internal_ref_dens_init

  !===========================================================================!

     subroutine internal_ref_dens_destroy(ref_dens)

       !-----------------------------------------------------------------------!
       ! lpl: Destroy DDEC_REF_DENS                                            !
       !-----------------------------------------------------------------------!
       ! InOut : ref_dens - Array of reference densities for ALL atoms         !
       !-----------------------------------------------------------------------!
       implicit none

       ! lpl: Input
       type(DDEC_REF_DENS), pointer, intent(inout) :: ref_dens(:)

       ! lpl: Local variables
       integer :: isp, it, iion

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Entering internal_ref_dens_destroy'

       do isp=par%num_species,1,-1
          ! lpl: Destroy radial grid interpolation first
          do it=ref_dens(isp)%max_interp_shell,1,-1
             deallocate(ref_dens(isp)%rad_grid(it)%coord,stat=ierr)
             call utils_dealloc_check('internal_ref_dens_destroy', &
                  'ref_dens(isp)%rad_grid(it)%coord',ierr)
          end do
          ! lpl: Destroy radial grid
          deallocate(ref_dens(isp)%rad_grid,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
               'ref_dens(isp)%rad_grid',ierr)

          if ( allocated(ref_dens(isp)%neutral_orb_info) ) then
             deallocate(ref_dens(isp)%neutral_orb_info,stat=ierr)
             call utils_dealloc_check('internal_ref_dens_init', &
                  'ref_dens(isp)%neutral_orb_info',ierr)
          end if
          deallocate(ref_dens(isp)%rcomp,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
               'ref_dens(isp)%rcomp',ierr)
          deallocate(ref_dens(isp)%refconf,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
               'ref_dens(isp)%refconf',ierr)
       end do

       ! lpl: Deallocate all reference densities that have been allocated
       !      according to record stored in 'ddec_ref_alloc_list%alloc'
       do it=ddec_ref_alloc_list%num,1,-1
          isp  = ddec_ref_alloc_list%alloc(it,1)
          iion = ddec_ref_alloc_list%alloc(it,2)
          deallocate(ref_dens(isp)%ion_dens(iion)%dens,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
               'ref_dens(isp)%ion_dens(iion)%dens',ierr)
       end do

       do isp=par%num_species,1,-1
          ! lpl: Nullify pointer to core reference density
          nullify(ref_dens(isp)%core_dens)

          deallocate(ref_dens(isp)%c3flag,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
               'ref_dens(isp)%c3flag',ierr)

          ! lpl: MOD#07 - Additional flag 'ref_dens%init_status' to mark init
          !      status of reference densities for each ion
          deallocate(ref_dens(isp)%init_status,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
               'ref_dens(isp)%init_status',ierr)

          deallocate(ref_dens(isp)%ion_dens,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
               'ref_dens(isp)%ion_dens',ierr)

          ! lpl: Destroy hypothetical average radius array and
          !      radial array
          deallocate(ref_dens(isp)%rad_shell_vol,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
                  'ref_dens(isp)%rad_shell_vol',ierr)
          deallocate(ref_dens(isp)%avg_rad,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
                  'ref_dens(isp)%avg_rad',ierr)
          deallocate(ref_dens(isp)%rad,stat=ierr)
          call utils_dealloc_check('internal_ref_dens_destroy', &
               'ref_dens(isp)%rad',ierr)
       end do

       deallocate(ref_dens,stat=ierr)
       call utils_dealloc_check('internal_ref_dens_destroy', &
            '=>ref_dens',ierr)

       ! lpl: Destroy allocation list
       deallocate(ddec_ref_alloc_list%alloc,stat=ierr)
       call utils_dealloc_check('internal_ref_dens_destroy', &
            'ddec_ref_alloc_list%alloc',ierr)

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Leaving internal_ref_dens_destroy'

     end subroutine internal_ref_dens_destroy

  !===========================================================================!

     ! lpl: COMMENT#03 - This routine was haphazardly written to deal with
     !      the unplanned and ever-expanding keyword options in the reference
     !      density input block '%block ddec_rcomp' - Feel free to clean this
     !      up, change it, or remove it completely. Since the 'SOLVE' feature
     !      (c.f. COMMENT#01) should no longer be supported (not merely
     !      deprecated, unless there is a particular need to do so), keywords
     !      pertaining to that should be removed. Originally, the aforementioned
     !      block accepted only a string of reals delimited by colons, but due
     !      to the hasty addition of the 'READ' option, other routines were not
     !      updated to read in anything more than a bunch of reals. This parser
     !      simply translates string keywords in various formats into columns
     !      of reals for processing by several other routines (especially the
     !      'ddec_solve_atom' subroutine (COMMENT#01, COMMENT#04), which is no
     !      longer supported. ASDF
     subroutine internal_parse_rcomp(ref_dens,min_ion_state,max_ion_state)

       use esdf, only: esdf_reduce, llength
       use rundat, only: pub_ddec_rcomp
       use utils, only: utils_read_check

       implicit none

       ! lpl: InOut
       type(DDEC_REF_DENS), intent(inout) :: ref_dens

       ! lpl: Input
       integer, intent(in) :: min_ion_state, max_ion_state

       ! lpl: Local variablees
       character(len=128) :: rcomp_read_conf
       character(len=4) :: rcomp_read_sp
       integer :: rcomp_read_ion

       integer :: delim_pos(128)
       integer :: delim_num, ilist
       integer :: iion, num_ions, iions(300)
       logical :: all_species, read_core
       character(len=1), parameter :: tab = char(9)
       character(len=1) :: delim_set(4)
       data delim_set(:) / ' ', tab, ':', '"' /

       if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: Entering internal_parse_rcomp'

       ! lpl: Read rcomp block if 'pub_ddec_rcomp' exists
       all_species = .false.
       if( allocated(pub_ddec_rcomp) ) then
          ! lpl: Parse over 'pub_ddec_rcomp' strings
          do ilist=lbound(pub_ddec_rcomp,1),ubound(pub_ddec_rcomp,1)

             ! lpl: Skip blank lines
!             if ( len_trim(pub_ddec_rcomp(ilist)) == 0 ) cycle
                !djc commented out to allow DDEC calculation
                !djc during geom opt / MD

             ! lpl: Parse data - format in 'pub_ddec_rcomp(:)' should be
             !      <species symbol> <ion charge> <config string> e.g.
             !      C -2  1.0 : 2.0 : 3.0 : 4.0 : 5.0 : conf=2s2 2p4
             call ddec_parse_string(delim_pos,delim_num, &
                  pub_ddec_rcomp(ilist),delim_set(1:2))
             if ( delim_num < 3 ) call utils_abort(" ERROR in &
                  &internal_parse_rcomp: Insufficient &
                  &number of parameters in rcomp string '"//&
                  trim(adjustl(pub_ddec_rcomp(ilist)))//"'")

             ! lpl: 1st column is species label
             read(pub_ddec_rcomp(ilist)(delim_pos(1):delim_pos(2)-2),*, &
                  iostat=ierr) rcomp_read_sp
             call utils_read_check('internal_parse_rcomp','rcomp_read_sp',ierr)
             rcomp_read_sp = adjustl(rcomp_read_sp)

             ! lpl: 2nd column is ionic charge unless 'ALL' or
             !      'CORE' is specified
             read(pub_ddec_rcomp(ilist)(delim_pos(2):delim_pos(3)-2),*, &
                  iostat=ierr) rcomp_read_ion

             ! lpl: If we fail to read in an integer, this could be an
             !      'ALL' or 'CORE' flag
             read_core = .false.
             if ( ierr /= 0 ) then
                read(pub_ddec_rcomp(ilist)(delim_pos(2):delim_pos(3)-2),*, &
                     iostat=ierr) rcomp_read_conf
                call utils_read_check('internal_parse_rcomp', &
                     'rcomp_read_ion',ierr)
                rcomp_read_conf = adjustl(rcomp_read_conf)
                if ( len_trim(rcomp_read_conf) > llength ) &
                     call utils_abort(' ERROR in internal_parse_rcomp: &
                          &len_trim(rcomp_read_conf) > llength')
                rcomp_read_conf(1:llength) = &
                     esdf_reduce(rcomp_read_conf(1:llength))
                ! lpl: Check that this is an 'ALL' or 'CORE' flag
                select case (trim(adjustl(rcomp_read_conf)))
                   case ('core')
                      read_core = .true.
                      all_species = .false.
                   case ('all')
                      all_species = .true.
                   case default
                      call utils_abort(' ERROR in internal_parse_string:&
                           & Unrecognized character flag in 2nd column of&
                           & rcomp block')
                end select

                if ( read_core .and. .not. include_coredens ) then
                   write(stdout,'(a)') ' WARNING: pub_ddec_read_core_dens&
                        & = .false. but core density was provided. Core&
                        & densities will be ignored.'
                   ! lpl: Skip this block line
                   cycle
                end if
             end if ! END if ( ierr /= 0 )

             ! lpl: 3rd column onwards is rcomp config
             rcomp_read_conf = &
                  trim( adjustl(pub_ddec_rcomp(ilist)(delim_pos(3):)) )

             ! lpl: If the species_id of this ref_dens matches the species_id
             !      of the rcomp string
             if( trim(rcomp_read_sp) == trim(adjustl( &
                  ref_dens%element%species_id )) ) then
                num_ions = -999
                ! lpl: If 'ALL' or 'CORE' keyword is not present then we're
                !      reading info for only one ion
                if ( .not. all_species .and. .not. read_core ) then
                   ! lpl: Ignore ionic species that are out of range
                   if ( rcomp_read_ion < min_ion_state .and. &
                        rcomp_read_ion > max_ion_state ) then
                      write(stdout,'(a,a4,1x,I4,a)') &
                           ' WARNING: Ionic index ', &
                           adjustl(ref_dens%element%species_id), &
                           rcomp_read_ion,' is out of range and rcomp info&
                           & for this ion will be ignored'
                      num_ions = 0
                   ! lpl: Else we are dealing only with one ion
                   else
                      num_ions = 1
                      iions(num_ions) = rcomp_read_ion
                   end if
                ! lpl: If this rcomp string instructs us to read all species
                else if ( all_species ) then
                   ! lpl: else provide a list of ions to be read
                   num_ions = 0
                   do it=min_ion_state,max_ion_state
                      num_ions = num_ions + 1
                      iions(num_ions) = it
                   end do
                   ! lpl: For the identified 'all species' flag, if
                   !      calculation requires core densities then add
                   !      one more ion to be read. If this fails then
                   !      no error will occur as core densities will be
                   !      checked later on
                   if ( include_coredens ) then
                      if ( iions(num_ions) /= ref_dens%eff_nuc_charge ) then
                         num_ions = num_ions + 1
                         iions(num_ions) = ref_dens%eff_nuc_charge
                      end if
                   end if
                ! lpl: If this rcomp string instructs us to read core density
                else if ( read_core ) then
                   num_ions = 1
                   iions(num_ions) = ref_dens%eff_nuc_charge
                end if

                ! lpl: Read rcomp config for ion of charge 'rcomp_read_ion'
                do iion=1,num_ions
                   ! lpl: Obtain ionic number from list
                   rcomp_read_ion = iions(iion)

                if (pub_debug_on_root) write(stdout,'(a,2(2x,I5),2(2x,L1),2x,I4)') &
                     ' READING: '//ref_dens%element%species_id, &
                     rcomp_read_ion, nint(ref_dens%rcomp(rcomp_read_ion,6)), &
                     read_core, all_species, &
                     nint(ref_dens%rcomp(rcomp_read_ion,6))

                   ! lpl: Check that this ionic sp has not been initialized
                   if( nint(ref_dens%rcomp(rcomp_read_ion,6)) > 0 ) then
                      ! lpl: If this has already been initialized previously
                      !      as a regular cationic density, overwrite config
                      !      with the core density specification
                      if ( read_core ) then
                         if ( nint(ref_dens%rcomp(rcomp_read_ion,6)) == 3 .or. &
                              nint(ref_dens%rcomp(rcomp_read_ion,6)) == 4 ) then
                            call utils_abort(' ERROR in internal_parse_rcomp: &
                                 &Multiple core density configurations &
                                 &declared for '//ref_dens%element%species_id)
!djc                         else
!                            write(stdout,'(a,/a,/a)') ' WARNING: Overlapping&
!                                 & rcomp specifications for', ' core and&
!                                 & cationic density for '//&
!                                 trim(ref_dens%element%species_id), &
!djc                                 ' The core specification will take precedence.'
                         end if
                         ! lpl: Set flag to uninitialized
                         ref_dens%rcomp(rcomp_read_ion,6) = -1
                      ! lpl: If we are not reading a core density
                      else
                         ! lpl: If this ionic species has already been init'd
                         !      as a core density (marked as '3' or '4'),
                         !      then skip parsing info for the current ion
                         if ( nint(ref_dens%rcomp(rcomp_read_ion,6)) >= 3 ) then
                            write(stdout,'(a,/a,/a)') ' WARNING: A cation&
                                 & configuration with the same', ' ionic charge&
                                 & as the core density of '//&
                                 trim(ref_dens%element%species_id), &
                                 ' has been found and will be ignored.'
                            cycle
                         ! lpl: Ohterwise abort with error
                         else
                            if (pub_debug_on_root) write(stdout, &
                                 '(a,2(2x,I5),2(2x,L1),2x,I4)') &
                                 'ERROR: '//ref_dens%element%species_id, &
                                 rcomp_read_ion, &
                                 nint(ref_dens%rcomp(rcomp_read_ion,6)), &
                                 read_core, all_species, &
                                 nint(ref_dens%rcomp(rcomp_read_ion,6))
                            call utils_abort(' ERROR: rcomp configuration&
                                 & for '//ref_dens%element%species_id//' has&
                                 & been multiply defined')
                         end if
                      end if
                   end if ! END if( nint(ref_dens%rcomp(rcomp_read_ion,6)) > 0 )

                   ! lpl: Parse the rcomp configuration string
                   delim_pos = -1
                   call ddec_parse_string(delim_pos,delim_num, &
                        rcomp_read_conf,delim_set(3:3))
                   ! lpl: By default all parameters must be -1
                   ref_dens%rcomp(rcomp_read_ion,:) = -1.0_DP

                   ! lpl: Read 1st parameter since this can be standalone, and
                   !     check whether we can read it as a real
                   read(rcomp_read_conf(delim_pos(1):delim_pos(2)-2),*, &
                           iostat=ierr) ref_dens%rcomp(rcomp_read_ion,1)
                   ! lpl: If we are unable to successfully read it as a real,
                   !      try reading the 1st field only as a string
                   if ( ierr /= 0 ) then
                      ! lpl: Read 1st field of the string - behavior might be
                      !      compiler-dependent since no format is specified?
                      rcomp_read_conf = adjustl(rcomp_read_conf)
                      read(rcomp_read_conf(1:len_trim(rcomp_read_conf)),*, &
                           iostat=ierr) ref_dens%refconf(rcomp_read_ion)
                      call utils_read_check('internal_parse_rcomp', &
                           'rcomp_read_conf #1',ierr)

                      ! lpl: Set the value of the extra 6th column of the
                      !      rcomp data array to 2 to indicate that this ion
                      !      is to be read not solved, and this ionic species
                      !      has been initialized
                      ref_dens%rcomp(rcomp_read_ion,6) = 2
                   else
                      ! lpl: Else if 1st parameter can be read as a real, read
                      !      columns 2 to 5 for real parameters
                      ! lpl: MOD#13 - fixed erronous limit 'delim_num'
                      do it=2,min(5,delim_num)
                         read(rcomp_read_conf(delim_pos(it):&
                              delim_pos(it+1)-2),*,iostat=ierr) &
                              ref_dens%rcomp(rcomp_read_ion,it)
                         call utils_read_check('internal_parse_rcomp', &
                              'rcomp_read_conf #2',ierr)
                      end do
                      ! lpl: Read 6th column for atomic config string and
                      !      ignore additional columns
                      ! lpl: MOD#09 - do not use quotes for config string
                      if ( delim_num >= 6 ) then
                         read(rcomp_read_conf(delim_pos(6):delim_pos(7)-2), &
                             '(a)',iostat=ierr) &
                             ref_dens%refconf(rcomp_read_ion)
                         call utils_read_check('internal_parse_rcomp', &
                              'rcomp_read_conf #3',ierr)
                         ref_dens%refconf(rcomp_read_ion) = &
                              adjustl(ref_dens%refconf(rcomp_read_ion))

                         ! lpl: MOD#09
                         if (pub_debug_on_root) write(stdout,'(a)') &
                              ' READ: solver conf string = "'//&
                              trim(ref_dens%refconf(rcomp_read_ion))//'"'

                      end if

                      ! lpl: Set the value of the extra 6th column of the rcomp
                      !      data array 'ref_dens%rcomp(:)' (not the 6th
                      !      column of the rcomp configuration portion of the
                      !      input string 'rcomp_read_conf') to 1 to indicate
                      !      that rcomp info for a regular ion has been
                      !      initialized for this species
                      ref_dens%rcomp(rcomp_read_ion,6) = 1
                   end if

                   ! lpl: If the current string instructs us to read or solve
                   !      for a core density, then the 6th column of rcomp
                   !      will be '3' for solve and '4' for read
                   if ( read_core ) ref_dens%rcomp(rcomp_read_ion,6) = &
                        ref_dens%rcomp(rcomp_read_ion,6) + 2
                end do ! END do rcomp_read_ion=min_ion_read,max_ion_read

                ! lpl: Empty string so we won't come back to it later
!                pub_ddec_rcomp(ilist) = ''
                !djc commented out to allow DDEC calculation
                !djc during geom opt / MD

             end if ! END if matching species_id found

          end do ! END loop over all pub_ddec_rcomp(it)

       end if ! END if( allocated(pub_ddec_rcomp) )

       if (pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: Leaving internal_parse_rcomp'

     end subroutine internal_parse_rcomp

  !===========================================================================!

     subroutine internal_update_ih_dens(interp_dens,ref_dens,atom_popn)

       !-----------------------------------------------------------------------!
       ! lpl: Linear interpolation between IH reference densities              !
       !-----------------------------------------------------------------------!
       ! Output: interp_dens - IH reference density interpolated  between      !
       !            densities of ions from floor/ceiling(atom_popn)            !
       ! Input : ref_dens - Reference density for one atom. Contains all       !
       !            generated ionic configurations                             !
       !       : atom_popn - Atomic population at which the IH density should  !
       !            be interpolated                                            !
       !-----------------------------------------------------------------------!

       ! lpl: Output
       real(kind=DP), intent(out) :: interp_dens(:)

       ! lpl: Input
       type(DDEC_REF_DENS), intent(in) ::  ref_dens
       real(kind=DP), intent(in) :: atom_popn

       integer :: ion_lower, ion_upper
       real(kind=DP) :: atom_charge

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Entering internal_update_ih_dens'

       ! lpl: MOD#05 - Compute atom_charge from 'ddec_nuclear_charge'
       !      instead of assuming 'atom_popn' contains val + core
!       atom_charge = ref_dens%element%atomic_number - atom_popn
       atom_charge = ref_dens%ddec_nuc_charge - atom_popn

       ion_lower = floor(atom_charge)
       ion_upper = ceiling(atom_charge)

       if(ion_upper > ref_dens%max_cation) then
          write(stdout,'(a)') &
               ' ERROR: ion_upper out of range in internal_update_ih_dens'
          write(stdout,'(a3,a,f14.7,a,I3)') ref_dens%element%symbol, &
               '  Q = ',atom_charge,'  ion_upper = ',ion_upper
          call utils_abort('')
       else if(ion_lower < ref_dens%max_anion) then
          write(stdout,'(a)') &
               ' ERROR: ion_lower out of range in internal_update_ih_dens'
          write(stdout,'(a3,a,f14.7,a,I3)') ref_dens%element%symbol, &
               '  Q = ',atom_charge,'  ion_upper = ',ion_lower
          call utils_abort('')
       end if

       ! lpl: Return error if lower ionic reference state does not exist.
       !      Possible expansion later to allow automatic generation of new
       !      densities on-the-fly as required
       ! lpl: MOD#07 - Changed checking allocation to checking init_status flag
       if ( .not. ref_dens%init_status(ion_lower) ) then
!       if( .not. allocated(ref_dens%ion_dens(ion_lower)%dens) ) then
          write(stdout,'(a)') &
               ' ERROR: ion_lower out of range in internal_update_ih_dens'
          write(stdout,'(a3,a,f14.7,a,I3)') ref_dens%element%symbol, &
               '  Q = ',atom_charge,'  ion_lower = ',ion_lower
          call utils_abort('')
       ! else ... allocate and initialize the required ionic reference density
       end if

       ! lpl: Return error if upper ionic reference state does not exist.
       ! lpl: MOD#07 - Changed checking allocation to checking init_status flag
       if ( .not. ref_dens%init_status(ion_upper) ) then
!       if( .not. allocated(ref_dens%ion_dens(ion_upper)%dens) ) then
          write(stdout,'(a)') &
               ' ERROR: ion_upper out of range in internal_update_ih_dens'
          write(stdout,'(a3,a,f14.7,a,I3)') ref_dens%element%symbol, &
               '  Q = ',atom_charge,'  ion_upper = ',ion_upper
          call utils_abort('')
       ! else ... allocate and initialize the required ionic reference density
       end if

       ! lpl: Linear-interpolate to get reference fractional density
       interp_dens = 0.0_DP
       interp_dens(1:ref_dens%npts) &
            = (real(ion_upper,kind=DP) - atom_charge)&
                 *ref_dens%ion_dens(ion_lower)%dens &
            + (atom_charge - real(ion_lower,kind=DP))&
                 *ref_dens%ion_dens(ion_upper)%dens

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Leaving internal_update_ih_dens'

     end subroutine internal_update_ih_dens

  !===========================================================================!

  !----------------------------------------------------------------------------!
  ! lpl: Miscellaneous routines                                                !
  !----------------------------------------------------------------------------!

     subroutine internal_multipole(dipole, quadrupole, aim_val, &
          local_box, rad_dens)

       implicit none

       ! Output
       real(kind=DP), intent(out) :: dipole(3)
       real(kind=DP), intent(out) :: quadrupole(5)

       ! Input
       type(DDEC_LOCAL_BOX), intent(in) :: local_box
       type(DDEC_RAD_DENS), intent(in) :: rad_dens
       real(kind=DP), intent(in) :: aim_val(:,:,:)

       ! Local variables
       type(POINT) :: local_atom_ctr, local_grid_coord, r_coord
       type(POINT) :: local_dipole
       real(kind=DP) :: local_quadrupole(5)
       integer :: iat, it1, it2, it3, ir

       if (pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Entering internal_multipole'
       if (pub_debug) call timer_clock('ddec_internal_multipole',1)

       if ( size(aim_val,1) /= local_box%n1 .or. &
            size(aim_val,2) /= local_box%n2 .or. &
            size(aim_val,3) /= local_box%n3 ) &
            call utils_abort(' ERROR in internal_multipole:&
                 & size(aim_val) /= local_box')

       iat = rad_dens%iat

       local_atom_ctr = mdl%elements(par%orig_atom(iat))%centre &
            - local_box%box_coord_origin(iat)

       local_dipole%x = 0.0_DP
       local_dipole%y = 0.0_DP
       local_dipole%z = 0.0_DP
       local_quadrupole = 0.0_DP

       do it1=1,local_box%n1
          do it2=1,local_box%n2
             do it3=1,local_box%n3

                ! lpl: Coord w.r.t. voxel 1,1,1 centre
                local_grid_coord &
                     = real(it1 - 1,kind=DP)*grid%da1 &
                     + real(it2 - 1,kind=DP)*grid%da2 &
                     + real(it3 - 1,kind=DP)*grid%da3

                ! lpl: r_coord w.r.t. atom center
                r_coord = local_grid_coord - local_atom_ctr

                ! lpl: Atom-centred radius index
                ir = aint(geometry_magnitude(r_coord)/rad_dens%dr) + 1

                if ( ir > rad_dens%npts ) cycle

                if ( aim_val(it1,it2,it3) < 0.0_DP ) &
                     call utils_abort(' ERROR in internal_multipole:&
                          & aim_val < 0')

                ! Sum dipole
                local_dipole = local_dipole + aim_val(it1,it2,it3)*r_coord

                ! Sum quadrupole
                ! Q_xy
                local_quadrupole(1) = local_quadrupole(1) &
                   + r_coord%x*r_coord%y*aim_val(it1,it2,it3)
                ! Q_xz
                local_quadrupole(2) = local_quadrupole(2) &
                   + r_coord%x*r_coord%z*aim_val(it1,it2,it3)
                ! Q_yz
                local_quadrupole(3) = local_quadrupole(3) &
                   + r_coord%y*r_coord%z*aim_val(it1,it2,it3)
                ! Q_x2-y2
                local_quadrupole(4) = local_quadrupole(4) &
                   + ( (r_coord%x)**2 - (r_coord%y)**2 )*aim_val(it1,it2,it3)
                ! Q_3z2-r2
                local_quadrupole(5) = local_quadrupole(5) &
                   + ( 3.0_DP*(r_coord%z)**2 &
                   - (r_coord.DOT.r_coord) )*aim_val(it1,it2,it3)

             end do
          end do
       end do

       dipole(1) = local_dipole%x
       dipole(2) = local_dipole%y
       dipole(3) = local_dipole%z
       quadrupole(:) = local_quadrupole(:)

       if (pub_debug) call timer_clock('ddec_internal_multipole',2)
       if ( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_multipole'

     end subroutine internal_multipole

  !===========================================================================!

     subroutine internal_moment(moment, aim_val, local_box, &
          rad_dens, moment_order)

       implicit none

       ! Output
       real(kind=DP), intent(out) :: moment

       ! Input
       type(DDEC_LOCAL_BOX), intent(in) :: local_box
       type(DDEC_RAD_DENS), intent(in) :: rad_dens
       real(kind=DP), intent(in) :: aim_val(:,:,:)
       integer, intent(in) :: moment_order

       ! Local variables
       type(POINT) :: local_atom_ctr, local_grid_coord
       real(kind=DP) :: local_moment, r_dist
       integer :: iat, it1, it2, it3, ir

       if (pub_debug_on_root)  write(stdout,'(a)') &
            ' DEBUG: Entering internal_moment'
       if (pub_debug) call timer_clock('ddec_internal_moment',1)

       if ( size(aim_val,1) /= local_box%n1 .or. &
            size(aim_val,2) /= local_box%n2 .or. &
            size(aim_val,3) /= local_box%n3 ) &
            call utils_abort(' ERROR in internal_moment:&
                 & size(aim_val) /= local_box')

       iat = rad_dens%iat

       local_atom_ctr = mdl%elements(par%orig_atom(iat))%centre &
            - local_box%box_coord_origin(iat)

       local_moment = 0.0_DP

       do it1=1,local_box%n1
          do it2=1,local_box%n2
             do it3=1,local_box%n3

                ! lpl: Coord w.r.t. voxel 1,1,1 centre
                local_grid_coord &
                     = real(it1 - 1,kind=DP)*grid%da1 &
                     + real(it2 - 1,kind=DP)*grid%da2 &
                     + real(it3 - 1,kind=DP)*grid%da3

                ! lpl: r_coord w.r.t. atom center
                r_dist = geometry_magnitude(local_atom_ctr - local_grid_coord)

                ! lpl: Atom-centred radius index
                ir = aint(r_dist/rad_dens%dr) + 1

                if ( ir > rad_dens%npts ) cycle

                if ( aim_val(it1,it2,it3) < 0.0_DP ) &
                     call utils_abort(' ERROR in internal_moment: aim_val < 0')

                local_moment = local_moment + &
                     (r_dist**moment_order)*aim_val(it1,it2,it3)

             end do
          end do
       end do

       moment = local_moment

       if (pub_debug) call timer_clock('ddec_internal_moment',2)
       if ( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: Leaving internal_moment'

     end subroutine internal_moment

  !===========================================================================!

     subroutine internal_eff_decay_exp(eff_decay_a, eff_decay_b, &
          eff_decay_rsquared, rad_dens_val, rad_dens, eff_decay_rmin)

       !-----------------------------------------------------------------------!
       ! lpl: Least-squareslog fit of f(r) = exp(a-b*r) to AIM density profile !
       !      n_A(r). Copypaste of CHARGEMOL's routine                         !
       !-----------------------------------------------------------------------!
       ! Output : eff_decay_a, eff_decay_b - Fitted parameters of f(r)         !
       ! Input  : rad_dens_val - Radial density profile to be fitted           !
       !        : rad_dens - DDEC_RAD_DENS object associated with              !
       !             'rad_dens_val'                                            !
       !        : eff_decay_rmin - Minimum r for wich n_A(r) is fitted to f(r) !
       !-----------------------------------------------------------------------!

       implicit none

       ! Output
       real(kind=DP), intent(out) :: eff_decay_a, eff_decay_b
       real(kind=DP), intent(out) :: eff_decay_rsquared

       ! Input
       real(kind=DP), intent(in) :: rad_dens_val(:)
       type(DDEC_RAD_DENS), intent(in) :: rad_dens

       real(kind=DP), intent(in) :: eff_decay_rmin

       ! Local variables
       real(kind=DP) :: f, r
       real(kind=DP) :: Sr, Srr, Sfr, Sf, Sff
       integer :: ir, irmin, num

       if( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: ddec_main: Entering internal_eff_decay_exp'

       if ( size(rad_dens_val) /= rad_dens%npts ) call utils_abort(' ERROR in&
            & internal_eff_decay_exp: size(rad_dens_val) /= rad_dens%npts')
       if ( eff_decay_rmin >= rad_dens%rcut ) call utils_abort(' ERROR in&
            & internal_eff_decay_exp: eff_decay_rmin >= rad_dens%rcut')

       iat = rad_dens%iat
       irmin = int( floor(eff_decay_rmin/rad_dens%dr) ) + 1

       Sr  = 0.0_DP
       Srr = 0.0_DP
       Sfr = 0.0_DP
       Sf  = 0.0_DP
       Sff = 0.0_DP
       num = 0

       do ir=irmin,rad_dens%npts
          ! lpl: Only fit if n_A(r) > threhsold. Follows CHARGEMOL
          !      by defining threhsold as zero_threshold^(3/2)
          if ( rad_dens_val(ir) >= pub_ddec_zero_thresh**1.5_DP ) then
             f = log(rad_dens_val(ir))
             r = rad_dens%rad(ir)

             Sr  = Sr + r
             Srr = Srr + r**2
             Sfr = Sfr + f*r
             Sf  = Sf + f
             Sff = Sff + f**2
             num = num + 1
          end if
       end do

       eff_decay_b = (num*Sfr - Sf*Sr)/(num*Srr - Sr*Sr)
       eff_decay_a = (Sf - eff_decay_b*Sr)/num
       eff_decay_rsquared = (num*Sfr - Sf*Sr)**2/&
            ((num*Srr - Sr*Sr)*(num*Sff - Sf*Sf))

       if( pub_debug_on_root ) write(stdout,'(a)') &
            ' DEBUG: ddec_main: Leaving internal_eff_decay_exp'

     end subroutine internal_eff_decay_exp

  !===========================================================================!

     subroutine internal_interp_rad_dens(interp_part_dens, orig_part_dens, &
          rad_dens, ref_dens, local_box_val, local_box)

       !-----------------------------------------------------------------------!
       ! lpl: Transfers density from local_box to radial grid for a specified  !
       !      number of shells from the atomic center by trilinear             !
       !      interpolation                                                    !
       !-----------------------------------------------------------------------!
       ! InOut : rad_dens - Radial density for one atom. The interpolated      !
       !            density will overwrite 'rad_dens%new_part_dens'. Also      !
       !            really an 'InOut'                                          !
       ! Input : ref_dens - Array of reference densities for ALL species.      !
       !            Contains the necessary information regarding the shells    !
       !            requiring interpolation                                    !
       !       : local_box_fine - Should already contain the deposited         !
       !            spherically-averaged AIM density in                        !
       !            'local_box_fine%total_pseudoval'                           !
       !-----------------------------------------------------------------------!

       implicit none

       ! Output
       real(kind=DP), intent(out) :: interp_part_dens(:)
       ! Input
       type(DDEC_RAD_DENS), intent(in) :: rad_dens
       type(DDEC_REF_DENS), intent(in) :: ref_dens(:)
       type(DDEC_LOCAL_BOX), intent(in) :: local_box
       real(kind=DP), intent(in) :: orig_part_dens(:)

       ! Buffer
       real(kind=DP), intent(inout) :: local_box_val(:,:,:)

       ! lpl: Local variables
       real(kind=DP) :: interp_dens
       type(POINT) :: local_atom_ctr, local_grid_coord
       integer :: r_index, it, iat, max_interp_shell

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Entering internal_interp_rad_dens'

       if ( size(local_box_val,1) /= local_box%n1 .or. &
            size(local_box_val,2) /= local_box%n2 .or. &
            size(local_box_val,3) /= local_box%n3 ) &
            call utils_abort(' ERROR in internal_interp_rad_dens:&
                 & size(local_box_val) /= local_box')

       if ( size(interp_part_dens) /= rad_dens%npts ) &
            call utils_abort(' ERROR in internal_interp_rad_dens:&
                 & size(interp_part_dens) /= rad_dens%npts')

       if ( size(orig_part_dens) /= rad_dens%npts ) &
            call utils_abort(' ERROR in internal_interp_rad_dens:&
                 & size(orig_part_dens) /= rad_dens%npts')

       iat = rad_dens%iat

       max_interp_shell = ref_dens( mdl%elements(par%orig_atom(iat) )% &
            species_number )%max_interp_shell

       ! lpl: Reset range of radial densities requiring interpolation
       interp_part_dens(1:max_interp_shell) = 0.0_DP

       ! lpl: Grid coord w.r.t. voxel 1,1,1 centre
       local_atom_ctr = mdl%elements(par%orig_atom(iat))%centre &
            - local_box%box_coord_origin(iat)

       ! lpl: Iterate over set of shells (of small radii) that
       !      have been marked as requring density interpolation
       !      (first shell up to index 'max_interp_shell' for each
       !      species) instead of the usual 'grid point binning'
       do r_index=1,max_interp_shell
          ! lpl: Iterate over the list of surface coordinates stored
          !      in 'rad_grid' for this radial shell 'r_index'
          do it=1,ref_dens( mdl%elements(par%orig_atom(iat) )%species_number &
               )%rad_grid(r_index)%npts
             ! lpl: Local Cartesian coordinate w.r.t. local_box origin
             local_grid_coord = ref_dens( mdl%elements(par%orig_atom(iat) )%&
                  species_number )%rad_grid(r_index)%coord(it) &
                  + local_atom_ctr
             ! lpl: Obtain value of densities at this surface coordinate via
             !      linear interpolation
             interp_dens = internal_box_interp(local_grid_coord, &
                  local_box_val)
             if ( interp_dens < 0.0_DP ) call utils_abort(' ERROR in &
                  &internal_interp_rad_dens - density < 0')

             interp_part_dens(r_index) = interp_dens + &
                  interp_part_dens(r_index)
          end do
       end do

       ! lpl: Get spherically-averaged densities for all radial grid values
       !      and count atomic population simultaneously
       do r_index=1,max_interp_shell
          if ( ref_dens( mdl%elements(par%orig_atom(iat) )%species_number &
               )%rad_grid(r_index)%npts > 0 ) &
               interp_part_dens(r_index) = &
                    interp_part_dens(r_index)/&
                    ( 1.0_DP*ref_dens( mdl%elements(par%orig_atom(iat) )&
                    %species_number )%rad_grid(r_index)%npts )
       end do

       ! lpl: Copy the uninterpolated densities into 'interp_part_dens' in the
       !      same way as how 'interp_part_dens' is updated after each
       !      iteration (Mostly done so 'internal_update_rad_dens' can be
       !      called straightforwardly later on)
       interp_part_dens(max_interp_shell+1:rad_dens%npts) &
            = orig_part_dens(max_interp_shell+1:rad_dens%npts)

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Leaving internal_interp_rad_dens'

     end subroutine internal_interp_rad_dens

  !===========================================================================!

     subroutine internal_radial_grid(rad_grid,sphere_rad)

       !-----------------------------------------------------------------------!
       ! lpl: Use golden spiral section to generate radial pts on sphere of    !
       !      radius 'sphere_rad'                                              !
       !-----------------------------------------------------------------------!
       ! InOut : rad_grid - Object which is just a list of arrays containing   !
       !            coordinates of the grid points on the surface of a sphere  !
       !            of radius 'sphere_rad'                                     !
       ! Input : sphere_rad - Radius of sphere                                 !
       !-----------------------------------------------------------------------!

       implicit none

       ! lpl: InOut
       type(RADIAL_GRID), intent(inout) :: rad_grid

       ! lpl: Input
       real(kind=DP), intent(in) :: sphere_rad

       ! lpl: Local variables
       real(kind=DP), parameter :: alpha = PI*( 3.0_DP - sqrt(5.0_DP) )
       real(kind=DP) :: dy, y, phi, r
       integer :: it, npts

       npts = rad_grid%npts
       dy = 1.0_DP/real(npts,kind=DP)

       do it=1,npts
          y = dy*real(2*it - 1,kind=DP) - 1.0_DP
          phi = alpha*it
          r = sqrt( 1.0_DP - y**2 )

          rad_grid%coord(it)%x = sphere_rad*r*cos(phi)
          rad_grid%coord(it)%y = sphere_rad*y
          rad_grid%coord(it)%z = sphere_rad*r*sin(phi)
       end do

     end subroutine internal_radial_grid

  !===========================================================================!

     real(kind=DP) function internal_box_interp(local_r, local_box_val)

       !-----------------------------------------------------------------------!
       ! lpl: Triinear interpolates value at 'local_r' based on grid calues in !
       !      'local_box_val'. 'local_r' is the Cartesian coordinates relative !
       !      to the origin (defined as index (1,1,1)) of 'local_box_val'.     !
       !      'local_box_val' is any generic 3D array (not necessarily a       !
       !      DDEC_LOCAL_BOX array) assumed to share the same GRID_INFO as the !
       !      supercell.                                                       !
       !-----------------------------------------------------------------------!
       ! Input : local_r - Coordinate requiring interpolation w.r.t. origin of !
       !            'local_box_val'.                                           !
       !       : local_box_val - 3D array. Origin ASSUMED to be 0,0,0 for      !
       !            index (1,1,1) and having the SAME lattice parameters as    !
       !            the simulation cell.                                       !
       !-----------------------------------------------------------------------!

       implicit none

       ! lpl: Inputs
       type(POINT), intent(in) :: local_r
       real(kind=DP), intent(in) :: local_box_val(:,:,:) ! 3D array

       ! lpl: Local variables
       integer :: n_000_1, n_000_2, n_000_3
       real(kind=DP) :: f_000_001, f_100_101, f_010_011, f_110_111 ! z
       real(kind=DP) :: f_00_01, f_10_11                           ! y
       real(kind=DP) :: frac_1, frac_2, frac_3

       ! lpl: Make sure 'local_box_val' is following convention
       if( (lbound(local_box_val,1) <= 0) .or. &
           (lbound(local_box_val,2) <= 0) .or. &
           (lbound(local_box_val,3) <= 0) ) call utils_abort(' ERROR in &
           &internal_box_interp: local_box_val lbound index is not (1,1,1)')

       ! lpl: Get array indices n_000_(1,2,3) where this 'local_r' is located.
       frac_1 = (local_r.DOT.grid_ua1)/grid_lda1
       frac_2 = (local_r.DOT.grid_ua2)/grid_lda2
       frac_3 = (local_r.DOT.grid_ua3)/grid_lda3

       ! lpl: lbound of 'local_box_val' is (1,1,1), not (0,0,0)
       !      This should be floor of real, since for interpolation, we require
       !      the lower and upper bound of the frac, and the upperbound is
       !      simply 'lower+1'
       n_000_1 = aint(frac_1) + 1
       n_000_2 = aint(frac_2) + 1
       n_000_3 = aint(frac_3) + 1

       ! lpl: Make sure 'local_r' is inside 'local_box_val'
       if ( (n_000_1 + 1 > size(local_box_val,1)) .or. &
            (n_000_2 + 1 > size(local_box_val,2)) .or. &
            (n_000_3 + 1 > size(local_box_val,3)) .or. &
            (n_000_1 < 1) .or. (n_000_2 < 1) .or. (n_000_3 < 1) ) &
            call utils_abort(' ERROR in internal_box_interp: &
               &local_r out of bounds')

       ! lpl: Fractional coord of local_r w.r.t. origin of the voxel
       frac_1 = frac_1 - real(n_000_1 - 1,kind=DP)
       frac_2 = frac_2 - real(n_000_2 - 1,kind=DP)
       frac_3 = frac_3 - real(n_000_3 - 1,kind=DP)

       if( frac_1 < 0.0_DP .or. frac_2 < 0.0_DP .or. frac_3 < 0.0_DP ) &
            call utils_abort(' ERROR in internal_box_interp: Fractional &
               &coordinate w.r.t. voxel < 0')

       ! lpl: Linear interpolate z, then y, then x
       f_000_001 = services_linear_interpolation(frac_3, &
            local_box_val(n_000_1,n_000_2,n_000_3), &
            local_box_val(n_000_1,n_000_2,n_000_3+1), &
            0.0_DP,1.0_DP)

       f_100_101 = services_linear_interpolation(frac_3, &
            local_box_val(n_000_1+1,n_000_2,n_000_3), &
            local_box_val(n_000_1+1,n_000_2,n_000_3+1), &
            0.0_DP,1.0_DP)

       f_010_011 = services_linear_interpolation(frac_3, &
            local_box_val(n_000_1,n_000_2+1,n_000_3), &
            local_box_val(n_000_1,n_000_2+1,n_000_3+1), &
            0.0_DP,1.0_DP)

       f_110_111 = services_linear_interpolation(frac_3, &
            local_box_val(n_000_1+1,n_000_2+1,n_000_3), &
            local_box_val(n_000_1+1,n_000_2+1,n_000_3+1), &
            0.0_DP,1.0_DP)

       f_00_01 = services_linear_interpolation(frac_2, &
            f_000_001,f_010_011,0.0_DP,1.0_DP)

       f_10_11 = services_linear_interpolation(frac_2, &
            f_100_101,f_110_111,0.0_DP,1.0_DP)

       internal_box_interp = services_linear_interpolation(frac_1, &
            f_00_01,f_10_11,0.0_DP,1.0_DP)

     end function internal_box_interp

  end subroutine ddec_main

  !===========================================================================!

  ! lpl: COMMENT#04 - This routine should no longer be used (c.f. COMMENT#01)

!  subroutine ddec_solve_atom(ref_dens,ionic_charge, &
!       input_target_radius,target_nfuncs,proc)
!
!    !------------------------------------------------------------------!
!    ! This routine solves the Schrodinger equation for a single ion,   !
!    ! and creates the corresponding starting radial functions for the  !
!    ! NGWFs.                                                           !
!    ! lpl: Plaigiarized from ngwfs_solve_atom and removed NGWF         !
!    !      generation to leave only density as return                  !
!    ! lpl: Input - 'current_element' with 'ionic_charge' solved using  !
!    !              Bessel basis of 'target_radius' and 'target_nfuncs' !
!    !      Output - ref_dens(:)                                        !
!    ! lpl: TODO - automatic ion generation w/o config specification    !
!    !      (right now only good for select range of ionic species)     !
!    ! lpl: Latest mod 18-07-2012                                       !
!    !------------------------------------------------------------------!
!    ! InOut : ref_dens - Reference density object for one atom         !
!    ! Input : ionic_charge - The ionic charge of the reference density !
!    !            to be generated                                       !
!    !       : input_target_radius - Radius for the atomic solver       !
!    !       : target_nfuncs - Number of orbitals for the atomic solver !
!    !       : proc - The proc to exectue this (this will likely be     !
!    !            removed and everything executed on root)              !
!    !------------------------------------------------------------------!
!
!    use atom, only: ATOM_TYPE, BASIS_TYPE, atom_create, atom_destroy, &
!         atom_create_basis, atom_destroy_basis, atom_solve, &
!         atom_write_orbitals, atom_get_lmax, atom_split_orbital, &
!         atom_polarise_orbital
!    use comms, only : comms_bcast, pub_my_proc_id, pub_on_root, &
!         pub_root_proc_id
!    use constants, only: DP, PI, stdout, VERBOSE
!    use ion, only: ELEMENT
!    use rundat, only: cutoff_energy, pub_rootname, &
!         pub_initial_dens_realspace, pub_cond_calculate, &
!         ! pub_ddec_refdens_init_mode, &
!         pub_output_detail, pub_ddec_write, pub_ddec_write_rcomp, &
!         pub_ddec_renormalize_refdens, pub_ddec_rcomp_maxit, &
!         pub_ddec_rcomp_econv, pub_ddec_ref_shell_mode, &
!         pub_ddec_rad_rcut
!    use utils, only: utils_alloc_check, utils_dealloc_check, &
!         utils_abort
!    use services, only: services_linear_interpolation, &
!         services_radial_integral_rmax, services_regular_integral
!
!    implicit none
!
!    ! lpl: InOut
!    type(DDEC_REF_DENS), intent(inout) :: ref_dens
!
!    ! lpl: Input
!    integer, intent(in) :: ionic_charge
!    real(kind=DP), intent(in) :: input_target_radius
!    integer, intent(in) :: target_nfuncs
!    integer, intent(in) :: proc
!
!    ! Local Variables
!    integer :: nfuncs
!    type(ATOM_TYPE) :: tatom
!    type(BASIS_TYPE) :: tbasis
!    character(len=128) :: report(400)
!    character(len=80) :: config
!    character(len=128) :: config_long
!    character(len=80) :: filename
!    character(len=20) :: tmp
!    real(kind=DP) :: cwidth, cwidth_l(0:5), cscale
!    real(kind=DP) :: rmax, rmax_l(0:5)
!    real(kind=DP) :: target_radius, rcomp_target_radius
!    integer :: nreport, ireport, pos
!    integer :: lmax
!    logical :: fix_derivs
!
!    real(kind=DP) :: r, rad_r1, rad_r2, dens_r1, dens_r2
!    integer :: it, rpt, num
!
!    ! Stuff for compensated charge
!    real(kind=DP), allocatable :: comp_pot(:)
!    real(kind=DP) :: rcomp, rcomp_min, rcomp_max, rcomp_target, rcomp_incr
!    real(kind=DP) :: comppot_energy, comppot_eq, comppot_nq, comppot_qq
!
!    real(kind=DP) :: dE, E_grad, c0, c1, c2, iE(0:2), ir(0:2)
!    real(kind=DP), parameter :: E_scale = 1000.0_DP
!
!    real(kind=DP) :: occ_diff
!    integer :: n_occ_diff, orb_it
!
!    integer :: icomp, old_it, new_it, E_sign
!    logical :: break_loop, auto_iterate, iterate_once
!
!    integer :: ierr
!
!       if ( pub_debug_on_root ) write(stdout,'(a)') &
!            ' DEBUG: Entering ddec_solve_atom'
!
!    ierr = 0
!    ref_dens%ion_dens(ionic_charge)%dens = 0.0_DP
!
!    ! lpl: Generate density on specified proc
!    if(pub_my_proc_id == proc) then
!
!       ! lpl: Make sure subroutines runs on root easier management
!       if ( pub_my_proc_id /= pub_root_proc_id ) call utils_abort(' ERROR: &
!            &ddec_solve_atom not called on root')
!
!       ! lpl: Read r_comp parameters
!       rcomp_target = ref_dens%rcomp(ionic_charge,1)
!       rcomp_min  = ref_dens%rcomp(ionic_charge,2)
!       rcomp_max  = ref_dens%rcomp(ionic_charge,3)
!       rcomp_incr = ref_dens%rcomp(ionic_charge,4)
!       rcomp_target_radius = ref_dens%rcomp(ionic_charge,5)
!
!       ! lpl: Override 'target_radius' if species-specific value exists
!       target_radius = -1.0_DP
!       if ( rcomp_target_radius > 0 ) then
!          target_radius = rcomp_target_radius
!       else
!          target_radius = input_target_radius
!       end if
!
!       ! lpl: Set up config string for atom_solve depending on ionic charge
!       !      This will usually be the minimal basis unless user has
!       !      manually specified their own electronic configuration
!       if ( ref_dens%refconf(ionic_charge) /= '' ) then
!          config = trim(adjustl( ref_dens%refconf(ionic_charge) ))
!          ! lpl: MOD#09 - Set this to -1 so that 'atom_get_lmax' will
!          !      determine it based on 'config'
!          nfuncs = -1
!       else
!          call internal_ion_config(config,nfuncs,ref_dens%element, &
!               ionic_charge)
!       end if
!       config_long = adjustl(config)
!
!       ! lpl: All parameters here hard-fixed to defaults for now.
!       !      This can be modified to read user-specified flags later,
!       !      ideally by calling 'ngwfs_solve_atom' from 'ngwfs_mod'.
!       ! Read the confining potential width
!       cwidth_l(:) = 3.0_DP ! Default
!       rmax_l(:) = target_radius
!       ! Read the confining potential height (set to zero for now)
!       ! cscale = 100.0_DP
!       cscale = 0.0_DP
!       fix_derivs = .false.
!
!       ! lpl: Set up the basis and the atom
!       ! lpl: Assuming target_nfuncs doesn't affect atom_solve and is only used
!       !      to count number of final initialized NGWFs, but for some reason is
!       !      used in atom_get_lmax (in atom_orb_range), target_nfuncs is set to
!       !      number of non-vacant orbitals in gs_occ
!#ifdef DEBUG
!       ! lpl: MOD#09
!       write(stdout,'(a)') ' SOLVE: solver config = "'//trim(config)//'"'
!#endif
!       call atom_get_lmax(lmax,ref_dens%element,nfuncs,config)
!       call atom_create_basis(tbasis,rmax_l,lmax,6.0_DP*cutoff_energy, &
!            fix_derivs)
!       ! This is done solely to extract occupancy and angular momentum info
!       ! for each orbital based on the config string without having to write
!       ! a separate code without presumtions regarding orbital ordering. This
!       ! 'tatom' will be destroyed in the 1st iteration of the rcomp loop
!       call atom_create(tatom,tbasis,ref_dens%element,nfuncs, &
!            config,cwidth_l,cscale)
!
!       write(stdout,'(/a)') '===============================================&
!            &====================='
!       write(stdout,'(a)') &
!            '                Reference Density Solver Parameters'
!       write(stdout,'(a)') '-----------------------------------------------&
!            &---------------------'
!       write(stdout,'(a)') ' Symbol  Element   Z    Charge    Radius'
!       write(stdout,'(2x,A4,4x,A2,6x,I3,4x,I3,6x,f7.4)') &
!            ref_dens%element%species_id, ref_dens%element%symbol, &
!            ref_dens%element%atomic_number, ionic_charge, target_radius
!       write(stdout,'(a)') '-----------------------------------------------&
!            &---------------------'
!       write(stdout,'(a)') 'Config = "'//trim(adjustl(config))//'"'
!       write(stdout,'(a)') '-----------------------------------------------&
!            &---------------------'
!
!       ! lpl: Set array 'ref_dens%neutral_orb_info(:,:)' to store info
!       !      of orbitals and their <r> for cationic rcomp calculation later on
!       if ( ionic_charge == 0 ) then
!
!          ! lpl: This will have been initialized as -1.0_DP - so if neutral
!          !      calculation did not start, or has been performed, later
!          !      calculations will know about it
!          if( all(ref_dens%neutral_orb_info == 0.0_DP) ) &
!               call utils_abort(' ERROR: ref_dens%neutral_orb_info should &
!                    &not have been initialized before neutral density &
!                    &calculation is performed.')
!          ref_dens%neutral_orb_info = 0.0_DP
!          ref_dens%neutral_norbs = tatom%norbs
!
!       else if ( ionic_charge > 0 ) then
!          if ( any(ref_dens%neutral_orb_info < 0.0_DP) ) &
!               call utils_abort(' ERROR: ref_dens%neutral_orb_info should &
!                  &have been filled with neutral species info before cation &
!                  &density calculation is performed.')
!       end if
!
!       ! lpl: Alter condition to perform uncompensated density calculations
!       !      for neutral ionic species or for ionic species with unspecified
!       !      'rcomp_target'
!       auto_iterate = .false.
!       iterate_once = .false.
!       if ( ionic_charge >= 0 ) then
!          ! lpl: For neutral species, set target rcomp < 0 as no compensation
!          !      charge is required
!          if ( ionic_charge == 0 ) then
!             rcomp_target = -1.0_DP
!          else
!             if ( rcomp_min > 0 .and. rcomp_max < 0 ) then
!                ! lpl: Determine rcomp automatically if rcomp_max set to < 0
!                !      (to be consistent with anion calculations where setting
!                !      rcomp_max < 0 results in automatic calculation)
!                !      For cations rcomp = average of <psi|r|psi> for orbitals
!                !      in the neutral species which become vacant in the cation
!                !      For automatic config specification, the neutral species
!                !      will have been generated using the ground state config as
!                !      in the 'gs_occ' array (accesible here from the function
!                !      'atom_gs_occ') for a minimal basis (whose value is stored
!                !      in the variable 'nfuncs'
!                ! lpl: Assume ordering of orbitals are identical save for
!                !      possibly extra orbitals in the neutral species.
!                !      Otherwise spit error
!                n_occ_diff = 0
!                rcomp_target = 0.0_DP
!
!                write(stdout,'(a)') &
!                     ' Orb#  L  Neutral<r>  NeutrOcc  CationOcc DiffOcc'
!
!                do orb_it=1,ref_dens%neutral_norbs
!                   occ_diff = 0
!                   if ( orb_it > tatom%norbs ) then
!                      occ_diff = ref_dens%neutral_orb_info(orb_it,1)
!                   else
!                      occ_diff = ref_dens%neutral_orb_info(orb_it,1) &
!                           - tatom%occ(orb_it)
!                      ! lpl: Check that angular momenta are the same
!                      if ( nint(ref_dens%neutral_orb_info(orb_it,3)) &
!                           /= tatom%orb_ang_mom(orb_it) ) &
!                           call utils_abort(' ERROR: Angular momentum mismatch&
!                              & when attempting to determine rcomp for&
!                              & cation in ddec_solve_atom')
!                   end if
!
!                   ! lpl: Might happen for transition metals - not sure what
!                   !      to do here. For now, the 'newly occupied' orbitals in
!                   !      the cation are ignored and we only consider orbitals
!                   !      vacated from the neutral species in determining rcomp
!                   if ( occ_diff < 0 ) then
!                      write(stdout,'(a,/a)') ' WARNING: Cation has orbitals&
!                           & which are unoccupied','in the neutral species.&
!                           & The extra orbitals in the cation will not be&
!                           & considered for rcomp'
!                      occ_diff = 0
!                   end if
!
!                   n_occ_diff = n_occ_diff + occ_diff
!                   rcomp_target = rcomp_target + &
!                        occ_diff*ref_dens%neutral_orb_info(orb_it,2)
!
!                   if ( orb_it > tatom%norbs ) then
!                      c0 = 0
!                   else
!                      c0 = tatom%occ(orb_it)
!                   end if
!                   write(stdout,'(1x,I3,3x,I1,2x,f10.5,3(2x,f8.5))') orb_it, &
!                        nint(ref_dens%neutral_orb_info(orb_it,3)), &
!                        ref_dens%neutral_orb_info(orb_it,2), &
!                        ref_dens%neutral_orb_info(orb_it,1), c0, occ_diff
!
!                end do
!
!                if ( n_occ_diff <= 0 ) then
!                   call utils_abort(' ERROR: n_occ_diff <= 0 in &
!                        &ddec_solve_atom')
!                else
!                   rcomp_target = rcomp_target/real(n_occ_diff,kind=DP)
!                end if
!
!                write(stdout,'(a)') '----------------------------------------&
!                     &----------------------------'
!             end if ! END if ( rcomp_target < 0 .and. rcomp_max < 0 )
!          end if ! END if ( ionic_charge == 0 ) else
!          ! Set initial rcomp
!          rcomp = rcomp_target
!          ! Flag to break loop after one iteration
!          iterate_once = .true.
!
!       ! If 'rcomp_max' and 'rcomp_min' were unspecified, only iterate once
!       else if ( rcomp_max < 0 .and. rcomp_min < 0 ) then
!          iterate_once = .true.
!          rcomp = rcomp_target
!       ! If 'rcomp_min' is specified but not 'rcomp_max' for anion, then this
!       ! is an automated calculation
!       else if ( rcomp_max < 0 .and. rcomp_min > 0 ) then
!          auto_iterate = .true.
!          rcomp_max = target_radius
!
!          if ( rcomp_incr <= 0.0_DP ) &
!               call utils_abort(' ERROR: rcomp_incr <= 0')
!          ! Initialize rcomp while accounting for min/max radius
!          if ( rcomp_target > 0 ) then
!             if ( rcomp_target < rcomp_min + rcomp_incr ) then
!                rcomp = rcomp_min
!             else if ( rcomp_target > target_radius - rcomp_incr ) then
!                rcomp = rcomp_target - 2.0_DP*rcomp_incr
!             else
!                rcomp = rcomp_target - rcomp_incr
!             end if
!          else
!             rcomp = rcomp_min
!          end if
!       ! lpl: Else, assume this is manual calculation
!       else
!          if ( rcomp_max*rcomp_min < 0 ) then
!             call utils_abort(' ERROR: rcomp_max or rcomp_min < 0')
!          else if ( rcomp_max < rcomp_min ) then
!             call utils_abort(' ERROR: rcomp_max < rcomp_min')
!          else if ( rcomp_incr <= 0.0_DP ) then
!             call utils_abort(' ERROR: rcomp_incr <= 0')
!          end if
!          rcomp = rcomp_min
!       end if
!
!       ! lpl: Terminate if manual atomic configuration was specified
!       !      but conditions are set to automatically determine rcomp
!       !      since stuff will break due to assumption that the solver
!       !      configuration follows strictly the aufbau principle and
!       !      the ground state for each ion.
!       ! lpl: MOD#13 - Allow the above
!!       if ( auto_iterate .and. &
!!            ref_dens%refconf(ionic_charge) /= '' ) call utils_abort(' ERROR:&
!!            & Manual reference density configuration cannot be used with&
!!            & automatic rcomp calculation')
!
!       ! lpl: Calculations over a range of rcomps
!       allocate(comp_pot(tbasis%npts),stat=ierr)
!       call utils_alloc_check('internal_ddec_solve_atom', &
!            'comp_pot',ierr)
!
!       ! Initialize 1/r potential
!       comp_pot(1:tbasis%npts) = 1.0_DP/tbasis%rad(1:tbasis%npts)
!
!#ifdef DEBUG
!       write(stdout,'(a)') ' iter rcomp       Qcomp       &
!             &Etotal                Edensity              Qencl&
!             &            dE              old_it new_it'
!#else
!       write(stdout,'(a)') &
!            ' iter rcomp       Qcomp       TotalE                dE'
!#endif
!
!       ! Calculate compensated density for (sets of) rcomp
!       icomp = 0
!       break_loop = .false.
!       dE = pub_ddec_rcomp_econv
!       old_it = 1
!       new_it = 2
!
!       do while ( .not. break_loop )
!
!          icomp = icomp + 1
!
!          ! Make last calculation the calculation for target rcomp
!          if ( .not. auto_iterate ) then
!             if ( iterate_once .or. rcomp > rcomp_max ) then
!                if ( .not. iterate_once ) rcomp = rcomp_target
!                break_loop = .true.
!                if ( rcomp < 0.0_DP ) write(stdout,'(/a)') &
!                     '   ***No compensated charge calculation for&
!                     & this configuration***'
!             end if
!          else
!             if ( icomp > pub_ddec_rcomp_maxit ) then
!                rcomp = rcomp_target
!                break_loop = .true.
!                write(stdout,'(a,2(/a))') ' WARNING: Failed to find rcomp with &
!                     &minimum energy.','Defaulting to user-specified rcomp. &
!                     &If this was not set, no compensation','charge will &
!                     &be used.'
!             end if
!          end if
!
!          ! lpl: Destroy and recreate atom and basis. Resetting would be better
!          !      but this is easier to code.
!          call atom_destroy(tatom)
!          call atom_create(tatom,tbasis,ref_dens%element,nfuncs, &
!               config,cwidth_l,cscale)
!          ! lpl: Set rcomp, rcomppot in newly-created atom and solve
!          !      with compensation potential
!          tatom%rcomp_extra_scf = .true.
!          if ( rcomp >= 0.0_DP ) then
!             tatom%include_rcomppot = .true.
!             tatom%rcomp = rcomp
!             rpt = int( rcomp/(tbasis%dr) )
!             tatom%rcomppot(rpt:tbasis%npts) = comp_pot(rpt:tbasis%npts)
!             tatom%rcomppot(1:rpt-1) = comp_pot(rpt)
!          else
!             tatom%include_rcomppot = .false.
!          end if
!
!          ! lpl: Set comppot details depending on whether ion is anion or cation
!          if (ionic_charge < 0) then
!             tatom%update_rcomppot = .true.
!          else if(ionic_charge >= 0) then
!             tatom%update_rcomppot = .false.
!             ! lpl: qcomp defined as having negatively-charged electrons
!             tatom%qcomp = -1.0_DP*real(ionic_charge,kind=DP)
!          end if
!
!          report = ''
!          call atom_solve(tatom,tbasis,report,nreport,config_long)
!          if (pub_output_detail >= VERBOSE) then
!             do ireport=1,nreport
!                write(stdout,'(a)') trim(report(ireport))
!             end do
!          end if
!
!          ! lpl: Calculate electrosatic energy contribution from compensation
!          !      shell
!          comppot_eq = 0.0_DP
!          comppot_nq = 0.0_DP
!          comppot_qq = 0.0_DP
!          if( tatom%include_rcomppot ) then
!             ! E_eq (electron-qcomp interaction energy)
!             tatom%work(:) = tatom%den(:)*tbasis%rad(:)**2
!             comppot_eq = services_regular_integral(tbasis%npts,tbasis%dr, &
!                  -1.0_DP*tatom%qcomp*tatom%rcomppot(:)*tatom%work(:))
!
!             ! E_nq (nuclear-qcomp interaction energy)
!             comppot_nq = &
!                  -1.0_DP*tatom%qcomp*tatom%locpspot_orig(tatom%rcomp_ipt)
!             ! E_qq (qcomp-qcomp interaction energy)
!             comppot_qq = 0.5_DP*(tatom%qcomp**2)/(tatom%rcomp_ipt*tbasis%dr)
!          end if
!
!          ! E_total_corrected
!          comppot_energy = tatom%total_energy + comppot_nq + comppot_qq
!
!          if ( auto_iterate ) then
!             iE(new_it) = E_scale*comppot_energy
!             ir(new_it) = tatom%rcomp
!             dE = ( iE(new_it) - iE(old_it) )/E_scale
!
!             if ( icomp > 1 ) then
!                if ( rcomp_incr < 0.5_DP*tbasis%dr .and. &
!                     abs(dE) > pub_ddec_rcomp_econv ) then
!                   call utils_abort(' ERROR: rcomp_incr < 0.5*tbasis%dr but &
!                        &E(rcomp) has not converged.')
!                else if ( abs(dE) < pub_ddec_rcomp_econv ) then
!                   break_loop = .true.
!                end if
!             end if
!          end if
!
!          if (break_loop) then
!             write(stdout,'(a)') ''
!
!             ! lpl: Write density profile from atomic solver
!             if (pub_ddec_write .and. pub_ddec_write_rcomp) then
!                tmp = ''
!                write(tmp,'(sp,I4.3)') ionic_charge
!                filename = ''
!                if ( rcomp >= 0.0_DP ) then
!                   write(filename,'(f10.5)') rcomp
!                   filename = '_rcomp_'//trim(adjustl(filename))
!                end if
!                call ddec_write_dens( (tatom%den)/FOUR_PI,tbasis%rad, &
!                     tbasis%npts, &
!                     'solve_'//trim(adjustl(ref_dens%element%symbol))// &
!                     '_'//trim(adjustl(tmp))//trim(adjustl(filename)) )
!             end if
!          end if
!
!#ifdef DEBUG
!          c0 = 0.0_DP
!          if ( tatom%include_rcomppot ) c0 = &
!               services_regular_integral(tatom%rcomp_ipt,tbasis%dr, &
!               tatom%work(1:tatom%rcomp_ipt))
!          write(stdout,'(1x,I3,2(2x,f10.5),2(2x,f20.10),&
!               &2(2x,f15.10),2(3x,I4))') icomp, tatom%rcomp,tatom%qcomp, &
!               comppot_energy,tatom%total_energy - comppot_eq, c0, &
!               dE,old_it,new_it
!#else
!          write(stdout,'(1x,I3,2(2x,f10.5),2x,f20.10,2x,f15.10)') &
!               icomp, tatom%rcomp, tatom%qcomp, comppot_energy, dE
!#endif
!
!          ! Update rcomp accordingly
!          if (auto_iterate .and. .not. break_loop) then
!             if ( icomp > 2 ) then
!                ! lpl: Parabolic fit
!                ! f(r) = c0(r-r1)(r-r2) + c1(r-r0)(r-r2) + c2(r-r0)(r-r1)
!                c0 = iE(0)/( (ir(0)-ir(1))*(ir(0)-ir(2)) )
!                c1 = iE(1)/( (ir(1)-ir(0))*(ir(1)-ir(2)) )
!                c2 = iE(2)/( (ir(2)-ir(0))*(ir(2)-ir(1)) )
!                ! f"(r) = 2(c0+c1+c2)
!                ! f'(r) = c0(2r-r1-r2) + c1(2r-r0-r2) + c2(2r-r0-r1)
!                E_grad = c0*(ir(1)-ir(2)) + c1*(2.0_DP*ir(1)-ir(0)-ir(2)) &
!                     + c2*(ir(1)-ir(0))
!#ifdef DEBUG
!                write(stdout,'(a)') '-----------------------------------------'
!                write(stdout,'(a)') '  old_it new_it E_grad          rcomp_incr'
!                write(stdout,'(1x,I4,6x,I4,3(2x,f14.7))') old_it, new_it, &
!                     E_grad, rcomp_incr
!                write(stdout,'(a)') '-----------------------------------------'
!                write(stdout,'(a)') ' iE(0)          ir(0)           '
!                write(stdout,'(2(2x,f14.7))') iE(0), ir(0)
!                write(stdout,'(2(2x,f14.7))') iE(1), ir(1)
!                write(stdout,'(2(2x,f14.7))') iE(2), ir(2)
!                write(stdout,'(a)') '-----------------------------------------'
!#endif
!
!                ! Only keep the sign
!                E_sign = sign(1.0_DP,E_grad)
!
!                ! If minimum (both sides > center), reduce step size
!                if ( iE(2) > iE(1) .and. iE(0) > iE(1) ) then
!                   rcomp_incr = 0.5_DP*rcomp_incr
!                   rcomp = ir(1) - real(E_sign,kind=DP)*rcomp_incr
!
!                   old_it = (1 + E_sign)/2
!                   iE(0) = iE(1-old_it)
!                   ir(0) = ir(1-old_it)
!                   iE(2) = iE(2-old_it)
!                   ir(2) = ir(2-old_it)
!                   new_it = 1
!                   old_it = 1 + E_sign
!                ! Else increment by dr in correct direction
!                else
!                   rcomp = ir(1) - 2.0_DP*real(E_sign,kind=DP)*rcomp_incr
!                   if (E_sign > 0) then
!                      iE(2) = iE(1)
!                      ir(2) = ir(1)
!                      iE(1) = iE(0)
!                      ir(1) = ir(0)
!                   else
!                      iE(0) = iE(1)
!                      ir(0) = ir(1)
!                      iE(1) = iE(2)
!                      ir(1) = ir(2)
!                   end if
!                   new_it = 1 - E_sign
!                   old_it = new_it + E_sign
!                end if
!             else
!                rcomp = rcomp + rcomp_incr
!                iE(0) = iE(1)
!                ir(0) = ir(1)
!                iE(1) = iE(2)
!                ir(1) = ir(2)
!                new_it = 2
!                old_it = 1
!             end if
!
!             ! Make sure rcomp does not exceed these bounds
!             if( rcomp < rcomp_min ) then
!                rcomp = rcomp_min
!             else if ( rcomp > target_radius ) then
!                rcomp = rcomp_max
!             end if
!          else
!             rcomp = rcomp + rcomp_incr
!          end if
!
!       end do ! END do while ( rcomp <= rcomp_max )
!
!       if ( tatom%include_rcomppot ) then
!          write(stdout,'(a)') '----------------------------------------&
!               &----------------------------'
!          write(stdout,'(a,f10.5,a)') ' rcomp = ', tatom%rcomp, ' bohr'
!       end if
!       if ( ionic_charge /= 0 ) write(stdout,'(a)') '===================&
!            &================================================='
!
!       ! lpl: Write density profile from atomic solver
!       if (pub_ddec_write) then
!          tmp = ''
!          write(tmp,'(sp,I4.3)') ionic_charge
!          call ddec_write_dens( (tatom%den)/FOUR_PI,tbasis%rad, &
!               tbasis%npts,'solve_'//trim(adjustl(ref_dens%element%symbol))// &
!               '_'//trim(adjustl(tmp)) )
!       end if
!
!       ! lpl: Copy the atom density into the rad_dens
!       ! lpl: Mod 09122013 - Following CHARGEMOL, such densities are 'binned'
!       !      per radial shell and not interpolated, so make option available
!       ref_dens%ion_dens(ionic_charge)%dens = 0.0_DP
!       select case ( pub_ddec_ref_shell_mode )
!          case (0) ! Should be default - interpolate value to each shell
!             call ddec_1d_interp(ref_dens%ion_dens(ionic_charge)%dens, &
!                  ref_dens%rad,tatom%den/FOUR_PI,tbasis%dr,0.0_DP,.true.)
!
!             ! lpl: Optional - renormalize densities for this radial grid
!             if ( pub_ddec_renormalize_refdens ) then
!                rad_r1 = 0.0_DP
!                rad_r2 = ref_dens%element%ion_charge - ionic_charge
!                do it=1,ref_dens%npts
!                   ! lpl: Use ref_dens%avg_rad to approximate summing in
!                   !      local_box
!                   rad_r1 = rad_r1 + &
!                        FOUR_PI*(ref_dens%avg_rad(it)**2)*ref_dens%dr* &
!                        ref_dens%ion_dens(ionic_charge)%dens(it)
!                end do
!                if ( rad_r1 > 0.0_DP ) then
!                   ref_dens%ion_dens(ionic_charge)%dens = &
!                   ref_dens%ion_dens(ionic_charge)%dens*(rad_r2/rad_r1)
!                else if ( rad_r1 < 0.0_DP ) then
!                   call utils_abort('ERROR: Somehow, the integrated&
!                        & density is < 0')
!                end if
!                if (pub_on_root) write(stdout,'(a,f14.7,a,f14.7)') &
!                     ' Norm = ', rad_r1, '     Renorm = ', rad_r2
!             end if
!          case (1) ! Average over all radial values for each shell
!             call ddec_shell_avg(ref_dens%ion_dens(ionic_charge)%dens, &
!                  ref_dens%dr,tatom%den/FOUR_PI,tbasis%dr)
!          case (2) ! As above but via Cartesian grid
!             call ddec_rad_to_box_to_rad(ref_dens%ion_dens(ionic_charge)&
!               %dens, ref_dens%dr,tatom%den/FOUR_PI,tbasis%dr, &
!               pub_ddec_rad_rcut, ref_dens%element%ion_charge - &
!               ionic_charge, pub_ddec_renormalize_refdens)
!          case (3) ! Copies values directly from input
!             call ddec_copy_rad(ref_dens%ion_dens(ionic_charge)%dens, &
!                  ref_dens%dr,tatom%den/FOUR_PI,tbasis%dr)
!          case default
!             call utils_abort(' ERROR in ddec_solve_atom: Invalid option&
!                  & for pub_ddec_ref_shell_mode')
!       end select
!
!       if ( allocated(comp_pot) ) then
!          deallocate(comp_pot,stat=ierr)
!          call utils_dealloc_check('internal_ddec_solve_atom', &
!               'comp_pot',ierr)
!       end if
!
!       ! lpl: Calculate expectation value of radius
!       if ( ionic_charge == 0 ) then
!          write(stdout,'(a)')    '-------------------------------------------&
!               &-------------------------'
!          write(stdout,'(a)') '                        Neutral Species Info'
!          write(stdout,'(a)')    '-------------------------------------------&
!               &-------------------------'
!          write(stdout,'(a)') ' Orbital  L  Occupancy   <r>         Norm'
!          do it=1,tatom%norbs
!
!             ! lpl:  <psi|r|psi> for real |psi>
!             tatom%work(:) = tatom%psi_r(:,it)**2 * tbasis%rad(:)**3
!             r = services_regular_integral(tbasis%npts,tbasis%dr,tatom%work(:))
!             ! lpl: <psi|psi>
!             tatom%work(:) = tatom%psi_r(:,it)**2 * tbasis%rad(:)**2
!             rad_r1 = services_regular_integral(tbasis%npts,tbasis%dr, &
!                  tatom%work(:))
!             ! lpl: Report <psi|r|psi>
!             write(stdout,'(3x,I3,4x,I1,3(2x,f10.5))') &
!                  it, tatom%orb_ang_mom(it), tatom%occ(it), r, rad_r1
!             ! lpl: Store in array
!             ref_dens%neutral_orb_info(it,1) = tatom%occ(it)
!             ref_dens%neutral_orb_info(it,2) = r
!             ref_dens%neutral_orb_info(it,3) = tatom%orb_ang_mom(it)
!
!#ifdef DEBUG
!             ! lpl: Write orbitals to files
!             if (pub_ddec_write .and. pub_ddec_write_rcomp) then
!                tmp = ''
!                write(tmp,'(sp,I4.3)') ionic_charge
!                filename = ''
!                write(filename,'(I3.3)') it
!                call ddec_write_dens( (tatom%psi_r(:,it)),tbasis%rad, &
!                     tbasis%npts, &
!                     'solve_'//trim(adjustl(ref_dens%element%symbol))// &
!                     '_'//trim(adjustl(tmp))//'_orb_'//trim(adjustl(filename)) )
!             end if
!#endif
!          end do
!
!          write(stdout,'(a)') '===========================================&
!               &========================='
!
!       end if ! END if ( ionic_charge == 0 )
!
!       ! lpl: Clean up the atom and basis
!       call atom_destroy(tatom)
!       call atom_destroy_basis(tbasis)
!
!       write(stdout,'(a)') ''
!
!    end if ! END if(pub_my_proc_id == proc)
!
!    call comms_bcast(proc,ref_dens%ion_dens(ionic_charge)%dens)
!
!#ifdef DEBUG
!       if ( pub_on_root ) write(stdout,'(a)') &
!            ' DEBUG: Leaving ddec_solve_atom'
!#endif
!
!   contains
!
!  !===========================================================================!
!
!      subroutine internal_ion_config(config_string,target_nfuncs, &
!         current_element,ionic_charge)
!
!        !------------------------------------------------------------------!
!        ! lpl: Generates configuration string and target_nfuncs based on   !
!        !      ionic charge for use in atomic solver. Might break for      !
!        !      heavy elements. Will likely generate superflous number of   !
!        !      vacant orbitals.                                            !
!        !------------------------------------------------------------------!
!        ! Output: config_string - String containing configuration for      !
!        !            atomic solver                                         !
!        !       : target_nfuncs - Number of orbitals for atomic solver     !
!        ! Input : current_element - The ELEMENT object of the current atom !
!        !       : ionic_charge - The ionic charge of the configuration     !
!        !------------------------------------------------------------------!
!
!        use atom, only: atom_gs_occ
!
!        implicit none
!
!        ! lpl: Output
!        character(len=80), intent(out) :: config_string
!        integer, intent(out) :: target_nfuncs
!
!        ! lpl: Input
!        type(ELEMENT), intent(in) :: current_element
!        integer, intent(in) :: ionic_charge
!
!        ! lpl: Internal variables
!        character(len=2) :: occ_str
!        character(len=6) :: term
!        character :: l_char
!        integer :: ion_num
!        integer :: n_left, nl_min, nl_max, nl_index_size, it
!        integer :: orb_n, orb_l
!
!#ifdef DEBUG
!       if ( pub_on_root ) write(stdout,'(a)') &
!            ' DEBUG: Entering internal_ion_config'
!#endif
!
!        nl_index_size = atom_gs_occ(2,0,'Size')
!
!        ! lpl: Assume that ionic configuration are simply atomic configurations
!        !      of neutral element with the same number of electrons. Does not
!        !      always work for transition metals or possibly heavy elements
!        ion_num = current_element%atomic_number - ionic_charge
!
!        if(ion_num < 0 .or. ion_num > atom_gs_occ(1,0,'Size')) &
!             call utils_abort(' ERROR: ion_num out of range in &
!                &internal_ion_config in ddec_mod')
!
!        ! lpl: Check that ion_num does not exceed ion_charge due to pseudized
!        !      core charges
!        if(ionic_charge > current_element%ion_charge) &
!             call utils_abort(' ERROR: ion_num > ion_charge in &
!                &internal_ion_config in ddec_mod')
!
!        ! lpl: Get num of pseudized nl shells based on floor of pseudized charge
!        n_left = int(real(current_element%atomic_number,kind=DP) &
!             - current_element%ion_charge)
!        nl_min = 1
!        do while(n_left > 0)
!           n_left = n_left - atom_gs_occ(ion_num,nl_min,'Occupancy')
!           nl_min = nl_min + 1
!        end do
!        if(nl_min > nl_index_size) call utils_abort(' ERROR: nl_min out of&
!             & range in internal_ion_config in ddec_mod')
!
!        ! lpl: Get max occpied nl_index for this ion_num
!        do nl_max=nl_index_size,nl_min,-1
!           ! lowest nl_max is minbas
!           if(atom_gs_occ(ion_num,nl_max,'Occupancy') /= 0 .or. &
!                atom_gs_occ(current_element%atomic_number,nl_max,&
!                'Occupancy') /= 0) exit
!        end do
!        if(nl_max < nl_min) nl_max = nl_min
!
!        config_string = ''
!        target_nfuncs = 0
!        do it=nl_min,nl_max
!
!           ! lpl: Get nl index in terms of config string format
!           write(occ_str,'(I2)') atom_gs_occ(0,it)
!           occ_str = adjustl(occ_str)
!           read(occ_str(1:1),'(I1)') orb_n
!           read(occ_str(2:2),'(I1)') orb_l
!
!           occ_str = ''
!           if ( ion_num /= 0 ) then
!              write(occ_str,'(I2)') atom_gs_occ(ion_num,it,'Occupancy')
!           else
!              write(occ_str,'(I2)') 0
!           end if
!
!           select case (orb_l)
!              case (0); l_char='s'
!              case (1); l_char='p'
!              case (2); l_char='d'
!              case (3); l_char='f'
!              case (4); l_char='g'
!              case (5); l_char='h'
!              case (6); l_char='i'
!              case (7); l_char='j'
!              case (8); l_char='k'
!              case (9); l_char='l'
!           end select
!
!           ! lpl: Accumulate number of required functions
!           target_nfuncs = target_nfuncs + (2*orb_l + 1)
!
!           term = ''
!           write(term,'(a,I2,A1,A2)') ' ',orb_n,l_char,adjustl(occ_str)
!           config_string = trim(adjustl(config_string))//term
!
!        end do
!
!        config_string = adjustl(config_string)
!
!#ifdef DEBUG
!       if ( pub_on_root ) write(stdout,'(a)') &
!            ' DEBUG: Leaving internal_ion_config'
!#endif
!
!      end subroutine internal_ion_config
!
!  end subroutine ddec_solve_atom
!
!  !===========================================================================!

  ! lpl: COMMENT#05 - This routine reads in a string obtained from the
  !      reference density block, finds the relevant ASCII reference density
  !      file, and reads in the appropriate reference density. Feel free to
  !      change, remove, or replace with something better.
  subroutine ddec_read_atom(ref_dens,ionic_charge,filename,reset)

    !--------------------------------------------------------------------------!
    ! This routine reads core and valence reference densities from a formatted !
    ! ASCII file. The radial cutoff is determined by the file itself.          !
    !--------------------------------------------------------------------------!
    ! InOut : ref_dens - Reference density object for ONE species              !
    ! Input : ionic_charge - ionic charge to be read                           !
    !         filename - File name containing reference densities for all      !
    !            ionic species                                                 !
    !       : reset - optional, tells the sobroutine to deallocate the static  !
    !            density buffer array 'dens_buffer'                            !
    !--------------------------------------------------------------------------!

    use comms, only: comms_bcast, pub_on_root
    use constants, only: ANGSTROM, DP, stdout
    use rundat, only: pub_ddec_zero_thresh, &
         pub_ddec_rad_rcut, pub_ddec_renormalize_refdens
    use utils, only: utils_abort, utils_alloc_check, utils_close_unit_check, &
         utils_dealloc_check, utils_open_unit_check, utils_read_check, &
         utils_unit

    implicit none

    ! InOut
    type(DDEC_REF_DENS), intent(inout) :: ref_dens

    ! Input
    integer, intent(in) :: ionic_charge
    character(len=*), intent(in) :: filename
    logical, optional, intent(in) :: reset

    ! Local variables
    type(KEY_PARAM_PAIR), target :: key_param(12)
    type(KEY_PARAM_PAIR_PTR), allocatable, target :: key_param_ptr(:)
    type(KEY_PARAM_PAIR), pointer :: temp_key_param_ptr

    integer :: it, num_keys_remaining
    integer :: file_unit
    logical :: key_found, file_exists
    character(len=1024) :: line_buffer

    integer :: atnum_read, ionic_charge_read
    integer :: eff_nuc_charge, electron_count
    integer :: ddec_npts, net_magnetization
    real(kind=DP) :: ddec_rcut, ddec_dr

    ! lpl: Unimportant variables, comment out if desired
    character(len=512) :: basis_set_str, xc_string
    character(len=512) :: prepared_by_str, prepared_date_str

    ! lpl: Static items
    real(kind=DP), allocatable, save :: dens_buffer(:)

    integer :: ierr

    if ( pub_debug_on_root )  write(stdout,'(a)') &
         ' DEBUG: Entering ddec_read_atom'

    ! lpl: Deallocate static array and reset static variable
    if ( present(reset) ) then
       if ( reset ) then

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Resetting ddec_read_atom'

          if ( allocated(dens_buffer) ) then
             deallocate(dens_buffer,stat=ierr)
             call utils_dealloc_check('ddec_read_atom','dens_buffer',ierr)
          end if
          return
       end if
    end if

    ! lpl: Make sure called from root
    if ( .not. pub_on_root ) call utils_abort(' ERROR: ddec_read_atom&
         & not called from root')

    ierr = 0
    file_exists = .false.

    ! lpl: Parameters to be read
    atnum_read = -1
    ionic_charge_read = -999
    eff_nuc_charge = -1
    electron_count = -1
    ddec_npts = -1
    net_magnetization = -1
    ddec_rcut = -1.0_DP
    ! lpl: Unimportant variables that could be dropped
    basis_set_str = ''
    xc_string = ''
    prepared_by_str = ''
    prepared_date_str = ''

    ! lpl: Relevant data keys
    do it=1,size(key_param)
       key_param(it)%param = ''
       key_param(it)%delimiter = '='
       key_param(it)%essential = .false.
       key_param(it)%initialized = .false.
    end do
    key_param(12)%delimiter = ''

    key_param( 1)%key = 'atomic number'
    key_param( 2)%key = 'effective nuclear charge'
    key_param( 3)%key = 'electron count'
    key_param( 4)%key = 'net charge'
    key_param( 5)%key = 'basis set'
    key_param( 6)%key = 'exchange correlation functional'
    key_param( 7)%key = 'cutoff radius (picometer)'
    key_param( 8)%key = 'number of shells'
    key_param( 9)%key = 'prepared by'
    key_param(10)%key = 'prepared date'
    key_param(11)%key = 'magnetization'
    key_param(12)%key = 'density in atomic units:'

    key_param( 1)%essential = .true.
    key_param( 2)%essential = .true.
    key_param( 3)%essential = .true.
    key_param( 4)%essential = .true.
    key_param( 7)%essential = .true.
    key_param( 8)%essential = .true.
    key_param(11)%essential = .true.
    key_param(12)%essential = .true.

    ! lpl: Pointer for sorting - since linked list is too much bulk
    allocate(key_param_ptr(size(key_param)),stat=ierr)
    call utils_alloc_check('ddec_read_atom','key_param_ptr',ierr)
    if ( size(key_param) /= size(key_param_ptr) ) &
         call utils_abort(' ERROR in ddec_read_atom: size(key_param) &
              &/= size(key_param_ptr)')
    do it=1,size(key_param_ptr)
       key_param_ptr(it)%ptr => key_param(it)
    end do

    ! lpl: Check if file exist and open (duh)
    inquire(file=trim(filename),exist=file_exists)
    if ( .not. file_exists ) call utils_abort(" ERROR in ddec_read_atom: &
         &file '"//trim(adjustl(filename))//"' not found")

    file_unit = utils_unit()
    open(file_unit,file=trim(filename),form='formatted',status='old', &
         position='rewind',iostat=ierr)
    call utils_open_unit_check('ddec_read_atom','file_unit',ierr)

    ierr = 0
    key_found = .false.
    num_keys_remaining = size(key_param_ptr)
    ! lpl: Find 1st (header) key
    do while( .not. key_found )

       line_buffer = ''
       read(file_unit,'(a)',iostat=ierr) line_buffer
       ! lpl: Usually fails if ionic species not found (EOF error). Do not
       !      abort if we fail to read in a reference density that should be
       !      empty
       ! lpl: MOD#07 - do not crash but return 'not initialized' flag instead
       if ( ierr /= 0 ) then
          if ( ref_dens%element%atomic_number /= &
               ref_dens%eff_nuc_charge ) then
             ref_dens%init_status(ionic_charge) = .false.
             close(file_unit,iostat=ierr)
             call utils_close_unit_check('ddec_read_atom','file_unit',ierr)
             return
!             call utils_read_check('ddec_read_atom','line_buffer #1',ierr)
          else
             write(stdout,'(a)') ' WARNING: No reference density found, but &
                  &atomic_number == eff_nuc_charge, so this reference density &
                  &will have no electrons.'
             exit
          end if
       end if

       ! lpl: Skip blank lines and those which are blank upon comment removal
       if ( internal_format_line(line_buffer,'!') <= 0 ) cycle
       call internal_read_parameter(key_param(1),line_buffer,key_found)

       ! lpl: If header key found test for correctness
       if ( key_found ) then
          read(key_param%param,*) atnum_read
          ! lpl: If correct header found search for 2nd key to determine
          !      correctness of subheader
          if ( ref_dens%element%atomic_number == atnum_read ) then

             ! lpl: Swap keys that have been read
             temp_key_param_ptr => key_param_ptr(1)%ptr
             key_param_ptr(1)%ptr => key_param_ptr(num_keys_remaining)%ptr
             key_param_ptr(num_keys_remaining)%ptr => temp_key_param_ptr

             ! lpl: Decrement number of keys left
             num_keys_remaining = num_keys_remaining - 1

             ! lpl: Read the rest of the keys to find 2nd subheader key
             do while ( key_found )

                ! lpl: Read line from file. No blank lines allowed, only
                !      comments
                read(file_unit,'(a)',iostat=ierr) line_buffer
                ! lpl: Error here is usually due to failure to parse keys
                call utils_read_check('ddec_read_atom','line_buffer #2',ierr)

                ! lpl: If there is a comment which upon removal results in a
                !      blank line
                if ( internal_format_line(line_buffer,'!') < 0 ) cycle

                ! lpl: Loop over until match is found
                do it=1,num_keys_remaining
                   key_found = .false.
                   call internal_read_parameter(key_param_ptr(it)%ptr, &
                        line_buffer,key_found)

                   ! lpl: If any key is found then inspect key for 2nd
                   !      subheader, or continue accumulating. Clause will exit
                   !      when the 2nd subheader is wrong, whereby 'key_found'
                   !      is changed to '.false.', or when the last key is
                   !      found, whereby 'key_found' is left unchanged.
                   if ( key_found ) then

                      ! lpl: Swap keys that have been read so we won't have to
                      !      iterate over the entire list every time. The
                      !      current ptr index corresponding to the parameter
                      !      currently observed will be of value
                      !      'num_keys_remaining + 1'
                      temp_key_param_ptr => key_param_ptr(it)%ptr
                      key_param_ptr(it)%ptr => &
                           key_param_ptr(num_keys_remaining)%ptr
                      key_param_ptr(num_keys_remaining)%ptr => &
                           temp_key_param_ptr

                      ! lpl: Decrement number of keys left
                      num_keys_remaining = num_keys_remaining - 1

                      ! lpl: If this subheader is the one we are looking for, or
                      !      the last one (i.e. there is a key that should be
                      !      placed last in the parameter file, which serves to
                      !      terminate the parameter block)
                      if ( key_param_match(key_param_ptr(num_keys_remaining+1)&
                           %ptr,key_param(4)) ) then
                         read(key_param_ptr(num_keys_remaining+1)%ptr%param,*) &
                              ionic_charge_read

                         ! lpl: If the ionic charge is not correct then exit the
                         !      loop that reads in all other parameters. Also
                         !      sets 'key_found = .false.' to cause outer loop
                         !      to reset all parameters read and continue
                         !      parsing file
                         if ( ionic_charge /= ionic_charge_read ) then
                           key_found = .false.
                           ! lpl: Use this flag to indicate wrong 2nd header was
                           !      read, and that the entire block should be
                           !      skipped - since this is the terminating key it
                           !      should never have been read and initialized
                           !      at this point anyway
                           key_param(12)%initialized = .false.
                           exit
                         end if
                      else if ( key_param_match(key_param_ptr&
                           (num_keys_remaining+1)%ptr,key_param(12)) ) then
                         ! lpl: If the 2nd header has not been found when
                         !      arriving at the terminating key
                         if ( .not. key_param(4)%initialized ) &
                              call utils_abort(' ERROR incomplete block')
                         key_found = .false.
                         exit
                      end if

                      ! lpl: Exit loop if key is found
                      exit
                   else
                      ! This is needed for the if clause below the loop
                      if ( it >= num_keys_remaining ) &
                           call utils_abort(" ERROR in ddec_read_atom: unknown&
                                & key '"//trim(adjustl(line_buffer))//"'")
                   end if ! END if ( key_found )

                end do ! END do it=1,num_keys_remaining

                ! lpl: This condition indicates that the incorrect 2nd header
                !      was read
                if ( .not. (key_found .or. key_param(12)%initialized) ) then
                   ! lpl: Reset keys and parameters that have been read
                   do it=1,size(key_param_ptr)
                      ! lpl: Reset pointer array
                      key_param_ptr(it)%ptr => key_param(it)
                      ! lpl: Reset key_param object (essential)
                      key_param_ptr(it)%ptr%param = ''
                      key_param_ptr(it)%ptr%initialized = .false.
                   end do
                   num_keys_remaining = size(key_param_ptr)
                end if

             end do ! END do while ( key_found ) (inner loop to read keys)
          else
             ! lpl: Stop reading block if 1st header is incorrect
             key_found = .false.
          end if ! END if ( ref_dens%element%atomic_number == atnum_read )

       end if ! END if ( key_found )

       ! lpl: If correct 1st and 2nd header were found, then validate
       !      correctness of data read from keys
       if ( key_param(4)%initialized .and. key_param(12)%initialized ) then

          ! lpl: Check that important keys have been read correctly
          do it=1,(size(key_param_ptr) - num_keys_remaining)
             ! lpl: All keys should have been initialized by this point
             if ( .not. key_param_ptr(it)%ptr%initialized .and. &
                  key_param_ptr(it)%ptr%essential) &
                  call utils_abort(' ERROR not all essential keys were &
                       &properly read')
          end do

          ! lpl: Read all keys into variables
          if ( key_param( 1)%initialized ) &
               read(key_param( 1)%param,*,iostat=ierr) atnum_read
          if ( key_param( 2)%initialized ) &
               read(key_param( 2)%param,*,iostat=ierr) eff_nuc_charge
          if ( key_param( 3)%initialized ) &
               read(key_param( 3)%param,*,iostat=ierr) electron_count
          if ( key_param( 4)%initialized ) &
               read(key_param( 4)%param,*,iostat=ierr) ionic_charge_read
          if ( key_param( 5)%initialized ) &
               read(key_param( 5)%param,*,iostat=ierr) basis_set_str
          if ( key_param( 6)%initialized ) &
               read(key_param( 6)%param,*,iostat=ierr) xc_string
          if ( key_param( 7)%initialized ) &
               read(key_param( 7)%param,*,iostat=ierr) ddec_rcut
          if ( key_param( 8)%initialized ) &
               read(key_param( 8)%param,*,iostat=ierr) ddec_npts
          if ( key_param( 9)%initialized ) &
               read(key_param( 9)%param,'(a)',iostat=ierr) prepared_by_str
          if ( key_param(10)%initialized ) &
               read(key_param(10)%param,'(a)',iostat=ierr) prepared_date_str
          if ( key_param(11)%initialized ) &
               read(key_param(11)%param,*,iostat=ierr) net_magnetization

          ! lpl: Error here indicates failure to convert read data
          call utils_read_check('ddec_read_atom','key_param%param',ierr)
          if ( min(atnum_read,eff_nuc_charge,electron_count, &
               nint(ddec_rcut),ddec_npts,net_magnetization) < 0 ) &
               call utils_abort(' ERROR in ddec_read_atom: Some parameters &
                    &were not properly initialized')

          ! lpl: Check important keys are correct
          if ( eff_nuc_charge > atnum_read ) &
               call utils_abort(' ERROR in ddec_read_atom: eff_nuc_charge &
                    &> atnum_read')
          if ( eff_nuc_charge - electron_count /= ionic_charge_read ) &
               call utils_abort(' ERROR in ddec_read_atom: charge is incorrect')

          if ( ddec_npts < ref_dens%npts ) write(stdout,'(a)') &
               ' WARNING: nshells < ref_dens%npts'
          ! lpl: ddec_rcut must be in picometers
          ddec_rcut = 0.01_DP*ANGSTROM*ddec_rcut ! Convert from pm to ang
          if ( ddec_rcut > ref_dens%rcut + pub_ddec_zero_thresh ) &
               call utils_abort(' ERROR in ddec_read_atom: &
                    &rcut < ref_dens%rcut')
          ddec_dr = ddec_rcut/ddec_npts
          if ( ddec_dr > ref_dens%dr + pub_ddec_zero_thresh ) &
               write(stdout,'(a)') ' WARNING: dr < ref_dens%dr'

          if (pub_debug) then
             write(stdout,'(4(1x,a,1x,e20.10))') 'rc1=', ddec_rcut, &
                  'rc2=', ref_dens%rcut, 'dr1=', ddec_dr, 'dr2=', ref_dens%dr
          end if

          ! lpl: Allocate array to read in density data
          if ( allocated(dens_buffer) ) then
             if ( size(dens_buffer) < ddec_npts ) then
                deallocate(dens_buffer,stat=ierr)
                call utils_dealloc_check('ddec_read_atom','dens_buffer',ierr)
                allocate(dens_buffer(5*ddec_npts),stat=ierr)
                call utils_alloc_check('ddec_read_atom','dens_buffer',ierr)
             end if
          else
             allocate(dens_buffer(5*ddec_npts),stat=ierr)
             call utils_alloc_check('ddec_read_atom','dens_buffer',ierr)
          end if

          ! lpl: Read density block - no blank lines except comments allowed
          ierr = 0
          dens_buffer = 0.0_DP
          do it=1,ddec_npts
             read(file_unit,*,iostat=ierr) line_buffer
             ! lpl: Error here indicates failure to read density data
             call utils_read_check('ddec_read_atom','line_buffer #3',ierr)
             if ( internal_format_line(line_buffer,'!') < 0 ) cycle
             read(line_buffer,*,iostat=ierr) dens_buffer(it)
             ! lpl: Error here indicates failure to convert read data
             call utils_read_check('ddec_read_atom','line_buffer #4',ierr)
          end do
          ! lpl: This terminates the outermost loop
          key_found = .true.
       else
          key_found = .false.
       end if
    end do ! do while ( .not. key_found )

    ! lpl: Close input file
    close(file_unit,iostat=ierr)
    call utils_close_unit_check('ddec_read_atom','file_unit',ierr)

    ! lpl: Either transfer density read from file into memory or set it
    !      to zero if atomic_number == ionic_charge
    if ( .not. key_found ) then
       if( ref_dens%element%atomic_number == ionic_charge ) then
          ! lpl: Set density to zero if atomic_number == ionic_charge
          !      and no such density exist in the input file
          ref_dens%ion_dens(ionic_charge)%dens = 0.0_DP
       else
          ! lpl: Report that we did not find the required density
          ref_dens%init_status(ionic_charge) = .false.
          close(file_unit,iostat=ierr)
          call utils_close_unit_check('ddec_read_atom','file_unit',ierr)
          return
       end if
    else
       ref_dens%ion_dens(ionic_charge)%dens = 0.0_DP
       select case ( 3 )
!       select case ( pub_ddec_ref_shell_mode )
          case (0) ! Should be default - interpolate value to each shell
             ! lpl: Transfer data into ref_dens by quartic interpolation
             call ddec_1d_interp(ref_dens%ion_dens(ionic_charge)%dens, &
                  ref_dens%rad,dens_buffer(1:ddec_npts),ddec_dr, &
                  0.5_DP*ddec_dr,.false.)
          case (1) ! Average over all radial values for each shell
             call ddec_shell_avg(ref_dens%ion_dens(ionic_charge)%dens, &
                  ref_dens%dr,dens_buffer(1:ddec_npts),ddec_dr)
          case (2) ! Deposit to box and spherically-averge for each shell
             if ( pub_ddec_renormalize_refdens ) then
                call ddec_rad_to_box_to_rad(ref_dens%ion_dens(ionic_charge)&
                     %dens, ref_dens%dr,dens_buffer(1:ddec_npts),ddec_dr, &
                     pub_ddec_rad_rcut, &
                     ref_dens%element%ion_charge - ionic_charge,.true.)
             else
                call ddec_rad_to_box_to_rad(ref_dens%ion_dens(ionic_charge)&
                     %dens, ref_dens%dr,dens_buffer(1:ddec_npts),ddec_dr, &
                     pub_ddec_rad_rcut, &
                     ref_dens%element%ion_charge - ionic_charge,.false.)
             end if
          case (3) ! Copies values directly from input
             call ddec_copy_rad(ref_dens%ion_dens(ionic_charge)%dens, &
                  ref_dens%dr,dens_buffer(1:ddec_npts),ddec_dr)
          case default
             call utils_abort(' ERROR in ddec_solve_atom: Invalid option&
                  & for pub_ddec_ref_shell_mode')
       end select

        ref_dens%init_status(ionic_charge) = .true.
    end if

       if ( pub_debug_on_root )  write(stdout,'(a)') &
            ' DEBUG: Leaving ddec_read_atom'

  contains

  !===========================================================================!

    subroutine internal_read_parameter(key_param,line_in,found)

      !------------------------------------------------------------------------!
      ! Reads key-parameter pair from an input 'line_in' string. Aborts        !
      ! program if attempting to read a variables flagged as initialized.      !
      !------------------------------------------------------------------------!
      ! InOut : key_param - key-parameter pair. Overwrites the %param string   !
      !            and sets the %initialized flag to '.true.' if key is found  !
      ! Input : line_in - input line string containing a '<key><delim><param>' !
      !            string                                                      !
      ! Output: found - optional, indicates whether key has been found. Does   !
      !            not care if %param is empty                                 !
      !------------------------------------------------------------------------!

      implicit none

      ! lpl: InOut
      type(KEY_PARAM_PAIR), intent(inout) :: key_param
      ! lpl: Output
      logical, optional, intent(out) :: found
      ! lpl: Input
      character(len=*), intent(in) :: line_in

      ! lpl: Local variables
      integer :: ipos, jpos, key_len
      logical :: read_param
      character(len=16) :: delimiter

      if ( len_trim(key_param%key) == 0 ) &
           call utils_abort(' ERROR: No parameter to be read')

      if ( len_trim((adjustl(key_param%delimiter))) <= 0 ) then
         read_param = .false.
      else
         read_param = .true.
         delimiter = trim(adjustl(key_param%delimiter))
      end if

      ! lpl: 'index' returns pos from 1 to length of string, 0
      !      if not found
      ipos = index( line_in,trim(key_param%key) )
      ! lpl: Return if key is not found
      if ( ipos == 0 ) then
         found = .false.
         return
      ! lpl: Read stuff if key matches
      else
         if ( read_param ) then
            key_len = len_trim( adjustl(key_param%key) )
            ipos = ipos + key_len
            ! lpl: Search specified delimiter string first
            jpos = index( line_in(ipos:), trim(adjustl(delimiter)) )
            if ( jpos == 0 ) then
               ! lpl: Assume space if specified delimiter is not found
               delimiter = ' '
               jpos = index( line_in(ipos+key_len:),trim(adjustl(delimiter)) )
               if ( jpos == 0 ) then
                  call utils_abort(' ERROR reading delimiter')
               end if
            end if

            ipos = ipos + jpos + len_trim(adjustl(delimiter))
            if ( ipos > len(line_in) ) &
                 call utils_abort(' ERROR ipos > size(line_in)')
            key_param%param = trim( adjustl(line_in(ipos:)) )
         end if
         found = .true.
         if ( key_param%initialized ) then
            call utils_abort(' ERROR in internal_read_parameter: attemptting &
                 &to overwrite an initialized parameter')
         else
            key_param%initialized = .true.
         end if
      end if

    end subroutine internal_read_parameter

    logical function key_param_match(first_key_param,second_key_param)

      !------------------------------------------------------------------------!
      ! Simple comparison function for keys in the KEY_PARAM_PAIR object.      !
      !------------------------------------------------------------------------!
      ! Input : first_key_param, second_key_param - pairs of KEY_PARAM_PAIR    !
      !            objects to be compared                                      !
      !------------------------------------------------------------------------!

      implicit none

      type(KEY_PARAM_PAIR), intent(in) :: first_key_param, second_key_param

      key_param_match = ( trim(adjustl(first_key_param%key)) == &
           trim(adjustl(second_key_param%key)) )

    end function key_param_match

    integer function internal_format_line(line_inout,comment_char)

      !------------------------------------------------------------------------!
      ! Simple routine to check and remove comment lines. Returns < -ipos with !
      ! ipos being the position of the comment character if found, else        !
      ! returns the trimmed length (>=0) of 'line_inout'                       !
      !------------------------------------------------------------------------!
      ! InOut : line_inout - line string                                       !
      ! Input : comment_char - single character comment symbol                 !
      !------------------------------------------------------------------------!

      implicit none

      character(len=*), intent(inout) :: line_inout
      character(len=1), intent(in) :: comment_char

      integer :: ipos

      ipos = scan(line_inout,comment_char)

      if ( ipos > 0 ) then
         line_inout = line_inout(1:ipos-1)
         internal_format_line = -1*ipos
      else
         internal_format_line = len_trim(line_inout)
      end if

    end function internal_format_line

  end subroutine ddec_read_atom

  !===========================================================================!

  subroutine ddec_write_dens(dens,r,npts,output_filename,proc)

    !--------------------------------------------------------------------------!
    ! lpl: Print radial profiles for rad_dens                                  !
    !--------------------------------------------------------------------------!
    ! Input : dens - Array containing density                                  !
    !       : r - Array containing radial positions for the density. If array  !
    !            dimension is r(1), this is 'dr' instead.                      !
    !       : npts - Number of points from 'dens' and 'r' to be written out    !
    !       : output_filename - Filename for density                           !
    !       : proc - Optional. Specifies which procs does the writing          !
    !--------------------------------------------------------------------------!

    use utils, only: utils_abort, utils_unit, utils_close_unit_check, &
         utils_open_unit_check
    use comms, only: pub_my_proc_id, pub_root_proc_id, pub_total_num_procs

    implicit none

    ! lpl: Inputs
    real(kind=DP), intent(in) :: dens(:), r(:)
    integer, intent(in) :: npts
    character(len=*), intent(in) :: output_filename
    integer, optional :: proc ! Optionally specify proc from which to write

    ! lpl: Local variables
    integer :: it, output, ierr, write_proc

    ierr = 0

    if( present(proc) ) then
       if ( proc >= pub_total_num_procs ) &
            call utils_abort(' ERROR: Invalid proc in ddec_write_dens')
       write_proc = proc
    else
       write_proc = pub_root_proc_id
    end if

    if ( (size(r) < size(dens)) .and. (size(r) /= 1) ) &
       call utils_abort(' ERROR: size(r) < size(dens) in ddec_write_dens')

    if(pub_my_proc_id == write_proc) then
       output = utils_unit()
       open(unit=output,form="formatted", &
            file=(trim(adjustl(output_filename))),action="write",iostat=ierr)
       call utils_open_unit_check('ddec_write_dens','output',ierr)
       write(output,'(a14,1x,a20)') 'r (bohr)','density (a.u.)'
       if ( size(r) == 1 ) then
          do it=1,npts
             write(output,'(f14.7,1x,e20.10,1x,I5)') &
                  (real(it,kind=DP) - 0.5_DP)*r, dens(it)
          end do
       else
          do it=1,npts
             write(output,'(f14.7,1x,e20.10,1x,I5)') r(it), dens(it)
          end do
       end if
       close(unit=output,iostat=ierr)
       call utils_close_unit_check('ddec_write_dens','output',ierr)
    end if

  end subroutine ddec_write_dens

  !===========================================================================!

  ! lpl: COMMENT#06 - Since this only works for nonperiodic systems, it is now
  !      removed to prevent unaware users from accidentally using it on
 !      periodic ones.
!  subroutine ddec_rmse(grid_fine,grid,elements,atom_charge,coulomb_grid)

    !--------------------------------------------------------------------------!
    ! lpl: Quick-and-dirty subroutine to calculate RMS V(r) difference between !
    !      atomic point charges and real density. Only works for nonperiodic   !
    !      systems. Will either be made to work for periodic systmes or        !
    !      be removed alltogether.                                             !
    !--------------------------------------------------------------------------!
    ! Output: coulomb_grid - Optional, stores electrostatic potential due to   !
    !            the point charges                                             !
    ! Input : grid_fine - Cartesian grid containing electrostatic V(r)         !
    !       : grid - GRID_INFO of the simulation supercell                     !
    !       : elements - ELEMENT of the current system                         !
    !       : atom_charge - Array containing atomic charges ordered according  !
    !            to iat                                                        !
    !--------------------------------------------------------------------------!

!     use cell_grid, only: GRID_INFO, cell_grid_deposit_box, &
!          cell_grid_extract_box
!     use comms, only: comms_barrier, comms_bcast, comms_reduce, pub_on_root, &
!          pub_my_proc_id, pub_root_proc_id, pub_total_num_procs
!     use constants, only: ANGSTROM, HARTREE_IN_EVS, PI, stdout
!     use geometry, only: geometry_magnitude, POINT, unit_vector, &
!          operator(.DOT.), operator(.CROSS.), operator(*), operator(+), &
!          operator(-)
!     use rundat, only: pub_ddec_rmse_vdW
!     use simulation_cell, only: pub_cell
!     use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check
!     use visual, only: visual_scalarfield
!     use ion, only: ELEMENT
!     use parallel_strategy, only: par=>pub_par


!     implicit none

!     real(kind=DP), parameter :: hartree_in_kcalmol = 627.509469_DP

!     ! lpl: Output
!     real(kind=DP), intent(out) :: coulomb_grid(:,:,:,:)

!     ! lpl: Input
!     real(kind=DP),intent(in) :: grid_fine(:,:,:,:)
!     type(GRID_INFO), intent(in) :: grid
!     type(ELEMENT), intent(in) :: elements(par%nat)
!     real(kind=DP), intent(in) :: atom_charge(par%nat)

!     ! lpl: Local variables
!     real(kind=DP), allocatable :: buffer(:,:,:)
!     integer, allocatable :: grid_label(:,:,:)

!     integer :: iproc, iat, isp, it1, it2, it3

!     type(POINT) :: grid_ua1, grid_ua2, grid_ua3
!     real(kind=DP) :: grid_lda1, grid_lda2, grid_lda3
!     real(kind=DP) :: sin_12, sin_23, sin_13, vdW_max
!     integer :: n1_half, n2_half, n3_half, l1_ctr, l2_ctr, l3_ctr
!     integer :: l1_start, l2_start, l3_start, l1_end, l2_end, l3_end

!     type(POINT) :: v_atom, v_grid
!     real(kind=DP) :: r_dist

!     integer :: n_coul
!     real(kind=DP) :: v_coul, rmse, dV, dV_disp, rmse0, dV0, dV_abs

!     real(kind=DP), allocatable :: vdW(:,:)
!     character(len=4) :: species_label

!     real(kind=DP) :: rp
!     integer :: p1, p2, p3

!     integer :: ierr

! if (pub_debug) then
!        if ( pub_on_root ) write(stdout,'(a)') &
!             ' DEBUG: Entering ddec_rmse'
! end if

!     allocate(buffer(grid%ld1,grid%ld2,grid%max_slabs12),stat=ierr)
!     call utils_alloc_check('ddec_rmse','buffer',ierr)
!     allocate(grid_label(grid%ld1,grid%ld2,grid%num_my_slabs12),stat=ierr)
!     call utils_alloc_check('ddec_rmse','grid_label',ierr)
!     allocate(vdW(par%num_species,2),stat=ierr)
!     call utils_alloc_check('ddec_rmse','vdW',ierr)

!     ! lpl: Read vdW radii block
!     if (pub_on_root) then
!        write(stdout,'(a)') '============================================'
!        write(stdout,'(a)') '          vdW Details for RMS V(r)'
!        write(stdout,'(a)') '--------------------------------------------'
!        write(stdout,'(a)') ' isp   label vdW_min         vdW_max'
!        write(stdout,'(a)') '--------------------------------------------'
!        do isp=1,par%num_species
!           read(pub_ddec_rmse_vdW(isp),*) species_label, vdW(isp,1), vdW(isp,2)
!           write(stdout,'(1x,I4,1x,a4,2(2x,f14.7))') isp, &
!                trim(adjustl(species_label)),vdW(isp,1),vdW(isp,2)
!        end do
!        write(stdout,'(a)') '--------------------------------------------'
!     end if
!     call comms_barrier
!     call comms_bcast(pub_root_proc_id, vdW(:,:))

!     grid_label = 0

!     ! lpl: Grid unit vectors
!     grid_ua1 = unit_vector(grid%da1)
!     grid_ua2 = unit_vector(grid%da2)
!     grid_ua3 = unit_vector(grid%da3)

!     grid_lda1 = geometry_magnitude(grid%da1)
!     grid_lda2 = geometry_magnitude(grid%da2)
!     grid_lda3 = geometry_magnitude(grid%da3)

!     ! lpl: Calculate minimum box dimension required based on
!     !      isa_rcut_max
!     sin_12 = sin(acos(grid_ua1.DOT.grid_ua2))
!     sin_23 = sin(acos(grid_ua2.DOT.grid_ua3))
!     sin_13 = sin(acos(grid_ua1.DOT.grid_ua3))

!     vdW_max = maxval(vdW(:,2))

!     ! lpl: Set 'box size' surrounding atom
!     n1_half = int(ceiling(max(vdW_max/sin_12,vdW_max/sin_13)/grid_lda1)) + 1
!     n2_half = int(ceiling(max(vdW_max/sin_12,vdW_max/sin_23)/grid_lda2)) + 1
!     n3_half = int(ceiling(max(vdW_max/sin_13,vdW_max/sin_23)/grid_lda3)) + 1

!     do iat=1,par%nat
!        isp = elements(par%orig_atom(iat))%species_number
!        ! lpl: Absolute Cartesian atom coord
!        v_atom = elements(par%orig_atom(iat))%centre
!        ! lpl: Get cell start index number
!        l1_ctr = 1 + (v_atom.DOT.grid_ua1)/grid_lda1
!        l2_ctr = 1 + (v_atom.DOT.grid_ua2)/grid_lda2
!        l3_ctr = 1 + (v_atom.DOT.grid_ua3)/grid_lda3

!        l1_start = l1_ctr - n1_half
!        l2_start = l2_ctr - n2_half
!        l3_start = l3_ctr - n3_half
!        if ( l1_start <= 0 .or. l2_start <= 0 .or. l3_start <= 0 ) &
!           call utils_abort(' ERROR: atom box out of bounds')
!        l3_start = max(l3_start - grid%first_slab12(pub_my_proc_id) + 1,1)

!        l1_end = l1_ctr + n1_half
!        l2_end = l2_ctr + n2_half
!        l3_end = l3_ctr + n3_half
!        if ( l1_end > grid%ld1 .or. l2_end > grid%ld2 &
!             .or. l3_end > grid%ld3 ) call utils_abort(' ERROR: atom box out &
!             &of bounds in ddec_rmse')
!        l3_end = min(l3_end - grid%first_slab12(pub_my_proc_id) + 1, &
!             grid%num_my_slabs12)

!        ! lpl: Iterate over local range
!        do it1=l1_start,l1_end
!           do it2=l2_start,l2_end
!              do it3=l3_start,l3_end
!                 if(it3 > l3_end) cycle
!                 v_grid = (it1-1)*grid%da1 + (it2-1)*grid%da2 &
!                      + (it3-1 + grid%first_slab12(pub_my_proc_id) - 1)*grid%da3
!                 r_dist = geometry_magnitude(v_grid - v_atom)
!                 ! lpl: Mark grid points that falls within ESP evaluation region
!                 !      while avoiding regions already considered too close to
!                 !      any atom
!                 if (r_dist >= vdW(isp,1) .and. r_dist < vdW(isp,2) &
!                      .and. grid_label(it1,it2,it3) >= 0) then
!                    grid_label(it1,it2,it3) = grid_label(it1,it2,it3) + 1
!                 ! lpl: Exclude regions too close to any atom (< vdW radius)
!                 else if (r_dist < vdW(isp,1)) then
!                    grid_label(it1,it2,it3) = -1
!                 end if
!              end do
!           end do
!        end do
!     end do
!     call comms_barrier

!     ! lpl: Calculate RMS V(r)
!     rmse    = 0.0_DP
!     dV_disp = 0.0_DP
!     dV_abs  = 0.0_DP
!     rmse0   = 0.0_DP
!     dV0     = 0.0_DP
!     n_coul  = 0
!     do it1=1,grid%ld1
!        do it2=1,grid%ld2
!           do it3=1,grid%num_my_slabs12
!              ! lpl: Calculate point charge potential
!              if (grid_label(it1,it2,it3) > 0) then
!                 v_grid = (it1-1)*grid%da1 + (it2-1)*grid%da2 &
!                     + (it3-1 + grid%first_slab12(pub_my_proc_id) - 1)*grid%da3
!                 v_coul = 0.0_DP
!                 n_coul = n_coul + 1

!                 do iat=1,par%nat
!                    v_atom = elements(par%orig_atom(iat))%centre
!                    r_dist = geometry_magnitude(v_grid - v_atom)
!                    ! lpl: Dumbest way to account for periodicity. Not actually
!                    !      used here
!                    do p1=0,0
!                       do p2=-0,0
!                          do p3=-0,0
!                             rp = geometry_magnitude(real(grid%n1*p1,kind=DP)*&
!                                  grid%da1 + real(grid%n2*p2,kind=DP)*grid%da2 &
!                                  + real(grid%n3*p3,kind=DP)*grid%da3)
!                             v_coul = v_coul + atom_charge(iat)/(r_dist + rp)
!                          end do
!                       end do
!                    end do
!                 end do

!                 ! lpl: Accumulate point charge V(r)
!                 coulomb_grid(it1,it2,it3,1) = v_coul

!                 dV = grid_fine(it1,it2,it3,1) + v_coul
!                 dV_disp = dV_disp + dV
!                 dV_abs  = dV_abs + abs(dV)
!                 rmse = rmse + dV*dV

!                 dV0 = dV0 + grid_fine(it1,it2,it3,1)
!                 rmse0 = rmse0 + grid_fine(it1,it2,it3,1)*&
!                      grid_fine(it1,it2,it3,1)
!              end if
!           end do
!        end do
!     end do
!     call comms_barrier

!     call comms_reduce('SUM',n_coul)
!     call comms_reduce('SUM',rmse)
!     call comms_reduce('SUM',dV_disp)
!     call comms_reduce('SUM',dV_abs)
!     call comms_reduce('SUM',rmse0)
!     call comms_reduce('SUM',dV0)

!     dV_disp = dV_disp/n_coul
!     dV_abs  = dV_abs/n_coul
!     rmse    = rmse/n_coul
!     dV0   = dV0/n_coul
!     rmse0 = rmse0/n_coul

!     if(pub_on_root) then
!        write(stdout,'(a)') '    RMS Electrostatic Potential Summary'
!        write(stdout,'(a)') '--------------------------------------------'
!        write(stdout,'(a,I10)')   ' npts  = ', n_coul
!        write(stdout,'(a,f14.7,a)') ' RMS   = ', hartree_in_kcalmol*sqrt(rmse), &
!             ' kcal/mol'
!        write(stdout,'(a,f14.7,a)') ' dVerr = ', hartree_in_kcalmol*dV_disp, &
!             ' kcal/mol'
!        write(stdout,'(a,f14.7,a)') ' dVabs = ', hartree_in_kcalmol*dV_abs, &
!             ' kcal/mol'
!        write(stdout,'(a,f14.7,a)') ' RMSE  = ', &
!             hartree_in_kcalmol*sqrt( rmse - dV_disp*dV_disp ),' kcal/mol'
!        write(stdout,'(a,f14.7,a)') ' RMSE0 = ', &
!             hartree_in_kcalmol*sqrt( rmse0 - dV0*dV0 ),' kcal/mol'
!        write(stdout,'(a)') '============================================'
!     end if

!     deallocate(grid_label,stat=ierr)
!     call utils_dealloc_check('ddec_rmse','grid_label',ierr)
!     deallocate(buffer,stat=ierr)
!     call utils_dealloc_check('ddec_rmse','buffer',ierr)
!     deallocate(vdW,stat=ierr)
!     call utils_dealloc_check('ddec_rmse','vdW',ierr)

! if (pub_debug) then
!        if ( pub_on_root ) write(stdout,'(a)') &
!             ' DEBUG: Leaving ddec_rmse'
! end if

!   end subroutine ddec_rmse

  !===========================================================================!

  subroutine ddec_parse_string(delim_pos,delim_num,input_string,delim_set)


    implicit none

    ! lpl: Output
    integer, intent(out) :: delim_pos(:)
    integer, intent(out) :: delim_num

    ! lpl: Input
    character(len=*), intent(in) :: input_string
    character(len=1), intent(in) :: delim_set(:)

    ! lpl: Local variables
    character(len=1) :: char_buff
    integer :: str_len, ipos, max_fields

    integer :: ierr

    ierr = 0

    ! lpl: Length of input string
    str_len = len(input_string)

    ! lpl: Max num of split fields found will depend on length
    !      of this array
    max_fields = size(delim_pos)

    ! lpl: Do nothing if either of these conditions is true
    if ( size(delim_set) <= 0 .or. max_fields <= 0 .or. &
         str_len <= 0 ) return

    ! lpl: Since this returns the starting index for the next hypotehtical
    !      field, which is at index 'i', the end of the previous field should
    !      be 'i-2' to account for the delimiter preceeding the next field
    delim_pos = str_len + 2
    ipos = 1
    delim_num = 0
    ! If 1st character is already data
    char_buff = input_string(ipos:ipos)
   if ( all( delim_set /= char_buff ) ) then
         delim_num = delim_num + 1
         delim_pos(delim_num) = ipos
    end if

    do while ( ipos < str_len )
       char_buff = input_string(ipos:ipos)
       ! If delimiter is found
       if ( any( delim_set == char_buff ) ) then
          ! Iterate until non-delimiter is found
          ipos = ipos + 1
          do while ( ipos <= str_len )
             char_buff = input_string(ipos:ipos)
             if ( all( delim_set /= char_buff ) ) exit
             ipos = ipos + 1
          end do

          ! Stop parsing if max_fields reached
          if ( delim_num >= max_fields ) exit
          ! Store the index of the start of the next non-delimeter
          delim_num = delim_num + 1
          delim_pos(delim_num) = ipos
       else
          ipos = ipos + 1
       end if
    end do

  end subroutine ddec_parse_string

  !===========================================================================!

  subroutine ddec_1d_interp(final_dens,final_grid,initial_dens,initial_dr, &
       initial_offset,use_services)

    !--------------------------------------------------------------------------!
    ! lpl: Wrapper for transferring radial density from one grid to another    !
    !      via density interpolation. Used in 'ddec_read_atom' and             !
    !      'ddec_solve_atom'. Input grid must be uniform and starts from       !
    !      'initial_offset', output grid can be anything                       !
    !--------------------------------------------------------------------------!
    ! Output: final_dens - radial density output                               !
    ! Input : initial_dens - radial density input. Must be on a uniform grid   !
    !            that starts from r=0 and with increment of 'initial_dr'       !
    !       : final_grid - radial coordinates for 'final_dens'. Can have       !
    !            variable spacing.                                             !
    !       : initial_dens - radial density input. Must be on a uniform grid   !
    !            that starts from r=0.5*initial_dr                             !
    !       : initial_dr - constant dr for 'initial_dens'                      !
    !--------------------------------------------------------------------------!

    use utils, only: utils_abort
    use services, only: services_1d_interpolation

    implicit none

    ! lpl: Output
    real(kind=DP), intent(out) :: final_dens(:)
    ! lpl: Input
    real(kind=DP), intent(in) :: final_grid(:)
    real(kind=DP), intent(in) :: initial_dens(:), initial_dr
    real(kind=DP), intent(in) :: initial_offset
    logical, intent(in) :: use_services

    ! lpl: Local variables
    real(kind=DP) :: ir, ir_offset
    integer :: it, final_npts, initial_npts

    initial_npts = size(initial_dens)
    final_npts = size(final_dens)
    final_dens = 0.0_DP

    ! lpl: Offset of r in the intial density
    ir_offset = initial_offset/initial_dr

    if ( ir_offset > final_grid(final_npts) ) call utils_abort(' ERROR &
         &in ddec_1d_interp: ir_offset > final_grid(final_npts)')
    if ( initial_offset < 0 ) call utils_abort(' ERROR in ddec_1d_interp: &
         &initial_rad_offset or initial_offset < 0')

    if ( use_services ) then
       if ( ir_offset > 1.0_DP ) call utils_abort(' ERROR in &
            &ddec_1d_interp: ir_offset > 1.0_DP')
       ! lpl: Perform regular interpolation assuming smoothness at r=0
       do it=1,final_npts
          ! lpl: Radial grid of angmom = 0, so f(r) = f(-r)
          !      Essentially ir = abs(final_grid(it)/initial_dr) - ir_offset
          ir = ( abs(final_grid(it)) - initial_offset )/(1.0_DP*initial_dr)
          if ( int(ceiling(ir)) > initial_npts ) exit
          final_dens(it) = services_1d_interpolation(initial_dens, &
               initial_npts,ir,0)
       end do
    else
       ! lpl: Perform regular interpolation without assuming smoothness
       !      at r=0. In reality, this is better approximated by an
       !      exponentially decaying function, but that is left for
       !      later....
       do it=1,final_npts
          ! lpl: Radial grid of angmom = 0, so f(r) = f(-r)
          !      Essentially ir = abs(final_grid(it)/initial_dr) - ir_offset
          ir = ( abs(final_grid(it)) - initial_offset )/(1.0_DP*initial_dr)
          if ( int(ceiling(ir)) > initial_npts ) exit
          final_dens(it) = internal_1d_interpolation(initial_dens, &
               initial_npts,ir,.true.)
       end do
    end if

    contains

  !===========================================================================!

      real(kind=DP) function internal_1d_interpolation(values,num_values, &
           &xx,snap)

        !===============================================================!
        ! Interpolation in one dimension. Asuumes x for values(1) is 0  !
        !---------------------------------------------------------------!
        ! Plagiarized from services_mod. Does not smooth function at 0  !
        !===============================================================!

        use constants, only: DP
        use utils, only: utils_abort

        implicit none

        ! Arguments
        integer, intent(in) ::  num_values
        real(kind=DP), intent(in) :: values(num_values)
        real(kind=DP), intent(in) :: xx
        logical, intent(in) :: snap ! Snap to grid

        ! Local Variables
        integer :: start_index
        real(kind=DP) :: x2
        real(kind=DP) :: off, f1, f2, f3, f4, t0, t1, t2, t3
        real(kind=DP), parameter :: snap_thresh = 5e-5_DP
        integer :: ioff

        if ( num_values < 4 ) call utils_abort(' ERROR in &
             &internal_1d_interpolation: num_values < 4')

        ! lpl: Obtain the central index and the offset from it
        !      If start_index < 1, make it 2, and we are now extrapolating
        !      If start index > (num_values-2), then set the central index
        !      x2 to correspond to (num_values-2), and extrapolate forward
        start_index = min( max(int(xx) + 1,2), num_values-2 )
        x2 = real(start_index-1,kind=DP)
        off = xx - x2
        if ( snap ) then
           ioff = nint(off)
           if ( abs(ioff - off) < snap_thresh ) off = ioff
        end if

        f1 = values(start_index-1)
        f2 = values(start_index)
        f3 = values(start_index+1)
        f4 = values(start_index+2)

        t0 = f2
        t1 = ((6.0_dp*f3)-(2.0_dp*f1)-(3.0_dp*f2)-f4)/6.0_dp
        t2 = (f1+f3-(2.0_dp*f2))/2.0_dp
        t3 = (f4-f1+(3.0_dp*(f2-f3)))/6.0_dp

        internal_1d_interpolation = t0+off*(t1+off*(t2+off*t3))

      end function internal_1d_interpolation

  end subroutine ddec_1d_interp

  !===========================================================================!

  subroutine ddec_zero_threshold(density,threshold,mode,warn)

    use utils, only: utils_abort

    implicit none

    real(kind=DP), intent(inout) :: density(:,:,:)
    real(kind=DP), intent(in) :: threshold
    character(len=*), intent(in) :: mode
    logical, intent(in) :: warn

    integer :: it1, it2, it3
    real(kind=DP) :: corrected_val

    select case (mode)
       case ('Zero')
          corrected_val = 0.0_DP
       case ('Threshold')
          corrected_val = threshold
       case default
          call utils_abort(' ERROR in ddec_zero_threshold: Invalid mode flag')
    end select

    do it1=lbound(density,1),ubound(density,1)
       do it2=lbound(density,2),ubound(density,2)
          do it3=lbound(density,3),ubound(density,3)
             if ( density(it1,it2,it3) < threshold ) then
                density(it1,it2,it3) = corrected_val
             end if
          end do
       end do
    end do

  end subroutine ddec_zero_threshold

  !===========================================================================!

  subroutine ddec_c3_dens(ref_dens,min_ion_state,max_ion_state)

    !===============================================================!
    ! lpl: Reshapes input 'c2' density to 'c3 density               !
    !---------------------------------------------------------------!
    !===============================================================!


    use comms, only: pub_on_root
    use constants, only: stdout
    use utils, only: utils_abort

    implicit none

    ! lpl: InOut
    type(DDEC_REF_DENS), intent(inout) :: ref_dens ! For one species
    integer, intent(in) :: min_ion_state, max_ion_state

    ! lpl: Local variables
    real(kind=DP) :: max_dens, sum_dens, &
         sum_dens_diff, dens_diff, min_dens_diff
    ! lpl: Threshold used in CHARGEMOL (decay not enforced for densities
    !      smaller than this)
    real(kind=DP), parameter :: min_dens_thresh = 1e-16_DP
    integer, parameter :: max_iit = 100
    integer :: iion_limit(-1:1)
    integer :: it, iit, iion, ion_sign
    logical :: tail

    if(pub_debug_on_root) write(stdout,'(a)') ' DEBUG: Entering ddec_c3_dens'

    ! lpl: We require neutral density as all monotonic decay across ions are
    !      performed with respect to the neutral density - even if we are
    !      generating extra non-neutral densities at a later time, where we
    !      would require previously-generated c3 densities
    if ( .not. allocated(ref_dens%ion_dens(0)%dens) ) &
         call utils_abort(' ERROR in ddec_c3_dens: Neutral species has not&
              & been initialized')

    iion_limit(-1) = min_ion_state ! Max anion
    iion_limit( 1) = max_ion_state ! Max cation

    ! lpl: Determine where to start enforcing decay
    iion_limit(0)  = max_ion_state - min_ion_state
    if ( abs(iion_limit(0)) /= abs(max_ion_state) + abs(min_ion_state) ) then
       if ( abs(max_ion_state) < abs(min_ion_state) ) &
            call utils_abort(' ERROR in ddec_c3_dens: Invalid ionic range&
                 & specified #1')
       iion_limit(0) = sign(min_ion_state,iion_limit(0))
    else
       if ( max_ion_state < min_ion_state ) call utils_abort(' ERROR in &
            &ddec_c3_dens: Invalid ionic range specified #2')
       iion_limit(0) = 0
    end if

    ! lpl: Enforce monotonic decay in the -ve -and +ve directions
    do ion_sign=-1,1,2

       ! lpl: Enforce monotonic decay and accumulate differences in
       !      integrated density
       do iion=iion_limit(0),iion_limit(ion_sign),ion_sign

          ! lpl: MOD#07 - if certain ionic range was not initialized, print
          !      warning and invalidate all other ions above/below
          if ( .not. ref_dens%init_status(iion) ) then
             if ( pub_on_root ) then
                write(stdout,'(a,sp,I4.3,a)') ' WARNING: refdens '//&
                     ref_dens%element%species_id, iion, ' was not&
                     & initialized.'
                select case (ion_sign)
                   case (-1); write(stdout,'(a)',advance='no') ' Min'
                   case ( 1); write(stdout,'(a)',advance='no') ' Max'
                   case default
                      call utils_abort(' ERROR in ddec_c3_dens: Invalid&
                           & ion_sign')
                end select
                write(stdout,'(a,sp,I4.3)') ' ion range possible is ', &
                     iion - ion_sign
             end if

             ! lpl: Invalide this ions and all others below/above it
             do it=iion,iion_limit(ion_sign),ion_sign
                ref_dens%init_status(it) = .false.
             end do

             ! lpl: Exit iion loop
             exit
          end if

          if ( .not. allocated(ref_dens%ion_dens(iion)%dens) ) &
               call utils_abort(' ERROR in ddec_c3_dens: Ionic density&
                    & not initialized')

          max_dens = min_dens_thresh
          sum_dens_diff = 0.0_DP  ! summed_oxidation_density in CHARGEMOL
          sum_dens = 0.0_DP       ! summed_delta_oxidation_density  "

          ! lpl: Integrate unaltered density
          do it=2,ref_dens%npts
             sum_dens = sum_dens + &
                  ref_dens%ion_dens(iion)%dens(it)*ref_dens%rad_shell_vol(it)
          end do

          if (pub_debug_on_root) write(stdout,'(a,1x,e20.10)') &
               ' sum_dens_0 =', sum_dens

          ! lpl: CHARGEMOL update - correct for density at cusp by adjusting
          !      so that its volume contribution is whaterver is left unsummed
          ref_dens%ion_dens(iion)%dens(1) = &
               max( (ref_dens%element%atomic_number - iion - sum_dens)/&
               ref_dens%rad_shell_vol(1), min_dens_thresh )
          sum_dens = sum_dens + &
               ref_dens%ion_dens(iion)%dens(1)*ref_dens%rad_shell_vol(1)

          if (pub_debug_on_root) write(stdout,'(a,1x,e20.10)') &
               ' sum_dens_1 =', sum_dens

          ! lpl: Enforce monotonic decay
          tail = .true.
          do it=ref_dens%npts,2,-1
             if ( ref_dens%ion_dens(iion)%dens(it) > max_dens ) then
                max_dens = ref_dens%ion_dens(iion)%dens(it)
                tail = .false.
             else
                ! lpl: Skip if we are still at the tail section where density
                !      is likely zero
                if ( tail ) cycle
                ! lpl: Integrate difference between altered and unaltered
                !      densities
                sum_dens_diff = sum_dens_diff + (max_dens - &
                     ref_dens%ion_dens(iion)%dens(it))*&
                     ref_dens%rad_shell_vol(it)
                ref_dens%ion_dens(iion)%dens(it) = max_dens
             end if
          end do

          ! lpl: Rescale altered density to correct charge
          ref_dens%ion_dens(iion)%dens = ref_dens%ion_dens(iion)%dens*&
               (sum_dens/(sum_dens + sum_dens_diff))

          ! lpl: Enforce monotonic decay between ionic densities
          if ( iion /= 0 ) then

             ! lpl: Make sure that the previous density is c3
             if ( .not. ref_dens%c3flag(iion-ion_sign) ) then
                call utils_abort(' ERROR in ddec_c3_dens: Attempting to&
                     & generate c3 density by building upon other non-c3&
                     & densities')
             else if ( ref_dens%c3flag(iion) ) then
                call utils_abort(' ERROR in ddec_c3_dens: Attempting to&
                     & generate c3 density by building upon a c3 density')
             end if

             do iit=1,max_iit
                sum_dens_diff = 0.0_DP
                min_dens_diff = 0.0_DP
                do it=ref_dens%npts,1,-1
                   dens_diff = ref_dens%ion_dens(iion-ion_sign)%dens(it) - &
                        ref_dens%ion_dens(iion)%dens(it)
                   min_dens_diff = max( ion_sign*dens_diff,min_dens_diff )
                   dens_diff = dens_diff - ion_sign*min_dens_diff
                   ref_dens%ion_dens(iion)%dens(it) = &
                        ref_dens%ion_dens(iion)%dens(it) + dens_diff

                   sum_dens_diff = sum_dens_diff + &
                        dens_diff*ref_dens%rad_shell_vol(it)

                end do

                ! lpl: Rescale altered density to correct charge
                ref_dens%ion_dens(iion)%dens = ref_dens%ion_dens(iion)%dens*&
                     (sum_dens/(sum_dens + sum_dens_diff))

                if ( sum_dens_diff == 0.0_DP ) exit
             end do ! END do iit=1,max_iit
          end if ! END if ( iion /= 0 )

          ! lpl: Mark this density as having been initialized as 'c3'
          ref_dens%c3flag(iion) = .true.

       end do ! END do iion=0,iion_limit(ion_sign),ion_sign

    end do ! END do ion_sign=-1,1,2

    if (pub_debug_on_root) write(stdout,'(a)') ' DEBUG: Leaving ddec_c3_dens'

  end subroutine ddec_c3_dens

  !===========================================================================!

  subroutine ddec_shell_avg(rad_dens_out, dr_out, rad_dens_in, dr_in)

    ! lpl: Averages over radial shells: rad_dens_in --> rad_dens_out

    use utils, only: utils_abort

    implicit none

    ! Output
    real(kind=DP), intent(out) :: rad_dens_out(:)

    ! Input
    real(kind=DP), intent(in) :: rad_dens_in(:)
    real(kind=DP), intent(in) :: dr_out, dr_in

    ! Local variables
    real(kind=DP), parameter :: zero_threshold = 1e-15
    real(kind=DP) :: iin_frac_start, iin_frac_end, shell_vol
    integer :: iin_start, iin_end, iin, iout
    logical :: term

    call utils_abort(' ERROR in ddec_shell_avg: This subroutine is&
         & untested and should not be used (yet)')

    if ( lbound(rad_dens_out,1) /= 1 .or. lbound(rad_dens_in,1) /= 1 ) &
         call utils_abort(' ERROR in ddec_shell_avg:&
              & lbound(rad_dens_out) /= 1 or lbound(rad_dens_in) /= 1')

    if ( dr_out <= 0.0_DP .or. dr_in <= 0.0_DP ) &
         call utils_abort(' ERROR in ddec_shell_avg: dr_out <= 0&
              & or dr_in <= 0')

    iin_frac_end = 0.0_DP
    iin_end = 0
    term = .false.

    rad_dens_out = 0.0_DP

    do iout=1,size(rad_dens_out)

       ! Continue from previous input bin position
       iin_start = iin_end
       iin_frac_start = iin_frac_end

       ! Sum remainder of current input bin that straddles boundary of
       ! output bin between 'iout-1' and 'iout'
       if ( iin_frac_start > zero_threshold ) then
          ! (ir+1)**3 - (ir+1-frac)**3
          shell_vol = iin_frac_start*( 3*(iin_start + 1)*(iin_start + 1 &
               - iin_frac_start) + iin_frac_start**2 )
          rad_dens_out(iout) = rad_dens_out(iout) + &
               shell_vol*rad_dens_in(iin_start + 1)
       else
          iin_frac_start = 0.0_DP
          ! (ir+1)**3 - ir**3
          shell_vol = 3*iin_start*(iin_start + 1) + 1
          rad_dens_out(iout) = rad_dens_out(iout) + &
               shell_vol*rad_dens_in(iin_start + 1)
       end if

       ! Find out where current bin 'iout' terminates relative to input bin
       iin_frac_end = 1.0_DP*iout*(dr_out/dr_in)
       ! Input bin index less 1 that straddles 'iout' boundary
       iin_end = aint(iin_frac_end)
       ! If input bin is the last, then set limit to 'iin_end'
       if ( iin_end + 1 > size(rad_dens_in) ) then
          iin_end = size(rad_dens_in) - 1
          term = .true.
       end if

       ! Sum non-boundary input bins that are completely within bin 'iout'
       do iin=iin_start+1,iin_end-1
          ! (ir+1)**3 - ir**3
          shell_vol = 3*iin*(iin + 1) + 1
          rad_dens_out(iout) = rad_dens_out(iout) + &
               shell_vol*rad_dens_in(iin + 1)
       end do

       ! Sum last input bin that straddles boundary of 'iout'
       iin_frac_end = iin_frac_end - iin_end
       if ( iin_frac_end > zero_threshold ) then
          ! (ir+eta)**3 - ir**3
          shell_vol = iin_frac_end*( 3*iin_end*(iin_frac_end + iin_end) + &
               iin_frac_end**2 )
          rad_dens_out(iout) = rad_dens_out(iout) + &
               shell_vol*rad_dens_in(iin_end + 1)
       else
          iin_frac_end = 0.0_DP
          ! ir**3 - (ir-1)**3
          shell_vol = 3*iin_end*(iin_end - 1) + 1
          rad_dens_out(iout) = rad_dens_out(iout) + &
               shell_vol*rad_dens_in(iin_end + 1)
          ! If this last bin does not straddle boundary of 'iout', then advance
          ! 'iin_end' to the next bin so 'iin_start' will not re-sum this bin
          iin_end = iin_end + 1
       end if

       ! Average over summed volume
       shell_vol = (iin_end + iin_frac_end)**3 - (iin_start + iin_frac_start)**3
       rad_dens_out(iout) = rad_dens_out(iout)/shell_vol

       if ( term ) exit

    end do

  end subroutine ddec_shell_avg

  !===========================================================================!

  subroutine ddec_rad_to_box_to_rad(rad_dens_out, dr_out, rad_dens_in, dr_in, &
       rcut_out, norm_out, norm)

    ! lpl: Perform spherical averaging by 'boxifying' i.e. depositing
    !      radial densities into Cartesian grid and accumulating them
    !      back into a radial density.

    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Output
    real(kind=DP), intent(out) :: rad_dens_out(:)

    ! Input
    real(kind=DP), intent(in) :: rad_dens_in(:)
    real(kind=DP), intent(in) :: dr_out, dr_in, rcut_out, norm_out
    logical, intent(in) :: norm

    ! Local variables
    real(kind=DP), allocatable :: grid_val(:,:,:)
    integer, allocatable :: rad_npts(:)
    real(kind=DP) :: vol, sum_val, grid_dr, r_dist
    integer :: ir, it1, it2, it3
    integer :: grid_n, grid_ctr, npts_in, npts_out

    integer :: ierr

    if (pub_debug_on_root) write(stdout,'(a)') ' DEBUG: Entering&
         & ddec_rad_to_box_to_rad'

    ierr = 0

    if ( lbound(rad_dens_out,1) /= 1 .or. lbound(rad_dens_in,1) /= 1 ) &
         call utils_abort(' ERROR in ddec_rad_to_box_to_rad:&
              & lbound(rad_dens_out) /= 1 or lbound(rad_dens_in) /= 1')

    if ( dr_out <= 0.0_DP .or. dr_in <= 0.0_DP ) &
         call utils_abort(' ERROR in ddec_rad_to_box_to_rad: dr_out <= 0&
              & or dr_in <= 0')

    npts_in = size(rad_dens_in)
    npts_out = size(rad_dens_out)

    ! Determine grid finesse based on dr_out
    vol = FOUR_PI*(dr_out**3)/3.0_DP
    grid_dr = 0.5_DP*vol**(1.0_DP/3.0_DP)
    grid_n  = nint(2.0_DP*rcut_out/grid_dr) + 2
    grid_ctr = grid_n/2

    allocate(rad_npts(npts_out),stat=ierr)
    call utils_alloc_check('ddec_rad_to_box_to_rad','rad_npts',ierr)
    allocate(grid_val(grid_n,grid_n,grid_n),stat=ierr)
    call utils_alloc_check('ddec_rad_to_box_to_rad','grid_val',ierr)

    ! Deposit from 'rad_dens_in'
    grid_val = 0.0_DP
    do it1=1,grid_n
       do it2=1,grid_n
          do it3=1,grid_n

             r_dist = grid_dr*sqrt( real((it1 - grid_ctr)**2 + &
                  (it2 - grid_ctr)**2 + (it3 - grid_ctr)**2, kind=DP) )
             ir = aint(r_dist/dr_in) + 1

             if ( ir > npts_in ) cycle

             grid_val(it1,it2,it3) = rad_dens_in(ir)

          end do
       end do
    end do

    ! Extract to 'rad_dens_out'
    rad_dens_out = 0.0_DP
    sum_val = 0.0_DP
    rad_npts = 0
    do it1=1,grid_n
       do it2=1,grid_n
          do it3=1,grid_n

             r_dist = grid_dr*sqrt( real((it1 - grid_ctr)**2 + &
                  (it2 - grid_ctr)**2 + (it3 - grid_ctr)**2, kind=DP) )

             ir = aint(r_dist/dr_out) + 1

             if ( ir > npts_out ) cycle

             rad_dens_out(ir) = rad_dens_out(ir) + grid_val(it1,it2,it3)
             sum_val = sum_val + grid_val(it1,it2,it3)
             rad_npts(ir) = rad_npts(ir) + 1

          end do
       end do
    end do
    sum_val = (grid_dr**3)*sum_val

    ! Spherically-average 'rad_dens_out'
    if ( norm ) then
       sum_val = norm_out/sum_val
    else
       sum_val = 1.0_DP
    end if

    do ir=1,npts_out
       if ( rad_npts(ir) > 0 ) then
          rad_dens_out(ir) = sum_val*(rad_dens_out(ir)/rad_npts(ir))
       else
          rad_dens_out(ir) = 0.0_DP
       end if
    end do

    if (pub_debug_on_root) write(stdout,'(a,1x,e20.10)') ' sum_val =', &
         grid_dr**3*sum(rad_dens_out*rad_npts)

    deallocate(grid_val,stat=ierr)
    call utils_dealloc_check('ddec_rad_to_box_to_rad','grid_val',ierr)
    deallocate(rad_npts,stat=ierr)
    call utils_dealloc_check('ddec_rad_to_box_to_rad','rad_npts',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') ' DEBUG: Leaving&
         & ddec_rad_to_box_to_rad'

  end subroutine ddec_rad_to_box_to_rad

  !===========================================================================!

  subroutine ddec_copy_rad(rad_dens_out, dr_out, rad_dens_in, dr_in)

    ! lpl: Null routine that copies values

    use utils, only: utils_abort

    implicit none

    ! Output
    real(kind=DP), intent(out) :: rad_dens_out(:)

    ! Input
    real(kind=DP), intent(in) :: rad_dens_in(:)
    real(kind=DP), intent(in) :: dr_out, dr_in

    ! Local variables
    real(kind=DP) :: zero_threshold = 1e-5_DP
    integer :: n_min

    if ( lbound(rad_dens_out,1) /= 1 .or. lbound(rad_dens_in,1) /= 1 ) &
         call utils_abort(' ERROR in ddec_copy_rad: lbound(rad_dens_out) /= 1&
              & or lbound(rad_dens_in) /= 1')

    if ( dr_out <= 0.0_DP .or. dr_in <= 0.0_DP ) &
          call utils_abort(' ERROR in ddec_copy_rad: dr_out <= 0 or dr_in <= 0')

    if ( abs(dr_out - dr_in) > zero_threshold ) &
         call utils_abort(' ERROR in ddec_copy_rad: dr_out /= dr_in')

    rad_dens_out = 0.0_DP
    n_min = min(size(rad_dens_out),size(rad_dens_in))
    rad_dens_out(1:n_min) = rad_dens_in(1:n_min)

  end subroutine ddec_copy_rad

  !===========================================================================!

  subroutine ddec_check_3D_array(array, n1, n2, n3)

    ! lpl: Dumb function to check size of input arrays since some array
    !      dimensions are stored in objects instead of integer values
    !      themselves when passed to 'ddec_main'

    use utils, only: utils_abort

    implicit none

    real(kind=DP), intent(in) :: array(:,:,:)
    integer, intent(in) :: n1, n2, n3

    if ( n1 < 0 .or. n2 < 0 .or. n3 < 0 )  call utils_abort(' ERROR in&
         & ddec_check_3D_array: n1 < 0 .or. n2 < 0 .or. n3 < 0')

    if ( size(array,1) /= n1 .or. size(array,2) /= n2 .or. &
         size(array,3) /= n3 ) then
       write(stdout,'(a,3(1x,I5))') ' Array dim:', &
            size(array,1), size(array,2), size(array,3)
       write(stdout,'(a,3(1x,I5))') ' Check dim:', n1, n2, n3
       call utils_abort(' ERROR in ddec_check_3D_array:&
            & dimensions are incorrect')
    end if

  end subroutine ddec_check_3D_array

  !===========================================================================!
  subroutine ddec_off_site_charges(optimum_posit, optimum_charges, &
       error_before, minerror, atom_ctr_point, atom_ctr_point_array,&
       remaining_charge, aim_volume, ierr, aim_val, local_box, element_type, &
       all_element_types, number, grid_spacing_x, grid_spacing_y, &
       grid_spacing_z, no_charges, vdw_radius )

   !--------------------------------------------------------------------------!
   ! aeaa : This subroutine finds the optimal positions and charges for a set !
   !        of point charges in order to recreate the electrostatic potential !
   !        of an atom. The number of charges tested goes from one to a       !
   !        maximum of three. The position of the charges is based on the     !
   !        number of bonds the atom has.                                     !
   !--------------------------------------------------------------------------!
   ! Output: optimum_posit - Optimum positions of added point charges         !
   !       : optimum_charges - Optimum charges of added point charges         !
   !       : error_before - The error before extra charges are added          !
   !       : minerror - gives the error after charges are added               !
   !       : no_charges - the number of additional charges needed             !
   ! Input : atom_ctr_point - the position of the center of the grid used     !
   !       : remaining_charge - the DDEC charge the has been found            !
   !       : aim_val - contains the charge density values for the atom        !
   !       : local_box - contain info on the grid of the box                  !
   !       : element_type - the element of the atom, needed for vdw radius    !
   !       : all_element_types - all the elements of the molecule             !
   !       : number - the atom number                                         !
   !       : grid_spacing - necessary to find charge from charge density      !
   !--------------------------------------------------------------------------!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: DP, PI, stdout
    use geometry, only: POINT
    use model_type, only: MODEL
    use parallel_strategy, only: par=>pub_par
    use utils, only:  utils_alloc_check, utils_dealloc_check, utils_abort
    use rundat, only: pub_ddec_aniso_error_thres, pub_ddec_aniso_max_dis,&
         pub_ddec_aniso_max_dis_halogen, pub_ddec_aniso_error_reduce

    implicit none

    ! aeaa: Output Variables
    integer, intent(out) :: no_charges

    real(kind=DP), intent(out) :: optimum_posit(3,5), optimum_charges(5)
    real(kind=DP), intent(out) :: minerror, error_before, vdw_radius


    ! aeaa: Input Variables
    integer, intent(inout) :: ierr
    integer, intent(in) :: number !Atom number assigned

    type(POINT), intent(in) :: atom_ctr_point
    type(DDEC_LOCAL_BOX), intent(in) :: local_box

    real(kind=DP), intent(inout) :: aim_volume
    real(kind=DP), intent(in) :: remaining_charge ! DDEC charge values
    real(kind=DP), intent(in) :: aim_val(:,:,:)
    real(kind=DP), intent(in) :: grid_spacing_x, grid_spacing_y, &
         grid_spacing_z
    type(POINT), intent(in) :: atom_ctr_point_array(:)

    ! aeaa: Local variables
    ! aeaa: Variables to define sampling grid use for potential
    integer, parameter :: no_div = 15, no_shells = 4
    integer :: length !Give array length 0f sampling grid
    integer :: i, j, n, it1, it2, it3, iat ! Indices
    integer :: no_bonded, index

    ! aeaa: Variables to define sampling grid use for potential
    real(kind=DP) :: free_volume, free_radius, start_radius, final_radius
    real(kind=DP) :: shell_radius
    real(kind=DP) :: error_thres !Above this value extra charges used
    real(kind=DP) :: error, error_reduce !States when extra charges added
    real(kind=DP) :: max_distance ! vdW radius for atom
    real(kind=DP) :: atom_ctr(3) !center point of grid
    real(kind=DP) :: bond_a(3), bond_b(3), bond_c(3) ! Neighbouring bonds
    real(kind=DP) :: dir_vector(3), nor_vector(3), per_vector(3), &
         sec_dir_vector(3) ! Direction vectors
    real(kind=DP) :: tmp(3), tmp_i(3), tmp_j(3)
    real(kind=DP) :: chargenucleus ! Valence electron value
    real(kind=DP) :: true_V(no_shells, no_div, no_div) ! Real potential
    ! aeaa: Sphere_positions gives coords of sphere
    real(kind=DP) :: sphere_positions(no_shells, 3, no_div, no_div)
    ! aeaa: Sample_positions give the spheres centered at nucleus point
    real(kind=DP) :: sample_positions(no_shells, 3, no_div, no_div)
    real(kind=DP), allocatable  :: density(:), charges(:) !For true potential
    real(kind=DP), allocatable  :: positions(:,:) ! For true  potential
    real(kind=DP), allocatable  :: neighbouring_vectors(:,:)
    real(kind=DP), allocatable :: atom_distances(:,:) ! array of atom distance
    real(kind=DP), dimension(no_div,no_div) :: xsphere, ysphere, zsphere

    ! aeaa: Element types defined
    character(len=2) :: element_type, all_element_types(:), number_str
    integer :: neighbour_element_types(3)

    ! aeaa: The positions and charges for the additional point charges
    real(kind=DP) ::  rcentermin(3), r1min(3),r2min(3), r3min(3), r4min(3), &
         temp_vector(3), qcentermin, q1min,q2min, q3min, q4min

    ! aeaa: Change type of atom center from point type
    atom_ctr = [atom_ctr_point%X,  atom_ctr_point%Y, atom_ctr_point%Z]

    ! aeaa: Change values from Bohr to Ang
    atom_ctr = atom_ctr * 0.529177_DP;
    aim_volume = aim_volume * ((0.529177_DP) ** 3)

    element_type = trim(element_type)

    ! aeaa: Free radius, free volume and max_distance set
    select case (element_type)
    case ('H');
       free_radius = 1.44_DP !Free radius - same as Louis Lee paper
       free_volume = 7.60_DP
       max_distance = pub_ddec_aniso_max_dis
    case ('B');
       free_radius = 2.04_DP
       free_volume = 46.7_DP
       max_distance = pub_ddec_aniso_max_dis
    case ('C');
       free_radius = 1.93_DP
       free_volume = 34.4_DP
       max_distance = pub_ddec_aniso_max_dis
    case ('N');
       free_radius = 1.83_DP
       free_volume = 25.9_DP
       max_distance = pub_ddec_aniso_max_dis
    case ('O');
       free_radius = 1.75_DP
       free_volume = 22.1_DP
       max_distance = pub_ddec_aniso_max_dis
    case ('F');
       free_radius = 1.68_DP
       free_volume = 18.2_DP
       max_distance = pub_ddec_aniso_max_dis_halogen
    case ('P');
       free_radius = 2.07_DP
       free_volume = 84.6_DP
       max_distance = pub_ddec_aniso_max_dis
    case ('S');
       free_radius = 2.02_DP
       free_volume = 75.2_DP
       max_distance = pub_ddec_aniso_max_dis
    case ('Cl');
       free_radius = 1.97_DP
       free_volume = 65.1_DP
       max_distance = pub_ddec_aniso_max_dis_halogen
    case ('Br');
       free_radius = 2.10_DP
       free_volume = 95.7_DP
       max_distance = pub_ddec_aniso_max_dis_halogen
    case default
        call utils_abort( 'ERROR in ddec_aniso, atom with unknown vdw &
             &radius is being used.' )
    end select

    vdw_radius = free_radius

    ! aeaa: radius of shells used for sampling potential
    start_radius = 1.4_DP * vdw_radius
    final_radius = 2.0_DP * vdw_radius

    length = local_box%n1 * local_box%n2 * local_box%n3

    allocate(density(length),stat=ierr)
    call utils_alloc_check('ddec_aniso', 'density', ierr)
    allocate(charges(length))
    call utils_alloc_check('ddec_aniso', 'charges', ierr)
    allocate(positions(length, 3))
    call utils_alloc_check('ddec_aniso', 'positions', ierr)
    allocate(atom_distances(size(atom_ctr_point_array), &
         size(atom_ctr_point_array)))
    call utils_alloc_check('ddec_aniso', 'atom_distances', ierr)

    ! aeaa: Positions of sample points set
    i = 1
    do it1=1,local_box%n1
       do it2=1,local_box%n2
          do it3=1,local_box%n3
             !Set the positions coords, change to Ang and
             !Multiply by the grid width(0.225)
             positions(i, 1)  = real(it1 -  1, kind=DP) * &
                  (0.529177_DP * grid_spacing_x)
             positions(i, 2)  = real(it2 -  1, kind=DP) * &
                  (0.529177_DP * grid_spacing_y)
             positions(i, 3)  = real(it3 -  1, kind=DP) * &
                  (0.529177_DP * grid_spacing_z)

             density(i) = aim_val(it1,it2,it3)
             i = i + 1

          end do
       end do
    end do

    ! aeaa: Charges and total charges calculated
    charges = - density * grid_spacing_x * grid_spacing_y * grid_spacing_z
    chargenucleus = - ( sum(charges) - remaining_charge )

    ! aeaa: Calculates coordinates of sphere of radius 1.0
    call internal_sphere_coordinates( xsphere, ysphere, zsphere, no_div)

    ! aeaa: The sample positions - a number of spheres at varying radii
    do n =  1, no_shells
       shell_radius = start_radius + (( n - 1 ) * &
            (final_radius - start_radius) ) / ( real(no_shells, kind=DP) - 1 )

       ! aeaa: Set true potential sample positions
       do i =  1,no_div
          do j = 1,no_div
             sphere_positions(n,1,i,j) =  xsphere(i,j) * shell_radius
             sphere_positions(n,2,i,j) =  ysphere(i,j) * shell_radius
             sphere_positions(n,3,i,j) =  zsphere(i,j) * shell_radius

             sample_positions(n,1,i,j) =  xsphere(i,j) * shell_radius &
                  + atom_ctr(1)
             sample_positions(n,2,i,j) =  ysphere(i,j) * shell_radius &
                  + atom_ctr(2)
             sample_positions(n,3,i,j) =  zsphere(i,j) * shell_radius &
                  + atom_ctr(3)
          end do
       end do
    end do

    ! aeaa: Finds the potential at the sampling points
    ! aeaa: Length gives the number of points used
    do n  = 1, no_shells
       do i =  1,no_div
          do j = 1,no_div
             true_V(n,i,j)=  internal_V(sample_positions(n,:,i,j), &
                  length, positions, charges);
             ! aeaa: Adds the potential from nucleus
             true_V(n,i,j)=  true_V(n,i,j) + internal_V &
                  (sample_positions(n,:,i,j), 1, atom_ctr, [chargenucleus])
          end do
       end do
    end do

    ! aeaa: Above this error extra point charges are needed
    error_thres = pub_ddec_aniso_error_thres
    ! aeaa: Set the improvement needed for extra charges to be added
    error_reduce = pub_ddec_aniso_error_reduce

    ! aeaa: The error if no off center charges are used
    error =  internal_calc_error(1, [0.0_DP, 0.0_DP, 0.0_DP], &
         [remaining_charge], no_div, true_V, sphere_positions, no_shells)

    no_charges = 0
    error_before = error
    optimum_posit = 0.0_DP
    optimum_charges = 0.0_DP

    ! aeaa: Finds the distance between the different atoms
    do i = 1,size(atom_ctr_point_array)
       do j = 1,size(atom_ctr_point_array)
          tmp_i = [atom_ctr_point_array(i)%X , atom_ctr_point_array(i)%Y, &
               atom_ctr_point_array(i)%Z]
          tmp_j = [atom_ctr_point_array(j)%X , atom_ctr_point_array(j)%Y, &
               atom_ctr_point_array(j)%Z]
          tmp = tmp_i - tmp_j
          atom_distances(i,j) = ( internal_norm2(tmp) ) * 0.529177_DP
       end do
    end do

    no_bonded = 0

    ! aeaa: Counts the number of neighbouring atoms
    do i = 1,size(atom_ctr_point_array)
       if ( (atom_distances(number, i) < free_radius) .and. &
            (number /= i) ) then
          no_bonded = no_bonded + 1
       end if
    end do

    allocate(neighbouring_vectors(no_bonded,3))
    call utils_alloc_check('ddec_aniso', 'neighbouring_vectors', ierr)

    ! aeaa: Sets neighbouring vectors of atom
    index = 1
    do i = 1,size(atom_ctr_point_array)
       if ( (atom_distances(number, i) < free_radius) .and. &
            (number /= i) ) then
          tmp_i = [atom_ctr_point_array(number)%X, &
               atom_ctr_point_array(number)%Y, atom_ctr_point_array(number)%Z]
          tmp_j = [atom_ctr_point_array(i)%X , atom_ctr_point_array(i)%Y, &
               atom_ctr_point_array(i)%Z]
          neighbouring_vectors(index,:) = tmp_j - tmp_i
          neighbour_element_types(index) = i
          index = index + 1
       end if
    end do

    minerror = error_before

    ! aeaa: Extra points are now added based on bond number
    if ( error_before > error_thres ) then

       ! aeaa: One bond is along the bond axis
       if (no_bonded == 1) then
          neighbouring_vectors = -neighbouring_vectors ! Reverse direction

          call internal_symmetric_1_1_extra_charges(error, rcentermin, r1min, &
               qcentermin, q1min, remaining_charge, max_distance, true_V, &
               sphere_positions, no_shells, no_div, neighbouring_vectors)

          ! aeaa: Only add extra charges if suitable drop in error
          if ( ( ( error + error_reduce) < error_before ) ) then
             no_charges = 1
             optimum_posit(:, 1) = rcentermin
             optimum_posit(:, 2) = r1min
             optimum_charges(1) = qcentermin
             optimum_charges(2) = q1min
             minerror = error
          end if

          ! aeaa: Only continues if above error threshold
          if (minerror > error_thres) then

             ! aeaa: Two point charges added, both along bond axis
             call internal_symmetric_1_2_extra_charges(error, rcentermin, &
                  r1min, r2min, qcentermin, q1min, q2min, remaining_charge, &
                  max_distance, true_V, sphere_positions, no_shells, no_div, &
                  neighbouring_vectors)

             ! aeaa: Added if lower error than previously
             if ( ( ( error + error_reduce) < error_before ) .and. &
                  ( ( error + error_reduce ) < minerror ) ) then
                optimum_posit(:, 1) = rcentermin
                optimum_posit(:, 2) = r1min
                optimum_posit(:, 3) = r2min
                optimum_charges(1) = qcentermin
                optimum_charges(2) = q1min
                optimum_charges(3) = q2min
                minerror = error
                no_charges = 2
             end if
          end if
       end if

       ! aeaa: Two bond case, bisector to bonds used
       if (no_bonded == 2) then

          bond_a = neighbouring_vectors(1,:)
          bond_a = bond_a / internal_norm2(bond_a)
          bond_b = neighbouring_vectors(2,:)
          bond_b = bond_b / internal_norm2(bond_b)

          ! aeaa: Bisector calculated
          dir_vector = - (bond_a/2 + bond_b/2)

          ! aeaa: Normal to plane
          nor_vector = internal_cross(bond_a, bond_b)
          nor_vector = nor_vector / internal_norm2(nor_vector)

          ! aeaa: In plane and normal to bisector
          per_vector = internal_cross(dir_vector, nor_vector)
          per_vector = per_vector / internal_norm2(per_vector)

          call internal_symmetric_2_1_extra_charges(error, rcentermin, r1min,&
               qcentermin, q1min, remaining_charge, max_distance, &
               true_V, sphere_positions, no_shells, no_div, dir_vector)

          if ( ( error + error_reduce ) < error_before ) then
             no_charges = 1
             optimum_posit(:, 1) = rcentermin
             optimum_posit(:, 2) = r1min
             optimum_charges(1) = qcentermin
             optimum_charges(2) = q1min
             minerror = error
          end if

          ! aeaa: Only continue if above error threshold
          if (minerror > error_thres) then

             ! aeaa: This first way puts points perp to the plane of ABC
             call internal_symmetric_2_2_extra_charges(error, rcentermin, &
                  r1min, r2min, qcentermin, q1min, q2min, remaining_charge, &
                  max_distance, true_V, sphere_positions, no_shells, no_div, &
                  dir_vector, nor_vector)

             ! aeaa: Added if lower error than previously
             if ( ( ( error + error_reduce) < error_before ) .and. &
                  ( ( error + error_reduce ) < minerror ) ) then
                optimum_posit(:, 1) = rcentermin
                optimum_posit(:, 2) = r1min
                optimum_posit(:, 3) = r2min
                optimum_charges(1) = qcentermin
                optimum_charges(2) = q1min
                optimum_charges(3) = q2min
                minerror = error
                no_charges = 2
             end if

             ! aeaa: This way puts points in the plane of ABC
             call internal_symmetric_2_2_extra_charges(error, rcentermin, &
                  r1min, r2min, qcentermin, q1min, q2min, remaining_charge, &
                   max_distance, true_V, sphere_positions, no_shells,&
                  no_div, dir_vector, per_vector)

             ! aeaa: Added if lower error than previously
             if ( ( ( error + error_reduce) < error_before )  .and. &
                  ( ( ( no_charges == 2 ) .and. ( error < minerror ) ) .or. &
                  ( ( no_charges == 1 ) .and. ( ( error + error_reduce ) &
                  < minerror ) ) ) ) then
                optimum_posit(:, 1) = rcentermin
                optimum_posit(:, 2) = r1min
                optimum_posit(:, 3) = r2min
                optimum_charges(1) = qcentermin
                optimum_charges(2) = q1min
                optimum_charges(3) = q2min
                minerror = error
                no_charges = 2
             end if

          end if

          if (minerror > error_thres) then

             call internal_symmetric_2_3_extra_charges(error, rcentermin, &
                  r1min, r2min, r3min, qcentermin, q1min, q2min, q3min, &
                  remaining_charge, max_distance, true_V, sphere_positions, &
                  no_shells, no_div, dir_vector, nor_vector)

             ! aeaa: Added if lower error than previously
             if ( ( ( error + error_reduce) < error_before ) .and. &
                  ( ( error + error_reduce ) < minerror ) ) then
                optimum_posit(:, 1) = rcentermin
                optimum_posit(:, 2) = r1min
                optimum_posit(:, 3) = r2min
                optimum_posit(:, 4) = r3min

                optimum_charges(1) = qcentermin
                optimum_charges(2) = q1min
                optimum_charges(3) = q2min
                optimum_charges(4) = q3min
                minerror = error
                no_charges = 3
             end if

             call internal_symmetric_2_3_extra_charges(error, rcentermin, &
                  r1min, r2min, r3min, qcentermin, q1min, q2min, q3min, &
                  remaining_charge, max_distance, true_V, sphere_positions, &
                  no_shells, no_div, dir_vector, per_vector)

             ! aeaa: Added if lower error than previously
             if ( ( ( error + error_reduce) < error_before ) .and. &
                ( ( ( no_charges == 3 ) .and. ( error < minerror ) ) &
                .or. ( ( error + error_reduce ) < minerror ) ) ) then

                optimum_posit(:, 1) = rcentermin
                optimum_posit(:, 2) = r1min
                optimum_posit(:, 3) = r2min
                optimum_posit(:, 4) = r3min

                optimum_charges(1) = qcentermin
                optimum_charges(2) = q1min
                optimum_charges(3) = q2min
                optimum_charges(4) = q3min
                minerror = error
                no_charges = 3
             end if

          end if

       end  if

       ! aeaa: Three bonded, same angle with all three bonds
       if (no_bonded == 3) then

          bond_a = neighbouring_vectors(1,:)
          bond_b = neighbouring_vectors(2,:)
          bond_c = neighbouring_vectors(3,:)

          bond_a = bond_a / internal_norm2(bond_a)
          bond_b = bond_b / internal_norm2(bond_b)
          bond_c = bond_c / internal_norm2(bond_c)

          ! aeaa: Direction vector is an equal angle from all bonds
          dir_vector = internal_cross( (bond_a - bond_b) , (bond_c - bond_b) )
          dir_vector = dir_vector / internal_norm2(dir_vector)

          ! aeaa: Checks correct direction
          if ( dot_product(dir_vector, bond_a) < 1 .and. &
               dot_product(dir_vector, bond_a) > 0 ) then
             dir_vector = - dir_vector
          end if

          sec_dir_vector = - dir_vector

          call internal_symmetric_1_1_extra_charges(error, rcentermin, r1min, &
               qcentermin, q1min, remaining_charge, max_distance, &
               true_V, sphere_positions, no_shells, no_div, dir_vector)

          ! aeaa: Added if lower error than previously
          if ( ( ( error + error_reduce ) < error_before ) ) then
             no_charges = 1
             optimum_posit(:, 1) = rcentermin
             optimum_posit(:, 2) = r1min
             optimum_charges(1) = qcentermin
             optimum_charges(2) = q1min
             minerror = error
          end if

          ! aeaa: H2 special case due to partioning, vector between both H
          if ( all_element_types(neighbour_element_types(1)) == 'H ' .and. &
               all_element_types(neighbour_element_types(2)) == 'H ' ) then
             sec_dir_vector = bond_a + bond_b
             sec_dir_vector = sec_dir_vector / internal_norm2(dir_vector)
          end if

          if ( all_element_types(neighbour_element_types(2)) == 'H ' .and. &
               all_element_types(neighbour_element_types(3)) == 'H ' ) then
             sec_dir_vector = bond_b + bond_c
             sec_dir_vector = sec_dir_vector / internal_norm2(dir_vector)
          end if

          if ( all_element_types(neighbour_element_types(1)) == 'H ' .and. &
               all_element_types(neighbour_element_types(3)) == 'H ' ) then
             sec_dir_vector = bond_a + bond_c
             sec_dir_vector = sec_dir_vector / internal_norm2(dir_vector)
          end if

          ! aeaa: Only continue if above error threshold
          if (minerror > error_thres) then

             call internal_symmetric_3_2_extra_charges(error, rcentermin, &
                  r1min, r2min, qcentermin, q1min, q2min, remaining_charge, &
                  max_distance, true_V, sphere_positions, no_shells, no_div, &
                  dir_vector, sec_dir_vector)

             ! aeaa: Added if lower error than previously
             if ( ( ( error + error_reduce) < error_before ) .and. &
                  ( ( error + error_reduce) < minerror ) ) then
                optimum_posit(:, 1) = rcentermin
                optimum_posit(:, 2) = r1min
                optimum_posit(:, 3) = r2min
                optimum_charges(1) = qcentermin
                optimum_charges(2) = q1min
                optimum_charges(3) = q2min
                minerror = error
                no_charges = 2
             end if

          end  if

       end if

   end if

!   aeaa: Uncomment below if you want to true potential printed to file
!
! if (error_before > error_thres) then
!      open(77, file = element_type(1:2) // number_str // &
!           '_true_V' )
!      write(77, *) grid_spacing_x
!      close(77)
!   end if

   ! aeaa: Deallocate arrays
   deallocate(density, stat=ierr)
   call utils_dealloc_check('ddec_aniso', 'density', ierr)
   deallocate(charges, stat=ierr)
   call utils_dealloc_check('ddec_aniso', 'charges', ierr)
   deallocate(positions, stat=ierr)
   call utils_dealloc_check('ddec_aniso', 'positions', ierr)
   deallocate(neighbouring_vectors, stat=ierr)
   call utils_dealloc_check('ddec_aniso', 'neighbouring_vectors', ierr)


  contains

    subroutine internal_sphere_coordinates(xsphere, ysphere, zsphere, no_div)

      !-----------------------------------------------------------------------!
      !aeaa: Calculates positions of grid points on a sphere                  !
      !    : These points then used as sampling positions for potential       !
      !-----------------------------------------------------------------------!
      ! Output : 'xsphere', x sphere coordinates                              !
      !        : 'ysphere', y sphere coordinates                              !
      !        : 'zsphere', z sphere coordinates                              !
      ! Input  : 'no_div', the number of divisions on the grid used           !
      !-----------------------------------------------------------------------!

      implicit none

      ! aeaa: Input variables
      integer, intent(in) :: no_div ! Number of divisions on the sampling grid

      ! aeaa: Output variables
      real(kind=DP), dimension(no_div, no_div), intent(out) :: xsphere, &
           ysphere, zsphere ! Sphere coordinates

      ! aeaa: Local variables
      real(kind=DP), dimension(no_div, no_div) ::  phi, theta
      real(kind=DP), parameter :: radius = 1.0
      integer :: i, j

      ! aeaa: Creates grid of evenly spaced phi and theta points
      do i = 1,no_div
         do j = 1,no_div
            phi(j,i) = PI * real( i - 1, kind=DP) / real( no_div - 1, kind=DP)
            theta(i,j) = 2.0_DP * PI * real( i - 1 , kind=DP) /  &
                 real( no_div - 1, kind=DP)
         end do
      end do

      ! aeaa: Finds the x,y,z coordinates
      do i = 1,no_div
         do j = 1,no_div
            xsphere(i,j) = radius * sin(phi(i,j)) * cos(theta(i,j))
            ysphere(i,j) = radius * sin(phi(i,j)) * sin(theta(i,j))
            zsphere(i,j) = radius * cos(phi(i,j))
         end do
      end do

    end subroutine internal_sphere_coordinates

  !===========================================================================!

    subroutine internal_symmetric_1_1_extra_charges(errorvar, rcentermin, &
         r1min, qcentermin, q1min, remaining_charge, max_distance, &
         truepotential, sample_positions, no_shells, no_div, dir_vector)

      !-----------------------------------------------------------------------!
      !    : Returns charge and position values that minimizes the error      !
      !    : In the electrostatic potential for the recreation of a           !
      !    : QM potential using one point charge for one bond. The direction  !
      !    : Is along the direction of the bond for one bond.                 !
      !-----------------------------------------------------------------------!
      ! Output : 'q1', the charge calculated for the extra point charge       !
      !        : 'r1', the positions calculated for the extra point charge    !
      !        : 'qcentermin', the charge calculated for the center point     !
      !        : 'rcentermin', the positions of the center (set to the origin)!
      ! Input  : 'max_distance', the maximum distance the charges can be      !
      !        : 'remaining_charge', the charge the extra charges must sum to !
      !        : 'truepotential', the QM potential found at the sampling pos  !
      !        : 'sample_positions',the shells of spherical sampling positions!
      !        : 'no_shells', the number of sampling shells used              !
      !        : 'no_div', the number of grid divisions in sampling pos       !
      !-----------------------------------------------------------------------!

      implicit none

      ! aeaa: Output variables
      real(kind=DP), intent(out) :: errorvar !minimum error found
      real(kind=DP), intent(out) :: rcentermin(3) !positions of center
      real(kind=DP), intent(out) :: r1min(3) !positions of extra charge
      ! aeaa: Charge on center and extra charge
      real(kind=DP), intent(out) :: qcentermin, q1min

      ! aeaa: Input variable
      integer, intent(in) :: no_div  !Number of grid divisions in sampling grid
      integer, intent(in) :: no_shells  !Number of shells for sampling points

      real(kind=DP), intent(in) :: remaining_charge !The charge from DDEC
      real(kind=DP), intent(in) :: max_distance !The max distance from nucleus
      real(kind=DP), intent(in) :: truepotential(no_shells, no_div, &
           no_div)
      ! aeaa: The coordinates of the sample positions
      real(kind=DP), intent(in) :: sample_positions(no_shells, 3, no_div, &
           no_div)

      ! aeaa: Direction Vector of extra points
      real(kind=DP), intent(inout) :: dir_vector(3)

      ! aeaa: Local variables
      integer :: i, j, n
      integer :: no_loops_charge = 200   !Total number of q1 values tried out

      ! aeaa: The charge values (q) and position values (r)
      real(kind=DP) ::  q1, r1(3), charges(2), positions(2,3), qcenter, &
           rcenter(3) = [0.0_DP, 0.0_DP, 0.0_DP]
      real(kind=DP) :: error
      real(kind=DP) :: minerror  !The minimum error must be initialised
      ! aeaa: The minimum charge value the point charge can take
      real(kind=DP) :: minsumcharge = -2.0_DP
      character(10) :: charge_string

      dir_vector = dir_vector / internal_norm2(dir_vector)

      rcenter = [0.0_DP, 0.0_DP, 0.0_DP]
      qcenter = 0.0_DP

      minerror = 10000.0_DP

      do i = 1,no_loops_charge

         q1 = minsumcharge - 2.0 * ( dble(minsumcharge) / &
              dble(no_loops_charge) ) * i
         qcenter = remaining_charge - q1
         r1 = [ 0.0_DP, 0.0_DP, 0.0_DP]

         do while ( internal_norm2(r1) <= max_distance )

            r1 =  r1 + (0.01_DP) * dir_vector

            error = 0
            charges = [ qcenter, q1 ]
            positions(1, :) = rcenter
            positions(2, :) = r1

            ! aeaa: Calculate error for this configuration
            error = internal_calc_error( 2, positions, charges, no_div, &
                 truepotential, sample_positions, no_shells)


            if ( error < minerror )  then
               ! aeaa: New values are accepted if the error is lower

               minerror = error

               rcentermin = rcenter
               r1min = r1

               qcentermin = qcenter
               q1min = q1

            end if

         end do
      end do

      errorvar = minerror

    end subroutine internal_symmetric_1_1_extra_charges

  !===========================================================================!

 subroutine internal_symmetric_1_2_extra_charges(errorvar, rcentermin, r1min, &
         r2min, qcentermin, q1min, q2min, remaining_charge, max_distance, &
         truepotential, sample_positions, no_shells, no_div, dir_vector)

      !-----------------------------------------------------------------------!
      !    : Returns charge and position values that minimizes the error      !
      !    : In the electrostatic potential for the recreation of a           !
      !    : QM potential using two point charges for one bond. The direction !
      !    : Is along the direction of the bond.                              !
      !-----------------------------------------------------------------------!
      ! Output : 'q1', the charge calculated for the extra point charge       !
      !        : 'r1', the positions calculated for the extra point charge    !
      !        : 'qcentermin', the charge calculated for the center point     !
      !        : 'rcentermin', the positions of the center (set to the origin)!
      ! Input  : 'max_distance', the maximum distance the charges can be      !
      !        : 'remaining_charge', the charge the extra charges must sum to !
      !        : 'truepotential', the QM potential found at the sampling pos  !
      !        : 'sample_positions',the shells of spherical sampling positions!
      !        : 'no_shells', the number of sampling shells used              !
      !        : 'no_div', the number of grid divisions in sampling pos       !
      !-----------------------------------------------------------------------!

      implicit none

      ! aeaa: Output variables
      real(kind=DP), intent(out) :: errorvar !minimum error found
      real(kind=DP), intent(out) :: rcentermin(3) !positions of centers
      real(kind=DP), intent(out) :: r1min(3), r2min(3) !positions of charges
      ! aeaa: Charge on center and extra charge
      real(kind=DP), intent(out) :: qcentermin, q1min, q2min

      ! aeaa: Input variable
      integer, intent(in) :: no_div  !Number of grid divisions in sampling grid
      integer, intent(in) :: no_shells  !Number of shells for sampling points

      real(kind=DP), intent(in) :: remaining_charge !The charge from DDEC
      real(kind=DP), intent(in) :: max_distance !The max distance from nucleus
      real(kind=DP), intent(in) :: truepotential(no_shells, no_div, &
           no_div)
      ! aeaa: The coordinates of the sample positions
      real(kind=DP), intent(in) :: sample_positions(no_shells, 3, no_div, &
           no_div)

      ! aeaa: Direction Vector along axis
      real(kind=DP), intent(inout) :: dir_vector(3)

      ! aeaa: Local variables
      integer :: i, j, ii, jj, n
      integer :: no_loops_charge = 40   !Total number of q1 values tried out

      ! aeaa: The charge values (q) and position values (r)
      real(kind=DP) ::  q1, q2, r1(3), charges(3), positions(3,3), qcenter, &
           rcenter(3) = [0.0_DP, 0.0_DP, 0.0_DP], r2(3)
      real(kind=DP) :: error
      real(kind=DP) :: minerror  !The minimum error must be initialised
      ! aeaa: The minimum charge value the point charge can take
      real(kind=DP) :: minsumcharge = -2.0_DP

      ! aeaa: Normalise vector
      dir_vector = dir_vector / internal_norm2(dir_vector)

      rcenter = [0.0_DP, 0.0_DP, 0.0_DP]
      qcenter = 0.0_DP

      minerror = 10000.0_DP

      do i = 1,no_loops_charge

         r1 = [ 0.0_DP, 0.0_DP, 0.0_DP]

         do while ( internal_norm2(r1) <= max_distance )

            r1 =  r1 + (0.05_DP) * dir_vector

            do ii =1,no_loops_charge

               r2 = [ 0.0_DP, 0.0_DP, 0.0_DP]

               do while ( internal_norm2(r2) <= max_distance )

                  r2 =  r2 + (0.05_DP) * dir_vector

                  q1 = minsumcharge - 2.0 * ( dble(minsumcharge) / &
                       dble(no_loops_charge) ) * dble(i)
                  q2 = minsumcharge - 2.0 * ( dble(minsumcharge) / &
                       dble(no_loops_charge) ) * dble(ii)
                  qcenter = remaining_charge - q1 - q2

                  error = 0
                  charges = [ qcenter, q1, q2 ]
                  positions(1, :) = rcenter
                  positions(2, :) = r1
                  positions(3, :) = r2

                  ! aeaa: Calculate error for this configuration
                  error = internal_calc_error( 3, positions, charges, no_div, &
                       truepotential, sample_positions, no_shells)

                  if ( error < minerror ) then
                     ! aeaa: New values are accepted if the error is lower

                     minerror = error

                     rcentermin = rcenter
                     r1min = r1
                     r2min = r2

                     qcentermin = qcenter
                     q1min = q1
                     q2min = q2

                  end if
               end do
            end do
         end do
      end do

      errorvar = minerror

    end subroutine internal_symmetric_1_2_extra_charges

  !===========================================================================!


    subroutine internal_symmetric_2_1_extra_charges(errorvar, rcentermin, &
         r1min, qcentermin, q1min, remaining_charge, max_distance, &
         truepotential, sample_positions, no_shells, no_div, dir_vector)

      !-----------------------------------------------------------------------!
      !    : Returns charge and position values that minimizes the error      !
      !    : In the electrostatic potential for the recreation of a           !
      !    : QM potential using one point charges for two bond. The direction !
      !    : Is along the direction that bisects the two bonds.               !
      !-----------------------------------------------------------------------!
      ! Output : 'q1', the charge calculated for the extra point charge       !
      !        : 'r1', the positions calculated for the extra point charge    !
      !        : 'qcentermin', the charge calculated for the center point     !
      !        : 'rcentermin', the positions of the center (set to the origin)!
      ! Input  : 'max_distance', the maximum distance the charges can be      !
      !        : 'remaining_charge', the charge the extra charges must sum to !
      !        : 'truepotential', the QM potential found at the sampling pos  !
      !        : 'sample_positions',the shells of spherical sampling positions!
      !        : 'no_shells', the number of sampling shells used              !
      !        : 'no_div', the number of grid divisions in sampling pos       !
      !-----------------------------------------------------------------------!

      ! aeaa: Function not the same as 1_1 as it allows negative lambda values

      implicit none

      !Output variables
      real(kind=DP), intent(out) :: errorvar !minimum error found
      real(kind=DP), intent(out) :: rcentermin(3) !positions of center
      real(kind=DP), intent(out) :: r1min(3) !positions of extra charge
      !charge on center and extra charge
      real(kind=DP), intent(out) :: qcentermin, q1min

      !Input variable
      integer, intent(in) :: no_div  !Number of grid divisions in sampling grid
      integer, intent(in) :: no_shells  !Number of shells for sampling points

      real(kind=DP), intent(in) :: remaining_charge !The charge from DDEC
      real(kind=DP), intent(in) :: max_distance !The max distance from nucleus
      real(kind=DP), intent(in) :: truepotential(no_shells, no_div, &
           no_div)
      !The coordinates of the sample positions
      real(kind=DP), intent(in) :: sample_positions(no_shells, 3, no_div, &
           no_div)

      !Direction Vector along axis
      real(kind=DP), intent(inout) :: dir_vector(3)

      ! Local variables
      integer :: i, j, n
      integer :: no_loops_charge = 200   !Total number of q1 values tried out

      !The charge values (q) and position values (r)
      real(kind=DP) ::  q1, r1(3), charges(2), positions(2,3), qcenter, &
           rcenter(3) = [0.0_DP, 0.0_DP, 0.0_DP]
      real(kind=DP) :: error
      real(kind=DP) :: minerror  !The minimum error must be initialised
      !The minimum charge value the point charge can take
      real(kind=DP) :: minsumcharge = -2.0_DP
      character(50) :: charge_string

      dir_vector = dir_vector / internal_norm2(dir_vector)

      rcenter = [0.0_DP, 0.0_DP, 0.0_DP]
      qcenter = 0.0_DP

      minerror = 10000.0_DP

      do i = 1,no_loops_charge

         q1 = minsumcharge - 2.0 * ( minsumcharge / dble(no_loops_charge) ) &
              * dble(i)
         qcenter = remaining_charge - q1

         r1 = - max_distance * dir_vector

         do while ( internal_norm2(r1) <= max_distance )

            r1 =  r1 + (0.01_DP) * dir_vector

            error = 0
            charges = [ qcenter, q1 ]
            positions(1, :) = rcenter
            positions(2, :) = r1

            ! aeaa: Calculate error for this configuration
            error = internal_calc_error( 2, positions, charges, no_div, &
                 truepotential, sample_positions, no_shells)

            if ( error < minerror ) then
               ! aeaa: New values are accepted if the error is lower

               minerror = error

               rcentermin = rcenter
               r1min = r1

               qcentermin = qcenter
               q1min = q1

            end if

         end do
      end do

      errorvar = minerror

    end subroutine internal_symmetric_2_1_extra_charges

  !===========================================================================!

    subroutine internal_symmetric_2_2_extra_charges(errorvar, rcentermin, &
         r1min, r2min, qcentermin, q1min, q2min, remaining_charge, &
         max_distance, truepotential, sample_positions, no_shells, no_div, &
         dir_vector, per_vector)

      !-----------------------------------------------------------------------!
      !    : Returns charge and position values that minimizes the error      !
      !    : In the electrostatic potential for the recreation of a           !
      !    : QM potential using two point charges for two bond. The direction !
      !    : Is along the direction of the sum of a vector that bisects the   !
      !    : The two bonds and either a vector perpendicular to the two bonds !
      !    : Or a vector in the plane and perpendicular to bisector.          !
      !-----------------------------------------------------------------------!
      ! Output : 'q1', the charge calculated for the extra point charge       !
      !        : 'r1', the positions calculated for the extra point charge    !
      !        : 'qcentermin', the charge calculated for the center point     !
      !        : 'rcentermin', the positions of the center (set to the origin)!
      ! Input  : 'max_distance', the maximum distance the charges can be      !
      !        : 'remaining_charge', the charge the extra charges must sum to !
      !        : 'truepotential', the QM potential found at the sampling pos  !
      !        : 'sample_positions',the shells of spherical sampling positions!
      !        : 'no_shells', the number of sampling shells used              !
      !        : 'no_div', the number of grid divisions in sampling pos       !
      !-----------------------------------------------------------------------!

      implicit none

      ! aeaa: Output variables
      real(kind=DP), intent(out) :: errorvar !minimum error found
      real(kind=DP), intent(out) :: rcentermin(3) !positions of center
      real(kind=DP), intent(out) :: r1min(3), r2min(3) !positions of charges
      ! aeaa: Charge on center and extra charge
      real(kind=DP), intent(out) :: qcentermin, q1min, q2min

      ! aeaa: Input variable
      integer, intent(in) :: no_div  !Number of grid divisions in sampling grid
      integer, intent(in) :: no_shells  !Number of shells for sampling points

      real(kind=DP), intent(in) :: remaining_charge !The charge from DDEC
      real(kind=DP), intent(in) :: max_distance !The max distance from nucleus
      real(kind=DP), intent(in) :: truepotential(no_shells, no_div, &
           no_div)
      ! aeaa: The coordinates of the sample positions
      real(kind=DP), intent(in) :: sample_positions(no_shells, 3, no_div, &
           no_div)

      ! aeaa: Direction Vector along axis
      real(kind=DP), intent(inout) :: dir_vector(3), per_vector(3)

      ! aeaa: Local variables
      integer :: i, j, ii, jj, n
      integer :: no_loops_charge = 40   !Total number of q1 values tried out

      ! aeaa: The charge values (q) and position values (r)
      real(kind=DP) ::  q1, q2, r1(3), r3(3), charges(3), positions(3,3), &
           qcenter, rcenter(3) = [0.0_DP, 0.0_DP, 0.0_DP], r2(3)
      real(kind=DP) :: perp_dis(3), para_dis(3)
      real(kind=DP) :: error
      real(kind=DP) :: minerror  !The minimum error must be initialised
      ! aeaa: The minimum charge value the point charge can take
      real(kind=DP) :: minsumcharge = -2.0_DP

      ! aeaa: Normalise vectors
      dir_vector = dir_vector / internal_norm2(dir_vector)

      rcenter = [0.0_DP, 0.0_DP, 0.0_DP]
      qcenter = 0.0_DP

      minerror = 10000.0_DP

      do i = 1,no_loops_charge

         perp_dis = [ 0.0_DP, 0.0_DP, 0.0_DP]

         do while ( internal_norm2(perp_dis) <= max_distance )

            perp_dis = perp_dis + (0.05_DP) * per_vector
            para_dis = [ 0.0_DP, 0.0_DP, 0.0_DP]

            do while ( internal_norm2(para_dis) <= max_distance )

               para_dis =  para_dis + (0.05_DP) * dir_vector

               q1 = minsumcharge - 2.0_DP * ( dble(minsumcharge) /  &
                    dble(no_loops_charge) ) * dble(i)

               ! aeaa: Want symettry in charges allocated
               q2 = q1
               qcenter = remaining_charge - q1 - q2

               r1 = para_dis + perp_dis
               r2 = para_dis - perp_dis

               error = 0
               charges = [ qcenter, q1, q2 ]
               positions(1, :) = rcenter
               positions(2, :) = r1
               positions(3, :) = r2

               ! aeaa: Calculate error for this configuration
               error = internal_calc_error( 3, positions, charges, no_div, &
                    truepotential, sample_positions, no_shells)

               if ( ( error < minerror ) .and. &
                    ( internal_norm2(r1) <= max_distance ))  then
                  ! aeaa: New values are accepted if the error is lower

                  minerror = error

                  rcentermin = rcenter
                  r1min = r1
                  r2min = r2

                  qcentermin = qcenter
                  q1min = q1
                  q2min = q2

               end if
            end do
         end do
      end do

      errorvar = minerror

    end subroutine internal_symmetric_2_2_extra_charges


    !===========================================================================!


    subroutine internal_symmetric_2_3_extra_charges(errorvar, rcentermin, &
         r1min, r2min, r3min, qcentermin, q1min, q2min, q3min, &
         remaining_charge, max_distance, truepotential, sample_positions, &
         no_shells, no_div, dir_vector, per_vector)

      !-----------------------------------------------------------------------!
      !    : Returns charge and position values that minimizes the error      !
      !    : In the electrostatic potential for the recreation of a           !
      !    : QM potential using 3 point charges for two bonds. The direction  !
      !    : Is along the direction of the sum of a vector that bisects the   !
      !    : The two bonds and either a vector perpendicular to the two bonds !
      !    : Or a vector in the plane and perpendicular to bisector.          !
      !    : The second vector is the negative of the first and the third is  !
      !    : Along the bisector.                                              !
      !-----------------------------------------------------------------------!
      ! Output : 'q1', the charge calculated for the extra point charge       !
      !        : 'r1', the positions calculated for the extra point charge    !
      !        : 'qcentermin', the charge calculated for the center point     !
      !        : 'rcentermin', the positions of the center (set to the origin)!
      ! Input  : 'max_distance', the maximum distance the charges can be      !
      !        : 'remaining_charge', the charge the extra charges must sum to !
      !        : 'truepotential', the QM potential found at the sampling pos  !
      !        : 'sample_positions',the shells of spherical sampling positions!
      !        : 'no_shells', the number of sampling shells used              !
      !        : 'no_div', the number of grid divisions in sampling pos       !
      !-----------------------------------------------------------------------!

      implicit none

      ! aeaa: Output variables
      real(kind=DP), intent(out) :: errorvar ! Minimum error found
      real(kind=DP), intent(out) :: rcentermin(3) ! Positions of center
      real(kind=DP), intent(out) :: r1min(3), r2min(3), r3min(3) ! Positions
      ! aeaa: charge on center and extra charge
      real(kind=DP), intent(out) :: qcentermin, q1min, q2min, q3min

      ! aeaa: Input variable
      integer, intent(in) :: no_div  !Number of grid divisions in sampling grid
      integer, intent(in) :: no_shells  !Number of shells for sampling points

      real(kind=DP), intent(in) :: remaining_charge !The charge from DDEC
      real(kind=DP), intent(in) :: max_distance !The max distance from nucleus
      real(kind=DP), intent(in) :: truepotential(no_shells, no_div, &
           no_div)
      ! aeaa: The coordinates of the sample positions
      real(kind=DP), intent(in) :: sample_positions(no_shells, 3, no_div, &
           no_div)
      ! aeaa: Direction Vector along axis and perpendicular
      real(kind=DP), intent(inout) :: dir_vector(3), per_vector(3)

      ! aeaa: Local variables
      integer :: i, j, ii, jj, n
      integer :: no_loops_charge = 40   !Total number of q1 values tried out

      ! aeaa: The charge values (q) and position values (r)
      real(kind=DP) ::  q1, q2, q3, r1(3), r2(3), r3(3),charges(4), &
           positions(4,3), qcenter, rcenter(3) = [0.0_DP, 0.0_DP, 0.0_DP]
      real(kind=DP) :: perp_dis(3), para_dis(3)
      real(kind=DP) :: error
      real(kind=DP) :: minerror  !The minimum error must be initialised
      ! aeaa: The minimum charge value the point charge can take
      real(kind=DP) :: minsumcharge = -2.0_DP
      real(kind=DP) :: a, b

      ! aeaa: Normalise vector
      dir_vector = dir_vector / internal_norm2(dir_vector)

      rcenter = [0.0_DP, 0.0_DP, 0.0_DP]
      qcenter = 0.0_DP

      minerror = 10000.0_DP

      do i = 1,no_loops_charge
         do ii = 1,no_loops_charge

            perp_dis = [ 0.0_DP, 0.0_DP, 0.0_DP]

            do while ( internal_norm2(perp_dis) <= max_distance )

               perp_dis = perp_dis + (0.05_DP) * per_vector
               para_dis = [ 0.0_DP, 0.0_DP, 0.0_DP]

               do while ( internal_norm2(para_dis) <= max_distance )

                  para_dis =  para_dis + (0.05_DP) * dir_vector

                  q1 = minsumcharge - 2.0_DP * ( dble(minsumcharge) &
                       / dble(no_loops_charge) ) * dble(i)
                  q2 = q1
                  q3 = minsumcharge - 2.0_DP * ( dble(minsumcharge) &
                       / dble(no_loops_charge) ) * dble(ii)
                  qcenter = remaining_charge - q1 - q2 - q3

                  r1 = para_dis + perp_dis
                  r2 = para_dis - perp_dis

                  a = internal_norm2(r1)
                  b = internal_norm2(r2)

                  ! aeaa: Third charge is restricted to same distance as
                  !       r2 and r3 from the center
                  r3 = (dir_vector * sqrt( (a*a) + (b*b) ))

                  error = 0
                  charges = [ qcenter, q1, q2, q3 ]
                  positions(1, :) = rcenter
                  positions(2, :) = r1
                  positions(3, :) = r2
                  positions(4, :) = r3

                  ! aeaa: Calculate error for this configuration
                  error = internal_calc_error( 4, positions, charges, no_div, &
                       truepotential, sample_positions, no_shells)

                  if ( ( error < minerror ) .and. &
                       ( internal_norm2(r1) <= max_distance ) &
                       .and. ( q1 < 2.0_DP ) )  then
                     ! aeaa: New values are accepted if the error is lower

                     minerror = error

                     rcentermin = rcenter
                     r1min = r1
                     r2min = r2
                     r3min = r3

                     qcentermin = qcenter
                     q1min = q1
                     q2min = q2
                     q3min = q3

                  end if

                  ! aeaa: With r3 on opposite side
                  error = 0
                  charges = [ qcenter, q1, q2, q3 ]
                  positions(1, :) = rcenter
                  positions(2, :) = r1
                  positions(3, :) = r2
                  positions(4, :) = -r3

                  ! aeaa: Calculate error for this configuration
                  error = internal_calc_error( 4, positions, charges, no_div, &
                       truepotential, sample_positions, no_shells)

                  if ( ( error < minerror ) .and. &
                       ( internal_norm2(r1) <= max_distance ) &
                       .and. ( q1 < 2.0_DP ) )  then
                     ! aeaa: New values are accepted if the error is lower

                     minerror = error

                     rcentermin = rcenter
                     r1min = r1
                     r2min = r2
                     r3min = -r3

                     qcentermin = qcenter
                     q1min = q1
                     q2min = q2
                     q3min = q3

                  end if

               end do
            end do
         end do
      end do
         errorvar = minerror

    end subroutine internal_symmetric_2_3_extra_charges

  !===========================================================================!

    subroutine internal_symmetric_3_2_extra_charges(errorvar, rcentermin, &
         r1min, r2min, qcentermin, q1min, q2min, remaining_charge, &
         max_distance, truepotential, sample_positions, no_shells, no_div, &
         dir_vector, sec_dir_vector)

      !-----------------------------------------------------------------------!
      !    : Returns charge and position values that minimizes the error      !
      !    : In the electrostatic potential for the recreation of a           !
      !    : QM potential using 2 point charges for 3 bonds. The direction    !
      !    : Is along the direction that makes the same angle with all three  !
      !    : Bonds. The second direction is in the negative directions of the !
      !    : First direction except for H2 groups when the second directions  !
      !    : Is in bisector of the bond to the H atoms.                       !
      !-----------------------------------------------------------------------!
      ! Output : 'q1', the charge calculated for the extra point charge       !
      !        : 'r1', the positions calculated for the extra point charge    !
      !        : 'qcentermin', the charge calculated for the center point     !
      !        : 'rcentermin', the positions of the center (set to the origin)!
      ! Input  : 'max_distance', the maximum distance the charges can be      !
      !        : 'remaining_charge', the charge the extra charges must sum to !
      !        : 'truepotential', the QM potential found at the sampling pos  !
      !        : 'sample_positions',the shells of spherical sampling positions!
      !        : 'no_shells', the number of sampling shells used              !
      !        : 'no_div', the number of grid divisions in sampling pos       !
      !-----------------------------------------------------------------------!

      implicit none

      ! aeaa: Output variables
      real(kind=DP), intent(out) :: errorvar !minimum error found
      real(kind=DP), intent(out) :: rcentermin(3) !positions of center
      real(kind=DP), intent(out) :: r1min(3), r2min(3) !positions of charges
      ! aeaa: charge on center and extra charge
      real(kind=DP), intent(out) :: qcentermin, q1min, q2min

      ! aeaa: Input variable
      integer, intent(in) :: no_div  !Number of grid divisions in sampling grid
      integer, intent(in) :: no_shells  !Number of shells for sampling points

      real(kind=DP), intent(in) :: remaining_charge !The charge from DDEC
      real(kind=DP), intent(in) :: max_distance !The max distance from nucleus
      real(kind=DP), intent(in) :: truepotential(no_shells, no_div, &
           no_div)
      ! aeaa: The coordinates of the sample positions
      real(kind=DP), intent(in) :: sample_positions(no_shells, 3, no_div, &
           no_div)
     real(kind=DP), intent(inout) :: dir_vector(3), sec_dir_vector(3)

      ! aeaa: Local variables
      integer :: i, j, ii, jj, n
      integer :: no_loops_charge = 40  !Total number of q1 values tried out

      ! aeaa: The charge values (q) and position values (r)
      real(kind=DP) ::  q1, q2, r1(3), charges(3), positions(3,3), qcenter, &
           rcenter(3) = [0.0_DP, 0.0_DP, 0.0_DP], r2(3)
      real(kind=DP) :: error
      real(kind=DP) :: minerror  !The minimum error must be initialised
      ! aeaa: The minimum charge value the point charge can take
      real(kind=DP) :: minsumcharge = -2.0_DP

      ! aeaa: Normalise vector
      dir_vector = dir_vector / internal_norm2(dir_vector)
      sec_dir_vector = sec_dir_vector / internal_norm2(sec_dir_vector)

      rcenter = [0.0_DP, 0.0_DP, 0.0_DP]
      qcenter = 0.0_DP

      minerror = 10000.0_DP

      do i = 1,no_loops_charge
         do ii =1,no_loops_charge
            r1 = [ 0.0_DP, 0.0_DP, 0.0_DP]

            do while ( internal_norm2(r1) <= max_distance )

               r1 =  r1 + (0.05_DP) * dir_vector
               ! aeaa: vector can be positive or negative
               r2 = - max_distance * sec_dir_vector

               do while ( internal_norm2(r2) <= max_distance )

                  q1 = minsumcharge - 2.0 * ( dble(minsumcharge) / &
                       dble(no_loops_charge) ) * dble(i)
                  q2 = minsumcharge - 2.0 * ( dble(minsumcharge) / &
                       dble(no_loops_charge) ) * dble(ii)
                  qcenter = remaining_charge - q1 - q2

                  r2 = r2 + (0.05_DP) * sec_dir_vector

                  error = 0
                  charges = [ qcenter, q1, q2 ]
                  positions(1, :) = rcenter
                  positions(2, :) = r1
                  positions(3, :) = r2

                  ! aeaa: Calculate error for this configuration
                  error = internal_calc_error( 3, positions, charges, no_div, &
                       truepotential, sample_positions, no_shells)

                  if ( ( error < minerror ) .and. &
                       ( internal_norm2(r1) < max_distance ) &
                       .and. ( q1 < 2.0_DP ) )  then
                     ! aeaa: New values are accepted if the error is lower
                     minerror = error

                     rcentermin = rcenter
                     r1min = r1
                     r2min = r2

                     qcentermin = qcenter
                     q1min = q1
                     q2min = q2

                  end if
               end do
            end do
         end do
      end do

      errorvar = minerror

    end subroutine internal_symmetric_3_2_extra_charges

  !===========================================================================!

    function internal_calc_error(no_charges, positions, charges, no_div, &
         truepotential, sample_positions, no_shells)
      !-----------------------------------------------------------------------!
      !aeaa: Calculates the error between the actual potential and the        !
      !    : potential from the point charges.                                !
      !-----------------------------------------------------------------------!
      ! Output : 'V', the potential found                                     !
      ! Input  : 'sample_positions', the position the potential is found at   !
      !        : 'no_charges', the number of point charges used               !
      !        : 'positions', the point charge positions                      !
      !        : 'charges', the point charge charge values                    !
      !-----------------------------------------------------------------------!

      implicit none

      ! aeaa: Output variables
      real(kind=DP) :: internal_calc_error !error value

      ! aeaa: Input variables
      integer, intent(in) :: no_charges, no_div, no_shells

      ! aeaa: Charge of point charges
      real(kind=DP), intent(in) :: charges(no_charges)
      ! aeaa: Positions of point charges
      real(kind=DP), intent(in) :: positions(no_charges, 3)
      ! aeaa: The true potential at the sample positions
      real(kind=DP), intent(in) :: truepotential(no_shells, no_div, no_div)
      ! aeaa: The positions the potential is found
      real(kind=DP), intent(in) :: sample_positions(no_shells,3,no_div, no_div)

      ! aeaa: Local variables
      integer :: nn, i,j ! Indices
      real(kind=DP) :: point_charge_V !Potential due to point charges

      internal_calc_error = 0.0_DP
      point_charge_V = 0.0_DP

      do nn = 1, no_shells
         do i = 1,no_div
            do j = 1,no_div
               ! aeaa: Point charge potential
               point_charge_V = internal_V( sample_positions(nn, :, i, j), &
                    no_charges, positions, charges )
               ! aeaa: Difference calculated
               internal_calc_error = internal_calc_error + &
                    abs( point_charge_V - truepotential(nn, i, j) )
            end do
         end do
      end do

      internal_calc_error  = ( internal_calc_error / &
           (no_div * no_div * no_shells) )

    end function internal_calc_error

  !===========================================================================!

    function internal_V(sample_positions, no_charges, positions, charges)
      !-----------------------------------------------------------------------!
      !    : Returns the value of the potential at sample postion due to the  !
      !    : Positions and charges that are input                             !
      !-----------------------------------------------------------------------!
      ! Output : 'V', the potential found                                     !
      ! Input  : 'sample_positions', the position the potential is found at   !
      !        : 'no_charges', the number of point charges used               !
      !        : 'positions', the point charge positions                      !
      !        : 'charges', the point charge charge values                    !
      !-----------------------------------------------------------------------!

      implicit none

      ! aeaa: Output variables
      real(kind=DP) :: internal_V

      ! aeaa: Input variables
      integer,  intent(in) :: no_charges

      ! aeaa: The posisitons the potential is sampled
      real(kind=DP), intent(in) :: sample_positions(3)
      ! aeaa: Charge of point charges
      real(kind=DP), intent(in) :: charges(no_charges)
      ! aeaa: Positions of point charges
      real(kind=DP), intent(in) :: positions(no_charges, 3)

      ! aeaa; Local variables
      integer :: ii

      internal_V = 0.0_DP

      do ii = 1,no_charges
         ! If statement as otherwise causes very large values
         ! This should never happen unless the max distance is altered
         if (  internal_norm2(positions(ii,:) - sample_positions) > 0.01_DP ) &
              then
            internal_V = internal_V + (charges(ii)) / &
                 internal_norm2(positions(ii,:) - sample_positions)
         end if
      end do

     internal_V = internal_V * 334.9_DP !change units to kcal/mol

    end function internal_V

  !===========================================================================!

    function internal_norm2(vector)

      implicit none

      real(kind=DP) :: internal_norm2
      real(kind=DP), intent(in) :: vector(3)

      internal_norm2 =  sqrt(sum(vector * vector))

    end function internal_norm2

  !===========================================================================!

    function internal_cross(vector_a, vector_b)

      implicit none

      real(kind=DP) :: internal_cross(3)
      real(kind=DP), intent(in) :: vector_a(3), vector_b(3)

      internal_cross = &
           [((vector_a(2)*vector_b(3)) -  (vector_a(3)*vector_b(2))), &
           -((vector_a(1)*vector_b(3)) - (vector_a(3)*vector_b(1))), &
           ((vector_a(1)*vector_b(2)) - (vector_a(2)*vector_b(1)))]

    end function internal_cross

  !===========================================================================!

  end subroutine ddec_off_site_charges
  !===========================================================================!

end module ddec
