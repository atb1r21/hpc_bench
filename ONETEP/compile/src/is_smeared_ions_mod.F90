! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                         I S_ S M E A R E D _ I O N S                        !
!=============================================================================!
! This module adds the smeared ion representation due to                      !
! [1] Scherlis et al, described in J. Chem. Phys. 124 (2006).                 !
!-----------------------------------------------------------------------------!
! Written by Jacek Dziedzic in May 2010, based largely on the abovementioned  !
! paper (Appendix A) and a similar module for CASTEP due to Hatem H Helal.    !
!-----------------------------------------------------------------------------!

module is_smeared_ions

  use constants, only: DP
  use utils, only: utils_trace_in, utils_trace_out

  implicit none
  private

  public :: smeared_ion_apply_Vloc_corr
  public :: smeared_ion_calc_E_self
  public :: smeared_ion_calc_E_smeared_OBC
  public :: smeared_ion_calc_E_smeared_PBC
  public :: smeared_ion_exit
  public :: smeared_ion_hartree
  public :: smeared_ion_density
  public :: smeared_ion_ii_forces
  public :: smeared_ion_force_correction
  public :: smeared_ion_initialise

  real(kind=DP), public, save :: smeared_ion_E_self
  real(kind=DP), public, save :: smeared_ion_E_smeared

  ! jd: Smeared-ion density
  real(kind=DP), public, allocatable, save :: rho_ion(:,:,:)
  ! JCW: Smeared ion electrostatic potential
  ! JCW: (Currently only used for fully periodic BCs)
  real(kind=DP), public, allocatable, save :: phi_ion(:,:,:)
  ! kkbd: Fully or partially masks the dielectric from exclusion regions defined
  !       by the is_dielectric_exclusions block. Also fully excludes the
  !       dielectric for gridpoints inside core regions, where we want to ensure
  !       eps is 1.0 and its derivative is 0.0 to prevent spurious cavities from
  !       forming.
  real(kind=DP), public, allocatable, save :: dielectric_mask(:,:,:)

  ! jd: Regions manually excluded from dielectric penetration (eps==1).
  !     idx1 -- region #
  !     idx2 -- 1: region type (1 = sphere, 2 = box, 3 = xcyl, 4 = ycyl,
  !                             5 = zcyl, 0 = last entry)
  !             if sphere, 2-4: sphere centre, 5: sphere radius, 6-7: unused
  !             if box, 2-4: lo-corner, 5-7: hi-corner
  !             if xcyl, 2-4: y-coord, z-coord, radius
  !             if ycyl, 2-4: x-coord, z-coord, radius
  !             if zcyl, 2-4: x-coord, y-coord, radius
  integer, public, parameter       :: n_max_exclusions = 10000
  real(kind=DP), public, save      :: dielectric_exclusions(n_max_exclusions,7)


contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_initialise(mdl,is_periodic)
    !=========================================================================!
    ! Initializes the smeared ions:                                           !
    ! - generates the smeared ion density,                                    !
    ! - calculates the self-interaction term (constant),                      !
    ! - calculates the nonself-interaction term (constant until the ions move)!
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   elements    (input): The array describing the ions.                   !
    !   is_periodic (input): Logical array indicating whether the grid is     !
    !                        periodic (T) or not (F) along each lattice       !
    !                        vector. This is used to determine how the        !
    !                        smeared ion-smeared ion interaction is computed. !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    ! Modified by James C. Womack to support periodic BCs, 05/2017            !
    ! Modified for embedding by Joseph Prentice, 09/2020                      !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use model_type, only: MODEL
    use utils, only: utils_abort, utils_alloc_check

    implicit none

    ! jd: Arguments
    type(MODEL), intent(in) :: mdl
    logical, intent(in)     :: is_periodic(3)

    ! Parameters
    character(len=*), parameter :: myself = "smeared_ion_initialise"

    ! jd: Internal variables
    integer :: ierr ! jd: Error flag
    integer :: isub ! jcap: region counter

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_initialise')

    ! jd: Allocate the array that will hold the smeared ion density
    allocate(rho_ion(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('smeared_ion_initialise','rho_ion',ierr)

    ! jd: Allocate the array that will hold the dielectric mask
    allocate(dielectric_mask(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('smeared_ion_initialise','dielectric_mask',ierr)

    if (all(is_periodic)) then
       ! JCW: If fully periodic BCs, then allocate the module global variable
       ! JCW: to store the potential due to the smeared ion density
       ! JCW: --> This will be used to calculate the smeared ion-smeared ion
       ! JCW:     interaction in PBCs and the correction to the local pseudopotential
       allocate(phi_ion(mdl%fine_grid%ld1, mdl%fine_grid%ld2, &
         mdl%fine_grid%max_slabs12),stat=ierr)
       call utils_alloc_check('smeared_ion_initialise','phi_ion',ierr)
    end if

    ! jd: Populate this array (and dielectric_mask)
    call smeared_ion_generate_density(mdl%elements,mdl%fine_grid,mdl%cell,is_periodic,mdl%nat)

    ! jd: Compute constant terms to energy
    ! jcap: loop over regions and sum up
    smeared_ion_E_self = 0.0_DP
    do isub=1, mdl%nsub
       smeared_ion_E_self = smeared_ion_E_self + &
            &smeared_ion_calc_E_self(mdl%regions(isub)%elements, mdl%regions(isub)%par)
    end do
     if(pub_on_root) then
       write(stdout,'(a,f0.6)') 'Smeared ion self correction:     ', smeared_ion_E_self
    end if
    ! jd: NB, this must be invalidated when ionic positions change!

    if (.not.any(is_periodic)) then
       ! JCW: Fully open BCs

       ! JCW: Compute the non-self-interaction energy of smeared ions using
       ! JCW: integrals over pairs of smeared ions, excluding self-interaction
       smeared_ion_E_smeared = smeared_ion_calc_E_smeared_OBC(mdl)

    else if (all(is_periodic)) then
       ! JCW: Fully periodic BCs

       ! JCW: Compute the potential due to the smeared ion density and store as
       ! JCW: a module global variable so it can be re-used to compute the
       ! JCW: correction to the local pseudopotential
       call smeared_ion_generate_potential_PBC(phi_ion,rho_ion,mdl%fine_grid)
       ! JCW: NB, as with rho_ion, this must be invalidated when ionic positions
       ! JCW: change

       ! JCW: Compute the non-self-interaction energy as the interaction of the
       ! JCW: full smeared ion density with the potential due to this density
       ! JCW: and then subtract the self-interaction energy
       smeared_ion_E_smeared = smeared_ion_calc_E_smeared_PBC(phi_ion,rho_ion,&
            smeared_ion_E_self,mdl%fine_grid)
    else
       ! JCW: Mixed periodic/open BCs
       ! JCW: Not currently supported
       call utils_abort("Error in "//myself//": Mixed periodic/open BCs are not &
            &currently supported for smeared ions.")
    end if

    call utils_trace_out('smeared_ion_initialise')

  end subroutine smeared_ion_initialise
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_exit
    !=========================================================================!
    ! Cleans up after the smeared ions:                                       !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! jd: Internal variables
    integer :: ierr ! jd: Error flag

    !------------------------------------------------------------------------

    ! jd: Deallocate the array that held the dielectric mask
    deallocate(dielectric_mask,stat=ierr)
    call utils_dealloc_check('smeared_ion_exit','dielectric_mask',ierr)

    ! jd: Deallocate the array that held the smeared ion density
    deallocate(rho_ion,stat=ierr)
    call utils_dealloc_check('smeared_ion_exit','rho_ion',ierr)

    if (allocated(phi_ion)) then
       ! JCW: Deallocate the array that held the smeared ion potential
       ! JCW: (Not always allocated)
       deallocate(phi_ion,stat=ierr)
       call utils_dealloc_check('smeared_ion_exit','phi_ion',ierr)
    end if

  end subroutine smeared_ion_exit
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function smeared_ion_calc_E_self(elements,par)
    !=========================================================================!
    ! Calculates the self-interaction energy of smeared ions.                 !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   elements(input): The array describing the ions.                       !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    ! Modified to remove pub_par by Joseph Prentice, September 2018           !
    !=========================================================================!

    use comms, only: comms_reduce, pub_on_root, pub_my_proc_id
    use constants, only: DP, pi, stdout
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_is_smeared_ion_width

    implicit none

    ! jd: Arguments
    type(PARAL_INFO), intent(in) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)

    ! jd: Internal variables
    integer :: my_first_at, my_last_at, I      ! jd: Indices and bounds
    real(kind=DP) :: inv_sigma_I               ! jd: Inverse of Gaussian spread
    real(kind=DP) :: Z_I                       ! jd: Core charge
    real(kind=DP) :: E_self                    ! jd: Accumulated quantity

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_calc_E_self')

    inv_sigma_I = 1.0_DP / pub_is_smeared_ion_width

    my_first_at = par%first_atom_on_proc(pub_my_proc_id)
    my_last_at = my_first_at + par%num_atoms_on_proc(pub_my_proc_id) - 1

    E_self = 0.0_DP

    ! jd: Loop over my atoms
    do I = my_first_at, my_last_at
       Z_I = real(elements(I)%ion_charge, kind=DP)
       E_self = E_self + Z_I*Z_I * inv_sigma_I
    end do

    ! jd: Apply prefactor and reduce
    E_self = E_self / (-sqrt(2.0_DP * pi))
    call comms_reduce('SUM', E_self)

    call utils_trace_out('smeared_ion_calc_E_self')

    smeared_ion_calc_E_self = E_self

  end function smeared_ion_calc_E_self
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function smeared_ion_calc_E_smeared_OBC(mdl)
    !=========================================================================!
    ! Calculates the non-self-interaction energy of smeared ions in fully     !
    ! open BCs.                                                               !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   mdl(input): The model for the system.                                 !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    ! Modified to remove pub_par by Joseph Prentice, September 2018           !
    ! Modified for embedding by Joseph Prentice, September 2020               !
    !=========================================================================!

    use comms, only: comms_reduce, pub_on_root, pub_my_proc_id
    use constants, only: DP, stdout
    use geometry, only: magnitude, OPERATOR(-)
    use model_type, only: MODEL
    use rundat, only: pub_is_smeared_ion_width
!$  use rundat, only: pub_threads_max
    use utils, only: utils_erf

    implicit none

    ! jd: Arguments
    type(MODEL), intent(in) :: mdl

    ! jd: Internal variables
    integer :: my_first_at, my_last_at, I, J         ! jd: Indices and bounds
    integer :: isub, global_i                        ! jcap: Region and global atom no. counters
    real(kind=DP) :: sigma_I, sigma_J, sigma_IJ      ! jd: Gaussian spreads
    real(kind=DP) :: Z_I, Z_J                        ! jd: Charges
    real(kind=DP) :: R_IJ                            ! jd: Distance betw. cores
    real(kind=DP) :: E_smeared                       ! jd: Accumulated quantity

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_calc_E_smeared_OBC')

    E_smeared = 0.0_DP

    ! jcap: loop over regions
    do isub = 1, mdl%nsub

       my_first_at = mdl%regions(isub)%par%first_atom_on_proc(pub_my_proc_id)
       my_last_at = my_first_at + mdl%regions(isub)%par%num_atoms_on_proc(pub_my_proc_id) - 1

       ! jd: Loop over my atoms
!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(I,global_i,Z_I,sigma_I,J,Z_J,R_IJ,sigma_J,sigma_IJ) &
!$OMP SHARED(mdl,my_first_at,my_last_at,isub,pub_is_smeared_ion_width, &
!$OMP      pub_threads_max) &
!$OMP REDUCTION(+:E_smeared)
       do I = my_first_at, my_last_at

          ! jcap: find the global atom number
          global_i = mdl%regions(isub)%elements(I)%global_atom_number

          Z_I = real(mdl%elements(global_i)%ion_charge, kind=DP)
          sigma_I = pub_is_smeared_ion_width

          ! jd: Loop over all other atoms
!$OMP DO
          do J = 1, mdl%nat

             if (J == global_i) cycle

             ! jd: Find charge and distance of other atom
             Z_J = real(mdl%elements(J)%ion_charge,kind=DP)
             R_IJ = magnitude(mdl%elements(global_i)%centre - mdl%elements(J)%centre)

             ! jd: Calculate sqrt(sigma_i^2 * sigma_j^2)
             sigma_J = pub_is_smeared_ion_width
             sigma_IJ = sqrt(sigma_I*sigma_I + sigma_J*sigma_J)

             ! jd: Add Gaussian smeared charge contribution to energy
             E_smeared = E_smeared + Z_i*Z_J / R_IJ * utils_erf(R_IJ / sigma_IJ)

          end do
!$OMP END DO

       end do
!$OMP END PARALLEL

    end do

    ! jd: E_smeared has a minus at the front and the sum runs over I<J.
    !     Since I/=J was excluded, we just have to divide by two...
    E_smeared = E_smeared * (-0.5_DP)
    call comms_reduce('SUM', E_smeared)

    if(pub_on_root) then
       write(stdout,'(a,f0.6)') 'Smeared ion non-self correction: ', E_smeared
    end if

    smeared_ion_calc_E_smeared_OBC = E_smeared

    call utils_trace_out('smeared_ion_calc_E_smeared_OBC')

  end function smeared_ion_calc_E_smeared_OBC
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  function smeared_ion_calc_E_smeared_PBC(phi_si,rho_si,E_self,grid) result(E_smeared)
    !=========================================================================!
    ! Calculates the non-self-interaction energy of smeared ions in fully     !
    ! periodic BCs, using the previously computed periodic electrostatic      !
    ! potential due to the smeared ions, and subtracting off the previously   !
    ! computed self-interaction energy.                                       !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   phi_si (input)  : The potential due to the smeared ions               !
    !   rho_si (input)  : The charge density of the smeared ions              !
    !   E_self (input)  : Self-interaction energy for smeared ions (i.e.      !
    !                     interaction of Gaussians with themselves but NOT    !
    !                     their periodic images                               !
    !   grid   (input)  : GRID_INFO instance describing the grids used to     !
    !                     store the charge and potential and also for the     !
    !                     grid used for integration                           !
    ! Result:                                                                 !
    !   E_smeared       : The non-self-interaction energy of the smeared ions !
    !                     in fully periodic BCs (i.e. interaction of          !
    !                     Gaussians with each other and their periodic images !
    !                     but not themselves)                                 !
    ! ----------------------------------------------------------------------- !
    ! Written by James C. Womack, 05/2017                                     !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: stdout
    use integrals, only: integrals_product_on_grid

    implicit none

    ! Result
    real(kind=DP) :: E_smeared

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in)   :: phi_si(grid%ld1,grid%ld2,grid%max_slabs12)
    real(kind=DP), intent(in)   :: rho_si(grid%ld1,grid%ld2,grid%max_slabs12)
    real(kind=DP), intent(in)   :: E_self

    ! Parameters
    character(len=*), parameter :: myself = "smeared_ion_calc_E_smeared_PBC"
    real(kind=DP), parameter    :: fac = -0.5_dp ! Avoid double counting of
                                                 ! ion-ion interactions, negative
                                                 ! as this is a correction

    call utils_trace_in(myself)

    ! JCW: Note that integrals_product_on_grid gives the integral over the
    ! JCW: entire grid (not just the region on the local MPI rank), so no
    ! JCW: additional reduction of E_smeared is required
    E_smeared = &
         fac * integrals_product_on_grid(grid, phi_si, rho_si) - E_self

    if(pub_on_root) then
       write(stdout,'(a,f0.6)') 'Smeared ion non-self correction: ', E_smeared
    end if

    call utils_trace_out(myself)

  end function smeared_ion_calc_E_smeared_PBC

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_apply_vloc_corr(elements, v_loc_on_fine, grid, is_periodic, par)
    !=========================================================================!
    ! Applies the smeared-ion correction to vloc. How this is calculated      !
    ! depends on the boundary conditions used.                                !
    !                                                                         !
    ! Fully open BCs:                                                         !
    ! - The correction is computed by summing over per-ion potentials.        !
    !                                                                         !
    ! Fully periodic BCs:                                                     !
    ! - The previously computed periodic potential due to the smeared ions    !
    !   is used to correct the local pseudopotential.                         !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   elements(input):      The array describing the ions.                  !
    !   v_loc_on_fine(inout): Vloc on the fine grid.                          !
    !   is_periodic (input) : Logical array indicating whether the grid is    !
    !                         periodic (T) or not (F) along each lattice      !
    !                         vector. This is used to determine how the       !
    !                         correction is computed.                         !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    ! Modified by James C. Womack to support periodic BCs, 05/2017            !
    ! Modified to remove pub_par by Joseph Prentice, 09/2020                  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id
    use constants, only: DP
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+), OPERATOR(*)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_is_smeared_ion_width
!$  use rundat, only: pub_threads_max
    use utils, only: utils_abort, utils_assert, utils_erf

    implicit none

    ! jd: Arguments
    type(PARAL_INFO), intent(in) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: v_loc_on_fine(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    logical, intent(in) :: is_periodic(3)

    ! jd: Local variables and parameters
    real(kind=DP), parameter :: sqrt_pi = 1.7724538509055_DP
    real(kind=DP) :: Z_I, d, sigma_I, term, accum
    integer :: ipt, i1, i2, i3, islab12, I
    type(POINT) :: R_I, r

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_apply_vloc_corr')

    if (.not.any(is_periodic)) then
       ! JCW: Fully open BCs
       ! JCW: Compute correction by summing over per-ion potentials.
       ! jd: Add the correction to v_loc_on_fine
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,i3,r,accum,I,R_I,Z_I,sigma_I,d,term) &
!$OMP SHARED(par,elements,grid,pub_is_smeared_ion_width, &
!$OMP      v_loc_on_fine,pub_my_proc_id,pub_threads_max)
       do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
          i1 = modulo(ipt-1,grid%n1) + 1
          i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
          islab12 = (ipt-(i2-1)*grid%n1-i1) / &
               (grid%n1*grid%n2) + 1
          i3 = grid%first_slab12(pub_my_proc_id) + islab12 - 1

          r = &
               real((i1-1),kind=DP) * grid%da1 + &
               real((i2-1),kind=DP) * grid%da2 + &
               real((i3-1),kind=DP) * grid%da3

          accum = 0.0_DP
          do I = 1, par%nat

             R_I = elements(I)%centre
             Z_I = real(elements(I)%ion_charge,kind=DP)
             sigma_I = pub_is_smeared_ion_width

             d = magnitude(r-R_I)

             ! jd: Calculate erf(d/sigma_I)/d, considering d==0 separately
             if (d < 1E-10_DP) then
                ! jd: lim_{d->0} erf(d/sigma_I)/d = 2/(sqrt(pi)*sigma_I)
                term = 2.0_DP/sigma_I/sqrt_pi
             else
                term = utils_erf(d/sigma_I) / d
             end if

             ! jd: It's faster to accumulate in a variable first
             accum = accum + Z_I * term

          end do ! over I

          ! jd: Subtract the whole thing from current point in v_loc_fine
          v_loc_on_fine(i1,i2,islab12) = v_loc_on_fine(i1,i2,islab12) + accum

       end do ! ipt
!$OMP END PARALLEL DO

    else if (all(is_periodic)) then
       ! JCW: Fully periodic BCs
       ! JCW: Use the previously computed periodic potential due to the
       ! JCW: smeared ions to correct the local pseudopotential.

       ! JCW: phi_ion must be allocated and initialised by smeared_ion_initialise
       call utils_assert(allocated(phi_ion),"Error in smeared_ion_apply_vloc_corr: &
            &For fully periodic BC correction of the local pseudopotential, phi_ion &
            &should be first allocated and initialised by smeared_ion_initialise.")

       v_loc_on_fine(:,:,:) = v_loc_on_fine(:,:,:) - phi_ion(:,:,:)

    else
       ! JCW: Mixed open/periodic BCs
       ! JCW: Not currently supported
       call utils_abort("Error in smeared_ion_apply_vloc_corr: &
            &Mixed BCs for non-cutoff-Coulomb calculations are not yet supported.")

    end if ! BCs

    call utils_trace_out('smeared_ion_apply_vloc_corr')

  end subroutine smeared_ion_apply_vloc_corr
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_generate_density(elements,grid,cell,is_periodic,nat)
    !=========================================================================!
    ! Generates the smeared-ion density, stores it into rho_ion. The          !
    ! way the density is constructed depends upon the boundary conditions.    !
    !                                                                         !
    ! Fully open BCs:                                                         !
    !   Smeared Gaussian charges only exist in the simulation cell. To avoid  !
    !   truncation of the charge, the Gaussians must be sufficiently far from !
    !   the cell boundaries that they do not bleed into adjacent cells.       !
    !                                                                         !
    ! Fully periodic BCs:                                                     !
    !   Smeared Gaussian charges are made approximately periodic by adding    !
    !   contributions from each Gaussian using a minimum-image convention.    !
    !   Only contributions from a single Gaussian are considered (not images) !
    !   images, but components of that Gaussian "wrap around" at cell         !
    !   boundaries.                                                           !
    !                                                                         !
    ! Note on periodic smeared ion representation:                            !
    ! Since the smeared Gaussian functions are simply a useful representation !
    ! of point charges for use in multigrid operations, it does not matter    !
    ! that they are not "truly" periodic in the sense that periodic images    !
    ! do not contribute. The smeared ion-smeared ion and smeared ion-electron !
    ! interactions are corrected for and replaced with the true ion-ion       !
    ! (Ewald) and ion-electron (local pseudopotential) interactions when      !
    ! evaluating the total energy.                                            !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   elements(input):      The array describing the ions.                  !
    !   is_periodic (input) : Logical array indicating whether the grid is    !
    !                         periodic (T) or not (F) along each lattice      !
    !                         vector. This is used to determine how the       !
    !                         smeared ion density is constructed.             !
    !   grid(input):          The fine grid.                                  !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    ! Modified by James C. Womack to support periodic BCs, 05/2017            !
    ! Smooth-walled exclusion regions added by Edward B. Linscott and         !
    ! Kevin K. B. Duff, March 2018.                                           !
    ! Modified for embedding and to remove pub_par by Joseph Prentice, 09/2020!
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id, pub_on_root
    use constants, only: pi, stdout, VERBOSE
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+), OPERATOR(*)
    use ion, only: ELEMENT
    use rundat, only: pub_is_smeared_ion_width, pub_is_core_width, &
          pub_output_detail
!$  use rundat, only: pub_threads_max
    use simulation_cell, only: CELL_INFO, minimum_image_distance_ortho
    use utils, only: utils_abort, utils_assert

    implicit none

    ! jd: Arguments
    integer, intent(in) :: nat
    type(ELEMENT), intent(in) :: elements(nat)
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    logical, intent(in) :: is_periodic(3)

    ! jd: Local variables
    integer :: islab12, i1, i2, i3, ipt, I
    type(POINT) :: r, R_I
    real(kind=DP) :: sigma_I, Z_I, d, accum, prefactor
    real(kind=DP) :: minx, maxx, miny, maxy, minz, maxz, dx, dy, dz

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_generate_density')

    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
       write(stdout,'(/a)',advance='no') 'Generating smeared ion density ...'
    end if
    prefactor = 1.0_DP/pi**1.5_DP
    rho_ion = 0.0_DP ! jd: Take care of padding between pt and ld
    dielectric_mask = 0.0_DP ! ebl: By default everywhere

    ! JCW: Check for supported BCs
    if (any(is_periodic).and.(.not.all(is_periodic))) then
       ! JCW: Mixed periodic/open BCs
       ! JCW: Not currently supported
       call utils_abort("Error in smeared_ion_generate_density: &
            &Mixed periodic/open BCs are not currently supported for smeared ions.")
    end if ! BCs

    ! jd: Generate smeared ion density on every point of fine grid
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,i3,r,accum,I,R_I,Z_I,sigma_I,d,minx,maxx,miny, &
!$OMP      maxy, minz, maxz, dx, dy, dz) &
!$OMP SHARED(nat,cell,elements,grid,is_periodic,pub_is_smeared_ion_width, &
!$OMP      pub_is_core_width,dielectric_mask,rho_ion,prefactor, &
!$OMP      pub_my_proc_id,pub_threads_max,dielectric_exclusions)
    do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / &
            (grid%n1*grid%n2) + 1
       i3 = grid%first_slab12(pub_my_proc_id) + islab12 - 1

       r = &
            real((i1-1),kind=DP) * grid%da1 + &
            real((i2-1),kind=DP) * grid%da2 + &
            real((i3-1),kind=DP) * grid%da3

       accum = 0.0_DP
       ! JCW: No need to check BCs are supported, since we have already checked this
       do I = 1, nat
          R_I = elements(I)%centre
          !ab: electronic convention, ionic charges are treated as negative
          Z_I = -elements(I)%ion_charge
          sigma_I = pub_is_smeared_ion_width
          if (.not.any(is_periodic)) then
             ! JCW: Fully open BCs
             ! JCW: Gaussians do not "wrap around" at cell boundaries, distance
             ! JCW: between points is simple magnitude of difference of vectors
             d = magnitude(r-R_I)
          else if (all(is_periodic)) then
             ! JCW: Fully periodic BCs
             ! JCW: Gaussians wrap around at cell boundaries using minimum image
             ! JCW: convention
             d = minimum_image_distance_ortho(r,R_I,cell)
          end if
          if(d <= pub_is_core_width) then
             dielectric_mask(i1,i2,islab12) = 1.0_DP
          end if
          ! jd: Calculate density value for d, accumulate
          accum = accum + Z_I / sigma_I**3 * exp(-d*d/(sigma_I*sigma_I))
       end do ! over I

       ! jd: Store
       rho_ion(i1,i2,islab12) = prefactor * accum

       ! jd: Handle exclusions from exclusion regions
       do I = 1, n_max_exclusions
          ! jd: --- no more exclusions ---
          if(dielectric_exclusions(I,1) == 0) then
             exit

          ! jd: @PBC@ will need to be supported here once we have solvation in PBC.
          ! JCW: Block execution with PBCs until implemented
          if (any(is_periodic)) call utils_abort("Error in &
               &smeared_ion_generate_density: dielectric_exclusions are not yet &
               &supported with PBCs.")

          ! jd: --- sphere ---
          else if(dielectric_exclusions(I,1) == 1) then
             R_I%x = dielectric_exclusions(I,2)
             R_I%y = dielectric_exclusions(I,3)
             R_I%z = dielectric_exclusions(I,4)
             d = magnitude(r-R_I)

             ! jd: Exclude what's in the sphere
             ! ebl: change to smeared sphere
             call smeared_ion_smooth_mask(dielectric_mask(i1,i2,islab12), &
                  d, dielectric_exclusions(I,5))

          ! jd: --- box ---
          else if(dielectric_exclusions(I,1) == 2) then
             minx = dielectric_exclusions(I,2)
             maxx = dielectric_exclusions(I,3)
             miny = dielectric_exclusions(I,4)
             maxy = dielectric_exclusions(I,5)
             minz = dielectric_exclusions(I,6)
             maxz = dielectric_exclusions(I,7)
             dx = max(minx-r%x, r%x-maxx)
             dy = max(miny-r%y, r%y-maxy)
             dz = max(minz-r%z, r%z-maxz)
             if(dx < 0.0_DP .and. dy < 0.0_DP .and. dz < 0.0_DP) then
                ! ebl: if inside box, give a negative distance
                d = max(dx, dy, dz)
             else
                if (dx < 0.0_DP) dx = 0.0_DP
                if (dy < 0.0_DP) dy = 0.0_DP
                if (dz < 0.0_DP) dz = 0.0_DP
                d = sqrt(dx**2 + dy**2 + dz**2)
             end if

             ! jd: Exclude what's in the box
             ! ebl: changed to smeared box
             call smeared_ion_smooth_mask(dielectric_mask(i1,i2,islab12), &
                  d, 0.0_DP)
          ! jd: --- xcyl ---
          else if(dielectric_exclusions(I,1) == 3) then
             R_I%x = r%x
             R_I%y = dielectric_exclusions(I,2)
             R_I%z = dielectric_exclusions(I,3)
             d = magnitude(r-R_I)

             ! jd: Exclude what's in the cylinder
             ! ebl: changed to smeared cylinder
             call smeared_ion_smooth_mask(dielectric_mask(i1,i2,islab12), &
                  d, dielectric_exclusions(I,4))
          else if(dielectric_exclusions(I,1) == 4) then
             R_I%x = dielectric_exclusions(I,2)
             R_I%y = r%y
             R_I%z = dielectric_exclusions(I,3)
             d = magnitude(r-R_I)

             ! jd: Exclude what's in the cylinder
             ! ebl: changed to smeared cylinder
             call smeared_ion_smooth_mask(dielectric_mask(i1,i2,islab12), &
                  d, dielectric_exclusions(I,4))
          else if(dielectric_exclusions(I,1) == 5) then
             R_I%x = dielectric_exclusions(I,2)
             R_I%y = dielectric_exclusions(I,3)
             R_I%z = r%z
             d = magnitude(r-R_I)

             ! jd: Exclude what's in the cylinder
             ! ebl: changed to smeared cylinder
             call smeared_ion_smooth_mask(dielectric_mask(i1,i2,islab12), &
                  d, dielectric_exclusions(I,4))
          else
             call utils_assert(.false., &
                  'Unrecognized dielectric exclusion type: ', &
                  dielectric_exclusions(I,1))
          end if

       end do
    end do ! ipt
!$OMP END PARALLEL DO

    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
       write(stdout,'(a)') ' done'
    end if

    call utils_trace_out('smeared_ion_generate_density')

  end subroutine smeared_ion_generate_density
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine smeared_ion_smooth_mask(mask, d, d0)
    !========================================================================!
    ! Smears the dielectric mask using a Fermi-Dirac function                !
    ! ---------------------------------------------------------------------- !
    ! Arguments:                                                             !
    !   mask (in/output) : a point in the dielectric mask                    !
    !   d    (input)     : the distance of that point from the centre of a   !
    !                      sphere/axis of a cylinder/edge of a box           !
    !                      (analogous to the energy)                         !
    !   d0   (input)     : the distance to the surface of the exclusion      !
    !                      region (analogous to the Fermi energy)            !
    ! ---------------------------------------------------------------------- !
    ! Written by Edward B. Linscott and Kevin K. B. Duff, 02/2018, partly    !
    ! based on edft_fermidirac_distr                                         !
    !========================================================================!

    use constants, only: SAFE_DIV_EPS
    use rundat, only: pub_is_dielectric_exclusions_smearing

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: mask
    real(kind=DP), intent(in   ) :: d
    real(kind=DP), intent(in   ) :: d0

    ! Local variables
    real(kind=DP) :: fermi_distr
    real(kind=DP) :: exp_arg, exp_val

    ! ebl: protecting against dividing by zero
    if(pub_is_dielectric_exclusions_smearing > SAFE_DIV_EPS) then
       ! ars: check for DP overflow
       exp_arg = (d - d0) / pub_is_dielectric_exclusions_smearing
       if(exp_arg < 700.0_DP) then ! exp(700) is just below max_DP
          exp_val = exp(exp_arg)
          fermi_distr = 1.0_DP / (1.0_DP + exp_val)
       else
          fermi_distr = 0.0_DP
       end if
    else
       if(d > d0) then
          fermi_distr = 0.0_DP
       else
          fermi_distr = 1.0_DP
       end if
    end if

    ! ebl: ensure mask stays within [0, 1]
    mask = min(1.0_DP, mask + fermi_distr)

  end subroutine smeared_ion_smooth_mask


  subroutine smeared_ion_generate_potential_PBC(phi_si,rho_si,grid)
    !=========================================================================!
    ! Generates the electrostatic potential due to smeared ion density in     !
    ! fully periodic BCs by solving the Poisson equation in reciprocal space. !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   phi_si (output) : The electrostatic potential on the real space grid  !
    !   rho_si (input)  : The smeared ion charge density on the grid          !
    !   grid   (input)  : GRID_INFO instance describing the grids to be used  !
    !                     to store the charge and potential and also for the  !
    !                     reciprocal space work array                         !
    ! ----------------------------------------------------------------------- !
    ! Written by James C. Womack, 05/2017, based on hartree_on_grid           !
    !=========================================================================!

    use constants, only: PI
    use cell_grid, only: GRID_INFO
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out) :: phi_si(grid%ld1,grid%ld2,grid%max_slabs12)
    real(kind=DP), intent(in)  :: rho_si(grid%ld1,grid%ld2,grid%max_slabs12)


    ! Parameters
    character(len=*), parameter :: myself = "smeared_ion_generate_potential_PBC"
    real(kind=DP), parameter :: fourpi = 4.0_DP*PI ! 4 pi

    ! Local variables
    integer :: ierr                                ! Error flag
    complex(kind=DP), allocatable :: zwork(:,:,:)  ! Workspace

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    ! Allocate workspace
    allocate(zwork(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check(myself,'zwork',ierr)

    ! Fourier transform smeared ion density to reciprocal space
    call fourier_apply_cell_forward(rho_si(:,:,:),zwork,grid)

    ! Compute electrostatic potential based on reciprocal space solution of the
    ! Poisson equation, i.e.
    !   \phi(G) = 4\pi * n(G) / G^2
    zwork = zwork * grid%coulomb_recip(:,:,:) * fourpi

    ! Fourier transform potential to real space
    call fourier_apply_cell_backward(phi_si(:,:,:),zwork,grid)

    ! Deallocate workspace
    deallocate(zwork,stat=ierr)
    call utils_dealloc_check(myself,'zwork',ierr)

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine smeared_ion_generate_potential_PBC

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_hartree(phi, rho_elec, grid, cell, E_Hartree, no_dump)
    !=========================================================================!
    ! Calculates the Hartree energy in the smeared-ion representation.        !
    ! This is the molecular Hartree energy as described in [1], Appendix.     !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   phi(input):      Molecular electrostatic potential on the fine grid.  !
    !   rho_elec(input): *Electronic* charge density on the fine grid.        !
    !                    (smeared-ion density is stored by the module.)       !
    !   E_Hartree(output), optional: Resultant energy.                        !
    !   no_dump, (opt)     intent=in,  logical parameter to avoid 3D dumps    !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use multigrid_methods, only: multigrid_calculate_hartree
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    ! jd: Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    ! jd: Electronic density in real space:
    real(kind=DP), intent(inout) :: rho_elec(grid%ld1, &
         grid%ld2,grid%max_slabs12)

    ! jd: Resultant molecular potential, from the MG
    real(kind=DP), intent(out) :: phi(grid%ld1,grid%ld2, &
         grid%max_slabs12)

    ! jd: Calculated Hartree energy, if required
    real(kind=DP), intent(out), optional :: E_Hartree
    logical, intent(in), optional :: no_dump

    ! jd: Local variables
    real(kind=DP), allocatable :: rho_tot(:,:,:)
    real(kind=DP) :: E_Hartree_loc

    integer :: ierr ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_hartree')

    ! jd: Allocate the array that will hold smeared ion density
    allocate(rho_tot(grid%ld1,grid%ld2, &
         grid%max_slabs12),stat=ierr)
    call utils_alloc_check('smeared_ion_hartree','rho_tot',ierr)

    ! jd: Prepare a molecular density
    rho_tot = rho_ion + rho_elec

    ! jd: Calculate Hartree potential and energy
    E_Hartree_loc = multigrid_calculate_hartree(phi, rho_tot, grid, cell, &
    no_dump=no_dump)
    if(present(E_Hartree)) E_Hartree = E_Hartree_loc

    ! jd: Clean up
    deallocate(rho_tot,stat=ierr)
    call utils_dealloc_check('smeared_ion_hartree','rho_tot',ierr)

    call utils_trace_out('smeared_ion_hartree')

  end subroutine smeared_ion_hartree
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_density(rho_tot, rho_elec, grid)
    !=========================================================================!
    ! Calculates the density in the smeared-ion representation.               !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   rho_elec(input): *Electronic* charge density on the fine grid.        !
    !                    (smeared-ion density is stored by the module.)       !
    ! ----------------------------------------------------------------------- !
    ! Modified from smeared_ion_hartree by Jolyon Aarons...2017               !
    !=========================================================================!

    use cell_grid, only: GRID_INFO

    ! jd: Arguments
    type(GRID_INFO), intent(in) :: grid
    ! jd: Electronic density in real space:
    real(kind=DP), intent(in) :: rho_elec(grid%ld1, &
         grid%ld2,grid%max_slabs12)

    ! jd: Resultant molecular density
    real(kind=DP), intent(out) :: rho_tot(grid%ld1,grid%ld2, &
         grid%max_slabs12)

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_density')

    ! jd: Prepare a molecular density
    rho_tot = rho_ion + rho_elec

    call utils_trace_out('smeared_ion_density')

  end subroutine smeared_ion_density
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_ii_forces(ii_forces,mdl)
    !=====================================================================!
    ! Calculates the direct smeared_ion-smeared_ion contribution to the   !
    ! energy. Replaces ewald_calculate_forces for finite systems that use !
    ! the smeared ion representation.                                     !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    !   elements (input): list of all the atoms in the system.            !
    !   ii_energy (output): all smeared_ion-smeared_ion forces.           !
    !---------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 09/2013 using cutoff_coulomb_ii_forces,  !
    ! written by Nicholas Hine, March 2008, as template.                  !
    ! Modified for embedding and to remove pub_par by Joseph Prentice,    !
    ! May 2018                                                            !
    ! Further modification for embedding by Joseph Prentice, Sep 2020     !
    !=====================================================================!

    use comms, only: comms_reduce, pub_my_proc_id
    use constants, only: SQRT_PI
    use geometry, only: POINT, OPERATOR(+), OPERATOR(-), OPERATOR(*), magnitude
    use model_type, only: MODEL
    use rundat, only: pub_is_smeared_ion_width
    use utils, only: utils_erf

    implicit none

    ! Arguments
    type(MODEL),intent(in) :: mdl
    real (kind=dp), dimension(1:3,1:mdl%nat), intent(out) :: ii_forces

    ! Locals
    integer       :: my_first_at, my_last_at  ! range of atoms on this proc
    integer       :: iatom, jatom             ! atom index counters
    integer       :: isub                     ! region counter
    integer       :: global_i                 ! global atom index counter
    real(kind=DP) :: qi, qj                   ! charge of atoms i and j
    real(kind=DP) :: fij(3)     ! force between smeared ions i and j
    type(POINT)   :: rijv       ! vector between particles i and j
    real(kind=DP) :: rij        ! |ri-rj|
    real(kind=DP) :: rij2       ! |ri-rj|^2
    real(kind=DP) :: rij3       ! |ri-rj|^3
    real(kind=DP), parameter :: two_over_sqrt_pi = 2.0_DP/SQRT_PI
    real(kind=DP) :: sigma_I, sigma_J, sigma_IJ ! smearing widths
    real(kind=DP) :: scalar_term

    ! -------------------------------------------------------------------------

    ii_forces(:,:) = 0.0_DP

    ! jcap: loop over regions
    do isub = 1, mdl%nsub

       my_first_at = mdl%regions(isub)%par%first_atom_on_proc(pub_my_proc_id)
       my_last_at = my_first_at + mdl%regions(isub)%par%num_atoms_on_proc(pub_my_proc_id) - 1

       ! Loop over my atoms
       do iatom = my_first_at, my_last_at

          ! jcap: find the global atom number
          global_i=mdl%regions(isub)%elements(iatom)%global_atom_number

          qi = real(mdl%elements(global_i)%ion_charge, kind=DP)

          sigma_I = pub_is_smeared_ion_width

          ! Loop over all other atoms j<i
          do jatom = 1, global_i - 1

             sigma_J = pub_is_smeared_ion_width ! might depend on J in the future

             sigma_IJ = sqrt(sigma_I*sigma_I + sigma_J*sigma_J)

             ! Find charge and distance of other atom
             qj = real(mdl%elements(jatom)%ion_charge,kind=DP)
             rijv = mdl%elements(global_i)%centre - mdl%elements(jatom)%centre
             rij  = magnitude(rijv)
             rij2 = rij*rij
             rij3 = rij2*rij

             ! Find force between atoms i and j
             scalar_term = (qi*qj/rij3) * &
                  (1.0_DP &
                  + two_over_sqrt_pi/sigma_IJ*exp(-rij2/(sigma_IJ*sigma_IJ)) * rij&
                  - utils_erf(rij/sigma_IJ))
             fij(1) = scalar_term * rijv%X
             fij(2) = scalar_term * rijv%Y
             fij(3) = scalar_term * rijv%Z

             ! Add contribution of jatom to force on iatom
             ii_forces(:,global_i) = ii_forces(:,global_i) + fij

             ! Add contribution of iatom to force on jatom
             ii_forces(:,jatom) = ii_forces(:,jatom) - fij

          enddo

       enddo

    end do

    call comms_reduce('SUM', ii_forces)

  end subroutine smeared_ion_ii_forces
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_force_correction(phi_vac, grid, elements, nat, &
       smeared_ion_forces)
    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the smeared ions interacting with other smeared ions and with      !
    ! electrons, in vacuum. These are forces due to terms (4B) and (4C) in    !
    ! solvation forces notes.                                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi_vac (in): The potential due to the electrons + smeared_ions,      !
    !                 obtained from a multigrid calculation in vacuum.        !
    !   grid (in): The integration grid.                                      !
    !   elements (in): The usual.                                             !
    !   smeared_ion_forces (out): The calculated forces on all smeared ions.  !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 09/2013 using code snippets from           !
    ! pseudo_local_calculate_forces.                                          !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id, comms_reduce
    use constants, only: PI, SQRT_PI, stdout
    use geometry, only: magnitude, POINT, operator(-), operator(*)
! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
    use geometry, only: add_points
#else
    use geometry, only: operator(+)
#endif
    use ion, only: ELEMENT
    use rundat, only: pub_forces_needed, pub_is_smeared_ion_width, &
         pub_debug_on_root
!$  use rundat, only: pub_threads_max
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)  :: grid
    real(kind=DP), intent(inout) :: phi_vac(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    integer, intent(in)          :: nat
    type(ELEMENT), intent(in)    :: elements(nat)
    real(kind=DP), intent(out)   :: smeared_ion_forces(1:3,nat)

    ! Local Variables
    integer       :: ipt, i1, i2, i3, islab12
    integer       :: I
    real(kind=DP) :: x
    real(kind=DP) :: scalar_term
    type(POINT)   :: r
    type(POINT)   :: R_I
    type(POINT)   :: r_minus_R_I
    real(kind=DP) :: sigma_I
    real(kind=DP) :: Z_I
!$  real(kind=DP) :: term
    real(kind=DP), parameter :: factor = 2.0_DP/(PI*SQRT_PI)

    ! -------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering smeared_ion_force_correction'

    ! Start timer
    call timer_clock('smeared_ion_force_correction',1)

    call utils_assert(pub_forces_needed, &
         'Internal error in smeared_ion_force_correction.')

    smeared_ion_forces(:,:) = 0.0_DP

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,i3,r,I,R_I,r_minus_R_I,x,scalar_term,term,Z_I, &
!$OMP      sigma_I) &
!$OMP SHARED(nat,elements,grid, pub_is_smeared_ion_width, &
!$OMP      pub_my_proc_id,pub_threads_max,phi_vac) &
!$OMP REDUCTION(+:smeared_ion_forces)
!$OMP DO

    ! jd: Loop over local grid points
    do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / &
            (grid%n1*grid%n2) + 1
       i3 = grid%first_slab12(pub_my_proc_id) + islab12 - 1

! jme: workaround for cray compiler bug 843178
#if _CRAYFTN
       r = add_points(add_points((i1-1)*grid%da1, (i2-1)*grid%da2), &
                      (i3-1)*grid%da3)
#else
       r = (i1-1) * grid%da1 + (i2-1) * grid%da2 + (i3-1) * grid%da3
#endif

       ! jd: Loop over all the ions
       do I = 1, nat

          R_I = elements(I)%centre
          r_minus_R_I = r - R_I
          x = magnitude(r_minus_R_I)
          Z_I = real(elements(I)%ion_charge,kind=DP)
          sigma_I = pub_is_smeared_ion_width

          scalar_term = factor * phi_vac(i1,i2,islab12) * Z_I &
               * sigma_I**(-5.0_DP) * exp(-x*x/(sigma_I*sigma_I))

          smeared_ion_forces(1,I) = smeared_ion_forces(1,I) &
               + scalar_term * r_minus_R_I%X
          smeared_ion_forces(2,I) = smeared_ion_forces(2,I) &
               + scalar_term * r_minus_R_I%Y
          smeared_ion_forces(3,I) = smeared_ion_forces(3,I) &
               + scalar_term * r_minus_R_I%Z

       end do ! atoms, I

    end do ! ipt
!$OMP END DO
!$OMP END PARALLEL

    ! jd: Apply grid weight and sum over grid slabs
    smeared_ion_forces = smeared_ion_forces * grid%weight
    call comms_reduce('SUM', smeared_ion_forces)

    ! Stop timer
    call timer_clock('smeared_ion_force_correction',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving smeared_ion_force_correction'

  end subroutine smeared_ion_force_correction
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module is_smeared_ions
