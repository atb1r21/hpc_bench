!=============================================================================!
!=============================================================================!
!                          MULTIGRID_METHODS MODULE                           !
!=============================================================================!
! This module defines the onetep interface to the multigrid solver and        !
! subroutines to calculate electrostatic quantities that are related to the   !
! solution of the Poisson equation.                                           !
!                                                                             !
! For more info on multigrid methods see the following useful references      !
!                                                                             !
! [1] 'A Multigrid Tutorial' Briggs, Henson and McCormick.                    !
!                                                                             !
! [2] TL Beck 'Real-space mesh techniques in density-functional theory'       !
!     Rev Mod Phys 72 pg 1041.                                                !
!                                                                             !
! Portions of this module use techniques or formulas described in             !
! [3] Scherlis, Fattebert, Gygi, Cococcioni and Marzari, J. Chem. Phys. 124,  !
!     074103 (2006).                                                          !
!-----------------------------------------------------------------------------!
! Written for castep by Hatem H Helal, starting in 10/2007                    !
! Adapted for onetep by Jacek Dziedzic, 04-05/2010.                           !
! Made parallel-ready by Jacek Dziedzic, 04-05/2010.                          !
!                                                                             !
! Implementation of defect correction solver by HH Helal 06/2010              !
! Described in some detail within the following references:                   !
! [4] U Trottenberg, C Oosterlee and A Schuller, "Multigrid", 2001.           !
! [5] S Schaffer, Math of Comp, 43, 89-115 (1984).                            !
!-----------------------------------------------------------------------------!

module multigrid_methods

  use constants, only: DP, stdout
  use finite_differences, only: FD_GRID_INFO
#ifdef HAVE_DL_MG
  use dl_mg  !! External dependency
#endif

  implicit none

  private

  !---------------------------------------------------------------------------!
  !                   P u b l i c   D a t a s t r u c t u r e s               !
  !---------------------------------------------------------------------------!

  ! jd: Geometry of the multigrid
  !     Mostly used privately, but is_solvation needs it to calculate FD
  !     gradients on the MG grid. Boltzmann solvation needs it to dimension
  !     steric_pot and gamma.
  type(FD_GRID_INFO), public, save :: mg

  ! jd: Last potential calculated through the multigrid. This is needed for
  !     the calculation of force correction due to smeared ions. During the
  !     force calculation phi is needed, and the total density is no longer
  !     available, so we resort to storing phi from the last multigrid_hartree
  !     calculation. This array is not allocated if smeared ions are not in use
  !     or if forces are not needed.
  real(kind=DP), public, allocatable, save :: phi_vac_for_forces(:,:,:)

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!

  public  ::  multigrid_initialise
  public  ::  multigrid_exit
  public  ::  multigrid_calculate_hartree
  public  ::  multigrid_electrostatic_energy
  public  ::  multigrid_prepare_bound_cond
  public  ::  multigrid_bc_for_dlmg

  !---------------------------------------------------------------------------!
  !                               P r i v a t e                               !
  !---------------------------------------------------------------------------!

  ! jd: Initialisation flags
  logical, save :: multigrid_initialised = .false.

  ! jd: Debug flag, set to true by DEVEL_CODE MG: ... :MG
  logical, save :: multigrid_debug = .false.

  ! JCW: Enable use of (original) defect correction code in ONETEP
  ! JCW: This is set in multigrid_initialise, and is .true. if the DEVEL_CODE
  ! JCW: MG:USE_ONETEP_DEFCO=T:MG is present in the input file.
  logical, save :: multigrid_use_onetep_defco = .false.

contains

  !---------------------------------------------------------------------------!
  !                      P u b l i c   s u b r o u t i n e s                  !
  !---------------------------------------------------------------------------!

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine multigrid_initialise(grid)
    !=========================================================================!
    ! Initialises the multigrid.                                              !
    ! Calculates the size of the multigrid, sanity-checks that the            !
    ! distribution makes sense, examines DEVEL_CODE.                          !
    ! There is no need to call this subroutine, it will be automatically      !
    ! called the first time you use any of multigrid_*() routines, but you can!
    ! call it explicitly, any number of times, if so desired.                 !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                  !
    ! Updated by James C. Womack, 02/2017 to reflect new DL_MG interface.     !
    !=========================================================================!

    use bibliography, only: bibliography_cite
    use cell_grid, only: GRID_INFO
    use constants, only: stdout
    use comms, only: pub_on_root, pub_root_proc_id, pub_total_num_procs, &
         comms_bcast, comms_allgather, comms_cart_create, pub_null_comm, &
         comms_rank
    use finite_differences, only: finite_difference_initialise, &
         finite_difference_set_geometry
    use geometry, only: operator(.DOT.)
    use rundat, only: pub_devel_code, pub_forces_needed, &
         pub_mg_granularity_power, pub_mg_continue_on_error, &
         pub_is_smeared_ion_rep, pub_multigrid_bc_is_periodic, pub_rootname
    use utils, only: utils_assert, utils_alloc_check, &
         utils_dealloc_check, utils_devel_code, utils_feature_not_supported, &
         utils_flush, utils_trace_in, utils_trace_out, utils_unit
#if defined (__INTEL_COMPILER) && defined (DEBUG)
!$  use constants, only: CRLF
!$  use omp_lib
#endif

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid

    ! jd: Local variables
    integer            :: multigrid_granularity
    character(len=200) :: mg_devel_code
    integer            :: start_pos, stop_pos, test_pos
    integer            :: bc(3) ! DL_MG-compatible BC specification
    integer            :: ierr
#ifdef HAVE_DL_MG
    integer              :: n(3)
    real(kind=DP)        :: d(3)
    integer              :: mgunit, mg_comm
    integer              :: idx_start(3), idx_end(3), mg_proc_id
    integer, allocatable :: zlim(:)
#endif
    character(len=*), parameter :: myself = 'multigrid_initialise'

#if defined (__INTEL_COMPILER) && defined (DEBUG)
    ! Used in workaround for possible compiler bug with Intel compiler in
    ! debug mode with > 1 thread
!$  integer              :: nthreads_save
#endif

    ! ------------------------------------------------------------------------

    if(multigrid_initialised) return

    call utils_trace_in(myself)
    call bibliography_cite('PBCs')

    ! jd: @warning: changes to this should be reflected in
    !     is_solvation_boltzmann::init_steric_pot
    multigrid_granularity = 2**pub_mg_granularity_power

    ! JCW: Check that grid is orthorhombic since this is required by
    ! JCW: the multigrid solver
    call utils_assert(&
         (grid%da1 .dot. grid%da2) == 0.0_DP .and. &
         (grid%da1 .dot. grid%da3) == 0.0_DP .and. &
         (grid%da2 .dot. grid%da3) == 0.0_DP, &
         'Error in multigrid_initialise: &
         &The multigrid solver only supports orthorhombic cells.')

    ! JCW: Derive DL_MG compatible BC array from pub_multigrid_bc_is_periodic
    call multigrid_bc_for_dlmg(bc,pub_multigrid_bc_is_periodic)

    ! jd: Initialise the MG grid
    call finite_difference_initialise
    mg = finite_difference_set_geometry(grid,multigrid_granularity,pub_multigrid_bc_is_periodic)

    ! jd: Check if in devel-code low-level debug mode
    if (pub_on_root) then
       mg_devel_code=pub_devel_code
       if (len_trim(mg_devel_code)>0) then
          start_pos=index(mg_devel_code,'MG:')
          stop_pos=index(mg_devel_code,':MG')
          if (stop_pos<=0) stop_pos=len_trim(mg_devel_code)
          if (start_pos>0) then
             test_pos=index(mg_devel_code,'DEBUG_DUMP')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                multigrid_debug = .true.
                write(stdout,'(a)') &
                     'NB: Multigrid debugging dumps will be produced.'
             end if
          end if
       end if
    end if

    ! JCW: Output some information regarding grid sizes and warn if
    ! JCW: grid truncation for multigrid is potentially dangerous
    if(pub_on_root) then
       write(stdout,'(a)',advance='no') 'ONETEP fine grid is '
       write(stdout,'(i0,a,i0,a,i0)',advance='no') &
            mg%pt1f,' x ',mg%pt2f,' x ',mg%pt3f
       write(stdout,'(a,F0.4,a,F0.4,a,F0.4,a)') ' gridpoints, ', &
            mg%pt1f*mg%d1f, ' x ',mg%pt2f*mg%d2f,' x ', mg%pt3f*mg%d3f,' bohr.'
       write(stdout,'(a)',advance='no') 'FD multigrid is     '
       write(stdout,'(i0,a,i0,a,i0)',advance='no') &
            mg%pq1f,' x ',mg%pq2f,' x ',mg%pq3f
       write(stdout,'(a,F0.4,a,F0.4,a,F0.4,a/)') ' gridpoints, ', &
            mg%pq1f*mg%d1f, ' x ',mg%pq2f*mg%d2f,' x ', mg%pq3f*mg%d3f,' bohr.'

       if(mg%pt1f-mg%pq1f>10 .or. mg%pt2f-mg%pq2f>10 .or. &
            mg%pt3f-mg%pq3f>10) then
          write(stdout,'(a)') 'WARNING: A considerable portion of the fine grid&
               & will be hidden from the'
          write(stdout,'(a)') 'multigrid calculation. If you''re sure that&
               & there is no charge density in the '
          write(stdout,'(a)') 'margin that constitues the above difference,&
               & you''re OK. If not, try to increase '

          if(mg%pt1f-mg%pq1f>10) then
             write(stdout,'(a,i0,a,i0,a)') 'the 1st dimension, ',mg%pt1f, &
                  ', to ', mg%pq1f+multigrid_granularity+1,','
          end if
          if(mg%pt2f-mg%pq2f>10) then
             write(stdout,'(a,i0,a,i0,a)') 'the 2nd dimension, ',mg%pt2f, &
                  ', to ', mg%pq2f+multigrid_granularity+1,','
          end if
          if(mg%pt3f-mg%pq3f>10) then
             write(stdout,'(a,i0,a,i0,a)') 'the 3rd dimension, ',mg%pt3f, &
                  ', to ', mg%pq3f+multigrid_granularity+1,','
          end if
          write(stdout,'(a)') 'and re-run.'
       end if

    end if

    call comms_bcast(pub_root_proc_id,multigrid_debug)

    call utils_flush

    ! la: generate a MPI topology and initialise dl_mg multigrid solver
#ifdef HAVE_DL_MG

    ! la: dl_mg needs a communicator with a topology and the index range of the
    ! la: global grid hold on this rank ( no halos, but do include the boundary layer )
    call comms_cart_create(mg_comm, pub_multigrid_bc_is_periodic, mg%num_slabs12)

    ! JCW: Note on communicators
    ! * The mg_comm Cartesian topology communicator features only MPI ranks
    !   with num_slabs12 > 0.
    ! * MPI ranks with num_slabs12 == 0  have mg_comm = MPI_COMM_NULL.
    ! * DL_MG checks whether mg_comm == MPI_COMM_NULL on entry, and exits
    !   without doing any work for both init and solver routines in this
    !   case.

    ! JCW: Check devel_code string for MG:USE_ONETEP_DEFCO=T:MG
    ! JCW: If present, use ONETEP's own defect correction code, rather than the
    ! JCW: defect correction code in DL_MG
    multigrid_use_onetep_defco = utils_devel_code(.false.,"MG","USE_ONETEP_DEFCO",pub_devel_code)
    if (pub_on_root.and.multigrid_use_onetep_defco) write(stdout,'(a)') &
       "MG Defect correction will be performed in ONETEP, if required"

    ! we need to find the start and end point in global grid coordinate for the local array
    ! x and y directions are trivial
    if ( mg_comm /= pub_null_comm) then
       call comms_rank(mg_proc_id, mg_comm)
       allocate( zlim(0:pub_total_num_procs-1), stat=ierr)
       call utils_alloc_check(myself,'zlim',ierr)
       call comms_allgather(zlim, mg%num_slabs12_pq, comm = mg_comm)

       idx_start = (/ 1,       1,       sum(zlim(0:mg_proc_id - 1)) + 1 /)
       idx_end   = (/ mg%pq1f, mg%pq2f, idx_start(3) + mg%num_slabs12_pq - 1 /)
       mgunit = utils_unit()
       deallocate(zlim,stat=ierr)
       call utils_dealloc_check(myself,'zlim',ierr)
       ! JCW: Set some local variables to be passed to DL_MG as arguments
       n  = [ mg%pq1f, mg%pq2f, mg%pq3f ]
       d  = [ SQRT(grid%da1%X**2+grid%da1%Y**2+grid%da1%Z**2), &
              SQRT(grid%da2%X**2+grid%da2%Y**2+grid%da2%Z**2), &
              SQRT(grid%da3%X**2+grid%da3%Y**2+grid%da3%Z**2) ]

    endif

#if defined (__INTEL_COMPILER) && defined (DEBUG)
    ! JCW: Workaround for possible compiler bug for Intel compiler in debug
    ! JCW: mode, with > 1 OpenMP thread. This occurs in DL_MG and is possible
    ! JCW: due to misreading of array range info in the mg struct.
    ! JCW: Since the issue only occurs with > 1 OpenMP thread, we temporarily
    ! JCW: change the number of threads while DL_MG is active:
!$  if (pub_on_root) then
!$     write(stdout,'(a)') &
!$          "MG WARNING: Temporarily reducing number of OpenMP threads to 1 to work around"//CRLF//&
!$          "            possible compiler bug when DL_MG is compiled with Intel compiler"//CRLF//&
!$          "            in debug mode."
!$  end if
!$  nthreads_save = omp_get_max_threads()
!$  call omp_set_num_threads(1)
    ! JCW: This may be the same possible compiler bug described in ONETEP
    ! JCW: CCPForge bug #1410.
#endif

    ! JCW: Identical initialization for DL_MG with ONETEP or DL_MG's defect
    ! JCW: correction code
    call dl_mg_init(n(1), n(2), n(3), d(1), d(2), d(3), bc, &
         idx_start, idx_end, mg_comm, mgunit, &
         trim(pub_rootname)//"_dl_mg_log.txt", &
         errors_return = pub_mg_continue_on_error)
#else
    call utils_feature_not_supported('DL_MG multigrid solver')
#endif

#if defined (__INTEL_COMPILER) && defined (DEBUG)
    ! JCW: Workaround for possible compiler bug for Intel compiler in debug
    ! JCW: mode, with > 1 OpenMP thread. This occurs in DL_MG and is possible
    ! JCW: due to misreading of array range info in the mg struct.
    ! JCW: Since the issue only occurs with > 1 OpenMP thread, we temporarily
    ! JCW: change the number of threads while DL_MG is active:
!$  if (pub_on_root) then
!$     write(stdout,'(a)') &
!$          "MG WARNING: Restoring number of OpenMP threads to previous value."
!$  end if
!$  call omp_set_num_threads(nthreads_save)
    ! JCW: This may be the same possible compiler bug described in ONETEP
    ! JCW: CCPForge bug #1410.
#endif

    ! jd: Allocate an array to hold phi_vac, required if smeared-ion correction
    !     to forces is needed.
    if(pub_forces_needed .and. pub_is_smeared_ion_rep) then
       allocate(phi_vac_for_forces(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check('multigrid_initialise','phi_vac_for_forces',ierr)
    end if

    multigrid_initialised = .true.

    call utils_trace_out(myself)

  end subroutine multigrid_initialise
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine multigrid_bc_for_dlmg(bc_dlmg,bc_is_periodic)
    !=========================================================================!
    ! Based on the values of the logical array bc_is_periodic, return an      !
    ! integer array bc_dlmg which contains the integer constants DL_MG uses   !
    ! to indicate the type of BCs to be used by the solver.                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  bc_dlmg        intent=out,      the BCs in DL_MG-compatible format     !
    !  bc_is_periodic intent=in ,      array of logicals indicating whether   !
    !                                  BCs are periodic in each direction     !
    !-------------------------------------------------------------------------!
    ! Written by James C. Womack, 05/2017                                     !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, dimension(3), intent(out) :: bc_dlmg
    logical, dimension(3), intent(in) :: bc_is_periodic

    ! Parameters
    character(len=*), parameter :: myself = "multigrid_bc_for_dlmg"

    ! Local variables
    integer :: ibc

    ! Derive DL_MG array from logical bc_is_periodic_loc array
#ifdef HAVE_DL_MG
    do ibc = 1, 3

       if (bc_is_periodic(ibc)) then
          bc_dlmg(ibc)        = DL_MG_BC_PERIODIC
       else
          bc_dlmg(ibc)        = DL_MG_BC_DIRICHLET
       end if

    end do
#else
    bc_dlmg(:) = -1
#endif

  end subroutine multigrid_bc_for_dlmg
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function multigrid_electrostatic_energy(phi,grid,eps)
    !=========================================================================!
    ! Calculates the electrostatic energy due to a potential phi in a dielec- !
    ! tric with permittivity eps, via (4) in [3], i.e.                        !
    ! E_es = 1/8pi \int eps(r) * (grad phi(r))^2 dr.                          !
    ! This formula is only valid when integrating over all space, in a finite !
    ! box this is only an approximation, rather accurate for neutral molecules!
    ! and quite crude for charged species.                                    !
    ! NB: This value is never used in the calculation, it is only displayed   !
    !     to the user, so that the discretization error can be assessed.      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi,               intent=in,  the potential                          !
    !   eps,               intent=in,  the dielectric permittivity            !
    !   ... both in the usual distributed representation.                     !
    ! Returns:                                                                !
    !   calculated electrostatic energy.                                      !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 13/12/2010.                                  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP, PI
    use finite_differences, only: finite_difference_mod_grad_sqd
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_finite_difference_order
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)           :: phi(mg%ld1f,mg%ld2f,mg%max_slabs12)
    type(GRID_INFO), intent(in)         :: grid
    real(kind=DP), intent(in), optional :: eps(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! jd: Local variables
    real(kind=DP), allocatable :: grad_phi_sqd(:,:,:) ! jd: (\nabla \phi)^2
    real(kind=DP) :: integral                         ! jd: see below
    integer :: ierr                                   ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('multigrid_electrostatic_energy')

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    allocate(grad_phi_sqd(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check('multigrid_electrostatic_energy','grad_phi_sqd',ierr)

    ! jd: Compute (\nabla \phi)^2
    call finite_difference_mod_grad_sqd(grad_phi_sqd, phi, &
         pub_finite_difference_order, mg, grid)

    ! jd: Compute \int (eps * (\nabla \phi)^2) dr
    if(present(eps)) then
       integral = integrals_product_on_grid(grid, &
            grad_phi_sqd,eps,mg%pq1f, mg%pq2f, mg%num_slabs12_pq)
    else
       integral = integrals_product_on_grid(grid, &
            grad_phi_sqd, m1 = mg%pq1f, m2 = mg%pq2f, m3 = mg%num_slabs12_pq)
    end if

    ! jd: Return 1/8pi times the above integral, cf (4) in [3].
    multigrid_electrostatic_energy = 1.0_DP/(8.0_DP * pi) * integral

    deallocate(grad_phi_sqd,stat=ierr)
    call utils_dealloc_check('multigrid_electrostatic_energy','grad_phi_sqd', &
         ierr)

    call utils_trace_out('multigrid_electrostatic_energy')

  end function multigrid_electrostatic_energy
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function multigrid_finite_box_correction(phi,grid,eps)
    !=========================================================================!
    ! Calculates the correction to the electrostatic energy returned by       !
    ! multigrid_electrostatic_energy() resulting from the finiteness of the   !
    ! box, i.e. the term -1/8pi * surfint eps(r) phi(r) grad_phi(r) dS.       !
    ! NB: This value is never used in the calculation, it is only displayed   !
    !     to the user, so that the discretization error can be assessed.      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi,               intent=in,  the potential                          !
    !   eps,               intent=in,  the dielectric permittivity            !
    !   ... both in the usual distributed representation.                     !
    ! Returns:                                                                !
    !   calculated correction.                                                !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 13/12/2010.                                  !
    ! Improved by Jacek Dziedzic, 09/06/2013 to handle 0 slabs per core.      !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id
    use constants, only: DP, PI
    use finite_differences, only: finite_difference_gradient
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_finite_difference_order
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in)         :: grid
    real(kind=DP), intent(in)           :: phi(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional :: eps(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! jd: Local variables
    real(kind=DP), allocatable :: grad_phi(:,:,:,:)       ! jd: grad phi
    real(kind=DP), allocatable :: phi_gradphi_ds(:,:,:)   ! jd: phi*grad_phi*dS
    real(kind=DP) :: integral                             ! jd: see below
    real(kind=DP) :: gradphi_ds                           ! jd: grad_phi * dS
    real(kind=DP) :: surf_weight_x, surf_weight_y, surf_weight_z ! jd: weights
    integer :: i1, i2, islab12, i3_global                 ! jd: indices
    integer :: ierr                                       ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('multigrid_finite_box_correction')

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    ! jd: Allocate work arrays
    allocate(grad_phi(mg%ld1f,mg%ld2f,mg%max_slabs12,3),stat=ierr)
    call utils_alloc_check('multigrid_finite_box_correction','grad_phi',ierr)
    allocate(phi_gradphi_ds(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check('multigrid_finite_box_correction', &
         'phi_gradphi_ds',ierr)

    ! jd: Compute grad \phi
    call finite_difference_gradient(grad_phi, phi, pub_finite_difference_order,&
         mg, grid)

    ! jd: Weights for dS -- the grid is orthorhombic, but not necessarily cubic
    surf_weight_x = grid%da2%y*grid%da3%z
    surf_weight_y = grid%da1%x*grid%da3%z
    surf_weight_z = grid%da1%x*grid%da2%y

    ! jd: Compute phi * grad_phi * dS on the surfaces of the box and
    !     a lot of zeros inside the box, will allow easy integration
    !     with the volume integrator.
    do islab12=1, mg%num_slabs12_pq
       i3_global = islab12 + grid%first_slab12(pub_my_proc_id)-1
       do i2=1, mg%pq2f
          do i1=1, mg%pq1f
             gradphi_ds = 0D0
             if(i1 == 1) gradphi_ds = &
                  gradphi_ds - grad_phi(i1,i2,islab12,1) * surf_weight_x
             if(i1 == mg%pq1f) gradphi_ds = &
                  gradphi_ds + grad_phi(i1,i2,islab12,1) * surf_weight_x
             if(i2 == 1) gradphi_ds = &
                  gradphi_ds - grad_phi(i1,i2,islab12,2) * surf_weight_y
             if(i2 == mg%pq2f) gradphi_ds = &
                  gradphi_ds + grad_phi(i1,i2,islab12,2) * surf_weight_y
             if(i3_global == 1) gradphi_ds = &
                  gradphi_ds - grad_phi(i1,i2,islab12,3) * surf_weight_z
             if(i3_global == mg%pq3f) gradphi_ds = &
                  gradphi_ds + grad_phi(i1,i2,islab12,3) * surf_weight_z

             phi_gradphi_ds(i1,i2,islab12) = phi(i1,i2,islab12) * gradphi_ds
          end do
       end do
    end do

    ! jd: Compute \int (eps * gphi * (grad \phi) dS)
    if(present(eps)) then
       integral = integrals_product_on_grid(grid, &
            phi_gradphi_ds,eps,mg%pq1f, mg%pq2f, mg%num_slabs12_pq)
    else
       integral = integrals_product_on_grid(grid, &
            phi_gradphi_ds, m1 = mg%pq1f, m2 = mg%pq2f, m3 = mg%num_slabs12_pq)
    end if

    ! jd: Take -1/8pi into account, also remove the dV weight that
    !     integrals_product_on_grid() adds, dS weights have already been used.
    multigrid_finite_box_correction = &
         -1.0_DP/(8.0_DP * pi) * integral / grid%weight

    ! jd: Clean up
    deallocate(phi_gradphi_ds,stat=ierr)
    call utils_dealloc_check('multigrid_finite_box_correction', &
         'phi_gradphi_ds', ierr)
    deallocate(grad_phi,stat=ierr)
    call utils_dealloc_check('multigrid_finite_box_correction','grad_phi', &
         ierr)

    call utils_trace_out('multigrid_finite_box_correction')

  end function multigrid_finite_box_correction
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine multigrid_prepare_bound_cond(bound, rho_in, grid, cell, uniform_eps, &
       bc_coarseness, bc_surface_coarseness)
    !=========================================================================!
    ! Calculates the Coulombic boundary conditions for a potential due to a   !
    ! density rho, on the surface points of the cell. Points in the bulk are  !
    ! zeroed, even though the MG solver never uses them. Coarse-graining of   !
    ! charge is used to keep the O(L^5) scaling in check.                     !
    ! In solvent an approximation is used in which it is assumed that the     !
    ! whole cell is filled with uniform dielectric, with a permittivity of    !
    ! uniform_eps. This will become inaccurate for solvents with low permitti-!
    ! vities. If uniform_eps is omitted, 1.0 is assumed.                      !
    !                                                                         !
    ! Also note that when the Boltzmann term is included (is_pbe /= NONE),    !
    ! Debye-length screening is employed in the BCs (unless is_pbe_bc_        !
    ! _debye_screening F). This is exact in the linearised formulation, and   !
    ! an approximation in the full formulation.                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   bound,       intent=out,  calculated Coulombic potential on the surf. !
    !   rho_in,      intent=in,   charge density generating the BCs           !
    !   ... both in the usual distributed representation.                     !
    !   uniform_eps, intent=in,   (optional) the uniform diel. permittivity   !
    !   bc_coarseness, intent=in, the coarseness of the BCs.                  !
    !   bc_surface_coarseness, intent=in, the surface coarseness of the BCs.  !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 13/12/2010.                                  !
    ! Revised by Jacek Dziedzic, 15/04/2011 to allow blocks spanning core     !
    !                                       boundaries.                       !
    ! Generalised by Jacek Dziedzic, 12/10/2014 to allow Debye screening.     !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, &
         cell_grid_extract_box, cell_grid_deposit_box
    use comms, only: pub_my_proc_id, pub_total_num_procs, &
         comms_bcast, comms_barrier, comms_reduce, pub_on_root, pub_root_proc_id
    use constants, only: VERBOSE, stdout
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(*), OPERATOR(+)
    use is_solvation_boltzmann, only: implicit_solvent_boltzmann_debye_length
    use rundat, only: pub_is_bc_threshold, pub_output_detail, &
         pub_is_pbe_bc_debye_screening, pub_debug_on_root, pub_charge, &
         pub_is_bc_allow_frac_charge
!$  use rundat, only: pub_threads_max
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_trace_in, utils_trace_out, utils_assert, utils_abort
    use visual, only: visual_scalarfield

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(out)  :: bound(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)   :: rho_in(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional :: uniform_eps
    integer, intent(in)         :: bc_coarseness
    integer, intent(in)         :: bc_surface_coarseness

    ! jd: Internal variables
    real(kind=DP) :: inv_uniform_eps   ! 1/uniform_eps
    integer :: i1, i2, i3, islab12     ! jd: Indices
    integer :: j_packed, n_packed, n_packed_plus, n_packed_minus, n_packed_max

    ! jd: Packed coarse-grained charges
    real(kind=DP), allocatable, dimension(:)      :: packed_q
    ! jd: ... and their locations
    real(kind=DP), allocatable, dimension(:,:)    :: packed_loc
    ! jd: Stores a one-block-thick layer of rho_in, extracted from other procs
    real(kind=DP), allocatable, dimension(:,:,:)  :: rho
    ! jd: Buffer for cell_grid_extract_box, cell_grid_deposit_box
    real(kind=DP), allocatable, dimension(:,:,:)  :: buffer
    ! jd: Value of charge per gridpoint below which all charges are ignored
    real(kind=DP) :: q_threshold
    real(kind=DP) :: q_i
    real(kind=DP) :: q_factor
    real(kind=DP) :: v_here
    type(POINT) :: r_i, r_j, r_0 ! jd: Various vectors
    type(POINT) :: r_of_block_plus, r_of_block_minus ! jd: "Centre of charge"
    integer :: n_packed_first, n_packed_last ! jd: CG'ed charges on this proc
    integer :: my_portion       ! jd: # of CG'ed charges on this proc
    integer :: typical_portion  ! jd: # of CG'ed charges on every proc but last
    integer :: tail             ! jd: Difference between the two, on last proc
    real(kind=DP) :: d          ! jd: Distance
    real(kind=DP) :: q_in_block_plus  ! jd: Positive charge in a CG'ed block
    real(kind=DP) :: q_in_block_minus ! jd: Negative charge in a CG'ed block
    real(kind=DP) :: q_sum      ! jd: Charge a in CG'ed block
    integer :: n_in_block_plus, n_in_block_minus ! jd: Relevant charge points in blk
    integer :: bsize, bheight   ! jd: Desired and actual block height
    integer :: nb1, nb2, nb3    ! jd: Number of blocks along x, y, z
    integer :: b1, b2, b3       ! jd: Block indices
    integer :: i1min, i1max, i2min, i2max, i3min, i3max ! jd: Bounds
    real(kind=DP) :: q ! jd: Current coarse-grained charge point
    integer :: n_points_plus, n_points_minus
    integer :: n_empty_blocks

    integer :: face    ! jd: Face of the cell under consideration (1-6).
    ! jd: Calculated quantity
    real(kind=DP), allocatable, target, dimension(:,:) :: face_pot
    ! jme: pointer used to remap bounds
    real(kind=DP), pointer, dimension(:,:,:) :: face_pot_pointer
    integer :: i1_pq, i1_ld, i2_pq, i2_ld, i3_pq ! jd: Bounds of generic indices
    type(POINT) :: i1_da, i2_da, i3_da ! Grid widths along each generic index
    ! jd: These refer to Cartesian (non-generic) representation
    integer :: width_x, width_y, width_z, start_x, start_y, start_z, ld_x, ld_y

    ! jd: The following refer to surface CG'ing (stage 3).
    integer :: prev1, prev2, next1, next2
    integer :: offset1, offset2
    real(kind=DP) :: t, u
    real(kind=DP) :: width1, width2
    real(kind=DP) :: val_prev_t_prev_u, val_prev_t_next_u, val_next_t_prev_u, &
         val_next_t_next_u

    ! jd: Varia
    integer :: ierr ! jd: Error flag
    character(len=*), parameter :: myself = 'multigrid_prepare_bound_cond'
    logical, save :: ringing_warning_issued = .false.

    ! ebl: handling fractional charge
    real(kind=DP) :: frac_charge

    !------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself//'.'
    call timer_clock(myself,1)
    call utils_trace_in(myself)

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    !     We don't test against 'multigrid_initialised', since we may have been
    !     called from the DMA module, which just wants the potential.
    call utils_assert(mg%ld1f >=0 .and. mg%ld2f >=0 .and. mg%max_slabs12 >=0, &
         trim(myself)//': The mg variable seems not to have been initialised.')

    if(multigrid_debug) then
       call visual_scalarfield(rho_in,grid,cell, &
             'Charge density','_rhoin_in_prepare_bound_cond')
    end if

    if(pub_is_pbe_bc_debye_screening) then
       call utils_assert(implicit_solvent_boltzmann_debye_length > 0.0_DP, &
            myself//': Illegal Debye length:', &
            implicit_solvent_boltzmann_debye_length)
    end if

    ! jme: initialize pointer status
    nullify(face_pot_pointer)

    bound = 0.0_DP

    ! jd: bc_coarseness of 0 means zero BC's
    if(bc_coarseness == 0) then
      call utils_trace_out(myself)
      call timer_clock(myself,2)
      return
    end if

    if(present(uniform_eps)) then
       inv_uniform_eps = 1.0_DP/uniform_eps
    else
       inv_uniform_eps = 1.0_DP
    end if

    ! jd: Determine the number of blocks along each axis, globally
    bsize = bc_coarseness
    nb3 = mg%pq3f/bsize
    nb2 = mg%pq2f/bsize
    nb1 = mg%pq1f/bsize

    ! jd: Count partial blocks as well
    if(mod(mg%pq3f,bsize) /= 0) nb3 = nb3 + 1
    if(mod(mg%pq2f,bsize) /= 0) nb2 = nb2 + 1
    if(mod(mg%pq1f,bsize) /= 0) nb1 = nb1 + 1

    ! jd: Allocate temporaries
    n_packed_max = 2*nb1*nb2*nb3
    ! JCW: Factor of 2 in n_packed_max because for each CG block we can have
    ! JCW: one positive and one negative charge stored in the packed_q and
    ! JCW: packed_loc arrays.
    allocate(packed_q(n_packed_max),stat=ierr)
    call utils_alloc_check(myself,'packed_q',ierr)
    allocate(packed_loc(3,n_packed_max),stat=ierr)
    call utils_alloc_check(myself,'packed_loc',ierr)
    allocate(rho(mg%ld1f,mg%ld2f,bsize),stat=ierr)
    call utils_alloc_check(myself,'rho',ierr)
    allocate(buffer(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'buffer',ierr)

    ! ---------------------------------------------------------------------
    ! STAGE 1
    ! jd: Go over rho_in in layers, one-block thick.
    !     Get the layer to root proc, regardless of who owns the slabs.
    !     On the root proc, go over all blocks in the layer, examine all
    !     entries in the block that are above a tiny threshold.
    !     Each entry in packed_q represents the total charge in the block.
    !     Each entry in packed_loc is the location of the centre-of-charge.
    !     A value of -1 in packed_loc signifies end of the array,
    !     beyond which the entries in packed_q are undefined
    !     This stage is very cheap (compared to stage 2), so we can get away
    !     with doing this on the root proc.
    ! ---------------------------------------------------------------------

    packed_loc = 1D33     ! jd: Garbage
    ! jd: All grid points with |charge| lower than that will be ignored.
    q_threshold = pub_is_bc_threshold
    n_packed = 0
    n_packed_plus = 0
    n_packed_minus = 0
    n_points_plus = 0
    n_points_minus = 0
    n_empty_blocks = 0
    q_sum = 0.0_DP

    ! jd: Go over layers, each one-block thick
    do b3 = 1, nb3

       ! jd: Height of the block in this layer
       !     (=bsize except in last global layer)
       if(b3 < nb3) then
          bheight = bsize
       else
          bheight = mod(mg%pq3f,bsize)
       end if
       ! jd: Extract bheight slabs, mg%pq1f x mg%pq2f wide into rho,
       !     regardless of who owns them
       rho = 0D0
       call cell_grid_extract_box(rho, buffer, rho_in, grid, &
            mg%ld1f, mg%ld2f, bheight, mg%ld1f, mg%ld2f, &
            1, 1, (b3-1)*bsize+1, &
            pub_on_root, .false.)

       ! jd: Go over blocks in the current layer, packing them
       if(pub_on_root) then

         i3min = (b3-1)*bsize+1
         i3max = b3*bsize
         if(i3max > mg%pq3f) i3max = mg%pq3f

         do b2=1, nb2
            i2min = (b2-1)*bsize+1
            i2max = b2*bsize
            if(i2max > mg%pq2f) i2max = mg%pq2f

            do b1=1, nb1
               i1min = (b1-1)*bsize+1
               i1max = b1*bsize
               if(i1max > mg%pq1f) i1max = mg%pq1f

               ! jd: Centre of the block
               r_0 = 0.5_DP * ( &
                    real((i1min-1),kind=DP) * grid%da1 + &
                    real((i2min-1),kind=DP) * grid%da2 + &
                    real((i3min-1),kind=DP) * grid%da3 + &
                    real((i1max-1),kind=DP) * grid%da1 + &
                    real((i2max-1),kind=DP) * grid%da2 + &
                    real((i3max-1),kind=DP) * grid%da3 )

               n_in_block_plus = 0
               n_in_block_minus = 0
               q_in_block_plus = 0.0_DP
               q_in_block_minus = 0.0_DP
               r_of_block_plus%X = 0.0_DP
               r_of_block_plus%Y = 0.0_DP
               r_of_block_plus%Z = 0.0_DP
               r_of_block_minus%X = 0.0_DP
               r_of_block_minus%Y = 0.0_DP
               r_of_block_minus%Z = 0.0_DP

               ! jd: Go over points in current block
               do islab12=1, bheight
                  i3 = i3min + islab12 - 1
                  do i2=i2min, i2max
                     do i1=i1min, i1max

!                        call utils_assert(i1<=mg%pq1f .and. i1>0, &
!                             'Bound check [1]')
!                        call utils_assert(i2<=mg%pq2f .and. i2>0, &
!                             'Bound check [2]')
!                        call utils_assert(i3<=mg%pq3f .and. i3>0, &
!                             'Bound check [3]')
                        q_i = rho(i1,i2,islab12) * grid%weight

                        if(q_i > q_threshold) then
                           r_i = &
                                real((i1-1),kind=DP) * grid%da1 + &
                                real((i2-1),kind=DP) * grid%da2 + &
                                real((i3-1),kind=DP) * grid%da3

                           q_in_block_plus = q_in_block_plus + q_i
                           r_of_block_plus = r_of_block_plus + q_i* (r_i-r_0)
                           n_in_block_plus = n_in_block_plus + 1
                        else if(q_i < -q_threshold) then
                           r_i = &
                                real((i1-1),kind=DP) * grid%da1 + &
                                real((i2-1),kind=DP) * grid%da2 + &
                                real((i3-1),kind=DP) * grid%da3

                           q_in_block_minus = q_in_block_minus + q_i
                           r_of_block_minus = r_of_block_minus + q_i* (r_i-r_0)
                           n_in_block_minus = n_in_block_minus + 1
                        end if

                     end do  ! }
                  end do     ! } over gridpoints in block
               end do        ! }

               ! jd: Store packed block (positive charges)
               if(n_in_block_plus /= 0) then
                  n_packed = n_packed + 1
                  n_packed_plus = n_packed_plus + 1
                  packed_q(n_packed) = q_in_block_plus
                  q_sum = q_sum + q_in_block_plus
                  r_of_block_plus = 1.0_DP/q_in_block_plus * r_of_block_plus + r_0
                  packed_loc(1,n_packed) = r_of_block_plus%X
                  packed_loc(2,n_packed) = r_of_block_plus%Y
                  packed_loc(3,n_packed) = r_of_block_plus%Z
                  n_points_plus = n_points_plus + n_in_block_plus
               end if

               ! jd: Store packed block (negative charges)
               if(n_in_block_minus /= 0) then
                  n_packed = n_packed + 1
                  n_packed_minus = n_packed_minus + 1
                  packed_q(n_packed) = q_in_block_minus
                  q_sum = q_sum + q_in_block_minus
                  r_of_block_minus = 1.0_DP/q_in_block_minus * r_of_block_minus + r_0
                  packed_loc(1,n_packed) = r_of_block_minus%X
                  packed_loc(2,n_packed) = r_of_block_minus%Y
                  packed_loc(3,n_packed) = r_of_block_minus%Z
                  n_points_minus = n_points_minus + n_in_block_minus
               end if

               if(n_in_block_minus == 0 .and. n_in_block_plus == 0) then
                  n_empty_blocks = n_empty_blocks + 1
               end if

            end do ! }
         end do    ! } over blocks
      end if ! if on root
    end do  ! over layers

    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
        write(stdout,'(a,i9,a,i11,a)') 'MG BCs: +ve charge blocks (gridpoints): ', &
             n_packed_plus, ' (', n_points_plus, ').'
        write(stdout,'(a,i9,a,i11,a)') 'MG BCs: -ve charge blocks (gridpoints): ', &
             n_packed_minus, ' (', n_points_minus, ').'
        write(stdout,'(a,i9,a,i11,a)') 'MG BCs:     Elided blocks (gridpoints): ', &
             n_empty_blocks, ' (', &
             mg%pq1f*mg%pq2f*mg%pq3f-n_points_plus-n_points_minus, ').'
        write(stdout,'(a,f8.4,a)') 'MG BCs: Computational effort saved: ', &
             (mg%pq1f*mg%pq2f*mg%pq3f-n_packed_plus-n_packed_minus) &
             /real(mg%pq1f*mg%pq2f*mg%pq3f,kind=DP)*100.0_DP,'%.'
        write(stdout,'(a,f20.7,a)') 'MG BCs: Charge after coarse-graining: ', &
             q_sum,'.'
    end if

    ! ebl: account for cases where we have a non-integer number of electrons
    frac_charge = abs(pub_charge - real(nint(pub_charge),kind=DP))

    if(pub_on_root .and. .not. pub_is_bc_allow_frac_charge) then
       call utils_assert(abs(abs(q_sum - real(nint(q_sum),kind=DP)) - &
            frac_charge) < 0.05_DP, 'The coarse-grained charge is suspiciously &
            &far from its expected value. First check that your charge is corre&
            &ctly localized within the FD grid. Allow up to 10 bohr of extra va&
            &cuum to account for the ringing, especially if fine_grid_scale is &
            &different from 2.0. If you get this error during a properties calc&
            &ulation, make sure that the calculation is well-converged first. I&
            &f none of the above apply, reduce is_bc_threshold (to, say, 1E-10)&
            & to make charge coarse-graining more accurate. If this doesn''t he&
            &lp, this error might indicate that you need to decrease ''is_bc_co&
            &arseness'' to improve accuracy. The obtained coarse-grained charge&
            & was', q_sum)

       call utils_assert(n_packed /= 0, 'Charge density is zero (or extremely &
            &low) everywhere! Or perhaps your is_bc_threshold is insanely big?')

       ! jd: Because we killed off all charge below the threshold, we've left
       !     out the ringing. Better adjust the CGed charge by multiplying it by
       !     q_desired / q_cged (a factor in the order of 1.00001) so that the
       !     the total charge agrees. This must happen on root only.
       q_factor = real(nint(q_sum),kind=DP) / q_sum
       if(q_factor /= 0.0_DP) then ! jd: But not when the molecule is neutral
          packed_q(1:n_packed) = packed_q(1:n_packed) * q_factor
       end if
    end if

    ! jd: These are no longer needed
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check(myself,'buffer',ierr)
    deallocate(rho,stat=ierr)
    call utils_dealloc_check(myself,'rho',ierr)

    ! jd: Broadcast the packed arrays to everyone
    call comms_bcast(pub_root_proc_id,packed_q)
    call comms_bcast(pub_root_proc_id,packed_loc)
    call comms_bcast(pub_root_proc_id,n_packed)

    ! ---------------------------------------------------------------------
    ! STAGE 2
    ! jd: Go over packed array and evaluate potential on coarse surface points.
    !     Then fill in the missing points by bilinear interpolation.
    !     Each processor deals with *all* the points on the boundary, but only
    !     a subset of the packed point charges. It doesn't make sense to use the
    !     usual slab distribution, because then the first and last procs have
    !     much more work -- they deal with top and bottom surfaces. Also the
    !     surface coarse-graining would become more involved (as surface blocks
    !     would cross processor boundaries).
    ! ---------------------------------------------------------------------

    ! jd: Determine how big our portion of the packed array is
    !     When n_packed is not divisible by n_procs, give the remainder to
    !     the last proc -- we have n_packed in the order of 10^3-10^6, so
    !     the extra tail is negligible.
    my_portion = n_packed / pub_total_num_procs
    typical_portion = my_portion
    tail = mod(n_packed,pub_total_num_procs)
    if(pub_my_proc_id == pub_total_num_procs-1) my_portion = my_portion + tail
    n_packed_first = pub_my_proc_id * typical_portion + 1
    n_packed_last = n_packed_first + my_portion - 1

    bsize = bc_surface_coarseness

    ! --------------------------------------------------
    ! --- Go over the 6 faces, fill them and deposit ---
    ! --------------------------------------------------
    ! jd: The face is covered by indices i1, i2.

    allocate(buffer(grid%ld1,grid%ld2, &
         grid%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'buffer',ierr)

    do face = 1, 6
       ! jd: Silence compiler warning
       i1_ld = 0; i1_pq = 0; i1_ld = 0; i2_ld = 0; i2_pq = 0; i3_pq = 0
       ! jd: Top or bottom face (XY): i1=x, i2=y, i3=z
       if(face == 1 .or. face == 2) then
          i1_ld = grid%ld1
          i1_pq = mg%pq1f
          i1_da = grid%da1
          i2_ld = grid%ld2
          i2_pq = mg%pq2f
          i2_da = grid%da2
          i3_pq = mg%pq3f
          i3_da = grid%da3
          ld_x = i1_ld
          ld_y = i2_ld
          width_x = i1_pq
          width_y = i2_pq
          width_z = 1
          start_x = 1
          start_y = 1
          if(face == 1 ) then
             start_z = 1           ! top
          else
             start_z = i3_pq       ! bottom
          end if
       end if
       ! jd: Left or right face (YZ): i1=y, i2=z, i3=x
       if(face == 3 .or. face == 4) then
          i1_ld = grid%ld2
          i1_pq = mg%pq2f
          i1_da = grid%da2
          i2_ld = grid%ld3
          i2_pq = mg%pq3f
          i2_da = grid%da3
          i3_pq = mg%pq1f
          i3_da = grid%da1
          ld_x = 1
          ld_y = i1_ld
          width_x = 1
          width_y = i1_pq
          width_z = i2_pq
          start_y = 1
          start_z = 1
          if(face == 3 ) then
             start_x = 1           ! left
          else
             start_x = i3_pq       ! right
          end if
       end if
       ! jd: Front or back (XZ): i1=x, i2=z, i3=y
       if(face == 5 .or. face == 6) then
          i1_ld = grid%ld1
          i1_pq = mg%pq1f
          i1_da = grid%da1
          i2_ld = grid%ld3
          i2_pq = mg%pq3f
          i2_da = grid%da3
          i3_pq = mg%pq2f
          i3_da = grid%da2
          ld_x = i1_ld
          ld_y = 1
          width_x = i1_pq
          width_y = 1
          width_z = i2_pq
          start_x = 1
          start_z = 1
          if(face == 5 ) then
             start_y = 1           ! front
          else
             start_y = i3_pq       ! back
          end if
       end if

       ! jd: Allocate and zero storage for this cell face
       allocate(face_pot(i1_ld,i2_ld),stat=ierr)
       call utils_alloc_check(myself,'face_pot',ierr)
       face_pot = 0D0

       ! jd: Low-index-value or high-index-value faces?
       if(face == 1 .or. face == 3 .or. face == 5) then
          i3 = 1     ! jd: Face location: top, left or front
       else
          i3 = i3_pq ! jd: Face location: bottom, right or back
       end if

       ! jd: Loop over points on current face.
       !     Everything's been cleared earlier, so no need to worry about
       !     padding margins.
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(i1,i2,r_i,r_j,v_here,j_packed,q,d) &
!$OMP SHARED(i3,i1_pq,i2_pq,bsize,i1_da,i2_da,i3_da,n_packed_first, &
!$OMP      n_packed_last,packed_loc,packed_q,ringing_warning_issued, &
!$OMP      inv_uniform_eps,face_pot,pub_threads_max,&
!$OMP      implicit_solvent_boltzmann_debye_length, &
!$OMP      pub_is_pbe_bc_debye_screening,stdout)
       do i2=1, i2_pq
          do i1=1, i1_pq
             ! jd: If coarse-graining surfaces, skip inner points
             if(bsize /=1 .and. &
                  ((mod(i1,bsize) /= 1 .and. i1 /= i1_pq) .or. &
                   (mod(i2,bsize) /= 1 .and. i2 /= i2_pq))) cycle

             r_i = &
                  real((i1-1),kind=DP) * i1_da + &
                  real((i2-1),kind=DP) * i2_da + &
                  real((i3-1),kind=DP) * i3_da

             v_here = 0.0_DP

             ! jd: Go over all my points in the packed array
             do j_packed = n_packed_first, n_packed_last
                r_j%X = packed_loc(1,j_packed)
                r_j%Y = packed_loc(2,j_packed)
                r_j%Z = packed_loc(3,j_packed)

                q = packed_q(j_packed)
                d = magnitude(r_i-r_j)

                ! jd: Packed charge points very close to the face indicate
                !     something's seriously wrong.
                if(abs(d) < 1D-2) then
                   if(.not. ringing_warning_issued) then
                      write(stdout,'(a)') 'WARNING: Non-zero charge at, or very&
                           & close to, one of the simulation cell faces.'
                      write(stdout,'(a,f0.5,a,f0.5,a,f0.5,a)') '     &
                           &    Charge at:     ',r_j%X,',',r_j%Y,',',r_j%Z,'.'
                      write(stdout,'(a,f0.5,a,f0.5,a,f0.5,a)') '     &
                           &    Point on face: ',r_i%X,',',r_i%Y,',',r_i%Z,'.'
                      write(stdout,'(a,e21.5)') &
                           'WARNING: Charge:                 ',q
                      write(stdout,'(a,e21.5)') &
                           'WARNING: Distance from box face: ',d
                   end if
                   if(abs(q) > 1D-7) then
                      call utils_abort('Non-zero charge density&
                           & detected on or very close to cell boundaries.&
                           & Check your cell bounds and charge localization.')
                   else
                      ! jd: Unless it's just the ringing, then ignore it
                      if(.not. ringing_warning_issued) then
                         write(stdout,'(a)') 'WARNING: Density ringing detected on the&
                              & face of the simulation cell.'
                         write(stdout,'(a)') '         This may lead to&
                              & inaccuracies in the charge coarse-graining&
                              & procedure.'
                         write(stdout,'(a)') '         Consider making your &
                              &simulation cell slightly larger.'
                         ringing_warning_issued = .true.
                      end if
                      ! jd: Avoid the singularity
                      q = 0D0
                      d = 1D0
                   end if
                end if
                if(.not. pub_is_pbe_bc_debye_screening) then
                   v_here = v_here + q / d
                else
                   v_here = v_here + q / d * &
                        exp(-d/implicit_solvent_boltzmann_debye_length)
                end if
             end do
             v_here = v_here * inv_uniform_eps
             face_pot(i1,i2) = v_here

          end do  ! } over points on this face
       end do     ! }
!$OMP END PARALLEL DO

       ! ---------------------------------------------------------------------
       ! STAGE 3
       ! jd: Now fill in the missing points on the surfaces by bilinear
       !     interpolation
       ! ---------------------------------------------------------------------

       if(bsize /= 1) then ! jd: No need to do this if all points filled already
          ! jd: Loop over points on current face.
          do i2=1, i2_pq
             do i1=1, i1_pq

                ! jd: Ignore points that *are* on the coarsened surface
                if(  (i1==i1_pq .or. mod(i1,bsize) == 1) .and. &
                     (i2==i2_pq .or. mod(i2,bsize) == 1)) then
                     if(face_pot(i1,i2) == 0D0) then
                        call utils_abort("Internal error [1] in "//myself)
                     end if
                     cycle
                end if

                r_i = &
                     real((i1-1),kind=DP) * i1_da + &
                     real((i2-1),kind=DP) * i2_da + &
                     real((i3-1),kind=DP) * i3_da

                call utils_assert(face_pot(i1,i2) == 0D0, &
                     'Internal error [2] in '//myself)

                ! jd: Our offset within a CG surface block (0..bsize-1)
                offset1 = mod(i1-1,bsize)
                offset2 = mod(i2-1,bsize)

                ! jd: Determine what the coarse-grained neighbours are
                prev1 = i1 - offset1
                prev2 = i2 - offset2
                next1 = i1 + (bsize-offset1)
                next2 = i2 + (bsize-offset2)
                if(prev1 < 1) prev1 = 1
                if(prev2 < 1) prev2 = 1
                if(next1 > i1_pq) next1 = i1_pq
                if(next2 > i2_pq) next2 = i2_pq

                ! jd: Determine the size of the bilinear interpolation block
                !     This is usually bsize, except close to pq boundaries
                width1=real(next1-prev1,kind=DP)
                width2=real(next2-prev2,kind=DP)

                ! jd: Take care of the corner case where prev==next.
                !     This can happen e.g. when i1==pq1, and prev is i1,
                !     but simultaneously next is i1, because it cannot lie
                !     outside ot the boundary. We then want 't' or 'u' to
                !     be zero, not NaN.
                if(width1==0.0_DP) width1=1.0_DP
                if(width2==0.0_DP) width2=1.0_DP

                ! jd: Pick right neighbours
                val_prev_t_prev_u = face_pot(prev1,prev2)
                val_prev_t_next_u = face_pot(prev1,next2)
                val_next_t_prev_u = face_pot(next1,prev2)
                val_next_t_next_u = face_pot(next1,next2)
                t=real(offset1,kind=DP)/width1
                u=real(offset2,kind=DP)/width2

                ! jd: Apply bilinear interpolation
                v_here = &
                     (1.0_DP-t)*(1.0_DP-u)*val_prev_t_prev_u + &
                     (1.0_DP-t)*u*val_prev_t_next_u + &
                     t*(1.0_DP-u)*val_next_t_prev_u + &
                     t*u*val_next_t_next_u

                face_pot(i1,i2) = v_here

             end do  ! } over points
          end do     ! } on the face
       end if ! jd: If surface coarse-graining

       ! jd: Careful not to double-count edges
       face_pot(1,:) = face_pot(1,:)/2.0_DP
       face_pot(:,1) = face_pot(:,1)/2.0_DP
       face_pot(i1_pq,:) = face_pot(i1_pq,:)/2.0_DP
       face_pot(:,i2_pq) = face_pot(:,i2_pq)/2.0_DP

       ! jd: Ditto for corners
       face_pot(1,1) = face_pot(1,1) * 4.0_DP/3.0_DP
       face_pot(i1_pq,1) = face_pot(i1_pq,1) * 4.0_DP/3.0_DP
       face_pot(1,i2_pq) = face_pot(1,i2_pq) * 4.0_DP/3.0_DP
       face_pot(i1_pq,i2_pq) = face_pot(i1_pq,i2_pq) * 4.0_DP/3.0_DP

       ! jd: Reduce contributions to this face over all procs to root
       call comms_reduce('SUM',face_pot,root=pub_root_proc_id)

       ! jme: use pointer to remap bounds (CCPForge issue #1742)
       face_pot_pointer(1:ld_x, 1:ld_y, 1:width_z) => face_pot(:,:)

       ! jd: Deposit this from root to the slab distrib'n on respective procs
       call cell_grid_deposit_box(bound, face_pot_pointer, buffer, grid, &
            width_x, width_y, width_z, ld_x, ld_y, &
            start_x, start_y, start_z, pub_on_root, .false.)

       nullify(face_pot_pointer) ! jme: prevent dangling pointer
       deallocate(face_pot,stat=ierr)
       call utils_dealloc_check(myself,'face_pot',ierr)

    end do ! over faces

    deallocate(buffer,stat=ierr)
    call utils_dealloc_check(myself,'buffer',ierr)

    ! jd: Sync to get accurate timings
    call comms_barrier

    ! jd: Clean up
    deallocate(packed_loc,stat=ierr)
    call utils_dealloc_check(myself,'packed_loc',ierr)
    deallocate(packed_q,stat=ierr)
    call utils_dealloc_check(myself,'packed_q',ierr)

    call utils_trace_out(myself)
    call timer_clock(myself,2)
    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving '//myself//'.'

  end subroutine multigrid_prepare_bound_cond
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function multigrid_calculate_hartree(phi, rho, grid, cell, &
       uniform_eps, eps_full, eps_half, solvation_terms, elements, no_dump)
    !=========================================================================!
    ! This is a back-end subroutine for the calculation of the Hartree energy !
    ! and Hartree potential using a multigrid solver. There are three,        !
    ! mutually-exclusive cases when this is in order:                         !
    !   a) A calculation without smeared ions, in vacuum.                     !
    !   b) A calculation, with smeared ions, in vacuum.                       !
    !   c) A calculation, with smeared ions, in solvent.                      !
    !                                                                         !
    ! In a) the Hartree energy and potential are only due to electrons, and   !
    !       rho is the electronic density.                                    !
    ! In b) and c) the Hartree energy and potential are due to the whole      !
    !       molecule, i.e. electrons and smeared ions, and rho is the total   !
    !       density.                                                          !
    ! In c) the dielectric permittivity must be supplied through eps.         !
    ! However, this is assumed to be a scalar value, we thus assume that in   !
    ! the solvated case the whole cell is filled with bulk-permittivity die-  !
    ! lectric.                                                                !
    !                                                                         !
    ! Note that because this subroutine cannot know whether the supplied      !
    ! density is electronic or total, it must be spin-unaware, i.e. if the    !
    ! system is spin-polarized, the electronic density must be summed across  !
    ! spins prior to calling this subroutine and the second component of the  !
    ! potential must be initialized from phi after the call.                  !
    !                                                                         !
    ! Open (Dirichlet) BCs and PBCs are supported.                            !
    ! Note that in PBCs the density is shifted such that it integrates to 0   !
    ! (the zero-mode is removed) -- unless the counterion model is used (then !
    ! the counterions take care of ensuring neutrality).                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi,               intent=out, the calculated Hartree potential       !
    !   rho,               intent=in,  the input density, electronic or total !
    !   grid,              intent=in,  the fine grid.                         !
    !   cell,              intent=in,  the simulation cell.                   !
    !   uniform_eps, (opt) intent=in,  the dielectric permittivity (scalar)   !
    !                                  to use for approximate BCs             !
    !   eps_full, (opt)    intent=in,  the dielectric permittivity (for all   !
    !                                  points)                                !
    !   eps_half, (opt)    intent=in,  the dielectric permittivity (for all   !
    !                                  points, shifted by half a grid point)  !
    ! The above optional arguments default to 1, or arrays of 1, if omitted.  !
    !   solvation_terms, (opt, inout), if passed, its components will be      !
    !                                  populated with a breakdown of energy   !
    !                                  terms relevant in solvation.           !
    !   elements, (opt)    intent=in,  elements array passed on to  multigrid_!
    !                                  polarization_data()                    !
    !   no_dump, (opt)     intent=in,  logical parameter to avoid printing any!
    !                                  output or 3D datafiles.                !
    ! Returns:                                                                !
    !   calculated Hartree potential due to rho.                              !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 10/12/2010.                                  !
    ! Updated by James C. Womack, 02/2017, to reflect new DL_MG interface, in !
    ! which the defect correction is performed within DL_MG                   !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root, pub_my_proc_id, comms_barrier
    use constants, only: CRLF, FOUR_PI, stdout, VERBOSE, garbage_real
    use integrals, only: integrals_product_on_grid
    use is_poisson, only: SOLVATION_ENERGY_TERMS, zero_solvation_terms
    use is_solvation_boltzmann, only: IS_PBE_MG_RETURN, boltzmann_ions, &
         n_boltzmann_ions, max_expcap_dl_mg, gamma, &
         implicit_solvent_boltzmann_main, &
         implicit_solvent_boltzmann_shifted_hartree, &
         implicit_solvent_boltzmann_mg_neutralisation_method
    use rundat, only: pub_finite_difference_order, &
         pub_forces_needed, pub_is_bc_coarseness, &
         pub_is_bc_surface_coarseness, pub_is_pbe, pub_is_pbe_exp_cap, &
         pub_is_pbe_temperature, pub_mg_pbe_use_fas, pub_is_smeared_ion_rep, &
         pub_multigrid_bc_is_periodic, pub_multigrid_bc_is_zero, pub_output_detail, &
         pub_mg_defco_fd_order, pub_mg_use_error_damping, &
         pub_mg_tol_res_rel, pub_mg_tol_res_abs, &
         pub_mg_tol_pot_rel, pub_mg_tol_pot_abs, &
         pub_mg_tol_newton_rel, pub_mg_tol_newton_abs, &
         pub_mg_tol_vcyc_rel, pub_mg_tol_vcyc_abs, &
         pub_mg_max_iters_defco, pub_mg_max_iters_vcycle, &
         pub_mg_max_iters_newton, pub_external_bc_from_cube, &
         pub_mg_vcyc_smoother_iter_pre, pub_mg_vcyc_smoother_iter_post, &
         pub_mg_max_res_ratio, pub_mg_tol_mu_rel, pub_mg_tol_mu_abs, &
         pub_is_pbe_neutralisation_scheme, pub_mg_continue_on_error, &
         pub_is_solvation_properties, pub_inner_loop_iteration, &
         pub_mg_tol_cg_res_rel, pub_mg_tol_cg_res_abs, pub_mg_max_iters_cg, &
         pub_mg_use_cg
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_feature_not_supported, utils_assert, utils_trace_in, &
         utils_trace_out
    use visual, only: visual_scalarfield_read
    use ion, only: ELEMENT

#if defined (__INTEL_COMPILER) && defined (DEBUG)
!$  use omp_lib
#endif

    ! jd: Arguments
    ! jd: Electronic or total density in real space.
    real(kind=DP), intent(in) :: rho(mg%ld1f,mg%ld2f,mg%max_slabs12)
    ! jd: Calculated potential due to rho (and, in PBE mode, Boltzmann ions)
    real(kind=DP), intent(out) :: phi(mg%ld1f,mg%ld2f,mg%max_slabs12)
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=DP), intent(in), optional :: uniform_eps
    real(kind=DP), intent(in), optional, target &
         :: eps_full(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional, target &
         :: eps_half(mg%ld1f,mg%ld2f,mg%max_slabs12,3)
    type(SOLVATION_ENERGY_TERMS), intent(inout), optional :: solvation_terms
    type(ELEMENT), intent(in), optional :: elements(:)
    logical, intent(in), optional :: no_dump

    ! jd: Boundary conditions
    real(kind=DP), allocatable :: bound(:,:,:), bound_ext(:,:,:)
    integer                    :: i1, i2, i3, loc_i3

    logical       :: pbc
    logical       :: is_linearised
    real(kind=DP) :: hartree_from_phi_formula
    real(kind=DP) :: finite_box_correction
    type(SOLVATION_ENERGY_TERMS) :: loc_solvation_terms
    type(IS_PBE_MG_RETURN) :: mg_return

    ! JCW: Local variables required to call DL_MG with defect correction
#ifdef HAVE_DL_MG
    real(kind=DP), dimension(:,:,:), pointer   :: eps_full_loc => NULL()
    real(kind=DP), dimension(:,:,:,:), pointer :: eps_half_loc => NULL()
    real(kind=DP)                              :: exp_cap
    character(len=DL_MG_MAX_ERROR_STRING)      :: err_msg
    integer                                    :: err_msg_length
#endif

    logical :: calculate_bcs
    logical :: create_dumps
    integer :: ierr ! jd: Error flag

    integer :: dlmg_nwt_iteration, dlmg_defco_iteration
    integer :: dlmg_nwt_ierror, dlmg_defco_ierror
    real(kind=DP) :: dlmg_nwt_solution_norm, dlmg_defco_solution_norm
    real(kind=DP) :: dlmg_nwt_residual_norm, dlmg_defco_residual_norm
    real(kind=DP) :: dlmg_nwt_target_residual_norm, dlmg_defco_target_residual_norm
    real(kind=DP) :: dlmg_defco_error_norm

    logical, save :: fd_order_discrepancy_warning_issued = .false.

    character(len=*), parameter :: myself = "multigrid_calculate_hartree"

#if defined (__INTEL_COMPILER) && defined (DEBUG)
    ! Used in workaround for possible compiler bug with Intel compiler in
    ! debug mode with > 1 thread
!$  integer              :: nthreads_save
#endif

    !------------------------------------------------------------------------

    call timer_clock(myself,1)
    call utils_trace_in('multigrid_calculate_hartree')

    if (present(no_dump)) then
       create_dumps = .not. no_dump
    else
       create_dumps = .true.
    end if

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    ! JCW: Issue warning about discrepancy in order of FD used for defect correction
    ! JCW: and computation of electrostatic energy and finite box correction if
    ! JCW: user has requested different values of pub_mg_defco_fd_order and
    ! JCW: pub_finite_difference_order
    if (pub_on_root .and. (.not. fd_order_discrepancy_warning_issued) .and. &
         (pub_finite_difference_order /= pub_mg_defco_fd_order) .and. create_dumps) then
       write(stdout,'(a,i2,a,i2,a)') "WARNING: There is a discrepancy in the &
            &finite difference order used to compute"//CRLF//&
            "the electrostatic energy and finite box correction, and the &
            &finite difference"//CRLF//&
            "order used in the defect correction procedure."//CRLF//&
            "FD order = ",pub_finite_difference_order,&
            " for electrostatic energy and finite box correction."//CRLF//&
            "FD order = ",pub_mg_defco_fd_order," for defect correction."//CRLF//&
            "The discretization error is computed using the electrostatic energy &
            &and finite "//CRLF//&
            "box correction and thus depends upon the FD order used. For a &
            &consistent"//CRLF//&
            "treatment of discretization error, it is recommended that you set &
            &these to the"//CRLF//&
            "same order, i.e. "//CRLF//&
            "  finite_difference_order == mg_defco_fd_order."//CRLF
       fd_order_discrepancy_warning_issued = .true.
    end if

    ! jd: Apolar terms might have been calculated in advance (calc_apolar_term),
    !     copy them over.
    if(present(solvation_terms)) then
       loc_solvation_terms%E_apolar_cavitation = &
            solvation_terms%E_apolar_cavitation
       loc_solvation_terms%E_apolar_disrep = solvation_terms%E_apolar_disrep
    end if

    ! jd: If the Boltzmann branch below does not initialise these, they'll be
    !     left at zero.
    call zero_solvation_terms(loc_solvation_terms, zero_apolar = .false., &
         zero_vacuum = .false.)

    if (multigrid_use_onetep_defco) then
       ! JCW: Use the (original) defect correction solver in ONETEP
       ! JCW: This is activated with the DEVEL_CODE value "MG:USE_ONETEP_DEFCO=T:MG"

       ! jd: Allocate the array that will hold the BCs for the MG calculation
       allocate(bound(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check('multigrid_calculate_hartree','bound',ierr)

       ! JCW: ONETEP defect correction only supports fully open (non-periodic)
       ! JCW: BCs
       if (any(pub_multigrid_bc_is_periodic)) then
          call utils_abort("Error in "//myself//": ONETEP-native defect &
               &correction code does not support periodic BCs.")
       end if

       ! JCW: Default is to calculate Dirichlet BCs
       calculate_bcs = .true.

       if (any(pub_multigrid_bc_is_zero)) then
          ! JCW: Boundary condition calculation only supports fully zero BCs or
          ! JCW: fully computed Dirichlet BCs (cannot mix computed and zero)
          call utils_assert(all(pub_multigrid_bc_is_zero),"Error in "//myself//": BCs must be &
               &all zero or all computed Dirichlet -- combinations for different &
               &directions are not currently supported.")
          ! JCW: Zero BCs, no need to compute
          calculate_bcs = .false.
       end if

       if (pub_external_bc_from_cube) then
          ! ebl: read in external boundary conditions from cube file
          allocate(bound_ext(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
          call utils_alloc_check('multigrid_calculate_hartree', &
               'bound_ext', ierr)

          ! implicitly assume no scaling is needed for the potential:
          ! require using atomic units
          call visual_scalarfield_read(bound_ext, grid, mg, cell, &
               '_POT_EXT_BC')

          ! ebl: zero everything but the boundaries
          do loc_i3 = 1, grid%max_slabs12
             i3 = loc_i3 + grid%first_slab12(pub_my_proc_id) - 1
             do i2 = 2, mg%pq2f - 1
                do i1 = 2, mg%pq1f - 1
                   if ((i3 /= 1) .and. (i3 /= mg%pq3f)) then
                      bound_ext(i1, i2, loc_i3) = 0.0_DP
                   end if
                end do
             end do
          end do

          ! ebl: The contents of the cube file are the negative of
          !      the boundary conditions due a sign difference between then
          !      internal quantity "pot" and the electrostatic potential
          bound(:,:,:) = - bound_ext(:,:,:)

          call comms_barrier

          deallocate(bound_ext,stat=ierr)
          call utils_dealloc_check('multigrid_calculate_hartree', &
               'bound_ext',ierr)
       else if (calculate_bcs) then
          ! jd: Determine the correct boundary conditions
          call multigrid_prepare_bound_cond(bound,rho,grid,cell,uniform_eps, &
               pub_is_bc_coarseness,pub_is_bc_surface_coarseness)
       else
          bound = 0.0_DP
       end if

       ! jd: Use the multigrid solver to get the Hartree potential
       call multigrid_defect_corr_solver(pot = phi, & ! output
            rho = rho, bound = bound, &               ! input
            eps_full = eps_full, eps_half = eps_half, grid = grid, cell = cell) ! input

       ! jd: Clean up
       deallocate(bound,stat=ierr)
       call utils_dealloc_check('multigrid_calculate_hartree','bound',ierr)

    else
       ! JCW: Use the newer defect correction solver in DL_MG (default)
#ifdef HAVE_DL_MG
       ! jd: If permittivity specified in 'eps_half', use it,
       !     else use unity everywhere, but must allocate the array.
       ! JCW: Do the same for eps_full, as it is required in the reorganized
       ! JCW: DL_MG interface
       if(present(eps_half)) then
          call utils_assert(present(eps_full),"Error in "//myself//": &
               &If eps_half is present, then eps_full must also be present.")
          eps_half_loc => eps_half
          eps_full_loc => eps_full
       else
          call utils_assert(.not.present(eps_full),"Error in "//myself//": &
               &If eps_half is NOT present, then eps_full must also NOT be present.")
          allocate(eps_half_loc(mg%ld1f,mg%ld2f,mg%max_slabs12,3),stat=ierr)
          call utils_alloc_check(myself,'=>eps_half_loc',ierr)
          allocate(eps_full_loc(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
          call utils_alloc_check(myself,'=>eps_full_loc',ierr)
          eps_half_loc = 1.0_DP
          eps_full_loc = 1.0_DP
       end if

       ! JCW: Default is to calculate Dirichlet BCs
       calculate_bcs = .true.

       if (any(pub_multigrid_bc_is_zero).or.any(pub_multigrid_bc_is_periodic)) then
          ! JCW: Boundary condition calculation only supports combinations of
          ! JCW: zero and periodic BCs (no need to compute BCs on cell faces) or
          ! JCW: fully computed Dirichlet BCs

          call utils_assert(all(pub_multigrid_bc_is_zero(:).or.pub_multigrid_bc_is_periodic(:)),&
               "Error in "//myself//": Mixing of computed Dirichlet BCs (open BC) &
               &and periodic/zero BCs is not currently supported. The BCs should &
               &either be all open, or some combination of zero and periodic.")
          ! JCW: Zero BCs, no need to compute Dirichlet BCs
          calculate_bcs = .false.

          ! JCW: In PBCs, we must subtract off the average density (G=0) component
          ! JCW: of the potential in reciprocal space to be consistent with the
          ! JCW: other electrostatic energy terms (Ewald).
          ! JCW: --> For fully periodic BCs, this is done automatically within
          ! JCW:     DL_MG, so the density and potential can be passed to
          ! JCW:     DL_MG unchanged.
          ! JCW: --> For mixed periodic/non-periodic BCs, this is not done within
          ! JCW:     DL_MG
       end if


       if (pub_external_bc_from_cube) then
          ! ebl: read in external boundary conditions from cube file
          allocate(bound_ext(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
          call utils_alloc_check('multigrid_calculate_hartree', &
               'bound_ext', ierr)

          ! implicitly assume no scaling is needed for the potential:
          ! require using atomic units
          call visual_scalarfield_read(bound_ext, grid, mg, cell, &
               '_POT_EXT_BC')

          ! ebl: zero everything but the boundaries
          do loc_i3 = 1, grid%num_my_slabs12
             i3 = loc_i3 + grid%first_slab12(pub_my_proc_id) - 1
             do i2 = 2, mg%pq2f - 1
                do i1 = 2, mg%pq1f - 1
                   if ((i3 /= 1) .and. (i3 /= mg%pq3f)) then
                      bound_ext(i1, i2, loc_i3) = 0.0_DP
                   end if
                end do
             end do
          end do

          ! ebl: The contents of the cube file are the negative of
          !      the boundary conditions due a sign difference between then
          !      internal quantity "phi" and the electrostatic potential
          phi(:,:,:) = - bound_ext(:,:,:)

          call comms_barrier

          deallocate(bound_ext,stat=ierr)
          call utils_dealloc_check('multigrid_calculate_hartree', &
               'bound_ext',ierr)
       else if (calculate_bcs) then
          ! JCW: Determine boundary conditions and place in potential array to
          ! JCW: be passed to DL_MG
          call multigrid_prepare_bound_cond(phi,rho,grid,cell,uniform_eps, &
               pub_is_bc_coarseness,pub_is_bc_surface_coarseness)
       else
          phi = 0.0_dp
       end if

       if (create_dumps) then
          call mg_debug_dump_input(grid, cell, rho, eps_full, eps_half, phi)
       end if

       if (pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a,i0)') &
               "MG DL_MG: Calling solver with defect correction order = ", &
               pub_mg_defco_fd_order
       end if

       if(pub_is_pbe_exp_cap == 0.0_DP) then
          exp_cap = max_expcap_dl_mg
       else
          exp_cap = pub_is_pbe_exp_cap
       end if

#if defined (__INTEL_COMPILER) && defined (DEBUG)
    ! JCW: Workaround for possible compiler bug for Intel compiler in debug
    ! JCW: mode, with > 1 OpenMP thread. This occurs in DL_MG and is possible
    ! JCW: due to misreading of array range info in the mg struct.
    ! JCW: Since the issue only occurs with > 1 OpenMP thread, we temporarily
    ! JCW: change the number of threads while DL_MG is active:
!$  if (pub_on_root) then
!$     write(stdout,'(a)') &
!$          "MG WARNING: Temporarily reducing number of OpenMP threads to 1 to work around"//CRLF//&
!$          "            possible compiler bug when DL_MG is compiled with Intel compiler"//CRLF//&
!$          "            in debug mode."
!$  end if
!$  nthreads_save = omp_get_max_threads()
!$  call omp_set_num_threads(1)
    ! JCW: This may be the same possible compiler bug described in ONETEP
    ! JCW: CCPForge bug #1410.
#endif

       ! jd: This is needed to later determine if we actually used counterion
       !     neutralisation. If not, the ratios will stay set to garbage_real.
       mg_return%used_neutralisation_ratios(1:n_boltzmann_ions) = garbage_real

       call timer_clock('dl_mg_solver_calls',1)
       ! JCW: New DL_MG interface exposes one name for dl_mg_solver
       ! JCW: and distinguishes whether to solve Poisson or Poisson-Boltzmann
       ! JCW: equation via arguments supplied.
       is_linearised = (pub_is_pbe == 'LINEARISED')
       if(pub_is_pbe == 'FULL' .or. pub_is_pbe == 'LINEARISED') then
          ! *************************
          ! *** Poisson-Boltzmann ***
          ! *************************
          call dl_mg_solver(eps_full_loc, eps_half_loc, -FOUR_PI, rho,&
               -FOUR_PI, pub_is_pbe_temperature, n_boltzmann_ions, &
               boltzmann_ions(1:n_boltzmann_ions)%concentration_user, &
               boltzmann_ions(1:n_boltzmann_ions)%charge, &
               phi, pub_mg_defco_fd_order, &
               tol_res_rel = pub_mg_tol_res_rel, &
               tol_res_abs = pub_mg_tol_res_abs, &
               tol_pot_rel = pub_mg_tol_pot_rel, &
               tol_pot_abs = pub_mg_tol_pot_abs, &
               tol_newton_res_rel = pub_mg_tol_newton_rel, &
               tol_newton_res_abs = pub_mg_tol_newton_abs, &
               tol_vcycle_res_rel = pub_mg_tol_vcyc_rel, &
               tol_vcycle_res_abs = pub_mg_tol_vcyc_abs, &
               tol_mu_rel = pub_mg_tol_mu_rel, &
               tol_mu_abs = pub_mg_tol_mu_abs, &
               tol_cg_res_rel = pub_mg_tol_cg_res_rel, &
               tol_cg_res_abs = pub_mg_tol_cg_res_abs, &
               max_iters_defco  = pub_mg_max_iters_defco, &
               max_iters_vcycle = pub_mg_max_iters_vcycle, &
               max_iters_newton = pub_mg_max_iters_newton, &
               max_iters_cg     = pub_mg_max_iters_cg, &
               use_cg = pub_mg_use_cg, &
               linearised=is_linearised, steric_weight=gamma, expcap=exp_cap, &
               use_fas=pub_mg_pbe_use_fas, &
               use_damping=pub_mg_use_error_damping, &

               neutralisation_method = &
               implicit_solvent_boltzmann_mg_neutralisation_method(mg,grid,rho), &
               neutralisation_ion_ratios = &
               boltzmann_ions(:)%neutralisation_ratio, &

               betamu_electrostatic = &
               mg_return%betamu_elec(1:n_boltzmann_ions), &
               used_ion_concentrations = &
               mg_return%used_ion_concentration(1:n_boltzmann_ions), &
               used_neutralisation_ratios = &
               mg_return%used_neutralisation_ratios(1:n_boltzmann_ions), &
               steric_weight_average = mg_return%steric_weight_average, &

               v_iterations = (/pub_mg_vcyc_smoother_iter_pre, &
               pub_mg_vcyc_smoother_iter_post/), &
               mg_max_conv_rate = pub_mg_max_res_ratio, &
               ierror=ierr)
       elseif(pub_is_pbe == 'NONE') then
          ! ********************
          ! *** Pure Poisson ***
          ! ********************
          call dl_mg_solver(eps_full_loc, eps_half_loc, -FOUR_PI, rho,&
               phi, pub_mg_defco_fd_order, &
               tol_res_rel = pub_mg_tol_res_rel, &
               tol_res_abs = pub_mg_tol_res_abs, &
               tol_pot_rel = pub_mg_tol_pot_rel, &
               tol_pot_abs = pub_mg_tol_pot_abs, &
               tol_vcycle_res_rel = pub_mg_tol_vcyc_rel, &
               tol_vcycle_res_abs = pub_mg_tol_vcyc_abs, &
               max_iters_defco  = pub_mg_max_iters_defco, &
               max_iters_vcycle = pub_mg_max_iters_vcycle, &
               max_iters_cg     = pub_mg_max_iters_cg, &
               use_cg = pub_mg_use_cg, &
               use_damping = pub_mg_use_error_damping, &
               v_iterations = (/pub_mg_vcyc_smoother_iter_pre, &
               pub_mg_vcyc_smoother_iter_post/), &
               mg_max_conv_rate = pub_mg_max_res_ratio, &
               ierror=ierr)
       else
          call utils_abort('Internal error. Unrecognized pub_is_pbe')
       endif
       call timer_clock('dl_mg_solver_calls',2)

       if (create_dumps) then
          call mg_debug_dump_output(grid, cell, phi, quit = .true., &
               mg_return = mg_return)
       end if

#if defined (__INTEL_COMPILER) && defined (DEBUG)
    ! JCW: Workaround for possible compiler bug for Intel compiler in debug
    ! JCW: mode, with > 1 OpenMP thread. This occurs in DL_MG and is possible
    ! JCW: due to misreading of array range info in the mg struct.
    ! JCW: Since the issue only occurs with > 1 OpenMP thread, we temporarily
    ! JCW: change the number of threads while DL_MG is active:
!$  if (pub_on_root) then
!$     write(stdout,'(a)') &
!$          "MG WARNING: Restoring number of OpenMP threads to previous value."
!$  end if
!$  call omp_set_num_threads(nthreads_save)
    ! JCW: This may be the same possible compiler bug described in ONETEP
    ! JCW: CCPForge bug #1410.
#endif

       if(ierr /= 0) then
          call dl_mg_error_string(ierr,err_msg,err_msg_length)

          if(pub_mg_continue_on_error) then
             if(pub_on_root) then
                write(stdout,'(a)') "MG DL_MG: ------------------------- &
                     &SERIOUS WARNING: -------------------------"
                write(stdout,'(a)') "MG DL_MG: Could not achieve convergence to&
                     & prescribed tolerance settings!"
                write(stdout,'(a)') "MG "//trim(err_msg)//'.'
                call dl_mg_solver_status(dlmg_nwt_iteration, &
                     dlmg_nwt_solution_norm, dlmg_nwt_residual_norm, &
                     dlmg_nwt_target_residual_norm, dlmg_nwt_ierror, &
                     dlmg_defco_iteration, dlmg_defco_solution_norm, &
                     dlmg_defco_residual_norm, dlmg_defco_error_norm, &
                     dlmg_defco_target_residual_norm, dlmg_defco_ierror)

                write(stdout,'(a,i3,a,i6,a,e10.3,a)') &
                     "MG DL_MG: Newton: error code:", dlmg_nwt_ierror, &
                     ", n_iter:", dlmg_nwt_iteration, ", sol_norm:", &
                     dlmg_nwt_solution_norm, "."
                write(stdout,'(a,e10.3,a,e10.3,a)') &
                     "MG DL_MG: Newton: res_norm:", dlmg_nwt_residual_norm, &
                     ", target_res_norm:", dlmg_nwt_target_residual_norm, &
                     "."
                write(stdout,'(a,i3,a,i6,a,e10.3,a)') &
                     "MG DL_MG: DEFCO:  error code:", dlmg_defco_ierror, &
                     ", n_iter:", dlmg_defco_iteration, ", sol_norm:", &
                     dlmg_defco_solution_norm, "."
                write(stdout,'(a,e10.3,a,e10.3,a)') &
                     "MG DL_MG: DEFCO:  res_norm:", dlmg_defco_residual_norm, &
                     ", target_res_norm:", dlmg_defco_target_residual_norm, "."
                write(stdout,'(a,e10.3,a)') &
                     "MG DL_MG: DEFCO:  err_norm:", dlmg_defco_error_norm, "."
                write(stdout,'(a)') "MG DL_MG: --------------------------&
                     &------------------------------------------"

             end if
          else
             call utils_abort("Error in the DL_MG multigrid solver: "//CRLF//&
                  trim(err_msg))
          end if
       end if

       if (pub_on_root .and. pub_output_detail >= VERBOSE) &
            write(stdout,'(a)') "MG DL_MG: Solver finished"

       ! Clean up temporary arrays / pointers
       call utils_assert(associated(eps_full_loc),"Error in "//myself//": &
            &Unexpected non-associated status for eps_full_loc pointer.")
       call utils_assert(associated(eps_half_loc),"Error in "//myself//": &
            &Unexpected non-associated status for eps_half_loc pointer.")
       if(.not. present(eps_full)) then
          deallocate(eps_full_loc,stat=ierr)
          call utils_dealloc_check(myself,'=>eps_full_loc',ierr)
       else
          nullify(eps_full_loc)
       end if
       if(.not. present(eps_half)) then
          deallocate(eps_half_loc,stat=ierr)
          call utils_dealloc_check(myself,'=>eps_half_loc',ierr)
       else
          nullify(eps_half_loc)
       end if
#else
       call utils_feature_not_supported('DL_MG multigrid solver')
#endif
    end if

    ! ab: produce polarization data in solvent
    if(present(elements) .and. ((pub_is_solvation_properties .and. &
         pub_inner_loop_iteration == -2) .or. multigrid_debug)) then
       call multigrid_polarization_data(grid, cell, elements, phi, eps_full)
    end if

    ! jd: If smeared ions are in use and forces are needed, store phi,
    !     it will be needed for the force calculation, where the density
    !     is no longer accessible.
    if(pub_forces_needed .and. pub_is_smeared_ion_rep) then
       phi_vac_for_forces = phi
    end if

    ! JCW: The finite box approximation is only present for non-periodic BCs
    if (.not.any(pub_multigrid_bc_is_periodic)) then
       ! JCW: Fully Open BCs
       ! jd: Calculate \surfint phi * eps * grad_phi dS, the correction
       finite_box_correction = multigrid_finite_box_correction(phi,grid,eps_full)
       pbc = .false.
    else if (all(pub_multigrid_bc_is_periodic)) then
       ! JCW: Fully periodic BCs (no finite box approximation)
       finite_box_correction = 0.0_DP
       pbc = .true.
    else
       ! JCW: Mixed open/periodic BCs
       ! TODO Currently calculation of the finite box correction assumes fully
       !      open BCs --> implement support for mixed BCs
       finite_box_correction = garbage_real
    end if

    ! jd: Calculate E_Hartree, the returned value
    ! jd: In OBCs or when counterions are used this is as usual.
    if(.not. pbc .or. &
         pub_is_pbe_neutralisation_scheme(1:11) == 'COUNTERIONS') then
       multigrid_calculate_hartree = &
            0.5_DP * integrals_product_on_grid(grid, &
            phi, rho, mg%pq1f, mg%pq2f, mg%num_slabs12_pq)
    else
       ! jd: In PBCs we need to shift the charge density so that it integrates
       !     to zero -- unless the counterion model is in use (then the
       !     counterions take care of neutralisation). Delegate this to the
       !     Boltzmann module. This is needed even in the absence of electrolyte
       !     (is_pbe NONE), because we may have a charged system (which invokes
       !     jellium).
       multigrid_calculate_hartree = &
            implicit_solvent_boltzmann_shifted_hartree(grid, mg, mg_return, &
            rho, phi)
    end if

    ! jd: Calculate 1/8pi \int eps (grad phi)^2, the other estimate
    hartree_from_phi_formula = multigrid_electrostatic_energy(phi,grid,eps_full)

    if(pub_is_pbe == 'NONE' .or. pub_is_pbe == 'LINEARISED') then
       if(create_dumps .and. pub_on_root .and. pub_output_detail >= VERBOSE) then
          write(stdout,'(a,f14.6)') &
               "MG Hartree energy from (\int phi*rho dr):           ", &
               multigrid_calculate_hartree
          write(stdout,'(a,f14.6)') &
               "MG Hartree energy from (\int eps (grad phi)^2 dr):  ", &
               hartree_from_phi_formula
          if (.not.any(pub_multigrid_bc_is_periodic)) then
             ! JCW: Fully open BCs, finite box approximation present and supported
             write(stdout,'(a,f14.6,a,f10.5,a)') &
                 "MG Discrepancy due to finite box and discretization:", &
                  hartree_from_phi_formula - multigrid_calculate_hartree, ' (', &
                  100.0_DP*(hartree_from_phi_formula - multigrid_calculate_hartree) /&
                  multigrid_calculate_hartree,'%)'
             write(stdout,'(a,f14.6)') &
                  "MG Correction to account for finite box:            ", &
                  finite_box_correction
             hartree_from_phi_formula = hartree_from_phi_formula+finite_box_correction
             write(stdout,'(a,f14.6,a,f10.5,a)') &
                  "MG Remaining discrepancy, due to discretization:    ", &
                  hartree_from_phi_formula - multigrid_calculate_hartree, ' (', &
                  100.0_DP*(hartree_from_phi_formula - multigrid_calculate_hartree) /&
                  multigrid_calculate_hartree,'%)'
          else if (all(pub_multigrid_bc_is_periodic)) then
             ! JCW: Fully periodic BCs -- finite box approximation not present
             ! JCW:
             ! JCW: The two forms used for the electrostatic energy are equivalent
             ! JCW: in PBCs: there is no finite box correction and the
             ! JCW: solution from DL_MG is the potential minus it's average, i.e.
             ! JCW: with the zero mode removed.
             ! JCW:
             ! JCW: In PBCs DL_MG solves
             ! JCW:   A[\epsilon]\phi' = f
             ! JCW: where
             ! JCW:   A[\epsilon]\phi' = \nabla.(\epsilon\nabla\phi')
             ! JCW:   f = -4\pi rho'
             ! JCW: and where \rho' is the charge density with compensating
             ! JCW: density added to ensure neutrality and \phi' is the
             ! JCW: electrostatic potential for \rho'.
             ! JCW: The Hartree energy (multigrid_calculate_hartree) is
             ! JCW:   (1/2) \int \phi'(r) \rho(r) dr
             ! JCW: integrating over the simulation cell, with \rho(r) the
             ! JCW: charge density without the compensating charge.
             ! jd 2020.02: ^^^ the above is no longer the case, we now exclude
             !                 the zero mode from the density in
             !                 multigrid_calculate_hartree()
             ! JCW: The alternative form (hartree_from_phi_formula) is
             ! JCW:   (1/8)\pi \int \epsilon(r) (\nabla \phi'(r))^2 dr
             ! JCW: which is equivalent to
             ! JCW:   (1/2) \int \phi'(r) \rho'(r) dr
             ! JCW: --> Note that this is the interaction of \rho' with \phi'
             ! JCW:     i.e. the charge density + compensating charge and the
             ! JCW:     corresponding potential.
             ! JCW: Since the zero mode (or G=0) component of the potential is
             ! JCW: removed by DL_MG under PBCs
             ! JCW:   (1/2) \int \phi'(r) \rho_cmp = 0
             ! JCW: where \rho_cmp is the uniform compensating density.
             ! JCW: --> Consequently, the two forms are formally equivalent and
             ! JCW:     any difference between the values is due to discretization
             ! JCW:
             ! JCW: In short:
             ! JCW: 1. The relationship between multigrid_calculate_hartree and
             ! JCW:    hartree_from_phi_formula may be obtained by integration by
             ! JCW:    parts. In PBCs, the surface term (finite_box_correction)
             ! JCW:    vanishes.
             ! JCW: 2. The lack of compensating charge in the density used
             ! JCW:    to compute multigrid_calculate_hartree is irrelevant,
             ! JCW:    since the removal of the G=0 component from the PBC
             ! JCW:    potential ensures that the average density in the cell
             ! JCW:    makes no contribution to the energy.
             write(stdout,'(a,f14.6,a,f10.5,a)') &
                  "MG Discrepancy due to discretization:               ", &
                  hartree_from_phi_formula - multigrid_calculate_hartree, ' (', &
                  100.0_DP*(hartree_from_phi_formula - multigrid_calculate_hartree) /&
                  multigrid_calculate_hartree,'%)'
             !write(stdout,'(a,f14.6)') &
             !     "MG No finite box approximation in periodic BCs"
          else
             ! JCW: Mixed open/periodic BCs
             ! TODO Currently calculation of the finite box correction assumes fully
             !      open BCs --> implement support for mixed BCs
             ! [ Do nothing ]
          end if ! BCs

       end if ! on root in VERBOSE
    end if ! FULL/LINEARISED/NONE

    if(pub_is_pbe == 'FULL' .or. pub_is_pbe == 'LINEARISED') then
       call implicit_solvent_boltzmann_main(loc_solvation_terms, mg, grid, &
            mg_return, phi, multigrid_calculate_hartree, pbc)
    end if

    if(present(solvation_terms)) then
       loc_solvation_terms%E_elec_fixed = multigrid_calculate_hartree
       solvation_terms = loc_solvation_terms
    end if

    call utils_trace_out('multigrid_calculate_hartree')
    call timer_clock(myself,2)

  end function multigrid_calculate_hartree
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine multigrid_exit
    !=========================================================================!
    ! Deinitialises the multigrid.                                            !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, September 2013.                              !
    !=========================================================================!

    use rundat, only: pub_forces_needed, pub_is_smeared_ion_rep
    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr

    ! -------------------------------------------------------------------------

#ifdef HAVE_DL_MG
    call dl_mg_free()
#endif

    ! jd: Deallocate an array to hold phi_vac, required if smeared-ion
    !     correction to forces is needed.
    if(pub_forces_needed .and. pub_is_smeared_ion_rep) then
       deallocate(phi_vac_for_forces,stat=ierr)
       call utils_dealloc_check('multigrid_exit','phi_vac_for_forces',ierr)
    end if

    ! JCW: Reset module variables
    multigrid_initialised = .false.
    multigrid_debug = .false.
    multigrid_use_onetep_defco = .false.

  end subroutine multigrid_exit
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !---------------------------------------------------------------------------!
  !                     P r i v a t e   s u b r o u t i n e s                 !
  !---------------------------------------------------------------------------!
  ! The following routines are used to perform the defect correction within   !
  ! ONETEP, rather than in DL_MG and represent the original implementation of !
  ! the defect correction (from which the DL_MG implementation was derived).  !
  !                                                                           !
  ! By default, the defect correction is performed within DL_MG, though the   !
  ! original ONETEP defect correction (which calls DL_MG only as a second-    !
  ! order solver) can be enabled by adding the DEVEL_CODE                     !
  ! "MG:USE_ONETEP_DEFCO=T:MG" to the ONETEP input file.                      !
  !---------------------------------------------------------------------------!

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine multigrid_defect_corr_solver(pot, & ! out
       rho, bound, eps_full, eps_half, grid, cell)     ! in
    !=========================================================================!
    ! For a given electronic density (rho) and dielectric functional (eps)    !
    ! this subroutine will compute the electrostatic potential (pot) which is !
    ! the solution to the Poisson equation:                                   !
    !                                                                         !
    !        -  \nabla \cdot [ eps \nabla (pot) ] = 4*pi*rho                  !
    !                                                                         !
    ! This subroutine uses the defect correction method with multigrid to     !
    ! obtain a solution to the Poisson equation which has high order accuracy !
    ! determined by the finite difference operators used.                     !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   pot,       intent=out,             the potential we seek              !
    !   rho,       intent=in,              the density                        !
    !                                                                         !
    !   bound,     intent=in,              the initial guess for the potential!
    !              optional                and its fixed boundary conditions  !
    !              If omitted, will default to 0.0 everywhere.                !
    !   eps_full,  intent=in,              the dielectric functional on       !
    !              optional                'full' grid points                 !
    !              If omitted, will default to 1.0 everywhere.                !
    !   eps_half,  intent=in,              the dielectric functional          !
    !              optional                interpolated to 'half' grids       !
    !              If omitted, will default to 1.0 everywhere.                !
    !-------------------------------------------------------------------------!
    ! Defect correction written for castep by HH Helal 29/04/2010.            !
    ! Mended for onetep by HH Helal 7/6/2010.                                 !
    ! Adapted to work in parallel by Jacek Dziedzic 11/06/2010.               !
    ! Updated by James C. Womack, 02/2017, to reflect new DL_MG interface.    !
    !=========================================================================!

    use constants, only: stdout, VERBOSE
    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root, comms_reduce
    use constants, only: FOUR_PI
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_mg_tol_pot_abs, pub_output_detail, &
         pub_mg_tol_res_abs, pub_mg_max_iters_defco, &
         pub_mg_defco_fd_order, pub_mg_use_error_damping, &
         pub_mg_max_iters_cg
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_flush, &
         utils_alloc_check, utils_dealloc_check, utils_sanity_check, &
         utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)                  :: pot(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)                   :: rho(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional         :: bound(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional, target :: eps_full(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional, target :: eps_half(mg%ld1f,mg%ld2f,mg%max_slabs12,3)
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell

    ! Internal variables
    integer, save :: total_calls = 0
    integer       :: ierr ! jd: Error flag
    integer       :: i1, i2, i3, islab12
    integer       :: owner_proc
    real(kind=DP) :: verbose_actual_y, verbose_actual_z
    real(kind=DP) :: res_norm, error_norm
    real(kind=DP) :: defect_norm, defect_norm_stpt, defect_norm_bt
    integer       :: damping_iterations
    integer       :: damping_iterations_total
    character(len=*), parameter :: myself = 'multigrid_defect_corr_solver'
    character(len=*), parameter :: explanation = "This could mean that the&
         & multigrid is too coarse to represent the quickly changing&
         & quantities (such as the permittivity). Setting a tighter value,&
         & such as 1E-5, for 'mg_tol_res_abs' might help, as&
         & might using 'mg_use_error_damping T'. The most robust solution&
         & will be making the grid finer through either 'fine_grid_scale' or&
         & by increasing the KE cutoff."

    ! Defect correction method variables
    real(kind=DP), allocatable :: defect(:,:,:)      ! current defect - total
    real(kind=DP), allocatable :: defect_stpt(:,:,:) ! - src and Poisson term
    real(kind=DP), allocatable :: defect_bt(:,:,:)   ! - Boltzmann term
    real(kind=DP), allocatable :: error(:,:,:)  ! hhh: curr. error approx.
    real(kind=DP), allocatable :: res_local(:,:,:)  ! jd: Workspace
    real(kind=DP), dimension(:,:,:), pointer   :: eps_full_loc => NULL()
    real(kind=DP), dimension(:,:,:,:), pointer :: eps_half_loc => NULL()
    integer                    :: defcorr_iter
    logical                    :: converged
    logical                    :: spokesman
    ! JCW: For output of norm of the potential when multigrid_debug is .true. :
    real(kind=DP)              :: pot_norm_before, pot_norm_after

    !------------------------------------------------------------------------

    call timer_clock(myself,1)
    call utils_trace_in(myself)

    spokesman = pub_on_root .and. pub_output_detail >= VERBOSE

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    total_calls = total_calls + 1

    ! jd: If permittivity specified in 'eps_half', use it,
    !     else use unity everywhere, but must allocate the array.
    ! JCW: Do the same for eps_full, as it is required in the reorganized
    ! JCW: DL_MG interface (as of 01/2017)
    if(present(eps_half)) then
       call utils_assert(present(eps_full),"Error in "//myself//": &
            &If eps_half is present, then eps_full must also be present.")
       eps_half_loc => eps_half
       eps_full_loc => eps_full
    else
       call utils_assert(.not.present(eps_full),"Error in "//myself//": &
            &If eps_half is NOT present, then eps_full must also NOT be present.")
       allocate(eps_half_loc(mg%ld1f,mg%ld2f,mg%max_slabs12,3),stat=ierr)
       call utils_alloc_check(myself,'=>eps_half_loc',ierr)
       allocate(eps_full_loc(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'=>eps_full_loc',ierr)
       eps_half_loc = 1.0_DP
       eps_full_loc = 1.0_DP
    end if

    ! jd: Sanity checks
    call utils_sanity_check(rho,'rho in '//myself)
    if(present(bound)) call utils_sanity_check(bound, 'bound in '//myself)
    call utils_sanity_check(eps_full_loc, '=>eps_full_loc in '//myself)
    call utils_sanity_check(eps_half_loc(:,:,:,1), '=>eps_half_loc/x in '//myself)
    call utils_sanity_check(eps_half_loc(:,:,:,2), '=>eps_half_loc/y in '//myself)
    call utils_sanity_check(eps_half_loc(:,:,:,3), '=>eps_half_loc/z in '//myself)

    call internal_multigrid_verbose_input
    call mg_debug_dump_input(grid, cell, rho, eps_full, eps_half, bound)

    if(spokesman) then
       write(stdout,'(a,i0,a)',advance='no') 'MG (Call #', total_calls, &
            ') Solving with 2nd order solver... '
       call utils_flush(stdout,.true.)
    end if

    ! jd: Allocate workspace for residual
    allocate(res_local(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'res_local',ierr)
    res_local = 0.0_DP ! jd: Takes care of max-num, ld-pt padding

    ! ************************************************************************
    ! *** CALL THE SOLVER (FOR PBE EQUATION PROPER)
    ! ************************************************************************
    call second_order_multigrid_solver(pot,res_local,rho,grid,cell,&
         eps_full_loc,eps_half_loc,bound)
    ! ************************************************************************

    if (multigrid_debug) then
       pot_norm_before = sum(bound(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq)**2)
       pot_norm_after  = sum(pot(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq)**2)
       call comms_reduce('SUM',pot_norm_before)
       call comms_reduce('SUM',pot_norm_after)
       pot_norm_before = sqrt(pot_norm_before)
       pot_norm_after  = sqrt(pot_norm_after)
       if(spokesman) then
          write(stdout,'(a,e23.15)') &
            "MG DEBUG: norm pot before solver call      = ", pot_norm_before
          write(stdout,'(a,e23.15)') &
            "MG DEBUG: norm pot after 1st solver call   = ", pot_norm_after
          call utils_flush(stdout,.true.)
       end if
    end if

    if(spokesman) then
       write(stdout,'(a)') 'done'
       call utils_flush(stdout,.true.)
    end if

    ! ************************************************************************
    ! *** DEFECT CORRECTION PROCEDURE
    ! ************************************************************************
    ! jd: If higher order was requested, employ defect correction
    if (pub_mg_defco_fd_order > 2) then

       damping_iterations_total = 0

       ! Allocate memory for local arrays needed for defect correction
       allocate(defect(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'defect',ierr)
       allocate(defect_stpt(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'defect_stpt',ierr)
       allocate(defect_bt(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'defect_bt',ierr)
       allocate(error(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'error',ierr)

       ! hhh: Error initialized to zero -this is fine for possibly arbitrary BCs
       !      set earlier by optional bound argument in multigrid_poisson_solver
       error(:,:,:) = 0.0_DP

       ! hhh: We now have all the inputs ready to start the defect correction
       !      iteration
       converged = .false.
       defcorr_iter = 0

       if(spokesman) then
          write(stdout,'(a,i0)') 'MG Starting defect correction with &
               &a discretization order of: ', pub_mg_defco_fd_order
          write(stdout,'(a)') "  iter |ST+PT defect|   |BT defect|  &
               &|tot defect|        |error|  |err.eqn res.|"
       end if


       ! **********************************************************************
       ! *** DEFECT CORRECTION LOOP
       ! **********************************************************************
       do while (.not. converged)

          defcorr_iter = defcorr_iter + 1

          if(spokesman) then
             write(stdout,'(a,i4,a)',advance='no') "MG", defcorr_iter,"  "
          end if

         ! hhh: Compute the defect using high order finite difference
         !      Level of accuracy is selected through the input param
         !      discretization_order
         !      defect = source + \nabla \cdot [\epsilon \grad potmg]
         call compute_current_defect(pub_mg_defco_fd_order, rho, &
              eps_full_loc, pot, grid, &
              defect, defect_norm, defect_stpt, defect_norm_stpt, &
              defect_bt, defect_norm_bt)

         call utils_sanity_check(defect,'defect in '//myself)

         if(spokesman) then
            write(stdout,'(e13.7,a,e13.7,a,e13.7,a)',advance='no') &
                 defect_norm_stpt," ",defect_norm_bt," ",defect_norm," -> "
         end if

         call internal_multigrid_verbose_input_defco

         ! hhh: Now that we have the high order defect we very approximately
         !      solve the error equation using second order multigrid method
         !      with a low convergence threshold
         !      - \nabla \cdot [epsilon \grad error] = defect
         defect = defect/(FOUR_PI) ! jd: (*) solver expects 'rho', not 'source'

         ! *****************************************************************
         ! *** CALL THE SOLVER (DEFECT EQUATION)
         ! *****************************************************************
         call second_order_multigrid_solver(error,res_local,defect,grid,cell,&
              eps_full_loc,eps_half_loc,&
              tolerance_override = pub_mg_tol_res_abs,&
              use_linear=.true., der_pot=pot)
         ! *****************************************************************

         ! jd: Compute error and error equation residual norms
         error_norm = sqrt(integrals_product_on_grid(grid, &
              error,error,mg%pq1f,mg%pq2f,mg%num_slabs12_pq))
         res_norm = sqrt(integrals_product_on_grid(grid, &
              res_local,res_local,mg%pq1f,mg%pq2f,mg%num_slabs12_pq))

         if(spokesman) then
            write(stdout,'(e13.7,a,e13.7)') error_norm," ",res_norm
         end if

         ! la: damp the error if needed
         damping_iterations = 0
         if (pub_mg_use_error_damping) then
            call internal_error_damping(error, pot, defect_norm)
            damping_iterations_total = &
                 damping_iterations_total + damping_iterations
         end if

         defect = defect*(FOUR_PI) ! jd: undo (*)

         ! Correct the current approximation to get the new approximation
         pot(:,:,:) = pot(:,:,:) + error(:,:,:)

         ! hhh: Compute error_norm and res_norm
         error_norm = sqrt(integrals_product_on_grid(grid, &
              error,error,mg%pq1f,mg%pq2f,mg%num_slabs12_pq))
         res_norm = sqrt(integrals_product_on_grid(grid, &
              res_local,res_local,mg%pq1f,mg%pq2f,mg%num_slabs12_pq))

         if(spokesman) then
            call utils_flush(stdout,.true.)
         end if

         ! jd: Sanity checks for runaway solutions
         call utils_assert(error_norm < 1D6, 'The error norm is absurdly &
              &high. '//trim(explanation),error_norm)
         call utils_assert(defect_norm < 1D6, 'The defect norm is absurdly &
              &high. '//trim(explanation),defect_norm)
         call utils_assert(res_norm < 1D6, 'The residual norm is absurdly &
              &high. '//trim(explanation),res_norm)

         call utils_assert(defcorr_iter < pub_mg_max_iters_defco, &
              'The multigrid defect-correction solver failed to converge within&
              & the prescribed number of iterations. '//trim(explanation))

         ! hhh: Evaluate stop criteria
         if (error_norm < pub_mg_tol_pot_abs) converged = .true.

       end do
       ! **********************************************************************
       ! *** END OF DEFECT CORRECTION LOOP
       ! **********************************************************************

       if(spokesman) then
          write(stdout,'(a)') 'MG Defect correct solution obtained.'
          write(stdout,'(a,e12.4)') 'MG  - final defect norm:            ', &
               defect_norm
          write(stdout,'(a,e12.4)') 'MG  - final error norm:             ', &
               error_norm
          write(stdout,'(a,i0)')    'MG  - defect-correction iterations: ', &
               defcorr_iter
          if(pub_mg_use_error_damping) then
             write(stdout,'(a,i0)') 'MG  - error-damping iterations:     ', &
                  damping_iterations_total
          end if
       end if

       ! jd: Clean up
       ! JCW: Important to clean up allocated pointers, since these
       ! JCW: are not automatically deallocated when they go out of
       ! JCW: scope (eps_full_loc, eps_half_loc)
       if(.not. present(eps_full)) then
          ! JCW: If allocated, deallocate memory referenced by pointer
          deallocate(eps_full_loc,stat=ierr)
          call utils_dealloc_check(myself,'=>eps_full_loc',ierr)
       else
          ! JCW: If not allocated, nullify pointer
          nullify(eps_full_loc)
       end if
       if(.not. present(eps_half)) then
          ! JCW: If allocated, deallocate memory referenced by pointer
          deallocate(eps_half_loc,stat=ierr)
          call utils_dealloc_check(myself,'=>eps_half_loc',ierr)
       else
          ! JCW: If not allocated, nullify pointer
          nullify(eps_half_loc)
       end if
       deallocate(error,stat=ierr)
       call utils_dealloc_check(myself,'error',ierr)
       deallocate(defect_bt,stat=ierr)
       call utils_dealloc_check(myself,'defect_bt',ierr)
       deallocate(defect_stpt,stat=ierr)
       call utils_dealloc_check(myself,'defect_stpt',ierr)
       deallocate(defect,stat=ierr)
       call utils_dealloc_check(myself,'defect',ierr)

    end if ! if defect correction used
    ! ************************************************************************
    ! *** END OF DEFECT CORRECTION LOOP
    ! ************************************************************************

    if (multigrid_debug) then
       pot_norm_after  = sum(pot(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq)**2)
       call comms_reduce('SUM',pot_norm_after)
       pot_norm_after  = sqrt(pot_norm_after)
       if(spokesman) then
          write(stdout,'(a,e23.15)') "MG DEBUG: norm pot after defect correction = ", pot_norm_after
          call utils_flush(stdout,.true.)
       end if
    end if

    call mg_debug_dump_output(grid, cell, pot, quit = .true., res = res_local)

    deallocate(res_local,stat=ierr)
    call utils_dealloc_check(myself,'res_local', ierr)

    call utils_trace_out(myself)
    call timer_clock(myself,2)

    return

  contains

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine internal_error_damping(error, pot, defect_norm)
      ! -----------------------------------------------------------------------
      ! Error damping facility for full PBE, by Lucian Anton
      ! Modified for ONETEP conventions by Jacek Dziedzic
      ! -----------------------------------------------------------------------

      implicit none

      real(kind=DP), intent(inout) :: error(:,:,:)
      real(kind=DP), intent(in)    :: pot(:,:,:)
      real(kind=DP)                :: defect_norm

      ! jd: Internal variables
      real(kind=DP), parameter     :: s_min = 1.0e-5_DP ! minimum damping
      real(kind=DP), parameter     :: q_damp =0.75_DP   ! damping factor

      real(kind=DP)                :: s
      real(kind=DP)                :: defect_norm_aux
      real(kind=DP)                :: defect_norm_aux_stpt
      real(kind=DP)                :: defect_norm_aux_bt
      logical                      :: found
      real(kind=DP), allocatable   :: pot_aux(:,:,:), defect_aux(:,:,:,:)
      integer                      :: ierr
      character(len=*), parameter  :: myself = 'internal_error_damping'

      ! ----------------------------------------------------------------------

      s = 1.0_DP
      damping_iterations = 0
      found = .false.

      allocate(pot_aux(mg%ld1f,mg%ld2f,mg%max_slabs12), stat=ierr)
      call utils_alloc_check(myself,'pot_aux',ierr)
      allocate(defect_aux(mg%ld1f,mg%ld2f,mg%max_slabs12,3), stat=ierr)
      call utils_alloc_check(myself,'defect_aux',ierr)

      do while ( s > s_min )

         pot_aux = pot + s * error

         call compute_current_defect(pub_mg_defco_fd_order, rho, &
              eps_full_loc, pot_aux, grid, &
              defect_aux(:,:,:,1), defect_norm_aux, &
              defect_aux(:,:,:,2), defect_norm_aux_stpt, &
              defect_aux(:,:,:,3), defect_norm_aux_bt)

         damping_iterations = damping_iterations + 1

         if(spokesman) then
            write(stdout,'(a,i4,a)',advance='no') "MG", defcorr_iter,"  "
            write(stdout,'(e13.7,a,e13.7,a,e13.7,a,i2,f11.8,a)') &
                 defect_norm_aux_stpt," ",defect_norm_aux_bt," ",&
                 defect_norm_aux,"        DAMPING (",damping_iterations,s,")"
         end if

         if (defect_norm_aux < defect_norm) then
            error = s * error
            found = .true.
            exit
         else
            s = q_damp * s
         endif

      end do

      call utils_assert(found, &
           'Multigrid defect-correction error damping failed.')

      deallocate(defect_aux, stat=ierr)
      call utils_dealloc_check(myself,'defect_aux',ierr)
      deallocate(pot_aux, stat=ierr)
      call utils_dealloc_check(myself,'pot_aux',ierr)

    end subroutine internal_error_damping
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine internal_multigrid_verbose_input
      ! -----------------------------------------------------------------------
      ! -----------------------------------------------------------------------

      use comms, only: pub_my_proc_id
      use is_solvation_boltzmann, only: steric_pot, gamma
      use rundat, only: pub_is_pbe, pub_is_multigrid_verbose, pub_rootname, &
           pub_is_multigrid_verbose_y, pub_is_multigrid_verbose_z
      use utils, only: utils_abort, utils_unit

      implicit none

      character(len=1024) :: dump_in_filename
      integer             :: dump_in_unit

      ! ----------------------------------------------------------------------

      ! jd: Write out multigrid input verbose info file
      if(pub_is_multigrid_verbose) then

         ! jd: Determine grid indices nearest to y, z, and proc on which z lives.
         i3 = int(pub_is_multigrid_verbose_z/mg%d3f) + 1
         i2 = int(pub_is_multigrid_verbose_y/mg%d2f) + 1
         verbose_actual_y = real((i2-1),kind=DP) * mg%d2f
         verbose_actual_z = real((i3-1),kind=DP) * mg%d3f
         if(i3<1) i3 = 1
         if(i2<1) i2 = 1
         if(i2 > grid%n2) i2 = grid%n2
         if(i3 > grid%n3) i3 = grid%n3
         owner_proc = grid%proc_slab12(i3)

         if(pub_my_proc_id == owner_proc) then

            write(dump_in_filename,'(a,i0,a)') &
                 trim(pub_rootname)//'_mg_input_', total_calls, '.txt'
            dump_in_unit = utils_unit()

            ! jd: If first call, delete file first.
            if(total_calls == 1) then
               open(unit=dump_in_unit, file=dump_in_filename, status="old", &
                    err=500)
500            continue ! no-op if file was not there to begin with
               close(unit=dump_in_unit, status="delete", err=500)
            end if

            ! jd: Open for append
            open(unit = dump_in_unit, file = dump_in_filename, &
                 action = "write", err = 1100, position = "append")

            ! jd: Write header.
            write(dump_in_unit,'(a,f10.3,a,f10.3,a,i0,a,i0,a,i0)') &
                 '# y: ', pub_is_multigrid_verbose_y,', z: ', &
                 pub_is_multigrid_verbose_z,' pq: ', mg%pq1f, ' pt: ', &
                 mg%pt1f, ' ld: ', mg%ld1f
            write(dump_in_unit,'(a)') &
                 '#   i1          x          rho      bound   eps_full &
                 &eps_half_x eps_half_y &
                 &eps_half_z steric_pot   gamma'

            islab12 = i3 - mg%my_first_slab12 + 1
            do i1 = 1, mg%ld1f
               write(dump_in_unit,'(i6,a,f10.3,a,e12.3,a)',advance='no') &
                    i1,' ',real((i1-1),kind=DP) * mg%d1f,' ',&
                    merge(&
                    rho(i1,i2,islab12),0D0,abs(rho(i1,i2,islab12))>1D-99&
                    ), &
                    ' '
               if(present(bound)) then
                  write(dump_in_unit,'(e10.3,a1)',advance='no') &
                       merge(&
                       bound(i1,i2,islab12),0D0,abs(bound(i1,i2,islab12))>1D-99&
                       ), &
                       ' '
               else
                  write(dump_in_unit,'(e10.3,a1)',advance='no') &
                       0D0,' '
               end if
               if(present(eps_full)) then
                  write(dump_in_unit,'(e10.3,a1)',advance='no') &
                       merge(&
                       eps_full(i1,i2,islab12),0D0,&
                       abs(eps_full(i1,i2,islab12))>1D-99&
                       ), &
                       ' '
               else
                  write(dump_in_unit,'(e10.3,a1)',advance='no') &
                       0D0,' '
               end if
               if(present(eps_half)) then
                  write(dump_in_unit,'(e10.3,a1,e10.3,a1,e10.3,a1)', &
                       advance='no') &
                       merge(&
                       eps_half(i1,i2,islab12,1),0D0,&
                       abs(eps_half(i1,i2,islab12,1))>1D-99&
                       ), &
                       ' ', &
                       merge(&
                       eps_half(i1,i2,islab12,2),0D0,&
                       abs(eps_half(i1,i2,islab12,2))>1D-99&
                       ), &
                       ' ', &
                       merge(&
                       eps_half(i1,i2,islab12,3),0D0,&
                       abs(eps_half(i1,i2,islab12,3))>1D-99&
                       ), &
                       ' '
               else
                  write(dump_in_unit,'(e10.3,a1,e10.3,a1,e10.3,a1)', &
                       advance='no') 0D0,' ',0D0,' ',0D0,' '
               end if
               if(pub_is_pbe /= 'NONE') then
                  write(dump_in_unit,'(e10.3,a1,e10.3,a1)',advance='no') &
                       merge(&
                       steric_pot(i1,i2,islab12),0D0,&
                       abs(steric_pot(i1,i2,islab12))>1D-99&
                       ), &
                       ' ', &
                       merge(&
                       gamma(i1,i2,islab12),0D0,&
                       abs(gamma(i1,i2,islab12))>1D-99&
                       ), &
                       ' '
               else
                  write(dump_in_unit,'(e10.3,a1,e10.3,a1)',advance='no') &
                       0D0,' ',0D0,' '
               end if

               write(dump_in_unit,'(a)') ''
            end do

            close(dump_in_unit, err = 2100)
         end if ! proc responsible for output
      end if ! multgrid_verbose

      return

1100  call utils_abort('Error during creation of file: '//dump_in_filename)
2100  call utils_abort('Error during closing of file: '//dump_in_filename)

    end subroutine internal_multigrid_verbose_input
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine internal_multigrid_verbose_input_defco
      ! -----------------------------------------------------------------------
      ! -----------------------------------------------------------------------

      use comms, only: pub_my_proc_id
      use is_solvation_boltzmann, only: steric_pot, gamma
      use rundat, only: pub_is_pbe, pub_is_multigrid_verbose, pub_rootname, &
           pub_is_multigrid_verbose_y, pub_is_multigrid_verbose_z
      use utils, only: utils_abort, utils_unit

      implicit none

      character(len=1024) :: dump_in_filename
      integer             :: dump_in_unit

      ! ----------------------------------------------------------------------

      ! jd: Write out multigrid input verbose info file
      if(pub_is_multigrid_verbose) then

         ! jd: Determine grid indices nearest to y, z, and proc on which z lives.
         i3 = int(pub_is_multigrid_verbose_z/mg%d3f) + 1
         i2 = int(pub_is_multigrid_verbose_y/mg%d2f) + 1
         verbose_actual_y = real((i2-1),kind=DP) * mg%d2f
         verbose_actual_z = real((i3-1),kind=DP) * mg%d3f
         if(i3<1) i3 = 1
         if(i2<1) i2 = 1
         if(i2 > grid%n2) i2 = grid%n2
         if(i3 > grid%n3) i3 = grid%n3
         owner_proc = grid%proc_slab12(i3)

         if(pub_my_proc_id == owner_proc) then

            write(dump_in_filename,'(a,i0,a,i0,a)') &
                 trim(pub_rootname)//'_mg_input_', total_calls, '_', &
                 defcorr_iter, '.txt'
            dump_in_unit = utils_unit()

            ! jd: If first call, delete file first.
            if(total_calls == 1) then
               open(unit=dump_in_unit, file=dump_in_filename, status="old", &
                    err=500)
500            continue ! no-op if file was not there to begin with
               close(unit=dump_in_unit, status="delete", err=500)
            end if

            ! jd: Open for append
            open(unit = dump_in_unit, file = dump_in_filename, &
                 action = "write", err = 3100, position = "append")

            ! jd: Write header.
            write(dump_in_unit,&
                 '(a,f10.3,a,f10.3,a,i0,a,i0,a,i0)')&
                 '# y: ', pub_is_multigrid_verbose_y,', z: ', &
                 pub_is_multigrid_verbose_z,' pq: ', mg%pq1f, ' pt: ', &
                 mg%pt1f, ' ld: ', mg%ld1f
            write(dump_in_unit,'(a)') &
                 '#   i1          x       defect  defect_stpt    &
                 &defect_bt      bound   eps_full &
                 &eps_half_x eps_half_y &
                 &eps_half_z steric_pot   gamma       der_pot'

            islab12 = i3 - mg%my_first_slab12 + 1
            do i1 = 1, mg%ld1f
               write(dump_in_unit,'(i6,a,f10.3,a,e12.3,a)',advance='no') &
                    i1,' ',real((i1-1),kind=DP) * mg%d1f,' ',&
                    merge(&
                    defect(i1,i2,islab12),0D0,abs(defect(i1,i2,islab12))>1D-99&
                    ), &
                    ' '
               write(dump_in_unit,'(e12.3,a)',advance='no') &
                    merge(&
                    defect_stpt(i1,i2,islab12),0D0,&
                    abs(defect_stpt(i1,i2,islab12))>1D-99&
                    ), &
                    ' '
               write(dump_in_unit,'(e12.3,a)',advance='no') &
                    merge(&
                    defect_bt(i1,i2,islab12),0D0,&
                    abs(defect_bt(i1,i2,islab12))>1D-99&
                    ), &
                    ' '
               if(present(bound)) then
                  write(dump_in_unit,'(e10.3,a1)',advance='no') &
                       merge(&
                       bound(i1,i2,islab12),0D0,abs(bound(i1,i2,islab12))>1D-99&
                       ), &
                       ' '
               else
                  write(dump_in_unit,'(e10.3,a1)',advance='no') &
                       0D0,' '
               end if
               if(present(eps_full)) then
                  write(dump_in_unit,'(e10.3,a1)',advance='no') &
                       merge(&
                       eps_full(i1,i2,islab12),0D0,&
                       abs(eps_full(i1,i2,islab12))>1D-99&
                       ), &
                       ' '
               else
                  write(dump_in_unit,'(e10.3,a1)',advance='no') &
                       0D0,' '
               end if
               if(present(eps_half)) then
                  write(dump_in_unit,'(e10.3,a1,e10.3,a1,e10.3,a1)', &
                       advance='no') &
                       merge(&
                       eps_half(i1,i2,islab12,1),0D0,&
                       abs(eps_half(i1,i2,islab12,1))>1D-99&
                       ), &
                       ' ', &
                       merge(&
                       eps_half(i1,i2,islab12,2),0D0,&
                       abs(eps_half(i1,i2,islab12,2))>1D-99&
                       ), &
                       ' ', &
                       merge(&
                       eps_half(i1,i2,islab12,3),0D0,&
                       abs(eps_half(i1,i2,islab12,3))>1D-99&
                       ), &
                       ' '
               else
                  write(dump_in_unit,'(e10.3,a1,e10.3,a1,e10.3,a1)', &
                       advance='no') 0D0,' ',0D0,' ',0D0,' '
               end if
               if(pub_is_pbe /= 'NONE') then
                  write(dump_in_unit,'(e10.3,a1,e10.3,a1,e10.3,a1)', &
                       advance='no') &
                       merge(&
                       steric_pot(i1,i2,islab12),0D0,&
                       abs(steric_pot(i1,i2,islab12))>1D-99&
                       ), &
                       ' ', &
                       merge(&
                       gamma(i1,i2,islab12),0D0,&
                       abs(gamma(i1,i2,islab12))>1D-99&
                       ), &
                       ' ', &
                       merge(&  ! 'pot' corresponds to argument der_pot
                       pot(i1,i2,islab12),0D0,&
                       abs(pot(i1,i2,islab12))>1D-99&
                       ), &
                       ' '
               else
                  write(dump_in_unit,'(e10.3,a1,e10.3,a1,e10.3,a1)',&
                       advance='no') &
                       0D0,' ',0D0,' ',0D0,' '
               end if

               write(dump_in_unit,'(a)') ''
            end do

            close(dump_in_unit, err = 4100)
         end if ! proc responsible for output
      end if ! multgrid_verbose

      return

3100  call utils_abort('Error during creation of file: '//dump_in_filename)
4100  call utils_abort('Error during closing of file: '//dump_in_filename)

    end subroutine internal_multigrid_verbose_input_defco
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  end subroutine multigrid_defect_corr_solver
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine second_order_multigrid_solver(pot, residual, rho, grid, cell, &
       eps_full, eps_half, bound, tolerance_override, use_linear, der_pot)
    !=========================================================================!
    ! For a given electronic density (rho) and dielectric functional (eps)    !
    ! this subroutine will compute the electrostatic potential (pot) which is !
    ! the solution to the Poisson equation:                                   !
    !                                                                         !
    !        -  \nabla \cdot [ eps \nabla (pot) ] = 4*pi*rho                  !
    !                                                                         !
    ! ... or the Poisson-Boltzmann equation if is_pbe is not 'NONE'.          !
    !                                                                         !
    ! This is the middleware subroutine that calls the DL_MG solver once, and !
    ! only takes care to prepare/fake optional values (eps if not given, pot  !
    ! from bound), do sanity checks and compute the residual.                 !
    !                                                                         !
    ! The obtained solution is to 2nd order. Use multigrid_defect_corr_solver !
    ! as the user-level subroutine.                                           !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   pot,               (output),     the potential we seek                !
    !                                                                         !
    !   residual,          (output)      useful for error analysis            !
    !                                                                         !
    !   rho,               (input),      the input electron density           !
    !                                                                         !
    !   grid,              (input),      while the solver operates on the mg  !
    !                                    grid, the 'grid' argument is needed  !
    !                                    in an integrals_mod call.            !
    !   eps_full           (input),      the dielectric functional on the     !
    !                                    grid                                 !
    !                                                                         !
    !   eps_half,          (input),      the dielectric functional            !
    !                                    interpolated to 'half' grids         !
    !                                                                         !
    !   bound,             (input)       initial guess for the potential      !
    !                      optional      and fixed boundary condition         !
    !                      Defaults to 0.0 everywhere if omitted              !
    !                                                                         !
    !   tolerance_override (input)       tolerance for stop criterion in the  !
    !                      optional      mg solver.                           !
    !                      Defaults to pub_mg_tol_pot_abs if omitted.         !
    !                                                                         !
    !   use_linear         (input, opt)  passed to the solver in PBE case.    !
    !   der_pot            (input, opt)  passed to the solver in PBE case.    !
    !                                                                         !
    ! All the arrays are assumed to be in the distributed slab representation.!
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal, 19/10/2007                                    !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    ! Updated for defect correction by HH Helal, 06/2010                      !
    ! Adapted for fully parallel operation by Jacek Dziedzic, 03/2011.        !
    ! Modified many times over 2011-2014 by Jacek Dziedzic.                   !
    ! Updated by James C. Womack, 02/2017, to reflect new DL_MG interface.    !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_barrier, pub_on_root
    use constants, only: FOUR_PI, CRLF, garbage_real
    use geometry, only: operator(.DOT.)
    use integrals, only: integrals_product_on_grid
    use is_solvation_boltzmann, only: boltzmann_ions, n_boltzmann_ions, gamma, &
         max_expcap_dl_mg, implicit_solvent_boltzmann_mg_neutralisation_method,&
         IS_PBE_MG_RETURN
    use rundat, only: pub_mg_tol_pot_abs, pub_is_pbe, &
         pub_is_pbe_exp_cap, pub_is_pbe_temperature, pub_mg_pbe_use_fas, &
         pub_mg_vcyc_smoother_iter_pre, pub_mg_vcyc_smoother_iter_post, &
         pub_mg_max_res_ratio, pub_mg_tol_mu_rel, pub_mg_tol_mu_abs, &
         pub_mg_tol_cg_res_rel, pub_mg_tol_cg_res_abs, &
         pub_mg_continue_on_error, pub_mg_max_iters_cg, pub_mg_use_cg
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_assert, utils_sanity_check, &
         utils_trace_in, utils_trace_out, utils_feature_not_supported

    implicit none

    ! --------------------------------------------------------------------------
    ! Output arguments
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12), intent(out) &
         :: pot
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12), intent(out) &
         :: residual

    ! Input arguments
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12), & ! [*]
         intent(in)                                                 :: rho
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12), &
         intent(in)                                                 :: eps_full
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12, 3), &
         intent(in)                                                 :: eps_half
    ! jd: The next two arrays are dimensioned like the first one @[*].
    !     For some reason declaring the optionals as assumed-shape works
    !     around ifort's bug that manifests in debug mode as the compiler
    !     thinking the optional der_pot is suddenly present in dl_mg_solver_pbe
    !     even though it is not passed (cf. ONETEP bug tracker #1553).
    !     Workaround due to jme. Bug present at least in ifort v16.0.
    real(kind=DP), dimension(:,:,:), intent(in), optional           :: bound
    real(kind=DP), dimension(:,:,:), intent(in), optional           :: der_pot
    type(GRID_INFO), intent(in)         :: grid
    type(CELL_INFO), intent(in)         :: cell
    real(kind=DP), intent(in), optional :: tolerance_override
    logical, intent(in), optional       :: use_linear

    ! --------------------------------------------------------------------------

    ! Internal variables
    real(kind=DP) :: res_norm   ! jd: Norm of the residual for the defect-corr'n
    integer :: ierr             ! jd: Error flag
    real(kind=DP) :: tolerance  ! jd: Requested multigrid tolerance
    logical :: use_linear_loc   ! JCW: Allows the value of pub_is_pbe to be
                                ! JCW: overridden
    real(kind=DP) :: exp_cap

    type(IS_PBE_MG_RETURN) :: mg_return

#ifdef HAVE_DL_MG
    character(len=DL_MG_MAX_ERROR_STRING) :: err_msg
    integer :: err_msg_length
#endif

    integer :: dlmg_nwt_iteration, dlmg_defco_iteration
    integer :: dlmg_nwt_ierror, dlmg_defco_ierror
    real(kind=DP) :: dlmg_nwt_solution_norm, dlmg_defco_solution_norm
    real(kind=DP) :: dlmg_nwt_residual_norm, dlmg_defco_residual_norm
    real(kind=DP) :: dlmg_nwt_target_residual_norm, dlmg_defco_target_residual_norm
    real(kind=DP) :: dlmg_defco_error_norm

    character(len=*), parameter :: myself = 'second_order_multigrid_solver'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier()    ! jd: Synchronize to get accurate timings for solver
    call timer_clock(myself,1)

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         trim(myself)//': Must call multigrid_initialise() first.')

    ! jd: Sanity checks on input
    call utils_sanity_check(rho,'rho in '//myself)
    call utils_sanity_check(eps_full, 'eps_full in '//myself)
    call utils_sanity_check(eps_half(:,:,:,1),'eps_half/x in '//myself)
    call utils_sanity_check(eps_half(:,:,:,2),'eps_half/y in '//myself)
    call utils_sanity_check(eps_half(:,:,:,3),'eps_half/z in '//myself)
    if(present(bound)) call utils_sanity_check(bound,'bound in '//myself)

    ! hhh: We have to be careful here as we implicitly assume our cell is
    !      orthorhombic for using the multigrid solver
    ! jd:  Explicitly check if we are safe
    call utils_assert(&
         (cell%a1_unit .dot. cell%a2_unit) == 0.0_DP .and. &
         (cell%a1_unit .dot. cell%a3_unit) == 0.0_DP .and. &
         (cell%a2_unit .dot. cell%a3_unit) == 0.0_DP, &
         'The multigrid solver only supports orthorhombic cells, sorry.')

    ! jd: If boundary conditions are specified in 'bound', use them,
    !     else use zero BC. Since bound is intent(in) and the solver expects
    !     the BCs in the pot array, let's use the pot array.
    !     First zero pot to take care of padding.
    pot = 0.0_DP
    if(present(bound)) then
       pot(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq) = &
            bound(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq)
    end if

    ! jd: Allow the number of iterations to be overridden, so that crude
    !     solutions are possible from within the defect_corr_solver
    if(present(tolerance_override)) then
       tolerance = tolerance_override
    else
       tolerance = pub_mg_tol_pot_abs
    end if

    if ( present(use_linear) ) then
       use_linear_loc = use_linear
    else
       use_linear_loc = .false.
    end if

    ! jd: Take care of padding in residual
    residual = 0.0_DP

    ! **********************************************************************
    ! jd: Call DL_MG
    ! **********************************************************************
#ifdef HAVE_DL_MG
    if(pub_is_pbe_exp_cap == 0.0_DP) then
       exp_cap = max_expcap_dl_mg
    else
       exp_cap = pub_is_pbe_exp_cap
    end if

    ! jd: This is needed to later determine if we actually used counterion
    !     neutralisation. If not, the ratios will stay set to garbage_real.
    mg_return%used_neutralisation_ratios(1:n_boltzmann_ions) = garbage_real

    call timer_clock('dl_mg_solver_calls',1)
    ! JCW: New DL_MG interface exposes one name for dl_mg_solver
    ! JCW: and distinguishes whether to solve Poisson or Poisson-Boltzmann
    ! JCW: equation via arguments supplied.
    if(pub_is_pbe == 'FULL') then
       ! **********************************************************************
       ! *** Poisson-Boltzmann: FULL or override to LINEARISED (use_linear_loc)
       ! **********************************************************************
       call dl_mg_solver(eps_full, eps_half,-FOUR_PI, rho,&
            -FOUR_PI, pub_is_pbe_temperature, n_boltzmann_ions, &
            boltzmann_ions(1:n_boltzmann_ions)%concentration_user, &
            boltzmann_ions(1:n_boltzmann_ions)%charge, &
            pot, 2,&
            tol_res_rel = 0.0_dp, &
            tol_res_abs = tolerance, &
            tol_mu_rel = pub_mg_tol_mu_rel, &
            tol_mu_abs = pub_mg_tol_mu_abs, &
            tol_cg_res_rel = pub_mg_tol_cg_res_rel, &
            tol_cg_res_abs = pub_mg_tol_cg_res_abs, &
            max_iters_cg = pub_mg_max_iters_cg, &
            use_cg = pub_mg_use_cg, &

            linearised=use_linear_loc, der_pot=der_pot, &
            steric_weight=gamma, expcap=exp_cap, &
            use_fas=pub_mg_pbe_use_fas, &

            neutralisation_method = &
            implicit_solvent_boltzmann_mg_neutralisation_method(mg,grid,rho), &
            neutralisation_ion_ratios = &
            boltzmann_ions(:)%neutralisation_ratio, &

            betamu_electrostatic = &
            mg_return%betamu_elec(1:n_boltzmann_ions), &
            used_ion_concentrations = &
            mg_return%used_ion_concentration(1:n_boltzmann_ions), &
            used_neutralisation_ratios = &
            mg_return%used_neutralisation_ratios(1:n_boltzmann_ions), &
            steric_weight_average = mg_return%steric_weight_average, &

            v_iterations = (/pub_mg_vcyc_smoother_iter_pre, &
            pub_mg_vcyc_smoother_iter_post/), &
            mg_max_conv_rate = pub_mg_max_res_ratio, &
            res=residual,ierror=ierr)
    elseif(pub_is_pbe == 'LINEARISED') then
       ! *************************************
       ! *** Poisson-Boltzmann: LINEARISED ***
       ! *************************************
       call dl_mg_solver(eps_full, eps_half,-FOUR_PI, rho,&
            -FOUR_PI, pub_is_pbe_temperature, n_boltzmann_ions, &
            boltzmann_ions(1:n_boltzmann_ions)%concentration_user, &
            boltzmann_ions(1:n_boltzmann_ions)%charge, &
            pot, 2, &
            tol_res_rel = 0.0_dp, &
            tol_res_abs = tolerance, &
            tol_mu_rel = pub_mg_tol_mu_rel, &
            tol_mu_abs = pub_mg_tol_mu_abs, &
            tol_cg_res_rel = pub_mg_tol_cg_res_rel, &
            tol_cg_res_abs = pub_mg_tol_cg_res_abs, &
            max_iters_cg = pub_mg_max_iters_cg, &
            use_cg = pub_mg_use_cg, &
            linearised=.true., &
            steric_weight=gamma,expcap=exp_cap, &
            use_fas=pub_mg_pbe_use_fas, &

            neutralisation_method = &
            implicit_solvent_boltzmann_mg_neutralisation_method(mg,grid,rho), &
            neutralisation_ion_ratios = &
            boltzmann_ions(:)%neutralisation_ratio, &

            betamu_electrostatic = &
            mg_return%betamu_elec(1:n_boltzmann_ions), &
            used_ion_concentrations = &
            mg_return%used_ion_concentration(1:n_boltzmann_ions), &
            used_neutralisation_ratios = &
            mg_return%used_neutralisation_ratios(1:n_boltzmann_ions), &
            steric_weight_average = mg_return%steric_weight_average, &

            v_iterations = (/pub_mg_vcyc_smoother_iter_pre, &
            pub_mg_vcyc_smoother_iter_post/), &
            mg_max_conv_rate = pub_mg_max_res_ratio, &
            res=residual,ierror=ierr)

    elseif(pub_is_pbe == 'NONE') then
       ! ********************
       ! *** Pure Poisson ***
       ! ********************
       call dl_mg_solver(eps_full, eps_half,-FOUR_PI, rho,&
            pot, 2, &
            tol_res_rel = 0.0_dp, &
            tol_res_abs = tolerance, &
            tol_cg_res_rel = pub_mg_tol_cg_res_rel, &
            tol_cg_res_abs = pub_mg_tol_cg_res_abs, &
            max_iters_cg = pub_mg_max_iters_cg, &
            use_cg = pub_mg_use_cg, &
            v_iterations = (/pub_mg_vcyc_smoother_iter_pre, &
            pub_mg_vcyc_smoother_iter_post/), &
            mg_max_conv_rate = pub_mg_max_res_ratio, &
            res=residual, ierror=ierr)
    else
       call utils_abort('Internal error. Unrecognized pub_is_pbe')
    endif
    call timer_clock('dl_mg_solver_calls',2)

    if(ierr /= 0) then
       call dl_mg_error_string(ierr,err_msg,err_msg_length)

       if(pub_mg_continue_on_error) then
          if(pub_on_root) then
             write(stdout,'(a)') "MG DL_MG: ------------------------- &
                  &SERIOUS WARNING: -------------------------"
             write(stdout,'(a)') "MG DL_MG: Could not achieve convergence to&
                  & prescribed tolerance settings!"
             write(stdout,'(a)') "MG "//trim(err_msg)//'.'
             call dl_mg_solver_status(dlmg_nwt_iteration, &
                  dlmg_nwt_solution_norm, dlmg_nwt_residual_norm, &
                  dlmg_nwt_target_residual_norm, dlmg_nwt_ierror, &
                  dlmg_defco_iteration, dlmg_defco_solution_norm, &
                  dlmg_defco_residual_norm, dlmg_defco_error_norm, &
                  dlmg_defco_target_residual_norm, dlmg_defco_ierror)

             write(stdout,'(a,i3,a,i6,a,e10.3,a)') &
                  "MG DL_MG: Newton: error code:", dlmg_nwt_ierror, &
                  ", n_iter:", dlmg_nwt_iteration, ", sol_norm:", &
                  dlmg_nwt_solution_norm, "."
             write(stdout,'(a,e10.3,a,e10.3,a)') &
                  "MG DL_MG: Newton: res_norm:", dlmg_nwt_residual_norm, &
                  ", target_res_norm:", dlmg_nwt_target_residual_norm, &
                  "."
             write(stdout,'(a,i3,a,i6,a,e10.3,a)') &
                  "MG DL_MG: DEFCO:  error code:", dlmg_defco_ierror, &
                  ", n_iter:", dlmg_defco_iteration, ", sol_norm:", &
                  dlmg_defco_solution_norm, "."
             write(stdout,'(a,e10.3,a,e10.3,a)') &
                  "MG DL_MG: DEFCO:  res_norm:", dlmg_defco_residual_norm, &
                  ", target_res_norm:", dlmg_defco_target_residual_norm, "."
             write(stdout,'(a,e10.3,a)') &
                  "MG DL_MG: DEFCO:  err_norm:", dlmg_defco_error_norm, "."
             write(stdout,'(a)') "MG DL_MG: --------------------------&
                  &------------------------------------------"

          end if
       else
          call utils_abort("Error in the DL_MG multigrid solver: "//CRLF//&
               trim(err_msg))
       end if
    end if
#else
    call utils_feature_not_supported('DL_MG multigrid solver')
#endif
    ! **********************************************************************
    ! **********************************************************************

    ! jd: Calculate the residual norm
    call utils_sanity_check(residual,'residual in '//myself)
    res_norm = sqrt(integrals_product_on_grid(grid,&
         residual,residual,mg%pq1f,mg%pq2f,mg%num_slabs12_pq))

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine second_order_multigrid_solver
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine compute_current_defect(order, rho, eps, v_i, grid, & ! in
    defect, defect_norm, defect_stpt, defect_norm_stpt, &         ! out
    defect_bt, defect_norm_bt)                                    ! out
    !=========================================================================!
    ! @updateme
    ! This subroutine will compute the current defect for the Poisson equation!
    ! in a dielectric medium.  The defect is defined as the amount that the   !
    ! current approximation fails to satisfy the equation or in other words:  !
    !                                                                         !
    !    defect = source + \nabla \cdot [eps \grad v_i],                      !
    ! where source = 4\pi \rho                                                !
    !                                                                         !
    ! The discretization of the operators is selected by the                  !
    ! discretization_order parameter.  This subroutine is intended for use    !
    ! within the defect correction method which will produce a higher         !
    ! accuracy solution than is possible with normal multigrid iteration.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !    order       (input): finite difference order to use                  !
    !    rho         (input): fixed charge density                            !
    !    eps         (input): the diel. functional on the full grid points    !
    !    v_i         (input): the current approximation of the potential      !
    !    grid        (input): the defect is computed on the mg grid, but the  !
    !                         host grid is needed here and there.             !
    !    defect      (output): defined as above                               !
    !    defect_norm (output): discrete L2 norm of the defect function        !
    !                                                                         !
    ! Arrays are assumed to be in the distributed slab representation.        !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal 30/04/2010                                     !
    ! Mended for onetep by HH Helal 7/6/2010                                  !
    ! Adapted for parallel operation by Jacek Dziedzic, 10-11/06/2010.        !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: FOUR_PI
    use finite_differences, only: finite_difference_gradient, &
         finite_difference_laplacian
    use integrals, only: integrals_product_on_grid
    use is_solvation_boltzmann, only: implicit_solvent_boltzmann_defect
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_sanity_check, utils_trace_in, utils_trace_out
    use rundat, only : pub_is_pbe

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in) :: grid
    integer,       intent(in)   :: order
    real(kind=DP), intent(in)   :: rho(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)   :: eps(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)   :: v_i(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(out)  :: defect(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(out)  :: defect_norm
    real(kind=DP), intent(out)  :: defect_stpt(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(out)  :: defect_norm_stpt
    real(kind=DP), intent(out)  :: defect_bt(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(out)  :: defect_norm_bt

    ! jd: Internal variables
    real(kind=DP), allocatable     :: lap_v(:,:,:)         ! nabla^2 v
    real(kind=DP), allocatable     :: grad_eps(:,:,:,:)    ! nabla eps
    real(kind=DP), allocatable     :: grad_v(:,:,:,:)      ! nabla v
    integer                        :: ierr                 ! jd: Error flag
    character(len=*), parameter    :: myself = 'compute_current_defect'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    call utils_sanity_check(rho,'rho in '//myself)
    call utils_sanity_check(eps,'eps in '//myself)
    call utils_sanity_check(v_i,'v_i in'//myself)

    ! Allocate memory for local variables
    allocate(lap_v(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'lap_v',ierr)
    allocate(grad_eps(mg%ld1f,mg%ld2f,mg%max_slabs12,3),stat=ierr)
    call utils_alloc_check(myself,'grad_eps',ierr)
    allocate(grad_v(mg%ld1f,mg%ld2f,mg%max_slabs12,3),stat=ierr)
    call utils_alloc_check(myself,'grad_v',ierr)

    ! Calculate the gradient and Laplacian of the potential
    call finite_difference_gradient(grad_v,v_i,order,mg,grid)
    call finite_difference_laplacian(lap_v,v_i,order,mg,grid)
    ! jd: Calculate the gradient of epsilon
    call finite_difference_gradient(grad_eps,eps,order,mg,grid)

    call utils_sanity_check(lap_v,'lap_v in '//myself)
    call utils_sanity_check(grad_v(:,:,:,1),'grad_v/x in '//myself)
    call utils_sanity_check(grad_v(:,:,:,2),'grad_v/y in '//myself)
    call utils_sanity_check(grad_v(:,:,:,3),'grad_v/z in '//myself)
    call utils_sanity_check(grad_eps(:,:,:,1),'grad_eps/x in '//myself)
    call utils_sanity_check(grad_eps(:,:,:,2),'grad_eps/y in '//myself)
    call utils_sanity_check(grad_eps(:,:,:,3),'grad_eps/z in '//myself)

    ! Now have everything we need to compute the current defect and defect_norm
    defect_stpt(:,:,:) = FOUR_PI * rho(:,:,:) + &
         eps(:,:,:) * lap_v(:,:,:) + &
         grad_eps(:,:,:,1) * grad_v(:,:,:,1) + &
         grad_eps(:,:,:,2) * grad_v(:,:,:,2) + &
         grad_eps(:,:,:,3) * grad_v(:,:,:,3)
    defect_norm_stpt = sqrt(integrals_product_on_grid(grid, &
         defect_stpt, defect_stpt, mg%pq1f, mg%pq2f, mg%num_slabs12_pq))

    if (pub_is_pbe /= 'NONE') then
       call implicit_solvent_boltzmann_defect(defect_bt, v_i, &
            mg%pq1f, mg%pq2f, mg%num_slabs12_pq)
       defect_norm_bt = sqrt(integrals_product_on_grid(grid, &
            defect_bt, defect_bt, mg%pq1f, mg%pq2f, mg%num_slabs12_pq))
       defect = defect_stpt + defect_bt
       defect_norm = sqrt(integrals_product_on_grid(grid, &
            defect, defect, mg%pq1f, mg%pq2f, mg%num_slabs12_pq))
    else
       defect_bt(:,:,:) = 0D0
       defect_norm_bt = 0D0
       defect = defect_stpt
       defect_norm = defect_norm_stpt
    end if

    ! jd: Clean up
    deallocate(grad_v,stat=ierr)
    call utils_dealloc_check(myself,'grad_v',ierr)
    deallocate(grad_eps,stat=ierr)
    call utils_dealloc_check(myself,'grad_eps',ierr)
    deallocate(lap_v,stat=ierr)
    call utils_dealloc_check(myself,'lap_v',ierr)

    call utils_trace_out(myself)

  end subroutine compute_current_defect
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine mg_debug_dump_input(grid, cell, rho, eps_full, eps_half, bound)
    !=========================================================================!
    ! Outputs 3D debug dumps of crucial multigrid input.                      !
    !-------------------------------------------------------------------------!
    ! Extracted by Jacek Dziedzic from internal_multigrid_debug_dump_input()  !
    ! in July 2019.                                                           !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: ANGSTROM
    use is_solvation_boltzmann, only: gamma
    use rundat, only: pub_is_pbe, pub_is_solvation_properties, &
         pub_inner_loop_iteration
    use simulation_cell, only: CELL_INFO
    use visual, only: visual_scalarfield

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in)         :: grid
    type(CELL_INFO), intent(in)         :: cell
    real(kind=DP), intent(in)           :: rho(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional :: bound(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional :: eps_full(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional :: eps_half(mg%ld1f,mg%ld2f,mg%max_slabs12,3)

    ! ----------------------------------------------------------------------

    ! jd: Produce debug dumps if asked to via devel code,
    !     or if we are in properties and asked to via keyword
    if((pub_is_solvation_properties .and. pub_inner_loop_iteration == -2) .or. &
         multigrid_debug) then
       call visual_scalarfield(rho, grid, cell, 'Total density (e/ang^3)', &
            '_in_rho', conversion_factor = ANGSTROM**3)
       ! converts to electrons -ve conv'n  ^
       if(pub_is_pbe /= 'NONE') then
          call visual_scalarfield(gamma,grid,cell, &
               'Accessibility (gamma)','_in_gamma')
       end if
       if(present(eps_full)) then
          call visual_scalarfield(eps_full,grid,cell, &
               'Permittivity','_in_epsfull')
       end if
       if(present(eps_half)) then
          call visual_scalarfield(eps_half(:,:,:,1),grid,cell, &
               'x-shifted permittivity', '_in_epshalfx')
          call visual_scalarfield(eps_half(:,:,:,2),grid,cell, &
               'y-shifted permittivity', '_in_epshalfy')
          call visual_scalarfield(eps_half(:,:,:,3),grid,cell, &
               'z-shifted permittivity', '_in_epshalfz')
       end if
       if(present(bound)) then
          call visual_scalarfield(bound,grid,cell, &
               'Boundary conditions','_in_bound')
       end if
    end if

  end subroutine mg_debug_dump_input
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine mg_debug_dump_output(grid, cell, pot, quit, mg_return, res)
    !=========================================================================!
    ! Outputs 3D debug dumps of crucial multigrid output.                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   grid, cell (in): ONETEP grid and cell, needed by visual_scalarfield().!
    !   pot (in): Potential to output.                                        !
    !   quit (in): If set to .true., utils_abort() will be called after the   !
    !              dump is produced.                                          !
    !   mg_return (in, opt): Details of DL_MG's solution -- needed in         !
    !                        Boltzmann solvation, omitted in pure Poisson.    !
    !   res (in, opt): Final equation residual -- used only in the old scheme,!
    !                  where defect correction is done within ONETEP.         !
    !-------------------------------------------------------------------------!
    ! Extracted by Jacek Dziedzic from internal_multigrid_debug_dump_output() !
    ! in July 2019.                                                           !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use is_solvation_boltzmann, only: IS_PBE_MG_RETURN, &
         implicit_solvent_boltzmann_debug_dump
    use rundat, only: pub_is_pbe, pub_is_solvation_properties, &
         pub_inner_loop_iteration
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_assert
    use visual, only: visual_scalarfield
    use constants, only: HARTREE_IN_EVS

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in)         :: grid
    type(CELL_INFO), intent(in)         :: cell
    real(kind=DP), intent(in)           :: pot(mg%ld1f,mg%ld2f,mg%max_slabs12)
    logical, intent(in)                 :: quit
    type(IS_PBE_MG_RETURN), intent(in), optional :: mg_return
    real(kind=DP), intent(in), optional :: res(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! ----------------------------------------------------------------------

    ! jd: Produce debug dumps if asked to via devel code,
    !     or if we are in properties and asked to via keyword
    if((pub_is_solvation_properties .and. pub_inner_loop_iteration == -2) .or. &
         multigrid_debug) then
       call visual_scalarfield(pot, grid, cell, &
            'Electrostatic potential (V)','_out_pot', &
            conversion_factor = -HARTREE_IN_EVS)



       if(present(res)) then
          call visual_scalarfield(res,grid,cell, &
               'Final error equation residual','_out_resid')
       end if

       ! jd: Debugs of Boltzmann part
       if(pub_is_pbe /= "NONE") then
          ! jd: betamu_elec and steric_weight_average are inaccessible in the
          !     obsoleted internal defect-correction procedure
          if(present(mg_return)) then
             call implicit_solvent_boltzmann_debug_dump(pot, mg_return, grid, &
                  mg, cell)
          end if
       end if

       ! jd: Quit after producing output, but only if invoked from the devcode.
       if(multigrid_debug .and. quit) then
          call utils_assert(.false.,'Intentionally terminating ONETEP after MG &
               &debug dumps have been produced.')
       end if

    end if

  end subroutine mg_debug_dump_output
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine multigrid_polarization_data(grid, cell, elements, phi, eps_full)
    !========================================================================!
    ! Produces polarization density, polarization potential, polarization    !
    ! energy and smeared ion polarization correction in solvent              !
    ! Reference: Andreussi (2012)                                            !
    !========================================================================!
    ! Written by Arihant Bhandari on 13 Oct, 2020                            !
    !========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: stdout, VERBOSE, ANGSTROM, HARTREE_IN_EVS
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_multigrid_bc_is_periodic, pub_output_detail, &
         pub_mg_defco_fd_order, pub_is_bc_allow_frac_charge, pub_is_pbe
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_trace_in, &
         utils_trace_out
    use visual, only: visual_scalarfield
    use ion, only: ELEMENT

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)     :: grid
    type(CELL_INFO), intent(in)     :: cell
    type(ELEMENT), intent(in)       :: elements(:)
    real(kind=DP), intent(in)       :: phi(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)       :: eps_full(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! variables
    real(kind=DP), allocatable  :: rho_pol(:,:,:)
    real(kind=DP), allocatable  :: phi_pol(:,:,:)
    real(kind=DP), allocatable  :: pot_dif(:,:,:)
    character(len=*), parameter :: myself = "multigrid_polarization_data"
    real(kind=DP)               :: E_pol, sgpce
    integer                     :: ierr
    logical                     :: backup_pub_is_bc_allow_frac_charge
    character(len=80)           :: orig_is_pbe

    call timer_clock(myself,1)
    call utils_trace_in(myself)

    ! ab: polarization density
    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
       write(stdout,'(a)',advance='no') 'Computing polarization density...'
    end if
    ! ab: allocate rho_pol
    allocate(rho_pol(mg%ld1f, mg%ld2f, mg%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'rho_pol',ierr)
    ! ab: compute polarization density
    call compute_rho_pol(pub_mg_defco_fd_order, eps_full, phi, grid, rho_pol)
    ! ab: calculation done
    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
        write(stdout,'(a)') ' done'
    end if
    ! ab: dump polarization density
    call visual_scalarfield(rho_pol, grid, cell, &
         'Polarization density (e/Ang**3)', '_rho_pol', &
         conversion_factor = ANGSTROM**3)

    ! ab: polarization potential
    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
       write(stdout,'(a)') 'Computing polarization potential using DL_MG...'
    end if
    ! ab: allocate phi_pol
    allocate(phi_pol(mg%ld1f, mg%ld2f, mg%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'phi_pol',ierr)

    ! jd: Temporarily disable bound_cond's check for integer charge
    backup_pub_is_bc_allow_frac_charge = pub_is_bc_allow_frac_charge
    pub_is_bc_allow_frac_charge = .true.

    ! ab: solve \nabla ^ 2 \phi_pol = -4 \pi \rho_pol
    ! which is basically a Poisson equation in vacuum.
    ! Storing original value of pub_is_pbe
    orig_is_pbe = pub_is_pbe
    pub_is_pbe = 'NONE'          !     in the vacuum calculation.
    E_pol = multigrid_calculate_hartree(phi_pol, rho_pol, grid, cell, &
         no_dump = .true.)
    ! ab: reverting pub_is_pbe back to original
    pub_is_pbe = orig_is_pbe
    pub_is_bc_allow_frac_charge = backup_pub_is_bc_allow_frac_charge

    ! ab: dump polarization potential
    call visual_scalarfield(phi_pol, grid, cell, &
         'Polarization potential (V)', '_phi_pol', &
         conversion_factor = -HARTREE_IN_EVS)

    ! ab: allocate the array that will hold the difference in electrostatic
    ! potential due to point ions and smeared ions.
    ! This is used to correct smeared ion polarization energy in solvent
    allocate(pot_dif(grid%ld1, grid%ld2, grid%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'pot_dif',ierr)
    call smeared_ion_point_ion_potential_difference(pot_dif, elements, &
         grid, cell, pub_multigrid_bc_is_periodic)

    sgpce = 0.5_DP * integrals_product_on_grid(grid, pot_dif, rho_pol)

    if(pub_on_root .and. pub_output_detail >= VERBOSE) then
       write(stdout,'(a,f14.6)') 'Polarization energy &
            &(1/2*\int\phi_pol\rho_pol):                  ', E_pol
       write(stdout,'(a,f14.6)') 'Smeared ion polariza&
            &tion corr. energy (1/2\int\pot_dif\rho_pol): ', sgpce
       write(stdout,'(a)')
    end if

    ! ab: Deallocate variables
    deallocate(pot_dif,stat=ierr)
    call utils_dealloc_check(myself,'pot_dif',ierr)
    deallocate(phi_pol,stat=ierr)
    call utils_dealloc_check(myself,'phi_pol',ierr)
    deallocate(rho_pol,stat=ierr)
    call utils_dealloc_check(myself,'rho_pol',ierr)

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine multigrid_polarization_data
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine compute_rho_pol(order, eps, v_i, grid, rho_pol)
    !=========================================================================!
    ! This subroutine computes the polarisation density for a Poisson equation!
    ! in a dielectric medium.                                                 !
    !                                                                         !
    !    rho_pol =  \nabla \cdot [(eps-1) \grad v_i]                          !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !    order       (input): finite difference order to use                  !
    !    eps         (input): the diel. functional on the full grid points    !
    !    v_i         (input): the current approximation of the potential      !
    !    grid        (input): the defect is computed on the mg grid, but the  !
    !                         host grid is needed here and there.             !
    !    rho_pol     (output): defined as above                               !
    !                                                                         !
    ! Arrays are assumed to be in the distributed slab representation.        !
    !-------------------------------------------------------------------------!
    ! Written by Arihant Bhandari on 30/09/2020 from compute_current_defect() !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: FOUR_PI
    use finite_differences, only: finite_difference_gradient, &
         finite_difference_laplacian
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_sanity_check, utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in) :: grid
    integer,       intent(in)   :: order
    real(kind=DP), intent(in)   :: eps(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)   :: v_i(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(out)  :: rho_pol(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! jd: Internal variables
    real(kind=DP), allocatable     :: lap_v(:,:,:)         ! nabla^2 v
    real(kind=DP), allocatable     :: grad_eps(:,:,:,:)    ! nabla eps
    real(kind=DP), allocatable     :: grad_v(:,:,:,:)      ! nabla v
    integer                        :: ierr                 ! jd: Error flag
    character(len=*), parameter    :: myself = 'compute_rho_pol'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)
    call utils_sanity_check(eps,'eps in '//myself)
    call utils_sanity_check(v_i,'v_i in'//myself)

    ! Allocate memory for local variables
    allocate(lap_v(grid%ld1,grid%ld2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'lap_v',ierr)
    allocate(grad_eps(grid%ld1,grid%ld2,grid%max_slabs12,3),stat=ierr)
    call utils_alloc_check(myself,'grad_eps',ierr)
    allocate(grad_v(grid%ld1,grid%ld2,grid%max_slabs12,3),stat=ierr)
    call utils_alloc_check(myself,'grad_v',ierr)

    ! Calculate the gradient and Laplacian of the potential
    call finite_difference_gradient(grad_v,v_i,order,mg,grid)
    call finite_difference_laplacian(lap_v,v_i,order,mg,grid)
    ! jd: Calculate the gradient of epsilon
    call finite_difference_gradient(grad_eps,eps,order,mg,grid)

    call utils_sanity_check(lap_v,'lap_v in '//myself)
    call utils_sanity_check(grad_v(:,:,:,1),'grad_v/x in '//myself)
    call utils_sanity_check(grad_v(:,:,:,2),'grad_v/y in '//myself)
    call utils_sanity_check(grad_v(:,:,:,3),'grad_v/z in '//myself)
    call utils_sanity_check(grad_eps(:,:,:,1),'grad_eps/x in '//myself)
    call utils_sanity_check(grad_eps(:,:,:,2),'grad_eps/y in '//myself)
    call utils_sanity_check(grad_eps(:,:,:,3),'grad_eps/z in '//myself)

    ! Now have everything we need to compute the rho_pol
    rho_pol(:,:,:) = 1.0_DP / FOUR_PI * ( &
         (eps(:,:,:)-1.0_DP) * lap_v(:,:,:) + &
         grad_eps(:,:,:,1) * grad_v(:,:,:,1) + &
         grad_eps(:,:,:,2) * grad_v(:,:,:,2) + &
         grad_eps(:,:,:,3) * grad_v(:,:,:,3))

    ! jd: Clean up
    deallocate(grad_v,stat=ierr)
    call utils_dealloc_check(myself,'grad_v',ierr)
    deallocate(grad_eps,stat=ierr)
    call utils_dealloc_check(myself,'grad_eps',ierr)
    deallocate(lap_v,stat=ierr)
    call utils_dealloc_check(myself,'lap_v',ierr)

    call utils_trace_out(myself)

  end subroutine compute_rho_pol
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_point_ion_potential_difference(phi_dif, &  ! out
       elements, grid, cell, is_periodic)                           ! in
  !===========================================================================!
  ! Calculate the difference in electrostatic potential due to point ions and !
  ! smeared ions:                                                             !
  ! phi_dif = \sum_I -Z_I/|r-R_I|*erfc(|r-R_I|/sigma).                        !
  !---------------------------------------------------------------------------!
  ! jd: NB. Conceptually this belongs in is_smeared_ions_mod, but this would  !
  !     lead to a circular dependency -- is_smeared_ions_mod depends on       !
  !     multigrid_methods_mod for multigrid_hartree().                        !
  !---------------------------------------------------------------------------!
  ! written by Arihant Bhandari in Sep 2020                                   !
  !===========================================================================!
    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id
    use constants, only: DP, safe_div_eps
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+), OPERATOR(*)
    use ion, only: ELEMENT
    use rundat, only: pub_is_smeared_ion_width
!$  use rundat, only: pub_threads_max
    use simulation_cell, only: CELL_INFO, minimum_image_distance_ortho
    use utils, only: utils_abort, utils_erfc, utils_trace_in, utils_trace_out

    implicit none

    type(ELEMENT), intent(in) :: elements(:)
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    logical, intent(in) :: is_periodic(3)
    real(kind=DP), intent(out) :: phi_dif(grid%ld1, grid%ld2, grid%max_slabs12)
    ! variables and parameters
    real(kind=DP) :: Z_I, d, sigma_I, term, accum
    integer :: ipt, i1, i2, i3, islab12, I
    type(POINT) :: R_I, r
    !------------------------------------------------------------------------
    call utils_trace_in('smeared_ion_point_ion_potential_difference')

    phi_dif=0.0_DP
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(ipt,i1,i2,islab12,i3,r,accum,I,R_I,Z_I,sigma_I,d,term) &
!$OMP SHARED(elements,cell,is_periodic,grid,pub_is_smeared_ion_width, &
!$OMP      phi_dif,pub_my_proc_id,pub_threads_max)
    do ipt=1,grid%num_my_slabs12*grid%n2*grid%n1
       i1 = modulo(ipt-1,grid%n1) + 1
       i2 = modulo((ipt-i1)/grid%n1,grid%n2) + 1
       islab12 = (ipt-(i2-1)*grid%n1-i1) / (grid%n1*grid%n2) + 1
       i3 = grid%first_slab12(pub_my_proc_id) + islab12 - 1

       r = &
            real((i1-1),kind=DP) * grid%da1 + &
            real((i2-1),kind=DP) * grid%da2 + &
            real((i3-1),kind=DP) * grid%da3

       accum = 0.0_DP
       do I = 1, size(elements)
          R_I = elements(I)%centre
          ! ab: electronic convention, ion charges treated -ve
          Z_I = -elements(I)%ion_charge
          sigma_I = pub_is_smeared_ion_width
          if (.not. any(is_periodic)) then ! all OBCs
             d = magnitude(r-R_I)
          else if(all(is_periodic)) then ! all PBCs
             d = minimum_image_distance_ortho(r,R_I,cell)
          else ! Mixed periodic/open BCs (Not currently supported)
             call utils_abort("Error in smeared_ion_generate_density: &
                  &Mixed periodic/open BCs are not currently supported &
                  &for smeared ions.")
          end if
          ! ab: avoid the singularity at d=0, put it as zero there
          ! as while calculating sgpce \rho_pol is zero there
          if (d < safe_div_eps) then
             term = 0.0_DP
          else
             term = utils_erfc(d/sigma_I) / d
          end if
          accum = accum + Z_I * term
       end do ! over I
       phi_dif(i1,i2,islab12) = accum
    end do ! ipt
!$OMP END PARALLEL DO

    call utils_trace_out('smeared_ion_point_ion_potential_difference')

  end subroutine smeared_ion_point_ion_potential_difference
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module multigrid_methods
