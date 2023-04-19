! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                  Density kernel utility module                 !
!                                                                !
! This module implements several utilities for the density       !
! kernel
!----------------------------------------------------------------!
! Written by Peter Haynes, 18/11/04                              !
! Modified by Nicholas Hine, December 2007 to reduce number of   !
! sparse algebra operations                                      !
! Modified by J.M. Escartin to remove its hidden global state by !
! using a kernel container (Fall 2014).                          !
!================================================================!

module kernel

  use constants, only: DP
  use rundat, only: pub_debug_on_root
  use sparse_array, only: SPAM3_ARRAY
  use sparse_embed, only: SPAM3_EMBED_ARRAY

  implicit none

  private

  public :: kernel_num_bound_states
  public :: kernel_bandgap
  public :: kernel_rms_err
  public :: kernel_occupancy_bounds
  public :: kernel_middle_occupancy
  public :: kernel_fix
  public :: kernel_normalise

  public :: kernel_rescale
  public :: kernel_rescale_spam3
  public :: kernel_rescale_fraglocalised_spam3

  public :: kernel_purify
  public :: kernel_commutator
  public :: kernel_init_core_ham
  public :: kernel_occ_check

  public :: kernel_create
  public :: kernel_destroy
  public :: kernel_workspace_create
  public :: kernel_workspace_destroy

  public :: kernel_workspace_invalidate
  public :: kernel_validate_ks
  public :: kernel_validate_ksk

  public :: kernel_copy
  public :: kernel_from_vecs_asc
  public :: kernel_from_vecs_asc2

  !public :: kernel_basis_transform
  public :: kernel_basis_update
  !public :: kernel_christoffel

  public :: kernel_build_subsystem

  ! Interfaces
  interface kernel_rescale
     module procedure kernel_rescale_int
     module procedure kernel_rescale_real
  end interface

  ! Kernel container
  type, public :: DKERN
     type(SPAM3_EMBED_ARRAY) :: kern
     type(SPAM3_EMBED_ARRAY) :: ks, ksk
     logical           :: workspace_created = .false.
     logical           :: ks_valid = .false.
     logical           :: ksk_valid = .false.
     ! rc2013: structure info for the SPAM3_ARRAYs (no embedding subscripts)
     character(len=10) :: core_struc(3)
  end type


contains

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_create(denskern, structure, is_cmplx)

     !=======================================================================!
     ! This subroutine creates the sparse matrices within the kernel of a    !
     ! DKERN object.                                                         !
     !=======================================================================!
     ! Arguments:                                                            !
     ! denskern (input/output) : Density kernel in DKERN format.             !
     ! structure (input)       : Structure of the sparse matrices.           !
     ! is_cmplx (input)        : (Optional) flag to create a complex kernel. !
     !=======================================================================!
     ! Written by J.M. Escartin, Fall 2014.                                  !
     ! Modified by Andrea Greco on 04/05/2015 to allow complex matrices      !
     ! Modified for embedding by Robert Charlton, July 2017.                 !
     !=======================================================================!

     use rundat, only: pub_num_spins, pub_num_kpoints
     use sparse_embed, only: sparse_embed_array_create
     implicit none

     ! Arguments
     type(DKERN), intent(inout)    :: denskern
     character(len=*), intent(in)  :: structure
     logical, intent(in), optional :: is_cmplx

     ! Local variables
     logical :: loc_cmplx

     ! agrecocmplx
     if (present(is_cmplx)) then
         loc_cmplx = is_cmplx
     else
         loc_cmplx = .false.
     end if

     ! rc2013: set core_struc i.e. strucure code without embedding subscripts
     denskern%core_struc(1) = trim(structure)

     call sparse_embed_array_create(denskern%kern, n_spins=pub_num_spins, &
          n_kpoints=pub_num_kpoints, structure=structure, iscmplx=loc_cmplx)

  end subroutine kernel_create

!------------------------------------------------------------------------------

  subroutine kernel_workspace_create(denskern, overlap)

     !=======================================================================!
     ! This subroutine creates the ks and ksk workspaces within a DKERN      !
     ! object.                                                               !
     !=======================================================================!
     ! Arguments:                                                            !
     ! denskern (input/output) : Density kernel in DKERN format.             !
     ! overlap (input)         : Overlap matrix in SPAM3 format.             !
     !=======================================================================!
     ! Written by Nick Hine, November 2007.                                  !
     ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
     ! Adapted for embedding by Robert Charlton, 23/07/2017.                 !
     !=======================================================================!

     use constants, only: stdout
     use rundat, only: pub_debug_on_root
     use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_create
     use utils, only: utils_assert, utils_alloc_check

     implicit none

     ! Arguments
     type(DKERN), intent(inout) :: denskern
     type(SPAM3_EMBED), intent(in)    :: overlap

     if (pub_debug_on_root) then
        write(stdout,'(a)') 'DEBUG: Entering kernel_workspace_create'
     end if

     call utils_assert(.not. denskern%workspace_created, &
          'Density kernel workspaces already created at entry to subroutine &
          &kernel_workspace_create')

     ! Initialise the sparse arrays
     call sparse_embed_array_create(denskern%ks, denskern%kern, overlap)
     call sparse_embed_array_create(denskern%ksk, denskern%ks, denskern%kern)

     denskern%workspace_created = .true.
     call kernel_workspace_invalidate(denskern)

     if (pub_debug_on_root) then
        write(stdout,'(a)') 'DEBUG: Leaving kernel_workspace_create'
     end if

  end subroutine kernel_workspace_create

!------------------------------------------------------------------------------

  subroutine kernel_destroy(denskern)
     !=======================================================================!
     ! This subroutine destroys a density kernel of type DKERN.              !
     ! Written by J.M. Escartin, Fall 2014.                                  !
     !=======================================================================!

     use sparse_embed, only: sparse_embed_array_destroy
     implicit none

     ! Arguments
     type(DKERN), intent(inout) :: denskern

     if (denskern%workspace_created) then
        call kernel_workspace_destroy(denskern)
     end if

     call sparse_embed_array_destroy(denskern%kern)

  end subroutine kernel_destroy

!------------------------------------------------------------------------------

  subroutine kernel_workspace_destroy(denskern)

     !=======================================================================!
     ! This subroutine destroys the ks and ksk workspaces within a DKERN     !
     ! object.                                                               !
     !=======================================================================!
     ! Argument:                                                             !
     ! denskern (input/output) : Density kernel in DKERN format.             !
     !=======================================================================!
     ! Written by Nick Hine, November 2007.                                  !
     ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
     ! Adapted for embedding by Robert Charlton, 23/07/2017.                 !
     !=======================================================================!

     use sparse_embed, only: sparse_embed_array_destroy
     implicit none

     ! Arguments
     type(DKERN), intent(inout) :: denskern

     call sparse_embed_array_destroy(denskern%ksk)
     call sparse_embed_array_destroy(denskern%ks)

     call kernel_workspace_invalidate(denskern)
     denskern%workspace_created = .false.

  end subroutine kernel_workspace_destroy

!------------------------------------------------------------------------------

  subroutine kernel_copy(dkern_dest, dkern_src)
     !=======================================================================!
     ! This subroutine copies a density kernel of type DKERN to another      !
     ! density kernel of type DKERN.                                         !
     !=======================================================================!
     ! Arguments:                                                            !
     ! dkern_src (input)         : Source density kernel.                    !
     ! dkern_dest (input/output) : Destionation density kernel.              !
     !=======================================================================!
     ! Written by J.M. Escartin, Fall 2014.                                  !
     ! Modified for embedding by Robert Charlton, 14/08/2017.                !
     !=======================================================================!

     use sparse_embed, only: sparse_embed_array_copy, sparse_embed_array_create, &
          sparse_embed_array_check
     use utils, only: utils_assert
     implicit none

     ! Arguments
     type(DKERN), intent(in)    :: dkern_src
     type(DKERN), intent(inout) :: dkern_dest

     ! Check consistency of source kernel.
     call utils_assert ( sparse_embed_array_check(dkern_src%kern) == 1 , &
          'Error in kernel_copy: source kernel does not contain valid SPAM3_EMBED_ARRAY data.')

     ! Check destination kernel.
     if ( sparse_embed_array_check(dkern_dest%kern) /= 1 ) then
        call sparse_embed_array_create(dkern_dest%kern, dkern_src%kern)
     end if

     ! Copy kernel content.
     call sparse_embed_array_copy(dkern_dest%kern, dkern_src%kern)

     ! Copy workspaces (if already created in the source density kernel).
     if ( dkern_src%workspace_created ) then

        ! If destination workspaces are not created, create them.
        if ( .not. dkern_dest%workspace_created ) then

           call sparse_embed_array_create(dkern_dest%ks, dkern_src%ks)
           call sparse_embed_array_create(dkern_dest%ksk, dkern_src%ksk)

           dkern_dest%workspace_created = .true.
        end if

        ! Copy contents and validity flags of workspaces.
        call sparse_embed_array_copy(dkern_dest%ks, dkern_src%ks)
        call sparse_embed_array_copy(dkern_dest%ksk, dkern_src%ksk)

        dkern_dest%ks_valid = dkern_src%ks_valid
        dkern_dest%ksk_valid = dkern_src%ksk_valid

     end if

  end subroutine kernel_copy

!------------------------------------------------------------------------------

  subroutine kernel_workspace_invalidate(denskern)

     !=======================================================================!
     ! This subroutine sets the flags of the workspaces within a DKERN       !
     ! object to invalid.                                                    !
     !=======================================================================!
     ! Argument:                                                             !
     ! denskern (input/output) : Density kernel in DKERN format.             !
     !=======================================================================!
     ! Written by Nick Hine, November 2007.                                  !
     ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
     ! Adapted for embedding by Robert Charlton, 23/07/2017.                 !
     !=======================================================================!

     implicit none
     type(DKERN), intent(inout) :: denskern

     denskern%ks_valid = .false.
     denskern%ksk_valid = .false.

  end subroutine kernel_workspace_invalidate

!------------------------------------------------------------------------------

  subroutine kernel_validate_ks(denskern, overlap)

     !=======================================================================!
     ! Given a DKERN object, this subroutine recalculates its ks workspace   !
     ! if required and sets its corresponding flag to valid.                 !
     !=======================================================================!
     ! Arguments:                                                            !
     ! denskern (input/output) : Density kernel in DKERN format              !
     ! overlap (input)         : Overlap matrix in SPAM3 format              !
     !=======================================================================!
     ! Written by Nick Hine, November 2007                                   !
     ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
     !=======================================================================!

     use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_product

     implicit none

     ! Arguments
     type(DKERN), intent(inout) :: denskern
     type(SPAM3_EMBED), intent(in) :: overlap

     if (.not. denskern%ks_valid) then
        call sparse_embed_array_product(denskern%ks, denskern%kern, overlap)
        denskern%ks_valid = .true.
     endif

  end subroutine kernel_validate_ks

!------------------------------------------------------------------------------

  subroutine kernel_validate_ksk(denskern, overlap)

     !=======================================================================!
     ! Given a DKERN object, this subroutine recalculates its ksk workspace  !
     ! if required and sets its corresponding flag to valid.                 !
     !=======================================================================!
     ! Arguments:                                                            !
     ! denskern (input/output) : Density kernel in DKERN format              !
     ! overlap (input)         : Overlap matrix in SPAM3 format              !
     !=======================================================================!
     ! Written by Nick Hine, November 2007                                   !
     ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
     !=======================================================================!

     use constants, only: stdout
     use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_product

     implicit none

     ! Arguments
     type(DKERN), intent(inout) :: denskern
     type(SPAM3_EMBED), intent(in)    :: overlap

     if (pub_debug_on_root) then
        write(stdout,'(a)') 'DEBUG: Entering kernel_validate_ksk.'
     end if

     if (.not. denskern%ksk_valid) then

        if (.not. denskern%ks_valid ) then
           call kernel_validate_ks(denskern, overlap)
        endif

        call sparse_embed_array_product(denskern%ksk, denskern%ks, denskern%kern)
        denskern%ksk_valid = .true.

     endif

     if (pub_debug_on_root) then
        write(stdout,'(a)') 'DEBUG: Leaving kernel_validate_ksk.'
     end if

  end subroutine kernel_validate_ksk

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  real(kind=DP) function kernel_rms_err(denskern, overlap)

    !=======================================================================!
    ! This subroutine returns the rms occupancy error of a given density    !
    ! kernel by calculating the value of the penalty functional.            !
    !                                                                       !
    ! See P.D. Haynes and M.C. Payne, Phys. Rev. B 59, 12173 (1999).        !
    !=======================================================================!
    ! Arguments:                                                            !
    ! denskern (input/output) : Density kernel in DKERN format              !
    ! overlap (input)         : Overlap matrix in SPAM3 format              !
    !=======================================================================!
    ! Written by Nicholas Hine, September 2009.                             !
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
    ! Modified for embedding by Robert Charlton, July 2017.                 !
    !=======================================================================!

    use constants, only: DP
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_copy, &
         sparse_embed_array_create, sparse_embed_array_destroy, &
         sparse_embed_array_num_rows, sparse_embed_array_product, &
         sparse_embed_array_trace, sparse_embed_array_scale
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(DKERN), intent(inout) :: denskern     ! Density kernel
    type(SPAM3_EMBED), intent(in) :: overlap               ! Overlap matrix

    ! Local Variables
    type(SPAM3_EMBED_ARRAY) :: ksc, kscksc
    real(kind=DP) :: pen
    real(kind=DP), allocatable :: trace_array(:,:)
    logical :: destroy_workspace
    integer :: ierr

    ! Create workspace if required
    destroy_workspace = .false.
    if (.not.denskern%workspace_created) then
       destroy_workspace = .true.
       call kernel_workspace_create(denskern, overlap)
    end if

    ! Initialise
    !! KPOINTS_DANGER
    pen = 0.0_DP

    ! Create temporary matrices
    call sparse_embed_array_create(ksc, denskern%ks)
    call sparse_embed_array_create(kscksc, ksc, ksc)

    ! ndmh: recalculate K.S if it has become invalid
    call kernel_validate_ks(denskern, overlap)

    ! Calculate ksc := 1 - K.S
    call sparse_embed_array_copy(ksc, denskern%ks)
    call sparse_embed_array_scale(ksc, -1.0_DP, 1.0_DP)

    ! Calculate P = Tr[KSKS(1-KS)(1-KS)]
    call sparse_embed_array_product(kscksc, denskern%ks, ksc)

    allocate(trace_array(kscksc%num_spins,kscksc%num_kpoints),stat=ierr)
    call utils_alloc_check('kernel_rms_err','trace_array',ierr)
    !! KPOINTS_DANGER
    call sparse_embed_array_trace(trace_array,kscksc,kscksc)
    pen = pen + SUM(trace_array)
    deallocate(trace_array,stat=ierr)
    call utils_dealloc_check('kernel_rms_err','trace_array',ierr)

    pen = pen / denskern%kern%num_spins
    kernel_rms_err = sqrt(abs(pen)/sparse_embed_array_num_rows(denskern%kern))

    call sparse_embed_array_destroy(kscksc)
    call sparse_embed_array_destroy(ksc)

    ! Destroy workspace if we created it on entry
    if (destroy_workspace) call kernel_workspace_destroy(denskern)


  end function kernel_rms_err

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_occupancy_bounds(max_occ, min_occ, denskern, overlap)

    !==========================================================================!
    ! This subroutine returns the bounds on the occupancies of a given density !
    ! kernel by calculating the maximum and minimum eigenvalues of KS.         !
    !==========================================================================!
    ! Arguments:                                                               !
    ! denskern (inout)        : Density kernel in DKERN format                 !
    ! overlap (input)         : Overlap matrix in SPAM3 format                 !
    !==========================================================================!
    ! Written by Nicholas Hine, September 2009.                                !
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.           !
    !==========================================================================!

    use constants, only: DP
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_copy, &
         sparse_embed_array_create, sparse_embed_array_destroy, &
         sparse_embed_array_extremal_eigenvalues, sparse_embed_array_scale

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: max_occ(:,:)              ! Max KS eigenvalue
    real(kind=DP), intent(out) :: min_occ(:,:)              ! Min KS eigenvalue
    type(DKERN), intent(inout) :: denskern                  ! Density kernel
    type(SPAM3_EMBED), intent(in) :: overlap                ! Overlap matrix

    ! Local Variables
    logical :: destroy_workspace

    ! Create workspace if required
    destroy_workspace = .false.
    if (.not.denskern%workspace_created) then
       destroy_workspace = .true.
       call kernel_workspace_create(denskern, overlap)
    end if

    ! Recalculate K.S if required
    call kernel_validate_ks(denskern, overlap)

    ! Find maximum eigenvalues of KS and 1-KS
    ! KPOINTS_DANGER
    call sparse_embed_array_extremal_eigenvalues(denskern%ks, overlap, &
         max_occ, 0.001_DP)
    call sparse_embed_array_scale(denskern%ks, -1.0_DP, 1.0_DP)
    call sparse_embed_array_extremal_eigenvalues(denskern%ks, overlap, &
         min_occ, 0.001_DP)
    min_occ = 1.0_DP - min_occ

    ! Restore ks to avoid recalculating it next time
    call sparse_embed_array_scale(denskern%ks, -1.0_DP, 1.0_DP)

    ! Destroy workspace if we created it on entry
    if (destroy_workspace) call kernel_workspace_destroy(denskern)

  end subroutine kernel_occupancy_bounds

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_middle_occupancy(mid_occ, denskern, overlap)

    !==========================================================================!
    ! This subroutine returns the the occupancies of a given density kernel    !
    ! that is closest to 0.5.                                                  !
    ! Note that this is considerably more expensive than computing the         !
    ! extremal occupancy bounds.                                               !
    !==========================================================================!
    ! Arguments:                                                               !
    ! denskern (inout)        : Density kernel in DKERN format                 !
    ! overlap (input)         : Overlap matrix in SPAM3 format                 !
    !==========================================================================!
    ! Written by Robert Bell, 25/07/14                                         !
    ! reusing code from kernel_occupancy bounds and palser_mano_kernel_optimise!
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.           !
    ! Modified for embedding by Robert Charlton, 13/09/2017.                   !
    !==========================================================================!

    use constants, only: DP
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_copy, &
        sparse_embed_array_create, sparse_embed_array_destroy, &
        sparse_embed_array_extremal_eigenvalues, sparse_embed_array_product, &
        sparse_embed_array_scale

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: mid_occ(:,:)         ! mid KS eigenvalue
    type(DKERN), intent(inout) :: denskern             ! Density kernel
    type(SPAM3_EMBED), intent(in) :: overlap                 ! Overlap matrix

    ! Local Variables
    type(SPAM3_EMBED_ARRAY) :: ksks ! (ks-0.5)(ks-0.5)
    logical :: destroy_workspace

    ! Create workspace if required
    destroy_workspace = .false.
    if (.not.denskern%workspace_created) then
       destroy_workspace = .true.
       call kernel_workspace_create(denskern, overlap)
    end if

    ! create ksks matrix
    call sparse_embed_array_create(ksks, denskern%ks, denskern%ks)

    ! Recalculate K.S if required
    call kernel_validate_ks(denskern, overlap)

    ! ndmh: find closest eigenvalue to 0.5
    call sparse_embed_array_scale(denskern%ks, -1.0_DP, 0.5_DP)
    call sparse_embed_array_product(ksks, denskern%ks, denskern%ks)
    call sparse_embed_array_scale(ksks, -1.0_DP, 0.0_DP)
    call sparse_embed_array_extremal_eigenvalues(ksks, overlap, mid_occ)
    mid_occ = sqrt(abs(mid_occ)) + 0.5_DP

    ! Restore ks to avoid recalculating it next time
    call sparse_embed_array_scale(denskern%ks, -1.0_DP, 0.5_DP)

    ! Destroy workspace if we created it on entry
    if (destroy_workspace) call kernel_workspace_destroy(denskern)
    call sparse_embed_array_destroy(ksks)

  end subroutine kernel_middle_occupancy

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_fix(denskern, overlap, inv_overlap, miniter, const_ne)

    !=======================================================================!
    ! This subroutine iteratively improves the idempotency of the density   !
    ! kernel (for a fixed set of NGWFs) using the penalty functional.       !
    !                                                                       !
    ! See P.D. Haynes and M.C. Payne, Phys. Rev. B 59, 12173 (1999).        !
    !                                                                       !
    !=======================================================================!
    ! Arguments:                                                            !
    ! denskern (inout)        : Density kernel in DKERN format              !
    ! overlap (input)         : Overlap matrix in SPAM3 format              !
    ! inv_overlap (input)     : Inverse overlap matrix in SPAM3 format      !
    ! miniter (input)         : Minimum number of iterations to do          !
    ! const_ne (input)        : Flag to conserve number of electrons        !
    !=======================================================================!
    ! Written by Peter Haynes, November 2004                                !
    ! Improvements by Chris-Kriton Skylaris, June 2009                      !
    ! Speed improvements by Nicholas Hine, November 2009                    !
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
    !=======================================================================!

    use comms, only: pub_on_root, comms_bcast, pub_root_proc_id
    use constants, only: stdout, DP, VERBOSE
    use rundat, only: pub_output_detail, pub_kerfix, pub_maxit_kernel_fix
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_extremal_eigenvalues, &
         sparse_embed_array_copy,  sparse_embed_array_num_rows, &
         sparse_embed_array_trace
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments: input/output
    type(DKERN), intent(inout)     :: denskern     ! Density kernel

    ! Arguments: input only
    type(SPAM3_EMBED), intent(in)        :: overlap      ! Overlap matrix
    type(SPAM3_EMBED), intent(in)        :: inv_overlap  ! Inverse overlap matrix
    integer, intent(in), optional  :: miniter      ! Min num iterations
    logical, intent(in), optional  :: const_ne     ! Flag for constant Ne

    ! Local variables
    integer :: maxiter
    integer :: iter,miniter_local
    integer :: nspins
    integer :: is, ik, ierr
    integer, allocatable :: n_occ(:,:)
    real(kind=DP), parameter :: minbound = -0.36602540378443864676_DP
    real(kind=DP), parameter :: maxbound = 1.36602540378443864676_DP
    real(kind=DP) :: spin_fac
    real(kind=DP) :: pen
    real(kind=DP) :: rms_err,old_rms_err
    real(kind=DP) :: steplen
    real(kind=DP), allocatable :: ne(:,:), trace_array(:,:)
    real(kind=DP), allocatable :: min_occ(:,:), max_occ(:,:)
    logical :: local_const_ne
    logical :: converged,temp_workspace
    type(SPAM3_EMBED_ARRAY) :: ksc, ksks, ksgs, gsgs
    type(SPAM3_EMBED_ARRAY) :: dir

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering kernel_fix'

    ! Start timer
    call timer_clock('kernel_fix',1)

    if (pub_on_root .and. pub_output_detail >= VERBOSE) &
         write(stdout,'(/a)') '=================== Penalty &
         &functional idempotency correction =================='

    ! Local copy of miniter
    if (present(miniter)) then
       miniter_local = miniter
    else
       miniter_local = 1
    end if

    ! Local copy of const_ne
    if (present(const_ne)) then
       local_const_ne = const_ne
    else
       local_const_ne = .true.
    end if

    ! Establish number of spins
    nspins = denskern%kern%num_spins
    spin_fac = 2.0_DP / nspins

    ! Allocate local arrays
    allocate(n_occ(nspins, denskern%kern%num_kpoints), stat=ierr)
    call utils_alloc_check('kernel_fix', 'n_occ', ierr)
    allocate(ne(nspins, denskern%kern%num_kpoints), stat=ierr)
    call utils_alloc_check('kernel_fix', 'ne', ierr)
    allocate(min_occ(nspins, denskern%kern%num_kpoints), stat=ierr)
    call utils_alloc_check('kernel_fix', 'min_occ', ierr)
    allocate(max_occ(nspins, denskern%kern%num_kpoints), stat=ierr)
    call utils_alloc_check('kernel_fix', 'max_occ', ierr)

    allocate(trace_array(denskern%ks%num_spins, denskern%ks%num_kpoints), &
         stat=ierr)
    call utils_alloc_check('kernel_fix', 'trace_array', ierr)

    ! cks: intialise variables
    maxiter = pub_maxit_kernel_fix
    ne     =-1.0_DP
    n_occ  = 0
    rms_err     =-100.0_DP
    old_rms_err =-100.0_DP

    ! cks: determine electron numbers
    call kernel_validate_ks(denskern, overlap)
    if (local_const_ne .and. (pub_kerfix == 2) ) then
       call sparse_embed_array_trace(trace_array,denskern%ks)
       n_occ = nint(trace_array)
    endif

    ! Allocate local workspace
    call internal_workspace_alloc

    ! Set convergence flag
    converged = .false.

    ! Set bounds
    min_occ = 0.0_DP
    max_occ = 1.0_DP

    ! Loop over iterations
    do iter=1,maxiter

       old_rms_err = rms_err

       ! Calculate penalty functional value and (contravariant) gradient
       call internal_penalty_value_grad

       ! cks: Check for convergence only after having done at least
       ! cks: one iteration. This prevents early exit in cases of little
       ! cks: or no denskern truncation which would otherwise hinder
       ! cks: LNV convergence.
       if (iter > 1) then
          ! Check for convergence (i.e. idempotency)
          converged = (pen/sparse_embed_array_num_rows(denskern%kern) < &
               epsilon(1.0_DP) * 100.0_DP)
       end if

       call comms_bcast(pub_root_proc_id,converged)
       if (converged) exit

       ! ndmh: check that rms_err is not stuck (usually indicates kernel
       ! ndmh: with occupation numbers that do not sum to Ne)
       !if ((iter>1).and.(rms_err>old_rms_err*0.9)) then
       !   if (pub_on_root) write(stdout,'(a,a\a)') 'WARNING in kernel_fix:', &
       !        'RMS occupancy error not reducing:', &
       !        'Possibly indicates wrong occupation number'
       !end if
       ! old_rms_err = rms_err

       ! Make sure gradient preserves electron number if required
       if (local_const_ne .and. (pub_kerfix ==1)) call internal_correct_gradient

       ! Take optimal length step in this direction
       call internal_optimal_step(steplen)

       ! Print out details
       if (pub_on_root .and. pub_output_detail >= VERBOSE) then
          do ik = 1, denskern%kern%num_kpoints
             if (denskern%kern%num_kpoints > 1) then
                write(stdout, '(a,i0)') 'k-point: ', ik
             end if
             write(stdout,'(a,i3,a,f9.6,a,f7.4,a,f12.4,a)') &
                  ' RMS occupancy error at iteration ',iter,': ', &
                  rms_err,' (step ',steplen,' Ne1=',ne(1,ik),')'
             if (nspins == 2) write(stdout,'(a,f12.4)') &
                  '                                      &
                  &                         Ne2=',ne(2,ik)
          end do
       endif


       ! cks: rescale density kernel to correct ne
       if (local_const_ne .and. (pub_kerfix ==2)) &
            call kernel_rescale(denskern, overlap, n_occ, .false., .true.)


       ! Recalculate K.S
       call kernel_validate_ks(denskern, overlap)

       ! Check bounds for convergence
       ! KPOINTS_DANGER
       call kernel_occupancy_bounds(max_occ, min_occ, denskern, overlap)
       if (minval(min_occ) > minbound .and. &
            maxval(max_occ) < maxbound .and. iter >= miniter_local) exit

    end do

    ! ndmh: evaluate kernel occupancy bounds and kernel RMS occupancy error
    ! ndmh: if this was not already done in the loop
    ! ndmh: 04/06/10 - also evaluate Ne here as this was not done in loop
    if (maxiter==0) then
       call kernel_occupancy_bounds(max_occ, min_occ, denskern, overlap)
       rms_err = kernel_rms_err(denskern, overlap)
       old_rms_err = rms_err
       call sparse_embed_array_trace(ne,denskern%ks)
    end if

    ! cks: fix occupancies by diagonalisation if it is a hopeless case
    ! cks: and pub_kerfix == 2
    if ( ((rms_err > 0.001_DP) .and. (rms_err > old_rms_err*0.9_DP) ) &
         .and. (pub_kerfix == 2)  )then

       if (pub_on_root .and. pub_output_detail >= VERBOSE) &
            write(stdout,'(a)') &
            '               Reset of occupancies applied'

       call kernel_reset_occupancies(denskern, &  ! input-output
            overlap, n_occ )

       ! cks: denskern has changed: invalidate ks and ksk workspaces
       ! ebl: removing invalidate (already done in reset_occupancies above)
       ! call kernel_workspace_invalidate(denskern)
    endif



    ! Print results
    ! ndmh: only recalculate pen if displayed value > 0
    if ((.not. converged).and.(rms_err>0.0000005_DP)) &
         call internal_penalty_value_grad

    if (pub_on_root .and. pub_output_detail >= VERBOSE) then

       do ik = 1, denskern%kern%num_kpoints
          if (denskern%kern%num_kpoints > 1) then
             write(stdout, '(a,i0)') 'k-point: ', ik
          end if

          write(stdout,'(a,f9.6,a,f12.4,a)') &
               '            Final RMS occupancy error   : ', &
               rms_err,' (Ne1=',ne(1,ik),')'

          if (nspins == 2) write(stdout,'(a,f12.4)') &
               '                                     &
               &                Ne2=',ne(2,ik)

          if (nspins == 1) then
             write(stdout,'(a,2(f8.4,a))') &
                  '               Final occupancy bounds   : [',min_occ(1,ik), &
                  ',',max_occ(1,ik),']'
          else
             do is=1,nspins
                write(stdout,'(a,i1,a,2(f7.3,a))') &
                     '            Final occupancy bounds spin ',is,': [', &
                     min_occ(is,ik),',',max_occ(is,ik),']'
             end do
          end if
       end do
    end if

    ! Deallocate local arrays
    deallocate(n_occ, stat=ierr)
    call utils_dealloc_check('kernel_fix', 'n_occ', ierr)
    deallocate(ne, stat=ierr)
    call utils_dealloc_check('kernel_fix', 'ne', ierr)
    deallocate(min_occ, stat=ierr)
    call utils_dealloc_check('kernel_fix', 'min_occ', ierr)
    deallocate(max_occ, stat=ierr)
    call utils_dealloc_check('kernel_fix', 'max_occ', ierr)
    deallocate(trace_array, stat=ierr)
    call utils_dealloc_check('kernel_fix', 'trace_array', ierr)

    ! Deallocate local workspace
    call internal_workspace_dealloc

    if (pub_on_root .and. pub_output_detail >= VERBOSE) &
         write(stdout,'(a/)') repeat('=',80)

    ! Stop timer
    call timer_clock('kernel_fix',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving kernel_fix'

    return

  contains

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine internal_workspace_alloc

      use sparse_embed, only: sparse_embed_array_create

      implicit none

      ! Allocate a temporary workspace if we need one
      if(.not.denskern%workspace_created)then
         call kernel_workspace_create(denskern, overlap)
         temp_workspace=.true.
      else
         temp_workspace=.false.
      endif

      call sparse_embed_array_create(ksc, denskern%ks)
      call sparse_embed_array_create(ksks, denskern%ks, denskern%ks)
      call sparse_embed_array_create(ksgs, ksks)
      call sparse_embed_array_create(gsgs, ksks)

      call sparse_embed_array_create(dir, denskern%kern)

    end subroutine internal_workspace_alloc

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine internal_workspace_dealloc

      use sparse_embed, only: sparse_embed_array_destroy

      implicit none

      ! Deallocate workspace
      call sparse_embed_array_destroy(dir)
      call sparse_embed_array_destroy(gsgs)
      call sparse_embed_array_destroy(ksgs)
      call sparse_embed_array_destroy(ksks)
      call sparse_embed_array_destroy(ksc)

      ! ndmh: destroy temporary workspace if we created one
      if(temp_workspace)then
         call kernel_workspace_destroy(denskern)
      endif

    end subroutine internal_workspace_dealloc

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !--------------------------------------------------------------------!
    ! The following subroutine calculates the penalty functional value:  !
    !                                                                    !
    !   P = Tr[KSKS(1-KS)(1-KS)]                                         !
    !                                                                    !
    ! and gradient:                                                      !
    !                                                                    !
    !   G = 2(1-KS)(1-2KS)K                                              !
    !                                                                    !
    !--------------------------------------------------------------------!

    subroutine internal_penalty_value_grad

      use sparse_embed, only: sparse_embed_array_trace, sparse_embed_array_product, &
           sparse_embed_array_axpy, sparse_embed_array_scale

      implicit none

      ! Initialise
      ne = 0.0_DP ; pen = 0.0_DP

      ! ndmh: calculate K.S and K.S.K again if they have become invalid
      call kernel_validate_ks(denskern, overlap)
      call kernel_validate_ksk(denskern, overlap)

      ! Calculate ne = Tr[KS]
      call sparse_embed_array_trace(ne,denskern%ks)

      ! Calculate ksc := 1 - K.S
      call sparse_embed_array_copy(ksc, denskern%ks)
      call sparse_embed_array_scale(ksc, -1.0_DP, 1.0_DP)

      ! Calculate P = Tr[KSKS(1-KS)(1-KS)]
      call sparse_embed_array_product(ksks, denskern%ks, ksc)
      ! KPOINTS_DANGER
      call sparse_embed_array_trace(trace_array,ksks,ksks)
      pen = pen + SUM(trace_array)

      ! ndmh: a faster way to calculate dir, which also doesn't need to use
      ! ndmh: an extra kskc matrix (saving memory)

      ! ndmh: Calculate G = 2(2.KS.KSK - 3.KSK + K) in dir
      call sparse_embed_array_product(dir, denskern%ks, denskern%ksk)
      call sparse_embed_array_scale(dir,2.0_DP)
      call sparse_embed_array_axpy(dir, denskern%ksk, -3.0_DP)
      call sparse_embed_array_axpy(dir, denskern%kern, 1.0_DP)
      call sparse_embed_array_scale(dir, 2.0_DP)

      ne = ne * spin_fac
      pen = pen / nspins
      rms_err = sqrt(abs(pen)/sparse_embed_array_num_rows(denskern%kern))

    end subroutine internal_penalty_value_grad

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !--------------------------------------------------------------------!
    ! The following subroutine ensures that the gradient preserves the   !
    ! correct number of electrons i.e.                                   !
    !                                                                    !
    !   Tr[DS] = 0                                                       !
    !                                                                    !
    !--------------------------------------------------------------------!
    ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009   !
    !--------------------------------------------------------------------!

    subroutine internal_correct_gradient

      use sparse_embed, only: sparse_embed_trace, sparse_embed_array_trace, &
           sparse_embed_array_axpy
      use utils, only: utils_abort

      implicit none

      real(kind=DP),dimension(denskern%kern%num_spins, denskern%kern%num_kpoints) &
           :: trds, fac
      real(kind=DP) :: trsinvs

      ! Calculate Tr(Sinv.S)
      call sparse_embed_trace(trsinvs, inv_overlap, overlap)
      if (abs(trsinvs) < epsilon(1.0_DP)*sparse_embed_array_num_rows(denskern%kern)) then
         call utils_abort('Error in internal_correct_gradient (kernel_fix): &
              &Tr(Sinv.S) = 0')
      end if


      ! Calculate Tr(D.S)
      call sparse_embed_array_trace(trds, dir, overlap)

      fac = - trds / trsinvs

      ! Correct gradient
      call sparse_embed_array_axpy(dir, inv_overlap, fac)

    end subroutine internal_correct_gradient

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !--------------------------------------------------------------------!
    ! The following subroutine takes an optimal step in the direction G, !
    ! whose length is calculated from the analytic quartic form of the   !
    ! penalty functional                                                 !
    !--------------------------------------------------------------------!

    subroutine internal_optimal_step(steplen)

      use sparse_embed, only: sparse_embed_array_is_dense, sparse_embed_array_axpy

      implicit none

      real(kind=DP), intent(out) :: steplen

      real(kind=DP) :: a(4)   ! quartic coefficients

      ! a(1) = 2 Tr[KS(1-KS)(1-2KS)GS]
      ! a(2) = Tr[(1-2KS)(1-2KS)GSGS] - 2 Tr[KSGS(1-KS)GS]
      ! a(3) = -2 Tr[(1-2KS)GSGSGS]
      ! a(4) = Tr[GSGSGSGS]

      ! ndmh: for gradient below threshold, set step to -0.5 if kernel is
      ! ndmh: dense (reliable and much faster than explicit evaluation of
      ! ndmh: actual coefficients in dense cases... may be ok for
      ! ndmh: much sparser kernels and higher thresholds).
      if ((rms_err < 0.000005_DP) .and. sparse_embed_array_is_dense(denskern%kern)) then

         steplen = -0.50000_DP

         ! ndmh: else calculate step explicitly
      else

         ! Calculate coefficients of quartic
         call internal_quartic_coeffs(a)

         ! Calculate optimal step length
         call internal_quartic_min(a,steplen)

      end if

      ! Take optimal step
      call sparse_embed_array_axpy(denskern%kern, dir, steplen)
      ! ndmh: denskern has changed: invalidate ks and ksk workspaces
      call kernel_workspace_invalidate(denskern)

    end subroutine internal_optimal_step

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine internal_quartic_coeffs(a)

      use sparse_embed, only: sparse_embed_array_trace, sparse_embed_array_product

      implicit none

      ! Arguments
      real(kind=DP), intent(out) :: a(4)

      real(kind=DP) :: tr_ksgs, tr_ksksgs, tr_ksgsgs, tr_ksksgsgs, tr_ksksksgs, &
           tr_gsgs, tr_gsgsgs, tr_gsgsgsgs, tr_ksgsgsgs, tr_ksgsksgs

      ! KPOINTS_DANGER

      ! ndmh: New way of finding coefficients involving far fewer
      ! ndmh: sparse_product calls (at cost of 1 extra matrix stored)

      ! Calculate G.S in ksc
      call sparse_embed_array_product(ksc, dir, overlap)

      ! Calculate KS.KS
      call sparse_embed_array_product(ksks, denskern%ks, denskern%ks)

      ! Calculate KS.GS
      call sparse_embed_array_product(ksgs, denskern%ks, ksc)

      ! Calculate GS.GS
      call sparse_embed_array_product(gsgs, ksc, ksc)

      ! Calculate traces involving the above
      call sparse_embed_array_trace(trace_array, ksgs)
      tr_ksgs     = sum(trace_array)
      call sparse_embed_array_trace(trace_array, gsgs)
      tr_gsgs     = sum(trace_array)
      call sparse_embed_array_trace(trace_array, denskern%ks, ksgs)
      tr_ksksgs   = sum(trace_array)
      call sparse_embed_array_trace(trace_array, ksgs, ksc)
      tr_ksgsgs   = sum(trace_array)
      call sparse_embed_array_trace(trace_array, gsgs, ksc)
      tr_gsgsgs   = sum(trace_array)
      call sparse_embed_array_trace(trace_array, ksks, ksgs)
      tr_ksksksgs = sum(trace_array)
      call sparse_embed_array_trace(trace_array, ksks, gsgs)
      tr_ksksgsgs = sum(trace_array)
      call sparse_embed_array_trace(trace_array, ksgs, ksgs)
      tr_ksgsksgs = sum(trace_array)
      call sparse_embed_array_trace(trace_array, ksgs, gsgs)
      tr_ksgsgsgs = sum(trace_array)
      call sparse_embed_array_trace(trace_array, gsgs, gsgs)
      tr_gsgsgsgs = sum(trace_array)

      ! a(1) = 2 Tr[KS(1-KS)(1-2KS)GS]
      a(1) = 2.0_DP*(tr_ksgs - 3.0_DP*tr_ksksgs + 2.0_DP*tr_ksksksgs)

      ! a(2) = Tr[(1-2KS)(1-2KS)GSGS] - 2 Tr[KSGS(1-KS)GS]
      a(2) = tr_gsgs - 4.0_DP*tr_ksgsgs + 4.0_DP*tr_ksksgsgs - &
           2.0_DP*(tr_ksgsgs - tr_ksgsksgs)

      ! a(3) = -2 Tr[(1-2KS)GSGSGS]
      a(3) = - 2.0_DP*(tr_gsgsgs - 2.0_DP*tr_ksgsgsgs)

      ! a(4) = Tr[GSGSGSGS]
      a(4) = tr_gsgsgsgs

      a = a / nspins

    end subroutine internal_quartic_coeffs

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !--------------------------------------------------------------------!
    ! The following subroutine finds the global minimum of               !
    !   a1*x + a2*x**2 + a3*x**3 + a4*x**4 = 0                           !
    !--------------------------------------------------------------------!

    subroutine internal_quartic_min(a,x)

      use constants, only: DP, PI
      use utils, only: utils_assert

      implicit none

      ! Arguments
      real(kind=DP), intent(in) :: a(4)
      real(kind=DP), intent(out) :: x

      ! Local variables
      integer :: i
      logical :: foundmin
      real(kind=DP), parameter :: THIRD = 1.0_DP / 3.0_DP
      real(kind=DP) :: aa,bb,cc
      real(kind=DP) :: q,r,theta,rtq,y,z
      real(kind=DP) :: x1(3),q1(3),qmin

      ! If a4 vanishes then the gradient must be zero in which case choose
      ! step length to correspond to purification transformation
      if (abs(a(4)) < epsilon(1.0_DP)) then
         x = -0.5_DP
         return
      end if

      ! Differentiate and write result as  x**3 + a*x**2 + b*x + c = 0
      aa = 0.75_DP * a(3) / a(4)
      bb = 0.5_DP * a(2) / a(4)
      cc = 0.25_DP * a(1) / a(4)

      ! Solve cubic equation
      q = (aa*aa - 3.0_DP*bb)/9.0_DP
      r = (2*aa*aa*aa - 9.0_DP*aa*bb + 27.0_DP*cc)/54.0_DP
      if (r*r < q*q*q) then   ! three real roots
         rtq = sqrt(q)
         theta = acos(r/(rtq*q))
         x1(1) = -2.0_DP * rtq * cos(theta*THIRD) - aa*THIRD
         x1(2) = -2.0_DP * rtq * cos((theta+2.0_DP*PI)*THIRD) - aa*THIRD
         x1(3) = -2.0_DP * rtq * cos((theta-2.0_DP*PI)*THIRD) - aa*THIRD
         foundmin = .false.
         qmin = 0.0_DP ! qoh: Initialise to prevent compiler warning
         do i=1,3
            if (a(2)+3.0_DP*x1(i)*(a(3)+2.0_DP*a(4)*x1(i)) > 0.0_DP) then
               q1(i) = x1(i)*(a(1)+x1(i)*(a(2)+x1(i)*(a(3)+x1(i)*a(4))))
               if (foundmin) then
                  if (q1(i) < qmin) then
                     x = x1(i)
                     qmin = q1(i)
                  end if
               else
                  x = x1(i)
                  qmin = q1(i)
                  foundmin = .true.
               end if
            end if
         end do
         call utils_assert(foundmin, 'Error in internal_quartic_min &
              &(kernel_mod.F90): no minimum found')
      else   ! only one root -> must be a minimum
         y = -sign(1.0_DP,r)*(abs(r)+sqrt(r*r-q*q*q))**THIRD
         if (abs(y) < epsilon(1.0_DP)) then
            z = 0.0_DP
         else
            z = q / y
         end if
         x = y + z - aa*THIRD
      end if

    end subroutine internal_quartic_min

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  end subroutine kernel_fix

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_normalise(denskern, overlap, inv_overlap, n_occ)

    !=======================================================================!
    ! This subroutine normalises the density kernel so that it satisfies a  !
    ! normalisation constraint, by adding an amount of the contravariant    !
    ! electron number gradient (inverse overlap matrix)                     !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes, November 2004                                !
    ! Modified to support spin polarisation, Peter Haynes, July 2006        !
    ! Fixed bug that occured when the number of up and down electrons       !
    ! differed, Chris-Kriton Skylaris, 17 June 2009.                        !
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
    !=======================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, VERBOSE, stdout, &
         EDA_POLSIMUL, EDA_POLFRAGLOC_DEVEL, EDA_CTFRAGLOC_DEVEL
    use rundat, only: pub_output_detail, pub_eda_mode, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_trace, &
         sparse_embed_array_axpy, sparse_embed_array_trace

    implicit none

    ! Arguments
    type(DKERN), intent(inout) :: denskern     ! Density kernel
    type(SPAM3_EMBED), intent(in)    :: overlap      ! Overlap matrix
    type(SPAM3_EMBED), intent(in)    :: inv_overlap  ! Inverse overlap matrix
    integer, intent(in)        :: n_occ(:,:)

    ! Local variables
    integer :: nspins                  ! Number of spins
    integer :: is, ik                  ! Spin and k-point counters.
    real(kind=DP), dimension(denskern%kern%num_spins, denskern%kern%num_kpoints) :: &
         corrected_ne,   & ! Desired electron number
         uncorrected_ne, & ! Initial uncorrected electron number
         step,           & ! Amount by which to correct kernel
         final_ne          ! Final electron number
    real(kind=DP) :: rate_of_change_ne ! Rate of change of electron number as
    ! gradient is added
    real(kind=DP) :: spin_fac          ! Spin polarisation factor

    ! mjsp: If this is a polarisation EDA calculation,
    ! then instead use the fragment--localised rescaling routine to avoid
    ! unintended delocalisations between the fragments:
    if ((pub_eda_mode == EDA_POLSIMUL) .or. &
        (pub_eda_mode == EDA_POLFRAGLOC_DEVEL) .or. &
        (pub_eda_mode == EDA_CTFRAGLOC_DEVEL)) then
      call kernel_rescale_fraglocalised_spam3(denskern%kern%m(:,PUB_1K), overlap, .false.)
      return
    end if

    ! Establish spin polarisation status
    nspins = denskern%kern%num_spins
    spin_fac = 2.0_DP / nspins

    ! Calculate rate of change of electron number
    call sparse_embed_trace(rate_of_change_ne, inv_overlap, overlap)

    ! Calculate current electron number
    call sparse_embed_array_trace(uncorrected_ne, denskern%kern, overlap)

    ! KPOINTS_DANGER
    if (pub_output_detail >= VERBOSE .and. pub_on_root) then
       do ik = 1, denskern%kern%num_kpoints
          if (denskern%kern%num_kpoints > 1) then
             write(stdout, '(3x,a,i0)') 'k-point ', ik
          end if
          do is = 1, nspins
             write(stdout,'(a,i1,a,f16.8)') '   Initial electron number for spin ', &
                  is,': ', uncorrected_ne(is, ik) * spin_fac
          end do
       end do
    end if

    ! Do the correction
    corrected_ne = real(n_occ, kind=DP)
    step = (corrected_ne - uncorrected_ne) / rate_of_change_ne

    call sparse_embed_array_axpy(denskern%kern, inv_overlap, step)

    ! ndmh: mark ks and ksk workspaces as invalid
    call kernel_workspace_invalidate(denskern)

    ! Check the correction if required
    if (pub_output_detail >= VERBOSE) then

       call sparse_embed_array_trace(final_ne, denskern%kern, overlap)

       final_ne = final_ne * spin_fac
       if (pub_on_root) then
          do ik = 1, denskern%kern%num_kpoints
             if (denskern%kern%num_kpoints > 1) then
                write(stdout, '(3x,a,i0)') 'k-point ', ik
             end if
             do is = 1, nspins
                write(stdout,'(a,i1,a,f16.8)') &
                     '     Final electron number for spin ', is, ': ', &
                     final_ne(is, ik)
             end do
          end do
       end if
    end if

  end subroutine kernel_normalise


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////


  subroutine kernel_rescale_spam3(denskern, overlap, corrected_ne, &
       silent, norm_fac)

    !=======================================================================!
    ! This subroutine rescales the density kernel so that it satisfies a    !
    ! normalisation constraint                                              !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes, November 2004                                !
    ! Spin polarised by Peter Haynes, July 2006                             !
    ! Fixed bug that occured when the number of up and down electrons       !
    ! differed, Chris-Kriton Skylaris, 5 May 2009.                          !
    ! Adapted for embedding structures by Robert Charlton, August 2017.     !
    ! Changed signature to take in a real t_occ instead of int n_occ - kkbd !
    !=======================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE, &
         EDA_POLSIMUL, EDA_POLFRAGLOC_DEVEL, &
         EDA_CTFRAGLOC_DEVEL
    use rundat, only: pub_output_detail, pub_eda_mode, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_scale, &
         sparse_embed_array_trace

    implicit none

    ! Arguments
    type(SPAM3_EMBED_ARRAY), intent(inout)   :: denskern     ! Density kernel
    type(SPAM3_EMBED), intent(in)      :: overlap      ! Overlap matrix
    real(kind=dp), intent(in)          :: corrected_ne(:,:)   ! pub_num_spins, kpts
    logical, intent(in),optional       :: silent       ! Whether to supress output
    real(kind=DP), optional            :: norm_fac(:,:)

    ! Local variables
    integer :: nspins               ! Number of spins
    integer :: is, ik               ! Spin and k-point counter
    real(kind=DP) :: spin_fac       ! Spin degeneracy factor
    real(kind=DP), dimension(denskern%num_spins, denskern%num_kpoints) :: &
         uncorrected_ne, &          ! Initial uncorrected electron number
         factor,         &          ! Scaling factor
         final_ne                   ! Final electron number
    logical       :: loc_silent     ! Whether to supress output

    ! ndmh: optional argument
    loc_silent = .false.
    if (present(silent)) then
       loc_silent = silent
    end if

    ! Determine number of spins
    nspins = denskern%num_spins
    spin_fac = 2.0_DP / nspins

    ! mjsp: If this is a polarisation EDA calculation,
    ! then instead use the fragment--localised rescaling routine to avoid
    ! unintended delocalisations between the fragments:
    if ((pub_eda_mode == EDA_POLSIMUL) .or. &
        (pub_eda_mode == EDA_POLFRAGLOC_DEVEL) .or. &
        (pub_eda_mode == EDA_CTFRAGLOC_DEVEL)) then
      call kernel_rescale_fraglocalised_spam3(denskern%m(:,PUB_1K), overlap, silent, norm_fac)
      return
    end if

    ! Calculate current electron number
    call sparse_embed_array_trace(uncorrected_ne, denskern, overlap)

    ! KPOINTS_DANGER
    if (pub_output_detail >= VERBOSE .and. pub_on_root .and. (.not.loc_silent)) then
       do ik = 1, denskern%num_kpoints
          if (denskern%num_kpoints > 1) then
             write(stdout, '(3x,a,i0)') 'k-point ', ik
          end if
          do is = 1, nspins
             write(stdout,'(a,i1,a,f16.8)') '    Initial Ne',&
                  is,': ', uncorrected_ne(is,ik)*spin_fac
          end do
       end do
    end if

    ! ndmh: protection against divide-by-zero in one spin channel
    where (uncorrected_ne == 0.0_DP)
       factor = 0.0_DP
    else where
       factor = corrected_ne / uncorrected_ne
    end where

    ! Rescale
    call sparse_embed_array_scale(denskern, factor)

    ! ndmh: store norm_fac if requested
    ! KPOINTS_DANGER
    if (present(norm_fac)) norm_fac = factor

    ! Check the correction if required
    if (pub_output_detail >= VERBOSE .and. (.not.loc_silent)) then

       call sparse_embed_array_trace(final_ne, denskern, overlap)
       if (pub_on_root) then
          do ik = 1, denskern%num_kpoints
             if (denskern%num_kpoints > 1) then
                write(stdout, '(3x,a,i0)') 'k-point ', ik
             end if
             do is = 1, nspins
                write(stdout,'(a,i1,a,f16.8)') &
                  '      Final Ne',is,': ', final_ne(is, ik) * spin_fac
             end do
          end do
       end if

    end if

  end subroutine kernel_rescale_spam3

!------------------------------------------------------------------------------

  subroutine kernel_rescale_fraglocalised_spam3(denskern, overlap, silent, norm_fac)

    !=======================================================================!
    ! This subroutine rescales the fragment blocks of the density kernel so !
    ! that each block satisfies a normalisation constraint.                 !
    ! The block-wise restriction prevents charge delocalisation between     !
    ! fragments from arising.                                               !
    !                                                                       !
    ! NOTE: the norm_fac value returned is the norm_fac value for the       !
    ! fragment block given by pub_frag_counter. This is because this        !
    ! routine is intended to be used to calculate fragment-specific         !
    ! polarisation contributions, and so norm_fac for all                   !
    ! other fragments (i.e. the fragments providing the polarizing field)   !
    ! will be equal to zero as these fragments remain frozen.               !
    ! TODO: The code does not function with                                 !
    ! pub_exact_lnv=.true. as elsewhere norm_fac is applied to the full     !
    ! denskern, which results in delocalisations between the fragments.     !
    !-----------------------------------------------------------------------!
    ! Written by Max Phipps April 2015.                                     !
    ! Modified from kernel_rescale_spam3                                    !
    ! Modified for embedding by Joseph Prentice, September 2018             !
    !=======================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, max_spins, VERBOSE
    use rundat, only : pub_output_detail, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_scale, &
         sparse_embed_trace, sparse_embed_copy, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product
    use rundat, only: pub_frag_iatm, pub_frag_counter, pub_frag_counter2
    use fragment_matrix_ops, only: fmo_copy_frag_block_spam_nxn
    use fragment_data, only: pub_frag_data

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: denskern(:)     ! Density kernel
    type(SPAM3_EMBED), intent(in)    :: overlap      ! Overlap matrix
    logical, intent(in),optional :: silent     ! Whether to supress output
    real(kind=DP), optional :: norm_fac(max_spins) ! TODO:currently hacked
                                                   ! do not rely on this!

    ! Local variables
    integer :: nspins                  ! Number of spins
    integer :: is                      ! Spin counter
    real(kind=DP) :: spin_fac          ! Spin degeneracy factor
    real(kind=DP) :: corrected_ne      ! Desired electron number
    real(kind=DP) :: uncorrected_ne    ! Initial uncorrected electron number
    real(kind=DP) :: factor            ! Scaling factor
    real(kind=DP) :: final_ne          ! Final electron number
    logical       :: loc_silent        ! Whether to supress output
    type(SPAM3_EMBED)   :: ks                ! K.S buffer
    type(SPAM3_EMBED)   :: ks_workspace      ! K.S workspace, used to calculate each
                                       ! fragment's Ne
    type(SPAM3_EMBED)   :: workspace_buf     ! workspace buffer
    type(SPAM3_EMBED)   :: k_workspace       ! K workspace, used for rescaling
    integer       :: it, it2           ! Fragment counter

    ! mjsp: optional argument
    loc_silent = .false.
    if (present(silent)) then
       loc_silent = silent
    end if

    ! mjsp: Determine number of spins
    nspins = size(denskern)
    spin_fac = 2.0_DP / nspins

    ! mjsp: initialise matrices
    call sparse_embed_create(ks, denskern(1))
    call sparse_embed_create(ks_workspace, denskern(1))
    call sparse_embed_create(workspace_buf, denskern(1))
    call sparse_embed_create(k_workspace, denskern(1))

    do is=1,nspins

       ! mjsp: ks = K.S
       call sparse_embed_product(ks, denskern(is), overlap)

       ! mjsp: loop the fragments
       do it=1,pub_frag_iatm(0)

          ! mjsp: if pub_frag_counter2 is used, then skip if at this fragment
          ! to avoid rescaling this fragment without its partner
          ! pub_frag_counter fragment:
          if (pub_frag_counter2 .eq. it) cycle


          ! mjsp: if pub_frag_counter2 is used, and pub_frag_counter=it
          if ((pub_frag_counter2 .ne. 0) .and. (pub_frag_counter .eq. it)) then

             ! mjsp: set optional second fragment ID in order to join fragments
             ! pub_frag_counter and pub_frag_counter2:
             it2 = pub_frag_counter2

             ! mjsp: include values of both fragments to combine these two fragments:
             corrected_ne = pub_frag_data(it)%rep%n_occ(is,PUB_1K) + &
                  pub_frag_data(it2)%rep%n_occ(is,PUB_1K)
             if (pub_output_detail >= VERBOSE .and. pub_on_root .and. (.not.loc_silent)) &
                  write(stdout,'(a,i2,a,i2)') ' -->Fragments: ', it, ' and ', it2

          else

             ! mjsp: if pub_frag_counter2 is unused, or this fragment is not to be
             ! combined with pub_frag_counter2...

             corrected_ne = pub_frag_data(it)%rep%n_occ(is,PUB_1K)
             if (pub_output_detail >= VERBOSE .and. pub_on_root .and. (.not.loc_silent)) &
                  write(stdout,'(a,i2)') ' -->Fragment: ', it

             ! mjsp: fragment ID2 is redundant:
             it2 = 0

          end if


          ! mjsp: copy ks to the ks_workspace
          call sparse_embed_copy(ks_workspace, ks)
          ! mjsp: copy k to the k_workspace
          call sparse_embed_copy(k_workspace, denskern(is))

          ! mjsp: ensure the workspace buffer is cleared
          call sparse_embed_scale(workspace_buf, 0.0_DP)

          ! mjsp: zero out all but the block that relates to this fragment
          call fmo_copy_frag_block_spam_nxn(workspace_buf, ks_workspace, it, &
               it2)
          call sparse_embed_copy(ks_workspace, workspace_buf)


          ! mjsp: Calculate current electron number for this fragment
          call sparse_embed_trace(uncorrected_ne, ks_workspace)

          if (pub_output_detail >= VERBOSE .and. pub_on_root .and. (.not.loc_silent)) &
               write(stdout,'(a,i1,a,f16.8)') '    Initial Ne',&
               is,': ', uncorrected_ne*spin_fac

          ! mjsp: protection against divide-by-zero in one spin channel
          if (uncorrected_ne == 0.0_DP) then
             factor = 0.0_DP
          else
             factor = corrected_ne / uncorrected_ne
          end if

          ! mjsp: Rescale the k_workspace to rescale the fragment block
          call sparse_embed_scale(k_workspace, factor)

          ! mjsp: Copy the fragment block back into denskern
          call fmo_copy_frag_block_spam_nxn(denskern(is), k_workspace, it, it2)

          ! mjsp: store norm_fac if requested
          ! mjsp: NOTE: This is not used and is untested.
          if ((it .eq. pub_frag_counter) &
            .and. present(norm_fac)) norm_fac(is) = factor

          ! mjsp: Check the correction if required
          if (pub_output_detail >= VERBOSE .and. (.not.loc_silent)) then

             ! mjsp: ensure the workspace buffer is cleared
             call sparse_embed_scale(workspace_buf, 0.0_DP)

             ! mjsp: isolate the fragment in the k_workspace
             call fmo_copy_frag_block_spam_nxn(workspace_buf, k_workspace, it, it2)
             call sparse_embed_copy(k_workspace, workspace_buf)

             ! mjsp: check N_e is correct
             call sparse_embed_trace(final_ne, k_workspace, overlap)
             if (pub_on_root) write(stdout,'(a,i1,a,f16.8)') &
                  '      Final Ne',is,': ',final_ne*spin_fac

          end if

       end do ! fragment loop

    end do  ! spin loop

    ! mjsp: Check the correction if required (final supermolecule)
    if (pub_output_detail >= VERBOSE .and. (.not.loc_silent)) then
       if (pub_on_root) write(stdout,'(a)') ' -->Supermolecule'

       do is=1,nspins
          call sparse_embed_trace(final_ne, denskern(is), overlap)
          if (pub_on_root) write(stdout,'(a,i1,a,f16.8)') &
               '      Final Ne',is,': ',final_ne*spin_fac
       enddo

    end if

    ! mjsp: cleanup
    call sparse_embed_destroy(ks)
    call sparse_embed_destroy(ks_workspace)
    call sparse_embed_destroy(workspace_buf)

  end subroutine kernel_rescale_fraglocalised_spam3

!------------------------------------------------------------------------------

  subroutine kernel_rescale_int(denskern, overlap, n_occ, can_rescale_ks, &
       silent, norm_fac)

    !=======================================================================!
    ! Wrapper to the new kernel_rescale which takes in real occupancies.    !
    ! Nothing is broken because n_occ was converted into a real when used   !
    ! anyway.                                                               !
    !-----------------------------------------------------------------------!
    ! Written by Kevin Duff, 2017.                                          !
    !=======================================================================!

    use constants, only: DP
    use sparse_embed, only: SPAM3_EMBED

    implicit none

    ! Arguments
    type(DKERN), intent(inout)   :: denskern       ! Density kernel
    type(SPAM3_EMBED), intent(in)      :: overlap        ! Overlap matrix
    integer, intent(in)          :: n_occ(:,:)
    logical, intent(in),optional :: can_rescale_ks
    logical, intent(in),optional :: silent         ! Whether to supress output
    real(kind=DP), optional      :: norm_fac(:,:)

    ! Internal
    real(kind=dp) :: t_occ(denskern%kern%num_spins, denskern%kern%num_kpoints)

    t_occ = real(n_occ,kind=DP)

    call kernel_rescale(denskern, overlap, t_occ, can_rescale_ks, silent, &
            norm_fac)

  end subroutine kernel_rescale_int

  subroutine kernel_rescale_real(denskern, overlap, t_occ, can_rescale_ks, &
       silent, norm_fac)

    !=======================================================================!
    ! This subroutine rescales the density kernel so that it satisfies a    !
    ! normalisation constraint                                              !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes, November 2004                                !
    ! Spin polarised by Peter Haynes, July 2006                             !
    ! Fixed bug that occured when the number of up and down electrons       !
    ! differed, Chris-Kriton Skylaris, 5 May 2009.                          !
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
    ! Accepts non-integer occupancies. Kevin Duff, 2017.                    !
    !=======================================================================!

    use constants, only: DP, EDA_POLSIMUL, EDA_POLFRAGLOC_DEVEL, &
         EDA_CTFRAGLOC_DEVEL, stdout
    use rundat, only: pub_eda_mode, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_scale

    implicit none

    ! Arguments
    type(DKERN), intent(inout)     :: denskern       ! Density kernel
    type(SPAM3_EMBED), intent(in)        :: overlap        ! Overlap matrix
    real(kind=dp), intent(in)            :: t_occ(:,:)
    logical, intent(in),optional         :: can_rescale_ks
    logical, intent(in),optional         :: silent         ! Whether to supress output
    real(kind=DP), optional              :: norm_fac(:,:)

    ! Local variables
    real(kind=DP) :: loc_norm_fac(denskern%kern%num_spins, denskern%kern%num_kpoints)

    ! mjsp: If this is a polarisation EDA calculation,
    ! mjsp: then instead use the fragment--localised rescaling routine
    if ((pub_eda_mode == EDA_POLSIMUL) .or. &
        (pub_eda_mode == EDA_POLFRAGLOC_DEVEL) .or. &
        (pub_eda_mode == EDA_CTFRAGLOC_DEVEL)) then

       ! fragment-wise rescaling:
       call kernel_rescale_fraglocalised_spam3(denskern%kern%m(:,PUB_1K), overlap, &
          silent, loc_norm_fac)

       ! ks cannot be directly rescaled as certain fragment blocks
       ! must remain frozen
       !if (present(can_rescale_ks) then
       !  if(denskern%ks_valid .and. can_rescale_ks) &
       !    call kernel_validate_ks(denskern, overlap)
       !end if

       call kernel_workspace_invalidate(denskern)

       if (present(norm_fac)) norm_fac(:,:) = loc_norm_fac(:,:)

    else

       if(pub_debug_on_root)write(stdout,*)&
                "DEBUG: KERNEL_RESCALE: rescaling denskern to ",t_occ(:,1)
       call kernel_rescale_spam3(denskern%kern, overlap, t_occ, &
          silent, loc_norm_fac)

       ! KPOINTS_DANGER

       ! ndmh: allow rescaling of ks rather than recalculation
       if (present(can_rescale_ks)) then
          if(denskern%ks_valid .and. can_rescale_ks) then

             call sparse_embed_array_scale(denskern%ks, loc_norm_fac)

             call kernel_workspace_invalidate(denskern)
             denskern%ks_valid=.true.
          else
             call kernel_workspace_invalidate(denskern)
          endif
       else
          call kernel_workspace_invalidate(denskern)
       endif

       if (present(norm_fac)) norm_fac(:,:) = loc_norm_fac(:,:)

    end if

  end subroutine kernel_rescale_real

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_purify(pur_denskern, denskern, overlap, inv_overlap, &
       n_occ, fixed_denskern, row, col)

    !=======================================================================!
    ! This subroutine accepts a density kernel and overlap matrix in        !
    ! SPAM3 format and returns a purified density kernel with the given     !
    ! sparsity pattern.  The subroutine is robust in the sense that there   !
    ! are no truncations in matrix multiplications.                         !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes on 19 November 2004                           !
    ! Based on earlier versions by Chris-Kriton Skylaris and Peter Haynes   !
    ! Modified by Chris-Kriton Skylaris on 28/03/2005 so that it conserves  !
    ! Ne during kernel_fix.                                                 !
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.        !
    ! Adapted for embedding by Robert Charlton, 24/07/17.                   !
    !-----------------------------------------------------------------------!

    use constants, only: stdout
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_copy, &
         sparse_embed_array_create, sparse_embed_array_destroy
    use timer, only: timer_clock

    implicit none

    ! Arguments:
    type(SPAM3_EMBED_ARRAY), intent(inout) :: pur_denskern ! Purified density kernel
    type(DKERN), intent(inout) :: denskern          ! Original density kernel
    type(SPAM3_EMBED), intent(in)    :: overlap           ! Overlap matrix
    type(SPAM3_EMBED), intent(in)    :: inv_overlap       ! Inverse overlap matrix
    integer, intent(in)              :: n_occ(:,:)        ! Number of occupied bands
    logical, intent(in), optional    :: fixed_denskern    ! Set this to true if
    ! denskern must be fixed
    ! rc2013: is this used anymore?
    integer, intent(in), optional    :: row    ! optional row specification
    integer, intent(in), optional    :: col    ! optional column specification

    ! Local variables
    integer :: nspins
    logical :: temp_workspace
    logical :: temp_denskern_copy
    logical :: loc_fixed_denskern
    type(SPAM3_EMBED_ARRAY) :: denskern_backup ! Backup of input denskern
    type(SPAM3_EMBED_ARRAY) :: tmp_purkern ! rc2013: temporary matrix if we
                                           ! need to extract a sub-block
    integer :: ierr
    ! in case it must be altered

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering kernel_purify'

    ! Start timer
    call timer_clock('kernel_purify',1)

    if (present(fixed_denskern)) then
       loc_fixed_denskern = fixed_denskern
    else
       loc_fixed_denskern = .false.
    end if

    ! Establish number of spins
    nspins = denskern%kern%num_spins

    ! Create temporary workspace if we need one
    if (.not.denskern%workspace_created) then
       temp_workspace=.true.
       call kernel_workspace_create(denskern, overlap)
    else
       temp_workspace=.false.
    endif
    temp_denskern_copy = .false.

    ! Check that purification will be stable
    call internal_check_denskern

    ! Now do the purification
    call internal_purify

    ! Restore input kernel if we altered it
    if (temp_denskern_copy.and.loc_fixed_denskern) then
       call sparse_embed_array_copy(denskern%kern, denskern_backup)
       call sparse_embed_array_destroy(denskern_backup)
    end if

    ! Destroy workspace if we created a temporary one
    if (temp_workspace) then
       call kernel_workspace_destroy(denskern)
    endif

    ! Stop timer
    call timer_clock('kernel_purify',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving kernel_purify'

  contains

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine internal_check_denskern

      use constants, only: DP
      use sparse_embed, only: sparse_embed_array_copy

      implicit none

      ! Local variables
      real(kind=DP), parameter :: maxocc = 1.36602540378443864676_DP
      real(kind=DP), parameter :: minocc = -0.36602540378443864676_DP
      real(kind=DP), dimension(denskern%kern%num_spins, denskern%kern%num_kpoints) :: &
           max_occ, min_occ

      if (pub_debug_on_root) then
         write(stdout,'(a)') 'DEBUG: Entering internal_check_denskernel (kernel_purify).'
      end if

      ! Calculate K.S product
      call kernel_validate_ks(denskern, overlap)

      ! Check bounds
      call kernel_occupancy_bounds(max_occ, min_occ, denskern, overlap)

      ! If the density kernel is stable
      if (.not. (maxval(max_occ) < maxocc .and. minval(min_occ) > minocc) ) then

         if (loc_fixed_denskern) then
            call sparse_embed_array_create(denskern_backup, denskern%kern)
            call sparse_embed_array_copy(denskern_backup, denskern%kern)
            temp_denskern_copy = .true.
         end if

         ! cks: It is massively important to normalise to Ne before the fix
         ! cks: to avoid obtaining wrong occupancies
         call kernel_rescale(denskern, overlap, n_occ, can_rescale_ks=.true.)

         ! Fix the kernel
         call kernel_fix(denskern, overlap, inv_overlap)

         ! Re-calculate product K.S
         call kernel_validate_ks(denskern, overlap)

      end if

      if (pub_debug_on_root) then
         write(stdout,'(a)') 'DEBUG: Leaving internal_check_denskernel (kernel_purify).'
      end if

    end subroutine internal_check_denskern

    subroutine internal_purify !TODO: maybe use pur_denskern ks's for this??

      use constants, only: DP
      use sparse_embed, only: sparse_embed_array_scale, &
          sparse_embed_array_product, sparse_embed_array_extract_sub

      implicit none

      if (pub_debug_on_root) then
         write(stdout,'(a)') 'DEBUG: Entering internal_purify (kernel_purify).'
      end if

      ! Calculate product K.S.K
      call kernel_validate_ksk(denskern, overlap)

      ! Put 3I - 2K.S in ks
      call sparse_embed_array_scale(denskern%ks,-2.0_DP,3.0_DP)

      ! Calculate product 3K.S.K - 2K.S.K.S.K in pur_denskern
      if(present(row) .or. present(col)) then
          ! rc2013: if we need a sub-matrix pass the density kernels into
          ! a temporary structure so we can extract the sub-blocks
          call sparse_embed_array_create(tmp_purkern, pur_denskern)
          call sparse_embed_array_product(tmp_purkern, denskern%ks, denskern%ksk)
          call sparse_embed_array_extract_sub(pur_denskern, tmp_purkern, row, col)
          call sparse_embed_array_destroy(tmp_purkern)
      else
          ! rc2013: otherwise just do a straightforward multiplication
          call sparse_embed_array_product(pur_denskern, denskern%ks, denskern%ksk)
      endif

      ! ndmh: ks was broken in calculating this, so invalidate it
      denskern%ks_valid = .false.

      if (pub_debug_on_root) then
         write(stdout,'(a)') 'DEBUG: Leaving internal_purify (kernel_purify).'
      end if

    end subroutine internal_purify

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  end subroutine kernel_purify

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  function kernel_num_bound_states(denskern,ham,overlap,&
       inv_overlap,num_ngwfs,num_occ,energy_range,energy_gap) &
       result(num_bound_states)
    !====================================================================!
    ! Function computes the number of bound states of the valence Ham-   !
    ! iltonian and returns it. For systems where the vacuum level is not !
    ! well defined, it is also possible to specify an energy range. In   !
    ! this case, the routine returns the number of unoccupied states in  !
    ! that energy range as measured from the highest occupied state.     !
    ! In case a energy_gap is specified, the code will, if possible,     !
    ! return the number of bound states plus the number of states to     !
    ! reach a gap of size energy_gap in the conduction space.            !
    ! Currently, the routine assumes spin_degeneracy.                    !
    !--------------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff, July 2015                              !
    ! Modified for embedding by Joseph Prentice, June 2018               !
    !====================================================================!

    ! agrecocmplx
    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_dealloc
    use rundat, only: pub_num_spins, pub_num_kpoints, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_copy, sparse_embed_axpy, &
         sparse_embed_extremal_eigenvalue, sparse_embed_product, &
         sparse_embed_outer_product
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    ! Arguments
    type(DKERN), intent(in) :: denskern       ! Density kernel
    type(SPAM3_EMBED), intent(in) :: ham(:)         ! hamiltonian
    type(SPAM3_EMBED), intent(in) :: overlap        ! overlap matrix
    type(SPAM3_EMBED), intent(in) :: inv_overlap    ! inverse overlap matrix
    integer, intent(in) :: num_ngwfs
    integer, intent(in) :: num_occ(:)
    real(kind=DP), optional, intent(in) :: energy_range ! energy range
                               ! of the conduction space manifold considered
    real(kind=DP), optional, intent(in) :: energy_gap ! if this is present,
      ! the code counts the number of conduction states beyond the energy
      ! range until a gap of size energy_gap in the conduction space many-
      ! fold is found.

    ! Result
    integer, dimension(denskern%kern%num_spins) :: num_bound_states

    ! local variables
    real(kind=DP) :: temp_energy_gap
    real(kind=DP) :: temp_val(denskern%kern%num_spins)
    real(kind=DP) :: max_eval(denskern%kern%num_spins)
    real(kind=DP) :: shift
    real(kind=DP) :: eval(denskern%kern%num_spins)
    real(kind=DP) :: prev_eval(denskern%kern%num_spins)
    ! agrecocmplx: use FUNCTIONS type to switch between
    ! real and complex eigenvectors
    type(FUNCTIONS), allocatable, dimension(:) :: evec
    type(SPAM3_EMBED) :: temp_proj_ham
    type(SPAM3_EMBED) :: cond_proj
    type(SPAM3_EMBED) :: SPS, SP, HPS
    type(SPAM3_EMBED), allocatable :: PS(:), S_invH(:)
    integer :: ierr
    integer :: num_extra_states
    integer :: is
    ! agrecocmplx
    logical :: loc_cmplx

    ! jme: KPOINTS_DANGER
    ! Some parts of this function enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Function kernel_num_bound_states not ready yet for more&
         & than one k-point.')

    num_extra_states=0
    temp_energy_gap=0.0_DP

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    ! allocate data
    ! agrecocmplx: allocate using appropriate routine
    allocate(evec(pub_num_spins),stat=ierr)
    call utils_alloc_check('kernel_num_bound_states','evec',ierr)
    do is=1, pub_num_spins
       call data_functions_alloc(evec(is),num_ngwfs,iscmplx=loc_cmplx)
    end do
    allocate(S_invH(pub_num_spins),stat=ierr)
    call utils_alloc_check('kernel_num_bound_states','S_invH',ierr)
    allocate(PS(pub_num_spins),stat=ierr)
    call utils_alloc_check('kernel_num_bound_states','PS',ierr)

    call sparse_embed_create(temp_proj_ham,ham(1))
    call sparse_embed_create(cond_proj,inv_overlap)
    call sparse_embed_create(SP,overlap,cond_proj)

    do is=1,pub_num_spins
       call sparse_embed_create(S_invH(is),inv_overlap,ham(is))
       call sparse_embed_create(PS(is),cond_proj,overlap)
    end do

    call sparse_embed_create(HPS,ham(1),PS(1))
    call sparse_embed_create(SPS,ham(1))

    do is=1,pub_num_spins
       ! start by computing the largest eigenvalue of the Hamiltonian
       ! this fixes the shift to which states are projected
       call sparse_embed_product(S_invH(is),inv_overlap,ham(is))
       call sparse_embed_extremal_eigenvalue(S_invH(is),overlap,temp_val(is), &
            tol=0.0001_DP)
       temp_val(is)=temp_val(is)+0.2_DP
    end do
    shift = maxval(temp_val(:))

    do is=1,pub_num_spins

       ! create a projected Hamiltonian, with all valence states
       ! projected to the shift value.
       call sparse_embed_product(PS(is),denskern%kern%m(is,PUB_1K),overlap)
       call sparse_embed_product(SP,overlap,denskern%kern%m(is,PUB_1K))
       call sparse_embed_product(SPS,SP,overlap)
       call sparse_embed_product(HPS,ham(is),PS(is))
       call sparse_embed_product(cond_proj,SP,HPS)
       call sparse_embed_axpy(cond_proj,SPS,-1.0_DP*shift)
       call sparse_embed_copy(temp_proj_ham,ham(is))
       call sparse_embed_axpy(temp_proj_ham,cond_proj,-1.0_DP)
       call sparse_embed_product(S_invH(is),inv_overlap,temp_proj_ham)

       ! agrecocmplx: modified call to use FUNCTIONS type array
       call sparse_embed_extremal_eigenvalue(S_invH(is),overlap,eval(is),&
            tol=0.00001_DP, min_val=.true.,evec=evec(is))

    end do

    ! set the maximum eval we want to find. Either this is
    ! zero (if we are interested in bound states of an isolated
    ! system, or it is the value of the HOMO+energy_range.
    temp_val(:) = kernel_bandgap(denskern,ham,overlap,inv_overlap)
    if(present(energy_range)) then
       ! temp_val = estimated bandgap
       ! HOMO = LUMO-band_gap
       max_eval(:) = eval(:) - temp_val(:) + energy_range
    else
       max_eval(:) = 0.0_DP
    endif

    ! initialise prev eval to HOMO
    if(present(energy_gap)) then
       prev_eval(:) = eval(:) - temp_val(:)
    endif

    ! if there are no bound states or states in the range, leave.
    ! else, start counting states.

    do is=1,pub_num_spins
       if (max_eval(is)<eval(is)) then
          num_bound_states(is) = 0
       else
          num_bound_states(is) = 1
          do while (eval(is)<max_eval(is) .and. num_bound_states(is) <&
                &num_ngwfs-num_occ(is))

             ! create cond_proj as the outer product of the eigenvector
             ! agrecocmplx: distinguish between real and complex case
             if (loc_cmplx) then
                call sparse_embed_outer_product(cond_proj,evec(is)%z(:),evec(is)%z(:))
             else
                call sparse_embed_outer_product(cond_proj,evec(is)%d(:),evec(is)%d(:))
             end if

             ! project out this value from the hamiltonian
             temp_val(is) = -(shift-eval(is))
             call sparse_embed_product(PS(is),cond_proj, overlap)
             call sparse_embed_axpy(S_invH(is),PS(is),-1.0_DP*temp_val(is))

             ! if we have specified an energy gap, set previous eigenvalue
             if(present(energy_gap)) then
                prev_eval(:)=eval(:)
             endif

             ! calculate smallest eigenvalue
             ! agrecocmplx
             call sparse_embed_extremal_eigenvalue(S_invH(is),overlap,eval(is),&
                  tol=0.00001_DP, min_val=.true.,evec=evec(is))

             if(eval(is)<max_eval(is)) then
               num_bound_states(is) = num_bound_states(is) + 1
             endif
          enddo
       endif

       ! successfully counted the number of bound states. Now keep counting
       ! until an energy gap of a given size is found in the conduction
       ! space manifold
       if(num_bound_states(is)>0 .and. present(energy_gap)) then
          ! set temp_energy_gap=eval-prev_eval
          temp_energy_gap=eval(is)-prev_eval(is)
          ! loop until gap if found or no more conduction states are left
          do while (temp_energy_gap<energy_gap .and. num_extra_states<&
             &num_ngwfs-num_occ(is)-num_bound_states(is))

             ! increase counter if condition is fulfilled.
             if(temp_energy_gap<energy_gap) then
                num_extra_states = num_extra_states+1
             endif
             prev_eval(is)=eval(is)

             ! construct cond_proj as the outer product of the eigenvector
             ! agrecocmplx
             if (loc_cmplx) then
                call sparse_embed_outer_product(cond_proj,evec(is)%z(:),evec(is)%z(:))
             else
                call sparse_embed_outer_product(cond_proj,evec(is)%d(:),evec(is)%d(:))
             end if

             ! project out this value from the hamiltonian
             temp_val(is) = -(shift-eval(is))
             call sparse_embed_product(PS(is),cond_proj,overlap)
             call sparse_embed_axpy(S_invH(is),PS(is),-1.0_DP*temp_val(is))

             ! calculate smallest eigenvalue of new ham
             ! agrecocmplx
             call sparse_embed_extremal_eigenvalue(S_invH(is),overlap,eval(is),&
                 tol=0.00001_DP,min_val=.true.,evec=evec(is))

             ! set new temp_energy_gap
             temp_energy_gap = eval(is)-prev_eval(is)
          enddo
          ! only add num_extra_states to kernel_num_bound_states if the
          ! combined number is less than the number of conduction states
          ! in the system
          if(num_bound_states(is)+num_extra_states<num_ngwfs-num_occ(is)) then
             num_bound_states(is) = num_bound_states(is)+num_extra_states
          endif
       endif
    end do

    call sparse_embed_destroy(SPS)
    call sparse_embed_destroy(HPS)

    do is=pub_num_spins,1,-1
       call sparse_embed_destroy(PS(is))
       call sparse_embed_destroy(S_invH(is))
    end do

    call sparse_embed_destroy(SP)
    call sparse_embed_destroy(cond_proj)
    call sparse_embed_destroy(temp_proj_ham)

    deallocate(PS,stat=ierr)
    call utils_dealloc_check('kernel_num_bound_states','PS',ierr)
    deallocate(S_invH,stat=ierr)
    call utils_dealloc_check('kernel_num_bound_states','S_invH',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    do is=1, pub_num_spins
       call data_functions_dealloc(evec(is))
    end do
    deallocate(evec, stat=ierr)
    call utils_dealloc_check('kernel_num_bound_states', 'evec',ierr)

  end function kernel_num_bound_states

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  function kernel_bandgap(denskern, ham, overlap, &
      inv_overlap)
    !====================================================================!
    ! Function computes the approximate bandgap of a given hamiltonian   !
    ! by projecting out all valence states and all conduction states     !
    ! respectively and computing the minimum/maximum eigenvalue.         !
    !--------------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff, July 2015                              !
    ! Edited for embedding by Robert Charlton, 13/09/2017.               !
    !====================================================================!

    use rundat, only: pub_num_spins, pub_num_kpoints, PUB_1K
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
        sparse_embed_destroy,sparse_embed_copy, sparse_embed_axpy, &
        sparse_embed_extremal_eigenvalue, sparse_embed_product
    use utils, only: utils_assert

    ! Arguments
    real(kind=DP)                 :: kernel_bandgap(pub_num_spins)
    type(DKERN), intent(in)       :: denskern       ! Density kernel
    type(SPAM3_EMBED), intent(in) :: ham(:)         ! hamiltonian
    type(SPAM3_EMBED), intent(in) :: overlap        ! overlap matrix
    type(SPAM3_EMBED), intent(in) :: inv_overlap    ! inverse overlap matrix

    ! local variables
    integer :: is, nspins
    type(SPAM3_EMBED) :: temp_proj_ham
    type(SPAM3_EMBED) :: cond_proj, cond_dkn
    type(SPAM3_EMBED) :: PS, SPS, SP, HPS, S_invH
    real(kind=DP) :: shift
    real(kind=DP) :: eval_HOMO, eval_LUMO

    ! jme: KPOINTS_DANGER
    ! Some parts of this function enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Function kernel_bandgap not ready yet for more than one k-point.')

    ! establish nspins and initialise temporary value for bandgap
    kernel_bandgap(:) = 0.0_DP
    nspins=denskern%kern%num_spins

    ! allocate data
    call sparse_embed_create(temp_proj_ham,ham(1))
    call sparse_embed_create(cond_proj,inv_overlap)
    call sparse_embed_create(cond_dkn,inv_overlap)
    call sparse_embed_create(PS,cond_proj,overlap)
    call sparse_embed_create(SP,overlap,cond_proj)
    call sparse_embed_create(HPS,ham(1),PS)
    call sparse_embed_create(SPS,overlap,PS)
    call sparse_embed_create(S_invH,inv_overlap,ham(1))

    ! loop over spin channels
    do is=1, nspins
       ! first compute largest eval of the hamiltonian
       call sparse_embed_product(S_invH,inv_overlap,ham(is))
       call sparse_embed_extremal_eigenvalue(S_invH,overlap,shift,tol=0.00001_DP)

       shift=shift+0.1_DP

       ! compute projected hamiltonian where all valence
       ! states are projected out
       call sparse_embed_product(PS,denskern%kern%m(is,PUB_1K),overlap)
       call sparse_embed_product(SP,overlap,denskern%kern%m(is,PUB_1K))
       call sparse_embed_product(SPS,SP,overlap)
       call sparse_embed_product(HPS,ham(is),PS)
       call sparse_embed_axpy(HPS,SPS,-1.0_DP*shift)
       call sparse_embed_product(cond_proj,SP,HPS)
       call sparse_embed_copy(temp_proj_ham,ham(is))
       call sparse_embed_axpy(temp_proj_ham,cond_proj,-1.0_DP)
       call sparse_embed_product(S_invH,inv_overlap,temp_proj_ham)

       ! calculate minimum eigenvalue of projected Hamiltonian
       call sparse_embed_extremal_eigenvalue(S_invH,overlap,eval_LUMO,&
            tol=0.00001_DP, min_val=.true.)

       ! set shift to minimum value of original hamiltonain

       call sparse_embed_product(S_invH,inv_overlap,ham(is))
       call sparse_embed_extremal_eigenvalue(S_invH,overlap,shift,min_val=.true.,&
            tol=0.00001_DP)

       shift=shift-0.1_DP

       ! compute projected hamiltonian where all conduction
       ! states are projected lower than the minimum eigenvalue

       call sparse_embed_copy(cond_dkn,inv_overlap)
       call sparse_embed_axpy(cond_dkn,denskern%kern%m(is,PUB_1K),-1.0_DP)
       call sparse_embed_product(PS,cond_dkn,overlap)
       call sparse_embed_product(SP,overlap,cond_dkn)
       call sparse_embed_product(SPS,SP,overlap)
       call sparse_embed_product(HPS,ham(is),PS)
       call sparse_embed_axpy(HPS,SPS,-1.0_DP*shift)
       call sparse_embed_product(cond_proj,SP,HPS)
       call sparse_embed_copy(temp_proj_ham,ham(is))
       call sparse_embed_axpy(temp_proj_ham,cond_proj,-1.0_DP)

       call sparse_embed_product(S_invH,inv_overlap,temp_proj_ham)
       ! calculate maximum eigenvalue of projected Hamiltonian ---> HOMO
       call sparse_embed_extremal_eigenvalue(S_invH,overlap,eval_HOMO,tol=0.00001_DP)

       kernel_bandgap(is) = eval_LUMO-eval_HOMO
    enddo

    ! deallocate data
    call sparse_embed_destroy(PS)
    call sparse_embed_destroy(SPS)
    call sparse_embed_destroy(SP)
    call sparse_embed_destroy(HPS)
    call sparse_embed_destroy(cond_dkn)
    call sparse_embed_destroy(cond_proj)
    call sparse_embed_destroy(temp_proj_ham)
    call sparse_embed_destroy(S_invH)

  end function kernel_bandgap

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  real(kind=DP) function kernel_commutator(denskern, ham, overlap, &
       inv_overlap)
    !====================================================================!
    ! Wrapper function that call kernel_invariant_commutator or          !
    ! kernel_rms_commutator dependent on which measure has been chosen   !
    !--------------------------------------------------------------------!
    ! Written by Robert Bell, July 2014                                  !
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.     !
    ! Modified for embedding by Robert Charlton, August 2017.            !
    !====================================================================!

    use rundat, only: pub_rms_kernel_measure
    use sparse_embed, only: SPAM3_EMBED

    implicit none

    ! Arguments
    type(DKERN), intent(inout) :: denskern   ! Density kernel
    type(SPAM3_EMBED), intent(in) :: ham(:)        ! Hamiltonian
    type(SPAM3_EMBED), intent(in) :: overlap       ! Overlap matrix
    type(SPAM3_EMBED), intent(in) :: inv_overlap   ! Inverse overlap matrix

    ! internal
    logical :: temp_workspace

    ! Create temporary workspace if we need one
    if (.not.denskern%workspace_created) then
       temp_workspace=.true.
       call kernel_workspace_create(denskern, overlap)
    else
       temp_workspace=.false.
    endif

    ! get commutator
    if (pub_rms_kernel_measure) then
       kernel_commutator = kernel_rms_commutator(denskern,ham,overlap)
    else
       kernel_commutator = kernel_invariant_commutator(denskern,ham,overlap,inv_overlap)
    endif

    ! destroy temporary workspace
    if (temp_workspace) then
       call kernel_workspace_destroy(denskern)
    endif


  end function kernel_commutator


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  real(kind=DP) function kernel_invariant_commutator(denskern,ham,overlap, &
       inv_overlap)

    !====================================================================!
    ! This function returns the invariant second moment (a generalised   !
    ! Frobenius norm) of the commutator between the Hamiltonial matrix   !
    ! and the density kernel, renormalised as to measure per-electron.   !
    ! ||[K,H]||_F = sqrt[ Tr[ [K,H]^dagger [K,H] ]] / N_elec             !
    !             = sqrt[ -Tr[ C^2 ]] / Tr[ KS ], C = (SKHS - HKS)S^^.   !
    ! For exact S^^S = 1 and idempotency KSK = K we have that            !
    ! ||[K,H]||_F = sqrt[ 2 Tr[ HKH (S^^-K) ]], so measuring the overlap !
    ! of valence and conduction manifolds in the Hamiltonian metric.     !
    !--------------------------------------------------------------------!
    ! Written by David O'Regan, July 2013, based on kernel_rms_commutator!
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.     !
    ! Modified for embedding by Robert Charlton, August 2017.            !
    !====================================================================!

    use constants, only: DP
    use sparse_embed, only: SPAM3_EMBED, &
        sparse_embed_create, sparse_embed_destroy, sparse_embed_product, &
        sparse_embed_transpose, sparse_embed_axpy, sparse_embed_trace

    implicit none

    ! Arguments
    type(DKERN), intent(inout) :: denskern   ! Density kernel
    type(SPAM3_EMBED), intent(in) :: ham(:)        ! Hamiltonian
    type(SPAM3_EMBED), intent(in) :: overlap       ! Overlap matrix
    type(SPAM3_EMBED), intent(in) :: inv_overlap   ! Inverse overlap matrix

    ! Local variables
    type(SPAM3_EMBED) :: mixed_commutator
    type(SPAM3_EMBED) :: skh, hks
    type(SPAM3_EMBED) :: sk, kh
    integer :: nspins
    integer :: is, ik
    real(kind=DP) :: trks, invariant

    ! Establish number of spins
    nspins = denskern%kern%num_spins

    ! ndmh: construct matrices with HKH structure
    ! ndmh: (a destroyed matrix is fine for structure code use)
    ! rc2013: EMBED_FIX! not for embedding matrices!
    call sparse_embed_create(kh,denskern%kern%m(1,1),ham(1))
    call sparse_embed_create(skh,ham(1),kh)
    call sparse_embed_destroy(kh)
    call sparse_embed_create(hks,skh)
    call sparse_embed_create(mixed_commutator,skh,inv_overlap)

    ! ndmh: create temporary for SK
    call sparse_embed_create(sk,overlap,denskern%kern%m(1,1))

    ! ndmh: re-evaluate ks if required
    call kernel_validate_ks(denskern, overlap)

    ! ddor: Zero the buffer
    kernel_invariant_commutator = 0.0_DP

    ! KPOINTS_DANGER
    do ik = 1, denskern%kern%num_kpoints
       do is=1,nspins

          ! ndmh: faster way of calculating commutator
          ! Calculate S.K
          call sparse_embed_transpose(sk, denskern%ks%m(is,ik))

          ! Calculate (S.K).H and transpose to get HKS
          call sparse_embed_product(skh,sk,ham(is))
          call sparse_embed_transpose(hks,skh)

          ! ddor: Calculate the doubly covariant commutator SKH-HKS
          call sparse_embed_axpy(skh,hks,-1.0_DP)

          ! ddor: Calculate the mixed-index commutator C = (SKH-HKS)
          call sparse_embed_product(mixed_commutator,skh,inv_overlap)

          ! ddor: Compute occupancy of this spin channel
          call sparse_embed_trace(trks, denskern%ks%m(is,ik))

          ! ddor: The invariant per-electron second moment of the commutator
          !     : is given by sqrt{ - Tr[C^2] - (Tr[C] )^2 } / N
          !     : where Tr[C]=0 by antisymmetry of the commutator

          ! ddor: Trace the squared mixed-index commutator and renormalise
          !     : to give -Tr[C^2] / N^2
          ! ddor: Allow for case of empty density kernel for one spin
          if (trks .ne. 0.0_DP) then
             call sparse_embed_trace(invariant, mixed_commutator, mixed_commutator)
             invariant = -1.0_DP * invariant / (trks)**2.0_DP
          else
             invariant = 0.0_DP
          endif

          ! ddor: Add the square of the commutator measure for this spin to buffer
          kernel_invariant_commutator = kernel_invariant_commutator + invariant

       end do
    end do

    ! ddor: Divide over spins for a per-electron measure, and take the norm
    kernel_invariant_commutator = sqrt(kernel_invariant_commutator / &
          real(nspins, kind=DP))

    ! Destroy structures
    call sparse_embed_destroy(sk)
    call sparse_embed_destroy(mixed_commutator)
    call sparse_embed_destroy(hks)
    call sparse_embed_destroy(skh)

  end function kernel_invariant_commutator

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  real(kind=DP) function kernel_rms_commutator(denskern, ham, overlap)

    !====================================================================!
    ! This function returns the rms value of the commutator between the  !
    ! Hamiltonial matrix and the density kernel: HKS-SKH.                !
    !--------------------------------------------------------------------!
    ! Written by Peter Haynes, November 2004                             !
    ! Based on an original version by Arash Mostofi, July 2002           !
    ! Modifications for faster calculation of commutator by Nicholas     !
    ! Hine in July 2009 and November 2010.                               !
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.     !
    ! Modified for embedding matrices by Robert Charlton, 16/08/2018.    !
    !====================================================================!

    use constants, only: DP
    use sparse_embed, only: SPAM3_EMBED, &
        sparse_embed_create, sparse_embed_destroy, sparse_embed_product, &
        sparse_embed_transpose, sparse_embed_axpy, sparse_embed_trace, &
        sparse_embed_rms_element, sparse_embed_array_num_rows, &
        sparse_embed_num_element
    use sparse_array, only: sparse_array_num_rows

    implicit none

    ! Arguments
    type(DKERN), intent(inout) :: denskern   ! Density kernel
    type(SPAM3_EMBED), intent(in) :: ham(:)        ! Hamiltonian
    type(SPAM3_EMBED), intent(in) :: overlap       ! Overlap matrix

    ! Local variables
    type(SPAM3_EMBED) :: skh,hks
    type(SPAM3_EMBED) :: sk,kh
    integer :: nspins
    integer :: is
    real(kind=DP) :: rms

    ! Establish number of spins
    nspins = denskern%kern%num_spins

    ! ndmh: faster way of calculating commutator

    ! ndmh: construct matrices with HKH structure
    ! ndmh: (a destroyed matrix is fine for structure code use)
    ! rc2013: not for embedding matrices!
    call sparse_embed_create(kh,denskern%kern%m(1,1),ham(1))
    call sparse_embed_create(skh,ham(1),kh)
    call sparse_embed_destroy(kh)
    call sparse_embed_create(hks,skh)

    ! ndmh: create temporary for SK
    call sparse_embed_create(sk,overlap,denskern%kern%m(1,1))

    ! ndmh: re-evaluate ks if required
    call kernel_validate_ks(denskern, overlap)

    kernel_rms_commutator = 0.0_DP
    do is=1,nspins

       ! ndmh: faster way of calculating commutator
       ! Calculate S.K
       call sparse_embed_transpose(sk, denskern%ks%m(is,1))

       ! Calculate (S.K).H and transpose to get HKS
       call sparse_embed_product(skh,sk,ham(is))
       call sparse_embed_transpose(hks,skh)

       ! Calculate the commutator SKH-HKS
       call sparse_embed_axpy(skh,hks,-1.0_DP)

       rms = sparse_embed_rms_element(skh)

       kernel_rms_commutator = kernel_rms_commutator + rms * rms

    end do

    kernel_rms_commutator = sqrt(kernel_rms_commutator * &
         sparse_embed_num_element(skh)/real(sparse_embed_array_num_rows(denskern%kern) &
         *nspins, kind=DP))

    ! Destroy structures
    call sparse_embed_destroy(sk)
    call sparse_embed_destroy(hks)
    call sparse_embed_destroy(skh)

  end function kernel_rms_commutator

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_init_core_ham(denskern, &  ! input-output
       coreham, overlap, n_occ)                ! input

    !======================================================================!
    ! This subroutine initialises the density matrix to the density        !
    ! that is obtained from the orbitals of the core-Hamiltonian matrix    !
    ! in row-indexed sparse matrix format.                                 !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 3/7/2001 for the ONES code.      !
    ! Adapted for ONETEP by Chris-Kriton Skylaris on 12/1/2004.            !
    ! Modified for parallel SPAM3 by Peter Haynes, July 2006               !
    ! Moved to kernel_mod by Nicholas Hine, August 2008                    !
    ! Reduced memory usage by re-using buffer, Nicholas Hine, October 2009 !
    ! Re-written for general dense matrices by Nicholas Hine, February 2010!
    ! Adapted for embedding by Joseph Prentice, May 2018                   !
    !======================================================================!

    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_convert, dense_destroy, &
         dense_eigensolve, dense_product
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_num_rows, &
         sparse_embed_scale
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none


    type(SPAM3_EMBED), intent(inout) :: denskern
    type(SPAM3_EMBED), intent(in)    :: coreham
    type(SPAM3_EMBED), intent(in)    :: overlap
    integer, intent(in)       :: n_occ


    ! cks: internal declarations
    integer :: ierr
    integer :: num
    type(DEM) :: coreham_dens, overlap_dens, denskern_dens
    type(DEM) :: eigenvecs_dens
    real(kind=DP), allocatable :: eigenvalues(:)
    ! agrecocmplx
    logical :: loc_cmplx
    character :: opB_loc

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering kernel_init_core_ham'

    ! Do nothing if no states are required
    if (n_occ==0) then
       call sparse_embed_scale(denskern,0.0_DP)
       return
    end if

    ! agrecocmplx
    loc_cmplx = denskern%iscmplx

    ! agrecocmplx
    if (loc_cmplx) then
       opB_loc = 'C'
    else
       opB_loc = 'T'
    end if

    ! ndmh: allocate storage for list of eigenvalues
    num = sparse_embed_num_rows(denskern)
    allocate(eigenvalues(num),stat=ierr)
    call utils_alloc_check('kernel_init_core_ham','eigenvalues',ierr)

    ! ndmh: create temporary dense matrices for eigenvecs, Ham and overlap
    ! agrecocmplx: initialise as complex if using complex NGWFs
    call dense_create(eigenvecs_dens,num,num,iscmplx=loc_cmplx)!n_occ)
    call dense_create(coreham_dens,num,num,iscmplx=loc_cmplx)
    call dense_create(overlap_dens,num,num,iscmplx=loc_cmplx)

    ! ndmh: convert sparse Ham and overlap to dense
    call dense_convert(coreham_dens,coreham)
    call dense_convert(overlap_dens,overlap)

    ! ndmh: solve the generalised eigenvalue problem.
    ! ndmh: On return from this subroutine eigenvecs_dens will contain
    !       the eigenvectors of coreham_dens
    call dense_eigensolve(n_occ,eigenvalues,coreham_dens,overlap_dens,1, &
         eigenvecs_dens)

    ! ndmh: remove temporary matrices that are no longer required
    call dense_destroy(overlap_dens)
    call dense_destroy(coreham_dens)

    ! ndmh: create density kernel from eigenvectors by summing
    ! ndmh: first n_occ eigenvectors multiplied by transposed matrix
    ! agrecocmplx: need to take hermitian of matrix B in complex case?
    ! I think so....
    call dense_create(denskern_dens,num,num,iscmplx=loc_cmplx)
    call dense_product(denskern_dens,eigenvecs_dens,eigenvecs_dens, &
         opA='N',opB=opB_loc,first_k=1,last_k=n_occ)

    ! ndmh: convert dense version of denskern back to sparse
    call dense_convert(denskern,denskern_dens)

    ! ndmh: deallocate temporary memory
    call dense_destroy(denskern_dens)
    call dense_destroy(eigenvecs_dens)
    deallocate(eigenvalues,stat=ierr)
    call utils_dealloc_check('kernel_init_core_ham','eigenvalues',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving kernel_init_core_ham'

  end subroutine kernel_init_core_ham

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_occ_check(denskern, kernel_reset, &  ! input-output
       ham, overlap, inv_overlap, n_occ)                       ! input

    !======================================================================!
    ! This subroutine compares an existing kernel to the canonical kernel  !
    ! obtained from a Palser-Manolopoulos purification for a supplied      !
    ! Hamiltonian. It takes the difference of the two matrices and finds   !
    ! the extremal eigenvalues of this difference: if they are greater     !
    ! than a predefined threshold, this indicates a flipped occupancy, so  !
    ! the old kernel is replaced with the new one. Strictly O(N) for       !
    ! truncated kernels.                                                   !
    !----------------------------------------------------------------------!
    ! Written by Nicholas Hine on 22/03/2013.                              !
    ! Modified for embedding by Joseph Prentice, June 2018                 !
    !======================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use palser_mano, only: palser_mano_kernel_optimise
    use rundat, only: pub_maxit_palser_mano
    use sparse_embed, only: SPAM3_EMBED, SPAM3_EMBED_ARRAY, sparse_embed_create,&
         sparse_embed_destroy, sparse_embed_copy, sparse_embed_axpy, &
         sparse_embed_product, sparse_embed_scale, sparse_embed_array_create, &
         sparse_embed_array_copy, sparse_embed_array_axpy, &
         sparse_embed_array_product, sparse_embed_array_extremal_eigenvalues, &
         sparse_embed_array_scale, sparse_embed_array_destroy

    implicit none

    ! Arguments
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    type(SPAM3_EMBED), intent(in)    :: ham(:)
    type(SPAM3_EMBED), intent(in)    :: overlap
    type(SPAM3_EMBED), intent(in)    :: inv_overlap
    integer, intent(in)        :: n_occ(:,:)
    logical, intent(inout)     :: kernel_reset

    ! Local Variables
    integer :: is, ik
    type(SPAM3_EMBED_ARRAY) :: denskern_new, diffkern, ks
    real(kind=DP), &
         dimension(denskern%num_spins, denskern%num_kpoints) :: u_bound, l_bound

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering kernel_occ_check'

    ! Allocate storage for canonical density kernel and difference
    call sparse_embed_array_create(denskern_new, denskern)
    call sparse_embed_array_create(diffkern, denskern)

    ! Perform Palser-Manolopoulos Canonical Purification to generate
    ! Canonical density kernel
    ! KPOINTS_DANGER
    do ik = 1, denskern%num_kpoints
       do is = 1, denskern%num_spins
          if (pub_maxit_palser_mano < 1) then
             ! Diagonalise the initial Hamiltonian
             call kernel_init_core_ham(denskern_new%m(is,ik), &
                  ham(is), overlap, n_occ(is,ik))
          else
             ! Palser-Manolopoulos Canonical Purification
             call palser_mano_kernel_optimise(denskern_new%m(is,ik), &
                  ham(is), overlap, inv_overlap, n_occ(is,ik), &
                  num_iter=pub_maxit_palser_mano)
          endif
       end do
    end do

    ! Form difference kernel between denskern and denskern_new
    call sparse_embed_array_copy(diffkern, denskern_new)
    call sparse_embed_array_axpy(diffkern, denskern, -1.0_DP)

    ! ndmh: create temporary matrix for KSKS structures
    call sparse_embed_array_create(ks, diffkern, overlap)
    ! ndmh: find upper bound of occupancies
    call sparse_embed_array_product(ks, diffkern, overlap)
    call sparse_embed_array_extremal_eigenvalues(ks, overlap, u_bound)
    ! ndmh: find lower bound of occupancies
    call sparse_embed_array_scale(ks, -1.0_DP)
    call sparse_embed_array_extremal_eigenvalues(ks, overlap, l_bound)
    l_bound = - l_bound

    ! ndmh: display and check occupancy bounds
    ! KPOINTS_DANGER: info on k-point and spin should probably be printed
    do ik = 1, denskern%num_kpoints
       do is = 1, denskern%num_spins
          if (pub_on_root) write(stdout,'(a,2f10.6)') &
               'Eigenvalues of difference kernel: ', u_bound(is,ik), l_bound(is,ik)
          if ((u_bound(is,ik) > 0.85_DP).or.(l_bound(is,ik) < -0.85_DP)) then
             if (pub_on_root) write(stdout,'(a)') 'Resetting kernel using &
                  &Palser-Manolopoulos canonically purified kernel guess'
             call sparse_embed_copy(denskern%m(is,ik), denskern_new%m(is,ik))
             kernel_reset = .true.
          end if
       end do
    end do

    ! Deallocate storage for canonical density kernel and difference
    call sparse_embed_array_destroy(ks)
    call sparse_embed_array_destroy(diffkern)
    call sparse_embed_array_destroy(denskern_new)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving kernel_occ_check'

  end subroutine kernel_occ_check

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////


  subroutine kernel_reset_occupancies(denskern, &  ! input-output
       overlap, n_occ)                ! input

    !====================================================================!
    ! This subroutine manually resets the occupancies of the density     !
    ! kernel using diagonalisation.                                      !
    !--------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2009.                     !
    ! Adapted for DKERN density kernels by J.M. Escartin, Fall 2014.     !
    ! Modified for embedding by Robert Charlton, 16/08/2018.             !
    !====================================================================!

    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_convert, dense_destroy, &
         dense_eigensolve, dense_product, dense_scale
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_array_num_rows
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none


    type(DKERN), intent(inout) :: denskern
    type(SPAM3_EMBED), intent(in)    :: overlap
    integer, intent(in)        :: n_occ(:,:)

    ! cks: internal declarations
    integer :: ierr
    integer :: num
    integer :: is, ik    ! spin and k-point counter
    integer :: nspins
    real(kind=DP), dimension(:), allocatable :: eigenvalues
    type(DEM) :: denskern_dens, overlap_dens, eigs_dens
    ! agrecocmplx
    logical :: loc_cmplx
    character :: opB_loc

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &kernel_reset_occupancies'

    ! agrecocmplx
    loc_cmplx = overlap%iscmplx

    ! agrecocmplx
    if (loc_cmplx) then
       opB_loc = 'C'
    else
       opB_loc = 'T'
    end if

    ! Establish number of spins
    nspins = denskern%kern%num_spins

    num = int( sparse_embed_array_num_rows(denskern%kern) )

    allocate(eigenvalues(num),stat=ierr)
    call utils_alloc_check('kernel_reset_occupancies','eigenvalues',ierr)

    ! ndmh: create temporary dense matrices
    ! agrecocmplx
    call dense_create(denskern_dens,num,num,iscmplx=loc_cmplx)
    call dense_create(overlap_dens,num,num,iscmplx=loc_cmplx)
    call dense_create(eigs_dens,num,num,iscmplx=loc_cmplx)

    do ik = 1, denskern%kern%num_kpoints

       ! cks: loop over spins
       do is=1,nspins

          ! ndmh: expand the density kernel to a dense matrix
          call dense_convert(denskern_dens,denskern%kern%m(is,ik))

          ! ndmh: multiply by -1 so that we get the highest occupancy eigenvalues
          call dense_scale(denskern_dens,-1.0_DP)

          ! ndmh: expand the overlap matrix to a dense matrix
          call dense_convert(overlap_dens,overlap)

          ! ndmh: solve the generalised eigenvalue problem
          call dense_eigensolve(n_occ(is,ik),eigenvalues,denskern_dens,overlap_dens, &
               2,eigs_dens)

          ! ndmh: construct a density kernel from the first n_occ(is,ik) eigenvectors
          ! agrecocmplx: need to take hermitian of matrix B in complex case?
          ! I think so....
          call dense_product(denskern_dens,eigs_dens,eigs_dens, &
               opB=opB_loc,first_k=1,last_k=n_occ(is,ik))

          ! ndmh: convert back to sparse matrix
          call dense_convert(denskern%kern%m(is,ik),denskern_dens)

       end do

    end do

    ! ndmh: clean up temporary matrices
    call dense_destroy(eigs_dens)
    call dense_destroy(overlap_dens)
    call dense_destroy(denskern_dens)

    deallocate(eigenvalues,stat=ierr)
    call utils_dealloc_check('kernel_reset_occupancies','eigenvalues',ierr)

    ! cks: denskern has changed: invalidate ks and ksk workspaces
    call kernel_workspace_invalidate(denskern)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &kernel_reset_occupancies'

  end subroutine kernel_reset_occupancies

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_from_vecs_asc(density_sq, &
       orbitals,num,n_occ)

    !================================================================!
    ! This subroutine constructs an idempotent density kernel from a !
    ! set of eigenvectors of the density kernel (assumes that the    !
    ! eigenvectors are ordered in ascending order of eigenvalue      !
    ! (occupancy).                                                   !
    !----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2009.                 !
    !================================================================!

    use constants, only: DP
    use utils, only: utils_assert

    implicit none

    integer, intent(in) :: num, n_occ
    real(kind=DP), intent(out) :: density_sq(num,num)
    real(kind=DP), intent(in) :: orbitals(num,num)

    ! cks: internal variable declarations
    integer :: row, col, occ_count
    real(kind=DP) :: el

    call utils_assert(n_occ <= num, 'Error in kernel_from_vecs_asc: &
            &number of occupied orbitals exceeds number of NGWFs')

    ! cks: initialise
    density_sq = 0.0_DP

    do col=1,num
       do row=1,num

          ! cks: loop over the last n_occ (highest occupancy) orbitals
          el = 0.0_DP
          do occ_count = num, num-n_occ +1, -1
             el = el + orbitals(row,occ_count) * orbitals(col,occ_count)
          end do

          density_sq(row,col) = el

       end do
    end do


  end subroutine kernel_from_vecs_asc

  subroutine kernel_from_vecs_asc2(density_sq, &
       orbitals,num,n_occ)

    !================================================================!
    ! This subroutine constructs an idempotent density kernel from a !
    ! set of eigenvectors of the density kernel (assumes that the    !
    ! eigenvectors are ordered in ascending order of eigenvalue      !
    ! (occupancy).                                                   !
    !----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2009.                 !
    !================================================================!

    use constants, only: DP
    use utils, only: utils_assert

    implicit none

    integer, intent(in) :: num, n_occ
    real(kind=DP), intent(out) :: density_sq(num,num)
    real(kind=DP), intent(in) :: orbitals(num,num)

    ! cks: internal variable declarations
    integer :: row, col, occ_count
    real(kind=DP) :: el

    call utils_assert(n_occ <= num, 'Error in kernel_from_vecs_asc2: &
            &number of occupied orbitals exceeds number of NGWFs')

    ! cks: initialise
    density_sq = 0.0_DP

    do col=1,num
       do row=1,num

          ! cks: loop over the last n_occ (highest occupancy) orbitals
          el = 0.0_DP
          do occ_count = 1, n_occ
             el = el + orbitals(row,occ_count) * orbitals(col,occ_count)
          end do

          density_sq(row,col) = el

       end do
    end do


  end subroutine kernel_from_vecs_asc2


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

#if 0
  subroutine kernel_christoffel(ngwf_basis, &
       new_ngwfs_on_grid, old_ngwfs_on_grid, new_denskern, &
       old_denskern, new_inv_overlap, old_inv_overlap, &
       new_sp_overlap, old_sp_overlap)

    !==============================================================!
    ! This routine updates the density kernel with terms from the  !
    ! Christoffel symbols which preserve the completeness of the   !
    ! basis.                                                       !
    !--------------------------------------------------------------!
    ! Written by David O'Regan in June 2010                        !
    ! Modifications by Nicholas Hine, November 2010.               !
    ! PAW added by David O'Regan 29/3/11                           !
    ! Moved from ngwf_cg_mod and modified accordingly by           !
    ! Simon Dubois 20/06/11                                        !
    !==============================================================!

    use augmentation, only: augmentation_overlap
    use datatypes, only: data_functions_copy, data_functions_axpy
    use bibliography, only: bibliography_cite
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_brappd_ketppd
    use rundat, only: pub_aug
    use sparse, only: SPAM3, sparse_create,  sparse_transpose, &
         sparse_product, sparse_destroy, sparse_axpy, sparse_copy, &
         sparse_scale
    use sparse_array, only: SPAM3_ARRAY, sparse_array_create, &
         sparse_array_transpose, sparse_array_product, sparse_array_destroy &
         sparse_array_axpy, sparse_scale
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: new_ngwfs_on_grid(ngwf_basis%size_on_grid)
    real(kind=DP), intent(in) :: old_ngwfs_on_grid(ngwf_basis%size_on_grid)
    type(SPAM3_ARRAY), intent(inout) :: new_denskern
    type(SPAM3_ARRAY), intent(in) :: old_denskern
    type(SPAM3), intent(in) :: new_inv_overlap
    type(SPAM3), intent(in) :: old_inv_overlap
    type(SPAM3), intent(in) :: new_sp_overlap
    type(SPAM3), intent(in) :: old_sp_overlap

    ! Local Variables
    type(SPAM3) :: step, tc_step, step_sp_overlap
    type(SPAM3_ARRAY) :: ktmp
    real(kind=DP), allocatable :: step_on_grid(:)
    logical, parameter :: using_old_inv_overlap = .true.
    integer :: ierr, is

    call bibliography_cite('CHRISTOFFEL')

    allocate(step_on_grid(ngwf_basis%size_on_grid),stat=ierr)
    call utils_alloc_check('kernel_christoffel',&
         'step_on_grid',ierr)

    ! Create temporary matrix structures
    step%structure = 'S'
    call sparse_create(step)
    call sparse_create(tc_step,step,new_inv_overlap)
    call sparse_array_create(ktmp, old_denskern)

    ! Calculate change |dphi> in NGWFs

    ! ndmh_pointerfun
    call data_functions_copy(step_on_grid,new_ngwfs_on_grid)
    call data_functions_axpy(step_on_grid,old_ngwfs_on_grid,-1.0_DP)
    !step_on_grid = new_ngwfs_on_grid - old_ngwfs_on_grid

    ! Calculate overlap of the change with the old NGWFs
    call integrals_brappd_ketppd(step, &
         step_on_grid,ngwf_basis,old_ngwfs_on_grid,ngwf_basis)

    ! ndmh: Apply the augmentated overlap operator to this matrix
    if (pub_aug) then

       ! ndmh: <proj_i|dphi_a> = <proj_i|new_phi_a> - <proj_i|old_phi_a>
       call sparse_create(step_sp_overlap,old_sp_overlap)
       call sparse_copy(step_sp_overlap,new_sp_overlap)
       call sparse_axpy(step_sp_overlap,old_sp_overlap,-1.0_DP)

       ! ndmh: Calculate the augmentation of the "step" matrix
       call augmentation_overlap(step,paw_sp,old_sp_overlap,step_sp_overlap)

       ! ndmh: Clean up <proj_i|dphi_a> (no longer needed)
       call sparse_destroy(step_sp_overlap)

    end if

    ! Raise an index in the step-NGWF overlap
    if (using_old_inv_overlap) then
       ! ddor: Ideally, the inverse overlap for the point at
       !       which the gradient is computed is used.
       ! ddor: Here, we have the old NGWFs, but translated
       !       to their new positions.
       call sparse_product(tc_step,step,old_inv_overlap)
    else
       call sparse_product(tc_step,step,new_inv_overlap)
    endif

    ! Compute -K <g|phi> S^^
    call sparse_array_product(new_denskern, old_denskern, tc_step)
    call sparse_array_scale(new_denskern, -1.0_DP)
    ! Transpose to get -S^^ <phi|g> K
    call sparse_array_transpose(ktmp, new_denskern)
    ! Add this to what we got before
    call sparse_array_axpy(new_denskern, ktmp, 1.0_DP)
    ! Add the previous kernel to the calculated change to get new kernel
    call sparse_array_axpy(new_denskern, old_denskern, 1.0_DP)
    enddo

    ! Destroy temporary matrix structures
    call sparse_array_destroy(ktmp)
    call sparse_destroy(tc_step)
    call sparse_destroy(step)

    deallocate(step_on_grid,stat=ierr)
    call utils_dealloc_check('kernel_christoffel',&
         &'step_on_grid',ierr)

  end subroutine kernel_christoffel

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_basis_transform(ngwf_basis, &
       new_ngwfs_on_grid, old_ngwfs_on_grid, new_denskern, &
       old_denskern, new_inv_overlap, new_sp_overlap, &
       old_sp_overlap)

    !==============================================================!
    ! Transforms the density kernel to its representation in terms !
    ! of new NGWFs                                                 !
    !--------------------------------------------------------------!
    ! Written by Nicholas Hine on 22/10/09                         !
    ! PAW added by David O'Regan 29/3/11                           !
    ! Moved from ngwf_cg_mod and modified accordingly by           !
    ! Simon Dubois 20/06/11                                        !
    !==============================================================!

    use augmentation, only: augmentation_overlap
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_brappd_ketppd
    use rundat, only: pub_aug
    use sparse, only: SPAM3, sparse_create,  sparse_transpose, &
         sparse_product, sparse_destroy, sparse_copy, sparse_axpy
    use sparse_array, only: sparse_array_create, sparse_array_product, &
         sparse_array_destroy

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: new_ngwfs_on_grid(ngwf_basis%size_on_grid)
    real(kind=DP), intent(in) :: old_ngwfs_on_grid(ngwf_basis%size_on_grid)
    type(SPAM3), intent(in) :: new_inv_overlap
    type(SPAM3), intent(in) :: new_sp_overlap
    type(SPAM3), intent(in) :: old_sp_overlap
    type(SPAM3_ARRAY), intent(inout) :: new_denskern
    type(SPAM3_ARRAY), intent(in) :: old_denskern

    ! Local Variables
    type(SPAM3) :: overlap_old_new
    type(SPAM3) :: sk,ks
    type(SPAM3_ARRAY) :: ksk
    integer :: is

    ! Create temporary matrix structures
    overlap_old_new%structure = 'S'
    call sparse_create(overlap_old_new)
    sk%structure = 'SK'
    call sparse_create(sk)
    ks%structure = 'KS'
    call sparse_create(ks)
    call sparse_array_create(ksk, old_denskern, structure='KSK')

    ! Calculate overlap of old NGWFs with new NGWFs
    call integrals_brappd_ketppd(overlap_old_new, &
         old_ngwfs_on_grid,ngwf_basis,new_ngwfs_on_grid,ngwf_basis)

    if (pub_aug) then
       ! ddor: Calculate the augmentation of the overlap matrix
       call augmentation_overlap(overlap_old_new,old_sp_overlap, &
            new_sp_overlap,paw_sp)
    end if

    ! Calculate sk_a^b = <f_a|f'_c>.(S^-1)^cb
    call sparse_product(sk,overlap_old_new,new_inv_overlap)
    ! Transpose to find ks^a_b = (S^-1)^ac.<f'_c|f_b>
    call sparse_transpose(ks,sk)

    ! Calculate density kernel transformed to new NGWF basis
    ! K'^ah = (S^-1)^ab.<f'_b|f_g> K^ge <f_e|f'_z>.(S^-1)^zh
    call sparse_array_product(ksk, ks, old_denskern)
    call sparse_array_product(new_denskern, ksk, sk)

    ! Destroy temporary matrix structures
    call sparse_array_destroy(ksk)
    call sparse_destroy(ks)
    call sparse_destroy(sk)
    call sparse_destroy(overlap_old_new)

  end subroutine kernel_basis_transform
#endif


  subroutine kernel_build_subsystem(kern_dst,kern_src)

    !=========================================================================!
    ! This subroutine copies the components of a DKERN type object      !
    ! to a subsystem matrix.                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  kern_src (inout) : Container for matrices to be copied into.           !
    !  kern_dst (inout) : Container for matrices to be copied from.           !
    !-------------------------------------------------------------------------!
    ! Adapted from existing code by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use constants, only: stdout
    USE sparse_embed, only : sparse_embed_array_create_sub,sparse_embed_array_destroy

    implicit none

    type(DKERN), intent(inout) :: kern_dst(:)
    type(DKERN), intent(in) :: kern_src

    ! Local variables
    integer :: isub

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering kernel_build_subsystem'


    do isub=1,kern_src%kern%m(1,1)%mrows

       kern_dst(isub)%core_struc(1)=trim(kern_src%core_struc(1))

       call sparse_embed_array_create_sub(kern_dst(isub)%kern, &
            kern_src%kern, row=isub, col=isub)

    end do

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving kernel_build_subsystem'


  end subroutine kernel_build_subsystem

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////


  subroutine kernel_basis_update(new_denskern, old_denskern, &
       new_overlap, old_overlap, new_inv_overlap)

    !==============================================================!
    ! Update the density kernel such that                          !
    ! K_new*S_new = K_old*S_old                                    !
    !--------------------------------------------------------------!
    ! Written by Simon Dubois on 20/06/11                          !
    !==============================================================!

    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_create,  &
         sparse_product, sparse_destroy, sparse_copy, sparse_axpy

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: new_denskern(pub_num_spins)
    type(SPAM3), intent(in) :: old_denskern(pub_num_spins)
    type(SPAM3), intent(in) :: new_overlap
    type(SPAM3), intent(in) :: old_overlap
    type(SPAM3), intent(in) :: new_inv_overlap

    ! Local Variables
    type(SPAM3) :: ds,ks,ksk
    integer :: is

    ! Create temporary matrix structures
    ds%structure = 'S'
    call sparse_create(ds)
    ks%structure = 'KS'
    call sparse_create(ks)
    ksk%structure = 'KSK'
    call sparse_create(ksk)

    ! Calculate ds
    call sparse_copy(ds,new_overlap)
    call sparse_axpy(ds,old_overlap,-1.0_DP)

    ! Calculate ks, ksk and new_denskern
    do is=1,pub_num_spins
       call sparse_product(ks,old_denskern(is),ds)
       call sparse_product(ksk,ks,new_inv_overlap)
       call sparse_copy(new_denskern(is),old_denskern(is))
       call sparse_axpy(new_denskern(is),ksk,-1.0_DP)
    end do

    ! Destroy temporary matrix structures
    call sparse_destroy(ksk)
    call sparse_destroy(ks)
    call sparse_destroy(ds)

  end subroutine kernel_basis_update


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////




end module kernel
