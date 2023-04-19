!================================================================!
!                                                                !
!          Spherical wave resolution of identity module          !
!                                                                !
! This module contains datastructures and subroutines for the    !
! resolution of identity approach based on truncated spherical   !
! waves.                                                         !
!                                                                !
! Everything that pertains to the metric matrices, spherical     !
! waves, spherical wave potentials (collectively called "SWOPs", !
! for "Spherical Wave Or Potential (thereof)") lives here.       !
! The only significant exception is the SWOP cache, which does   !
! belong here conceptually, but is part of the sw_expansion      !
! module for reasons of efficiency.                              !
!                                                                !
! Everything that depends on NGWFs (expansion of their products) !
! should go into sw_expansion_mod.                               !
!                                                                !
!----------------------------------------------------------------!
! - The original module (as hf_exchange_mod) was written by      !
!   Quintin Hill in 2008/9 with supervision by Chris-Kriton      !
!   Skylaris, then extensively modified by Jacek Dziedzic in     !
!   2012 and 2013.                                               !
! - sw_expansion_mod was then extracted in April 2014 from       !
!   hf_exchange_mod by Jacek Dziedzic, in preparation for DMA.   !
! - This module was then extracted in February 2015 from         !
!   sw_expansion_mod by Jacek Dziedzic.                          !
! - The approach to metric matrix calculation was reworked by    !
!   Jacek Dziedzic in Feb-Apr 2016 to support restarts, dynamic  !
!   load balancing, proximity sorting, and more fine-grained OMP.!
!----------------------------------------------------------------!
!================================================================!

! TODO:
! - Periodicity broken by MIC not used in Chebyshevs.
! - Periodicity broken in multiple other places, look for "@MIC"
! - Choose direction on whether to support non-identical NGWF radii *within a
!   single SW_RI* (likely not). If not, remove 'init_centre' and associated
!   cruft. Places where we rely on the identical-radii assumption are denoted
!   with @IDENTICAL_RADII. There might be some where I'm oblivious to it.
! - 'Parent' SW_RI's are not supported, I can't even remember what the purpose
!   of this might have been.
! - Enable overlap metric (O) to be evaluated using 2Dn-1Da scheme --- currently
!   this is only implemented for the electrostatic metric (V). The string
!   @2DN_1DA_OVERLAP indicates relevant locations in the module.

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

module sw_resolution_of_identity

  use constants, only: DP, CRLF, stdout, max_rundat_swri_block_entries
  use geometry, only: POINT
  use hash_table, only: HT_HASH_TABLE
  use neighbour_list, only: NL_NEIGHBOUR_LIST
  use sparse, only: SPAM3
  use utils, only: utils_trace_in, utils_trace_out

  implicit none

  private

  ! --------------------------------
  ! --- Major public subroutines ---
  ! --------------------------------

  public :: swri_init_module_stage_0
  public :: swri_init_module_stage_1
  public :: swri_init_module_stage_2
  public :: swri_cleanup_module

  public :: swri_init_container_stage_1
  public :: swri_init_container_stage_2
  public :: swri_destroy_container
  public :: swri_get_handle_to

  ! --------------------------------
  ! --- Minor public subroutines ---
  ! --------------------------------

  public :: swri_expansion_centres
  public :: swri_num_sws_on_centre
  public :: swri_init_centre
  public :: swri_extract_matrix_2
  public :: swri_build_matrix_2
  public :: swri_make_wmatrix_2
  public :: swri_sw_to_ablmq
  public :: swri_sph_bess_solidharm_int
  public :: swri_obtain_swops_in_ppd
  public :: swri_obtain_swops_in_ppd_ptr
  public :: swri_obtain_swops_in_ppd_ptr_fast
  public :: swri_obtain_swops_in_ppd_set_omp
  public :: swri_obtain_swops_at_point
  public :: swri_swop_calc_all_in_ppd
  public :: swri_calc_sws_at_point
  public :: swri_calc_swpots_at_point
  public :: swri_calc_swops_at_point
  public :: swri_qualities_equal

  ! ----------------------------------
  ! --- Public types and variables ---
  ! ----------------------------------

  public :: BESSEL
  public :: ATOM_CENTRE
  public :: PBC_IMAGE_INFO
  public :: SW_RI
  public :: SW_QUALITY

  public :: swri_library, swri_library_size
  public :: swri_rundat_block_data ! Accessed from rundat_blocks
  ! --------------------------------------------------------------------------
  ! Parameters of a Bessel function
  ! --------------------------------------------------------------------------
  type BESSEL
     integer        :: lval   ! Quantum number l
     integer        :: qidx   ! Ordinal number of bessel zero corresp. to qval
     real(kind=DP)  :: qval   ! Set so that j_l(qa) = 0
     real(kind=DP)  :: aval   ! Cut off radius of spherical Bessel function
     real(kind=DP)  :: farpotint ! Radial potential integral when r > a
     real(kind=DP)  :: nearpotint ! Constant part of radial potential when r < a
  end type BESSEL
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  ! An atom centre in an SW expansion
  ! --------------------------------------------------------------------------
  type ATOM_CENTRE
     integer       :: species_number
     integer       :: global_idx
     real(kind=DP) :: radius         ! jd: Localization radius of the NGWFs
     type(POINT)   :: incell         ! Centre of spherical wave in cell
  end type ATOM_CENTRE

  ! --------------------------------------------------------------------------
  ! Information on a PBC image (for debugging only)
  ! --------------------------------------------------------------------------
  type PBC_IMAGE_INFO
     integer       :: image    ! jd: Integer -13..13 specifying which image
     type(POINT)   :: curpoint ! jd: Point where the SWOP was calculated
     type(POINT)   :: centre   ! jd: Where the atom centre was *originally*
     type(POINT)   :: disp     ! jd: Final obtained MIC displacement vector
     type(POINT)   :: disp0    ! jd: Original (pre-MIC) displacement vector
     real(kind=DP) :: rad      ! jd: Final obtained radius from image-of-centre to curpoint
     real(kind=DP) :: rad0     ! jd: Original distance from centre to curpoint
  end type PBC_IMAGE_INFO

  ! --------------------------------------------------------------------------
  ! Specification of auxiliary basis set quality
  ! --------------------------------------------------------------------------
  type SW_QUALITY
     integer       :: min_l                     ! jd: Any changes to this
     integer       :: max_l                     !     type should be reflected
     integer       :: max_q                     !     in swri_qualities_equal().
     integer       :: num_sws_per_centre        !
     integer       :: max_sws_per_centre        !
  end type SW_QUALITY

  ! --------------------------------------------------------------------------
  ! jd: Spherical-wave based resolution of identity.
  !     Encapsulates the metric matrices, Bessel info, index xlat tables.
  !     This does *not* depend on the NGWFs or DKN.
  !     This does depend on the positions of the ions, on the NGWF radii,
  !     on max l and q of the auxiliary basis.
  ! --------------------------------------------------------------------------
  type SW_RI
     logical           :: created = .false.
     logical           :: initialised = .false.
     logical           :: full_metric_matrices_destroyed = .false.
     character(len=32) :: swri_name = '(*I have not been set up yet!*) '
     type(SW_QUALITY)  :: quality
     integer     :: cheb_intervals  ! ] <- quality of the Chebyshev
     integer     :: cheb_order      ! ]    interpolation
     logical     :: store_blocks_on_disk
     logical     :: has_metric(2)           ! 1:V, 2:O
     type(SPAM3) :: full_metric_matrices(2) ! 1:V, 2:O
     integer     :: int_scheme     ! JCW: integral evaluation scheme (see below)
     integer     :: natoms

     ! jd: Bessel and SW book-keeping
     integer                    :: n_xbatches ! number of xbatches
     integer                    :: xbatchsize ! min(num_sws_per_centre,swri_batchsize_loc)
     integer                    :: max_num_bessels ! over all species
     type(BESSEL), allocatable  :: sphbessels(:,:) ! (species:bessel_idx)
     real(kind=DP), allocatable :: radtable(:) ! radius for each species
     integer, allocatable       :: sw_idx_to_bessel_idx(:) ! xlat sw->bessel
     integer, allocatable       :: sw_idx_to_m(:)          ! xlat sw->m
     integer, allocatable       :: sw_idx_to_l(:)          ! xlat sw->l
     integer, allocatable       :: sw_idx_to_q_idx(:)      ! xlat sw->q_idx
     ! JCW: Number of allowed q-values for each l: only SWs with spherical
     ! JCW: Bessel components with E = 0.5*q^2 < pub_cutoff_energy (psinc KE
     ! JCW: cutoff) are included in SWRI. This is per-l because the spherical
     ! JCW: Bessels for each l have different zeroes (q-values).
     integer, allocatable       :: num_q_per_l(:,:) ! (1:num_species,l_min:lmax)
     ! jd: Neighbour list for s-neighbours
     type(NL_NEIGHBOUR_LIST)    :: s_atoms_nl
     ! jd: Radial lookup for Bessels
     real(kind=DP), allocatable :: rad_lookup(:,:)

     ! jd: Hash tables and their file output control

     ! jd: Stores partial atomblocks while they are merged from tiles.
     !     Blocks that are fully merged are transferred to metric matrices,
     !     following consistency-forcing. This HT is thus empty after stage 2
     !     is complete.
     !     #1: atoma, #2: atomb, #3: xbatch, #4: SW_V or SW_O.
     type(HT_HASH_TABLE)        :: tiles_ht

     ! jd: SWOP caches, independent of DKN, NGWFs.
     type(HT_HASH_TABLE) :: swops_in_ppds_ht
     type(HT_HASH_TABLE) :: swops_at_points_ht

     ! jd: Periodic image map -- temporary functionality for debugging.
     !     Will be removed soon.
     integer, allocatable :: image_map(:,:,:) ! n_pts, n_ppds, nat

     integer                    :: hash_table_info_unit
     character(len=256)         :: hash_table_info_filename

  end type SW_RI
  ! --------------------------------------------------------------------------

  ! ---------------------------------------
  ! --- Public variables and parameters ---
  ! ---------------------------------------

  ! NB: Making the library un-protected and a target allows use association
  !     e.g. s_atoms_nl => swri_library(dma_swri_h)%s_atoms_nl
  type(SW_RI), allocatable, save, target :: swri_library(:)

  ! -----------------------------------------------------------------------
  ! ---                       P R I V A T E                             ---
  ! -----------------------------------------------------------------------

  ! -------------
  ! --- Types ---
  ! -------------

  type RUNDAT_BLOCK_ENTRY
     character(len=32) :: entry_name
     type(SW_QUALITY)  :: quality
     logical           :: has_metric(2)
     integer           :: cheb_intervals
     integer           :: cheb_order
     logical           :: read_from_file
     logical           :: write_to_file
     logical           :: print
     logical           :: disassemble
     logical           :: erase
     logical           :: quit
     integer           :: int_scheme
     character(len=32) :: parent_name
  end type RUNDAT_BLOCK_ENTRY

  ! --------------------------------------------------------------------------
  ! Describes a single off-centre atomblock
  ! --------------------------------------------------------------------------
  type ATOMBLOCK_DESCRIPTOR
     integer :: atoma
     integer :: atomb
     integer :: orig_atoma
     integer :: orig_atomb
     integer :: atoma_sparse_owner
  end type ATOMBLOCK_DESCRIPTOR

  ! --------------------------------------------------------------------------
  ! Describes a single off-centre atomblock tile
  ! --------------------------------------------------------------------------
  type TILE_DESCRIPTOR
     integer :: atomblock
     integer :: xbatch
     integer :: worker_proc
  end type TILE_DESCRIPTOR

  ! -----------------------------
  ! --- Module-wide variables ---
  ! -----------------------------

  integer, save, protected :: swri_library_size = -1
  type(RUNDAT_BLOCK_ENTRY), save, protected, &
       dimension(max_rundat_swri_block_entries) :: swri_rundat_block_data

  ! ------------------------------
  ! --- Module-wide parameters ---
  ! ------------------------------

  logical, parameter :: force_metric_matrix_consistency = .true.

  integer, parameter :: TILE_UNASSIGNED = -1
  integer, parameter :: TILE_COMPLETED = -2

  integer, parameter :: TAG_TILE_REQUEST = 1001
  integer, parameter :: TAG_ATOMBLOCK_ID = 1002
  integer, parameter :: TAG_XBATCH_ID = 1003
  integer, parameter :: TAG_TILE_ID = 1004
  integer, parameter :: TAG_TILE_DATA = 1005

  ! JCW: Integral evaluation scheme variants, to set value of
  ! JCW: RUNDAT_BLOCK_ENTRY%int_scheme (and subsequently SW_RI%int_scheme).
  ! JCW: This controls the integration scheme used for off-site metric matrix
  ! JCW: elements. Same-centre elements are always evaluated analytically.
  ! ----------------
  ! Undefined scheme
  ! ----------------
  ! If RUNDAT_BLOCK_ENTRY%int_scheme has this value, then we will select an
  ! appropriate default scheme during stage 1 of the initialization.
  integer, parameter :: SCHEME_UNDEFINED = -1
  ! --------------------------
  ! 3-D Chebyshev (3Dc) scheme
  ! --------------------------
  ! This is the original scheme described in section II.D.1 of
  !   J. Dziedzic et al.,J. Chem. Phys. 139, 214103 (2013)
  ! and is based on the expansion of truncated SWs and SWpots (in ONETEP's
  ! Cartesian coordinate system) in terms of Chebyshev polynomials in a 3-D
  ! sphere. The integral is evaluated by analytic integration over the
  ! products of the SWs/SWpots (depending on the metric) expanded in the
  ! Chebyshev polynomials.
  integer, parameter :: SCHEME_3DC       = 0
  ! ----------------------------------
  ! 2-D numerical, 1-D analytic scheme
  ! ----------------------------------
  ! In this scheme, off-site atomblocks of the V-matrix are evaluated in
  ! an "integral frame" coordinate system, in which the z-axis is aligned
  ! along the vector between the two atoms. The real-spherical harmonics (RSH)
  ! component of the SWs/SWpots are aligned along this z-axis.
  ! In this coordinate system, the metric matrix integrals can be separated into
  ! a 2-D numerical integral over the radial and polar coordinates of a
  ! spherical polar coordinate system centred on the atom in the "home sphere"
  ! and a 1-D analytic integral over the azimuthal coordinate.
  integer, parameter :: SCHEME_2DN_1DA       = 1

  ! Default integral evaluation scheme for each type of metric
  ! @2DN_1DA_OVERLAP
  integer, parameter :: SCHEME_DEFAULT(2)    = [ SCHEME_2DN_1DA, &
                                                 SCHEME_3DC ] ! 1:V, 2:O

  ! Min and max allowed scheme numbers (increment if more schemes are added)
  integer, parameter :: SCHEME_MIN           = 0
  integer, parameter :: SCHEME_MAX           = 1

  ! Per-scheme default max tile thickness for evaluation of atomblock
  ! (can be overidden by setting swri_cheb_batchsize keyword)
  integer, parameter :: DEFAULT_XBATCHSIZE(SCHEME_MIN:SCHEME_MAX) &
                                             = [ 12, HUGE(int(1)) ]
                                               ! 0: 3Dc
                                               ! 1: 2Dn-Da
  ! HUGE(int(1)) is the largest integer representable for the default integer
  ! kind. Since xbatchsize is set for each SW_RI instance in
  ! swri_init_container_stage_1 as
  !   swri%xbatchsize = min(n_sws, swri_cheb_batchsize_loc)
  ! this means that for 2Dn-1Da, xbatchsize defaults to n_sws (i.e. no tiling).

  ! Performance testing for the 2Dn-1Da scheme with variety of system sizes
  ! indicates that using smaller batches has no performance advantage, and
  ! for larger system sizes and numbers of MPI processes (e.g. ~1000 atoms,
  ! 100 MPI processes) can result in a significant load imbalance. This appears
  ! to be because worker processes in the task farm decomposition are able to
  ! evaluate tiles more rapidly than the master process can dispatch them.
  ! Consequently, the majority of overall wall time is spent by MPI processes
  ! waiting for tiles to process.

  ! Maximising tile size for the 2Dn-1Da scheme vastly reduces the load
  ! imbalance, by reducing the number of requests from workers to the
  ! master process.

  ! Strings containing human-readable integral evaluation scheme names
  ! Note that an array constructor with variable name character strings
  ! must contain a type specification, i.e.
  ! SCHEME_UNDEFINED = -1
  ! SCHEME_3DC = 0
  ! SCHEME_2DN_1DA = 1
  character(len=32), parameter :: NAME_SCHEME(SCHEME_UNDEFINED:SCHEME_MAX) = [ &
        character(len=32) :: &
        "Undefined", &
        "3D Chebyshev", &
        "2D numerical-1D analytic" ]

  character(len=32), parameter :: hash_table_info_rootname = &
       'hash_table_swri'

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------------!
  ! ****             H I G H - L E V E L    R O U T I N E S              **** !
  !---------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_init_module_stage_0(nbuf, lbuf, qbuf, mbuf, cibuf, cobuf, &
       rwbuf, n_entries)

    !==========================================================================!
    ! Performs stage 0 of the initialisation of the swri module.               !
    !--------------------------------------------------------------------------!
    ! This is called from rundat_blocks and takes care of the lowest-level     !
    ! set-up of the module, it translates the raw data passed from             !
    ! rundat_blocks and stores them in entries of swri_rundat_block_data.      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   See internal_swri() in rundat_blocks.                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 18/02/2015.                                 !
    ! Modified by James C. Womack to add support for flag controlling          !
    ! metric matrix integation scheme, 2018.                                   !
    !==========================================================================!

    use constants, only: SW_V, SW_O
    use rundat, only: pub_ions_will_move
    use utils, only: utils_assert, utils_int_to_str

    implicit none

    ! Parameters
    integer, parameter :: nflags = 10

    ! jd: Arguments
    character(len=32), intent(in) :: nbuf(max_rundat_swri_block_entries)
    integer, intent(in)           :: lbuf(max_rundat_swri_block_entries)
    integer, intent(in)           :: qbuf(max_rundat_swri_block_entries)
    character(len=2), intent(in)  :: mbuf(max_rundat_swri_block_entries)
    integer, intent(in)           :: cibuf(max_rundat_swri_block_entries)
    integer, intent(in)           :: cobuf(max_rundat_swri_block_entries)
    character(len=nflags), intent(in) :: rwbuf(max_rundat_swri_block_entries)
    integer, intent(in)           :: n_entries

    ! jd: Local variables
    integer :: e, c

    ! -------------------------------------------------------------------------

    swri_library_size = n_entries
    do e = 1, n_entries
       swri_rundat_block_data(e)%read_from_file = .false.
       swri_rundat_block_data(e)%write_to_file = .false.
       swri_rundat_block_data(e)%print = .false.
       swri_rundat_block_data(e)%disassemble = .false.
       swri_rundat_block_data(e)%erase = .false.
       swri_rundat_block_data(e)%quit = .false.
       swri_rundat_block_data(e)%has_metric(:) = .false.
       swri_rundat_block_data(e)%entry_name = nbuf(e)
       swri_rundat_block_data(e)%quality%min_l = 0
       swri_rundat_block_data(e)%quality%max_l = lbuf(e)
       swri_rundat_block_data(e)%quality%max_q = qbuf(e)
       swri_rundat_block_data(e)%quality%num_sws_per_centre = -1 ! these will be
       swri_rundat_block_data(e)%quality%max_sws_per_centre = -1 ! filled later
                                                                 ! in stage 1
       swri_rundat_block_data(e)%cheb_intervals = cibuf(e)
       swri_rundat_block_data(e)%cheb_order = cobuf(e)
       swri_rundat_block_data(e)%parent_name= '(@currently unsupported@)       '

       swri_rundat_block_data(e)%int_scheme = SCHEME_UNDEFINED

       do c = 1, 2
          if(mbuf(e)(c:c) == 'V' .or. mbuf(e)(c:c) == 'v') then
             swri_rundat_block_data(e)%has_metric(SW_V) = .true.
          end if
          if(mbuf(e)(c:c) == 'O' .or. mbuf(e)(c:c) == 'o') then
             swri_rundat_block_data(e)%has_metric(SW_O) = .true.
          end if
       end do

       call utils_assert(any(swri_rundat_block_data(e)%has_metric), &
            'block %swri entry #'//trim(utils_int_to_str(e))//' (['//&
            trim(swri_rundat_block_data(e)%entry_name)//']) does not define &
            &any metric. The fourth column must be one of {V,O,VO,OV}')

       do c = 1, nflags
          if(rwbuf(e)(c:c) == 'R') then
             swri_rundat_block_data(e)%read_from_file = .true.
             call utils_assert(.not. pub_ions_will_move, &
                  'block %swri entry #'//trim(utils_int_to_str(e))//' (['//&
                  trim(swri_rundat_block_data(e)%entry_name)//']) asks for the&
                  & metric matrix/matrices to be read from file, but TASK&
                  & specifies the ions will move (which requires re-calculating&
                  & the metric matrix/matrices). I will not let you do that.')
          end if
          if(rwbuf(e)(c:c) == 'W') then
             swri_rundat_block_data(e)%write_to_file = .true.
          end if
          if(rwbuf(e)(c:c) == 'P') then
             swri_rundat_block_data(e)%print = .true.
          end if
          if(rwbuf(e)(c:c) == 'D') then
             swri_rundat_block_data(e)%disassemble = .true.
          end if
          if(rwbuf(e)(c:c) == 'E') then
             swri_rundat_block_data(e)%erase = .true.
          end if
          if(rwbuf(e)(c:c) == 'Q') then
             swri_rundat_block_data(e)%quit = .true.
          end if
          if(rwbuf(e)(c:c) == '2') then
             swri_rundat_block_data(e)%int_scheme = SCHEME_2DN_1DA
          end if
          if(rwbuf(e)(c:c) == '3') then
             swri_rundat_block_data(e)%int_scheme = SCHEME_3DC
          end if
       end do
    end do

  end subroutine swri_init_module_stage_0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_init_module_stage_1(elements, cell, uni_tightbox, par)
    !==========================================================================!
    ! Performs stage 1 of the initialisation of the swri module.               !
    !--------------------------------------------------------------------------!
    ! Reads SW_RI entries from swri_rundat_block_data and performs stage 1     !
    ! initialisation for all SW_RI objects, putting them in swri_library.      !
    ! In stage 1:                                                              !
    ! - Bessels are initialised,                                               !
    ! - Persistent caches (currently none) are initialised,                    !
    ! - Metric matrices are initialised (blocking schemes, parallel strategy), !
    !   but *not* filled,                                                      !
    ! - If spherical harmonic rotations are required to to evaluate the metric !
    !   matrix, a rundat variable is set accordingly.                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   (input): Self explanatory.                                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 18/02/2015.                                 !
    ! Updated by James C. Womack to detect need for sph harm rotations.        !
    !==========================================================================!

    use bibliography, only: bibliography_cite
    use fft_box, only: FFTBOX_INFO
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: rundat_set_pub_use_sph_harm_rot_true
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_flushed_string_output

    implicit none

    ! jd: Arguments
    type(PARAL_INFO), intent(inout) :: par
    type(ELEMENT), intent(in)       :: elements(par%nat)
    type(CELL_INFO), intent(in)     :: cell
    type(FFTBOX_INFO), intent(in)   :: uni_tightbox

    ! jd: Local variables
    integer :: swri_h
    integer :: ierr
    character(len=*), parameter :: myself = 'swri_init_module_stage_1'

    ! -------------------------------------------------------------------------

    call bibliography_cite('HFX')

    call utils_flushed_string_output(&
         CRLF//'SWRI: Initialising module (stage 1)... '//CRLF)

    allocate(swri_library(swri_library_size),stat=ierr)
    call utils_alloc_check(myself,'swri_library',ierr)

    do swri_h = 1, swri_library_size
       call swri_init_container_stage_1(swri_library(swri_h), &   ! output
            swri_rundat_block_data(swri_h)%entry_name, &          ! input
            swri_rundat_block_data(swri_h)%quality, &             ! input
            swri_rundat_block_data(swri_h)%has_metric, &          ! input
            swri_rundat_block_data(swri_h)%int_scheme, &          ! input
            swri_rundat_block_data(swri_h)%cheb_intervals, &      ! input
            swri_rundat_block_data(swri_h)%cheb_order, &          ! input
            elements, uni_tightbox, cell, par)                    ! input
    end do

    ! JCW: Update rundat variable pub_use_sph_harm_rot to .true. if
    ! JCW: spherical harmonic rotations are required for metric matrix
    ! JCW: evaluation (e.g. in 2Dn-1Da scheme). This will ensure that
    ! JCW: the sph_harm_rotation module is initialized.
    if (any(swri_library(:)%int_scheme == SCHEME_2DN_1DA)) then
       call rundat_set_pub_use_sph_harm_rot_true
    end if

    call utils_flushed_string_output(&
         CRLF//'SWRI: Stage 1 completed.'//CRLF)

  end subroutine swri_init_module_stage_1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_init_module_stage_2(elements, overlap, cell)
    !==========================================================================!
    ! Performs stage 2 of the initialisation of the swri module.               !
    !--------------------------------------------------------------------------!
    ! Essentially calls stage 2 initialisation for all SW_RIs in the library.  !
    ! In stage 2 the metric matrices are filled (by Chebyshev interpolation or !
    ! by reading from files).                                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 18/02/2015.                                 !
    !==========================================================================!

    use ion, only: ELEMENT
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3
    use utils, only: utils_flushed_string_output

    implicit none

    ! jd: Arguments
    type(ELEMENT), intent(in)   :: elements(:)
    type(SPAM3), intent(in)     :: overlap
    type(CELL_INFO), intent(in) :: cell

    ! jd: Local variables
    integer :: swri_h

    ! -------------------------------------------------------------------------

    call utils_flushed_string_output(&
         CRLF//'SWRI: Initializing module (stage 2)... '//CRLF)

    do swri_h = 1, swri_library_size
       call swri_init_container_stage_2(swri_library(swri_h), & ! output
            elements, overlap, cell, &                          ! input
            swri_rundat_block_data(swri_h)%read_from_file, &    ! input
            swri_rundat_block_data(swri_h)%write_to_file, &     ! input
            swri_rundat_block_data(swri_h)%print, &             ! input
            swri_rundat_block_data(swri_h)%disassemble, &       ! input
            swri_rundat_block_data(swri_h)%erase, &             ! input
            swri_rundat_block_data(swri_h)%quit)                ! input
    end do

    call utils_flushed_string_output('SWRI: Done.'//CRLF)

  end subroutine swri_init_module_stage_2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_cleanup_module
    !==========================================================================!
    ! This subroutine cleans up after the swri module.                         !
    !--------------------------------------------------------------------------!
    ! Calls destroy_container() on every SW_RI object in the library, then     !
    ! deallocates the library itself.                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 18/02/2015.                                 !
    !==========================================================================!

    use utils, only: utils_dealloc_check, utils_flushed_string_output

    implicit none

    ! jd: Local variables
    integer :: swri_h
    integer :: ierr
    character(len=*), parameter :: myself = 'swri_cleanup_module'

    ! -------------------------------------------------------------------------

    call utils_flushed_string_output('SWRI: Clean-up... ')

    do swri_h = 1, swri_library_size
       call swri_destroy_container(swri_library(swri_h))
    end do

    deallocate(swri_library,stat=ierr)
    call utils_dealloc_check(myself,'swri_library',ierr, &
         allocated_in = 'swri_init_module_stage_1')
    call utils_flushed_string_output('done.'//CRLF)

  end subroutine swri_cleanup_module

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_init_container_stage_1(swri, &             ! output
     swri_name, quality, has_metric, int_scheme, &           ! input
     cheb_intervals, cheb_order, &                           ! input
     elements, uni_tightbox, cell, par)                      ! external input
    !==========================================================================!
    ! This subroutine creates an SW_RI container.                              !
    ! - Bessels are initialised,                                               !
    ! - Persistent caches (currently none) are initialised,                    !
    ! - Metric matrices are initialised (blocking schemes, parallel strategy), !
    !   but *not* filled.                                                      !
    ! - swri's quality's max_sws_per_centre and num_sws_per_centre are set.    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   - swri (out): Contains the intialised SW_RI container on ouput.        !
    !   - swri_name, quality, has_metric, int_scheme, cheb_intervals,          !
    !     cheb_order (all in): describe the container.                         !
    !   - elements, uni_tightbox, par (all in): self-explanatory.              !
    !   - cell (in): cell%n_pts is needed for dimensioning PPD data.           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 17/04/2012.                                 !
    ! Modified by Jacek Dziedzic on 08/04/2014.                                !
    ! Reworked by Jacek Dziedzic on 12/02/2015.                                !
    ! Modified by James C. Womack to add support for multiple off-site metric  !
    ! matrix element integation schemes, 2018.                                 !
    !==========================================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: CRLF, SW_V, SW_O, SW_LETTERS, garbage_int
    use fft_box, only: FFTBOX_INFO
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_swri_cheb_batchsize, pub_swri_verbose, pub_use_hfx, &
         pub_use_activehfx
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort, utils_assert, utils_alloc_check, &
         utils_flushed_string_output, utils_int_to_str

    implicit none

    ! jd: Arguments
    type(SW_RI),       intent(out)   :: swri
    character(len=*),  intent(in)    :: swri_name
    type(SW_QUALITY),  intent(in)    :: quality
    logical,           intent(in)    :: has_metric(2)
    integer,           intent(in)    :: int_scheme
    integer,           intent(in)    :: cheb_intervals
    integer,           intent(in)    :: cheb_order
    type(PARAL_INFO),  intent(inout) :: par
    type(ELEMENT),     intent(in)    :: elements(par%nat)
    type(FFTBOX_INFO), intent(in)    :: uni_tightbox
    type(CELL_INFO),   intent(in)    :: cell

    ! jd: Local variables
    integer       :: ierr
    real(kind=DP) :: rad
    integer       :: species
    integer       :: q_idx
    integer       :: sw_idx
    integer       :: bessel_idx
    integer       :: n_bessels, n_sws
    integer       :: old_l, l
    integer       :: m
    integer       :: V_or_O
    integer       :: swri_cheb_batchsize_loc

    character(len=*), parameter :: myself = 'swri_init_container_stage_1'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)

    call utils_flushed_string_output(&
         'SWRI: ['//trim(swri_name)//'] Initialising (stage 1)...'//CRLF)

    swri%swri_name = swri_name
    swri%quality = quality
    swri%has_metric(:) = has_metric(:)
    swri%cheb_intervals = cheb_intervals
    swri%cheb_order = cheb_order
    swri%natoms = par%nat

    allocate(swri%image_map(cell%n_pts,cell%n_ppds,par%nat),stat=ierr)
    call utils_alloc_check(myself,'image_map',ierr)
    swri%image_map(:,:,:) = garbage_int

    ! JCW: Set integration scheme
    ! JCW: If explicitly set in swri rundat block data, then respect this choice
    ! JCW: If not set in swri rundat block data, then select default
    if (int_scheme == SCHEME_UNDEFINED) then
       ! JCW: Set default
       ! @2DN_1DA_OVERLAP
       ! JCW: At present, only a single integration scheme can be selected for
       ! JCW: all metric matrices to be evaluated. Since not all schemes support
       ! JCW: all metric matrix types, some logic is required to select the
       ! JCW: appropriate default based on the metric matrices requested
       select case (count(swri%has_metric))
       case (1)
         ! JCW: If only one of SW_V and SW_O are selected, then simply use
         ! JCW: default for that metric matrix type
         do V_or_O = SW_V, SW_O
            if (swri%has_metric(V_or_O)) then
               swri%int_scheme = SCHEME_DEFAULT(V_or_O)
               call utils_flushed_string_output(&
                    'SWRI: ['//trim(swri_name)//'] - Using default '//&
                    trim(SW_LETTERS(V_or_O))//' metric matrix evaluation &
                    &scheme ('//trim(NAME_SCHEME(SCHEME_DEFAULT(V_or_O)))//&
                    ').'//CRLF)
               exit
            end if
         end do
       case (2)
         ! JCW: If both SW_V and SW_O are selected, then use the default
         ! JCW: for SW_O, since SW_O only supports the 3Dc scheme
         swri%int_scheme = SCHEME_DEFAULT(SW_O)
         call utils_flushed_string_output(&
              'SWRI: ['//trim(swri_name)//'] - Using default '//&
              trim(SW_LETTERS(SW_O))//' metric matrix evaluation &
              &scheme ('//trim(NAME_SCHEME(SCHEME_DEFAULT(SW_O)))//&
              ').'//CRLF)
         ! JCW ... and emit a warning that informs the user that a slower
         ! JCW: default scheme has been selected for evaluating SW_V
         call utils_flushed_string_output(&
              'SWRI: [WARNING] Evaluating the '//trim(SW_LETTERS(SW_O))//&
              ' matrix is currently only possible using the '//&
              trim(NAME_SCHEME(SCHEME_3DC))//' scheme.'//CRLF)
         call utils_flushed_string_output(&
              'SWRI: [WARNING] Only one metric matrix evaluation scheme can &
              &be used per SWRI specification, so the '//&
              trim(SW_LETTERS(SW_V))//' matrix will also be evaluated using &
              &this scheme.'//CRLF)
         call utils_flushed_string_output(&
             'SWRI: [WARNING] The '//trim(NAME_SCHEME(SCHEME_2DN_1DA))&
             //' scheme is significantly less computationally intensive for &
             &evaluating the '//trim(SW_LETTERS(SW_V))//' matrix.'//CRLF)
         call utils_flushed_string_output(&
             'SWRI: [WARNING] Consider switching to evaluating the '//&
             trim(SW_LETTERS(SW_V))//' matrix only, if possible.'//CRLF)
       case default
         ! JCW: Something has gone wrong if we end up here...
         call utils_abort('Error in '//myself//': Unexpected number of metric &
              &matrix types in SWRI specification. Expect 1 or 2, but got',&
              opt_int_to_print1=count(swri%has_metric))
       end select
    else
       ! JCW: Set user-selected scheme
       ! JCW: First check that the selected scheme is in allowed range
       call utils_assert(int_scheme>=SCHEME_MIN.and.int_scheme<=SCHEME_MAX,&
            'Error in '//myself//': Unexpected integration scheme &
            &specification, ',int_scheme)
       ! JCW: If in allowed range, state which scheme is being used.
       call utils_flushed_string_output(&
         'SWRI: ['//trim(swri_name)//'] - Using user-selected metric matrix &
            &evaluation scheme ('//trim(NAME_SCHEME(int_scheme))//').'//CRLF)
       swri%int_scheme = int_scheme
    end if

    call utils_flushed_string_output('SWRI: ['//trim(swri_name)//&
         '] - Initialising Bessels (l_max='//&
         trim(utils_int_to_str(quality%max_l))//')... ')

    ! jd: Initialise other internal 'constants'
    swri%quality%max_sws_per_centre = quality%max_q * &
         ((1+quality%max_l-quality%min_l) * (1+quality%max_l+quality%min_l))

    ! @jd: This is not entirely right -- we shouldn't be going over all species
    !      but only over those in the SWRI.
    ! qoh: Allocate arrays for Bessel and Spherical Wave data
    allocate(swri%radtable(par%num_species),stat=ierr)
    call utils_alloc_check(myself,'swri%radtable',ierr)

    ! qoh: Initialise radtable, max_num_bessels
    ! JCW: Also initialize number of q-values per l
    allocate(swri%num_q_per_l(par%num_species,quality%min_l:quality%max_l),&
         stat=ierr)
    call utils_alloc_check(myself,'swri%num_q_per_l',ierr)
    call swri_num_sph_functions(swri%radtable, swri%max_num_bessels, & ! fills
         swri%quality, par, num_q_per_l=swri%num_q_per_l)
    allocate(swri%sphbessels(par%num_species,swri%max_num_bessels),stat=ierr)
    call utils_alloc_check(myself,'swri%sphbessels',ierr)
    call swri_sph_bessels_init(swri,uni_tightbox, par%num_species)

    ! jd: Now that we have radtable, we can ensure the cell is big enough
    !     if PBCs are in effect
    call internal_check_cell_size(swri,cell)

    if (pub_swri_verbose.and.pub_on_root) then
       ! JCW: Tell us the number of allowed q-values for each species
       write(stdout,"(a)") "SWRI: Number of allowed q-values per species"
       write(stdout,"(a)") "SWRI: (One line per species, values listed per-l)"
       write(stdout,"(a5,1x,a7,1x,a12)") "SWRI:","Species","Num q-values"
       write(stdout,"(a5,1x,a7,1x,a12)") "SWRI:","-------","------------"
       do species = 1, par%num_species
          write(stdout,"(a5,1x,i7,1x)",advance="no") "SWRI:", species
          do l = quality%min_l, quality%max_l
             write(stdout,"(i0,1x)",advance="no") swri%num_q_per_l(species,l)
          end do
          write(stdout,*) ! End record
       end do
       write(stdout,"(a,i0)") "SWRI: User requested qmax: ", quality%max_q
    end if

    ! Check for unsupported states when using 2Dn-1Da scheme
    if (swri%int_scheme == SCHEME_2DN_1DA) then
       ! @2DN_1DA_OVERLAP
       ! Support for overlap metrix not (yet) implemented using 2Dn-1Da scheme
       call utils_assert(.not.swri%has_metric(SW_O),"Error in "//myself//": &
            &Current implementation of 2Dn-1Da integration scheme does not &
            &support the overlap metric, sorry.")

       ! Check some assumptions used by the sph_harm_rotation module when
       ! applying RSH rotation matrices to atomblocks
       ! * swri%quality%min_l == 0
       call utils_assert(swri%quality%min_l == 0,"Error in "//myself//": &
            &The 2Dn-1Da integration scheme requires that that min_l == 0 in &
            &the SWRI.")

       ! @IDENTICAL_RADII
       ! * swri%num_q_per_l(species,0:swri%quality%max_l) is identical for all
       !   species
       do species = 1, par%num_species
          call utils_assert(&
               all(swri%num_q_per_l(species,0:swri%quality%max_l)&
               ==swri%num_q_per_l(1,0:swri%quality%max_l)),&
               "Error in "//myself//": &
               &The 2Dn-1Da integration scheme requires that the number of &
               &q-values per l-value is identical for all species.")
       end do
    end if

    ! @IDENTICAL_RADII
    call utils_assert(minval(swri%radtable) == maxval(swri%radtable), &
         "Current implementation of spherical-wave resolution of identity &
         &only works if all NGWF radii are identical, sorry.")
    rad = swri%radtable(1)
    species = 1 ! @IDENTICAL_RADII

    ! jd: Count bessels, spherical waves
    n_sws = 0
    n_bessels = 0
    do bessel_idx = 1, swri%max_num_bessels
       l = swri%sphbessels(species,bessel_idx)%lval
       if(l == -1) exit
       n_bessels = n_bessels + 1
       n_sws = n_sws + (2*l+1) ! jd: account for m = -l..l
    end do
    call utils_assert(n_bessels > 0,'n_bessels must be positive')
    call utils_assert(n_sws > 0,'n_sws must be positive')

    swri%quality%num_sws_per_centre = n_sws

    if (pub_swri_cheb_batchsize > 0) then
       ! Use user-specified batch size (keyword has been set)
       swri_cheb_batchsize_loc = pub_swri_cheb_batchsize
    else
       ! Use per-integration scheme default batch size (keyword not set)
       swri_cheb_batchsize_loc = DEFAULT_XBATCHSIZE(swri%int_scheme)
    end if

    swri%xbatchsize = min(n_sws, swri_cheb_batchsize_loc)

    ! jd: Fill sw_idx_to_bessel_idx and sw_idx_to_m
    allocate(swri%sw_idx_to_bessel_idx(swri%quality%max_sws_per_centre),stat=ierr)
    call utils_alloc_check(myself,'swri%sw_idx_to_bessel_idx',ierr)
    allocate(swri%sw_idx_to_m(swri%quality%max_sws_per_centre),stat=ierr)
    call utils_alloc_check(myself,'swri%sw_idx_to_m',ierr)
    allocate(swri%sw_idx_to_l(swri%quality%max_sws_per_centre),stat=ierr)
    call utils_alloc_check(myself,'swri%sw_idx_to_l',ierr)
    allocate(swri%sw_idx_to_q_idx(swri%quality%max_sws_per_centre),stat=ierr)
    call utils_alloc_check(myself,'swri%sw_idx_to_q_idx',ierr)
    swri%sw_idx_to_bessel_idx(:) = -1
    swri%sw_idx_to_m(:) = -1
    swri%sw_idx_to_l(:) = -1
    swri%sw_idx_to_q_idx(:) = -1
    sw_idx = 1
    q_idx = 1
    l = -1
    do bessel_idx = 1, n_bessels
       old_l = l
       l = swri%sphbessels(species,bessel_idx)%lval
       if(old_l /= l) q_idx = 1
       do m = -l, l
          swri%sw_idx_to_bessel_idx(sw_idx) = bessel_idx
          swri%sw_idx_to_m(sw_idx) = m
          swri%sw_idx_to_l(sw_idx) = l
          swri%sw_idx_to_q_idx(sw_idx) = q_idx
          sw_idx = sw_idx + 1
       end do
       q_idx = q_idx + 1
    end do

    ! jd: Prepare interpolation of SW radial components. Only used in HFx.
    if(pub_use_hfx .or. pub_use_activehfx) then
       call utils_flushed_string_output('done.'//CRLF//&
            'SWRI: ['//trim(swri_name)//'] - Initialising radial lookup... ')
       call swri_init_sw_radial_lookup(swri, cell)
    end if

    ! jd: Calculate number of xbatches from the batch size
    swri%n_xbatches = swri%quality%num_sws_per_centre / swri%xbatchsize
    if(mod(swri%quality%num_sws_per_centre,swri%xbatchsize) /= 0) then
       swri%n_xbatches  = swri%n_xbatches + 1
    end if

    call comms_barrier ! jd: Otherwise we wait on the flush below, giving the
                       !     impression that the cache init takes a long time

    call utils_flushed_string_output('done.'//CRLF//&
         'SWRI: ['//trim(swri_name)//'] - SW batch size (tile width) &
         &set to '//trim(utils_int_to_str(swri%xbatchsize))//'.'//CRLF)

    call utils_flushed_string_output(&
         'SWRI: ['//trim(swri_name)//&
         '] - Initialising persistent cache structures... ')

    ! jd: Initialise persistent caches
    call swri_init_persistent_caches(swri, cell%n_pts)

    call utils_flushed_string_output('done.'//CRLF//&
         'SWRI: ['//trim(swri_name)//'] - Initialising metric matrix... ')

    ! jd: Initialise (but not fill) the metric matrix
    call swri_init_metric_matrix(swri,elements,par)

    call utils_flushed_string_output('done.'//CRLF//&
         'SWRI: ['//trim(swri_name)//'] Done.')

    swri%created = .true.
    swri%full_metric_matrices_destroyed = .false.

    call utils_trace_out(myself)

    return

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_check_cell_size(swri,cell)
      !========================================================================!
      ! Ensures that the cell is suitably big for SWRI in PBCs, if PBCs are    !
      ! in effect. The condition for this is that the cell dimension along any !
      ! axis is at least 4 rNGWF. In OBCs nothing is checked -- the usual      !
      ! requirement that the cell dimensions along any axis is at lest 2 rNGWF !
      ! in that case will be trapped by parallel strategy.                     !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2022.09.                                  !
      !========================================================================!

      use geometry, only: magnitude
      use rundat, only: pub_hfx_bc_is_periodic
      use simulation_cell, only: CELL_INFO
      use utils, only: utils_abort, utils_real_to_str

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(in)     :: swri
      type(CELL_INFO), intent(in) :: cell

      ! jd: Local variables
      real(kind=DP) :: rNGWF
      character(len=*), parameter :: myself = 'swri::internal_check_cell_size'

      ! -----------------------------------------------------------------------

      if(.not. any(pub_hfx_bc_is_periodic)) then
         ! jd: OBC
         return
      else if(all(pub_hfx_bc_is_periodic)) then
         ! jd: PBC
         rNGWF = swri%radtable(1) ! @IDENTICAL_RADII
         if(magnitude(cell%a1) < 4.0_DP * rNGWF .or. &
              magnitude(cell%a2) < 4.0_DP * rNGWF .or. &
              magnitude(cell%a3) < 4.0_DP * rNGWF) then
            call utils_abort('Spherical-wave resolution of identity &
                 &in PBCs (hfx_pbc P P P) requires'//CRLF//'all cell dimensions&
                 & to be at least 4 times the NGWF radius (in contrast to'//&
                 CRLF//'2 times the (largest) NGWF radius for standard ONETEP&
                 & operation.'//CRLF//'If you intend to use SWRI (needed for &
                 &HFx, hybrid functionals, DMA, polarisable'//CRLF//&
                 'embedding), you will either have to make the cell larger, &
                 &the NGWFs smaller,'//CRLF//'or choose open boundary &
                 &conditions (hfx_pbc O O O).'//CRLF//'The NGWF radius was '//&
                 trim(utils_real_to_str(rNGWF,'f12.3'))//' a0, making the &
                 &shortest possible cell vector '//&
                 trim(utils_real_to_str(4.0_DP * rNGWF,'f12.3'))//' a0.'//&
                 CRLF//'Length of cell vector a1 was '//&
                 trim(utils_real_to_str(magnitude(cell%a1),'f12.3'))//' a0.'//&
                 CRLF//'Length of cell vector a2 was '//&
                 trim(utils_real_to_str(magnitude(cell%a2),'f12.3'))//' a0.'//&
                 CRLF//'Length of cell vector a3 was '//&
                 trim(utils_real_to_str(magnitude(cell%a3),'f12.3'))//' a0.')
         end if
      else
         call utils_abort(myself//': Mixed BCs are not currently supported.')
      end if

    end subroutine internal_check_cell_size

  end subroutine swri_init_container_stage_1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_init_container_stage_2(swri, &                       ! in/out
       elements, overlap, cell, read_from_file, write_to_file, &       ! input
       print, disassemble, erase, quit)                                ! input
    !==========================================================================!
    ! This subroutine completes the intialisation of an SW_RI container:       !
    ! - It sparse-creates the metric matrices.                                 !
    ! - It populates the S-neighbour list of the SW_RI.                        !
    ! - If fills the metric matrices of an SW_RI container.                    !
    ! The remainder of this description deals with the filling.                !
    !                                                                          !
    ! If 'read_from_file' is T, the entire matrices are read from binary files.!
    ! Otherwise:                                                               !
    !  - On-site elements, for which analytical formulas are available, are    !
    !    calculated.                                                           !
    !  - A list of off-site atomblocks is constructed based on the atoms that  !
    !    belong to the SWRI ('members'). We need blocks between members *and*  !
    !    between members and S-neighbours of members (see 'subsystem case'     !
    !    below for more detail). Due to symmetry considerations, only one      !
    !    triangle of the matrices is needed. Because load is balanced dynami-  !
    !    cally, the checkerboard pattern is not used anymore.                  !
    !  - The list of off-site atomblocks is pruned to remove atomblocks, which !
    !    can be succesfully read from files. Atomblock files have the advantage!
    !    of being transferable between two calculations which share some of the!
    !    atomic positions. The result is a list of atomblocks that need to be  !
    !    computed.                                                             !
    !  - Off-site atomblocks are computed through Chebyshev interpolation,     !
    !    with functions expanded in 3-D or 2-D, depending on the selected      !
    !    metric matrix integration scheme. For 2-D Chebyshev interpolation     !
    !    the atomblock is evaluated in the "integral coordinate frame" where   !
    !    the z-axis is aligned along the interatomic coordinate. This allows   !
    !    the integral to be evaluated as a product of a 2-D numerical and 1-D  !
    !    analytic integral (2Dn-1Da scheme).                                   !
    !  - Each block is divided into tiles (slabs), whose thickness depends on  !
    !    swri_cheb_batchsize. A task farm decomposition is used, with the root !
    !    proc dispatching tiles to workers upon request. The root proc is also !
    !    a worker, it only probes for pending comms more often to avoid convoy !
    !    effects. Once all tiles of an atomblock are computed, the atomblock is!
    !    written out to a file, and is inserted into the SPAM3 metric matrix.  !
    !    The atomblock files are not erased by default.                        !
    !  - Once all atomblocks are completed, the metric matrix(/matrices)       !
    !    is/are ready.                                                         !
    !  - For 2-D Chebyshev expansion (2Dn-1Da scheme), a real-spherical        !
    !    harmonic (RSH) rotation matrix is computed and applied to the         !
    !    atomblock to evaluate the atomblock in terms of SWs in ONETEP's       !
    !    Cartesian coordinate system.                                          !
    !    If 'write_to_file' is T, the metric matrix(/matrices) is/are written  !
    !    to a file.                                                            !
    !  - If 'print' is T, the metric matrix(/matrices) is/are printed out in   !
    !    text format to stdout. This can be useful for validation.             !
    !  - If 'disassemble' is T, the metric matrix(/matrices) is/are disassem-  !
    !    bled into atomblocks, and each of these is written to a file. This is !
    !    useful for obtaining the atomblock representation when having only    !
    !    the entire-matrix representation.                                     !
    !  - If 'erase' is T, an attempt is made to erase the atomblock files.     !
    !    No error checking is performed.                                       !
    !  - If 'quit' is T, ONETEP quits gracefully at this point. This is useful !
    !    if only doing disassembles, for instance.                             !
    !                                                                          !
    ! As completed tiles are accumulated, they are temporarily kept in a hash  !
    ! table on the proc that is the sparse-owner of the atomblock's atom-a.    !
    ! Once all tiles of a block are ready, they are merged into an atomblock,  !
    ! and the hash table entries are evicted. The memory overhead of keeping   !
    ! the hash table is bounded and very small. The RSH rotation matrix is     !
    ! applied in the 2Dn-1Da integration scheme immediately after the          !
    ! (integral coordinate frame) atomblock has been merged.                   !
    !                                                                          !
    ! OpenMP parallelism is contained within two compute-intensive subroutines:!
    ! - the Chebyshev interpolation internal_tile_by_chebs), and               !
    ! - the calculation of the SWs in the home sphere (internal_home_sws).     !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swri (in/out): The container to initialise, assumed to have undergone  !
    !                  stage-1 initialisation already.                         !
    !   elements (in): Needed for ionic positions.                             !
    !   overlap (in):  Needed only to establish S-neighbours for the ket side, !
    !                  because we only need V/O matrix elements for NGWFs that !
    !                  overlap.                                                !
    !   cell (in):     Simulation cell. Needed for PBCs.                       !
    !   read_from_file(in): Whether to read the metric matrices from files.    !
    !   write_to_file(in): Whether to write the metric matrices to files.      !
    !   print(in): Whether to print the metric matrices to stdout.             !
    !   disassemble(in): Whether to disassemble the metric matrices to atomblk !
    !                    files.                                                !
    !   erase(in): Whether to erase atomblock files.                           !
    !   quit(in): Whether to quit after having performed the operations above. !
    !--------------------------------------------------------------------------!
    ! Subsystem case:                                                          !
    !   When the SW_RI is only defined for a subsystem:                        !
    !   - on-site blocks are needed for atoms in the subsystem *and* for their !
    !     neighbours.                                                          !
    !   - off-site blocks are needed between every pair at least one atom of   !
    !     which belongs to the subsystem.                                      !
    !--------------------------------------------------------------------------!
    ! Timing metric matrix evaluation:                                         !
    !   If the execution time to evaluate the metric matrix is desired, then   !
    !   then devel_code SW:QUIT_WITH_TIMER=T:SW can be used, which will ensure !
    !   that timing data is output when the calculation is aborted immediately !
    !   following evaluation of the metric matrix (i.e. "Q" flag is used in    !
    !   the SWRI block).                                                       !
    !--------------------------------------------------------------------------!
    ! - Written by Quintin Hill in late 2009 and early 2010.                   !
    ! - Severely modified by Jacek Dziedzic in 2012.                           !
    ! - Parallelisation reworked by Jacek Dziedzic in November 2013.           !
    ! - Adapted to be used in a scenario where V is not needed (DMA w/o HFx)   !
    !   by Jacek Dziedzic in April 2014.                                       !
    ! - Generalised for SW_RI by Jacek Dziedzic in February 2015.              !
    ! - Generalised to scenarios where the SW_RI does not encompass the entire !
    !   system by Jacek Dziedzic in April 2015.                                !
    ! - Generalised to add support for multiple off-site metric matrix element !
    !   evaluation schemes by James Womack, 2018.                              !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: CRLF, SW_O, SW_V
    use ion, only: ELEMENT
    use neighbour_list, only: neighbour_list_init_from_sparse
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_use_activeswx, pub_ngwf_regions_ngroups, pub_devel_code
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_create, sparse_show_matrix, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_dealloc_check, utils_devel_code, &
         utils_flushed_string_output

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout)    :: swri
    type(ELEMENT), intent(in)     :: elements(:)
    type(SPAM3), intent(in)       :: overlap
    type(CELL_INFO), intent(in)   :: cell
    logical, intent(in)           :: read_from_file
    logical, intent(in)           :: write_to_file
    logical, intent(in)           :: print
    logical, intent(in)           :: disassemble
    logical, intent(in)           :: erase
    logical, intent(in)           :: quit

    ! jd: Local Variables
    integer :: swri_handle
    integer :: ierr
    logical :: quit_with_timer
    type(ATOMBLOCK_DESCRIPTOR), allocatable :: atomblocklist(:)
    character(len=*), parameter :: myself = 'swri_init_container_stage_2'
    type(PARAL_INFO), pointer   :: par
    character(len=30) :: matstruc

    ! -----------------------------------------------------------------------

    if (swri%initialised) return

    call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
         '] Initialising (stage 2)... '//CRLF)

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    swri_handle = swri_get_handle_to(swri%swri_name)
    swri%store_blocks_on_disk = .not. erase

    ! jd: Sparse-create metric matrices
    if(swri%has_metric(SW_V)) then
       ! rc2013: if we're only using SWRI in active region set indices correctly
       if(pub_use_activeswx .and. pub_ngwf_regions_ngroups .gt. 1) then
          matstruc   = 'Vs11'
       else
          matstruc   = 'Vs'
       end if
       swri%full_metric_matrices(SW_V)%structure = matstruc
       call sparse_create(swri%full_metric_matrices(SW_V))
    end if
    if(swri%has_metric(SW_O)) then
       ! rc2013: if we're only using SWRI in active region set indices correctly
       if(pub_use_activeswx .and. pub_ngwf_regions_ngroups .gt. 1) then
          matstruc   = 'Ss11'
       else
          matstruc   = 'Ss'
       end if
       swri%full_metric_matrices(SW_O)%structure = matstruc
       call sparse_create(swri%full_metric_matrices(SW_O))
    end if

    ! rc2013: get the parallel strategy from the diagonal matrix
    call sparse_get_par(par, overlap, 'R')

    ! jd: Initialise the S-neighbour list
    call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
         '] - ')
    call neighbour_list_init_from_sparse(swri%s_atoms_nl, 'SWRI_S_ATOMS', &
         overlap)

    ! jd: Read metric matrices from file?
    if(read_from_file) then
       if(swri%has_metric(SW_V)) then
          call internal_read_matrix(swri%full_metric_matrices(SW_V), &
               swri%swri_name, trim(swri%swri_name)//'.vmatrix')
       end if
       if(swri%has_metric(SW_O)) then
          call internal_read_matrix(swri%full_metric_matrices(SW_O), &
               swri%swri_name, trim(swri%swri_name)//'.omatrix')
       end if
    end if

    ! jd: Generate the list of atomblocks needed, unless individual atomblocks
    !     will never be needed. This allocates and populates atomblocklist.
    !     This also assembles the blocks that are successfully found on disk,
    !     (unless read_from_file is .true.).
    if(.not. read_from_file .or. disassemble .or. erase) then
       call internal_gen_atomblock_list(swri, elements, cell, atomblocklist, &
            .not. read_from_file, erase)
    end if

    ! ---------------------------------------------------------------------
    ! --- The main branch.                                              ---
    ! --- Calculate same-centre blocks analytically.                    ---
    ! --- Look for atomblock files, assemble them into metric matrices. ---
    ! --- Atomblocks not found in files are calculated on demand.       ---
    ! ---------------------------------------------------------------------
    if(.not. read_from_file) then
       call internal_calc_metric_matrices_sc(swri, elements)
       if(size(atomblocklist) == 0) then
          call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
               '] - Not calculating off-site elements, all have been read &
               &from disk.'//CRLF)
          ! jd: When the entire matrix has been read from blocks, be sure to
          !     fill the rest by symmetry, as this usually happens in
          !     internal_calc_metric_matrices_oc(), which is now not called.
          call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
              '] - Filling the remainder by symmetry (path 2).'//CRLF)
          call internal_oc_fill_by_symmetry(swri)
       else
          call internal_calc_metric_matrices_oc(swri, atomblocklist, elements, &
               cell)
       end if
    end if

    ! jd: Write metric matrices to a file, if requested
    if(write_to_file) then
       if(swri%has_metric(SW_V)) then
          call internal_write_matrix(swri%full_metric_matrices(SW_V), &
               swri%swri_name, trim(swri%swri_name)//'.vmatrix')
       end if
       if(swri%has_metric(SW_O)) then
          call internal_write_matrix(swri%full_metric_matrices(SW_O), &
               swri%swri_name, trim(swri%swri_name)//'.omatrix')
       end if
    end if

    ! jd: Print metric matrices, if requested
    if(print) then
       if(swri%has_metric(SW_V)) then
          if(pub_on_root) then
             print *, ""
             print *, "-VMAT:-=============================@"
          end if
          call sparse_show_matrix(swri%full_metric_matrices(SW_V))
          if (pub_on_root) then
             print *, ""
             print *, "-------------------------------------"
          end if
       end if
       if(swri%has_metric(SW_O)) then
          if(pub_on_root) then
             print *, ""
             print *, "-OMAT:-=============================@"
          end if
          call sparse_show_matrix(swri%full_metric_matrices(SW_O))
          if (pub_on_root) then
             print *, ""
             print *, "-------------------------------------"
          end if
       end if
    end if ! needs to be printed?

    ! jd: Disassemble metric matrices to atomblocks, if requested
    if(disassemble) then
       call internal_disassemble(swri,atomblocklist)
    end if ! disassemble

    ! jd: Deallocate atomblocklist
    if(.not. read_from_file .or. disassemble .or. erase) then
       deallocate(atomblocklist,stat=ierr)
       call utils_dealloc_check(myself,'atomblocklist',ierr)
    end if

    if(quit) then
       ! JCW: If devel_code SW:QUIT_WITH_TIMER=T:SW is used, then call
       ! JCW: timer_clock routine to output timing information gathered
       ! JCW: so far before quitting following metric matrix evaluation.
       quit_with_timer = utils_devel_code(.false.,'SW','QUIT_WITH_TIMER',&
            pub_devel_code)
       if (quit_with_timer) then
          ! Finish timing current routine (swri_init_container_stage_2)
          call timer_clock(myself,2)
          ! Output timing data for ONETEP run up to this point
          call timer_clock(myself,3)
       end if

       call utils_abort('Intentional abort on account of ['//&
            trim(swri%swri_name)//']''s Q flag.')
    end if

    swri%initialised = .true.

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_calc_metric_matrices_sc(swri, elements)
      !========================================================================!
      ! Calculates the same-centre (on-site) elements of the metric matrices.  !
      ! These are calculated analytically.                                     !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri(in/out): The SW_RI whose on-site metric matrix elements are to  !
      !                 be filled. Results are returned here.                  !
      !   elements(in): Needed to determine which atoms are within the SW_RI.  !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      ! Modified to accept par as input by Robert Charlton, 13/09/2018.        !
      !========================================================================!

      use comms, only: pub_my_proc_id
      use constants, only: SW_O
      use ion, only: ELEMENT
      use rundat, only: pub_devel_code
      use sparse, only: sparse_put_block
      use utils, only: utils_devel_code, utils_flushed_string_output, &
           utils_alloc_check, utils_dealloc_check
      use timer, only: timer_clock

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(inout) :: swri
      type(ELEMENT), intent(in)  :: elements(par%nat)

      ! jd: Local variables
      integer :: local_a, atoma, orig_atoma ! atom indices
      integer :: b_idx, atomb, orig_atomb   ! atom indices
      type(ATOM_CENTRE) :: centrea, centreb ! Centre data
      integer :: swri_handle
      logical :: has_neighbour_in_swri
      logical :: debug_show_skipped
      real(kind=DP), allocatable :: vatomblock(:,:) ! Atomblock of V matrix
      real(kind=DP), allocatable :: oatomblock(:,:) ! Atomblock of O matrix
      character(len=*), parameter :: myself='internal_calc_metric_matrices_sc'

      ! -----------------------------------------------------------------------

      call timer_clock('swri_metric_mat_sc',1)

      if(swri%has_metric(SW_V) .and. swri%has_metric(SW_O)) then
         call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
              '] - Calculating on-site elements of metric matrices V, O... ')
      end if
      if(swri%has_metric(SW_V) .and. .not. swri%has_metric(SW_O)) then
         call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
              '] - Calculating on-site elements of the metric matrix V... ')
      end if
      if(swri%has_metric(SW_O) .and. .not. swri%has_metric(SW_V)) then
         call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
             '] - Calculating on-site elements of the metric matrix O... ')
      end if

      swri_handle = swri_get_handle_to(swri%swri_name)
      debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
           pub_devel_code, no_bcast=.true., no_warn=.true.)

      ! jd: Allocate workspace atomblock
      if(swri%has_metric(SW_V)) then
         allocate(vatomblock(&
              swri%quality%max_sws_per_centre,swri%quality%max_sws_per_centre),&
              stat=ierr)
         call utils_alloc_check(myself,'vatomblock',ierr)
      end if
      if(swri%has_metric(SW_O)) then
         allocate(oatomblock(&
              swri%quality%max_sws_per_centre,swri%quality%max_sws_per_centre),&
              stat=ierr)
         call utils_alloc_check(myself,'oatomblock',ierr)
      end if

      ! jd: Over all atoms local to this proc
      loop_A_SC: do local_a=1, par%num_atoms_on_proc(pub_my_proc_id)
         atoma = par%first_atom_on_proc(pub_my_proc_id) + local_a -1
         orig_atoma = par%orig_atom(atoma)

         ! jd: Atoms that are neither represented in this SWRI *nor* are their
         !     S-neighbours should be ignored. Remember, we *do* need
         !     same-centre elements for neighbours!
         if(.not. elements(orig_atoma)%in_swri(swri_handle)) then
            ! jd: Atom A is not in SWRI. Check neighbours.
            has_neighbour_in_swri = .false.
            do b_idx = swri%s_atoms_nl%first_idx(atoma), &
                 swri%s_atoms_nl%last_idx(atoma)
               atomb = swri%s_atoms_nl%neighbours(b_idx)
               orig_atomb = par%orig_atom(atomb)
               if(elements(orig_atomb)%in_swri(swri_handle)) then
                  has_neighbour_in_swri = .true.
                  exit
               end if
            end do
            if(.not. has_neighbour_in_swri) then
               if(debug_show_skipped) then
                  write(stdout,'(a,i0,a,i0,a,a,a)') &
                  'Skipping (SC) orig/SFC: ', orig_atoma, '/', atoma, ' (', &
                  trim(elements(orig_atoma)%species_id),')'
               end if
               cycle
            end if
         end if

         call swri_init_centre(centrea, par, elements, atoma)

         ! jd: Calculate atomblock
         if(swri%has_metric(SW_V)) then
            call swri_vmatrixblock_sc(swri, vatomblock, centrea)
         end if
         if(swri%has_metric(SW_O)) then
            call swri_omatrixblock_sc(swri, oatomblock, centrea)
         end if

         ! qoh: Insert atomblock into sparse matrix structure.
         ! jd: NB 'atoma' must belong to this proc
         if(swri%has_metric(SW_V)) call sparse_put_block(vatomblock, &
              swri%full_metric_matrices(SW_V), atoma, atoma)
         if(swri%has_metric(SW_O)) call sparse_put_block(oatomblock, &
              swri%full_metric_matrices(SW_O), atoma, atoma)

      end do loop_A_SC

      if(swri%has_metric(SW_O)) then
         deallocate(oatomblock, stat=ierr)
         call utils_dealloc_check(myself,'oatomblock',ierr)
      end if
      if(swri%has_metric(SW_V)) then
         deallocate(vatomblock, stat=ierr)
         call utils_dealloc_check(myself,'vatomblock',ierr)
      end if

      call utils_flushed_string_output('done.'//CRLF)

      call timer_clock('swri_metric_mat_sc',2)

    end subroutine internal_calc_metric_matrices_sc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_calc_metric_matrices_oc(swri, atomblocklist, &
         elements, cell)
      !========================================================================!
      ! Calculates the off-site (off-centre) elements of the metric matrices   !
      ! by Chebyshev interpolation.                                            !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri(in/out): The SW_RI whose off-site metric matrix elements are to !
      !                 be filled. Results are returned here.                  !
      !   atomblocklist(in/out): The list of atomblocks that need to be        !
      !                          processed. This argument is modified to tag   !
      !                          blocks that have been completed (sparse_owner !
      !   elements(in): Needed to determine which atoms are within the SW_RI.  !
      !   cell(in): Simulation cell. Needed for PBCs.                          !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04, from bits and pieces of previous !
      ! versions.                                                              !
      ! Modified by James C. Womack to support multiple metric matrix element  !
      ! integral evaluation schemes, 2018.                                     !
      !========================================================================!

      use chebyshev_rep, only: CHEB_NODES, CHEB_COEFFS, &
           cheb_alloc_coeffs, cheb_dealloc_coeffs, cheb_gen_nodes_piecewise_nD
      use comms, only: pub_my_proc_id, comms_barrier, pub_total_num_procs, &
           pub_on_root
      use constants, only: SW_O, SW_V, LONG, SIZEOF_DOUBLE, SIZEOF_INT
      use hash_table, only: hash_table_init, hash_table_free
      use ion, only: ELEMENT
      use rundat, only: pub_threads_max
      use simulation_cell, only: CELL_INFO
      use utils, only: utils_flushed_string_output, utils_abort, &
           utils_alloc_check, utils_dealloc_check, utils_report_memory_estimate
      use timer, only: timer_clock

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(inout)                             :: swri
      type(ATOMBLOCK_DESCRIPTOR), intent(inout), allocatable :: atomblocklist(:)
      type(ELEMENT), intent(in)                              :: elements(:)
      type(CELL_INFO), intent(in)                            :: cell

      ! jd: Local variables
      ! --- SWs, SWpots and their Chebyshev expansions ---
      integer :: num_sws, max_sws
      type(CHEB_NODES)   :: nodes_template
      ! coeffs of a batch of SWpots (b,x) [remote]
      type(CHEB_COEFFS), allocatable :: swpot_batch_coeffs(:) ! sw_idx
      ! coeffs of a batch of SWs (b,x) [remote]
      type(CHEB_COEFFS), allocatable :: swrem_batch_coeffs(:) ! sw_idx
      ! coeffs of a batch of SWs (a,y) local to this atom [home]
      type(CHEB_COEFFS), allocatable :: swhome_all_coeffs(:) ! sw_idx
      real(kind=DP), allocatable :: swop_values(:,:) ! indexed by xproc_idx, tid

      ! --- Atomblocks, tiles ---
      real(kind=DP), allocatable :: vtile(:,:) ! Tile of V matrix
      real(kind=DP), allocatable :: otile(:,:) ! Tile of O matrix
      type(TILE_DESCRIPTOR), allocatable :: tilelist(:)
      integer :: n_atomblocks, n_tiles
      integer :: atomblock_to_do
      integer :: xbatch_to_do

      ! --- Varia ---
      logical :: all_done
      logical :: all_tiles_dispatched
      logical, allocatable :: workers_done(:)
      integer :: n_dim ! Number of dimensions of integration domain (2 or 3)
      logical :: half_domain ! If true, cut integration domain in half, i.e.
                             ! n_dim == 2, half-disc; n_dim == 3, hemisphere
      integer :: i
      integer :: ierr
      integer :: max_n_atomblocks
      integer, parameter :: maxslots_tiles = 100000
      integer, parameter :: max_tiles = 100000000 ! essentially allow unlimited growth
      ! Used to force promotion of ints in memory estimate i/o/2 avoid wraparound
      integer(kind=LONG), parameter :: L = 1_LONG
      character(len=*), parameter :: myself='internal_calc_metric_matrices_oc'
      integer :: nat


      ! -----------------------------------------------------------------------

      call timer_clock('swri_metric_mat_oc',1)

      ! rc2013: get the number of atoms from the size of elements list
      nat = size(elements)

      if(swri%has_metric(SW_V) .and. swri%has_metric(SW_O)) then
         call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
              '] - Calculating off-site elements of metric matrices V, O.'//CRLF)
      end if
      if(swri%has_metric(SW_V) .and. .not. swri%has_metric(SW_O)) then
         call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
              '] - Calculating off-site elements of the metric matrix V.'//CRLF)
      end if
      if(swri%has_metric(SW_O) .and. .not. swri%has_metric(SW_V)) then
         call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
             '] - Calculating off-site elements of the metric matrix O.'//CRLF)
      end if

      max_sws = swri%quality%max_sws_per_centre
      num_sws = swri%quality%num_sws_per_centre

      ! jd: Absolute maximum # of atomblocks (fully dense scenario). Used only
      !     to dimension arrays of handles or indices, so no worry about O(N^2).
      max_n_atomblocks = nat * (nat-1) / 2

      ! -----------------------------------------------------------------------
      ! JCW: Select number of dimensions in which we need to generate nodes
      ! JCW: and do this prior to memory estimate, so we can correctly estimate
      ! JCW: the sizes of quantities based on dimensions of quantities expanded
      ! JCW: in Chebyshev polynomials
      ! JCW: 3Dc scheme:     3 (sphere)
      ! JCW: 2Dn-1Da scheme: 2 (half-disc)
      ! -----------------------------------------------------------------------
      select case (swri%int_scheme)
         case (SCHEME_3DC)
            ! -----------------------------------------------------------------
            ! JCW: Integration in a 3D sphere
            ! -----------------------------------------------------------------
            n_dim = 3
            half_domain = .false.
         case (SCHEME_2DN_1DA)
            ! -----------------------------------------------------------------
            ! JCW: Integration on a 2D half-disc
            ! -----------------------------------------------------------------
            n_dim = 2
            half_domain = .true.
         case default
            call utils_abort('Error in '//myself//': &
                 &Unexpected integration scheme specification, ',swri%int_scheme)
      end select

      call utils_report_memory_estimate(&
           'Spherical wave res. of identity Chebyshev engine',&
           (/&
           'workspace atomblock          ',&
           'vtile                        ',&
           'otile                        ',&
           'recvtile                     ',&
           'workers_done                 ',&
           'atomblocklist                ',&
           'tilelist (root only)         ',&
           'tile hash table (upper bound)',&
           'nodes_template               ',&
           'swop_values                  ',&
           'swy_values                   ',&
           'swhome_all_coeffs            ',&
           'swpot_batch_coeffs           ',&
           'swrem_batch_coeffs           ',&
           'all_atomblocks (indices)     ',&
           'my_atomblocks (indices)      ',&
           'num_q_per_l                  ',&
           'atomblock rotation matrix    '&
           /),&
           (/&
           ! workspace atomblock
           L*max_sws*max_sws*&
           SIZEOF_DOUBLE,&
           ! vtile
           merge(L*swri%xbatchsize*max_sws*&
           SIZEOF_DOUBLE,&
           0_LONG,swri%has_metric(SW_V)),&
           ! otile
           merge(L*swri%xbatchsize*max_sws*&
           SIZEOF_DOUBLE,&
           0_LONG,swri%has_metric(SW_O)),&
           ! recvtile
           L*swri%xbatchsize*max_sws*&
           SIZEOF_DOUBLE,&
           ! workers_done
           L*pub_total_num_procs*&
           SIZEOF_INT,&
           ! atomblocklist (5 ints per entry)
           L*size(atomblocklist)*5*&
           SIZEOF_INT,&
           ! tilelist (n_tiles = n_atomblocks * n_xbatches), 3 ints per entry
           L*size(atomblocklist)*swri%n_xbatches*3*&
           SIZEOF_INT,&
           ! tile hash table. At most every worker can stall the merging of an
           ! almost full atomblock (by working on its last tile).
           L*max_sws*max_sws*pub_total_num_procs*&
           SIZEOF_DOUBLE,&
           ! nodes_template
           L*(swri%cheb_intervals*swri%cheb_order)**n_dim*SIZEOF_DOUBLE,&
           ! swop_values
           L*(swri%cheb_intervals*swri%cheb_order)**n_dim*pub_threads_max*&
           SIZEOF_DOUBLE,&
           ! swy_values
           L*(swri%cheb_intervals*swri%cheb_order)**n_dim*pub_threads_max*&
           SIZEOF_DOUBLE,&
           ! swhome_all_coeffs
           L*num_sws*&
           (swri%cheb_intervals*swri%cheb_order)**n_dim*SIZEOF_DOUBLE,&
           ! swpot_batch_coeffs
           merge(L*swri%xbatchsize*&
           (swri%cheb_intervals*swri%cheb_order)**n_dim*SIZEOF_DOUBLE,&
           0_LONG,swri%has_metric(SW_V)),&
           ! swrem_batch_coeffs
           merge(L*swri%xbatchsize*&
           (swri%cheb_intervals*swri%cheb_order)**n_dim*SIZEOF_DOUBLE,&
           0_LONG,swri%has_metric(SW_O)),&
           ! all_atomblocks
           L*3*max_n_atomblocks*SIZEOF_INT,&
           ! my_atomblocks
           L*3*max_n_atomblocks*SIZEOF_INT,&
           ! swri%num_q_per_l (get actual size, since already allocated)
           L*size(swri%num_q_per_l)*SIZEOF_INT,&
           ! per-atomblock Rmat (for 2Dn-1Da scheme), allocated and evaluated
           ! within sph_harm_rot_apply_RSH_Rmat_to_atomblock
           ! (unlike the workspace atomblock, this is only allocated to be
           ! num_sws*num_sws, i.e. does not include space for Bessels which have
           ! been excluded for having a q corresponding to an E greater than the
           ! psinc cutoff
           merge(L*num_sws*num_sws*SIZEOF_DOUBLE,&
           0_LONG,swri%int_scheme==SCHEME_2DN_1DA)&
           /))

      ! jd: Initialise the HT for storing partially merged atomblocks.
      call hash_table_init(swri%tiles_ht, 'TILES', 4, maxslots_tiles, &
           max_tiles, swri%hash_table_info_unit)

      ! jd: Allocate workspace tile for V and/or O.
      if(swri%has_metric(SW_V)) then
         allocate(vtile(swri%xbatchsize,max_sws),stat=ierr)
         call utils_alloc_check(myself,'vtile',ierr)
      end if
      if(swri%has_metric(SW_O)) then
         allocate(otile(swri%xbatchsize,max_sws),stat=ierr)
         call utils_alloc_check(myself,'otile',ierr)
      end if


      ! ---------------------------------------------------------------------
      ! JCW: Generate Chebyshev nodes for selected scheme
      ! ---------------------------------------------------------------------
      call timer_clock('swri_metric_mat_cheb_gen_nodes',1)
      call cheb_gen_nodes_piecewise_nD(nodes_template, n_dim, swri%radtable(1), &
           swri%cheb_intervals, swri%cheb_order, half_domain = half_domain)
      ! @IDENTICAL_RADII --------( swri%radtable(1) )----------------^
      call timer_clock('swri_metric_mat_cheb_gen_nodes',2)

      ! ---------------------------------------------------------------------
      ! jd: Allocate the big arrays:
      ! ---------------------------------------------------------------------
      !     - swop_values: Values of SWOPx (remote) [threaded]
      !     - swpot_batch_coeffs: Cheb coeffs of a batch of SWpotx's
      !     - swrem_batch_coeffs: Cheb coeffs of a batch of SWx's
      !     - swhome_all_coeffs: Cheb coeffs of *all* SWy's (for efficiency)
      call timer_clock('swri_metric_mat_allocs',1)

      allocate(swop_values(nodes_template%n_points_total,pub_threads_max), &
           stat=ierr)
      call utils_alloc_check(myself,'swop_values',ierr)
      allocate(swhome_all_coeffs(num_sws), stat=ierr)
      call utils_alloc_check(myself,'swhome_all_coeffs',ierr)
      if(swri%has_metric(SW_V)) then
         allocate(swpot_batch_coeffs(swri%xbatchsize), stat=ierr)
         call utils_alloc_check(myself,'swpot_batch_coeffs',ierr)
      end if
      if(swri%has_metric(SW_O)) then
         allocate(swrem_batch_coeffs(swri%xbatchsize), stat=ierr)
         call utils_alloc_check(myself,'swrem_batch_coeffs',ierr)
      end if

      ! ---------------------------------------------------------------------
      ! jd: Allocate memory for coefficients
      ! ---------------------------------------------------------------------
      do i = 1, num_sws
         call cheb_alloc_coeffs(swhome_all_coeffs(i), n_dim, &
              swri%cheb_intervals, swri%cheb_order)
      end do
      do i = 1, swri%xbatchsize
         if(swri%has_metric(SW_V)) then
            call cheb_alloc_coeffs(swpot_batch_coeffs(i), n_dim, &
                 swri%cheb_intervals, swri%cheb_order)
         end if
         if(swri%has_metric(SW_O)) then
            call cheb_alloc_coeffs(swrem_batch_coeffs(i), n_dim, &
                 swri%cheb_intervals, swri%cheb_order)
         end if
      end do
      call timer_clock('swri_metric_mat_allocs',2)

      ! ---------------------------------------------------------------------
      !  Prepare coeffs for home SWs
      ! ---------------------------------------------------------------------
      call internal_calc_home_sws(swhome_all_coeffs,nodes_template,swri)

      ! ---------------------------------------------------------------------
      !  Divide the atomblocks into tiles
      ! ---------------------------------------------------------------------
      if(pub_on_root) then
         call internal_gen_tile_list(swri, atomblocklist, tilelist, n_tiles)
      end if

      allocate(workers_done(0:pub_total_num_procs-1),stat=ierr)
      call utils_alloc_check(myself,'workers_done',ierr)

      ! *********************************************************************
      ! ---------------------------------------------------------------------
      !  Workers:
      !  - ask master for a tile to process
      !  - process a tile
      !  - send the tile to the sparse-owner of its atomblock
      !  ... in the meantime they merge received tiles into atomblocks that
      !      they sparse-own
      !  Master:
      !  - responds to tile requests and assigns tiles to workers
      !  - is also a worker
      ! ---------------------------------------------------------------------
      ! *********************************************************************
      call timer_clock('swri_metric_mat_main_loop',1)
      workers_done(:) = .false.
      all_done = .false.
      all_tiles_dispatched = .false.
      do while (.not. all_done)

         if(.not. all_tiles_dispatched) then
            ! jd: --- There are still tiles to dispatch ---
            !         so it's crucial that master keeps
            !         satisfying requests.

            ! jd: Request a tile.
            !     [Merging takes place during waits]
            !     [Master satisfies requests]
            call internal_request_tile(atomblock_to_do, xbatch_to_do, & ! out
                 workers_done, swri, &                                  ! inout
                 atomblocklist, tilelist, n_tiles)                      ! in

            if(.not. (atomblock_to_do == -1 .and. xbatch_to_do == -1)) then
               ! jd: Process a tile
               !     [No comms here] -> except master keeps dispatching
               call internal_tile_by_chebs(vtile, otile, &
                    swhome_all_coeffs, &
                    swpot_batch_coeffs, swrem_batch_coeffs, &
                    nodes_template, swop_values, elements, &
                    atomblocklist(atomblock_to_do)%atoma, &
                    atomblocklist(atomblock_to_do)%atomb, xbatch_to_do, cell, &
                    swri, &
                    workers_done, atomblocklist, tilelist, n_tiles)
#ifdef MPI
               ! jd: Send to owner
               !     [Merging takes place during waits]
               call timer_clock('swri_metric_mat_wait_send',1)
               call internal_send_tile_to_owner(swri, cell, atomblocklist, &
                    vtile, otile, atomblock_to_do, xbatch_to_do, &
                    atomblocklist(atomblock_to_do)%atoma_sparse_owner, max_sws)
               call timer_clock('swri_metric_mat_wait_send',2)
#else
               call internal_merge_serialmode(swri, cell, atomblocklist, &
                    atomblock_to_do, xbatch_to_do, vtile, otile)
#endif
               ! jd: On master we must manually track completed tiles, as it
               !     doesn't request tiles the way other workers do.
               if(pub_on_root) then
                  tilelist(swri%n_xbatches * (atomblock_to_do-1) &
                       + xbatch_to_do)%worker_proc = TILE_COMPLETED
               end if

            else
               ! jd: Master signalled there are no more tiles to dispatch.
               all_tiles_dispatched = .true.
            end if
         else
            ! jd: --- All tiles dispatched ---
            !         Only merging left

            ! Master keeps answering requests, since, if it happens that it
            ! dealt the last tile to itself, only it knows that all tiles have
            ! been dispatched, workers will still come once with a request
            ! before they figure it out. So master must be the last to leave
            ! the loop.
            if(pub_on_root) then
               call internal_dispatch_tile(workers_done, tilelist, n_tiles)
            end if
            call internal_merge_attempt(swri, cell, atomblocklist)

            all_done = .not. (any(atomblocklist(:)%atoma_sparse_owner == &
                 pub_my_proc_id))

            ! Master can only leave after everyone's done
            if(pub_on_root) then
               if(.not. all(workers_done)) all_done = .false.
            end if

         end if
      end do
      call timer_clock('swri_metric_mat_main_loop',2)

      call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
          '] - Filling the remainder by symmetry.'//CRLF)

      deallocate(workers_done,stat=ierr)
      call utils_dealloc_check(myself,'workers_done',ierr)

      ! jd: Free hash table
      call hash_table_free(swri%tiles_ht)

      ! ---------------------------------------------------------------------
      ! jd: Deallocate arrays, Cheb coefficients, block- and tile-lists
      ! ---------------------------------------------------------------------
      if(pub_on_root) then
         deallocate(tilelist, stat=ierr)
         call utils_dealloc_check(myself,'tilelist',ierr)
      end if

      deallocate(swop_values,stat=ierr)
      call utils_dealloc_check(myself,'swop_values',ierr)

      if(swri%has_metric(SW_O)) then
         do i = 1, swri%xbatchsize
            call cheb_dealloc_coeffs(swrem_batch_coeffs(i))
         end do
      end if
      if(swri%has_metric(SW_V)) then
         do i = 1, swri%xbatchsize
            call cheb_dealloc_coeffs(swpot_batch_coeffs(i))
         end do
      end if

      do i = 1, swri%quality%num_sws_per_centre
         call cheb_dealloc_coeffs(swhome_all_coeffs(i))
      end do

      if(swri%has_metric(SW_O)) then
         deallocate(swrem_batch_coeffs, stat=ierr)
         call utils_dealloc_check(myself,'swrem_batch_coeffs',ierr)
         deallocate(otile, stat=ierr)
         call utils_dealloc_check(myself,'otile',ierr)
      end if
      if(swri%has_metric(SW_V)) then
         deallocate(swpot_batch_coeffs, stat=ierr)
         call utils_dealloc_check(myself,'swpot_batch_coeffs',ierr)
         deallocate(vtile, stat=ierr)
         call utils_dealloc_check(myself,'vtile',ierr)
      end if
      deallocate(swhome_all_coeffs, stat=ierr)
      call utils_dealloc_check(myself,'swhome_all_sph_coeffs',ierr)

      call timer_clock('swri_metric_mat_oc_work_imbalance',1)
      call comms_barrier
      call timer_clock('swri_metric_mat_oc_work_imbalance',2)

      ! ---------------------------------------------------------------------
      ! jd: Fill the second triangle of metric matrices by symmetry.
      ! ---------------------------------------------------------------------
      call internal_oc_fill_by_symmetry(swri)

      ! ---------------------------------------------------------------------
      ! jd: And we're done
      ! ---------------------------------------------------------------------
      call timer_clock('swri_metric_mat_oc_work_imbalance',1)
      call comms_barrier
      call timer_clock('swri_metric_mat_oc_work_imbalance',2)

      call timer_clock('swri_metric_mat_oc',2)

    end subroutine internal_calc_metric_matrices_oc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_oc_fill_by_symmetry(swri)
      !========================================================================!
      ! Fills the second triangle of metric matrices by symmetry.              !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri(in/out): The SW_RI whose off-site metric matrix elements are to !
      !                 be filled. Results are returned here.                  !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04, by extraction from               !
      ! internal_calc_metric_matrices_oc.                                      !
      !========================================================================!

      use comms, only: pub_my_proc_id
      use constants, only: SW_V, SW_O
      use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_transpose, &
           sparse_destroy, sparse_put_block
      use timer, only: timer_clock
      use utils, only: utils_dealloc_check

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(inout)                          :: swri

      ! jd: Local variables
      real(kind=DP), allocatable :: zero_atomblock(:,:) ! For zeroing diagonal
      type(SPAM3) :: mat_transpose  ! Workspace for transposed metric matrix
      integer :: local_a, atoma
      integer :: ierr
      character(len=*), parameter :: myself = 'internal_oc_fill_by_symmetry'

      ! -----------------------------------------------------------------------

      call timer_clock('swri_metric_mat_symmetry',1)

      allocate(zero_atomblock(&
           swri%quality%max_sws_per_centre,swri%quality%max_sws_per_centre), &
           stat=ierr)
      call utils_dealloc_check(myself,'zero_atomblock',ierr)
      zero_atomblock = 0.0_DP
      if(swri%has_metric(SW_V)) then
         call sparse_create(mat_transpose,swri%full_metric_matrices(SW_V))
         call sparse_transpose(mat_transpose,swri%full_metric_matrices(SW_V))
         ! qoh: Zero diagonal atomblocks in transpose
         do local_a=1, par%num_atoms_on_proc(pub_my_proc_id)
            atoma = par%first_atom_on_proc(pub_my_proc_id) + local_a -1
            ! jd: NB: No ignoring of atoms should happen here
            call sparse_put_block(zero_atomblock, mat_transpose, atoma, atoma)
         end do
         call sparse_axpy(swri%full_metric_matrices(SW_V),mat_transpose,1.0_DP)
         call sparse_destroy(mat_transpose)
      end if
      if(swri%has_metric(SW_O)) then
         call sparse_create(mat_transpose,swri%full_metric_matrices(SW_O))
         call sparse_transpose(mat_transpose,swri%full_metric_matrices(SW_O))
         ! qoh: Zero diagonal atomblocks in transpose
         do local_a=1, par%num_atoms_on_proc(pub_my_proc_id)
            atoma = par%first_atom_on_proc(pub_my_proc_id) + local_a -1
            ! jd: NB: No ignoring of atoms should happen here
            call sparse_put_block(zero_atomblock, mat_transpose, atoma, atoma)
         end do
         call sparse_axpy(swri%full_metric_matrices(SW_O),mat_transpose,1.0_DP)
         call sparse_destroy(mat_transpose)
      end if
      deallocate(zero_atomblock, stat=ierr)
      call utils_dealloc_check(myself,'zero_atomblock',ierr)

      call timer_clock('swri_metric_mat_symmetry',2)

    end subroutine internal_oc_fill_by_symmetry

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_merge_serialmode(swri, cell, atomblocklist, &
         atomblock_to_do, xbatch_to_do, vtile, otile)
      !========================================================================!
      ! A non-MPI version of the merge-attempt subroutine.                     !
      ! In the absence of MPI the usual scenario of sending completed tiles    !
      ! between subroutines through MPI messages does not apply. This is the   !
      ! serial mode version.                                                   !
      ! This subroutine adds the tiles to the hash table, and if all tiles for !
      ! a block are completed it calls internal_merge_my_tiles() and marks the !
      ! atomblock as completed.                                                !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      !========================================================================!

      use constants, only: SW_V, SW_O
      use hash_table, only: hash_table_add, hash_table_probe
      use simulation_cell, only: CELL_INFO

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(inout)                             :: swri
      type(ATOMBLOCK_DESCRIPTOR), intent(inout), allocatable :: atomblocklist(:)
      type(CELL_INFO), intent(in)                            :: cell
      integer :: atomblock_to_do, xbatch_to_do
      real(kind=DP), allocatable :: vtile(:,:) ! Tile of V matrix
      real(kind=DP), allocatable :: otile(:,:) ! Tile of O matrix

      ! jd: Local variables
      integer :: tile_length
      integer :: atoma, atomb
      integer :: V_or_O
      integer :: cur_xbatch
      logical :: have_entire_atomblock

      ! -----------------------------------------------------------------------

      tile_length = swri%xbatchsize * swri%quality%max_sws_per_centre
      atoma = atomblocklist(atomblock_to_do)%atoma
      atomb = atomblocklist(atomblock_to_do)%atomb

      ! jd: Store tile in hash table.
      if(swri%has_metric(SW_V)) then
         call hash_table_add(swri%tiles_ht, reshape(vtile,(/tile_length/)), &
              tile_length, atoma, atomb, xbatch_to_do, SW_V, &
              overfill_strategy = 'F')
      end if
      if(swri%has_metric(SW_O)) then
         call hash_table_add(swri%tiles_ht, reshape(otile,(/tile_length/)), &
              tile_length, atoma, atomb, xbatch_to_do, SW_O, &
              overfill_strategy = 'F')
      end if

      ! jd: See if we now have all tiles for this atomblock.
      do V_or_O = SW_V, SW_O
         have_entire_atomblock = .true.
         do cur_xbatch = 1, swri%n_xbatches
            if(-1 == hash_table_probe(swri%tiles_ht, atoma, atomb, &
                 cur_xbatch, V_or_O)) then
               have_entire_atomblock = .false.
               exit
            end if
         end do

         ! jd --- Merging ---
         if(have_entire_atomblock) then
            call internal_merge_my_tiles(swri, cell, V_or_O, atoma, atomb)
            ! jd: Mark block as done
            if(swri%has_metric(SW_V) .and. swri%has_metric(SW_O)) then
               ! jd: Corner case, where we have both V and O matrices. Then
               !     only mark as done once the latter is done.
               if(V_or_O == SW_O) then
                  atomblocklist(atomblock_to_do)%atoma_sparse_owner = -1
               end if
            else
               atomblocklist(atomblock_to_do)%atoma_sparse_owner = -1
            end if
         end if
      end do

    end subroutine internal_merge_serialmode

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_send_tile_to_owner(swri, cell, atomblocklist, &  ! inout
        vtile, otile, atomblock, xbatch, owner_proc, max_sws)            ! in
      !========================================================================!
      ! Sends completed tile (or tiles if doing both the electrostatic and     !
      ! overlap matrix) to the proc that sparse-owns the atomblock.            !
      !                                                                        !
      ! This is a first attempt and is suboptimal. Since the workspace in      !
      ! vtile and otile is needed for subsequent tiles, we must currently mimic!
      ! a blocking send. At least we keep merging tiles while busy-waiting on  !
      ! comms_test, but waits are unavoidable. Calls to the merge-attempt      !
      ! routine are interleaved with computation in internal_tile_by_chebs(),  !
      ! so the waits should never be too long.                                 !
      !                                                                        !
      ! A better approach would be to construct a queue of pending sends, with !
      ! an array of vtiles and otiles in transfer. Keeping it simple for now.  !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri(in/out): As merging takes place during busy waits, we need to   !
      !                 pass the SW_RI -- it holds the metric matrices.        !
      !   cell(in): Simulation cell. Needed for PBCs.                          !
      !   atomblockslist(in/out): Needs to be 'out', because the merge routine !
      !                           tags atomblocks as done when they're done.   !
      !                           Needed for tile merging.                     !
      !   vtile, otile(in): Tiles that are to be sent. Depending on which met- !
      !                     rics are active in the SW_RI object one of these   !
      !                     arrays might not be referenced.                    !
      !   atomblock(in): The atomblock to which this tile belongs.             !
      !   xbatch(in): The x-batch of the atomblock which this tile represents. !
      !   owner_proc(in): The proc where we want the tile to be sent to.       !
      !   max_sws(in): The maximum number of SWs in the SW_RI. Needed for      !
      !                dimensioning.                                           !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      !========================================================================!

      use comms, only: comms_send, comms_test

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(inout)                             :: swri
      type(CELL_INFO), intent(in)                            :: cell
      type(ATOMBLOCK_DESCRIPTOR), intent(inout), allocatable :: atomblocklist(:)
      real(kind=DP), intent(in), allocatable                 :: vtile(:,:)
      real(kind=DP), intent(in), allocatable                 :: otile(:,:)
      integer, intent(in)                                    :: atomblock
      integer, intent(in)                                    :: xbatch
      integer, intent(in)                                    :: owner_proc
      integer, intent(in)                                    :: max_sws

      ! jd: Local variables
      integer :: tile_length
      integer :: tile_id
      integer :: send_handle1, send_handle2
      logical :: send1_completed, send2_completed

      ! -----------------------------------------------------------------------
#ifdef MPI
      tile_length = swri%xbatchsize * max_sws
      tile_id = swri%n_xbatches * (atomblock-1) + xbatch

      if(swri%has_metric(SW_V)) then
         call comms_send(owner_proc, tile_id, &
              length = 1, tag = TAG_TILE_ID, &
              return_handle = send_handle1, add_to_stack = .false.)
         call comms_send(owner_proc, vtile, &
              length = tile_length, tag = TAG_TILE_DATA, &
              return_handle = send_handle2, add_to_stack = .false.)
         send1_completed = .false.
         send2_completed = .false.
         do while (.not. (send1_completed .and. send2_completed))
            call internal_merge_attempt(swri, cell, atomblocklist)
            call comms_test(send1_completed, send_handle1)
            call comms_test(send2_completed, send_handle2)
         end do
      end if

      if(swri%has_metric(SW_O)) then
         tile_id = -tile_id ! jd: Sign distinguishes O tiles from V tiles.
         call comms_send(owner_proc, tile_id, &
              length = 1, tag = TAG_TILE_ID, &
              return_handle = send_handle1, add_to_stack = .false.)
         call comms_send(owner_proc, otile, &
              length = tile_length, tag = TAG_TILE_DATA, &
              return_handle = send_handle2, add_to_stack = .false.)
         send1_completed = .false.
         send2_completed = .false.
         do while (.not. (send1_completed .and. send2_completed))
            call internal_merge_attempt(swri, cell, atomblocklist)
            call comms_test(send1_completed, send_handle1)
            call comms_test(send2_completed, send_handle2)
         end do
      end if
#endif
    end subroutine internal_send_tile_to_owner

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_tile_by_chebs(vtile, otile, &
         swhome_all_coeffs, swpot_batch_coeffs, swrem_batch_coeffs,&
         nodes_template, swop_values, elements, atoma, atomb, xbatch, cell, &
         swri, &
         workers_done, atomblocklist, tilelist, n_tiles)
      !========================================================================!
      ! Calculates a tile of an atomblock of the V matrix (and/or the O matrix)!
      ! by Chebyshev interpolation. This is only for off-site atomblocks, as   !
      ! on-site (same-centre) blocks can be evaluated analytically.            !
      !                                                                        !
      ! OpenMP threading is used over batch elements.                          !
      !                                                                        !
      ! Evaluation of off-site metric matrix elements may be evaluated using   !
      ! different schemes, depending on the value of swri%int_scheme.          !
      !                                                                        !
      ! Available integral evaluation schemes:                                 !
      !                                                                        !
      ! 3Dc (swri%int_scheme == SCHEME_3DC)                                    !
      ! ===================================                                    !
      ! The original scheme, 3-D Chebyshev (3Dc), evaluates the matrix elements!
      ! by integration over the product of home and remote SWOPs expanded in   !
      ! in 3-D, in the localization sphere of the home SWOP. The integration   !
      ! is performed in a Cartesian coordinate system with axes parallel with  !
      ! ONETEP's global Cartesian coordinate system.                           !
      !                                                                        !
      ! This is the scheme described in section II.D.1                         !
      !   J. Dziedzic et al., J. Chem. Phys. 139, 214103 (2013)                !
      !                                                                        !
      ! 2Dn-1Da (swri%int_scheme == SCHEME_2DN_1DA)                            !
      ! ===========================================                            !
      ! The newer scheme, 2-D numerical, 1-D analytic (2Dn-1Da), evaluates     !
      ! matrix elements in an "integral coordinate frame", in which the z-axis !
      ! is aligned along the vector between the two atoms. The real-spherical  !
      ! harmonic (RSH) components of home and remote SWOPs are aligned along   !
      ! this z-axis, which enables the metric matrix integrals to be separated !
      ! into a 2-D numerical integral over the radial and polar coordinates    !
      ! of a spherical polar coordinate system centred on the atom in the      !
      ! home sphere and a 1-D analytic integral over the azimuthal coordinate. !
      !                                                                        !
      ! The 2-D numerical integral is performed by integration over an         !
      ! home and remote SWOPRTs expanded in Chebyshev polynomials on a 2-D     !
      ! half-disc within the home sphere, in a plane parallel to the integral  !
      ! coordinate frame z-axis (a "SWOPRT" is the (r,theta)-dependent part of !
      ! a SWOP). The matrix elements are obtained by taking the product of the !
      ! 2-D numerical integral with the 1-D analytic integral (this is done    !
      ! within this subroutine).                                               !
      !                                                                        !
      ! While the separation of the integral into 2-D and 1-D components       !
      ! occurs in spherical polar coordinates, the 2-D integration over        !
      ! SWOPRTs expanded in Chebyshev polynomials is performed in Cartesian    !
      ! coordinates (within the half-disc). This change of coordinates         !
      ! is accounted for by Jacobian determinants included in the expansion    !
      ! of the SWRT in the home sphere (half-disc), see documentation for      !
      ! swri_eval_swrt_at_nodes_on_halfdisc_fast and                           !
      ! swri_eval_swrt_at_nodes_on_halfdisc.                                   !
      !                                                                        !
      ! Atomblock tiles are stored in the integral coordinate frame, and are   !
      ! transformed to the ONETEP coordinate frame (with SWOPs aligned         !
      ! parallel to the global Cartesian coordinate system z-axis) following   !
      ! assembly into a full atomblock in internal_merge_my_tiles by           !
      ! application of a RSH rotation matrix.                                  !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   vtile, otile (in/out): Atomblock tiles to be filled. They are always !
      !                          first overwritten with garbage (this is cheap !
      !                          compared to the cost of the interpolation) to !
      !                          aid with debugging.                           !
      !   swhome_all_coeffs (in):     Chebyshev expansion coefficients for the !
      !                               home sphere (the local SW expanded).     !
      !                               This holds all coefficients.             !
      !   swpot_batch_coeffs (in/out):     Workspace, allocated in parent, for !
      !                                    holding expansion coefficients of   !
      !                                    the remote SW potential (for the V  !
      !                                    metric.                             !
      !                                    This holds an x-batch of coeffs.    !
      !   swrem_batch_coeffs (in/out):     Workspace, allocated in parent, for !
      !                                    holding expansion coefficients of   !
      !                                    remote SW (for the O metric).       !
      !                                    This holds an x-batch of coeffs.    !
      !   nodes_template (in):     Template for the positions of Chebyshev     !
      !                            nodes, passed from parent.                  !
      !   swop_values (inout): Workspace, allocated in parent, for holding the !
      !                        SWOP values at Chebyshev nodes.                 !
      !   elements (in): Needed, in swri_init_centre(), to get positions of the!
      !                  atomic cores.                                         !
      !   atoma, atomb (in): Indices of the atoms making up the atomblock.     !
      !   xbatch (in): Number of x-batch to process.                           !
      !   cell (in): Simulation cell. Needed for PBCs.                         !
      !   swri (in/out): The SW_RI in which we're operating. Must be in/out,   !
      !                  because the interleaved comms perform tile merging.   !
      !                                                                        !
      !   The following arguments are only needed for the interleaved comms.   !
      !   ------------------------------------------------------------------   !
      !   workers_done (in/out): Keeps track of which workers have completed   !
      !                          their share of atomblocks.                    !
      !   atomblocklist (in/out): Details about the atomblocks to compute.     !
      !   tilelist (in/out): Details about the tiles to compute.               !
      !   n_tiles (in/out): Number of elements in tilelist.                    !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Caveats:                                                               !
      !   - Arguments are explicitly passed from the parent routine to work    !
      !     around a PGI/OMP bug whereby arrays in SHARED clauses in internal  !
      !     procedures are spuriously deallocated (TPR 21390).                 !
      !   - All tiles are dimensioned as (swri%xbatchsize, max_sws). What is   !
      !     filled by Chebyshevs here is limited to a chunk that can be smaller!
      !     in either or both of these dimensions. For the first dimension this!
      !     is because the last batch can be smaller (if num_sws is not divisi-!
      !     ble by the batch size. For the second dimension, we only fill up   !
      !     to num_sws). Since we are responsible for initialising the entire  !
      !     array, the rest will contain the 'garbage_real' constant (GR).     !
      !                                                                        !
      !     |<------------------ max_sws ------------------>|                  !
      !     |< ------------ num_sws ------------>|                             !
      !  -  -------------------------------------------------                  !
      !  ^  |                                    .          | TILE 1           !
      !  |n |                 DATA               .    GR    |                  !
      !  |u |                                    .          |                  !
      !  |m -------------------------------------|-----------                  !
      !  |_ |                                    .          | TILE 2           !
      !  |s |                 DATA               .    GR    |                  !
      !  |w |                                    .          |                  !
      !  |s -------------------------------------|-----------                  !
      !  v  |                 DATA               .    GR    | TILE 3           !
      !  -  |....................................+..........|                  !
      !     |                  GR                .    GR    |                  !
      !     -------------------------------------------------                  !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      ! Generalised to add support for multiple off-site metric matrix element !
      ! evaluation schemes by James Womack, 2018.                              !
      !========================================================================!

      use comms, only: pub_on_root
      use constants, only: garbage_real
      use chebyshev_rep, only: CHEB_COEFFS, CHEB_NODES, &
           cheb_int_product_piecewise_nD, cheb_expansion_piecewise_nD
      use geometry, only: POINT, operator(-), geometry_magnitude
      use rundat, only: pub_hfx_bc_is_periodic, pub_threads_max
      use simulation_cell, only: CELL_INFO, minimum_image_displacement
      use timer, only: timer_clock
      use utils, only: utils_abort, utils_assert

      implicit none

      ! jd: Arguments
      real(kind=DP), allocatable, intent(inout)     :: vtile(:,:)
      real(kind=DP), allocatable, intent(inout)     :: otile(:,:)
      type(CHEB_COEFFS), allocatable, intent(in)    :: swhome_all_coeffs(:)
      type(CHEB_COEFFS), allocatable, intent(inout) :: swpot_batch_coeffs(:)
      type(CHEB_COEFFS), allocatable, intent(inout) :: swrem_batch_coeffs(:)
      type(CHEB_NODES), intent(in)       :: nodes_template
      real(kind=DP), intent(inout)       :: swop_values(:,:) ! (xproc_idx,thread)
      type(ELEMENT), intent(in)          :: elements(par%nat)
      integer, intent(in)                :: atoma
      integer, intent(in)                :: atomb
      integer, intent(in)                :: xbatch
      type(CELL_INFO), intent(in)        :: cell
      type(SW_RI), intent(inout)         :: swri

      ! jd: Arguments needed for comms
      logical, intent(inout), allocatable                    :: workers_done(:)
      type(ATOMBLOCK_DESCRIPTOR), intent(inout), allocatable :: atomblocklist(:)
      type(TILE_DESCRIPTOR), intent(inout), allocatable      :: tilelist(:)
      integer, intent(in)                                    :: n_tiles

      ! jd: Local variables
      type(ATOM_CENTRE) :: centrea
      type(ATOM_CENTRE) :: centreb
      integer :: sbidx
      integer :: lval
      integer :: mval
      integer :: swx, swy, n_swxy, swxy
      integer :: mvalx, lvalx, xsbidx
      type(POINT) :: disp, disp_int_frame
      real(kind=DP) :: length_disp
      integer :: swx_first, swx_last, swx_in_batch
      integer :: besselx
      integer :: lx, mx, my
      real(kind=DP) :: qx, ax
      logical :: otile_is_zero ! .true., iff SWs x and y do not overlap
      real(kind=DP) :: int_1Da ! value of 1-D analytic int in 2Dn-1Da scheme
      logical :: int_1Da_is_zero ! .true. in 2Dn-1Da scheme where 1-D
                                 ! analytic integral is zero (no need to do
                                 ! 2-D numerical integration in this case)
!$    integer, external :: omp_get_thread_num
      integer :: tid
      character(len=*), parameter :: myself = 'internal_tile_by_chebs'

      ! -----------------------------------------------------------------------

      call utils_trace_in(myself)
      call timer_clock(myself,1)

      if(swri%has_metric(SW_V)) vtile(:,:) = garbage_real
      if(swri%has_metric(SW_O)) otile(:,:) = garbage_real

      call swri_init_centre(centrea, par, elements, atoma)
      call swri_init_centre(centreb, par, elements, atomb)

      call utils_assert(atoma /= atomb,'Internal error [1] in '//trim(myself))

      if(.not. any(pub_hfx_bc_is_periodic)) then
         ! jd: OBC
         disp = centreb%incell - centrea%incell
         length_disp = geometry_magnitude(disp)
      else if(all(pub_hfx_bc_is_periodic)) then
         ! jd: PBC
         call minimum_image_displacement(length_disp, disp, &
              centreb%incell, centrea%incell, cell)
      else
         call utils_abort(myself//': Mixed BCs are not currently supported.')
      end if

      ! JCW: Compute displacement in integral coordinate frame if needed
      if (swri%int_scheme == SCHEME_2DN_1DA) then
         ! JCW: Displacement in the integral frame (for 2Dn-1Da scheme) is
         ! JCW: along z-axis.
         ! JCW: centreb is swx, i.e. the SWOP in non-home sphere
         ! JCW: centrea is swy, i.e. the SWOP in home sphere
         ! JCW: centreb - centrea is the position of the centre of swx relative
         ! JCW: in a coordinate system centred on swy, so the displacement in the
         ! JCW: integral frame (with z-axis pointing from swy to swx) is
         ! JCW: ( 0, 0, |disp| )
         disp_int_frame%X = 0.0_DP
         disp_int_frame%Y = 0.0_DP
         disp_int_frame%Z = length_disp
      end if

      if(length_disp > centreb%radius + centrea%radius) then
         ! jd: These two SWs don't overlap
         otile_is_zero = .true.
      else
         otile_is_zero = .false.
      end if

      swx_first = (xbatch-1) * swri%xbatchsize + 1
      swx_last = swx_first + swri%xbatchsize - 1
      if(swx_last > swri%quality%num_sws_per_centre) then
         swx_last=swri%quality%num_sws_per_centre
      end if

      ! -----------------------------------------------------------
      ! Prepare coeffs for remote SWOPs
      ! -----------------------------------------------------------

      ! jd: Loop over SWOPs in this batch
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(swx, swx_in_batch, besselx, mx, lx, qx, ax, tid) &
!$OMP SHARED(swx_first, swx_last, nodes_template, disp, disp_int_frame, &
!$OMP      swri, otile_is_zero, &
!$OMP      centreb, swpot_batch_coeffs, &
!$OMP      swrem_batch_coeffs, swop_values, &
!$OMP      pub_threads_max, workers_done, atomblocklist, tilelist, n_tiles, &
!$OMP      pub_on_root, cell)

      do swx = swx_first, swx_last

         swx_in_batch = swx - swx_first + 1

         besselx = swri%sw_idx_to_bessel_idx(swx)
         mx = swri%sw_idx_to_m(swx)

         lx = swri%sphbessels(centreb%species_number,besselx)%lval
         qx = swri%sphbessels(centreb%species_number,besselx)%qval
         ax = swri%sphbessels(centreb%species_number,besselx)%aval

         tid = 0
!$       tid = omp_get_thread_num();

         ! jd: --------------------------------------------------------------
         !     Interleave with comms to minimise waits in:
         !     - internal_request_tile()
         !     - internal_send_tile_to_owner().
         ! jd: --------------------------------------------------------------
         if(tid == 0) then
            if(pub_on_root) then
               call internal_dispatch_tile(workers_done, tilelist, n_tiles)
            end if
            call internal_merge_attempt(swri, cell, atomblocklist)
         end if
         ! jd: The above is susceptible to a corner case where on some or all
         !     procs thread 0 is unlucky and does not get any loop iterations
         !     scheduled to it (eg. there are more threads than the batch size),
         !     thus never getting to do the comms. But if this is the case, the
         !     wait here is very short anyway, and we'll have a chance later,
         !     so no deadlocks.

         ! JCW: Evaluate and expand SWOP (or SWOPRT) from x acting on y
         select case (swri%int_scheme)
            case (SCHEME_3DC)
               ! JCW: Evaluate and expand SWOP on nodes in sphere
               if(swri%has_metric(SW_V)) then
                  call timer_clock('swri_metric_mat_eval_swpot',1)
                  call swri_eval_swpot_at_nodes(swri, &    ! in
                       swop_values(:,tid+1), &              ! out
                       nodes_template, lx, mx, disp, &
                       centreb%species_number, besselx)
                  call timer_clock('swri_metric_mat_eval_swpot',2)
                  call timer_clock('swri_metric_mat_cheb_expand',1)
                  call cheb_expansion_piecewise_nD( &
                       swpot_batch_coeffs(swx_in_batch), &
                       nodes_template, swop_values(:,tid+1))
                  call timer_clock('swri_metric_mat_cheb_expand',2)
               end if

               ! jd: Evaluate and expand the SW from x acting on y
               if(swri%has_metric(SW_O) .and. .not. otile_is_zero) then
                  call timer_clock('swri_metric_mat_eval_sw_remote',1)
                  call swri_eval_sw_at_nodes(swop_values(:,tid+1), & ! out
                       cell, nodes_template, lx, mx, qx, ax, disp)   ! in
                  call timer_clock('swri_metric_mat_eval_sw_remote',2)
                  call timer_clock('swri_metric_mat_cheb_expand',1)
                  call cheb_expansion_piecewise_nD( &
                       swrem_batch_coeffs(swx_in_batch), &
                       nodes_template, swop_values(:,tid+1))
                  call timer_clock('swri_metric_mat_cheb_expand',2)
               end if

            case (SCHEME_2DN_1DA)
               ! JCW: Evaluate and expand SWOPRT (SWOP (r,theta)-dependent part)
               ! JCW: at nodes on a 2-D half-disc in the "integral frame"

               if(swri%has_metric(SW_V)) then
                  call timer_clock('swri_metric_mat_eval_swpot',1)
                  call swri_eval_swpotrt_at_nodes_on_halfdisc(swri, &    ! in
                       swop_values(:,tid+1), &                       ! out
                       nodes_template, lx, mx, disp_int_frame, &
                       centreb%species_number, besselx)
                  call timer_clock('swri_metric_mat_eval_swpot',2)
                  call timer_clock('swri_metric_mat_cheb_expand',1)
                  call cheb_expansion_piecewise_nD( &
                       swpot_batch_coeffs(swx_in_batch), &
                       nodes_template, swop_values(:,tid+1))
                  call timer_clock('swri_metric_mat_cheb_expand',2)
               end if

               ! @2DN_1DA_OVERLAP
               if (swri%int_scheme == SCHEME_2DN_1DA) then
                  call utils_assert(.not.swri%has_metric(SW_O),&
                       "Current implementation of 2Dn-1Da integration scheme &
                       &does not support the overlap metric, sorry.")
               end if

               !! jd: Evaluate and expand the SW from x acting on y
               !if(swri%has_metric(SW_O) .and. .not. otile_is_zero) then
               !   call timer_clock('swri_metric_mat_eval_sw_remote',1)
               !   call swri_eval_sw_at_nodes(swop_values(:,tid+1), & ! out
               !        nodes_template, lx, mx, qx, ax, disp)     ! in
               !   call timer_clock('swri_metric_mat_eval_sw_remote',2)
               !   call timer_clock('swri_metric_mat_cheb_expand',1)
               !   call cheb_expansion_piecewise_nD( &
               !        swrem_batch_coeffs(swx_in_batch), &
               !        nodes_template, swop_values(:,tid+1))
               !   call timer_clock('swri_metric_mat_cheb_expand',2)
               !end if
            case default
               ! JCW: Unknown integration scheme
               call utils_abort('Error in '//myself//': &
                    &Unexpected integration scheme specification, ',&
                    swri%int_scheme)
         end select
      end do
!$OMP END PARALLEL DO

      ! ------------------------------------------
      ! Perform Chebyshev integration of products
      ! ------------------------------------------
      ! jd: This used to be a double loop structure, with swx going over a
      !     batch in the outer loop, and swy going over all num_sws in the
      !     inner loop. I rewrote it to use a single loop in 2016.04 to
      !     - better load-balance threads (x-batches tend to be rather small,
      !       and we typically do just one triangle of the entire matrix, so
      !       different swx's have different cost),
      !     - do the comms polling more often, to reduce wait times.

      n_swxy = (swx_last-swx_first+1)         * swri%quality%num_sws_per_centre
      !        ^ outer loop over swx            ^ inner loop over swy
      !        ^ == xbatches, except in the last batch

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(swxy, swx, swx_in_batch, swy, int_1Da, int_1Da_is_zero, &
!$OMP      mx, my, tid) &
!$OMP SHARED(n_swxy, swx_first, swx_last, nodes_template, &
!$OMP      swpot_batch_coeffs, otile_is_zero, otile, vtile, &
!$OMP      swhome_all_coeffs, swrem_batch_coeffs, cell, &
!$OMP      swri, pub_threads_max, workers_done, atomblocklist, tilelist, &
!$OMP      n_tiles, pub_on_root, xbatch, atoma, atomb,centrea,centreb)
      do swxy = 1, n_swxy

         tid = 0
!$       tid = omp_get_thread_num();

         ! jd: --------------------------------------------------------------
         !     Interleave with comms to minimise waits in:
         !     - internal_request_tile()
         !     - internal_send_tile_to_owner().
         ! jd: --------------------------------------------------------------
         if(tid == 0) then
            if(pub_on_root) then
               call internal_dispatch_tile(workers_done, tilelist, n_tiles)
            end if
            call internal_merge_attempt(swri, cell, atomblocklist)
         end if
         ! jd: The above is technically susceptible to a corner case where on
         !     some or all procs thread 0 is unlucky and does not get any loop
         !     iterations scheduled to it (eg. there are more threads than
         !     n_swxy), thus never getting to do the comms.
         !     Since n_swxy is rather large (~1000 typically), this is unlikely
         !     to happen in practice. If it does, that means that the wait here
         !     is very short anyway, and we'll have a chance later (outside of
         !     this subroutine), so no deadlocks -- only very occasional hiccups
         !     in very rarely used scenarios.

         ! jd: Index arithmetics to get swx, swy, swx_in_batch from swxy
         swx = swxy / swri%quality%num_sws_per_centre + swx_first
         swy = mod(swxy, swri%quality%num_sws_per_centre)
         if(swy == 0) then
            swy = swri%quality%num_sws_per_centre
            swx = swx - 1
         end if
         call utils_assert(swx >= swx_first .and. swx <= swx_last, &
              myself//': Broken index arithmetics [1]', swx, swx_last, swxy, xbatch)
         call utils_assert(swy >=1 .and. &
              swy <= swri%quality%num_sws_per_centre, &
              myself//': Broken index arithmetics [2]', swy, &
              swri%quality%num_sws_per_centre, swxy, xbatch)
         swx_in_batch = swx - swx_first + 1

         ! jd: If consistency is forced, only half of the block is needed
         if(force_metric_matrix_consistency .and. swx > swy) then
            if(swri%has_metric(SW_V)) vtile(swx_in_batch,swy) = 0.0_DP
            if(swri%has_metric(SW_O)) otile(swx_in_batch,swy) = 0.0_DP
            cycle
         end if

         call timer_clock('swri_metric_mat_cheb_int_product',1)
         ! JCW: Evaluate matrix element using selected integration scheme
         select case (swri%int_scheme)
            case (SCHEME_3DC)
               ! JCW: 3-D integration over a sphere in ONETEP's Cartesian
               ! JCW: Cartesian coordinate system

               ! jd: In the order of 30ms per call
               if(swri%has_metric(SW_V)) then
                  vtile(swx_in_batch,swy) = cheb_int_product_piecewise_nD( &
                       swhome_all_coeffs(swy), &
                       swpot_batch_coeffs(swx_in_batch), &
                       nodes_template)
               end if

               if(swri%has_metric(SW_O)) then
                  if(otile_is_zero) then
                     otile(swx_in_batch,swy) = 0.0_DP
                  else
                     otile(swx_in_batch,swy) = cheb_int_product_piecewise_nD( &
                          swhome_all_coeffs(swy), &
                          swrem_batch_coeffs(swx_in_batch), &
                          nodes_template)
                  end if
               end if
            case (SCHEME_2DN_1DA)
               ! JCW: 2-D integration over a half-disc in an "integral frame"
               ! JCW: coordinate system, in which the z-axis is aligned with the
               ! JCW: vector between the two atoms in the atomblock. The 2-D
               ! JCW: half-disc is in a plane aligned parallel to the z-axis
               ! JCW: and the third coordinate is integrated analytically.

               ! JCW: 1-D analytic integral depends on m-values for SWOP x and
               ! JCW: SW y
               mx = swri%sw_idx_to_m(swx)
               my = swri%sw_idx_to_m(swy)

               if(swri%has_metric(SW_V)) then

                  call internal_eval_int_1Da(int_1Da,int_1Da_is_zero,mx,my)

                  if(int_1Da_is_zero) then
                     ! No need to do 2-D numerical integration
                     vtile(swx_in_batch,swy) = 0.0_DP
                  else
                     ! Do 2-D numerical integration and evaluate product
                     vtile(swx_in_batch,swy) = cheb_int_product_piecewise_nD( &
                          swhome_all_coeffs(swy), &
                          swpot_batch_coeffs(swx_in_batch), &
                          nodes_template) * int_1Da
                  end if
               end if

               ! @2DN_1DA_OVERLAP
               if (swri%int_scheme == SCHEME_2DN_1DA) then
                  call utils_assert(.not.swri%has_metric(SW_O),&
                       "Current implementation of 2Dn-1Da integration scheme &
                       &does not support the overlap metric, sorry.")
               end if

               !if(swri%has_metric(SW_O)) then

               !   if(otile_is_zero) then
               !      otile(swx_in_batch,swy) = 0.0_DP
               !   else
               !      otile(swx_in_batch,swy) = cheb_int_product_piecewise_nD( &
               !           swhome_all_coeffs(swy), &
               !           swrem_batch_coeffs(swx_in_batch), &
               !           nodes_template)
               !   end if
               !end if
            case default
               ! JCW: Unknown integration scheme
               call utils_abort('Error in '//myself//': &
                    &Unexpected integration scheme specification, ',&
                    swri%int_scheme)
         end select
         call timer_clock('swri_metric_mat_cheb_int_product',2)

      end do
!$OMP END PARALLEL DO

      call timer_clock(myself,2)
      call utils_trace_out(myself)

    end subroutine internal_tile_by_chebs

    subroutine internal_eval_int_1Da(int_1Da,is_zero,m1,m2)
      !========================================================================!
      ! Evaluate the 1-D analytic integral over the phi-dependent part of the  !
      ! integrand in the 2Dn-1Da scheme. This is the integral over the phi-    !
      ! dependent part of the product of two real spherical harmonics (RSHs),  !
      ! i.e.                                                                   !
      !   \int_{0}^{2\pi} \mathrm{d}\phi f_{m1}(\phi) g_{m2}(\phi)             !
      ! where f_{m}(\phi), g_{m}(\phi) are                                     !
      !   sin(|m|\phi) for m < 0                                               !
      !   cos(|m|\phi) for m > 0, and                                          !
      !   1            for m == 0.                                             !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   int_1Da (out) : result of integration                                !
      !   is_zero (out) : true if integral evaluates to zero, false otherwise  !
      !   m1, m2  (in)  : m indexes numbers of RSHs                            !
      !------------------------------------------------------------------------!
      ! Written by James C. Womack, 2018.                                      !
      !========================================================================!

      use constants, only: PI

      implicit none

      ! Arguments
      real(kind=DP), intent(out) :: int_1Da
      logical, intent(out)       :: is_zero
      integer, intent(in)        :: m1
      integer, intent(in)        :: m2

      if (m1 /= m2) then
         int_1Da = 0.0_DP
         is_zero = .true.
      else if (m1 == 0) then
         int_1Da = 2*PI
         is_zero = .false.
      else
         int_1Da = PI
         is_zero = .false.
      end if

    end subroutine internal_eval_int_1Da

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_request_tile(atomblock_to_do, xbatch_to_do, &
         workers_done, swri, atomblocklist, tilelist, n_tiles)
      !========================================================================!
      ! On non-root procs (workers): Requests, from master, a tile to process. !
      !                              Merges tiles during busy waits.           !
      ! On root proc (master): Dispatches subsequent tiles to workers, and to  !
      !                        itself, shows a progress indicator. Merges tiles!
      !                        once every time this is called. Keeps track of  !
      !                        which workers (including itself) completed.     !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   atomblock_to_do (out): The index of the tile's atomblock, received   !
      !                          from master is returned here. If there are    !
      !                          no more tiles, -1 is returned.                !
      !   xbatch_to_do (out): Ditto for the x-batch.                           !
      !   workers_done (in/out): Only relevant on master. Keeps track of which !
      !                          workers were given '-1' and so are done.      !
      !   swri (in/out): Needed for merging completed tiles.                   !
      !   atomblocklist (in/out): Needed for merging and keeping track of which!
      !                           atomblocks are completed.                    !
      !   tilelist (in/out): Only relevant on master. Keeps track of which     !
      !                      tiles need to be processed and which ones are     !
      !                      complete.                                         !
      !   n_tiles (in): The number of elements in tilelist.                    !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      !========================================================================!

      use comms, only: pub_root_proc_id, pub_my_proc_id, comms_send, &
           comms_irecv, comms_test
      use constants, only: NORMAL
      use rundat, only: pub_output_detail
      use utils, only: utils_int_to_str, utils_flushed_string_output

      implicit none

      ! jd: Arguments
      integer, intent(out)                                   :: atomblock_to_do
      integer, intent(out)                                   :: xbatch_to_do
      logical, intent(inout), allocatable                    :: workers_done(:)
      type(SW_RI), intent(inout)                             :: swri
      type(ATOMBLOCK_DESCRIPTOR), intent(inout), allocatable :: atomblocklist(:)
      type(TILE_DESCRIPTOR), intent(inout), allocatable      :: tilelist(:)
      integer, intent(in)                                    :: n_tiles

      ! jd: Local variables
      integer :: send_handle, recv_handle1, recv_handle2, recv_handle
      logical :: send_completed
      logical :: recv1_completed, recv2_completed
      integer :: tile
      integer :: n_tiles_completed
      integer, save :: n_tiles_completed_prev = -1 ! harmless state, only for
                                                   ! progress indicator
      character(len=*), parameter :: myself = 'internal_request_tile'

      ! -----------------------------------------------------------------------

      ! ----------------------------------------------
      ! --- Non-master workers ask master for tile ---
      ! ----------------------------------------------
      if(.not. pub_on_root) then
         ! jd: If a worker, request a tile. Deal with tile merging during comms.
         call timer_clock('swri_metric_mat_wait_req_worker',1)
         call comms_send(pub_root_proc_id, pub_my_proc_id, &
              length = 1, tag = TAG_TILE_REQUEST, &
              return_handle = send_handle, add_to_stack = .false.)
         send_completed = .false.
         ! jd: Deal with tile merging during the busy-wait
         do while (.not. send_completed)
            call internal_merge_attempt(swri, cell, atomblocklist)
            call comms_test(send_completed, send_handle)
         end do
         ! jd: Post receives for the tile
         call comms_irecv(pub_root_proc_id, atomblock_to_do, &
              length = 1, tag = TAG_ATOMBLOCK_ID, &
              handle = recv_handle1)
         call comms_irecv(pub_root_proc_id, xbatch_to_do, &
              length = 1, tag = TAG_XBATCH_ID, &
              handle = recv_handle2)
         ! jd: Deal with tile merging during the busy-wait
         recv1_completed = .false.
         recv2_completed = .false.
         do while (.not. (recv1_completed .and. recv2_completed))
            call internal_merge_attempt(swri, cell, atomblocklist)
            call comms_test(recv1_completed, recv_handle1)
            call comms_test(recv2_completed, recv_handle2)
         end do
         ! jd: Got the answer. Return.
         call timer_clock('swri_metric_mat_wait_req_worker',2)
      else
         ! ----------------------------------------------------
         ! --- Master dispatches tiles and merges its share ---
         ! ----------------------------------------------------
         ! jd: On master we deal with requests one by one. We exit only if there
         !     are no pending requests OR maximum pending sends reached.
         !     Theoretically there can be more pending sends than 2*num_workers
         !     at any one time (!) -- this would happen if a worker is very
         !     quick in dealing with its tile and manages to come back twice
         !     within a single do-while.
         call timer_clock('swri_metric_mat_wait_req_master',1)
         call internal_dispatch_tile(workers_done, tilelist, n_tiles)

         ! -------------------------------------------------------------
         ! --- Master is also a worker, but self-dispatch is simpler ---
         ! -------------------------------------------------------------
         ! jd: Master must deal a tile to itself too
         call internal_next_tile(tile,tilelist,n_tiles)
         xbatch_to_do = tilelist(tile)%xbatch
         atomblock_to_do = tilelist(tile)%atomblock
         tilelist(tile)%worker_proc = pub_root_proc_id
         ! jd: If tiles exhausted, returns -1, which is fine, the master-worker
         !     will know we're done..
         if(tile == -1) then
            workers_done(pub_root_proc_id) = .true.
         end if

         ! jd: Master must merge too.
         call internal_merge_attempt(swri, cell, atomblocklist)

         ! jd: Show a progress indicator. Towards the end of the run, when
         !     there are no more tiles to dispatch and master has finished its
         !     tile, it visits this subroutine much more often, as it's the
         !     only thing it has left to do. Hence we want to supress
         !     the output if nothing changed.
         n_tiles_completed = count(tilelist(1:)%worker_proc == TILE_COMPLETED)
         if(n_tiles_completed /= n_tiles_completed_prev) then
            if(pub_output_detail >= NORMAL) then
               call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
                    ']   Tiles completed: '//trim(utils_int_to_str(&
                    n_tiles_completed,'i8'))//&
                    ', left: '//trim(utils_int_to_str(&
                    count(tilelist(1:)%worker_proc == TILE_UNASSIGNED),'i8'))//&
                    ', in progress: '//trim(utils_int_to_str(&
                    count(tilelist(1:)%worker_proc >= 0),'i5'))//'.'//CRLF,.true.)
            end if
         end if

         n_tiles_completed_prev = n_tiles_completed
         call timer_clock('swri_metric_mat_wait_req_master',2)
      end if

    end subroutine internal_request_tile

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_dispatch_tile(workers_done, tilelist, n_tiles)
      !========================================================================!
      ! To be executed on the root proc only.                                  !
      ! Scans the tilelist for the first tile that has not been dispatched,    !
      ! and dispatches it -- ie. comms_sends the tile id to the proc that      !
      ! requested a tile. Also marks the previous tile that worker has been    !
      ! working on as completed.                                               !
      ! This is repeated for all pending request, or until their number exceeds!
      ! max_pending_sends. Then the master proc comms_waits for the sends to   !
      ! complete. Each send uses a different buffer, but the buffer can be     !
      ! modified outside this routine, that's why we wait.                     !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   workers_done (in/out): Keeps track which workers were given '-1' as  !
      !                          an indicator that there are no more tiles.    !
      !   tilelist (in/out): Keeps track of which tiles need to be processed   !
      !                      and which ones are complete.                      !
      !   n_tiles (in): The number of elements in tilelist.                    !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      !========================================================================!

      use comms, only: pub_total_num_procs, pub_any_source, comms_send, &
           comms_wait, comms_probe, comms_recv
      use timer, only: timer_clock
      use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

      implicit none

      ! jd: Arguments
      logical, intent(inout), allocatable                 :: workers_done(:)
      type(TILE_DESCRIPTOR), intent(inout), allocatable   :: tilelist(:)
      integer, intent(in)                                 :: n_tiles

      ! jd: Local variables
      integer :: worker
      integer :: tile
      logical :: request_pending
      integer :: pending_send, n_pending_sends, max_pending_sends
      integer, allocatable :: send_handles(:)
      character(len=*), parameter :: myself = 'internal_dispatch_tile'

      ! -----------------------------------------------------------------------

      call utils_assert(pub_on_root,myself//': Can only be called on root')

      call timer_clock('swri_metric_mat_dispatch',1)

      ! jd: Maximum number of pending sents (used to dimension handle arrays).
      !     Can be set to something larger (or smaller), but must be even.
      max_pending_sends = pub_total_num_procs * 2

      allocate(send_handles(max_pending_sends),stat=ierr)
      call utils_alloc_check(myself,'send_handles',ierr)

      request_pending = .true.
      n_pending_sends = 0
      ! -----------------------------------------------------------------------
      ! --- Loop until all requests are answered, or we exhaust handle buffers.
      ! -----------------------------------------------------------------------
      do while(request_pending .and. n_pending_sends < max_pending_sends)
         call comms_probe(request_pending, pub_any_source, TAG_TILE_REQUEST)
         if(request_pending) then
            call comms_recv(pub_any_source, worker, &
                 length = 1, tag = TAG_TILE_REQUEST)
            ! jd: Got a request, handle it.
            ! jd: Also, this means that whatever tile this worker has been
            !     working on until now is complete. Remember that.
            do tile = 1, n_tiles
               if(tilelist(tile)%worker_proc == worker) then
                  tilelist(tile)%worker_proc = TILE_COMPLETED
                  exit
               end if
            end do
            ! jd: Find the tile to dispatch.
            call internal_next_tile(tile,tilelist,n_tiles)
            ! jd: If tiles exhausted, returns -1, which is fine, the worker
            !     will know we're done.
            tilelist(tile)%worker_proc = worker
            call comms_send(worker, tilelist(tile)%atomblock, &
                 length = 1, tag = TAG_ATOMBLOCK_ID, &
                 return_handle = send_handles(n_pending_sends+1), &
                 add_to_stack = .false.)
            n_pending_sends = n_pending_sends + 1
            call utils_assert(n_pending_sends < max_pending_sends, myself//&
                 ': Too many pending sends or max_pending_sends not even', &
                 n_pending_sends)
            call comms_send(worker, tilelist(tile)%xbatch, &
                 length = 1, tag = TAG_XBATCH_ID, &
                 return_handle = send_handles(n_pending_sends+1), &
                 add_to_stack = .false.)
            n_pending_sends = n_pending_sends + 1
            ! jd: NB, no assert here, n_pending_sends examined in the do clause.
            ! jd: Keep track of which workers have been notified that there are
            !     no more tiles
            if(tile == -1) then
               workers_done(worker) = .true.
            end if
            ! jd: At this stage it's OK to deal with another request, since
            !     each request uses a different element of tilelist(:) as
            !     send buffer. We'll comms_wait() for all of these later.
         end if
      end do
      ! jd: There are no pending requests (or handle buffers full). However,
      !     we must wait for our dispatches (sends) to actually complete before
      !     exiting -- the send buffers could potentially be volatile outside
      !     this subroutine.
      do pending_send = 1, n_pending_sends
         call comms_wait(send_handles(pending_send))
      end do

      deallocate(send_handles,stat=ierr)
      call utils_dealloc_check(myself,'send_handles',ierr)

      call timer_clock('swri_metric_mat_dispatch',2)

    end subroutine internal_dispatch_tile

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_next_tile(tile_idx, tilelist, n_tiles)
      !========================================================================!
      ! To be executed on the root proc only.                                  !
      ! Scans the tilelist for the first tile that has not been assigned, and  !
      ! returns its index. If all tiles have been assigned, returns -1.        !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   tile_idx (out): Here the tile index is returned.                     !
      !   tilelist (in/out): Keeps track of which tiles need to be processed   !
      !                      and which ones are complete.                      !
      !   n_tiles (in): The number of elements in tilelist.                    !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      !========================================================================!

      use comms, only: pub_on_root
      use utils, only: utils_assert

      implicit none

      ! jd: Arguments
      integer, intent(out)                           :: tile_idx
      type(TILE_DESCRIPTOR), intent(in), allocatable :: tilelist(:)
      integer, intent(in)                            :: n_tiles

      ! -----------------------------------------------------------------------

      call utils_assert(pub_on_root, &
           myself//': The tile list is only defined on root')

      do tile_idx = 1, n_tiles
         if(tilelist(tile_idx)%worker_proc == TILE_UNASSIGNED) return
      end do

      ! jd: If we get here, there are no unassigned tiles. Return -1.
      tile_idx = -1

    end subroutine internal_next_tile

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_merge_attempt(swri, cell, atomblocklist)
      !========================================================================!
      ! Probes if there are any tiles waiting to be sent to this proc.         !
      ! If not, returns immediately. If yes, processes the tile and checks for !
      ! more. Processing a tile consists in:                                   !
      ! - Blocking-recv'ing it.                                                !
      ! - Storing it in the TILES hash table.                                  !
      ! - Checking if all tiles for an atomblock have been received.           !
      ! If so,                                                                 !
      ! - internal_merge_my_tiles() is called, where the tiles are merged into !
      !   an atomblock and evicted from the hash table. The atomblock is added !
      !   to the SPAM3 metric matrix structures.                               !
      ! - The atomblock is marked as completed in the atomblocklist.           !
      !                                                                        !
      ! This subroutine is called any time we busy-wait for comms, and is also !
      ! interleaved with computation to minimise wait times.                   !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri (in/out): The SW_RI within which we operate. It's modified, be- !
      !                  cause it hosts the hash table and the metric matrices.!
      !   cell (in): Simulation cell. Needed for PBCs.                         !
      !   atomblocklist (in/out): Needed for merging and keeping track of which!
      !                           atomblocks are completed.                    !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      !========================================================================!

      use comms, only: pub_my_proc_id, pub_any_source, comms_probe, comms_recv
      use hash_table, only: hash_table_add, hash_table_probe
      use simulation_cell, only: CELL_INFO
      use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(inout)                             :: swri
      type(CELL_INFO), intent(in)                            :: cell
      type(ATOMBLOCK_DESCRIPTOR), intent(inout), allocatable :: atomblocklist(:)

      ! jd: Local variables
      logical :: request_pending
      integer :: actual_source
      integer :: tile_length
      integer :: tile_id, atomblock_id, xbatch
      integer :: n_xbatches
      integer :: atoma, atomb, atoma_sparse_owner
      integer :: cur_xbatch
      integer :: V_or_O
      logical :: have_entire_atomblock
      real(kind=DP), allocatable :: recvtile(:)
      character(len=*), parameter :: myself = 'internal_merge_attempt'

      ! -----------------------------------------------------------------------

      ! jd: Do not put timer calls here, as there will be overheads.
      !     This subroutine is called millions of times and if there are no
      !     requests pending, it exits immediately.

      ! jd: See if there any tiles sent to us
      request_pending = .true.
      do while(request_pending)

         call comms_probe(request_pending, pub_any_source, TAG_TILE_ID, &
              status_source = actual_source)

         ! jd: Exit immediately if there are no requests pending, otherwise...
         if(request_pending) then

            n_xbatches = swri%n_xbatches
            tile_length = swri%xbatchsize * swri%quality%max_sws_per_centre

            allocate(recvtile(tile_length),stat=ierr)
            call utils_alloc_check(myself,'recvtile',ierr)

            ! jd: Blocking receive a pending tile. 'actual_source' ensures
            !     message overtaking between different senders does not harm us.
            call comms_recv(actual_source, tile_id, &
                 length = 1, tag = TAG_TILE_ID)
            call comms_recv(actual_source, recvtile, &
                 length = tile_length, tag = TAG_TILE_DATA)
            if(tile_id > 0) then
               V_or_O = SW_V
            else
               V_or_O = SW_O
               tile_id = -tile_id ! jd: Negative values signify O tiles
            end if
            ! jd: Decode atomblock id and xbatch from tile_id.
            atomblock_id = tile_id / n_xbatches + 1
            xbatch = mod(tile_id,n_xbatches)
            if(xbatch == 0) then
               xbatch = n_xbatches
               atomblock_id = atomblock_id - 1
            end if
            atoma = atomblocklist(atomblock_id)%atoma
            atomb = atomblocklist(atomblock_id)%atomb
            atoma_sparse_owner = atomblocklist(atomblock_id)%atoma_sparse_owner

            call utils_assert(atoma_sparse_owner /= -1, &
                 myself//': Received a tile for a block that''s completed', &
                 pub_my_proc_id, atomblock_id, tile_id)
            call utils_assert(atoma_sparse_owner == pub_my_proc_id, &
                 myself//': Received someone else''s tile or decode arithmetics&
                 & broken', atoma_sparse_owner, pub_my_proc_id)

            ! jd: Store tile in hash table.
            call hash_table_add(swri%tiles_ht, recvtile, tile_length, &
                 atoma, atomb, xbatch, V_or_O, overfill_strategy = 'F')

            ! jd: See if we now have all tiles for this atomblock.
            have_entire_atomblock = .true.
            do cur_xbatch = 1, n_xbatches
               if(-1 == hash_table_probe(swri%tiles_ht, &
                    atoma, atomb, cur_xbatch, V_or_O)) then
                  have_entire_atomblock = .false.
                  exit
               end if
            end do

            ! jd --- Merging ---
            if(have_entire_atomblock) then
               call internal_merge_my_tiles(swri,cell,V_or_O,atoma,atomb)
               ! jd: Mark block as done
               if(swri%has_metric(SW_V) .and. swri%has_metric(SW_O)) then
                  ! jd: Corner case, where we have both V and O matrices. Then
                  !     only mark as done once the latter is done.
                  if(V_or_O == SW_O) then
                     atomblocklist(atomblock_id)%atoma_sparse_owner = -1
                  end if
               else
                  ! jd: Default case
                  atomblocklist(atomblock_id)%atoma_sparse_owner = -1
               end if
            end if

            deallocate(recvtile,stat=ierr)
            call utils_dealloc_check(myself,'recvtile',ierr)

         end if ! have pending request
      end do ! jd: Loop until no more pending requests

    end subroutine internal_merge_attempt

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_merge_my_tiles(swri, cell, V_or_O, atoma, atomb)
      !========================================================================!
      ! A local operation in which all tiles of an atomblock are extracted from!
      ! the TILES hash table, merged into an atomblock, the remaining triangle !
      ! of the atomblock is filled by symmetry, and the atomblock is inserted  !
      ! into the SPAM3 metric matrix. This subroutine also stores the block on !
      ! disk, if swri%store_blocks_on_disk is T.                               !
      !                                                                        !
      ! If the 2Dn-1Da metric matrix evaluation scheme is in use, the          !
      ! assembled atomblock is for SWs in the integral coordinate frame. A     !
      ! real spherical harmonic (RSH) rotation matrix is constructed and       !
      ! applied to the atomblock to yield the atomblock in the ONETEP          !
      ! coordinate frame (i.e. with SWs aligned with Cartesian coordinate      !
      ! axes parallel to ONETEP's global Cartesian coordinate system).         !
      !                                                                        !
      ! The RSH rotation matrix is applied prior to insertion of the atomblock !
      ! into the SPAM3 metric matrix or storing the block on disk, so the      !
      ! result of this routine is the same regardless of the metric matrix     !
      ! evaluation scheme used.                                                !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri (in/out): The SW_RI within which we operate. It's modified, be- !
      !                  cause it hosts the hash table and the metric matrices.!
      !   cell (in): Simulation cell. Needed for PBCs.                         !
      !   V_or_O (in): Specifies whether we're dealing with tiles and blocks of!
      !                the electrostatic or overlap metric matrix.             !
      !   atoma, atomb (in): Define the atomblock. The proc on which this      !
      !                      subroutine is called is the sparse-owner of atoma.!
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      ! Modified by James C. Womack to support 2Dn-1Da metric matrix element   !
      ! evaluation scheme, 2018                                                !
      !========================================================================!

      use constants, only: garbage_real, SW_LETTERS
      use geometry, only: POINT, unit_vector,operator(-)
      use hash_table, only: hash_table_lookup, hash_table_list, &
           hash_table_remove
      use rundat, only: pub_hfx_bc_is_periodic
      use simulation_cell, only: CELL_INFO, minimum_image_displacement
      use sparse, only: sparse_put_block
      use sph_harm_rotation, only: &
           sph_harm_rot_get_active_Euler_ang_z2z, &
           sph_harm_rot_apply_RSH_Rmat_to_atomblock
      use timer, only: timer_clock
      use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
           utils_assert


      implicit none

      ! jd: Arguments
      type(SW_RI), intent(inout)  :: swri
      type(CELL_INFO), intent(in) :: cell
      integer, intent(in)         :: V_or_O
      integer, intent(in)         :: atoma, atomb

      ! jd: Local variables
      integer :: max_sws, num_sws
      integer :: n_xbatches, cur_xbatch
      integer :: tile_length, stored_tile_length
      integer :: xbeg, xend
      integer :: lvalx, mvalx, lval, mval, xsbidx, sbidx, swx, swy, species
      integer :: actual_xbatchsize
      real(kind=DP) :: dummy
      real(kind=DP) :: beta_ang, gamma_ang ! Euler rotation angles
      type(ATOM_CENTRE) :: centrea ! centre of atoma
      type(ATOM_CENTRE) :: centreb ! centre of atomb
      type(POINT) :: z_vec ! unit vector in z-direction of "integral" frame
      type(POINT) :: disp
      real(kind=DP), allocatable :: atomblock(:,:)
      real(kind=DP), allocatable :: tile_1D(:)
      real(kind=DP), allocatable :: tile(:,:)
      character(len=*), parameter :: myself = 'internal_merge_my_tiles'

      ! -----------------------------------------------------------------------

      max_sws = swri%quality%max_sws_per_centre
      num_sws = swri%quality%num_sws_per_centre
      n_xbatches = swri%n_xbatches
      tile_length = swri%xbatchsize * max_sws

      ! jd: Allocate workspace atomblock and tile
      allocate(atomblock(max_sws,max_sws),stat=ierr)
      call utils_alloc_check(myself,'atomblock',ierr)

      allocate(tile_1D(tile_length),stat=ierr)
      call utils_alloc_check(myself,'tile_1D',ierr)
      allocate(tile(swri%xbatchsize, max_sws),stat=ierr)
      call utils_alloc_check(myself,'tile',ierr)

      atomblock(:,:) = garbage_real ! Will help detect races or uninitialised
                                    ! elements, if any. Also, takes care of
                                    ! padding between num_sws and max_sws.

      ! jd: Transfer the tiles from the HT into the atomblock
      do cur_xbatch = 1, n_xbatches
         xbeg = (cur_xbatch-1) * swri%xbatchsize + 1
         xend = xbeg + swri%xbatchsize - 1
         if(xend > num_sws) then                 ! Handle last batch, which may
            xend = num_sws                       ! ..need to be truncated to
            actual_xbatchsize = xend - xbeg + 1  ! ..the size of the atomblock.
         else
            actual_xbatchsize = swri%xbatchsize
         end if
         call hash_table_lookup(tile_1D, stored_tile_length, &
              swri%tiles_ht, atoma, atomb, cur_xbatch, V_or_O)
         call utils_assert(tile_length == stored_tile_length, &
              myself//': tiles_ht misbehaves', tile_length, &
              stored_tile_length)
         tile = reshape(tile_1D,(/swri%xbatchsize,max_sws/))
         atomblock(xbeg:xend,1:num_sws) = tile(1:actual_xbatchsize,1:num_sws)
      end do

      ! jd: Forced consistency of V matrix and O matrix, i.e. filling of the
      !     second triangle by symmetry.
      call timer_clock('swri_metric_mat_consistency',1)
      if (force_metric_matrix_consistency) then
         swy = 0
         species = 1 ! @IDENTICAL_RADII
         do sbidx = 1, swri%max_num_bessels
            lval = swri%sphbessels(species,sbidx)%lval
            if (lval == -1) exit
            do mval=-lval,lval
               swy = swy + 1
               swx = 0
               do xsbidx = 1, swri%max_num_bessels
                  lvalx = swri%sphbessels(species,xsbidx)%lval
                  if (lvalx == -1) exit
                  do mvalx=-lvalx,lvalx
                     swx = swx + 1
                     if (swy > swx) then
                        atomblock(swy,swx) = &
                             (-1)**(lval+lvalx)*atomblock(swx,swy)
                     end if
                  end do
               end do
            end do
         end do
      end if
      call timer_clock('swri_metric_mat_consistency',2)

      ! -----------------------------------------------------------------------
      ! JCW: If swri%int_scheme is such that the atomblock has been
      ! JCW: evaluated in the "integral frame" (e.g. SCHEME_2DN_1DA)
      ! JCW: then the z-axis in for integrals in this atomblock is aligned
      ! JCW: along the vector between the two atoms associated with the
      ! JCW: atomblock.
      ! JCW:
      ! JCW: We must apply a rotation to the real spherical harmonic (RSH)
      ! JCW: component of the SWs in the matrix elements to express our metric
      ! JCW: matrix in terms of SWs in the "ONETEP frame" (with z-direction
      ! JCW: of RSHs along the z-direction of the Cartesian coordinate system
      ! JCW: of the simulation cell).
      ! -----------------------------------------------------------------------
      select case (swri%int_scheme)
        case (SCHEME_3DC)
           ! ------------------------------------------------------------------
           ! JCW: Integration in a 3D sphere
           ! ------------------------------------------------------------------
           ! Nothing to do here, since the metric matrix elements are evaluated
           ! in the ONETEP frame
           continue

        case (SCHEME_2DN_1DA)
           ! ------------------------------------------------------------------
           ! JCW: Integration on a 2D half-disc with 1D analytic integration
           ! ------------------------------------------------------------------
           ! 0. Determine the number of spherical Bessels for each l-value.
           !    This is not always swri%quality%max_q, as certain q,l
           !    combinations are excluded because the q value corresponds to
           !    a spherical Bessel kinetic energy greater than the psinc
           !    KE cutoff
           !    [ This is computed in swri_num_sph_functions and the result is
           !    placed in swri%num_q_per_l when swri_init_container_stage_1
           !    is called ]
           !
           call timer_clock('swri_metric_mat_RSH_rotation',1)
           ! 1. Determine Euler rotation angles for active rotation about fixed
           !    coordinate axes which brings a vector aligned along the z-axis
           !    of the integral frame into alignment with the z-axis of the
           !    Cartesian coordinate system of the simulation cell (ONETEP
           !    frame, providing the fixed coordinate axes of rotation),
           !    i.e. 'z to z' rotation)
           !    * Get z-direction basis vector of integral frame expressed in
           !      ONETEP frame (fixed coord system of rotation)
           call swri_init_centre(centrea, par, elements, atoma)
           call swri_init_centre(centreb, par, elements, atomb)
           !      @optimize: swri_init_centre sets values we don't need
           !                 (only need position in cell of atoma, atomb)

           if(.not. any(pub_hfx_bc_is_periodic)) then
              ! jd: OBC
              disp = centreb%incell - centrea%incell
           else if(all(pub_hfx_bc_is_periodic)) then
              ! jd: PBC (@?)
              call minimum_image_displacement(dummy, disp, &
                   centreb%incell, centrea%incell, cell)
           else
              call utils_abort(myself//': Mixed BCs are not currently supported.')
           end if

           z_vec = unit_vector(disp)

           !      atoma and atomb are from an ATOMBLOCK_DESCRIPTOR, as in
           !      internal_tile_by_chebs, so we can assume the same conventions
           !      for which atom refers to home and non-home sphere SWOP:
           !      - centreb is swx, i.e. the SWOP in non-home sphere
           !      - centrea is swy, i.e. the SWOP in home sphere

           !    * Get active Euler rotation angles for z-to-z rotation
           call sph_harm_rot_get_active_Euler_ang_z2z(beta_ang,gamma_ang,z_vec)

           ! 2. Apply RSH rotation to integral frame metric matrix atomblock
           !    (i.e. R* M R) to yield metric matrix atomblock in ONETEP frame
           !    Perform rotation of integral coordinate frame atomblock into
           !    ONETEP coordinate frame
           call sph_harm_rot_apply_RSH_Rmat_to_atomblock(atomblock,&
                swri%quality%max_l,swri%num_q_per_l(1,0:swri%quality%max_l),&
                0.0_DP,beta_ang,gamma_ang)
           ! @IDENTICAL_RADII is assumed here, since the number of q-values per
           ! l-value (and in general, the number of SWs per atomic centre) is
           ! related to the radius of the SWs. Here we assume that the number of
           ! q-values per l-value is identical for all species, but this is only
           ! true if species have identical radii (see swri_num_sph_functions).
           ! If the size of the SW basis did vary per-species, then the
           ! sph_harm_rot_apply_RSH_Rmat_to_atomblock routine would need further
           ! modification to support non-square atomblocks and to accept the
           ! full swri%num_q_per_l array, rather than simply the first entry.
           call timer_clock('swri_metric_mat_RSH_rotation',2)
        case default
           call utils_abort('Error in '//myself//': &
                &Unexpected integration scheme specification, ',&
                swri%int_scheme)
      end select

      ! jd: Store block on disk, unless we've just erased blocks and don't want
      !     to regenerate them.
      if(swri%store_blocks_on_disk) then
         call internal_write_block(SW_LETTERS(V_or_O), &
              atomblock, atoma, atomb)
      end if

      ! jd: Store into the full metric matrix
      call sparse_put_block(atomblock, &
           swri%full_metric_matrices(V_OR_O), atomb, atoma)

      ! jd: Deallocate workspace atomblock and tile
      deallocate(tile,stat=ierr)
      call utils_dealloc_check(myself,'tile',ierr)
      deallocate(tile_1D,stat=ierr)
      call utils_dealloc_check(myself,'tile_1D',ierr)
      deallocate(atomblock,stat=ierr)
      call utils_dealloc_check(myself,'atomblock',ierr)

      ! jd: Evict the tiles from the hash table
      call hash_table_list(swri%tiles_ht,0)
      do cur_xbatch = 1, n_xbatches
         call hash_table_remove(swri%tiles_ht, atoma, atomb, cur_xbatch, V_or_O)
      end do

    end subroutine internal_merge_my_tiles

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_gen_tile_list(swri, atomblocklist, tilelist, n_tiles)
      !========================================================================!
      ! To be executed on the root proc only.                                  !
      ! Allocates and fills the tilelist, which is a structure for keeping     !
      ! track of which tiles have been assigned to which procs, and which ones !
      ! are completed.                                                         !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri (in): The SW_RI in which we are operating.                      !
      !   atomblocklist (in): The list of atomblocks to be divided into tiles. !
      !   tilelist (out): This is allocated here and contains the tile list    !
      !                   on exit.                                             !
      !   n_tiles (out): The number of elements in tilelist, set on exit.      !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      !========================================================================!

      use constants, only: garbage_int
      use utils, only: utils_alloc_check, utils_assert, utils_int_to_str, &
           utils_flushed_string_output

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(in)                             :: swri
      type(ATOMBLOCK_DESCRIPTOR), intent(in), allocatable :: atomblocklist(:)
      type(TILE_DESCRIPTOR), intent(out), allocatable     :: tilelist(:)
      integer, intent(out)                                :: n_tiles

      ! jd: Local variables
      integer :: num_sws
      integer :: atomblock, n_atomblocks
      integer :: xbatch, n_xbatches
      integer :: tile
      integer :: ierr
      character(len=*), parameter :: myself = 'internal_gen_tile_list'

      ! -----------------------------------------------------------------------

      num_sws = swri%quality%num_sws_per_centre
      n_xbatches = swri%n_xbatches
      n_atomblocks = size(atomblocklist)
      n_tiles = n_atomblocks * n_xbatches

      call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
           ']   '//trim(utils_int_to_str(n_tiles))//' tiles to process ('//&
           trim(utils_int_to_str(n_atomblocks))//' atomblocks x '//&
           trim(utils_int_to_str(n_xbatches))//' batches)'//CRLF,.true.)

      allocate(tilelist(-1:n_tiles), stat=ierr)
      call utils_alloc_check(myself,'tilelist',ierr)

      ! jd: Partition atomblocks into tiles.
      xbatch = 1
      atomblock = 1
      do tile = 1, n_tiles

         tilelist(tile)%atomblock = atomblock
         tilelist(tile)%xbatch = xbatch
         tilelist(tile)%worker_proc = TILE_UNASSIGNED

         xbatch = xbatch + 1
         if(xbatch > n_xbatches) then
            xbatch = 1
            atomblock = atomblock + 1
         end if
      end do

      call utils_assert(atomblock == n_atomblocks + 1, &
           myself//': Atomblock count mismatch')
      call utils_assert(xbatch == 1, &
           myself//': Xbatch count mismatch')

      ! jd: Special tile -1 is used to signal completion
      tilelist(-1)%atomblock = -1
      tilelist(-1)%xbatch = -1
      tilelist(-1)%worker_proc = TILE_COMPLETED

      ! jd: Tile 0 is unused, will help debugging if something goes south
      tilelist(0)%atomblock = garbage_int
      tilelist(0)%xbatch = garbage_int
      tilelist(0)%worker_proc = garbage_int

    end subroutine internal_gen_tile_list

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_calc_home_sws(swhome_all_coeffs, & ! out
         nodes_template, swri)                                 ! in
      !========================================================================!
      ! Performs a Chebyshev expansion, in the home sphere, of all SWs.        !
      ! The home sphere is a sphere with a radius of the common NGWF radius,   !
      ! centred on the origin of the SWs. This expansion is thus valid for     !
      ! *any* atom in the system, provided the current approach of identical   !
      ! NGWF radii is retained.                                                !
      !                                                                        !
      ! This calculation takes only a short while, compared at least Nat expan-!
      ! sions that will be needed later for SWOPs. Thus it is only OMP paralle-!
      ! lised, and is replicated across all MPI ranks.                         !
      !                                                                        !
      ! The functions to be expanded depend on the metric matrix scheme used:  !
      !                                                                        !
      ! 3Dc (swri%int_scheme == SCHEME_3DC)                                    !
      ! ===================================                                    !
      ! Full SW functions evaluated at Chebyshev nodes in 3-D Cartesian        !
      ! coordinate system within the home sphere.                              !
      !                                                                        !
      ! 2Dn-1Da (swri%int_scheme == SCHEME_2DN_1DA)                            !
      ! ===========================================                            !
      ! (r,theta)-dependent component of SW functions evaluated at Chebyshev   !
      ! notes on a half-disc in the Cartesian yz plane. These are multiplied   !
      ! by Jacobian determinants associated with the coordinate                !
      ! transformations necessary to express the metric matrix elements in a   !
      ! form separable into 2-D and 1-D components AND evaluate the 2-D        !
      ! integral in Cartesian coordinates (for details, see documentation for  !
      ! swri_eval_swrt_at_nodes_on_halfdisc_fast).                             !
      !------------------------------------------------------------------------!
      ! *** Important considerations (3Dc scheme): ***                         !
      ! Keeping all SWs in the Chebyshev representation requires a lot of RAM. !
      ! However, since this expansion is used ubiquitously, the alternatives   !
      ! are worse. The alternatives are:                                       !
      ! - Re-calculating them as needed: takes significant time, becoming a    !
      !   bottleneck by a factor of about 3 (tried that).                      !
      ! - Having each proc own a batch of home SWs ('y-batching'). This is     !
      !   also slower due to comms times, even when comms can be initiated     !
      !   instantly. In one tested case evaluation of one SWy cost 25ms, and   !
      !   MPI transmission of the same thing cost 11ms.                        !
      ! - Batching over the home SWs, ie. having all procs calculate the same  !
      !   batch of SWs, processing them, and moving to the next batch. This    !
      !   has the huge disadvantage of then having to recalculate SWOPs (or    !
      !   remote SWs, in the case of overlap metric) each time. Also, the      !
      !   Chebyshev product routines benefit from using large batches. This    !
      !   has not been tried, but looks like the worst alternative.            !
      !                                                                        !
      ! The leading assumption is that the memory strain can be relieved by    !
      ! using OpenMP and having 4-6 MPI ranks per node, or, in the worst-case  !
      ! scenario, by undersubscribing nodes.                                   !
      !                                                                        !
      ! NOTE: For the newer 2Dn-1Da metric matrix evaluation scheme, memory    !
      ! costs are significantly reduced, since we are now storing SWRTs        !
      ! (r,theta)-component of SWs in 2-D, rather than 3-D as above. As a      !
      ! consequence, the above memory concerns are likely no-longer relevant   !
      ! when using 2Dn-1Da to evaluate the metric matrix.                      !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swhome_all_coeffs (in/out):     Here the Chebyshev expansion coeffi- !
      !                                   cients are returned. Intent is in/out!
      !                                   because this structure needs to be   !
      !                                   allocated with cheb_alloc_coeff() in !
      !                                   advance.                             !
      !   nodes_template (in): The Chebyshev nodes where the expansion         !
      !                        takes place.                                    !
      !   swri (in): The SW_RI within which we operate.                        !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04 based on much older subroutines.  !
      ! Modified by James C. Womack to support 2Dn-1Da metric matrix element   !
      ! evaluation scheme, 2018                                                !
      !========================================================================!

      use chebyshev_rep, only: CHEB_NODES, CHEB_COEFFS, &
           cheb_expansion_piecewise_nD
      use comms, only: comms_barrier
      use rundat, only: pub_threads_max
      use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check
      use timer, only: timer_clock

      implicit none

      ! jd: Arguments
      type(CHEB_COEFFS), intent(inout) :: swhome_all_coeffs(:) ! sw_idx
      type(CHEB_NODES), intent(in)     :: nodes_template
      type(SW_RI), intent(in)          :: swri

      ! jd: Local variables
      real(kind=DP), allocatable :: swy_values(:,:) ! 2nd index is tid
      integer :: ly, my, bessely, swy
      real(kind=DP) :: ay, qy
      integer :: tid
      integer :: ierr
      integer :: species
!$    integer, external :: omp_get_thread_num
      character(len=*), parameter :: myself = 'internal_calc_home_sws'

      ! -----------------------------------------------------------------------

      call timer_clock('swri_metric_mat_home_sws',1)

      allocate(swy_values(nodes_template%n_points_total,pub_threads_max), &
           stat=ierr)
      call utils_alloc_check(myself,'swy_values',ierr)

      ! jd: Loop over home SWs (y)
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(swy, bessely, my, ly, qy, ay, tid) &
!$OMP SHARED(nodes_template, swy_values, swri, species, &
!$OMP      swhome_all_coeffs, pub_threads_max)
      do swy = 1, swri%quality%num_sws_per_centre

         bessely = swri%sw_idx_to_bessel_idx(swy)
         my = swri%sw_idx_to_m(swy)

         species=1 ! @IDENTICAL_RADII
         ly = swri%sphbessels(species,bessely)%lval
         qy = swri%sphbessels(species,bessely)%qval
         ay = swri%sphbessels(species,bessely)%aval

         tid = 0
!$       tid = omp_get_thread_num()

         ! --------------------------------------------------------------------
         ! JCW: Select method for computing function at nodes based on scheme,
         ! JCW: since the function evaluated at the nodes depends upon this.
         ! JCW:
         ! JCW: 3Dc scheme:     3 (sphere)
         ! JCW: 2Dn-1Da scheme: 2 (half-disc)
         ! --------------------------------------------------------------------
         select case (swri%int_scheme)
            case (SCHEME_3DC)
               ! --------------------------------------------------------------
               ! JCW: 3-D Chebyshev: evaluate full spherical wave at nodes
               ! --------------------------------------------------------------
               ! jd: Eval home SW y
               call swri_eval_sw_at_nodes_fast(swy_values(:,tid+1), &
                    nodes_template, ly, my, qy)
            case (SCHEME_2DN_1DA)
               ! --------------------------------------------------------------
               ! JCW: 2-D numerical, 1-D analytic: evaluate component of
               ! JCW: spherical wave dependent on (r, theta)
               ! JCW: * Azimuthal coordinate will be integrated analytically
               ! JCW:   and separate to the 2-D numerical integration on a
               ! JCW:   half-disc
               ! JCW: * "swrt" refers to the (r,theta)-dependent part of the
               ! JCW:   full SW
               ! --------------------------------------------------------------
                 call swri_eval_swrt_at_nodes_on_halfdisc_fast(swy_values(:,tid+1), &
                      nodes_template, ly, my, qy)
              case default
                 call utils_abort('Error in '//myself//': &
                      &Unexpected integration scheme specification, ',&
                      swri%int_scheme)
         end select

         ! jd: Perform Chebyshev expansion for home SW y
         call cheb_expansion_piecewise_nD(swhome_all_coeffs(swy), &
              nodes_template, swy_values(:,tid+1))

      end do ! home SWs
!$OMP END PARALLEL DO

      deallocate(swy_values,stat=ierr)
      call utils_dealloc_check(myself,'swy_values',ierr)

      call comms_barrier
      call timer_clock('swri_metric_mat_home_sws',2)

    end subroutine internal_calc_home_sws

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_gen_atomblock_list(swri, elements, cell, atomblocklist,&
        read_blocks_from_disk, erase_blocks_from_disk)
      !========================================================================!
      ! Generates a list of metric matrix atomblocks that need to be processed.!
      ! 'Processed' almost always means 'computed', except for the 'RD'        !
      ! scenario, where the metric matrix is read from disk, and the atomblock !
      ! list is only used to determine the blocks to disassembly.              !
      !
      ! The atomblock list is global (contains blocks between atoms on all     !
      ! procs) and uses after-space-filling-curve numbering for atoms. It is   !
      ! constructed as follows. For every atom A local to this proc all        !
      ! S-neighbours are determined (atoms B). If neither belongs to the SWRI, !
      ! the pair is discarded. This is because we need atomblocks between all  !
      ! atoms within the SWRI *and* those where at least one is within. All    !
      ! pairs where orig_atomb > orig_atoma are discarded. This is because we  !
      ! only need to fill one triangle of the metric matrix due to symmetry.   !
      ! The blocks on the diagonal are not needed, since they are found analy- !
      ! tically. Note the use of orig numbering in this one case. This is be-  !
      ! cause we don't want the triangle definitions to change depending on    !
      ! the SFC -- we'd like atomblocks on disk to be transferable to calcula- !
      ! tions on subsets of atoms. Unless 'R' is used, we then make an attempt !
      ! to read each atomblock from a file. If this succeeds, the atomblock is !
      ! *removed* from the list (as it does not need to be computed). If 'E'   !
      ! has been specified, an attempt is made to erase the block at this      !
      ! stage. Up till now each proc has handled their own share of blocks, ie.!
      ! those, whose atoma it sparse-owns. Subsequently the results are        !
      ! gathered across procs, so the atomblocklist is replicated. Finally,    !
      ! the atomblocks are heapsorted by proximity to a user-specified point.  !
      ! The larger of the two distances (atoma-refpt, atomb-refpt) is used as  !
      ! a criterion.                                                           !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri (in/out): The SW_RI within which we operate. Intent is in/out,  !
      !                  because we put blocks into the metric matrices if     !
      !                  they can be read from disk.                           !
      !   elements (in): Needed to determine whether an atom belongs to SW_RI. !
      !   cell (in): Simulation cell. Needed for PBCs.                         !
      !   atomblocklist (out): Is allocated here and results are returned here.!
      !   read_blocks_from_disk (in): If T, an attempt is made to read each    !
      !                               needed atomblock from disk and insert it !
      !                               into the metric matrix.                  !
      !   erase_blocks_from_disk (in): If T, an attempt is made to erase each  !
      !                                needed atomblock from disk.             !
      !                                                                        !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04 based on much older subroutines.  !
      !========================================================================!

      use comms, only: pub_my_proc_id, comms_allgather, comms_reduce
      use geometry, only: geometry_DISTANCE, POINT
      use ion, only: ELEMENT
      use rundat, only: pub_devel_code, pub_swri_proximity_sort_point, &
           pub_hfx_bc_is_periodic
      use simulation_cell, only: CELL_INFO, minimum_image_distance
      use sparse, only: sparse_put_block
      use utils, only: utils_devel_code, utils_int_to_str, &
           utils_alloc_check, utils_dealloc_check, utils_flushed_string_output,&
           utils_point_to_str, utils_heapsort

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(inout)  :: swri
      type(ATOMBLOCK_DESCRIPTOR), intent(out), allocatable :: atomblocklist(:)
      type(ELEMENT), intent(in)   :: elements(par%nat)
      type(CELL_INFO), intent(in) :: cell
      logical, intent(in)         :: read_blocks_from_disk
      logical, intent(in)         :: erase_blocks_from_disk

      ! jd: Internal variables
      integer :: local_a, atoma, orig_atoma ! atom indices
      integer :: b_idx, atomb, orig_atomb   ! atom indices
      integer :: max_n_atomblocks
      integer :: n_atomblocks, n_atomblocks_read, n_my_atomblocks
      integer :: i_blk
      integer :: error
      logical :: debug_show_skipped
      logical :: block_must_be_calculated
      type(POINT) :: refpt
      real(kind=DP) :: d1, d2
      integer, allocatable :: my_atomblocks(:,:)
      integer, allocatable :: all_atomblocks(:,:)
      real(kind=DP), allocatable :: distances(:)
      integer, allocatable :: distances_key(:)
      type(ATOMBLOCK_DESCRIPTOR), allocatable :: unsorted_atomblocklist(:)
      real(kind=DP), allocatable :: atomblock(:,:)
      character(len=*), parameter :: myself='internal_gen_atomblock_list'

      ! -----------------------------------------------------------------------

      call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
           '] - Generating and communicating atomblock list.'//CRLF)

      debug_show_skipped = utils_devel_code(.false.,'SW','SHOW_SKIPPED', &
           pub_devel_code, no_bcast=.true., no_warn=.true.)

      ! jd: Absolute maximum # of atomblocks (fully dense scenario). Used only
      !     to dimension arrays of handles or indices, so no worry about O(N^2).
      max_n_atomblocks = par%nat * (par%nat-1) / 2

      ! jd: Go over my ('sparse-my') blocks, remember which they are and who
      !     they belong to.
      allocate(all_atomblocks(3,max_n_atomblocks),stat=ierr)
      call utils_alloc_check(myself,'all_atomblocks',ierr)
      allocate(my_atomblocks(3,max_n_atomblocks),stat=ierr)
      call utils_alloc_check(myself,'my_atomblocks',ierr)

      ! jd: Allocate workspace atomblock
      allocate(atomblock(&
           swri%quality%max_sws_per_centre,swri%quality%max_sws_per_centre), &
           stat=ierr)
      call utils_alloc_check(myself,'atomblock',ierr)

      ! jd: First prepare portion for atoms A local to this proc
      n_atomblocks = 0
      n_atomblocks_read = 0
      loop_A1:                                                                &
      do local_a=1, par%num_atoms_on_proc(pub_my_proc_id)
         atoma = par%first_atom_on_proc(pub_my_proc_id) + local_a -1
         orig_atoma = par%orig_atom(atoma)

         ! jd: Loop over atoms B that are S-neighbours with A
         loop_B1:                                                             &
         do b_idx = swri%s_atoms_nl%first_idx(atoma), &
              swri%s_atoms_nl%last_idx(atoma)
            atomb = swri%s_atoms_nl%neighbours(b_idx)
            orig_atomb = par%orig_atom(atomb)

            ! jd: Only pick one triangle, the rest will be filled by symmetry.
            !     Use original numbering.
            if (orig_atoma .lt. orig_atomb) then

               ! jd: Eliminate blocks that do not involve at least one atom
               !     that is inside our SWRI
               if(.not. elements(orig_atoma)%in_swri(swri_handle) .and. &
                    .not. elements(orig_atomb)%in_swri(swri_handle)) then
                  if(debug_show_skipped) then
                     write(stdout,'(a,i0,a,i0,a,a,a,i0,a,i0,a,a,a)') &
                     'Skipping atomblock (off-diag) orig/SFC: ', orig_atoma, &
                     '/', atoma, ' (', trim(elements(orig_atoma)%species_id), &
                     ') : ', orig_atomb, '/', atomb, ' (', &
                     trim(elements(orig_atomb)%species_id),') '
                  end if
                  cycle
               end if

               block_must_be_calculated = .true.

               if(read_blocks_from_disk) then
                  ! jd: This block is needed, check if it can be read
                  block_must_be_calculated = .false.
                  if(swri%has_metric(SW_V)) then
                     call internal_read_block('V',atomblock,atoma,atomb,error)
                     if(error == 0) then
                        call sparse_put_block(atomblock, &
                             swri%full_metric_matrices(SW_V), atomb, atoma)
                     else
                        block_must_be_calculated = .true.
                     end if
                  end if
                  if(swri%has_metric(SW_O)) then
                     call internal_read_block('O',atomblock,atoma,atomb,error)
                     if(error == 0) then
                        call sparse_put_block(atomblock, &
                             swri%full_metric_matrices(SW_O), atomb, atoma)
                     else
                        block_must_be_calculated = .true.
                     end if
                  end if
               end if

               ! jd: Erase blocks from disk, if user asked to. Each proc erases
               !     the files whose atoma it owns.
               if(erase_blocks_from_disk) then
                  call internal_erase_atomblock_file(swri, atoma, atomb)
               end if

               ! jd: The block is needed and can't be read, add it to the list
               !     of blocks that must be calculated.
               if(block_must_be_calculated) then
                  n_atomblocks = n_atomblocks + 1
                  my_atomblocks(1,n_atomblocks) = atoma
                  my_atomblocks(2,n_atomblocks) = atomb
                  my_atomblocks(3,n_atomblocks) = pub_my_proc_id
               else
                  n_atomblocks_read = n_atomblocks_read + 1
               end if

            end if ! atomblock is needed

         end do loop_B1
      end do loop_A1

      n_my_atomblocks = n_atomblocks

      call comms_reduce('SUM',n_atomblocks)
      call comms_reduce('SUM',n_atomblocks_read)

      ! jd: Allgather results, now everyone knows what atomblocks need to
      !     be taken care of. If there are no blocks (system with 1 atom,
      !     or otherwise no overlapping atoms), do not do comms, as all_atomblocks
      !     and my_atomblocks have one of the extents == 0, meaning loc() returns
      !     the same address for both, and aliasing check in allgather complains.
      if(n_atomblocks /= 0) then
         call comms_allgather(all_atomblocks, my_atomblocks, 3*n_my_atomblocks, &
              (/-1/), (/-1/))
      end if

      ! jd: Convert all_atomblocks() to atomblocklist, now that we know the size
      allocate(atomblocklist(n_atomblocks),stat=ierr)
      call utils_alloc_check(myself,'atomblocklist',ierr) ! < dealloc in parent
      allocate(unsorted_atomblocklist(n_atomblocks),stat=ierr)
      call utils_alloc_check(myself,'unsorted_atomblocklist',ierr)

      allocate(distances(n_atomblocks),stat=ierr)
      call utils_alloc_check(myself,'distances',ierr)
      allocate(distances_key(n_atomblocks),stat=ierr)
      call utils_alloc_check(myself,'distances_key',ierr)

      ! jd: Sort blocks by the larger of the two distances from refpt
      refpt = POINT(pub_swri_proximity_sort_point(1), &
           pub_swri_proximity_sort_point(2), pub_swri_proximity_sort_point(3))
      do i_blk = 1, n_atomblocks
         unsorted_atomblocklist(i_blk)%atoma = all_atomblocks(1,i_blk)
         unsorted_atomblocklist(i_blk)%atomb = all_atomblocks(2,i_blk)
         unsorted_atomblocklist(i_blk)%orig_atoma = &
              par%orig_atom(unsorted_atomblocklist(i_blk)%atoma)
         unsorted_atomblocklist(i_blk)%orig_atomb = &
              par%orig_atom(unsorted_atomblocklist(i_blk)%atomb)
         unsorted_atomblocklist(i_blk)%atoma_sparse_owner = &
              all_atomblocks(3,i_blk)

         if(.not. any(pub_hfx_bc_is_periodic)) then
            ! jd: OBC
            d1 = geometry_DISTANCE(elements(unsorted_atomblocklist(i_blk)&
                 %orig_atoma)%centre, refpt)
            d2 = geometry_DISTANCE(elements(unsorted_atomblocklist(i_blk)&
                 %orig_atomb)%centre, refpt)
         else if(all(pub_hfx_bc_is_periodic)) then
            ! jd: PBC
            d1 = minimum_image_distance(elements(unsorted_atomblocklist(i_blk)&
                 %orig_atoma)%centre, refpt, cell)
            d2 = minimum_image_distance(elements(unsorted_atomblocklist(i_blk)&
                 %orig_atomb)%centre, refpt, cell)
         else
            call utils_abort(myself//': Mixed BCs are not currently supported.')
         end if

         distances(i_blk) = max(d1,d2)
      end do

      ! jd: Clean up
      deallocate(my_atomblocks,stat=ierr)
      call utils_dealloc_check(myself,'my_atomblocks',ierr)
      deallocate(all_atomblocks,stat=ierr)
      call utils_dealloc_check(myself,'all_atomblocks',ierr)
      deallocate(atomblock,stat=ierr)
      call utils_dealloc_check(myself,'atomblock',ierr)

      if(read_blocks_from_disk) then
         call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
              ']   '//trim(utils_int_to_str(n_atomblocks_read))//&
              ' atomblocks read from disk.'//CRLF)
      else
         call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
              ']   Atomblocks not read from disk.'//CRLF)
      end if
      call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
           ']   '//trim(utils_int_to_str(n_atomblocks))&
           //' atomblocks to process.'//CRLF)

      call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
           ']   Sorting atomblocks according to proximitiy to ('//&
           trim(utils_point_to_str(refpt,'f12.6'))//').'//CRLF)

      ! jd: Sort atomblocks based on proximity to user-defined point
      call utils_heapsort(n_atomblocks, distances, distances_key)

      ! jd: Prepare the sorted version of atomblocklist
      do i_blk = 1, n_atomblocks
         atomblocklist(i_blk) = unsorted_atomblocklist(distances_key(i_blk))
      end do

      deallocate(distances_key,stat=ierr)
      call utils_dealloc_check(myself,'distances_key',ierr)
      deallocate(distances,stat=ierr)
      call utils_dealloc_check(myself,'distances',ierr)
      deallocate(unsorted_atomblocklist,stat=ierr)
      call utils_dealloc_check(myself,'unsorted_atomblocklist',ierr)

    end subroutine internal_gen_atomblock_list

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_read_matrix(matrix, swriname, matrixname)
      !========================================================================!
      ! Reads a SPAM3 matrix from a file. Aborts if the file cannot be opened. !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   matrix (inout):  The SPAM3 matrix to read into.                      !
      !                    The argument is (inout), because the matrix needs   !
      !                    have been initialized first (sparsity pattern, etc.)!
      !   swriname (in):   The name of the swri to display in the message.     !
      !   matrixname (in): The filename suffix.                                !
      !========================================================================!

      use comms, only: pub_on_root
      use rundat, only: pub_rootname
      use utils, only: utils_abort, utils_flushed_string_output
      use sparse, only: sparse_read

      implicit none

      ! Arguments
      type(SPAM3), intent(inout)     :: matrix
      character(len=*), intent(in)   :: swriname
      character(len=*), intent(in)   :: matrixname

      ! Local variables
      character(len=512) :: filename
      logical :: fileexists
      character(len=*), parameter :: myself = 'internal_read_matrix'

      ! --------------------------------------------------------------------

      write(filename,'(2a)') trim(pub_rootname),'.'//trim(matrixname)
#ifdef HDF5
      filename=trim(filename)//'.h5'
#endif
      call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
           '] - Reading '//trim(matrixname)//' from file "'//&
           trim(filename)//'"...')

      ! Check that the file exists
      ! ndmh: only root proc needs to be able to see the file
      if (pub_on_root) then
         inquire(file=filename,exist=fileexists)
      else
         fileexists = .true.
      end if

      if (fileexists) then
         ! Read the matrix from this file
         call sparse_read(matrix,trim(filename))
      else
         call utils_abort('File "'//trim(filename)//'" not found.')
      end if

      if (pub_on_root) write(stdout,'(a)') 'done'

    end subroutine internal_read_matrix

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_write_matrix(matrix, swriname, matrixname)
      !========================================================================!
      ! Writes a SPAM3 matrix to a file.                                       !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   matrix (in):     The SPAM3 matrix to write out.                      !
      !   swriname (in):   The name of the swri to display in the message.     !
      !   matrixname (in): The filename suffix.                                !
      !========================================================================!

      use comms, only: pub_on_root
      use rundat, only: pub_rootname
      use sparse, only: sparse_write
      use utils, only: utils_flushed_string_output

      implicit none

      ! Arguments
      type(SPAM3), intent(in)        :: matrix
      character(len=*), intent(in)   :: swriname
      character(len=*), intent(in)   :: matrixname

      ! Local variables
      character(len=512) :: filename
      character(len=*), parameter :: myself = 'internal_write_matrix'

      ! --------------------------------------------------------------------

      write(filename,'(2a)') trim(pub_rootname),'.'//trim(matrixname)
#ifdef HDF5
      filename=trim(filename)//'.h5'
#endif
      call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
           '] - Writing '//trim(matrixname)//' to file "'//&
           trim(filename)//'"...')

      call sparse_write(matrix,trim(filename))

      if (pub_on_root) then
         write(stdout,'(a)') 'done'
      end if

    end subroutine internal_write_matrix

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_write_block(matchar, atomblock, atoma, atomb)
      !========================================================================!
      ! Writes a V/O matrix atomblock to a file.                               !
      ! The silent assumption is every proc can do I/O.                        !
      ! Errors abort.                                                          !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   matchar (in):    'V' or 'O'.                                         !
      !   atomblock (in):  The atomblock data to write.                        !
      !   atoma, atomb (in): Indices to the block in SFC global scheme.        !
      !========================================================================!

      use utils, only: utils_unit, utils_assert

      implicit none

      ! jd: Arguments
      character, intent(in) :: matchar
      real(kind=DP), intent(in) :: atomblock(swri%quality%max_sws_per_centre, &
          swri%quality%max_sws_per_centre)
      integer, intent(in) :: atoma, atomb

      ! jd: Local variables
      integer :: blk_unit
      character(len=1024) :: filename
      integer :: ierr

      ! --------------------------------------------------------------------

      call swri_pos_to_filename(filename, &
           elements(par%orig_atom(atoma))%centre, &
           elements(par%orig_atom(atomb))%centre, matchar, swri%swri_name)

      blk_unit = utils_unit()
      open(blk_unit, file=filename, form='unformatted', &
               action='write', iostat=ierr, status='replace')
      call utils_assert(ierr == 0, &
           'Error in internal_write_block (sw_resolution_of_identity_mod): &
           &creating file "'//trim(filename)//'" failed with code ', ierr)

      write(blk_unit, iostat=ierr) atomblock
      call utils_assert(ierr == 0, &
           'Error in internal_write_block (sw_resolution_of_identity_mod): &
           &writing file "'//trim(filename)//'" failed with code ', ierr)

      close(blk_unit, iostat=ierr)
      call utils_assert(ierr == 0, &
           'Error in internal_write_block (sw_resolution_of_identity_mod): &
           &closing file "'//trim(filename)//'" failed with code ', ierr)

    end subroutine internal_write_block

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_read_block(matchar, atomblock, atoma, atomb, error)
      !========================================================================!
      ! Reads a V/O matrix atomblock from a file.                              !
      ! The silent assumption is every proc can do I/O.                        !
      ! If 'error' is absent, aborts on error.                                 !
      ! If 'error' is present, returns 0 if succeeded, or >0 if failed.        !
      ! Tentative list of error meanings:                                      !
      ! 1 - error opening file.                                                !
      ! 2 - error reading file.                                                !
      ! 3 - error closing file.                                                !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   matchar (in):    'V' or 'O'.                                         !
      !   atomblock (out):  The atomblock to read into.                        !
      !                     If an error occured, garbage_real is returned in   !
      !                     all elements.                                      !
      !   atoma, atomb (in): Indices to the block in SFC global scheme.        !
      !   error (in, opt): See above.                                          !
      !========================================================================!

      use constants, only: garbage_real
      use rundat, only: pub_swri_verbose
      use utils, only: utils_unit, utils_abort, utils_assert

      implicit none

      ! jd: Arguments
      character, intent(in) :: matchar
      real(kind=DP), intent(out) :: atomblock(swri%quality%max_sws_per_centre, &
          swri%quality%max_sws_per_centre)
      integer, intent(in) :: atoma, atomb
      integer, intent(out), optional :: error

      ! jd: Local variables
      integer :: blk_unit
      character(len=1024) :: filename
      integer :: ierr

      ! --------------------------------------------------------------------

      call swri_pos_to_filename(filename, &
           elements(par%orig_atom(atoma))%centre, &
           elements(par%orig_atom(atomb))%centre, matchar, swri%swri_name)

      ! Open
      blk_unit = utils_unit()
      open(blk_unit, file=filename, form='unformatted', &
               action='read', iostat=ierr, status='old')
      if(ierr /= 0) then
         if(present(error)) then
            error = 1
            if(pub_swri_verbose) then
               write(stdout,'(a,i0,a,a)') &
                    '! Open error (', ierr,') for ', trim(filename)
            end if
            atomblock(:,:) = garbage_real
            return
         else
            call utils_abort(&
                 'Error in internal_read_block (sw_resolution_of_identity_mod): &
                 &opening file "'//trim(filename)//'" failed with code ', ierr)
         end if
      end if

      ! Read
      read(blk_unit, iostat=ierr) atomblock
      if(ierr /= 0) then
         if(present(error)) then
            error = 2
            if(pub_swri_verbose) then
               write(stdout,'(a,i0,a,a)') &
                    '! Read error (', ierr,') for ', trim(filename)
            end if
            atomblock(:,:) = garbage_real
            return
         else
            call utils_abort(&
                 'Error in internal_read_block (sw_resolution_of_identity_mod): &
                 &reading file "'//trim(filename)//'" failed with code ', ierr)
         end if
      end if

      ! Close
      close(blk_unit, iostat=ierr)
      if(ierr /= 0) then
         if(present(error)) then
            error = 3
            if(pub_swri_verbose) then
               write(stdout,'(a,i0,a,a)') &
                    '! Close error (', ierr,') for ', filename
            end if
            atomblock(:,:) = garbage_real
            return
         else
            call utils_assert(ierr == 0, &
                 'Error in internal_read_block (sw_resolution_of_identity_mod): &
                 &closing file "'//trim(filename)//'" failed with code ', ierr)
         end if
      end if

      if(present(error)) error = 0
      if(pub_swri_verbose) then
         write(stdout,'(a,a)') &
              '! Success for ', trim(filename)
      end if

    end subroutine internal_read_block

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_disassemble(swri, atomblocklist)
      !========================================================================!
      ! 'Disassembles' the metric matrix (or matrices) to atomblock files.     !
      ! Each procs writes out the blocks whose atoma it sparse-owns, so we     !
      ! assume each proc can do I/O.                                           !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri (in): The SW_RI within which we operate, containing the metric  !
      !              matrices.                                                 !
      !   atomblockslist (in): Describes the atomblocks of interest.           !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04 based on older subroutines.       !
      !========================================================================!

      use comms, only: pub_my_proc_id
      use sparse, only: sparse_get_block
      use utils, only: utils_alloc_check, utils_dealloc_check, &
           utils_flushed_string_output

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(in) :: swri
      type(ATOMBLOCK_DESCRIPTOR), intent(in), allocatable :: atomblocklist(:)

      ! jd: Local variables
      real(kind=DP), allocatable :: atomblock(:,:)
      integer :: i_atomblock
      integer :: atoma, atomb, atoma_sparse_owner
      character(len=*), parameter :: myself = 'internal_disassemble'

      ! -----------------------------------------------------------------------

      call utils_flushed_string_output('SWRI: ['//trim(swri%swri_name)//&
           '] - Disassembling metric matrix to files...')

      ! jd: Allocate workspace atomblock
      allocate(atomblock(&
           swri%quality%max_sws_per_centre,swri%quality%max_sws_per_centre),&
           stat=ierr)
      call utils_alloc_check(myself,'atomblock',ierr)

      ! jd: Loop over atomblocks in the list, write out an atomblock if I am
      !     its sparse owner.
      do i_atomblock = 1, size(atomblocklist)
         atoma = atomblocklist(i_atomblock)%atoma
         atomb = atomblocklist(i_atomblock)%atomb
         atoma_sparse_owner = atomblocklist(i_atomblock)%atoma_sparse_owner

         if(pub_my_proc_id == atoma_sparse_owner) then
            if(swri%has_metric(SW_V)) then
               call sparse_get_block(atomblock, &
                    swri%full_metric_matrices(SW_V), atomb, atoma)
               call internal_write_block('V', atomblock, atoma, atomb)
            end if
            if(swri%has_metric(SW_O)) then
               call sparse_get_block(atomblock, &
                    swri%full_metric_matrices(SW_O), atomb, atoma)
               call internal_write_block('O', atomblock, atoma, atomb)
            end if
         end if
       end do

       deallocate(atomblock, stat=ierr)
       call utils_dealloc_check(myself,'atomblock',ierr)

       call utils_flushed_string_output('done.'//CRLF)

    end subroutine internal_disassemble

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_erase_atomblock_file(swri, atoma, atomb)
      !========================================================================!
      ! Makes an attempt to erase an atomblock file from disk. Errors are      !
      ! ignored.
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   swri (in): The SW_RI within which we operate.                        !
      !   atoma, atomb: Atoms in the atomblock.                                !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in 2016.04.                                  !
      !========================================================================!

      use utils, only: utils_erase_file

      implicit none

      ! jd: Arguments
      type(SW_RI), intent(in) :: swri
      integer, intent(in) :: atoma, atomb

      ! jd: Local variables
      character(len=1024) :: filename
      character(len=*), parameter :: myself = 'internal_erase_atomblock_file'

      ! -----------------------------------------------------------------------

      if(swri%has_metric(SW_V)) then
         call swri_pos_to_filename(filename, &
              elements(par%orig_atom(atoma))%centre, &
              elements(par%orig_atom(atomb))%centre,'V', swri%swri_name)
         call utils_erase_file(filename)
      end if
      if(swri%has_metric(SW_O)) then
         call swri_pos_to_filename(filename, &
              elements(par%orig_atom(atoma))%centre, &
              elements(par%orig_atom(atomb))%centre,'O', swri%swri_name)
         call utils_erase_file(filename)
      end if

    end subroutine internal_erase_atomblock_file

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine swri_init_container_stage_2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_destroy_container(swri)
    !==========================================================================!
    ! Destroys an SW_RI container.                                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 18/06/2012.                                 !
    ! Modified by Jacek Dziedzic on 08/04/2014.                                !
    !==========================================================================!

    use constants, only: SW_V, SW_O
    use neighbour_list, only: neighbour_list_free
    use rundat, only: pub_use_hfx, pub_use_activehfx
    use sparse, only: sparse_destroy
    use utils, only: utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout)  :: swri

    ! jd: Local variables
    integer :: ierr
    character(len=*), parameter :: myself = 'swri_destroy_container'

    ! -------------------------------------------------------------------------

    if(.not. swri%created) return

    deallocate(swri%image_map,stat=ierr)
    call utils_dealloc_check(myself,'image_map',ierr)

    ! jd: Destroy metric matrices, unless they were destroyed earlier.
    !     HFx will destroy theirs on their own once they are re-distributed to
    !     its native parallel scheme. This has not been implemented yet.
    if(.not. swri%full_metric_matrices_destroyed) then
       if (swri%has_metric(SW_V)) then
          call sparse_destroy(swri%full_metric_matrices(SW_V))
       end if
       if (swri%has_metric(SW_O)) then
          call sparse_destroy(swri%full_metric_matrices(SW_O))
       end if
    end if
    swri%full_metric_matrices_destroyed = .true.

    ! jd: Free persistent caches
    call swri_free_persistent_caches(swri)

    ! jd: Free persistent neighbour lists
    call neighbour_list_free(swri%s_atoms_nl)

    ! jd: Clean up structures
    deallocate(swri%radtable,stat=ierr)
    call utils_dealloc_check(myself,'swri%radtable',ierr)
    deallocate(swri%sphbessels,stat=ierr)
    call utils_dealloc_check(myself,'swri%sphbessels',ierr)
    deallocate(swri%sw_idx_to_bessel_idx,stat=ierr)
    call utils_dealloc_check(myself,'swri%sw_idx_to_bessel_idx',ierr)
    deallocate(swri%sw_idx_to_m,stat=ierr)
    call utils_dealloc_check(myself,'swri%sw_idx_to_m',ierr)
    deallocate(swri%sw_idx_to_l,stat=ierr)
    call utils_dealloc_check(myself,'swri%sw_idx_to_l',ierr)
    deallocate(swri%sw_idx_to_q_idx,stat=ierr)
    call utils_dealloc_check(myself,'swri%sw_idx_to_q_idx',ierr)
    deallocate(swri%num_q_per_l,stat=ierr)
    call utils_dealloc_check(myself,'swri%num_q_per_l',ierr)

    ! jd: Destroy interpolation of SW radial components. Only used in HFx.
    if(pub_use_hfx .or. pub_use_activehfx) then
       deallocate(swri%rad_lookup,stat=ierr)
       call utils_dealloc_check(myself,'swri%rad_lookup',ierr)
    end if

    swri%swri_name = '[destroyed]'

    swri%created = .false.
    swri%initialised = .false.

  end subroutine swri_destroy_container

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function swri_get_handle_to(swri_name)
    !==========================================================================!
    ! Returns an integer handle to an SW_RI specified by name.                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2015.                                 !
    !==========================================================================!

    use utils, only: utils_abort, utils_assert

    implicit none

    ! jd: Arguments
    character(len=*), intent(in) :: swri_name

    ! jd: Local variables
    integer :: swri_h
    character(len=*), parameter :: myself = 'swri_get_handle_to'

    ! -------------------------------------------------------------------------

    ! jd: Trap attempts to call between stage_0 and stage_1, where our internal
    !     state is somewhat inconsistent: swri_library_size is set, but the
    !     library itself has not been allocated
    call utils_assert(allocated(swri_library), &
         myself//': Cannot be used before swri_init_container_stage_1')

    do swri_h = 1, swri_library_size
       if(trim(swri_name) == trim(swri_library(swri_h)%swri_name)) then
          swri_get_handle_to = swri_h
          return
       end if
    end do

    call utils_abort(myself//&
         ': ['//trim(swri_name)//'] not found in swri library.')

  end function swri_get_handle_to


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------------!
  ! ****               L O W - L E V E L    R O U T I N E S              **** !
  !---------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine swri_vmatrixblock_sc(swri,vatomblock,centrex)
    !==========================================================================!
    ! This subroutine calculates a same centre V matrix atomblock analytically.!
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in March 2010.                                   !
    ! Adapted to SW_RI by Jacek Dziedzic in February 2015.                     !
    !==========================================================================!

    use constants, only: PI
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(in)       :: swri
    real(kind=DP), intent(out)    :: vatomblock(:,:)
    type(atom_centre), intent(in) :: centrex

    ! jd: Local variables
    integer :: swx, swy       ! x and y spherical wave index
    integer :: species_number
    integer :: xlval, ylval   ! x and y l value
    integer :: xsbidx, ysbidx ! x and y spherical bessel index
    integer :: xmval, ymval   ! x and y m value
    real(kind=DP) :: aval
    real(kind=DP) :: qval
    real(kind=DP) :: qa
    real(kind=DP) :: swxpart
    real(kind=DP) :: diagel
    character(len=*), parameter :: myself = 'swri_vmatrixblock_sc'

    ! ------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    species_number = centrex%species_number

    vatomblock = 0.0_DP
    swx = 0
    aval = swri%sphbessels(species_number,1)%aval
    xsb: do xsbidx = 1, swri%max_num_bessels
       qval = swri%sphbessels(species_number,xsbidx)%qval
       xlval = swri%sphbessels(species_number,xsbidx)%lval
       qa = qval * aval

       ! jd: This likely realizes (5.7.3) from qoh's thesis

       select case (xlval)
       case(0)
          swxpart = -4.0_DP * PI * aval * cos(qa) / (qval**2)
          diagel = 4.0_DP*PI * aval/(2.0_DP*(qval**4)) &
               +swxpart * swri%sphbessels(species_number,xsbidx)%nearpotint
       case(1)
          swxpart = -4.0_DP*PI/(3.0_DP*qval**4) * sin(qa)*(qa*qa)
          diagel = 4.0_DP*PI/(qval**6) *(qa*qa - sin(qa)**2)/&
               (2.0_DP*aval)&
               + swxpart * swri%sphbessels(species_number,xsbidx)%nearpotint
       case(2)
          swxpart = 4.0_DP*PI/5.0_DP*&
               swri%sphbessels(species_number,xsbidx)%farpotint
          diagel = PI /(qa**3 * qval**5) * (qa*qa*(2.0_DP*qa*qa - 12.0_DP) + &
               (sin(qa))**2 *(qa*qa*(qa*qa*2.0_DP/3.0_DP + 2.0_DP) +12.0_DP))&
               + swxpart * swri%sphbessels(species_number,xsbidx)%nearpotint
       case(3)
          swxpart = 4.0_DP*PI/7.0_DP*&
               swri%sphbessels(species_number,xsbidx)%farpotint
          diagel = 4.0_DP*PI*((2.0_DP*(-45.0_DP - 15.0_DP*qa**2 - &
               6.0_DP*qa**4 + qa**6) + 6.0_DP*(15.0_DP - 25.0_DP*qa**2 + &
               2.0_DP*qa**4)*cos(2.0_DP*qa) + qa*(180.0_DP - 60.0_DP*qa**2 + &
               qa**4)*sin(2.0_DP*qa)))/ (4.0_DP*qa**5*qval**5) + &
               swxpart * swri%sphbessels(species_number,xsbidx)%nearpotint
       case(4)
          swxpart = 4.0_DP*PI/9.0_DP*&
               swri%sphbessels(species_number,xsbidx)%farpotint
          diagel = 4.0_DP*PI*((-3150.0_DP - 630.0_DP*qa**2 &
               - 90.0_DP*qa**4 - 20.0_DP*qa**6 + 2.0_DP*qa**8) &
               + 10.0_DP*(315.0_DP - 567.0_DP*qa**2 + 93.0_DP*qa**4 &
               - 2.0_DP*qa**6)*cos(2.0_DP*qa) &
               + qa*(6300.0_DP - 2940.0_DP*qa**2 + &
               180.0_DP*qa**4 - qa**6)*sin(2.0_DP*qa))/ (4.0_DP*qa**7*qval**5)&
               + swxpart * swri%sphbessels(species_number,xsbidx)%nearpotint
       case default
          swxpart = 0.0_DP
          diagel = 0.0_DP
          call utils_abort(trim(myself)//': Value of xlval too big.')
          exit
       end select
       xm: do xmval = -xlval,xlval
          swx = swx + 1
          vatomblock(swx,swx) = diagel
          swy = 0
          ysb: do ysbidx= 1, swri%max_num_bessels
             ylval = swri%sphbessels(species_number,ysbidx)%lval
             ym: do ymval= -ylval,ylval
                swy = swy + 1
                if (xlval /= ylval .or. xmval /= ymval .or. swx == swy) cycle
                vatomblock(swx,swy) =  swxpart * &
                     swri%sphbessels(species_number,ysbidx)%nearpotint
             end do ym ! jd: m's on y
          end do ysb ! jd: Bessels on y
       end do xm ! jd: m's on x
    end do xsb ! jd: Bessels on x

    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swri_vmatrixblock_sc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_omatrixblock_sc(swri, oatomblock, centrex)
    !==========================================================================!
    ! This subroutine calculates a same centre O matrix atomblock analytically.!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2012 basing on a similar subroutine     !
    ! for the vmatrix, by Quintin Hill.                                        !
    ! Adapted to SW_RI by Jacek Dziedzic in February 2015.                     !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(in)       :: swri
    real(kind=DP), intent(out)    :: oatomblock(:,:)
    type(atom_centre), intent(in) :: centrex

    ! jd: Local variables
    integer :: swx
    integer :: species_number
    integer :: xlval
    integer :: xsbidx
    integer :: xmval
    real(kind=DP) :: aval
    real(kind=DP) :: qval
    real(kind=DP) :: qa
    real(kind=DP) :: ovlp
    character(len=*), parameter :: myself = 'swri_omatrixblock_sc'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)

    species_number = centrex%species_number

    oatomblock = 0.0_DP
    swx = 0
    aval = swri%sphbessels(species_number,1)%aval

    do xsbidx = 1, swri%max_num_bessels
       qval = swri%sphbessels(species_number,xsbidx)%qval
       xlval = swri%sphbessels(species_number,xsbidx)%lval
       qa = qval * aval

       ! jd: Calculated in Mathematica. This is the overlap of two SWs on the
       !     same centre. The angular part is
       !     \int sin(phi) (Z_lm(phi,theta))^2 dphi dtheta and is 1.
       !     The radial part is \int r^2(j_l(qr))^2 dr. This is what this
       !     subroutine calculates.
       !     Results verified against Mathematica by comparing calculated
       !     values to 12 digits for l=0..4 and first two q's.
       select case (xlval)
       case(0)
          ovlp = (2.0_DP*aval*qval - sin(2.0_DP*qa)) / (4.0_DP * qval**3.0_DP)
       case(1)
          ovlp = (-2.0_DP + 2.0_DP*qa*qa + 2.0_DP*cos(2.0_DP*qa) + &
               qa*sin(2.0_DP*qa)) / (4.0_DP * aval * qval**4.0_DP)
       case(2)
          ovlp = (-6.0_DP - 6.0_DP*qa*qa + 2.0_DP*qa**4.0_DP + &
               (6.0_DP - 6.0_DP*qa*qa)*cos(2.0_DP*qa) + qa*(12.0_DP - qa*qa) * &
               sin(2.0_DP*qa))/(4.0_DP*aval**3.0_DP*qval**6.0_DP)
       case(3)
          ovlp = (2.0_DP*(-45.0_DP - 15.0_DP*qa*qa - 6.0_DP*qa**4.0_DP + &
               qa**6.0_DP) + 6.0_DP*(15.0_DP - 25.0_DP*qa*qa + &
               2.0_DP*qa**4.0_DP)*cos(2.0_DP*qa) + &
               qa*(180.0_DP - 60.0_DP *qa*qa + qa**4.0_DP) * sin(2.0_DP*qa)) / &
               (4.0_DP * aval**5.0_DP * qval**8.0_DP)
       case(4)
          ovlp = -((3150.0_DP + 630.0_DP*qa*qa + 90.0_DP*qa**4.0_DP + &
               20.0_DP*qa**6.0_DP - 2.0_DP*qa**8.0_DP + 10.0_DP*(-315.0_DP + &
               567.0_DP*qa*qa - 93.0_DP*qa**4.0_DP + 2.0_DP*qa**6.0_DP) * &
               cos(2.0_DP*qa) + qa*(-6300.0_DP + 2940.0_DP*qa*qa - 180.0_DP* &
               qa**4.0_DP + qa**6.0_DP) * sin (2.0_DP*qa)) / &
               (4.0_DP*aval**7.0_DP*qval**10.0_DP))
       case default
          ovlp = 0.0_DP
          exit
       end select

       do xmval = -xlval,xlval
          swx = swx + 1
          oatomblock(swx,swx) = ovlp
       end do ! jd: m's on x

    end do ! jd: Bessels on x

    call utils_trace_out(myself)

  end subroutine swri_omatrixblock_sc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_expansion_centres(swri, centres_out, atoms_out, &
       num_sws_in_expansion_out, centres_in, atoms_in, swexside_filter, &
       want_sorted)
    !==========================================================================!
    ! Determines the unique centres (up to 2) for an SW expansion on two       !
    ! centres. The centres and atom indices are sorted, unless 'want_sorted'   !
    ! is provided and is .false. (useful in the absence of Bb/Cc symmetry,     !
    ! e.g. with mixed NGWFs).                                                  !
    ! For single-centre expansion the second value will be -1. The number of   !
    ! SWs in the expansion is returned in num_sws_in_expansion_out. This is    !
    ! filtered down to a given swex quality by specifying a swexside_filter.   !
    !--------------------------------------------------------------------------!
    ! This subroutine is unnecessarily complicated, because it used to support !
    ! larger numbers of centres. Might make sense to simplify it.              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2012 to support n centres.            !
    ! Modified by Jacek Dziedzic in June 2012 to support 2 centres only.       !
    ! Modified by Jacek Dziedzic in February 2015 to support swexside_filter.  !
    ! Modified by Jacek Dziedzic in July 2018 to allow centre sorting to be    !
    ! disabled.                                                                !
    !==========================================================================!

    use geometry, only: geometry_distance

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(in)        :: swri
    type(ATOM_CENTRE), intent(out) :: centres_out(2)
    integer, intent(out)           :: atoms_out(2)
    integer, intent(out)           :: num_sws_in_expansion_out
    type(atom_centre), intent(in)  :: centres_in(2)
    integer, intent(in)            :: atoms_in(2)
    type(SW_QUALITY), intent(in)   :: swexside_filter
    logical, intent(in)            :: want_sorted

    ! jd: Local variables
    integer           :: sbidx
    integer           :: qidx
    integer           :: lval
    integer           :: i, j
    logical           :: is_unique
    integer           :: n_out
    integer           :: min_atom_loc
    type(atom_centre) :: centres_temp(2)
    integer           :: atoms_temp(2)
    character(len=*), parameter :: myself = 'swri_expansion_centres'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)

    ! jd: Ensure *all* centres_out are initialized to something, since this is
    !     an intent(out) argument. Below we only fill *up to* 2 elements with
    !     meaningful data. The remaining elements of this array will never be
    !     accessed anyway, as the corresponding elements in atoms_out are -1,
    !     but to placate the compiler we'll fill all elements of centres_out
    centres_temp(:) = centres_in(:)

    atoms_temp(:) = 1999999999

    num_sws_in_expansion_out = 0
    n_out = 1

    do i = 1, 2

       is_unique = .true.
       do j = 1, i-1
          if(.not. geometry_distance( &
               centres_in(i)%incell, &
               centres_in(j)%incell) > epsilon(1.0_DP)) is_unique=.false.
       end do

       if(is_unique) then

          centres_temp(n_out) = centres_in(i)
          atoms_temp(n_out) = atoms_in(i)

          do sbidx=1, swri%max_num_bessels
             lval = swri%sphbessels(centres_in(i)%species_number,sbidx)%lval
             if (lval == -1) exit
             if (lval < swexside_filter%min_l) cycle
             if (lval > swexside_filter%max_l) cycle

             qidx = swri%sphbessels(centres_in(i)%species_number,sbidx)%qidx
             if (qidx > swexside_filter%max_q) cycle

             num_sws_in_expansion_out = num_sws_in_expansion_out + 2*lval + 1
          end do

          n_out = n_out + 1

       end if

    end do

    ! jd: Sort the atoms and centres in ascending order so that the result
    !     never depends on how input centres are permuted.
    if(want_sorted) then
       do i = 1, 2
          min_atom_loc = minloc(atoms_temp,1)
          if(atoms_temp(min_atom_loc) == 1999999999) then
             atoms_out(i) = -1
          else
             atoms_out(i) = atoms_temp(min_atom_loc)
          end if
          centres_out(i) = centres_temp(min_atom_loc)
          atoms_temp(min_atom_loc) = 2000000000
       end do
    else
       ! jd: Here we do not sort, but if there's only one atom, we make sure
       !     it's first, and -1 is second
       if(atoms_temp(1) /= 1999999999 .and. atoms_temp(2) /= 1999999999) then
          ! jd: Two-centre expansion
          atoms_out(1:2) = atoms_temp(1:2)
          centres_out(1:2) = centres_temp(1:2)
       end if
       if(atoms_temp(1) /= 1999999999 .and. atoms_temp(2) == 1999999999) then
          ! jd: One-centre expansion, 2nd atom absent
          atoms_out(1) = atoms_temp(1)
          atoms_out(2) = -1
          centres_out(1:2) = centres_temp(1:2)
       end if
       if(atoms_temp(1) == 1999999999 .and. atoms_temp(2) /= 1999999999) then
          ! jd: One-centre expansion, 1st atom absent
          atoms_out(1) = atoms_temp(2)
          atoms_out(2) = -1
          centres_out(1) = centres_temp(2)
          centres_out(2) = centres_temp(1)
       end if
    end if

    call utils_trace_out(myself)

  end subroutine swri_expansion_centres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_init_centre(atomcentre, par, elements, atom)
    !==========================================================================!
    ! This subroutine initialises centre for the current atom.                 !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in September 2009                                !
    ! Modified to run on multiple cores on 13 October 2009 by Quintin Hill.    !
    ! Extensively prunned by Jacek Dziedzic in 2014/11.                        !
    !==========================================================================!

    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO

    implicit none

    type(atom_centre), intent(out) :: atomcentre
    type(PARAL_INFO), intent(in)   :: par
    type(ELEMENT),    intent(in)   :: elements(par%nat)
    integer,          intent(in)   :: atom

    atomcentre%species_number = elements(par%orig_atom(atom))%species_number
    atomcentre%global_idx = atom
    atomcentre%incell = elements(par%orig_atom(atom))%centre
    atomcentre%radius = elements(par%orig_atom(atom))%radius

  end subroutine swri_init_centre

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_sw_to_ablmq(swri, sw, a, b, l, m, q)
    !==========================================================================!
    ! Decomposes a swriside SW index into a Bessel number (b), and quantum     !
    ! numbers l, m. Also returns Bessel 'q' and 'a'.                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2016.                               !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(in)    :: swri
    integer, intent(in)        :: sw
    real(kind=DP), intent(out) :: a
    integer, intent(out)       :: b
    integer, intent(out)       :: l
    integer, intent(out)       :: m
    real(kind=DP), intent(out) :: q

    ! jd: Local variables
    integer :: species

    ! -------------------------------------------------------------------------

    species = 1 ! @IDENTICAL_RADII
    b = swri%sw_idx_to_bessel_idx(sw)
    l = swri%sphbessels(species,b)%lval
    q = swri%sphbessels(species,b)%qval
    a = swri%sphbessels(species,b)%aval
    m = swri%sw_idx_to_m(sw)

  end subroutine swri_sw_to_ablmq

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! @2DN_1DA_OVERLAP
  ! TODO In order to support 2Dn-1Da scheme with overlap metric, a variant of
  !      this routine for swri_eval_swrt_at_nodes which supports evaluation of
  !      SWRTs in a half-disc (2D) should be created (in analogy to the SWpot
  !      version swri_eval_swpotrt_at_nodes_on_halfdisc). This routine should
  !      be called swri_eval_swrt_at_nodes_on_halfdisc.
  subroutine swri_eval_sw_at_nodes(values, cell, nodes, l, m, q, a, disp)
    !==========================================================================!
    ! Calculates the density of a SW at Chebyshev nodes of a sphere.           !
    ! If the optional argument 'disp' is omitted, the centre of the SW is at   !
    ! the origin, otherwise the SW is centred at 'disp'.                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   values (out): The returned values of the SW at Chebyshev nodes.        !
    !   cell (in):    The simulation cell. Needed for PBCs.                    !
    !   nodes (in):   The Chebyshev nodes where the SW is to be evaluated.     !
    !   l, m, q (in): Indices of the source SW.                                !
    !   a (in):       The extent of the source SW.                             !
    !   disp (in, optional): The displacement vector to the source SW.         !
    !                        Assumed to be zero if omitted.                    !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !   For nodes which are outside the localization sphere of the source SW   !
    !   this subroutine will return zero without calculating the Bessel or the !
    !   harmonic, but it still takes a while to process these points. Note,    !
    !   however, that if the SW is so far away that there is *no* overlap      !
    !   whatsoever, it's better to elide the call to this subroutine entirely. !
    !   This is exactly what internal_vomatrixblock_cheb() does.               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2012.                               !
    ! Extended to PBCs by Jacek Dziedzic in July 2022.                         !
    !==========================================================================!

    use chebyshev_rep, only: CHEB_NODES
    use geometry, only: geometry_magnitude, unit_vector
    use rundat, only: pub_hfx_bc_is_periodic
    use simulation_cell, only: CELL_INFO, minimum_image_displacement
    use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm_unit
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)   :: values(:) ! Values for all xnodes of this SW
    type(CELL_INFO), intent(in)  :: cell
    type(CHEB_NODES), intent(in) :: nodes ! Chebyshev nodes
    integer, intent(in)          :: l, m             ! Quantum numbers of SW
    real(kind=DP), intent(in)    :: q, a       ! Bessel q and a of the SW
    type(POINT), intent(in), optional :: disp ! Vector to the centre of the SW

    ! jd: Local variables
    integer       :: xi, yi, zi
    integer       :: n_stripes(3)
    integer       :: loc
    real(kind=DP) :: r
    real(kind=DP) :: sphbess
    type(POINT)   :: pos, unit_pos, SW_disp, nodept
    character(len=*), parameter :: myself = 'swri_eval_sw_at_nodes'

    ! --------------------------------------------------------------------------


    call utils_trace_in(myself)
    call timer_clock(myself,1)

    if(present(disp)) then
       SW_disp = disp
    else
       SW_disp%x = 0.0_DP
       SW_disp%y = 0.0_DP
       SW_disp%z = 0.0_DP
    end if

    n_stripes(1:3) = nodes%ranges%n_stripes(1:3)

    if(.not. any(pub_hfx_bc_is_periodic) .or. .not. present(disp)) then
       ! jd: OBC; or center at zero, meaning PBCs don't matter -- entire sphere
       !          is guaranteed to fit in cell
       do zi = 1, n_stripes(3)
          pos%z = nodes%znodes(zi) - SW_disp%z

          do yi = 1, n_stripes(2)
             pos%y = nodes%ynodes(yi,zi) - SW_disp%y

             do xi = 1, n_stripes(1)
                pos%x = nodes%xnodes(xi,yi,zi) - SW_disp%x

                loc = xi + (yi-1) * n_stripes(1) + (zi-1) * n_stripes(1)*n_stripes(2)

                unit_pos = unit_vector(pos) ! @optimize: unit_vector calcs magnitude already
                r = geometry_magnitude(pos)

                if( r <= a ) then  ! jd: Within the sphere of the SW
                  sphbess = sw_bessel_fast(l,r*r*q*q) * (-1.0_DP)**l
                  values(loc) = sphbess * &
                       sw_real_sph_harm_unit(unit_pos%x,unit_pos%y,unit_pos%z,l,m)
                else               ! jd: Outside the sphere of the SW
                   values(loc) = 0.0_DP
                end if

             end do
          end do
       end do

    else if(all(pub_hfx_bc_is_periodic)) then
       ! jd: PBC
       do zi = 1, n_stripes(3)
          nodept%z = nodes%znodes(zi)

          do yi = 1, n_stripes(2)
             nodept%y = nodes%ynodes(yi,zi)

             do xi = 1, n_stripes(1)
                nodept%x = nodes%xnodes(xi,yi,zi)

                loc = xi + (yi-1) * n_stripes(1) + (zi-1) * n_stripes(1)*n_stripes(2)

                call minimum_image_displacement(r, pos, nodept, SW_disp, cell)

                unit_pos = unit_vector(pos)

                if( r <= a ) then  ! jd: Within the sphere of the SW
                  sphbess = sw_bessel_fast(l,r*r*q*q) * (-1.0_DP)**l
                  values(loc) = sphbess * &
                       sw_real_sph_harm_unit(unit_pos%x,unit_pos%y,unit_pos%z,l,m)
                else               ! jd: Outside the sphere of the SW
                   values(loc) = 0.0_DP
                end if

             end do
          end do
       end do

    else
       call utils_abort(myself//': Mixed BCs are not currently supported.')
    end if


    call timer_clock(myself,2)
    call utils_trace_out(myself)

  end subroutine swri_eval_sw_at_nodes

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_eval_sw_at_nodes_fast(values, nodes, l, m, q)
    !==========================================================================!
    ! Calculates the density of a SW at Chebyshev nodes of a sphere.           !
    ! Differs from swri_eval_sw_at_nodes() in that it does not support 'disp', !
    ! which can be used to make it more efficient.                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   values (out): The returned values of the SW at Chebyshev nodes.        !
    !   nodes (in):   The Chebyshev nodes where the SW is to be evaluated.     !
    !   l, m, q (in): Indices of the source SW.                                !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !   Since this operates on the home sphere, it can be assumed that the     !
    !   nodes lie inside the localisation sphere. This is not asserted for     !
    !   to squeeze out extra efficiency.                                       !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2013, based on swx_eval_at_nodes.  !
    !==========================================================================!

    use chebyshev_rep, only: CHEB_NODES
    use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm_unit

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: values(:) ! Values for all xnodes of this SW
    type(CHEB_NODES), intent(in) :: nodes ! Chebyshev nodes
    integer, intent(in) :: l, m             ! Quantum numbers of SW
    real(kind=DP), intent(in) :: q          ! Bessel q of the SW

    ! jd: Local variables
    integer       :: xi, yi, zi
    integer       :: n_stripes(3)
    integer       :: loc
    real(kind=DP) :: r, inv_r
    real(kind=DP) :: sphbess
    real(kind=DP) :: q2
    real(kind=DP) :: l_factor
    real(kind=DP) :: x,y,z,z2,y2_plus_z2
    real(kind=DP) :: unit_x, unit_y, unit_z
    integer       :: yz_offset, z_offset
!    character(len=*), parameter :: myself = 'swri_eval_sw_at_nodes_fast'

    ! --------------------------------------------------------------------------

    n_stripes(1:3) = nodes%ranges%n_stripes(1:3)
    q2 = q*q
    l_factor = (-1.0_DP)**l

    do zi = 1, n_stripes(3)
       z = nodes%znodes(zi)
       z2 = z*z
       z_offset = (zi-1) * n_stripes(1)*n_stripes(2)

       do yi = 1, n_stripes(2)
          y = nodes%ynodes(yi,zi)
          y2_plus_z2 = y*y + z2
          yz_offset = (yi-1) * n_stripes(1) + z_offset

          do xi = 1, n_stripes(1)
             x = nodes%xnodes(xi,yi,zi)
             loc = xi + yz_offset

             r = sqrt(x*x+y2_plus_z2) ! @MIC? Probably not
             inv_r = 1.0_DP / r
             unit_x = inv_r * x
             unit_y = inv_r * y
             unit_z = inv_r * z

             sphbess = sw_bessel_fast(l,r*r*q2) * l_factor
             values(loc) = sphbess * &
                  sw_real_sph_harm_unit(unit_x, unit_y, unit_z, l,m)

          end do
       end do
    end do

  end subroutine swri_eval_sw_at_nodes_fast

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_eval_swrt_at_nodes_on_halfdisc_fast(values, nodes, l, m, q)
    !==========================================================================!
    ! Calculates the density of the (r,theta)-dependent parts of a SW at       !
    ! Chebyshev nodes on a 2-D half-disc in the Cartesian yz plane and returns !
    ! the product of this with a function derived from Jacobian determinants   !
    ! for the coordinate transfomations required to (i) separate the 1-D       !
    ! analytic integral from the 2-D numerical integral; and (ii) express the  !
    ! 2-D numerical integral in Cartesian coordinates.                         !
    !                                                                          !
    ! A half-disc is used, since the function is destined to be used in an     !
    ! integration over Cartesian coordinates corresponding to the polar        !
    ! coordinate ranges r in [0,a], theta in [0,pi], i.e. half a disc.         !
    !                                                                          !
    ! The values returned by the routine are not simply the SWRT, but          !
    !   J1(y,z) * J2(y,z) * SWRT(r(y,z),theta(y,z))                            !
    ! where                                                                    !
    !   * J1(y,z) is the Jacobian determinant for change of coordinates from   !
    !     2-D polar (r,theta) to 2-D Cartesian (y,z).                          !
    !   * J2(y,z) is the Jacobian determinant for change of coordinates from   !
    !     3-D Cartesian (x,y,z) to 3-D spherical polar (r,theta,phi) in        !
    !     terms of Cartesian coordinates (y,z).                                !
    !   * SWRT is the (r,theta)-dependent part of the spherical wave in        !
    !     Cartesian coordinates (see below)                                    !
    ! J1 and J2 are necessary because the the separation into 2-D numerical    !
    ! and 1-D analytic integral is done in spherical polar coordinates, i.e.   !
    ! the volume element is J2 dr dtheta dphi. The 1-D integral over phi can   !
    ! be evaluated analytically, while the remaining 2-D integral is done      !
    ! numerically. We perform this 2-D integral in Cartesian coordinates, so   !
    ! a further coordinate transformation is required, converting our volume   !
    ! element dr dtheta to J1 dx dy.                                           !
    !                                                                          !
    ! SWRT is the (r,theta)-dependent part of the spherical wave, i.e.         !
    !   Full SW:  j_l(q*r) * Z1(l,m,theta) * Z2(m,phi)                         !
    !   SWRT:     j_l(q*r) * Z1(l,m,theta)                                     !
    ! where j_l is a spherical Bessel function and Z1, Z2 components of a      !
    ! real-spherical harmonic (RSH).                                           !
    !                                                                          !
    ! The RSH is decomposed as                                                 !
    !   Z_l,m = Z1(l,m,theta) * Z2(m,phi)                                      !
    ! with                                                                     !
    !   Z1 = N_l,|m| * P_{l,|m|}(cos(theta))                                   !
    ! and                                                                      !
    !         { 1,               m == 0                                        !
    !   Z2 =  { sin(abs(m)*phi), m < 0                                         !
    !         { cos(  m   *phi), m > 0                                         !
    ! where N_l,m is a normalization constant and P_{l,m} is an associated     !
    ! Legendre function (or, more precisely, a Ferrers function, see           !
    ! https://dlmf.nist.gov/14.7.E10).                                         !
    !                                                                          !
    ! Differs from swri_eval_swrt_at_nodes_on_halfdisc() in that it does not   !
    ! support 'disp', which can be used to make it more efficient.             !
    ! (@2DN_1DA_OVERLAP This routine has yet to be implemented, since it is    !
    ! not required to evaluate the electrostatic metric. Instead, the          !
    ! routine swri_eval_swpotrt_at_nodes_on_halfdisc is used to evaluate the   !
    ! remote SWPOTRT at nodes on the halfdisc.)                                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   values (out): The returned values of the SWRT at Chebyshev nodes.      !
    !   nodes (in):   The Chebyshev nodes where the SWRT is to be evaluated.   !
    !   l, m, q (in): Indices of the source SW.                                !
    !--------------------------------------------------------------------------!
    ! Notes:                                                                   !
    !   Since this operates on the home sphere, it can be assumed that the     !
    !   nodes lie inside the localisation sphere. This is not asserted for     !
    !   to squeeze out extra efficiency.                                       !
    !--------------------------------------------------------------------------!
    ! Written by James C. Womack in June 2018.                                 !
    ! Derived from swri_eval_sw_at_nodes_fast by Jacek Dziedzic, Nov 2013.     !
    !==========================================================================!

    use chebyshev_rep, only: CHEB_NODES
    use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm_theta_part_unit

    implicit none

    ! Parameters
    !character(len=*), parameter :: myself = &
    !    "swri_eval_swrt_at_nodes_on_halfdisc_fast"

    ! Arguments
    real(kind=DP), intent(out) :: values(:) ! Values for all xnodes of this SW
    type(CHEB_NODES), intent(in) :: nodes   ! Chebyshev nodes
    integer, intent(in) :: l, m             ! Quantum numbers of SW
    real(kind=DP), intent(in) :: q          ! Bessel q of the SW

    ! Local variables
    integer       :: yi, zi
    integer       :: n_stripes(3)
    integer       :: loc
    real(kind=DP) :: r, inv_r
    real(kind=DP) :: sphbess
    real(kind=DP) :: q2
    real(kind=DP) :: l_factor
    real(kind=DP) :: y,z,z2,y2_plus_z2
    real(kind=DP) :: unit_x, unit_y, unit_z
    integer       :: z_offset
    real(kind=DP) :: j1_j2 ! Product of J1 and J2 (see documentation above)

    ! JCW: n_stripes(1) == 1 for integration on a disc
    n_stripes(1:3) = nodes%ranges%n_stripes(1:3)
    q2 = q*q
    l_factor = (-1.0_DP)**l

    ! JCW: evaluate (r,theta)-dependent part of real-spherical harmonic in
    ! JCW: yz plane (x = 0)
    unit_x = 0.0_DP
    do zi = 1, n_stripes(3)
       z = nodes%znodes(zi)
       z2 = z*z
       z_offset = (zi-1) * n_stripes(2)

       do yi = 1, n_stripes(2)
          y = nodes%ynodes(yi,zi)
          y2_plus_z2 = y*y + z2

          loc = yi + z_offset

          r = sqrt(y2_plus_z2) ! @MIC? Likely not
          inv_r = 1.0_DP / r
          unit_y = inv_r * y
          unit_z = inv_r * z

          ! JCW: J1(y,z) is the Jacobian determinant for change of coordinates
          ! JCW: from 2-D polar (r,theta) to 2-D Cartesian (y,z), i.e.
          ! JCW:    dr dtheta = J1 dy dz
          ! JCW: J1 = 1/r = 1/sqrt(y^2+z^2), which can be easily seen from the
          ! JCW: well-known Jacobian determinant for Cartesian to polar
          ! JCW: transformation:
          ! JCW:    dy dz = r dr dtheta
          ! JCW:
          ! JCW: J2(y,z) is the Jacobian determinant for change of coordinates
          ! JCW: from 3-D Cartesian (x,y,z) to 3-D spherical polar
          ! JCW: (r,theta,phi), i.e.
          ! JCW:    dx dy dz = r^2 sin(theta) dr dtheta dphi
          ! JCW:
          ! JCW: Thus,
          ! JCW:    J1 * J2 = r sin(theta)
          ! JCW: which is, in Cartesian coords
          ! JCW:    J1 * J2 = sqrt(y^2 + z^2) * y/sqrt(y^2 + z^2)
          ! JCW:            = y
          j1_j2 = y

          sphbess = sw_bessel_fast(l,r*r*q2) * l_factor
          values(loc) = sphbess * &
               sw_real_sph_harm_theta_part_unit(unit_x, unit_y, unit_z, l,m) * &
               j1_j2

       end do
    end do

  end subroutine swri_eval_swrt_at_nodes_on_halfdisc_fast

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_eval_swpot_at_nodes(swri, values, nodes, l, m, disp, &
       species_number, bessel_idx)
    !==========================================================================!
    ! Calculates the potential of a SW at Chebyshev nodes of a sphere.         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swri (in):    The swri container (needed for sphbessels).
    !   values (out): The returned values of the SWpot at Chebyshev nodes.     !
    !   nodes (in):   The Chebyshev nodes where the SWpot is to be evaluated.  !
    !   l, m (in):    Indices of the source SW.                                !
    !   disp (in):    The displacement vector to the source SW.                !
    !   species_number (in): Species of the atom generating the SWpot.         !
    !   bessel_idx (in): The number of the Bessel corresponding to l, q.       !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2012.                               !
    !==========================================================================!

    use chebyshev_rep, only: CHEB_NODES
    use constants, only: PI
    use geometry, only: geometry_magnitude, unit_vector, operator(-)
    use rundat, only: pub_debug
    use simulation_cell, only: minimum_image_displacement
    use spherical_wave, only: sw_real_sph_harm_unit
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(in)        :: swri
    real(kind=DP), intent(out)     :: values(:) ! Returned values
    type(CHEB_NODES), intent(in)   :: nodes     ! Chebyshev nodes
    integer, intent(in)            :: l, m      ! Which SWpot
    type(POINT), intent(in)        :: disp      ! Displacement to SWpot
    integer, intent(in)            :: species_number  ! Species and
    integer, intent(in)            :: bessel_idx      ! Bessel index of SWpot

    ! jd: Local variables
    integer       :: xi, yi, zi
    integer       :: n_stripes(3)
    integer       :: val_idx
    real(kind=DP) :: r
    real(kind=DP) :: bessint
    type(POINT)   :: pos, vec, unit_vec
    real(kind=DP) :: factor
    character(len=*), parameter :: myself = 'swri_eval_swpot_at_nodes'

    ! --------------------------------------------------------------------------

    call utils_trace_in(myself)

    n_stripes(1:3) = nodes%ranges%n_stripes(1:3)

    factor = 4.0_DP * PI / real(((2*l+1) * (-1)**l),kind=DP)

    do zi = 1, n_stripes(3)

       pos%z = nodes%znodes(zi)

       do yi = 1, n_stripes(2)

          pos%y = nodes%ynodes(yi,zi)

          do xi = 1, n_stripes(1)

             pos%x = nodes%xnodes(xi,yi,zi)


             vec = pos - disp
             r = geometry_magnitude(vec)

#if 0
             ! jd: No MICs between centre and point, only between centres
             !     Leaving this code for a while, but should be removed
             !     eventually.
             if(.not. pbc) then
                ! jd: OBC
                vec = pos - disp
                r = geometry_magnitude(vec)
             else
                ! jd: PBC
                call minimum_image_displacement(r, vec, pos, disp, cell)
             end if
#endif

             unit_vec = unit_vector(vec)

             if(pub_debug) call utils_assert(r /= 0.0_DP, &
                  myself//': r is unexpectedly zero')

             val_idx = xi + (yi-1) * n_stripes(1) + (zi-1) * n_stripes(1)*n_stripes(2)

             bessint = swri_sph_bess_pot_int(r,1.0_DP/r, &
                 swri%sphbessels(species_number, bessel_idx))

             values(val_idx) = bessint * factor * &
                  sw_real_sph_harm_unit(unit_vec%x, unit_vec%y, unit_vec%z, l,m)

          end do
       end do
    end do

    call utils_trace_out(myself)

  end subroutine swri_eval_swpot_at_nodes

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_eval_swpotrt_at_nodes_on_halfdisc(swri, values, nodes, l, m, &
       disp, species_number, bessel_idx)
    !==========================================================================!
    ! Calculates the the (r,theta)-dependent parts of a SWpot at Chebyshev     !
    ! nodes on a 2-D disc in the Cartesian yz plane.                           !
    !                                                                          !
    ! Unlike swri_eval_swrt_at_nodes_on_halfdisc_fast, the values returned by  !
    ! this routine are simply the potential of the SW. The function J1 * J2    !
    ! from the coordinate transformations involved in the 2Dn-1Da scheme       !
    ! has already been applied to the SWRT evaluated in that routine.          !
    !                                                                          !
    ! SWpotRT is the (r,theta)-dependent part of the Coulombic potential of a  !
    ! spherical wave, i.e.                                                     !
    !   Full SWpot:  u(l,q,r) * Z1(l,m,theta) * Z2(m,phi)                      !
    !   SWpotRT:     u(l,q,r) * Z1(l,m,theta)                                  !
    ! where u(l,q,r) is r-dependent part of the spherical wave potential       !
    ! and Z1, Z2 are components of a real-spherical harmonic (RSH).            !
    !                                                                          !
    ! The RSH is decomposed as                                                 !
    !   Z_l,m = Z1(l,m,theta) * Z2(m,phi)                                      !
    ! with                                                                     !
    !   Z1 = N_l,|m| * P_{l,|m|}(cos(theta))                                   !
    ! and                                                                      !
    !         { 1,               m == 0                                        !
    !   Z2 =  { sin(abs(m)*phi), m < 0                                         !
    !         { cos(  m   *phi), m > 0                                         !
    ! where N_l,m is a normalization constant and P_{l,m} is an associated     !
    ! Legendre function (or, more precisely, a Ferrers function, see           !
    ! https://dlmf.nist.gov/14.7.E10).                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swri (in):    The swri container (needed for sphbessels).              !
    !   values (out): The returned values of the SWpotRT at Chebyshev nodes.   !
    !   nodes (in):   The Chebyshev nodes where the SWpotRT is to be evaluated.!
    !   l, m (in):    Indices of the source SW.                                !
    !   disp (in):    The displacement vector to the source SW (in the         !
    !                 integral frame, i.e. with z-axis aligned along           !
    !                 internuclear vector.                                     !
    !   species_number (in): Species of the atom generating the SWpotRT.       !
    !   bessel_idx (in): The number of the Bessel corresponding to l, q.       !
    !--------------------------------------------------------------------------!
    ! Written by James C. Womack in June 2018.                                 !
    ! Derived from swri_eval_swpot_at_nodes by Jacek Dziedzic, Jan 2012.       !
    !==========================================================================!

    use chebyshev_rep, only: CHEB_NODES
    use constants, only: PI
    use geometry, only: geometry_magnitude, unit_vector, operator(-)
    use rundat, only: pub_debug
    use spherical_wave, only: sw_real_sph_harm_theta_part_unit
    use utils, only: utils_assert

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = &
         'swri_eval_swpotrt_at_nodes_on_halfdisc'

    ! Arguments
    type(SW_RI), intent(in)        :: swri
    real(kind=DP), intent(out)     :: values(:) ! Returned values
    type(CHEB_NODES), intent(in)   :: nodes     ! Chebyshev nodes
    integer, intent(in)            :: l, m      ! Which SWpot
    type(POINT), intent(in)        :: disp      ! Displacement to SWpot in
                                                ! integral frame
    integer, intent(in)            :: species_number  ! Species and
    integer, intent(in)            :: bessel_idx      ! Bessel index of SWpot

    ! Local variables
    integer       :: yi, zi
    integer       :: n_stripes(3)
    integer       :: val_idx
    real(kind=DP) :: r
    real(kind=DP) :: bessint
    type(POINT)   :: pos, vec, unit_vec
    real(kind=DP) :: factor

    ! --------------------------------------------------------------------------

    call utils_trace_in(myself)

    ! JCW: n_stripes(1) == 1 for integration on a disc
    n_stripes(1:3) = nodes%ranges%n_stripes(1:3)

    factor = 4.0_DP * PI / real(((2*l+1) * (-1)**l),kind=DP)

    ! JCW: evaluate (r,theta)-dependent part of real-spherical harmonic in
    ! JCW: yz plane (x = 0)
    pos%x = 0.0_DP

    do zi = 1, n_stripes(3)

       pos%z = nodes%znodes(zi)

       do yi = 1, n_stripes(2)

          pos%y = nodes%ynodes(yi,zi)


          ! JCW: Vector from the centre of this SWpotRT to the current node
          ! JCW: in the home sphere (i.e. the half-disc portion of this which
          ! JCW: we are integrating over), i.e.
          ! JCW:   vec = pos - disp
          vec = pos - disp ! @MIC -- likely not, as we only apply it between
                           ! centres, not between centre and point.
                           ! Can also move some coords out of loop

          unit_vec = unit_vector(vec) ! @optimize: unit_vector calculates magnitude already
          r = geometry_magnitude(vec)

          if(pub_debug) call utils_assert(r /= 0.0_DP, &
               myself//': r is unexpectedly zero')

          val_idx = yi + (zi-1) * n_stripes(2)

          bessint = swri_sph_bess_pot_int(r,1.0_DP/r, &
              swri%sphbessels(species_number, bessel_idx))

          values(val_idx) = bessint * factor * &
               sw_real_sph_harm_theta_part_unit(&
               unit_vec%x, unit_vec%y, unit_vec%z, l,m)

       end do
    end do

    call utils_trace_out(myself)

  end subroutine swri_eval_swpotrt_at_nodes_on_halfdisc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_init_metric_matrix(swri,elements,par)
    !==========================================================================!
    ! This subroutine initialises the structure used to hold all atomblocks of !
    ! the metric matrices.                                                     !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill on 29/10/2010.                                   !
    ! Modified to take par as input by Robert Charlton, September 2018.        !
    !==========================================================================!

    use comms, only: pub_total_num_procs
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO, parallel_strategy_distr_funcs
    use sparse, only: sparse_init_blocking_scheme, BLKS_SW
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout) :: swri
    type(PARAL_INFO), intent(inout) :: par
    type(ELEMENT), intent(in)  :: elements(par%nat)

    ! jd: Local variables
    integer, allocatable :: first_func_on_proc(:)
    integer, allocatable :: num_funcs_on_proc(:)
    integer, allocatable :: first_func_on_atom(:)
    integer, allocatable :: num_funcs_on_atom(:)
    integer, allocatable :: proc_of_func(:)
    integer, allocatable :: atom_of_func(:)
    integer, allocatable :: num_sw(:)
    integer, allocatable :: nfuncs_orig(:)
    integer :: max_funcs_on_proc
    integer :: max_funcs_on_atom
    integer :: ierr      ! Error flag
    integer :: eidx      ! Element idx
    integer :: total_sw  ! Total number of spherical waves
    character(len=*), parameter :: myself = 'swri_init_metric_matrix'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)

    allocate(num_sw(par%num_species), stat=ierr)
    call utils_alloc_check(myself,'num_sw',ierr)
    allocate(nfuncs_orig(par%nat), stat=ierr)
    call utils_alloc_check(myself,'nfuncs_orig',ierr)

    call swri_num_sph_functions(swri%radtable, swri%max_num_bessels, & ! fills
         swri%quality, par, &                                          ! input
         num_sw = num_sw)                                              ! fills

    do eidx = 1, par%nat
       nfuncs_orig(eidx) = num_sw(elements(eidx)%species_number)
    end do
    total_sw = sum(nfuncs_orig)
    par%num_sw = total_sw

    allocate(num_funcs_on_proc(0:pub_total_num_procs-1), stat=ierr)
    call utils_alloc_check(myself,'num_funcs_on_proc',ierr)
    allocate(first_func_on_proc(0:pub_total_num_procs),stat=ierr)
    call utils_alloc_check(myself,'first_func_on_proc',ierr)
    allocate(num_funcs_on_atom(1:par%nat),stat=ierr)
    call utils_alloc_check(myself,'num_funcs_on_atom',ierr)
    allocate(first_func_on_atom(1:par%nat),stat=ierr)
    call utils_alloc_check(myself,'first_func_on_atom',ierr)
    allocate(proc_of_func(1:total_sw),stat=ierr)
    call utils_alloc_check(myself,'proc_of_func',ierr)
    allocate(atom_of_func(1:total_sw),stat=ierr)
    call utils_alloc_check(myself,'atom_of_func',ierr)

    call parallel_strategy_distr_funcs(total_sw, nfuncs_orig, par%orig_atom, &
         par, first_func_on_proc, num_funcs_on_proc, first_func_on_atom, &
         num_funcs_on_atom, proc_of_func, atom_of_func, max_funcs_on_proc, &
         max_funcs_on_atom)

    ! rc2013: pass par explicitly for sparse initialisation
    call sparse_init_blocking_scheme(BLKS_SW,total_sw,num_funcs_on_proc, &
         num_funcs_on_atom, first_func_on_proc, first_func_on_atom, &
         atom_of_func, proc_of_func, par)

    deallocate(num_funcs_on_proc, stat=ierr)
    call utils_dealloc_check(myself,'num_funcs_on_proc',ierr)
    deallocate(first_func_on_proc,stat=ierr)
    call utils_dealloc_check(myself,'first_func_on_proc',ierr)
    deallocate(num_funcs_on_atom,stat=ierr)
    call utils_dealloc_check(myself,'num_funcs_on_atom',ierr)
    deallocate(first_func_on_atom,stat=ierr)
    call utils_dealloc_check(myself,'first_func_on_atom',ierr)
    deallocate(proc_of_func,stat=ierr)
    call utils_dealloc_check(myself,'proc_of_func',ierr)
    deallocate(atom_of_func,stat=ierr)
    call utils_dealloc_check(myself,'atom_of_func',ierr)
    deallocate(num_sw, stat=ierr)
    call utils_dealloc_check(myself,'num_sw',ierr)
    deallocate(nfuncs_orig, stat=ierr)
    call utils_dealloc_check(myself,'nfuncs_orig',ierr)

    call utils_trace_out(myself)

  end subroutine swri_init_metric_matrix

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_extract_matrix_2(outmatrix, swri, &
       num_sws_in_expansion, expansion_centres, expansion_atoms, &
       v_or_o, swexside_filter)
    !==========================================================================!
    ! This subroutine constructs the V or O matrix for up to 2 atoms, returning!
    ! the result in outmatrix.                                                 !
    ! The atoms passed in expansion_atoms, and corresponding centres, must be  !
    ! unique (this is taken care of by swri_expansion_centres).                !
    ! Elements of expansion_atoms() equal to -1 indicate ignored atoms.        !
    ! If there are ignored atoms, the returned matrix will be smaller than 2x2 !
    ! atomblocks.                                                              !
    !--------------------------------------------------------------------------!
    ! This subroutine is unnecessarily complicated, as it evolved from the     !
    ! 4-centre version. Also, it honours Quintin's generalisation of # of SWs  !
    ! per centre depending on centre, at least swriside. It would be good to   !
    ! rewrite this cleanly one day.                                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   outmatrix (out): The constructed V or O matrix (up to 2x2 atomblocks). !
    !   swri (in): The SW_RI containing the full metric matrices from which to !
    !              extract.                                                    !
    !   num_sws_in_expansion (in): Dimensions the output matrix.               !
    !   expansion_centres (in) } These come from                               !
    !   expansion_atoms (in)   } swri_expansion_centres.                       !
    !   v_or_o (in): Specifies whether V or O matrix blocks are extracted.     !
    !   swexside_filter (in): Allows applying a filter to remove some SWs if   !
    !                         a smaller auxiliary basis set is needed SWEXside.!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2012 basing on Quintin Hill's      !
    ! subroutine for a pair of atoms.                                          !
    ! Modified to go back to 2 atoms in expansion.                             !
    ! Modified in February 2015 to support swexside_filter.                    !
    ! Modified in August 2018 to calculate and print eigenvalues.              !
    ! Modified to extract par from SW_RI by Robert Charlton, September 2018.   !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: SW_O, SW_V
    use linalg, only: linalg_dsyev_lt
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_swri_print_eigenvalues
    use sparse, only: sparse_get_block, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort

    implicit none

    ! jd: Arguments
    integer, intent(in)           :: num_sws_in_expansion
    real(kind=DP), intent(out)    :: outmatrix( &
         num_sws_in_expansion, num_sws_in_expansion)
    type(SW_RI), intent(in)       :: swri
    type(ATOM_CENTRE), intent(in) :: expansion_centres(2)
    integer, intent(in)           :: expansion_atoms(2)
    integer, intent(in)           :: v_or_o
    type(SW_QUALITY), intent(in)  :: swexside_filter

    ! jd: Local variables
    real(kind=DP) :: swriside_outmatrix( &
         2*swri%quality%num_sws_per_centre, 2*swri%quality%max_sws_per_centre)
    integer :: offset_to(2) ! jd: Offsets to blocks 1-4 in the matrix
    integer :: last_of(2)   ! jd: Last element of blocks 1-4 in the matrix
    integer :: n_atoms
    integer :: i, j
    ! jd: For debugging purposes only
    real(kind=DP) :: eigenvecs(num_sws_in_expansion, num_sws_in_expansion)
    real(kind=DP) :: eigenvals(num_sws_in_expansion)
    character(len=16) :: formt

    character(len=*), parameter :: myself = 'swri_extract_matrix_2'
    type(PARAL_INFO), pointer  :: par

    ! ------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    outmatrix = 0.0_DP

    ! rc2013: get the parallel strategy from the diagonal metric matrix
    if(swri%has_metric(SW_V)) then
       call sparse_get_par(par, swri%full_metric_matrices(SW_V))
    else if(swri%has_metric(SW_O)) then
       call sparse_get_par(par, swri%full_metric_matrices(SW_O))
    else
       call utils_abort('Error in swri_extract_matrix_2: metric not set, &
            &no parallel strategy to extract.')
    end if

    ! jd: Determine how many atoms out of 2 we have
    do i = 1, 2
       if(expansion_atoms(i) == -1) exit
    end do
    n_atoms = i-1

    call utils_assert(n_atoms > 0 .and. n_atoms <= 2, &
         'Internal error [1] in '//trim(myself))

    ! jd: Determine the swriside offsets in the out V or O matrix for each of
    !     the atoms
    offset_to(1) = 1
    last_of(1) = swri_num_sws_on_centre(swri,expansion_centres(1))
    do i = 2, n_atoms
       offset_to(i) = offset_to(i-1) + &
            swri_num_sws_on_centre(swri,expansion_centres(i-1))
       last_of(i) = offset_to(i) + &
            swri_num_sws_on_centre(swri,expansion_centres(i)) - 1
    end do

    ! jd: Single-centre expansion
    if(n_atoms == 1) then
       call utils_assert(par%proc_of_atom(par%orig_atom(expansion_atoms(1))) ==&
            pub_my_proc_id, 'Internal error [2] in '//trim(myself))
       call sparse_get_block(swriside_outmatrix(&
            offset_to(1):last_of(1), &
            offset_to(1):last_of(1)), &
            swri%full_metric_matrices(v_or_o), &
            expansion_atoms(1), expansion_atoms(1))

       call swri_atomblock_filter_to_swex(outmatrix, &
            swriside_outmatrix(&
            offset_to(1):last_of(1), offset_to(1):last_of(1)), &
            swexside_filter, swri)

       if(all(outmatrix(:,:) == 0.0_DP)) call utils_abort(&
            'Single-centre matrixblock all zero in '//trim(myself),&
            expansion_atoms(1))

       goto 900

    end if

    ! jd: Two-centre expansion
    call utils_assert(par%proc_of_atom(par%orig_atom(expansion_atoms(1))) == &
         pub_my_proc_id .or. &
         par%proc_of_atom(par%orig_atom(expansion_atoms(2))) == &
         pub_my_proc_id, 'Internal error [3] in '//trim(myself))

    do j = 1, n_atoms
       do i = 1, n_atoms
          if(par%proc_of_atom(par%orig_atom(expansion_atoms(j))) == &
               pub_my_proc_id) then
             ! jd: if I own j-th atom, I can get block directly ,regardless of i.
             call sparse_get_block( &
                  swriside_outmatrix(offset_to(i):last_of(i),&
                  offset_to(j):last_of(j)), &
                  swri%full_metric_matrices(v_or_o), &
                  expansion_atoms(i), expansion_atoms(j))
          else
             if(par%proc_of_atom(par%orig_atom(expansion_atoms(i))) == &
                  pub_my_proc_id) then
                ! jd: so I don't own j-th, but if I own i-th, I can use symmetry
                call sparse_get_block( &
                     swriside_outmatrix(offset_to(i):last_of(i),&
                     offset_to(j):last_of(j)), &
                     swri%full_metric_matrices(v_or_o), &
                     expansion_atoms(j), expansion_atoms(i))
                swriside_outmatrix(offset_to(i):last_of(i),&
                     offset_to(j):last_of(j)) = &
                     transpose(swriside_outmatrix(offset_to(i):last_of(i), &
                     offset_to(j):last_of(j)))
             else
                ! jd: I don't own i-th or j-th, so this is a same-centre remote
                !     element. Re-calculate it.
                call utils_assert(i==j,'Internal error [4] in '//trim(myself))
                if(v_or_o == SW_V) then
                   call swri_vmatrixblock_sc(swri, &
                        swriside_outmatrix(offset_to(i):last_of(i),&
                        offset_to(j):last_of(j)), expansion_centres(i))
                else
                   if(v_or_o == SW_O) then
                      call swri_omatrixblock_sc(swri, &
                           swriside_outmatrix(offset_to(i):last_of(i),&
                           offset_to(j):last_of(j)), expansion_centres(i))
                   else
                      call utils_abort('Unrecognized v_or_o in '//trim(myself))
                   end if
                end if
             end if
          end if

          call utils_assert(last_of(j)-offset_to(j) == last_of(i)-offset_to(i),&
               'Internal error [5] in '//trim(myself))

          call swri_atomblock_filter_to_swex(&
               outmatrix(&
               (i-1)*swexside_filter%num_sws_per_centre+1: &
               i*swexside_filter%num_sws_per_centre, &
               (j-1)*swexside_filter%num_sws_per_centre+1: &
               j*swexside_filter%num_sws_per_centre), &
               swriside_outmatrix(&
               offset_to(i):last_of(i), offset_to(j):last_of(j)), &
               swexside_filter, swri)

          if(all(outmatrix(&
               (i-1)*swexside_filter%num_sws_per_centre+1: &
               i*swexside_filter%num_sws_per_centre, &
               (j-1)*swexside_filter%num_sws_per_centre+1: &
               j*swexside_filter%num_sws_per_centre) &
               == 0.0_DP)) call utils_abort(&
               'Cross-centre matrixblock all zero in '//trim(myself), &
               expansion_atoms(i),expansion_atoms(j))
       end do
    end do

900 continue

    ! jd: Calculate the eigenvalues of the returned matrix to get a sense
    !     of linear dependence of auxiliary basis functions, if asked to.
    if(pub_swri_print_eigenvalues) then
       eigenvecs = outmatrix
       call linalg_dsyev_lt(eigenvecs, eigenvals, num_sws_in_expansion)
       write(stdout, '(a, i10, i10)', advance='no') 'SWRI: EIGVALS: ', &
            expansion_atoms(:)
       write(formt, '(a,i0,a)') '(', num_sws_in_expansion, 'e12.4)'
       write(stdout, formt) eigenvals
    end if

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine swri_extract_matrix_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_build_matrix_2(outmatrix, packed_metric_matrix_blocks_ht, &
       swri, num_sws_in_expansion, expansion_centres, expansion_atoms, v_or_o)
    !==========================================================================!
    ! This subroutine constructs the V or O matrix for up to 2 atoms, returning!
    ! the result in outmatrix.                                                 !
    ! The atoms passed in expansion_atoms, and corresponding centres, must be  !
    ! unique (this is taken care of by swri_expansion_centres).                !
    ! Elements of expansion_atoms() equal to -1 indicate ignored atoms.        !
    ! If there are ignored atoms, the returned matrix will be smaller than 2x2 !
    ! atomblocks.                                                              !
    !--------------------------------------------------------------------------!
    ! In contrast to swri_extract_matrix_2(), this subroutine does not use the !
    ! SPAM3 metric matrix. This is because in the new parallelisation scheme   !
    ! we might be interested in matrix blocks that we neither sparse-own, nor  !
    ! sparse-transpose-own. This used to be the case in the old scheme.        !
    ! What we do is we redistribute the necessary atomblocks across MPI ranks  !
    ! and store them in a hash table (packed_metric_matrix_blocks_ht). This    !
    ! subroutine locates them in this hash table.                              !
    !--------------------------------------------------------------------------!
    ! This subroutine is unnecessarily complicated, as it evolved from the     !
    ! 4-centre version. Also, it honours Quintin's generalisation of # of SWs  !
    ! per centre depending on centre, at least swriside. It would be good to   !
    ! rewrite this cleanly one day.                                            !
    !--------------------------------------------------------------------------!
    ! This subroutine does *not* allow for swri-to-swex filtering.             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   outmatrix (out): The constructed V or O matrix (up to 2x2 atomblocks). !
    !   packed_metric_matrix_blocks_ht (in): The HT that contains all my       !
    !                                        atomblocks of full V or O.        !
    !   swri (in): The SW_RI in which we're operating.                         !
    !   num_sws_in_expansion (in): Dimensions the output matrix.               !
    !   expansion_centres (in) } These come from                               !
    !   expansion_atoms (in)   } swri_expansion_centres.                       !
    !   v_or_o (in): Specifies if we want the V or O matrix.                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2019 using swri_extract_matrix_2() as !
    ! a template.                                                              !
    !==========================================================================!

    use constants, only: SW_O, SW_V
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup_nocount
    use linalg, only: linalg_dsyev_lt
    use rundat, only: pub_swri_print_eigenvalues, pub_debug
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort

    implicit none

    ! jd: Arguments
    integer, intent(in)                     :: num_sws_in_expansion
    type(HT_HASH_TABLE), intent(in), target :: packed_metric_matrix_blocks_ht
    real(kind=DP), intent(out)              :: outmatrix( &
         num_sws_in_expansion, num_sws_in_expansion)
    type(SW_RI), intent(in)                 :: swri
    type(ATOM_CENTRE), intent(in)           :: expansion_centres(2)
    integer, intent(in)                     :: expansion_atoms(2)
    integer, intent(in)                     :: v_or_o

    ! jd: Local variables
    integer :: offset_to(2) ! jd: Offsets to blocks 1-4 in the matrix
    integer :: last_of(2)   ! jd: Last element of blocks 1-4 in the matrix
    integer :: n_atoms
    integer :: i, j
    integer :: ndata
    logical :: filled_already
    real(kind=DP) :: outmatrix1D(num_sws_in_expansion*num_sws_in_expansion)
    !                ^ could be four times too big, but this is not costly
    ! jd: For debugging purposes only
    real(kind=DP) :: eigenvecs(num_sws_in_expansion, num_sws_in_expansion)
    real(kind=DP) :: eigenvals(num_sws_in_expansion)
    character(len=16) :: formt

    character(len=*), parameter :: myself = 'swri_build_matrix_2'

    ! ------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    outmatrix = 0.0_DP

    ! jd: Determine how many atoms out of 2 we have
    do i = 1, 2
       if(expansion_atoms(i) == -1) exit
    end do
    n_atoms = i-1

    call utils_assert(n_atoms > 0 .and. n_atoms <= 2, &
         'Internal error [1] in '//trim(myself))

    ! jd: Determine the swriside offsets in the out V or O matrix for each of
    !     the atoms
    offset_to(1) = 1
    last_of(1) = swri_num_sws_on_centre(swri,expansion_centres(1))
    do i = 2, n_atoms
       offset_to(i) = offset_to(i-1) + &
            swri_num_sws_on_centre(swri,expansion_centres(i-1))
       last_of(i) = offset_to(i) + &
            swri_num_sws_on_centre(swri,expansion_centres(i)) - 1
    end do

    ! jd: Single-centre expansion
    if(n_atoms == 1) then
       call hash_table_lookup_nocount(outmatrix1D, ndata, &
            packed_metric_matrix_blocks_ht, &
            expansion_atoms(1), expansion_atoms(1))
       call utils_assert(num_sws_in_expansion * num_sws_in_expansion == &
            ndata, 'Internal error [2] in '//trim(myself), &
            num_sws_in_expansion * num_sws_in_expansion, ndata)
       outmatrix(offset_to(1):last_of(1),offset_to(1):last_of(1)) = &
            reshape(outmatrix1D(1:ndata),&
            (/last_of(1)-offset_to(1)+1,last_of(1)-offset_to(1)+1/))

       if(pub_debug) then
          if(all(outmatrix(:,:) == 0.0_DP)) call utils_abort(&
               'Single-centre matrixblock all zero in '//trim(myself),&
               expansion_atoms(1))
       end if

       goto 900

    end if

    ! jd: Two-centre expansion
    do j = 1, n_atoms
       do i = 1, n_atoms

          call hash_table_lookup_nocount(outmatrix1D, ndata, &
               packed_metric_matrix_blocks_ht, &
               expansion_atoms(i), expansion_atoms(j))
          filled_already = .false.
          if(ndata == -1) then
             ! jd: Not found, look for transpose
             call hash_table_lookup_nocount(outmatrix1D, ndata, &
                  packed_metric_matrix_blocks_ht, &
                  expansion_atoms(j), expansion_atoms(i))
             if(ndata /= -1) then
                outmatrix1D(1:ndata) = reshape(transpose(reshape(&
                     outmatrix1D,(/last_of(i)-offset_to(i)+1,&
                     last_of(j)-offset_to(j)+1/))),(/ndata/))
             else ! jd: Still not found, could be i==j block of 2-atom
                if(expansion_atoms(i) == expansion_atoms(j)) then
                   if(v_or_o == SW_V) then
                      call swri_vmatrixblock_sc(swri, &
                           outmatrix(offset_to(i):last_of(i),&
                           offset_to(j):last_of(j)), expansion_centres(i))
                      filled_already = .true.
                   else
                      if(v_or_o == SW_O) then
                         call swri_omatrixblock_sc(swri, &
                              outmatrix(offset_to(i):last_of(i),&
                              offset_to(j):last_of(j)), expansion_centres(i))
                         filled_already = .true.
                      else
                         call utils_abort('Unrecognized v_or_o in '//trim(myself))
                      end if
                   end if
                else
                   call utils_abort(myself//&
                        ': Internal error [5]. Metric matrix block missing ', &
                        expansion_atoms(i), expansion_atoms(j))
                end if
             end if
          end if
          if(.not. filled_already) then
             call utils_assert(num_sws_in_expansion * num_sws_in_expansion / 4 == &
                  ndata, 'Internal error [3] in '//trim(myself), &
                  num_sws_in_expansion * num_sws_in_expansion / 4 , ndata, &
                  expansion_atoms(i), expansion_atoms(j))

             outmatrix(offset_to(i):last_of(i),offset_to(j):last_of(j)) = &
                  reshape(outmatrix1D(1:ndata),&
                  (/last_of(i)-offset_to(i)+1,last_of(j)-offset_to(j)+1/))
          end if
       end do
    end do

    if(pub_debug) then
       if(all(outmatrix(:,:) == 0.0_DP)) call utils_abort(&
            'Two-centre matrixblock all zero in '//trim(myself),&
            expansion_atoms(1),expansion_atoms(2))
    end if

900 continue

    ! jd: Calculate the eigenvalues of the returned matrix to get a sense
    !     of linear dependence of auxiliary basis functions, if asked to.
    if(pub_swri_print_eigenvalues) then
       eigenvecs = outmatrix
       call linalg_dsyev_lt(eigenvecs, eigenvals, num_sws_in_expansion)
       write(stdout, '(a, i10, i10)', advance='no') 'SWRI: EIGVALS: ', &
            expansion_atoms(:)
       write(formt, '(a,i0,a)') '(', num_sws_in_expansion, 'e12.4)'
       write(stdout, formt) eigenvals
    end if

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine swri_build_matrix_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_make_wmatrix_2(wmatrix, &           ! out
       vmatrix, omatrix, num_sws_in_expansion)        ! in
    !==========================================================================!
    ! This subroutine constructs the W matrix for up to 2 atoms.               !
    ! The W matrix is only used in the "OVO" expansion, which is obsolete.     !
    ! It should still work, except perhaps in HFx gradient.                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   wmatrix (out): The constructed W matrix (see below for definition).    !
    !   vmatrix, omatrix (in): Matrix blocks from which W is constructed.      !
    !   num_sws_in_expansion (in): Used to dimension arrays.                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2012.                                 !
    ! Modified for 2 atoms in July 2012.                                       !
    !==========================================================================!

    use comms, only: pub_on_root
    use linalg, only: linalg_invert_sym_matrix, linalg_dgemm_serial_square
    use rundat, only: pub_swri_overlap_indirect, pub_swri_improve_inverse
    use timer, only: timer_clock

    implicit none

    integer, intent(in)           :: num_sws_in_expansion
    real(kind=DP), intent(out)    :: wmatrix( &
         num_sws_in_expansion, num_sws_in_expansion)
    real(kind=DP), intent(in)     :: vmatrix( &
         num_sws_in_expansion, num_sws_in_expansion)
    real(kind=DP), intent(in)     :: omatrix( &
         num_sws_in_expansion, num_sws_in_expansion)

    ! Local variables
    real(kind=DP) :: work1(num_sws_in_expansion, num_sws_in_expansion)
    real(kind=DP) :: work2(num_sws_in_expansion, num_sws_in_expansion)
    real(kind=DP) :: work3(num_sws_in_expansion, num_sws_in_expansion)

    integer       :: row, col
    real(kind=DP) :: cond

    character(len=*), parameter :: myself = 'swri_make_wmatrix_2'

    ! ------------------------------------------------------------------------

    call utils_trace_in(myself)
    call timer_clock(myself,1)

    ! jd: work1 := V^-1 (iff overlap_indirect T)
    !     work1 := O^-1 (iff overlap_indirect F)
    if(pub_swri_overlap_indirect) then
       work1 = vmatrix
    else
       work1 = omatrix
    end if
    call linalg_invert_sym_matrix(work1,num_sws_in_expansion,cond)
    write(stdout,'(a,e8.2)') 'Condition number for matrix inversion: ',cond

!    call utils_dump_array2D_to_file(work1,'work1')

    if(pub_swri_improve_inverse) then
       if(pub_on_root) then
          write(stdout,'(a)') 'Improving inverse with Hotelling.'
       end if
       ! jd: Let's denote by M:
       !     - the matrix V (iff overlap_indirect T) or
       !     - the matrix O (iff overlap_indirect F).
       !     Now, let's improve the calculated M^-1 by using the relation
       !     better M^1 := 2M^-1 - M^-1.M.M^-1
       !     - First, work2 := M.M^-1 = M.work1
       if(pub_swri_overlap_indirect) then
          call linalg_dgemm_serial_square(work2,vmatrix,work1,&
               num_sws_in_expansion)
       else
          call linalg_dgemm_serial_square(work2,omatrix,work1,&
               num_sws_in_expansion)
       end if
       !     - Then, work3 := M^-1.M.M^-1 = work1.work2
       call linalg_dgemm_serial_square(work3,work1,work2,num_sws_in_expansion)
       !     - Finally, work1 := 2M^-1 - M^-1.M.M^-1 = 2*work1 - work3
       !       And this is is our improved inverse
       do row = 1, num_sws_in_expansion
          do col = 1, num_sws_in_expansion
             work1(col,row) = 2.0_DP * work1(col,row) - work3(col,row)
          end do
       end do
       ! jd: work2 will now be re-used for something else, work3 is not needed

    end if

    ! jd: work2 = V^-1.O (iff overlap_indirect T)
    !     work2 = V.O^-1 (iff overlap_indirect F)
    if(pub_swri_overlap_indirect) then
       call linalg_dgemm_serial_square(work2,work1,omatrix,num_sws_in_expansion)
    else
       call linalg_dgemm_serial_square(work2,vmatrix,work1,num_sws_in_expansion)
    end if

    ! jd: wmatrix = O.V^-1.O    (iff overlap_indirect T)
    !     wmatrix = O^-1.V.O^-1 (iff overlap_indirect F)
    if(pub_swri_overlap_indirect) then
       call linalg_dgemm_serial_square(wmatrix, omatrix, work2, &
            num_sws_in_expansion)
    else
       call linalg_dgemm_serial_square(wmatrix, work1, work2, &
            num_sws_in_expansion)
    end if

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine swri_make_wmatrix_2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function swri_sph_bess_pot_int(rad,invrad,sphbessel)
    !==========================================================================!
    ! This function calculates a spherical Bessel potential integral at a      !
    ! point.                                                                   !
    !--------------------------------------------------------------------------!
    ! Split from internal_swpot_point by Quintin Hill on 29/04/2010.           !
    !==========================================================================!

    implicit none

    real(kind=DP)             :: swri_sph_bess_pot_int
    type(bessel),  intent(in) :: sphbessel ! Spherical Bessel  parameters
    real(kind=DP), intent(in) :: rad       ! Radius from centre
    real(kind=DP), intent(in) :: invrad    ! 1/rad

    real(kind=DP) :: invq
    real(kind=DP) :: q, qr

    if (rad > sphbessel%aval) then
       swri_sph_bess_pot_int = invrad**(sphbessel%lval+1) * &
            sphbessel%farpotint ! jd: (5.6.4, downstairs)
    else

       ! jd: Formulas below verified with Mathematica
       invq = 1.0_DP / sphbessel%qval
       q = sphbessel%qval
       qr = q*rad
       select case(sphbessel%lval)
       case (0)
          swri_sph_bess_pot_int = invrad * invq**3 * sin(qr) &
               + sphbessel%nearpotint
       case(1)
          swri_sph_bess_pot_int = 3.0_DP * invrad * invq**3 * &
               (sin(qr)*invq*invrad - cos(qr)) &
               + rad*sphbessel%nearpotint
       case(2)
          swri_sph_bess_pot_int = invq**5 * invrad**3 *&
               ((15.0_DP - 5.0_DP*qr*qr) * sin(qr) - 15.0_DP*qr*cos(qr))  &
               +rad*rad*sphbessel%nearpotint

       case(3)
          swri_sph_bess_pot_int = invq**6 * invrad**4 * &
               ((105.0_DP-42.0_DP*qr*qr) * sin(qr) &
               - (105.0_DP - 7.0_DP*qr*qr)*qr*cos(qr)) &
               +rad**3 * sphbessel%nearpotint
       case(4)
          swri_sph_bess_pot_int = invq**7 * invrad**5 * 9.0_DP *&
               ((10.0_DP*qr*qr-105.0_DP)*qr*cos(qr) + &
               (qr*qr*(qr*qr-45.0_DP) + 105.0_DP)*sin(qr)) &
               + rad**4 * sphbessel%nearpotint
       case default
          swri_sph_bess_pot_int = 0.0_DP
       end select
    end if

  end function swri_sph_bess_pot_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function swri_sph_bess_solidharm_int(l, q, a) ! in
    !==========================================================================!
    ! This function computes Int(r^(l+2)j_l(qr))dr for l=0..4.                 !
    !--------------------------------------------------------------------------!
    ! Rewritten from scratch with Mathematica's help by Jacek Dziedzic         !
    ! in Jun 2015 and extended to l=4. Noved to swri by Jacek Dziedzic in      !
    ! Oct 2016.                                                                !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    integer, intent(in)         :: l
    real(kind=DP), intent(in)   :: q
    real(kind=DP), intent(in)   :: a

    ! jd: Local variables
    character(len=*), parameter :: myself = 'swri_sph_bess_solidharm_int'

    ! -------------------------------------------------------------------------

    select case(l)
    case(0)
       swri_sph_bess_solidharm_int = &
           sin(q*a)/(q**3.0_DP)-a*cos(q*a)/(q**2.0_DP)
    case(1)
       swri_sph_bess_solidharm_int = -3.0_DP*a*cos(q*a)/(q**3.0_DP)-&
          (q**2.0_DP*a**2.0_DP-3.0_DP)*sin(q*a)/(q**4.0_DP)
    case(2)
       swri_sph_bess_solidharm_int = (((q*a)**3.0_DP-15.0_DP*q*a)*&
            cos(q*a)+(15.0_DP-6.0_DP*(q*a)**2)*sin(q*a))/(q**5.0_DP)
    case(3)
       swri_sph_bess_solidharm_int = (5.0_DP*q*a*(-21.0_DP+2.0_DP*&
            q**2.0_DP*&
            a**2.0_DP)*cos(q*a)+(105.0_DP-45.0_DP*q**2.0_DP*&
            a**2.0_DP+ q**4.0_DP*a**4.0_DP)*sin(q*a))/q**6.0_DP

    case(4)
       swri_sph_bess_solidharm_int = (-a*q*(945.0_DP - 105.0_DP * a*a * q*q + &
            a**4.0_DP * q**4.0_DP) * cos(a*q) + 15.0_DP * (63.0_DP - 28.0_DP * &
            a*a * q*q + a**4.0_DP * q**4.0_DP) * sin(a*q))/q**7.0_DP

    case default
       call utils_abort(trim(myself)//': Unsupported l',l)
    end select

  end function swri_sph_bess_solidharm_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_obtain_swops_in_ppd_set_omp(swri, &                  ! inout
       swops_in_ppd_set, n_sws, &                                      ! out
       expansion_atoms, expansion_centres, ppd_indices_in_set, &       ! in
       n_ppds_in_set, sw_or_swpot, cell, swex_quality)                 ! in
    !==========================================================================!
    ! Not used currently anywhere in the code.                                 !
    !                                                                          !
    ! Obtains the density or potential of all spherical waves coming from      !
    ! src_centre for all points of a ppd. Whether it's the density or the      !
    ! potential that is generated depends on sw_or_swpot.                      !
    !                                                                          !
    ! If possible, the swop is retrieved from swri's hash table cache. In case !
    ! of a cache miss, swri_swop_calc_all_in_ppd() is used to calculate it and !
    ! it is added to the cache. Caching works on unfiltered (swriside) SWOPs.  !
    ! When filtering is in use (swex_quality is lower than swri%quality),      !
    ! we still calculate and cache swriside SWOPs, only filter them for output.!
    !                                                                          !
    ! This version uses OMP internally and HT stashes to improve performance.  !
    ! If filtering is not needed, swri_obtain_swops_in_ppd_ptr_fast() is       !
    ! expected to perform better, provided the SWOP cache is populated in      !
    ! advance.                                                                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swri (inout):       Needed for Bessel bookkeeping and for caching SWOPs!
    !                       behind the scenes.                                 !
    !   swops_in_ppd_set (out): Output buffer, at least                        !
    !                           (n_ppds_in_set * n_pts * n_sws) big.           !
    !   n_sws (out):        The number of SWs written to output.               !
    !                       *** probably broken, since it's PRIVATE and not    !
    !                           OMP-reduced.                                   !
    !   expansion_centres (in) } These come from                               !
    !   expansion_atoms (in)   } swri_expansion_centres.                       !
    !   ppd_indices_in_set (in): Array of PPD for which we want the SWOPs.     !
    !   n_ppds_in_set (in): The number of elements in the above.               !
    !   sw_or_swpot (in):   If 'S' density (SWs) is generated, if 'P', pot.    !
    !   cell (in):          The usual.                                         !
    !   swexside_quality (in): Quality to which the SWOPs are potentially      !
    !                          downfiltered.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2017.                                !
    !==========================================================================!

    use hash_table, only: hash_table_lookup, hash_table_stash_prepare, &
         hash_table_stash, hash_table_stash_commit, hash_table_stash_free
    use rundat, only: pub_ht_stash_size, pub_threads_max
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(SW_RI),       intent(inout)  :: swri
    real(kind=DP),     intent(out)    :: swops_in_ppd_set(:,:,:) ! [*]
    integer,           intent(out)    :: n_sws
    integer,           intent(in)     :: expansion_atoms(2)
    type(ATOM_CENTRE), intent(in)     :: expansion_centres(2)
    integer,           intent(in)     :: ppd_indices_in_set(:)
    integer,           intent(in)     :: n_ppds_in_set
    character,         intent(in)     :: sw_or_swpot
    type(CELL_INFO),   intent(in)     :: cell
    type(SW_QUALITY),  intent(in)     :: swex_quality

    ! jd: Local variables
    integer :: size_of_cached_data
    integer :: num_swri_sws
    logical :: filtering_needed
    integer :: src_offset, dest_offset
    integer :: lval, qidx, sbidx
    integer :: ec, n_centres
    integer :: points_in_bessel
    integer :: ppd_idx, ppd_in_set
    integer :: src_atom
    type(ATOM_CENTRE) :: src_centre
    real(kind=DP)     :: swri_swops_in_ppd_workspace(&
         cell%n_pts * swri%quality%num_sws_per_centre,2,n_ppds_in_set) ! [*]
    real(kind=DP), allocatable :: swops_stash(:)

    ! Fast index:           points in PPD } 1st
    ! Second-fastest index: SW #          }              [*]
    ! Second-slowest:       centre        2nd
    ! Slow index:           PPD #         3rd

    ! -------------------------------------------------------------------------

    filtering_needed = .not. (swri%quality%max_l == swex_quality%max_l .and. &
         swri%quality%max_q == swex_quality%max_q)

    if(expansion_atoms(2) == -1) then
       n_centres = 1
    else
       n_centres = 2
    end if

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP FIRSTPRIVATE(swops_stash) &
!$OMP PRIVATE(ppd_in_set, ppd_idx, src_atom, src_centre, size_of_cached_data, &
!$OMP      num_swri_sws, points_in_bessel, lval, qidx, n_sws, src_offset, &
!$OMP      dest_offset) &
!$OMP SHARED(n_centres, ppd_indices_in_set, n_ppds_in_set, filtering_needed, &
!$OMP      expansion_centres, expansion_atoms, swri, pub_ht_stash_size, &
!$OMP      swri_swops_in_ppd_workspace, sw_or_swpot, swops_in_ppd_set, &
!$OMP      swex_quality, cell)
    ! jd: Have each thread prepare its own stash
    call hash_table_stash_prepare(swops_stash, pub_ht_stash_size, &
         overwrite=.false., overfill_strategy = 'N')

    ! --------------------------------------------------------------------------
    ! jd: Loop over PPDs in set                                             PPD
    ! --------------------------------------------------------------------------
!$OMP DO SCHEDULE(DYNAMIC)
    ppd_loop:                                                                 &
    do ppd_in_set = 1, n_ppds_in_set
       ppd_idx = ppd_indices_in_set(ppd_in_set)

       ! -----------------------------------------------------------------------
       ! jd: Loop over expansion centres                                    EEE
       ! -----------------------------------------------------------------------
       two_centres:                                                           &
       do ec = 1, n_centres
          src_atom = expansion_atoms(ec)
          src_centre = expansion_centres(ec)

          ! jd: Look up SWOPs in PPD in cache.
          !     Owing to stashes, no OMP critical needed
          if(filtering_needed) then
             ! jd: Will have to filter, look up into workspace
             call hash_table_lookup(&
                  swri_swops_in_ppd_workspace(:,ec,ppd_in_set), &
                  size_of_cached_data, swri%swops_in_ppds_ht, &
                  key1 = src_atom, key2 = ppd_idx, key3 = ichar(sw_or_swpot))
          else
             ! jd: No filtering, look up straight into result
             call hash_table_lookup(&
                  swops_in_ppd_set(:,ec,ppd_in_set), &
                  size_of_cached_data, swri%swops_in_ppds_ht, &
                  key1 = src_atom, key2 = ppd_idx, key3 = ichar(sw_or_swpot))
          end if

          if(size_of_cached_data /= -1) then
             ! Cache HIT
             num_swri_sws = size_of_cached_data / cell%n_pts
          else
             ! Cache MISS
             ! - calculate the SWOPs for this PPD. Don't filter at this stage,
             !   as we want to stash the swriside version.
             if(filtering_needed) then
                ! jd: Will have to filter, calculate into workspace
                call swri_swop_calc_all_in_ppd(&
                     swri_swops_in_ppd_workspace(:,ec,ppd_in_set), &      ! out
                     num_swri_sws, &                                      ! out
                     src_centre, ppd_idx, sw_or_swpot, cell, &            ! in
                     swri, swri%quality)                                  ! in
             else
                ! jd: No filtering, calculate straight into result
                call swri_swop_calc_all_in_ppd(&
                     swops_in_ppd_set(:,ec,ppd_in_set), num_swri_sws, &    ! out
                     src_centre, ppd_idx, sw_or_swpot, cell, &             ! out
                     swri, swri%quality)                                   ! in
             end if
             ! - add these to stash

             if(filtering_needed) then
                call hash_table_stash(swops_stash, swri%swops_in_ppds_ht, &
                     swri_swops_in_ppd_workspace(:,ec,ppd_in_set), &
                     num_swri_sws * cell%n_pts, &
                     key1 = src_atom, key2 = ppd_idx, key3 = ichar(sw_or_swpot))
             else
                call hash_table_stash(swops_stash, swri%swops_in_ppds_ht, &
                     swops_in_ppd_set(:,ec,ppd_in_set), &
                     num_swri_sws * cell%n_pts, &
                     key1 = src_atom, key2 = ppd_idx, key3 = ichar(sw_or_swpot))
             end if
          end if

          if(num_swri_sws /= swri%quality%num_sws_per_centre) then
             call utils_abort('swri_obtain_swops_in_ppd_set_omp: Logic error --&
                  & disagreement about num_swri_sws',num_swri_sws, &
                  swri%quality%num_sws_per_centre)
          end if

          ! jd: Filter the SWOPs for this PPD/centre, if needed
          dest_offset = 1
          src_offset = 1
          if(filtering_needed) then
             n_sws = 0
             do sbidx = 1, swri%max_num_bessels
                lval = swri%sphbessels(1,sbidx)%lval ! @IDENTICAL_RADII
                qidx = swri%sphbessels(1,sbidx)%qidx ! @IDENTICAL_RADII
                points_in_bessel = (2*lval+1)*cell%n_pts

                if (lval == -1) exit
                if (lval > swex_quality%max_l) then
                   ! skip to the next Bessel, i.e. by 2l+1 PPDs
                   src_offset = src_offset + points_in_bessel
                   cycle
                end if
                if (qidx > swex_quality%max_q) then
                   ! skip to the next Bessel, i.e. by 2l+1 PPDs
                   src_offset = src_offset + points_in_bessel
                   cycle
                end if

                ! jd: Copy over the next Bessel, i.e. 2l+1 PPDs
                swops_in_ppd_set(&
                     dest_offset:dest_offset+points_in_bessel-1, &
                     ec, ppd_in_set) = &
                     swri_swops_in_ppd_workspace(&
                     src_offset:src_offset+points_in_bessel-1, &
                     ec, ppd_in_set)
                dest_offset = dest_offset + points_in_bessel
                src_offset = src_offset + points_in_bessel
                n_sws = n_sws + 2*lval+1

             end do
          else
             n_sws = num_swri_sws
          end if

       end do two_centres

    end do ppd_loop
!$OMP END DO

    ! jd: Merge each thread's stash into the shared hash table now that
    !     it's not being accessed anymore. Destroy stashes.
!$OMP CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)
    call hash_table_stash_commit(swops_stash, swri%swops_in_ppds_ht)
!$OMP END CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)

    call hash_table_stash_free(swops_stash)

!$OMP END PARALLEL

  end subroutine swri_obtain_swops_in_ppd_set_omp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_obtain_swops_in_ppd(swri, &                            ! inout
       swops_in_ppd, n_sws, &                                            ! out
       src_centre, src_atom, ppd_idx, sw_or_swpot, cell, swex_quality, & ! in
       which_images)
    !==========================================================================!
    ! Obtains the density or potential of all spherical waves coming from      !
    ! src_centre for all points of a ppd. Whether it's the density or the      !
    ! potential that is generated depends on sw_or_swpot.                      !
    !                                                                          !
    ! If possible, the swop is retrieved from swri's hash table cache. In case !
    ! of a cache miss, swri_swop_calc_all_in_ppd() is used to calculate it and !
    ! it is added to the cache. Caching works on unfiltered (swriside) SWOPs.  !
    ! When filtering is in use (swex_quality is lower than swri%quality),      !
    ! we still calculate and cache swriside SWOPs, only filter them for output.!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swri (inout):       Needed for Bessel bookkeeping and for caching SWOPs!
    !                       behind the scenes.                                 !
    !   swops_in_ppd (out): Output buffer, at least (n_pts * n_sws) big.       !
    !   n_sws (out):        The number of SWs written to output.               !
    !   src_centre (in):    Centre generating the SW densities/potentials.     !
    !   src_atom (in):      Atom generating the SW densities/potentials.       !
    !   ppd_idx (in):       The index of the PPD, needed to find its position. !
    !   sw_or_swpot (in):   If 'S' density (SWs) is generated, if 'P', pot.    !
    !   cell (in):          The usual.                                         !
    !   swexside_quality (in): Quality to which the SWOPs are potentially      !
    !                          downfiltered.                                   !
    !   which_images (out, opt): Debugging facility.
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2016.                               !
    !==========================================================================!

    use hash_table, only: hash_table_lookup, hash_table_add
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(SW_RI),       intent(inout)  :: swri
    real(kind=DP),     intent(out)    :: swops_in_ppd(:)
    integer,           intent(out)    :: n_sws
    type(ATOM_CENTRE), intent(in)     :: src_centre
    integer,           intent(in)     :: src_atom
    integer,           intent(in)     :: ppd_idx
    character,         intent(in)     :: sw_or_swpot
    type(CELL_INFO),   intent(in)     :: cell
    type(SW_QUALITY),  intent(in)     :: swex_quality
    type(PBC_IMAGE_INFO), intent(out), optional :: which_images(cell%n_pts)

    ! jd: Local variables
    integer :: size_of_cached_data
    integer :: num_swri_sws
    logical :: filtering_needed
    integer :: src_offset, dest_offset
    integer :: lval, qidx, sbidx
    integer :: points_in_bessel
    real(kind=DP) :: swri_swops_in_ppd_workspace(cell%n_pts * &
         swri%quality%num_sws_per_centre)

    ! -------------------------------------------------------------------------

    ! *** Concurrency ***
    ! We are called from OMP regions.
    ! Two issues of importance here (both addressed).
    ! 1) hash_table_adds cannot happen concurrently, hence the critical.
    ! 2) hash_table_adds would normally invalidate all pointers
    !    (swops_in_ppd here) and indeed, it breaks subtly when suffi-
    !    ciently many threads are used and the table is made to cleanup
    !    often. We prevent this by using 'F' overfill strategy and performing
    !    any cleanup after the OMP region. This might overfill the ht
    !    just a tiny bit. The caller then calls hash_table_cleanup_at_will()
    !    outside the OMP loop.
    !
    ! *** Pointers ***
    ! Before this was moved to swri, we had two copies of this -- one in
    ! hf_exchange, and another in sw_expansion. They looked up by pointer,
    ! or, in case of cache misses, pointed to a local workspace. The workspace
    ! was allocated in the caller. Now with swexside filtering done on the fly,
    ! this is not that easy -- there is no hash table element to point to, as
    ! the hash table contains swriside SWOPs. We could return by pointer if
    ! there is no filtering, but workspace would have to be allocated in
    ! caller -- the workspace here is volatile. @optimizeme@

    filtering_needed = .not. (swri%quality%min_l == swex_quality%min_l .and. &
         swri%quality%max_l == swex_quality%max_l .and. &
         swri%quality%max_q == swex_quality%max_q)

!$OMP CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)
    if(filtering_needed) then
       ! jd: Will have to filter, look up into workspace
       call hash_table_lookup(swri_swops_in_ppd_workspace, size_of_cached_data,&
            swri%swops_in_ppds_ht, key1 = src_atom, key2 = ppd_idx, &
            key3 = ichar(sw_or_swpot))
    else
       ! jd: No filtering, look up straight into result
       call hash_table_lookup(swops_in_ppd, size_of_cached_data,&
            swri%swops_in_ppds_ht, key1 = src_atom, key2 = ppd_idx, &
            key3 = ichar(sw_or_swpot))
    end if
!$OMP END CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)

    if(size_of_cached_data /= -1) then
       ! Cache HIT
       num_swri_sws = size_of_cached_data / cell%n_pts
    else
       ! Cache MISS
       ! - calculate the SWOPs for this PPD. Don't filter at this stage, as
       !   we want to cache the swriside version.
       if(filtering_needed) then
          ! jd: Will have to filter, calculate into workspace
          call swri_swop_calc_all_in_ppd(swri_swops_in_ppd_workspace, &    ! out
               num_swri_sws, &                                             ! out
               src_centre, ppd_idx, sw_or_swpot, cell, swri, swri%quality, &! in
               which_images)
       else
          ! jd: No filtering, calculate straight into result
          call swri_swop_calc_all_in_ppd(swops_in_ppd, num_swri_sws, &     ! out
               src_centre, ppd_idx, sw_or_swpot, cell, swri, swri%quality, &! in
               which_images)
       end if
       ! - add these to the cache
!$OMP CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)
       if(filtering_needed) then
          call hash_table_add(swri%swops_in_ppds_ht, &
               swri_swops_in_ppd_workspace, &
               num_swri_sws * cell%n_pts, key1 = src_atom, key2 = ppd_idx, &
               key3 = ichar(sw_or_swpot), overfill_strategy = 'F')
       else
          call hash_table_add(swri%swops_in_ppds_ht, swops_in_ppd, &
               num_swri_sws * cell%n_pts, key1 = src_atom, key2 = ppd_idx, &
               key3 = ichar(sw_or_swpot), overfill_strategy = 'F')
       end if
!$OMP END CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)
    end if

    if(num_swri_sws /= swri%quality%num_sws_per_centre) then
       call utils_abort('swri_obtain_swops_in_ppd: Logic error -- disagree&
            &ment about num_swri_sws',num_swri_sws, &
            swri%quality%num_sws_per_centre)
    end if

    ! jd: Filter the obtained SWOPs, if needed
    dest_offset = 1
    src_offset = 1
    if(filtering_needed) then
       n_sws = 0
       do sbidx = 1, swri%max_num_bessels
          lval = swri%sphbessels(1,sbidx)%lval ! @IDENTICAL_RADII
          qidx = swri%sphbessels(1,sbidx)%qidx ! @IDENTICAL_RADII
          points_in_bessel = (2*lval+1)*cell%n_pts

          if (lval == -1) exit
          if (lval < swex_quality%min_l .or. lval > swex_quality%max_l) then
             ! skip to the next Bessel, i.e. by 2l+1 PPDs
             src_offset = src_offset + points_in_bessel
             cycle
          end if
          if (qidx > swex_quality%max_q) then
             ! skip to the next Bessel, i.e. by 2l+1 PPDs
             src_offset = src_offset + points_in_bessel
             cycle
          end if

          ! jd: Copy over the next Bessel, i.e. 2l+1 PPDs
          swops_in_ppd(dest_offset:dest_offset+points_in_bessel-1) = &
               swri_swops_in_ppd_workspace(src_offset:src_offset+points_in_bessel-1)
          dest_offset = dest_offset + points_in_bessel
          src_offset = src_offset + points_in_bessel
          n_sws = n_sws + 2*lval+1

       end do
    else
       n_sws = num_swri_sws
    end if

  end subroutine swri_obtain_swops_in_ppd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_obtain_swops_in_ppd_ptr(swri, &                      ! inout
       swops_in_ppd_ptr, n_sws, &                                      ! out
       src_centre, src_atom, ppd_idx, sw_or_swpot, cell, swex_quality,&! in
       dont_lookup, use_fast)                                          ! in opt
    !==========================================================================!
    ! Obtains the density or potential of all spherical waves coming from      !
    ! src_centre for all points of a ppd. Whether it's the density or the      !
    ! potential that is generated depends on sw_or_swpot.                      !
    !                                                                          !
    ! The result is returned by pointer to improve efficiency. The pointer     !
    ! points to the internals of the hash table used to cache the SWOPs.       !
    ! For this to happen two conditions must be satisfied:                     !
    ! - no filtering (as we store unfiltered SWOPs),                           !
    ! - the SWOP must first be added to the cache and must fit.                !
    !                                                                          !
    ! The first condition is checked and we abort if it is not satisfied.      !
    ! The second condition is taken care of by using overfill_strategy = 'F'.  !
    !                                                                          !
    ! If possible, the swop is retrieved from swri's hash table cache. In case !
    ! of a cache miss, swri_swop_calc_all_in_ppd() is used to calculate it and !
    ! it is added to the cache. Caching works on unfiltered (swriside) SWOPs.  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swri (inout):       Needed for Bessel bookkeeping and for caching SWOPs!
    !                       behind the scenes.                                 !
    !   swops_in_ppd_ptr (out): Pointer to result, at least (n_pts*n_sws) big. !
    !   n_sws (out):        The number of SWs in the buffer pointed to.        !
    !   src_centre (in):    Centre generating the SW densities/potentials.     !
    !   src_atom (in):      Atom generating the SW densities/potentials.       !
    !   ppd_idx (in):       The index of the PPD, needed to find its position. !
    !   sw_or_swpot (in):   If 'S' density (SWs) is generated, if 'P', pot.    !
    !   cell (in):          The usual.                                         !
    !   swexside_quality (in): Expected quality of SWOPs. If different from the!
    !                          quality in SWRI, we abort.                      !
    !   dont_lookup (in, opt): If present and .true., SWOPs will not be looked !
    !                          up in the cache and assumed absent.             !
    !   use_fast (in, opt): If present and .true., the subroutine              !
    !                       swri_swop_calc_all_in_ppd_fast() will be used for  !
    !                       calculating SWOPs.                                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2019.                                 !
    !==========================================================================!

    use constants, only: garbage_int
    use hash_table, only: hash_table_lookup_ptr, hash_table_add
    use rundat, only: pub_hfx_bc_is_periodic, pub_hfx_debug
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(SW_RI),       intent(inout)  :: swri
    real(kind=DP),     intent(out), pointer :: swops_in_ppd_ptr(:)
    integer,           intent(out)    :: n_sws
    type(ATOM_CENTRE), intent(in)     :: src_centre
    integer,           intent(in)     :: src_atom
    integer,           intent(in)     :: ppd_idx
    character,         intent(in)     :: sw_or_swpot
    type(CELL_INFO),   intent(in)     :: cell
    type(SW_QUALITY),  intent(in)     :: swex_quality
    logical, optional, intent(in)     :: dont_lookup
    logical, optional, intent(in)     :: use_fast

    ! jd: Local variables
    integer :: size_of_cached_data
    integer :: num_swri_sws
    integer, save :: sample = 0
!$OMP THREADPRIVATE(sample)
    logical :: filtering_needed
    logical :: loc_dont_lookup
    logical :: loc_use_fast
    logical :: pbc
    real(kind=DP) :: swri_swops_in_ppd_workspace(cell%n_pts * &
         swri%quality%num_sws_per_centre)
    type(PBC_IMAGE_INFO) :: which_images(cell%n_pts)
    integer :: ii, ipt
    integer, parameter :: sample_period = 100
    character(len=*), parameter :: myself = 'swri_obtain_swops_in_ppd_ptr'

    ! -------------------------------------------------------------------------

    ! *** Concurrency ***
    ! We are called from OMP regions.
    ! Two issues of importance here (both addressed).
    ! 1) hash_table_adds cannot happen concurrently, hence the critical.
    ! 2) hash_table_adds would normally invalidate all pointers
    !    (swops_in_ppd here) and indeed, it breaks subtly when suffi-
    !    ciently many threads are used and the table is made to cleanup
    !    often. We prevent this by using 'F' overfill strategy and performing
    !    any cleanup after the OMP region. This might overfill the ht
    !    just a tiny bit. The caller then calls hash_table_cleanup_at_will()
    !    outside the OMP loop.
    !

    filtering_needed = .not. (swri%quality%min_l == swex_quality%min_l .and. &
         swri%quality%max_l == swex_quality%max_l .and. &
         swri%quality%max_q == swex_quality%max_q)

    if(filtering_needed) then
       call utils_abort(myself//' does not support swexside filtering.')
    end if

    if(present(dont_lookup)) then
       loc_dont_lookup = dont_lookup
    else
       loc_dont_lookup = .false.
    end if

    if(present(use_fast)) then
       loc_use_fast = use_fast
    else
       loc_use_fast = .false.
    end if

    if(loc_dont_lookup) then
       size_of_cached_data = -1
    else
!$OMP CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)
       call hash_table_lookup_ptr(swops_in_ppd_ptr, size_of_cached_data,&
            swri%swops_in_ppds_ht, key1 = src_atom, key2 = ppd_idx, &
            key3 = ichar(sw_or_swpot))
!$OMP END CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)
    end if

    if(size_of_cached_data /= -1) then
       ! Cache HIT
       num_swri_sws = size_of_cached_data / cell%n_pts
    else
       ! Cache MISS
       ! - calculate the SWOPs for this PPD.
       if(loc_use_fast) then

          if(.not. any(pub_hfx_bc_is_periodic)) then
             pbc = .false.
          else if(all(pub_hfx_bc_is_periodic)) then
             pbc = .true.
          else
             call utils_abort(myself//': Mixed BCs are not currently supported.')
          end if

          if(.not. pub_hfx_debug) then
             call swri_swop_calc_all_in_ppd_fast(swri_swops_in_ppd_workspace,&! out
                  num_swri_sws, &                                             ! out
                  src_centre, ppd_idx, sw_or_swpot, cell, swri, pbc, &        ! in
                  sample, sample_period)                                      ! in
          else
             call swri_swop_calc_all_in_ppd_fast(swri_swops_in_ppd_workspace,&! out
                  num_swri_sws, &                                             ! out
                  src_centre, ppd_idx, sw_or_swpot, cell, swri, pbc, &        ! in
                  sample, sample_period, which_images)                        ! in

             do ipt = 1, cell%n_pts
               if(swri%image_map(ipt,ppd_idx,src_centre%global_idx) == garbage_int) then
                  call utils_abort(myself//': The impossible happened: ', &
                       swri%image_map(ipt,ppd_idx,src_centre%global_idx), &
                       ipt, ppd_idx, src_centre%global_idx)
               end if
               swri%image_map(ipt,ppd_idx,src_centre%global_idx) = which_images(ipt)%image
             end do
          end if
          sample = sample + 1
       else
          call swri_swop_calc_all_in_ppd(swri_swops_in_ppd_workspace, &    ! out
               num_swri_sws, &                                             ! out
               src_centre, ppd_idx, sw_or_swpot, cell, swri, swri%quality) ! in
       end if
       ! - add these to the cache
!$OMP CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)
       call hash_table_add(swri%swops_in_ppds_ht, &
            swri_swops_in_ppd_workspace, &
            num_swri_sws * cell%n_pts, key1 = src_atom, key2 = ppd_idx, &
            key3 = ichar(sw_or_swpot), overfill_strategy = 'F')
!$OMP END CRITICAL(SEC_SWOPS_IN_PPDS_HT_ACCESS)
       ! - return pointer to HT data, not the temp workspace
       call hash_table_lookup_ptr(swops_in_ppd_ptr, size_of_cached_data,&
            swri%swops_in_ppds_ht, key1 = src_atom, key2 = ppd_idx, &
            key3 = ichar(sw_or_swpot))
       num_swri_sws = size_of_cached_data / cell%n_pts
    end if

    if(num_swri_sws /= swri%quality%num_sws_per_centre) then
       call utils_abort(myself//': Logic error -- disagreement about &
            &num_swri_sws', num_swri_sws, swri%quality%num_sws_per_centre)
    end if

    n_sws = num_swri_sws

  end subroutine swri_obtain_swops_in_ppd_ptr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_obtain_swops_in_ppd_ptr_fast(swops_in_ppd_ptr, n_sws,& ! out
       src_centre, src_atom, ppd_idx, sw_or_swpot, cell, swri, &         ! in
       swops_in_ppd_workspace, which_images)                             ! inout
    !==========================================================================!
    ! Obtains the density or potential of all spherical waves coming from      !
    ! src_centre for all points of a PPD. Whether it's the density or the      !
    ! potential that is generated depends on sw_or_swpot.                      !
    !                                                                          !
    ! This version is only used in HFx. DMA and polemb use non-fast versions.  !
    ! This version does not support swri-to-swex filtering.                    !
    ! This version relies on the SWOP cache, but only reads from it, it does   !
    ! not update or clean up the SWOP cache. This allows us to forgo CRITICAL  !
    ! regions.                                                                 !
    !                                                                          !
    ! The result is returned by pointer (swops_in_ppds_ptr) to improve         !
    ! efficiency. If the SWOPs are found in the cache, the pointer will point  !
    ! to the internals of the hash table used to cache the SWOPs. If the SWOPs !
    ! are *not* found in the cache, things are more tricky, since we cannot    !
    ! and do not want to add them here. In this case the data is returned in   !
    ! swops_in_ppd_workspace (provided by caller) and the pointer points there.!
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swops_in_ppd_ptr (out): Pointer to result, at least (n_pts*n_sws) big. !
    !   n_sws (out):        The number of SWs in the buffer pointed to.        !
    !   src_centre (in):    Centre generating the SW densities/potentials.     !
    !   src_atom (in):      Atom generating the SW densities/potentials.       !
    !   ppd_idx (in):       The index of the PPD, needed to find its position. !
    !   sw_or_swpot (in):   If 'S' density (SWs) is generated, if 'P', pot.    !
    !   cell (in):          The usual.                                         !
    !   swri (in):          Needed for Bessel bookkeeping and for caching SWOPs!
    !                       behind the scenes.                                 !
    !   swops_in_ppd_workspace (inout): Workspace provided by caller, at least !
    !                                   n_pts * swri%quality%num_sws_per_centre!
    !                                   elements big.                          !
    !   which_images (out, opt): Debugging facility.
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2019.                                !
    !==========================================================================!

    use constants, only: garbage_real
    use geometry, only: POINT
    use hash_table, only: hash_table_lookup_ptr_nocount
    use rundat, only: pub_hfx_bc_is_periodic, pub_hfx_debug
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(CELL_INFO), intent(in)          :: cell
    type(SW_RI), intent(in)              :: swri
    real(kind=DP), intent(out), pointer  :: swops_in_ppd_ptr(:)
    integer, intent(out)                 :: n_sws
    type(ATOM_CENTRE), intent(in)        :: src_centre
    integer, intent(in)                  :: src_atom
    integer, intent(in)                  :: ppd_idx
    character, intent(in)                :: sw_or_swpot
    real(kind=DP), intent(inout), target :: swops_in_ppd_workspace(&
         cell%n_pts * swri%quality%num_sws_per_centre)
    type(PBC_IMAGE_INFO), intent(out), optional  :: which_images(cell%n_pts)

    ! jd: Local variables
    integer :: size_of_cached_data
    logical :: pbc
    integer, save :: sample = 0
!$OMP THREADPRIVATE(sample)
    integer, parameter :: sample_period = 100
    character(len=*), parameter :: myself = 'swri_obtain_swops_in_ppd_ptr_fast'
    ! -------------------------------------------------------------------------

    ! jd: No critical needed, because no adds happen concurrently.
    call hash_table_lookup_ptr_nocount(swops_in_ppd_ptr, size_of_cached_data,&
         swri%swops_in_ppds_ht, key1 = src_atom, key2 = ppd_idx, &
         key3 = ichar(sw_or_swpot))

    if(size_of_cached_data /= -1) then
       ! Cache HIT
       n_sws = size_of_cached_data / cell%n_pts
       if(pub_hfx_debug) then
          if(present(which_images)) then
             ! jd: For cached SWOPs we only have image indices from the image_map,
             !     but not the other quantities.
             which_images(:)%image = swri%image_map(:,ppd_idx,src_centre%global_idx)
             which_images(:)%centre = POINT(garbage_real,garbage_real,garbage_real)
             which_images(:)%curpoint = which_images(:)%centre
             which_images(:)%disp = which_images(:)%centre
             which_images(:)%disp0 = which_images(:)%centre
             which_images(:)%rad = garbage_real
             which_images(:)%rad0 = garbage_real
          else
             call utils_abort(myself//&
                  ': Optional parameter is required when hfx_debug is T.')
          end if
       end if
    else

       if(.not. any(pub_hfx_bc_is_periodic)) then
          pbc = .false.
       else if(all(pub_hfx_bc_is_periodic)) then
          pbc = .true.
       else
          call utils_abort(myself//': Mixed BCs are not currently supported.')
       end if

       ! Cache MISS
       ! - calculate the SWOPs for this PPD, goes into workspace
       call swri_swop_calc_all_in_ppd_fast(swops_in_ppd_workspace, n_sws, & ! out
            src_centre, ppd_idx, sw_or_swpot, cell, swri, pbc, &            ! in
            sample, sample_period, which_images)                            ! in
       sample = sample + 1
       n_sws = swri%quality%num_sws_per_centre
       ! - return pointer to the temp workspace, which caller supplied
       swops_in_ppd_ptr => swops_in_ppd_workspace
    end if

    if(n_sws /= swri%quality%num_sws_per_centre) then
       call utils_abort(myself//': Logic error -- disagreement about &
            &n_sws', n_sws, swri%quality%num_sws_per_centre)
    end if

  end subroutine swri_obtain_swops_in_ppd_ptr_fast

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_swop_calc_all_in_ppd(swops_in_ppd, n_sws, &           ! out
       src_centre, ppd_idx, sw_or_swpot, cell, swri, swexside_filter, & ! in
       which_images)
    !==========================================================================!
    ! Generates the density or potential of all spherical waves coming from    !
    ! src_centre for all points of a ppd. Whether it's the density or the      !
    ! potential that is generated depends on sw_or_swpot.                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swops_in_ppd (out): Output buffer, at least (n_pts * n_sws) big.       !
    !   n_sws (out):        The number of SWs written to output.               !
    !   src_centre (in):    Centre generating the SW densities/potentials.     !
    !   ppd_idx (in):       The index of the PPD, needed to find its position. !
    !   sw_or_swpot (in):   If 'S' density (SWs) is generated, if 'P', pot.    !
    !   cell (in):          The usual.                                         !
    !   swri (in):          Needed to dimension a local buffer and for Bessel  !
    !                       book-keeping in the internal.                      !
    !   swexside_filter (in): Used in the internal to figure out n_sws.        !
    !   which_images (out, opt): Debugging facility.
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2012 using pot_all_in_tb as a         !
    ! template. Rewritten in June 2012 to use PPDs. Updated in February 2015   !
    ! to use swri and to employ filter.                                        !
    !==========================================================================!

    use basis, only: basis_ppd_location
    use geometry, only: POINT, operator(*), operator(-), operator(+)
    use rundat, only: pub_swri_swop_smoothing, pub_hfx_bc_is_periodic
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP),     intent(out) :: swops_in_ppd(:)
    integer,           intent(out) :: n_sws
    type(ATOM_CENTRE), intent(in)  :: src_centre
    integer,           intent(in)  :: ppd_idx
    character,         intent(in)  :: sw_or_swpot
    type(CELL_INFO),   intent(in)  :: cell
    type(SW_RI),       intent(in)  :: swri
    type(SW_QUALITY),  intent(in)  :: swexside_filter
    type(PBC_IMAGE_INFO), intent(out), optional :: which_images(cell%n_pts)

    ! jd: Local variables
    character(len=*), parameter :: myself = 'swri_swop_calc_all_in_ppd'
    type(POINT)   :: dest_corner1
    type(POINT)   :: curpoint
    type(POINT)   :: a1, a2, a3
    integer       :: d1idx, d2idx, d3idx
    integer       :: offs_to_point
    integer       :: sw
    real(kind=DP) :: buffer(swri%quality%max_sws_per_centre)
    logical       :: pbc

    ! jd: Variables and parameters required for SWOP smoothing
    real(kind=DP) :: buffer2(swri%quality%max_sws_per_centre)
    integer, parameter :: n_stencil_pts = 27
    integer       :: stencil_pt
    integer       :: i_pt_in_ppd
    real(kind=DP), parameter :: C0 = 1.0_DP/8.0_DP
    real(kind=DP), parameter :: C1 = 0.5_DP/8.0_DP
    real(kind=DP), parameter :: C2 = 0.25_DP/8.0_DP
    real(kind=DP), parameter :: C3 = 0.125_DP/8.0_DP
    real(kind=DP) :: stencil(n_stencil_pts,4)
    DATA stencil(:,1) / 0D0, &
        0.5_DP, -0.5_DP, 0D0, 0D0, 0D0, 0D0, &
        0.5_DP, 0.5_DP, -0.5_DP, -0.5_DP, &
        0.5_DP, 0.5_DP, -0.5_DP, -0.5_DP, &
        0D0, 0D0, 0D0, 0D0, &
        0.5_DP, 0.5_DP, 0.5_DP, 0.5_DP, &
        -0.5_DP, -0.5_DP, -0.5_DP, -0.5_DP /
    DATA stencil(:,2) / 0D0, &
        0D0, 0D0, 0.5_DP, -0.5_DP, 0D0, 0D0, &
        0.5_DP, -0.5_DP, 0.5_DP, -0.5_DP, &
        0D0, 0D0, 0D0, 0D0, &
        0.5_DP, 0.5_DP, -0.5_DP, -0.5_DP, &
        0.5_DP, 0.5_DP, -0.5_DP, -0.5_DP, &
        0.5_DP, 0.5_DP, -0.5_DP, -0.5_DP /
    DATA stencil(:,3) / 0D0, &
        0D0, 0D0, 0D0, 0D0, 0.5_DP, -0.5_DP, &
        0D0, 0D0, 0D0, 0D0, &
        0.5_DP, -0.5_DP, 0.5_DP, -0.5_DP, &
        0.5_DP, -0.5_DP, 0.5_DP, -0.5_DP, &
        0.5_DP, -0.5_DP, 0.5_DP, -0.5_DP, &
        0.5_DP, -0.5_DP, 0.5_DP, -0.5_DP /
    DATA stencil(:,4) / C0, &
        C1, C1, C1, C1, C1, C1, &
        C2, C2, C2, C2, &
        C2, C2, C2, C2, &
        C2, C2, C2, C2, &
        C3, C3, C3, C3, &
        C3, C3, C3, C3 /

    ! ------------------------------------------------------------------------

    call timer_clock(myself,1)
    ! don't utils_trace_in/out this subroutine -- it's called from OMP contexts

!    n_sws = swri_num_sws_on_centre(swri,src_centre)
!    -> now uses swexside, set by internal_*

    dest_corner1 = basis_ppd_location(ppd_idx,cell)

    a1 = cell%d1 * cell%a1_unit
    a2 = cell%d2 * cell%a2_unit
    a3 = cell%d3 * cell%a3_unit

    curpoint = dest_corner1 - a1 - a2 -a3
    offs_to_point = 1

    if(.not. any(pub_hfx_bc_is_periodic)) then
       pbc = .false.
    else if(all(pub_hfx_bc_is_periodic)) then
       pbc = .true.
    else
       call utils_abort(myself//': Mixed BCs are not currently supported.')
    end if

    ! NB: The ordering of swops_in_ppd wrt cache-friendliness is crucial here.
    !     It's much faster to calculate all SWs for a point at once, hence
    !     internal_* calculates all SWs at once. However, later on it's crucial
    !     that the data is organized with the SW being the major index and the
    !     PPD point being the minor index. Thus the cache-unfriendly loop over
    !     SWs here that re-orders the layout. It's much, much cheaper to pay
    !     the cache penalty here, when generating the SWs which are going to be
    !     cached a moment later than to pay it many times when accessing the
    !     data fetched from the cache.
    i_pt_in_ppd = 1
    do d3idx = 1, cell%n_pt3
       curpoint = curpoint + a3
       do d2idx = 1, cell%n_pt2
          curpoint = curpoint + a2
          do d1idx = 1, cell%n_pt1
             curpoint = curpoint + a1

             if(sw_or_swpot == 'S') then
                if(.not. pub_swri_swop_smoothing) then
                   call internal_sws_at_point(buffer, &
                        curpoint, src_centre, pbc) ! also sets n_sws to swexside
                else
                   call utils_abort(myself//': swri_swop_smoothing not &
                        &implemented for overlap metric.')
                end if
             else
                ! jd: If no smoothing, just calculate the SWpots
                if(.not. pub_swri_swop_smoothing) then
                   ! jd: which_images is an optional, so we can pass it around
                   !     if absent, but not index things out of it. Hence the if.
                   if(present(which_images)) then
                      call internal_swpots_at_point(buffer, &
                           curpoint, src_centre, pbc, which_images(i_pt_in_ppd))
                           ! also sets n_sws to swexside
                   else
                      call internal_swpots_at_point(buffer, &
                           curpoint, src_centre, pbc)
                           ! also sets n_sws to swexside
                   end if
                else
                   ! jd: With smooting use a 27-point stencil to include
                   !     neighbouring points in a weighted way
                   buffer = 0.0_DP
                   do stencil_pt = 1, n_stencil_pts
                      call internal_swpots_at_point(buffer2, &
                          curpoint + & ! @MIC
                          stencil(stencil_pt,1)*a1 + &
                          stencil(stencil_pt,2)*a2 + &
                          stencil(stencil_pt,3)*a3,  &
                          src_centre, pbc, which_images(i_pt_in_ppd))
                      buffer = buffer + stencil(stencil_pt,4) * buffer2
                   end do
                end if
             end if

             do sw = 1, n_sws
                swops_in_ppd((sw-1)*cell%n_pts + offs_to_point) = buffer(sw)
             end do

             offs_to_point = offs_to_point+1
             i_pt_in_ppd = i_pt_in_ppd + 1

          end do
          curpoint = curpoint &
               - real(cell%n_pt1, kind=DP) * a1
       end do
       curpoint = curpoint &
            - real(cell%n_pt2,kind=DP) * a2
    end do

    call timer_clock(myself,2)

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_sws_at_point(sw_val, curpoint, centrex, pbc)

      ! Updated in February 2015 by Jacek Dziedzic to use swex filter.
      ! Updated in 2016 by Jacek Dziedzic, filter moved elsewhere.

      use constants, only: SAFE_DIV_EPS, SQRT_PI
      use geometry, only: POINT, operator(*), operator(-), geometry_MAGNITUDE
      use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm_unit
      use simulation_cell, only: minimum_image_displacement

      implicit none

      ! jd: Arguments
      real(kind=DP),     intent(out) :: sw_val(swri%quality%max_sws_per_centre)
      type(POINT),       intent(in)  :: curpoint
      type(ATOM_CENTRE), intent(in)  :: centrex
      logical,           intent(in)  :: pbc

      ! jd: Local variables
      real(kind=DP) :: rad
      integer       :: lval
      integer       :: qidx
      integer       :: mval
      real(kind=DP) :: qval
      integer       :: sbidx
      real(kind=DP) :: sphbess
      type(POINT)   :: disp, unit_disp
      integer       :: sw_idx

      !@simplify loops for the case of 0

      ! ----------------------------------------------------------------------

      sw_val = 0.0_DP

      if(pbc) then
         ! jd: PBC version
         call minimum_image_displacement(rad, disp, curpoint, centrex%incell, cell)
      else
         ! jd: OBC version
         disp = curpoint - centrex%incell
         rad = geometry_MAGNITUDE(disp)
      end if

      if(rad > SAFE_DIV_EPS) then
         unit_disp = 1.0_DP / rad * disp
      else
         ! jd: Corner case of r=0. All SWs are 0, except for l=0, m=0.
         sw_idx = 1
         do sbidx = 1, swri%max_num_bessels
            lval = swri%sphbessels(centrex%species_number,sbidx)%lval
            qidx = swri%sphbessels(centrex%species_number,sbidx)%qidx

            if (lval == -1) exit
            if (lval < swexside_filter%min_l) cycle
            if (lval > swexside_filter%max_l) cycle
            if (qidx > swexside_filter%max_q) cycle

            do mval=-lval,lval
               sw_val(sw_idx) = merge(1.0_DP/(2.0_DP * SQRT_PI), 0.0_DP, &
                    lval == 0 .and. mval == 0)
               sw_idx = sw_idx + 1
            end do

         end do

         n_sws = sw_idx - 1 ! (in parent)

         return

      end if

      ! --- main case ---

      sw_idx = 1
      do sbidx = 1, swri%max_num_bessels
         lval = swri%sphbessels(centrex%species_number,sbidx)%lval
         qidx = swri%sphbessels(centrex%species_number,sbidx)%qidx
         qval = swri%sphbessels(centrex%species_number,sbidx)%qval

         if (lval == -1) exit
         if (lval < swexside_filter%min_l) cycle
         if (lval > swexside_filter%max_l) cycle
         if (qidx > swexside_filter%max_q) cycle

         ! jd: If outside the localization radius, then all SWs are trivially
         !     zero. This happens for PPDs of overlapping NGWFs, for points
         !     in PPDs that are > radius.
         if(rad > centrex%radius) then
            do mval=-lval,lval
               sw_val(sw_idx) = 0.0_DP
               sw_idx = sw_idx + 1
            end do
         else
            sphbess = sw_bessel_fast(lval, rad*rad*qval*qval) * (-1.0_DP)**lval
            do mval=-lval,lval
               sw_val(sw_idx) = sphbess * &
                    sw_real_sph_harm_unit(unit_disp%x,unit_disp%y,unit_disp%z,&
                    lval,mval)
               sw_idx = sw_idx + 1
            end do
         end if

      end do

      n_sws = sw_idx - 1 ! (in parent)

    end subroutine internal_sws_at_point

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_swpots_at_point(sw_pot, curpoint, centrex, pbc, &
         which_image)

      use constants, only: PI, SAFE_DIV_EPS
      use geometry, only: point, operator(*), operator(-), geometry_MAGNITUDE
      use rundat, only: pub_hfx_debug
      use spherical_wave, only: sw_real_sph_harm_unit
      use simulation_cell, only: minimum_image_displacement

      implicit none

      ! jd: Arguments
      real(kind=DP),     intent(out) :: sw_pot(swri%quality%max_sws_per_centre)
      type(POINT),       intent(in)  :: curpoint
      type(ATOM_CENTRE), intent(in)  :: centrex
      logical,           intent(in)  :: pbc
      type(PBC_IMAGE_INFO), intent(out), optional :: which_image

      ! jd: Local variables
      real(kind=DP) :: invrad, rad
      integer :: lval
      integer :: qidx
      integer :: species_number
      integer :: sbidx, mval
      real(kind=DP) :: factor, bessint
      type(POINT) :: disp, unit_disp
      real(kind=DP) :: common_factor
      integer :: sw_idx

      ! -----------------------------------------------------------------------

      sw_pot = 0.0_DP
      species_number = centrex%species_number
      common_factor = 4.0_DP * PI

      if(pbc) then
         ! jd: PBC version
         if(.not. pub_hfx_debug) then
            ! jd: Usual PBC version
            call minimum_image_displacement(rad, disp, &
                 curpoint, centrex%incell, cell)
         else
            ! jd: Debugged PBC version
            if(present(which_image)) then
               call minimum_image_displacement(rad, disp, &
                    curpoint, centrex%incell, cell, which_image%image)
               which_image%curpoint = curpoint
               which_image%centre = centrex%incell
               which_image%disp = disp
               which_image%rad = rad
               which_image%disp0 = curpoint - centrex%incell
               which_image%rad0 = geometry_MAGNITUDE(which_image%disp0)
            else
               call utils_abort('internal_swpots_at_point: &
                    &Optional argument missing')
            end if
         end if
      else
         ! jd: OBC version
         disp = curpoint - centrex%incell
         rad = geometry_MAGNITUDE(disp)
         if(present(which_image)) then
            which_image%image = 0
            which_image%curpoint = curpoint
            which_image%centre = centrex%incell
            which_image%disp = disp
            which_image%rad = rad
            which_image%disp0 = disp
            which_image%rad0 = rad
         end if
      end if

      sw_idx = 1
      do sbidx = 1, swri%max_num_bessels
         lval = swri%sphbessels(species_number,sbidx)%lval
         qidx = swri%sphbessels(species_number,sbidx)%qidx

         if (lval == -1) exit
         if (lval < swexside_filter%min_l) cycle
         if (lval > swexside_filter%max_l) cycle
         if (qidx > swexside_filter%max_q) cycle

         ! qoh: At rad = 0 only consider l=0 (we have 0/0 for other values of l)
         ! jd:  Mathematica seems to indicate all swpots evaluate to 0 for l>0
         !      when rad is 0, so it should be safe just to ignore them, only
         ! [*]  making sure to fast forward to the right n_sws and store 0's.

         if(rad > SAFE_DIV_EPS) then
            ! jd: Usual case
            invrad = 1.0_DP / rad
            unit_disp = invrad * disp
            bessint = swri_sph_bess_pot_int(rad,invrad,&
                 swri%sphbessels(species_number,sbidx))
            factor = common_factor / real(((2*lval+1) * (-1)**lval),kind=DP)
            do mval=-lval,lval
               sw_pot(sw_idx) = bessint * factor * &
                    sw_real_sph_harm_unit(unit_disp%x,unit_disp%y,unit_disp%z,&
                    lval,mval)
               sw_idx = sw_idx + 1
            end do
         else
            ! jd: Corner case of SWpot centre on gridpoint where we evaluate
            if(lval == 0) then
               ! z_00 = 0.5_DP/sqrt(pi) = 0.282094791773878_DP
               bessint = 1.0_DP / (swri%sphbessels(species_number,sbidx)%qval**2) &
                    + swri%sphbessels(species_number,sbidx)%nearpotint
            else
               bessint = 0.0_DP ! [*]
            end if
            do mval=-lval,lval
               sw_pot(sw_idx) = bessint * common_factor * 0.282094791773878_DP
               sw_idx = sw_idx + 1
            end do
         end if
      end do

      n_sws = sw_idx - 1 ! (in parent)

    end subroutine internal_swpots_at_point

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine swri_swop_calc_all_in_ppd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_swop_calc_all_in_ppd_fast(swops_in_ppd, n_sws, &    ! out
       src_centre, ppd_idx, sw_or_swpot, cell, swri, pbc, &           ! in
       sample, sample_period, which_images)                           ! in
    !==========================================================================!
    ! Generates the density or potential of all spherical waves coming from    !
    ! src_centre for all points of a ppd. Whether it's the density or the      !
    ! potential that is generated depends on sw_or_swpot.                      !
    ! This routine prioritises speed over readability. It doesn't work with    !
    ! swri-to-swex filtering.                                                  !
    ! Overlap metric (calculation of SWs) has not been optimised as carefully  !
    ! as the electrostatic metric (SWpots).                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   swops_in_ppd (out): Output buffer, at least (n_pts * n_sws) big.       !
    !   n_sws (out):        The number of SWs written to output.               !
    !   src_centre (in):    Centre generating the SW densities/potentials.     !
    !   ppd_idx (in):       The index of the PPD, needed to find its position. !
    !   sw_or_swpot (in):   If 'S' density (SWs) is generated, if 'P', pot.    !
    !   cell (in):          The usual.                                         !
    !   swri (in):          Needed to dimension a local buffer and for Bessel  !
    !                       book-keeping in the internal.                      !
    !   pbc (in):           Pass .true. in PBCs, .false. otherwise.            !
    !   sample (in):       } Helper arguments used to only time the subroutine !
    !   sample_period (in):} every sample_period calls.                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2019, based on                      !
    ! swri_swop_calc_all_in_ppd() whose history is as follows:                 !
    ! Written by Jacek Dziedzic in April 2012 using pot_all_in_tb as a         !
    ! template. Rewritten in June 2012 to use PPDs. Updated in February 2015   !
    ! to use swri and to employ filter.                                        !
    !==========================================================================!

    use basis, only: basis_ppd_location
    use geometry, only: POINT, operator(*), operator(-), operator(+)
    use rundat, only: pub_hfx_debug
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock

    implicit none

    ! jd: Arguments
    real(kind=DP),     intent(out) :: swops_in_ppd(:)
    integer,           intent(out) :: n_sws
    type(atom_centre), intent(in)  :: src_centre
    integer,           intent(in)  :: ppd_idx
    character,         intent(in)  :: sw_or_swpot
    type(CELL_INFO),   intent(in)  :: cell
    type(SW_RI),       intent(in)  :: swri
    logical,           intent(in)  :: pbc
    integer,           intent(in)  :: sample
    integer,           intent(in)  :: sample_period
    type(PBC_IMAGE_INFO), intent(out), optional :: which_images(cell%n_pts)

    ! jd: Local variables
    type(POINT)   :: dest_corner1
    type(POINT)   :: curpoint
    type(POINT)   :: a1, a2, a3
    integer       :: d1idx, d2idx, d3idx
    integer       :: offs_to_point
    integer       :: sw
    integer       :: max_lval
    integer       :: i_pt_in_ppd
    real(kind=DP) :: buffer(swri%quality%max_sws_per_centre)
    character(len=*), parameter :: myself = 'swri_swop_calc_all_in_ppd_fast'

    ! ------------------------------------------------------------------------

    if(mod(sample,sample_period) == 0) then
       call timer_clock(myself,1,multiplier = sample_period)
    end if
    ! don't utils_trace_in/out this subroutine -- it's called from OMP contexts

    n_sws = swri_num_sws_on_centre(swri,src_centre)

    dest_corner1 = basis_ppd_location(ppd_idx,cell)

    a1 = cell%d1 * cell%a1_unit
    a2 = cell%d2 * cell%a2_unit
    a3 = cell%d3 * cell%a3_unit

    max_lval = maxval(swri%sphbessels(src_centre%species_number,:)%lval)
    offs_to_point = 1

    curpoint = dest_corner1 - a1 - a2 -a3

    ! NB: The ordering of swops_in_ppd wrt cache-friendliness is crucial here.
    !     It's much faster to calculate all SWs for a point at once, hence
    !     internal_* calculates all SWs at once. However, later on it's crucial
    !     that the data is organized with the SW being the major index and the
    !     PPD point being the minor index. Thus the cache-unfriendly loop over
    !     SWs here that re-orders the layout. It's much, much cheaper to pay
    !     the cache penalty here, when generating the SWs which are going to be
    !     cached a moment later than to pay it many times when accessing the
    !     data fetched from the cache.
    i_pt_in_ppd = 1
    do d3idx = 1, cell%n_pt3
       curpoint = curpoint + a3
       do d2idx = 1, cell%n_pt2
          curpoint = curpoint + a2
          do d1idx = 1, cell%n_pt1
             curpoint = curpoint + a1

             if(sw_or_swpot == 'P') then
                if(pub_hfx_debug) then
                   call internal_swpots_at_point_fast(buffer, curpoint, &
                        src_centre, max_lval, pbc, which_images(i_pt_in_ppd))
                else
                   call internal_swpots_at_point_fast(buffer, curpoint, &
                        src_centre, max_lval, pbc)
                end if
             else
                call internal_sws_at_point(buffer, curpoint, src_centre, pbc)
             end if

             do sw = 1, n_sws
                swops_in_ppd((sw-1)*cell%n_pts + offs_to_point) = buffer(sw)
             end do

             offs_to_point = offs_to_point+1
             i_pt_in_ppd = i_pt_in_ppd + 1

          end do
          curpoint = curpoint - real(cell%n_pt1, kind=DP) * a1
       end do
       curpoint = curpoint - real(cell%n_pt2,kind=DP) * a2
    end do

    if(mod(sample,sample_period) == 0) then
       call timer_clock(myself,2,multiplier = sample_period)
    end if

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_sws_at_point(sw_val, curpoint, centrex, pbc)
    ! @optimizeme similarly to swpots_at_point below

      use constants, only: SAFE_DIV_EPS, SQRT_PI
      use geometry, only: point, operator(*), operator(-), geometry_MAGNITUDE
      use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm_unit
      use simulation_cell, only: minimum_image_displacement

      implicit none

      ! jd: Arguments
      real(kind=DP),     intent(out) :: sw_val(swri%quality%max_sws_per_centre)
      type(POINT),       intent(in)  :: curpoint
      type(ATOM_CENTRE), intent(in)  :: centrex
      logical,           intent(in)  :: pbc

      ! jd: Local variables
      real(kind=DP) :: rad
      real(kind=DP) :: radsq
      integer       :: lval
      integer       :: qidx
      integer       :: mval
      real(kind=DP) :: qval
      integer       :: sbidx
      real(kind=DP) :: sphbess
      type(POINT)   :: disp, unit_disp
      integer       :: sw_idx

      !@simplify loops for the case of 0

      ! ----------------------------------------------------------------------

      sw_val = 0.0_DP

      if(pbc) then
         ! jd: PBC version
         call minimum_image_displacement(rad, disp, curpoint, centrex%incell, cell)
      else
         ! jd: OBC version
         disp = curpoint - centrex%incell
         rad = geometry_MAGNITUDE(disp)
      end if

      if(rad > SAFE_DIV_EPS) then
         unit_disp = 1.0_DP / rad * disp
      else
         ! jd: Corner case of r=0. All SWs are 0, except for l=0, m=0.
         sw_idx = 1
         do sbidx = 1, swri%max_num_bessels
            lval = swri%sphbessels(centrex%species_number,sbidx)%lval
            if (lval == -1) exit
            qidx = swri%sphbessels(centrex%species_number,sbidx)%qidx
            if (qidx > swri%quality%max_q) cycle

            do mval=-lval,lval
               sw_val(sw_idx) = merge(1.0_DP/(2.0_DP * SQRT_PI), 0.0_DP, &
                    lval == 0 .and. mval == 0)
               sw_idx = sw_idx + 1
            end do

         end do

         return

      end if

      ! --- main case ---

      radsq = rad*rad

      sw_idx = 1
      do sbidx = 1, swri%max_num_bessels
         lval = swri%sphbessels(centrex%species_number,sbidx)%lval
         if (lval == -1) exit
         qidx = swri%sphbessels(centrex%species_number,sbidx)%qidx
         if (qidx > swri%quality%max_q) cycle
         qval = swri%sphbessels(centrex%species_number,sbidx)%qval

         ! jd: If outside the localization radius, then all SWs are trivially
         !     zero. This happens for PPDs of overlapping NGWFs, for points
         !     in PPDs that are > radius.
         if(rad > centrex%radius) then
            sw_val(sw_idx:sw_idx+(2*lval)) = 0.0_DP ! simulates loop over
            sw_idx = sw_idx + 2*lval+1              ! m=-lval to lval
         else
            sphbess = sw_bessel_fast(lval, radsq*qval*qval)
            if(lval == 1 .or. lval == 3) sphbess = -sphbess
            do mval=-lval,lval
               sw_val(sw_idx) = sphbess * &
                    sw_real_sph_harm_unit(unit_disp%x,unit_disp%y,unit_disp%z,&
                    lval,mval)
               sw_idx = sw_idx + 1
            end do
         end if

      end do

    end subroutine internal_sws_at_point

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_swpots_at_point_fast(sw_pot, &
         curpoint, centrex, max_lval, pbc, which_image)

      use constants, only: PI, SAFE_DIV_EPS
      use geometry, only: point, operator(*), operator(-), geometry_MAGNITUDE
      use rundat, only: pub_hfx_bessel_rad_nptsx
      use spherical_wave, only: sw_real_sph_harm_unit_all_l_all_m
      use simulation_cell, only: minimum_image_displacement
      use utils, only: utils_abort

      implicit none

      ! jd: Arguments
      real(kind=DP), intent(out)     :: sw_pot(swri%quality%max_sws_per_centre)
      type(POINT), intent(in)        :: curpoint
      type(ATOM_CENTRE), intent(in)  :: centrex
      integer, intent(in)            :: max_lval
      logical, intent(in)            :: pbc
      type(PBC_IMAGE_INFO), intent(out), optional :: which_image

      ! jd: Local variables
      real(kind=DP) :: invrad, rad
      integer :: lval
      integer :: qidx
      integer :: species_number
      integer :: sbidx, mval
      real(kind=DP) :: factor, bessint
      type(POINT) :: disp, unit_disp
      integer :: sw_idx
      integer :: xpt_1, xpt_2
      real(kind=DP) :: xmax
      real(kind=DP) :: frac, y1, y2, hx, rad1
      real(kind=DP) :: sph_harm_buff(25) ! jd: sufficient for l up to 4
      integer, parameter :: loffset(0:4) = (/1,2,5,10,17/) ! jd: ditto
      real(kind=DP), parameter :: lfactors(0:4) = &
           (/4.0_DP*PI/real(2*0+1,kind=DP), &
           -4.0_DP*PI/real(2*1+1,kind=DP), &
           4.0_DP*PI/real(2*2+1,kind=DP), &
           -4.0_DP*PI/real(2*3+1,kind=DP), &
           4.0_DP*PI/real(2*4+1,kind=DP)/)

      ! -----------------------------------------------------------------------

      species_number = centrex%species_number

      xmax = geometry_MAGNITUDE(cell%a1 + cell%a2 + cell%a3)

      if(pbc) then
         ! -----------------
         ! jd: PBC version
         ! -----------------
         if(.not. pub_hfx_debug) then
            ! jd: Usual PBC version
            call minimum_image_displacement(rad, disp, &
                 curpoint, centrex%incell, cell)
         else
            ! jd: Debugged PBC version
            if(present(which_image)) then
               call minimum_image_displacement(rad, disp, &
                    curpoint, centrex%incell, cell, which_image%image)
               which_image%curpoint = curpoint
               which_image%centre = centrex%incell
               which_image%disp = disp
               which_image%rad = rad
               which_image%disp0 = curpoint - centrex%incell
               which_image%rad0 = geometry_MAGNITUDE(which_image%disp0)
            else
               call utils_abort('internal_swpots_at_point_fast: &
                    &Optional argument missing')
            end if
         end if

      else
         ! -----------------
         ! jd: OBC version
         ! -----------------
         disp = curpoint - centrex%incell
         rad = geometry_MAGNITUDE(disp)
         if(present(which_image)) then
            which_image%image = 0
            which_image%curpoint = curpoint
            which_image%centre = centrex%incell
            which_image%disp = disp
            which_image%rad = rad
            which_image%disp0 = disp
            which_image%rad0 = rad
         end if
      end if

      if(rad > SAFE_DIV_EPS) then
         ! jd: Usual case -- point not practically exactly on atom centre
         invrad = 1.0_DP / rad
         unit_disp = invrad * disp

         ! jd: Calculate all the spherical harmonics -- all l, all m, but
         !     no dependence on q.
         call sw_real_sph_harm_unit_all_l_all_m(sph_harm_buff, &
              unit_disp%x, unit_disp%y, unit_disp%z, max_lval)

         xpt_1 = int(rad / xmax * real(pub_hfx_bessel_rad_nptsx-1,kind=DP)) + 1
         xpt_2 = xpt_1 + 1
         if(xpt_2 > pub_hfx_bessel_rad_nptsx) then
            xpt_2 = xpt_1
         end if
         hx = xmax/(pub_hfx_bessel_rad_nptsx-1)

         rad1 = (xpt_1-1) * hx
         frac = (rad-rad1)/hx

         sw_idx = 1
         do sbidx = 1, swri%max_num_bessels
            lval = swri%sphbessels(species_number,sbidx)%lval
            if (lval == -1) exit ! jd: signals the end

            ! jd: Calculate Bessel integral for this {l,q}
            ! Different options have been commented out

! 1. Direct calculation
!            bessint = swri_sph_bess_pot_int(rad, invrad, &
!                 swri%sphbessels(species_number,sbidx))

! 2. Quartic interpolation -- proves 50% slower (!)
!            bessint = services_1d_interpolation(swri%rad_lookup(:,sbidx), &
!                 pub_hfx_bessel_rad_nptsx, rad / xmax * &
!                 real(pub_hfx_bessel_rad_nptsx-1,kind=DP),0) ! -1, sic!

! 3. Zero-order interpolation -- not nearly accurate enough, unless npts_x is
!    increased to unreasonably large values, then is much slower than direct
!    calculation, possibly because the giant lookup is passed by value somewhere
!              bessint = swri%rad_lookup(xpt,sbidx)

! 4. Simple 1st order interpolation -- makes SWOP calcs about 21% faster,
!    meaning about 50% cut in walltime of radial part.

            y1 = swri%rad_lookup(xpt_1,sbidx)
            y2 = swri%rad_lookup(xpt_2,sbidx)
            bessint = y1 + frac * (y2-y1)

            ! jd: Do the entire 'm' loop in one go
            sw_pot(sw_idx:sw_idx+2*lval) = bessint * lfactors(lval) * &
                 sph_harm_buff(loffset(lval):loffset(lval)+2*lval)
            sw_idx = sw_idx + 2*lval + 1
         end do

      else
         ! jd: Corner case of SWpot centre on gridpoint where we evaluate
         sw_idx = 1
         do sbidx = 1, swri%max_num_bessels
            lval = swri%sphbessels(species_number,sbidx)%lval
            if (lval == -1) exit ! jd: signals the end

            if(lval == 0) then
               ! z_00 = 0.5_DP/sqrt(pi) = 0.282094791773878_DP
               bessint = 1.0_DP / (swri%sphbessels(species_number,sbidx)%qval**2) &
                    + swri%sphbessels(species_number,sbidx)%nearpotint
            else
               bessint = 0.0_DP ! [*]
            end if
            do mval=-lval,lval
               sw_pot(sw_idx) = bessint * 4.0_DP * PI * 0.282094791773878_DP
               sw_idx = sw_idx + 1
            end do
         end do
      end if

    end subroutine internal_swpots_at_point_fast

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine swri_swop_calc_all_in_ppd_fast

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_distribute_atomblocks(atomblock_first, atomblock_count, &
       n_atomblocks)
    !==========================================================================!
    ! Distributes metric matrix atomblocks across MPI ranks.                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2013.                              !
    !==========================================================================!

    use comms, only: pub_total_num_procs
    use rundat, only: pub_debug_on_root
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(out) :: atomblock_first(0:pub_total_num_procs-1)
    integer, intent(out) :: atomblock_count(0:pub_total_num_procs-1)
    integer, intent(in)  :: n_atomblocks

    ! Local variables
    integer :: portion, remainder
    integer :: p
    integer :: k
    real(kind=DP) :: delta

    character(len=*), parameter :: myself = 'swri_distribute_atomblocks'

    ! -------------------------------------------------------------------------

    portion = n_atomblocks / pub_total_num_procs
    remainder = n_atomblocks - portion * pub_total_num_procs

    ! Deal out 'portion' to every proc
    atomblock_count(:) = portion

    ! Deal out the remainder to every delta'th processor, if necessary
    if(remainder /= 0) then
       delta = real(pub_total_num_procs,kind=DP) / real(remainder,kind=DP)
       do k = 1, remainder
          p = nint(real(k-1,kind=DP) * delta)
          atomblock_count(p) = atomblock_count(p) + 1
       end do
    end if

    call utils_assert(sum(atomblock_count) == n_atomblocks, &
         'Sanity check failed in '//myself)

    ! Construct sw_first from sw_count
    do p = 0, pub_total_num_procs - 1
       if(p == 0) then
          atomblock_first(p) = 1
       else
          atomblock_first(p) = atomblock_first(p-1) + atomblock_count(p-1)
       end if

       if(atomblock_count(p) == 0) atomblock_first(p) = -1

       if(pub_debug_on_root) then
          if(atomblock_count(p) > 0) then
             write(stdout,'(a,i4,a,i4,a,i4,a)') '      Proc ',p, &
                  ': atomblocks ', atomblock_first(p), ' to ', &
                  atomblock_first(p)+atomblock_count(p)-1,'.'
          else
             write(stdout,'(a,i4,a)') 'Proc ',p,': no atomblocks.'
          end if
       end if

    end do

  end subroutine swri_distribute_atomblocks

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_num_sph_functions(radtable, max_num_bessels, & ! output
       quality, par, &                                           ! input
       num_sw, &                                                 ! opt output
       num_q_per_l )                                             ! opt output
    !=========================================================================!
    ! This subroutine fills radtable (later stored as member of SW_RI),       !
    ! calculates max_num_bessels (also stored) and num_sw (not stored, only   !
    ! used briefly in swri_init_metric_matrix()).                             !
    !-------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2009.                                        !
    ! Subsequently modified by Jacek Dziedzic between 2012 and 2015.          !
    ! Modified by James C. Womack to optionally output num_q_per_l, 11/2018   !
    !=========================================================================!

    use comms, only: comms_reduce, pub_my_proc_id, pub_on_root
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_cutoff_energy
    use spherical_wave, only: sw_bessel_zeros_init, sw_bessel_zeros

    implicit none

    ! jd: Arguments
    !                                 NGWF radii for all species
    type(PARAL_INFO), intent(in)   :: par
    real(kind=DP), intent(out)     :: radtable(par%num_species)
    integer, intent(out)           :: max_num_bessels
    type(SW_QUALITY)               :: quality
    integer, intent(out), optional :: num_sw(par%num_species)
    ! JCW: Maximum allowed qvalues for each l value (on each species)
    integer, intent(out), optional :: num_q_per_l(par%num_species,&
                                      quality%min_l:quality%max_l)

    ! jd: Local variables
    integer :: lval ! Quantum number l
    integer :: qidx ! Zero index
    integer :: ridx ! Radius index
    integer :: eidx ! Element index
    integer :: nbessels  ! Number of spherical bessels (per species)
    integer :: nsphwaves ! Number of spherical waves (per species)
    integer :: num_q     ! Number of allowed q-values (per l, per species)
    real(kind=DP) :: maxzero ! Maximum acceptable Bessel zero
    integer, save :: lmax_warnings = 0 ! So that warnings are not repeated
    character(len=*), parameter :: myself = 'swri_num_sph_functions'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)

    max_num_bessels = 0

    call sw_bessel_zeros_init(quality%max_q, quality%max_l)

    radtable = -1.0_DP
    do ridx = 1, par%num_species
       do eidx=1,par%num_atoms_on_proc(pub_my_proc_id)
          if (par%elements_on_proc(eidx)%species_number == ridx) then
             radtable(ridx)=par%elements_on_proc(eidx)%radius
             exit
          end if
       end do
    end do

    !qoh: Get radii of all species
    call comms_reduce('MAX',radtable,par%num_species)

    do ridx = 1, par%num_species
       nbessels = 0
       nsphwaves = 0

       ! qoh: qmax = sqrt(2.0_DP*pub_cutoff_energy)
       maxzero = sqrt(2.0_DP*pub_cutoff_energy)*radtable(ridx)


       do lval = quality%min_l, quality%max_l
          ! JCW: num_q is number of allowed spherical Bessel q-values
          ! JCW: for each l-value, for each species
          num_q = 0
          do qidx = 1, quality%max_q
             if (sw_bessel_zeros(qidx,lval) .lt. maxzero) then
                nbessels = nbessels + 1
                nsphwaves=nsphwaves + 1 + 2*lval
                num_q = num_q + 1
             else
                if (pub_on_root .and. lmax_warnings >= 0) then
                   write(stdout,'(a,i0,a,i0,a,i0,a)') &
                        'WARNING: Spherical waves with l=',lval,' and qi=', &
                        qidx,' will be ignored for species #',ridx,'.'
                   lmax_warnings = 1
                end if
             end if
          end do
          if (present(num_q_per_l)) num_q_per_l(ridx,lval) = num_q
       end do
       max_num_bessels = max(nbessels,max_num_bessels)
       if (present(num_sw)) num_sw(ridx) = nsphwaves
    end do

    ! jd: Ensure that once warnings were issued, they will not be repeated
    lmax_warnings = -abs(lmax_warnings)

    call utils_trace_out(myself)

  end subroutine swri_num_sph_functions

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_sph_bessels_init(swri, uni_tightbox, num_species)
    !==========================================================================!
    ! This subroutine initialises the sphbessels component of an SW_RI.        !
    ! It defines the spherical Bessel functions and calculates parts of the    !
    ! potential integrals.                                                     !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2009.                                         !
    ! Adapted to SW_RI representation by Jacek Dziedzic in February 2015.      !
    !==========================================================================!

    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_cutoff_energy
    use spherical_wave, only: sw_bessel_zeros, sw_init

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout)   :: swri
    type(FFTBOX_INFO),intent(in) :: uni_tightbox
    integer, intent(in)          :: num_species

    ! jd: Local variables
    integer :: lval  ! Quantum number l
    integer :: qidx  ! zero index
    integer :: ridx  ! Radius index
    integer :: sbidx ! spherical bessel idx
    real(kind=DP) :: maxzero ! maximum acceptable Bessel zero
    real(kind=DP) :: rmaxsq  ! Maximum arguement of Bessel function squared
    real(kind=DP) :: gmaxsq  ! Maximum square magnitude of G vector
    integer :: n1half, n2half, n3half
    character(len=*), parameter :: myself = 'swri_sph_bessels_init'

    ! -------------------------------------------------------------------------

    call utils_trace_in(myself)

    swri%sphbessels%lval = -1

    do ridx=1, num_species

       ! qoh: qmax = sqrt(2.0_DP*pub_cutoff_energy)
       maxzero = sqrt(2.0_DP*pub_cutoff_energy)*swri%radtable(ridx)
       sbidx = 0
       do lval = swri%quality%min_l, swri%quality%max_l
          do qidx = 1, swri%quality%max_q
             if ( sw_bessel_zeros(qidx,lval) .lt. maxzero) then
                sbidx = sbidx + 1
                swri%sphbessels(ridx,sbidx)%lval = lval
                swri%sphbessels(ridx,sbidx)%qidx = qidx
                swri%sphbessels(ridx,sbidx)%qval = sw_bessel_zeros(qidx,lval)&
                     /swri%radtable(ridx)
                swri%sphbessels(ridx,sbidx)%aval = swri%radtable(ridx)
                call internal_const_pot
             end if
          end do ! over Bessel zeros
       end do ! over l
    end do ! over species

    ! qoh: Need to set rmax as default is too small
    n1half = uni_tightbox%total_pt1/2+1
    n2half = uni_tightbox%total_pt2/2+1
    n3half = uni_tightbox%total_pt3/2+1

    gmaxsq = 2.0_DP*uni_tightbox%recip_grid(5,n1half,&
         n2half,n3half)
    rmaxsq = gmaxsq*(maxval(swri%radtable)**2)

    ! qoh: Initialise SW module
    call sw_init(swri%quality%max_l, swri%quality%max_q, rmaxsq)

    call utils_trace_out(myself)

  contains

    subroutine internal_const_pot

      ! qoh: Calculate the constant part of the spherical Bessel potential
      ! qoh: integrals.

      use utils, only: utils_abort

      ! jd: Formulas verified with Mathematica

      implicit none
      real(kind=DP) :: qa
      real(kind=DP) :: qval

      qa = swri%sphbessels(ridx,sbidx)%qval * swri%sphbessels(ridx,sbidx)%aval
      qval = swri%sphbessels(ridx,sbidx)%qval

      select case (lval)
      case(0)
         swri%sphbessels(ridx,sbidx)%farpotint = &
              (sin(qa) -qa*cos(qa))/ (qval**3) ! jd: OK
         swri%sphbessels(ridx,sbidx)%nearpotint = &
              - cos(qa) / (qval*qval) ! jd: OK
      case(1)
         swri%sphbessels(ridx,sbidx)%farpotint = ((3.0_DP- qa*qa )*sin(qa) &
              - 3.0_DP*qa*cos(qa))/(qval**4) ! jd: OK
         swri%sphbessels(ridx,sbidx)%nearpotint = - sin(qa) /(qval*qa) ! jd: OK
      case(2)
         swri%sphbessels(ridx,sbidx)%farpotint = (qa*(qa*qa-15.0_DP)*cos(qa) - &
              (6.0_DP*qa*qa-15.0_DP)*sin(qa))/(qval**5) ! jd: OK
         swri%sphbessels(ridx,sbidx)%nearpotint = &
              (qa*cos(qa) - sin(qa))/(qa**3) ! jd: OK
      case(3)
         swri%sphbessels(ridx,sbidx)%farpotint = &
              ((qa**4 -45.0_DP*qa*qa + 105.0_DP) &
              * sin(qa) + qa*(10.0_DP*qa*qa - 105.0_DP)*cos(qa))/(qval**6) ! jd: OK
         swri%sphbessels(ridx,sbidx)%nearpotint = &
              qval*(3.0_DP * qa * cos (qa) +&
              (qa*qa-3.0_DP)*sin(qa))/(qa**5) ! jd: OK
      case(4)
         swri%sphbessels(ridx,sbidx)%farpotint = &
              ((15.0_DP * qa**4 - 420.0_DP*qa*qa&
              + 945.0_DP) *sin(qa) &
              -qa*(qa**4 - 105.0_DP*qa*qa + 945.0_DP)*cos(qa))/ (qval**7) ! jd: OK
         swri%sphbessels(ridx,sbidx)%nearpotint = &
              qval*qval*(qa*(15.0_DP - qa*qa)&
              *cos(qa) + (6.0_DP * qa*qa - 15.0_DP) * sin(qa)) / (qa**7) ! jd: OK

      case default
         call utils_abort('internal_const_pot: Angular momenta >4 unsupported')

      end select

    end subroutine internal_const_pot

  end subroutine swri_sph_bessels_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function swri_num_sws_on_centre(swri,centrea)
    !==========================================================================!
    ! This subroutine calculates the size of the set of spherical waves        !
    ! centred on centrea.                                                      !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in September 2009                                !
    ! Modified by Jacek Dziedzic in June 2012.                                 !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(in)       :: swri
    type(atom_centre), intent(in) :: centrea

    ! jd: Local variables
    integer :: sbidx
    integer :: lval

    ! -------------------------------------------------------------------------

    swri_num_sws_on_centre = 0

    do sbidx=1, swri%max_num_bessels
       lval = swri%sphbessels(centrea%species_number,sbidx)%lval
       if (lval == -1) exit
       swri_num_sws_on_centre = swri_num_sws_on_centre + 2*lval + 1
    end do

  end function swri_num_sws_on_centre

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_atomblock_filter_to_swex(swexside_atomblock, & ! out
       swriside_atomblock, swexside_filter, swri)   ! in
    !==========================================================================!
    ! Extracts a swexside atomblock from a swriside atomblock.                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 19/02/2015.                                 !
    !==========================================================================!

    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)   :: swexside_atomblock(:,:)
    ! jd: ^ Both indices run from 1 to swex_num_sws_per_centre
    !       Using assumed-shape to enforce passing by descriptor and avoid copy
    real(kind=DP), intent(in)    :: swriside_atomblock(:,:)
    ! jd: ^ Both indices run from 1 to swri_num_sws_per_centre
    !       Using assumed-shape to enforce passing by descriptor and avoid copy
    type(SW_QUALITY), intent(in) :: swexside_filter
    type(SW_RI), intent(in)      :: swri

    ! jd: Local variables
    integer :: swex_sw_a, swex_sw_b
    integer :: swri_sw_a, swri_sw_b
    integer :: bessel_a, bessel_b
    integer :: lval_a, lval_b
    integer :: qidx_a, qidx_b
    integer :: species_a, species_b
    character(len=*), parameter :: myself = 'swri_atomblock_filter_to_swex'

    ! -------------------------------------------------------------------------

    ! jd: Short-circuit for the case of no conversion
    if(swexside_filter%max_l == swri%quality%max_l .and. &
         swexside_filter%max_q == swri%quality%max_q .and. &
         swexside_filter%min_l == swri%quality%min_l) then
       swexside_atomblock(:,:) = swriside_atomblock(:,:)
       return
    end if

    species_a = 1 ! @IDENTICAL_RADII
    species_b = 1 ! @IDENTICAL_RADII

    ! jd: Standard case of having to convert
    swex_sw_a = 1
    do swri_sw_a = 1, swri%quality%num_sws_per_centre
       bessel_a = swri%sw_idx_to_bessel_idx(swri_sw_a)
       lval_a = swri%sphbessels(species_a,bessel_a)%lval
       qidx_a = swri%sphbessels(species_a,bessel_a)%qidx

       if (lval_a == -1) exit
       if (lval_a < swexside_filter%min_l) cycle
       if (lval_a > swexside_filter%max_l) cycle
       if (qidx_a > swexside_filter%max_q) cycle

       swex_sw_b = 1
       do swri_sw_b = 1, swri%quality%num_sws_per_centre
          bessel_b = swri%sw_idx_to_bessel_idx(swri_sw_b)
          lval_b = swri%sphbessels(species_b,bessel_b)%lval
          qidx_b = swri%sphbessels(species_b,bessel_b)%qidx

          if (lval_b == -1) exit
          if (lval_b < swexside_filter%min_l) cycle
          if (lval_b > swexside_filter%max_l) cycle
          if (qidx_b > swexside_filter%max_q) cycle

          swexside_atomblock(swex_sw_a, swex_sw_b) = &
               swriside_atomblock(swri_sw_a, swri_sw_b)

          swex_sw_b = swex_sw_b + 1

       end do

       swex_sw_a = swex_sw_a + 1

    end do

    swex_sw_a = swex_sw_a - 1
    swex_sw_b = swex_sw_b - 1

    call utils_assert(swex_sw_a == swexside_filter%num_sws_per_centre, &
         myself//': Internal error [1]', swex_sw_a, &
         swexside_filter%num_sws_per_centre)
    call utils_assert(swex_sw_b == swexside_filter%num_sws_per_centre, &
         myself//': Internal error [2]', swex_sw_b, &
         swexside_filter%num_sws_per_centre)

  end subroutine swri_atomblock_filter_to_swex

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_init_persistent_caches(swri, cell_n_pts)
    !==================================================================!
    ! Initialises the cache structures that need to persist across     !
    ! calls. These are the SWOP caches, currently.                     !
    !------------------------------------------------------------------!
    !==================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: SIZEOF_DOUBLE
    use hash_table, only: hash_table_init, hash_table_calc_cache_size
    use rundat, only: pub_cache_limit_for_swops, pub_cache_limit_for_swops2, &
         pub_swx_dbl_grid, pub_pol_emb_dbl_grid
    use utils, only: utils_unit, utils_abort

    implicit none

    ! Arguments
    type(SW_RI), intent(inout) :: swri
    integer, intent(in)        :: cell_n_pts

    ! jd: Local variables
    integer :: max_cached_swops_in_ppds
    integer :: max_cached_swops_at_points

    ! -------------------------------------------------------------------------

    ! jd: Open a file for debugging hash-tables, if desired
    swri%hash_table_info_unit = utils_unit()
    write(swri%hash_table_info_filename,'(a,a,a,a,i0,a)') &
         trim(hash_table_info_rootname), '_', trim(swri%swri_name), &
         '_proc_', pub_my_proc_id, '.log'
    open(swri%hash_table_info_unit, file=swri%hash_table_info_filename, err=10)

    max_cached_swops_in_ppds = hash_table_calc_cache_size( &
         pub_cache_limit_for_swops, &
         swri%quality%max_sws_per_centre * cell_n_pts * SIZEOF_DOUBLE)

    max_cached_swops_at_points = hash_table_calc_cache_size( &
         pub_cache_limit_for_swops2, &
         swri%quality%max_sws_per_centre * SIZEOF_DOUBLE)

    ! jd: This is costly, let's use fewer slots to ensure we fill
    !     practically all of them.
    call hash_table_init(swri%swops_in_ppds_ht, 'SWOPS_IN_PPDS', 3, &
         max_cached_swops_in_ppds/4, max_cached_swops_in_ppds, &
         swri%hash_table_info_unit)

    if(pub_swx_dbl_grid .or. pub_pol_emb_dbl_grid) then
       call hash_table_init(swri%swops_at_points_ht, 'SWOPS_AT_POINTS', 5, &
            max_cached_swops_at_points, max_cached_swops_at_points, &
            swri%hash_table_info_unit)
    end if

    return

10  call utils_abort('Could not create file: '//&
         trim(swri%hash_table_info_filename))

  end subroutine swri_init_persistent_caches

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_free_persistent_caches(swri)
    !==================================================================!
    ! Frees the persistent caches.                                     !
    ! These are the SWOPS caches currently.                            !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   swri (inout): The swri whose persistent caches are to be freed.!
    !------------------------------------------------------------------!
    !==================================================================!

    use hash_table, only: hash_table_free
    use rundat, only: pub_swx_dbl_grid, pub_pol_emb_dbl_grid

    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout) :: swri

    ! -------------------------------------------------------------------------

    if(pub_swx_dbl_grid .or. pub_pol_emb_dbl_grid) then
       call hash_table_free(swri%swops_at_points_ht)
    end if
    call hash_table_free(swri%swops_in_ppds_ht)

    close(swri%hash_table_info_unit, err=10)

    return

10  call utils_abort('Could not close file: '//swri%hash_table_info_filename)

  end subroutine swri_free_persistent_caches

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_calc_sws_at_point(swri, cell, sw_val, curpoint, centrex, &
       swexside_filter, pbc)
    !==================================================================!
    ! Calculates all spherical waves at a point. Normally the PPD      !
    ! version is much more efficient. This version is useful on the    !
    ! double grid.                                                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2017.                       !
    !==================================================================!

    use constants, only: SAFE_DIV_EPS, SQRT_PI
    use geometry, only: point, operator(*), operator(-), geometry_MAGNITUDE
    use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm_unit
    use simulation_cell, only: minimum_image_displacement, CELL_INFO

    implicit none

    ! jd: Arguments
    type(SW_RI),       intent(in)  :: swri
    type(SW_QUALITY),  intent(in)  :: swexside_filter
    type(CELL_INFO),   intent(in)  :: cell
    real(kind=DP),     intent(out) :: sw_val(swexside_filter%num_sws_per_centre)
    type(POINT),       intent(in)  :: curpoint
    type(ATOM_CENTRE), intent(in)  :: centrex
    logical,           intent(in)  :: pbc

    ! jd: Local variables
    real(kind=DP) :: rad
    integer       :: lval
    integer       :: qidx
    integer       :: mval
    real(kind=DP) :: qval
    integer       :: sbidx
    real(kind=DP) :: sphbess
    type(POINT)   :: disp, unit_disp
    integer       :: sw_idx

    !@simplify loops for the case of 0

    ! ----------------------------------------------------------------------

    if(pbc) then
       ! jd: PBC version
       call minimum_image_displacement(rad, disp, curpoint, centrex%incell, cell)
    else
       ! jd: OBC version
       disp = curpoint - centrex%incell
       rad = geometry_MAGNITUDE(disp)
    end if

    ! jd: If outside the localization radius, then all SWs are trivially
    !     zero. This happens in PPDs of overlapping NGWFs, for points
    !     in PPDs that are > radius.
    if(rad > centrex%radius) then
       sw_val = 0.0_DP
       return
    end if

    if(rad > SAFE_DIV_EPS) then
       unit_disp = 1.0_DP / rad * disp
    else
       ! jd: Corner case of r=0. All SWs are 0, except for l=0, m=0.
       sw_idx = 1
       do sbidx = 1, swri%max_num_bessels
          lval = swri%sphbessels(centrex%species_number,sbidx)%lval
          if (lval == -1) exit
          qidx = swri%sphbessels(centrex%species_number,sbidx)%qidx
          if (qidx > swri%quality%max_q) cycle

          do mval=-lval,lval
             sw_val(sw_idx) = merge(1.0_DP/(2.0_DP * SQRT_PI), 0.0_DP, &
                  lval == 0 .and. mval == 0)
             sw_idx = sw_idx + 1
          end do

       end do

       return

    end if

    ! --- main case ---

    sw_idx = 1
    do sbidx = 1, swri%max_num_bessels
       lval = swri%sphbessels(centrex%species_number,sbidx)%lval
       qidx = swri%sphbessels(centrex%species_number,sbidx)%qidx
       qval = swri%sphbessels(centrex%species_number,sbidx)%qval

       if (lval == -1) exit
       if (lval < swexside_filter%min_l) cycle
       if (lval > swexside_filter%max_l) cycle
       if (qidx > swexside_filter%max_q) cycle

       sphbess = sw_bessel_fast(lval, rad*rad*qval*qval) * (-1.0_DP)**lval
       do mval = -lval, lval
          sw_val(sw_idx) = sphbess * &
               sw_real_sph_harm_unit(unit_disp%x, unit_disp%y, unit_disp%z, &
               lval, mval)
          sw_idx = sw_idx + 1
       end do
    end do

  end subroutine swri_calc_sws_at_point

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_calc_swpots_at_point(swri, cell, sw_pot, curpoint, centrex, &
       swexside_filter, pbc, rcut)
    !==================================================================!
    ! Calculates all spherical wave potentials at a point. Normally the!
    ! PPD version is much more efficient. This version is useful on the!
    ! double grid.                                                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2017.                       !
    !==================================================================!

    use constants, only: PI, SAFE_DIV_EPS
    use geometry, only: point, operator(*), operator(-), geometry_magnitude
    use spherical_wave, only: sw_real_sph_harm_unit
    use simulation_cell, only: minimum_image_displacement, CELL_INFO

    implicit none

    ! jd: Arguments
    type(SW_RI),       intent(in)       :: swri
    type(SW_QUALITY),  intent(in)       :: swexside_filter
    type(CELL_INFO),   intent(in)       :: cell
    real(kind=DP),     intent(out)      :: sw_pot(swexside_filter%num_sws_per_centre)
    type(POINT),       intent(in)       :: curpoint
    type(ATOM_CENTRE), intent(in)       :: centrex
    logical,           intent(in)       :: pbc
    real(kind=DP), intent(in), optional :: rcut


    ! jd: Local variables
    real(kind=DP) :: invrad, rad
    integer :: lval
    integer :: qidx
    integer :: species_number
    integer :: sbidx, mval
    real(kind=DP) :: factor, bessint
    type(POINT) :: disp, unit_disp
    real(kind=DP) :: common_factor
    integer :: sw_idx

    ! ----------------------------------------------------------------------

    species_number = 1 ! @IDENTICAL_RADII
    common_factor = 4.0_DP * PI

    if(pbc) then
       ! jd: PBC version
       call minimum_image_displacement(rad, disp, curpoint, centrex%incell, cell)
    else
       ! jd: OBC version
       disp = curpoint - centrex%incell
       rad = geometry_MAGNITUDE(disp)
    end if

    if(present(rcut)) then
       if(rad > rcut) then
          sw_pot = 0.0_DP
          return
       end if
    end if

! @opt:
! sw_real_sph_harm_unit does not depend on q
! bessint does not depend on m, there is hidden dependence on l. Depends on q, r.
! factor does not depend on m, q
! Instead of caching SWpots in the sphere, should we maybe keep interpolations
! or lookups for each part?
! a) rad-independent part: factor * harm_unit for all l, m, for all points in
!        sphere around centre. Indexing over points will make unit_disp and MIC
!        calculation redundant, except when generating the lookup.
!        Cost: 130k-300k points in sphere (with fine grid)
!              9lm's in DMA, 25lm's in HFx
!        Total: 8-60MB per centre, meh
!    ^ angular look-up
! b) rad part: 2l's in DMA, 5l's in HFx.
!              15q's in DMA, 20q's in HFx
!              r interpolated every 0.001a0, over box diagonal, say 80a0
!              18-60MB... and shared between al centres, global array.

!
! Low-level routine -- swri_obtain_swops_at_point().
! Then PPD routine and fine-grid routine can use this.

! Could even cache all (potentially) MIC'ed distances in a sphere

    sw_idx = 1
    do sbidx = 1, swri%max_num_bessels
       lval = swri%sphbessels(species_number,sbidx)%lval
       qidx = swri%sphbessels(species_number,sbidx)%qidx

       if (lval == -1) exit
       if (lval < swexside_filter%min_l) cycle
       if (lval > swexside_filter%max_l) cycle
       if (qidx > swexside_filter%max_q) cycle

       ! qoh: At rad = 0 only consider l=0 (we have 0/0 for other values of l)
       ! jd:  Mathematica seems to indicate all swpots evaluate to 0 for l>0
       !      when rad is 0, so it should be safe just to ignore them, only
       ! [*]  making sure to fast forward to the right n_sws and store 0's.

       if(rad > SAFE_DIV_EPS) then
          ! jd: Usual case
          invrad = 1.0_DP / rad
          unit_disp = invrad * disp
          bessint = swri_sph_bess_pot_int(rad, invrad, &
               swri%sphbessels(species_number, sbidx))
          factor = common_factor / real(((2*lval+1) * (-1)**lval),kind=DP)
          do mval = -lval, lval
             sw_pot(sw_idx) = bessint * factor * &
                  sw_real_sph_harm_unit(unit_disp%x, unit_disp%y, unit_disp%z, &
                  lval, mval)
             sw_idx = sw_idx + 1
          end do
       else
          ! jd: Corner case of SWpot centre on gridpoint where we evaluate
          if(lval == 0) then
             ! z_00 = 0.5_DP/sqrt(pi) = 0.282094791773878_DP
             bessint = 1.0_DP / (swri%sphbessels(species_number,sbidx)%qval**2) &
                  + swri%sphbessels(species_number,sbidx)%nearpotint
          else
             bessint = 0.0_DP ! [*]
          end if
          do mval = -lval, lval
             sw_pot(sw_idx) = bessint * common_factor * 0.282094791773878_DP
             sw_idx = sw_idx + 1
          end do
       end if
    end do

  end subroutine swri_calc_swpots_at_point

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_calc_swops_at_point(swri, cell, swop_val, curpoint, centrex, &
       swexside_filter, pbc, sw_or_swpot, rcut)
    !==================================================================!
    ! A simple wrapper around calc_sws_at_point() and                  !
    ! calc_swpots_at_point().                                          !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2017.                       !
    !==================================================================!

    use geometry, only: POINT
    use utils, only: utils_assert
    use simulation_cell, only: CELL_INFO

    implicit none

    ! jd: Arguments
    type(SW_RI),       intent(in)           :: swri
    type(SW_QUALITY),  intent(in)           :: swexside_filter
    type(CELL_INFO),   intent(in)           :: cell
    real(kind=DP),     intent(out)          :: swop_val(&
         swexside_filter%num_sws_per_centre)
    type(POINT),       intent(in)           :: curpoint
    type(ATOM_CENTRE), intent(in)           :: centrex
    character,         intent(in)           :: sw_or_swpot
    logical,           intent(in)           :: pbc
    real(kind=DP),     intent(in), optional :: rcut


    ! jd: Local variables

    ! ----------------------------------------------------------------------

    if(sw_or_swpot == 'S') then
       call swri_calc_sws_at_point(swri, cell, swop_val, curpoint, &
            centrex, swexside_filter, pbc)
    elseif(sw_or_swpot == 'P') then
       call swri_calc_swpots_at_point(swri, cell, swop_val, curpoint, &
            centrex, swexside_filter, pbc, rcut)
    else
       call utils_assert(.false., &
            'swri_calc_swops_at_point: Unrecognized sw_or_swpot:', sw_or_swpot)
    end if

  end subroutine swri_calc_swops_at_point

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_obtain_swops_at_point(swri, cell, swops_at_point, curpoint, &
       src_centre, src_atom, sw_or_swpot, swexside_filter, stash, rcut)
    !==================================================================!
    ! Retrieves from cache, or if absent, recalculates SWOPs at a      !
    ! given point. This should only be used for calculations on the    !
    ! fine grid, as on the coarse grid the PPD version is much more    !
    ! efficient.                                                       !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !------------------------------------------------------------------!
    ! Caveats:                                                         !
    !   The caller should call hash_table_cleanup_at_will() on the     !
    !   hash table holding cached SWOPs after calling this routine to  !
    !   ensure user-set memory requirements are not exceeded.          !
    !==================================================================!

    use constants, only: PI
    use geometry, only: POINT, operator(.DOT.)
    use hash_table, only: hash_table_lookup, hash_table_stash
    use rundat, only: pub_debug, pub_hfx_bc_is_periodic
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort, utils_assert, utils_point_to_str, &
         utils_real_to_str, utils_int_to_str

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout)   :: swri
    type(SW_QUALITY), intent(in) :: swexside_filter
    type(CELL_INFO), intent(in)  :: cell
    real(kind=DP), intent(out)   :: swops_at_point(&
         swexside_filter%num_sws_per_centre)
    type(POINT), intent(in)      :: curpoint
    type(ATOM_CENTRE), intent(in):: src_centre
    integer, intent(in)          :: src_atom
    character, intent(in)        :: sw_or_swpot
    real(kind=DP), intent(inout) :: stash(:)
    real(kind=DP), intent(in), optional :: rcut

    ! jd: Local variables
    logical :: filtering_needed
    logical :: pbc
    integer :: points_in_bessel
    integer :: num_swri_sws
    integer :: src_offset, dest_offset
    integer :: lval, sbidx, qidx
    integer :: size_of_cached_data
    real(kind=DP) :: fgp1, fgp2, fgp3
    integer :: igp1, igp2, igp3
    real(kind=DP), parameter :: invpi = 1.0_DP/PI
    real(kind=DP), parameter :: eps = 1D-9
    real(kind=DP) :: swri_swops_at_point_workspace(&
         swri%quality%num_sws_per_centre)

    ! ----------------------------------------------------------------------

    filtering_needed = .not. (swri%quality%min_l == swexside_filter%min_l .and. &
         swri%quality%max_l == swexside_filter%max_l .and. &
         swri%quality%max_q == swexside_filter%max_q)

    if(.not. any(pub_hfx_bc_is_periodic)) then
       pbc = .false.
    else if(all(pub_hfx_bc_is_periodic)) then
       pbc = .true.
    else
       call utils_abort('swri_obtain_swops_at_point: &
            &Mixed BCs are not currently supported.')
    end if

    ! jd: Figure out the fine grid point indices so that we have an integer
    !     key for the hash table. 1/pi and not 1/2pi accounts for double grid.
    fgp1 = invpi * real(cell%total_pt1,kind=DP) * (curpoint.DOT.cell%b1)
    fgp2 = invpi * real(cell%total_pt2,kind=DP) * (curpoint.DOT.cell%b2)
    fgp3 = invpi * real(cell%total_pt3,kind=DP) * (curpoint.DOT.cell%b3)
    igp1 = nint(fgp1)
    igp2 = nint(fgp2)
    igp3 = nint(fgp3)
    if(pub_debug) then ! jd: Ensure this does indeed lie on a (fine) grid point
       if(abs(fgp1 - real(igp1, kind=DP)) > eps .or. &
            abs(fgp2 - real(igp2, kind=DP)) > eps .or. &
            abs(fgp3 - real(igp3, kind=DP)) > eps) then
          call utils_assert(.false., &
               'swri_obtain_swops_at_point: curpoint does not look like a &
               &fine grid point', 'point: '//&
               trim(utils_point_to_str(curpoint,'f20.9'))//&
               ', fgp1: '//trim(utils_real_to_str(fgp1))//&
               ', fgp2: '//trim(utils_real_to_str(fgp2))//&
               ', fgp3: '//trim(utils_real_to_str(fgp3))//&
               ', igp1: '//trim(utils_int_to_str(igp1))//&
               ', igp2: '//trim(utils_int_to_str(igp2))//&
               ', igp3: '//trim(utils_int_to_str(igp3)))
       end if
    end if
! --- jd: No OMP CRITICAL needed -- using stash to effect writes ---

    if(filtering_needed) then
       ! jd: Will have to filter, look up into workspace
       call hash_table_lookup(swri_swops_at_point_workspace, size_of_cached_data,&
            swri%swops_at_points_ht, key1 = src_atom, &
            key2 = igp1, key3 = igp2, key4 = igp3, key5 = ichar(sw_or_swpot))
    else
       ! jd: No filtering, look up straight into result
       call hash_table_lookup(swops_at_point, size_of_cached_data, &
            swri%swops_at_points_ht, key1 = src_atom, &
            key2 = igp1, key3 = igp2, key4 = igp3, key5 = ichar(sw_or_swpot))
    end if

! --- jd: No OMP CRITICAL needed -- using stash to effect writes ---

    if(size_of_cached_data /= -1) then
       ! Cache HIT
       num_swri_sws = size_of_cached_data
    else
       ! Cache MISS
       ! - calculate the SWOPs for this PPD. Don't filter at this stage, as
       !   we want to cache the swriside version.

       if(filtering_needed) then
          ! jd: Will have to filter, calculate into workspace
          call swri_calc_swops_at_point(swri, cell, &                      ! in
               swri_swops_at_point_workspace, &                            ! out
               curpoint, src_centre, swri%quality, pbc, sw_or_swpot, rcut) ! in
       else
          ! jd: No filtering, calculate straight into result
          call swri_calc_swops_at_point(swri, cell, &                      ! in
               swops_at_point, &                                           ! out
               curpoint, src_centre, swri%quality, pbc, sw_or_swpot, rcut) ! in
       end if
       num_swri_sws = swri%quality%num_sws_per_centre
       ! - add these to the cache
! --- jd: No OMP CRITICAL needed -- hash_table_stash() takes care of that ---
       if(filtering_needed) then
          call hash_table_stash(stash, swri%swops_at_points_ht, &
               swri_swops_at_point_workspace, &
               num_swri_sws, key1 = src_atom, key2 = igp1, key3 = igp2, &
               key4 = igp3, key5 = ichar(sw_or_swpot))
       else
          call hash_table_stash(stash, swri%swops_at_points_ht, swops_at_point,&
               num_swri_sws, key1 = src_atom, key2 = igp1, key3 = igp2, &
               key4 = igp3, key5 = ichar(sw_or_swpot))
       end if
! --- jd: No OMP CRITICAL needed -- hash_table_stash() takes care of that ---
    end if

    if(num_swri_sws /= swri%quality%num_sws_per_centre) then
       call utils_abort('swri_obtain_swops_at_point: Logic error -- disagree&
            &ment about num_swri_sws',num_swri_sws, &
            swri%quality%num_sws_per_centre)
    end if

    ! jd: Filter the obtained SWOPs, if needed
    dest_offset = 1
    src_offset = 1
    if(filtering_needed) then
       do sbidx = 1, swri%max_num_bessels
          lval = swri%sphbessels(1,sbidx)%lval ! @IDENTICAL_RADII
          qidx = swri%sphbessels(1,sbidx)%qidx ! @IDENTICAL_RADII
          points_in_bessel = (2*lval+1)

          if (lval == -1) exit
          if (lval < swexside_filter%min_l .or. lval > swexside_filter%max_l) then
             ! skip to the next Bessel, i.e. by 2l+1 PPDs
             src_offset = src_offset + points_in_bessel
             cycle
          end if
          if (qidx > swexside_filter%max_q) then
             ! skip to the next Bessel, i.e. by 2l+1 PPDs
             src_offset = src_offset + points_in_bessel
             cycle
          end if

          ! jd: Copy over the next Bessel, i.e. 2l+1 PPDs
          swops_at_point(dest_offset:dest_offset+points_in_bessel-1) = &
               swri_swops_at_point_workspace(src_offset:src_offset+points_in_bessel-1)
          dest_offset = dest_offset + points_in_bessel
          src_offset = src_offset + points_in_bessel

       end do
    end if

  end subroutine swri_obtain_swops_at_point

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function swri_qualities_equal(quality1, quality2)
    !==================================================================!
    ! Returns .true. if two SW_QUALITY objects are equal (that is,     !
    ! all their components are identical). Otherwise returns .false..  !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   quality1 (in): The first SW_QUALITY object.                    !
    !   quality2 (in): The second SW_QUALITY object.                   !
    !------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2019.                        !
    !==================================================================!

    implicit none

    ! jd: Arguments
    type(SW_QUALITY), intent(in) :: quality1
    type(SW_QUALITY), intent(in) :: quality2

    ! -------------------------------------------------------------------------

    swri_qualities_equal = &
         (quality1%min_l == quality2%min_l .and. &
         quality1%max_l == quality2%max_l .and. &
         quality1%max_q == quality2%max_q .and. &
         quality1%num_sws_per_centre == quality2%num_sws_per_centre .and. &
         quality1%max_sws_per_centre == quality2%max_sws_per_centre)

  end function swri_qualities_equal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_init_sw_radial_lookup(swri, cell)
    !==================================================================!
    ! Initialises the radial lookup table for the radial component of  !
    ! SWpots. This works similarly to the openbc localpseudo lookup,   !
    ! except it doesn't use services_1D_interpolation() but linear     !
    ! interpolation between points for efficiency.                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   swri (inout): Holds the rad_lookup(:,:) array that is          !
    !                 populated here.                                  !
    !   cell (in): Simulation cell, to determine the maximum distance  !
    !              for the lookup (cell diagonal).                     !
    !------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2020.                      !
    !==================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_allgather
    use constants, only: SAFE_DIV_EPS
    use geometry, only: magnitude, operator(+), operator(-), operator(*)
    use rundat, only: pub_hfx_bessel_rad_nptsx, pub_threads_max, pub_rootname, &
         pub_debug_on_root
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_unit

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(inout)  :: swri
    type(CELL_INFO), intent(in) :: cell

    ! jd: Local variables
    integer :: npts_x                    ! jd: Size of the radial realspace grid
    integer :: npts_x_per_core           ! jd: One core's share of the above
    integer :: my_npts_x                 ! jd: One core's share of the above
    integer :: npts_x_rem                ! ndmh: remainder of points/cores
    integer :: species, x_rad_pt         ! jd: Indices
    real(kind=DP) :: hx                  ! jd: Gridpoint separation
    real(kind=DP) :: xmax                ! jd: Last x stored in lookups
    real(kind=DP) :: rad                 ! jd: Current x
    real(kind=DP) :: invrad              ! jd: 1/x if safe to divide
    integer :: x1, x2                    ! jd: min and max x for this core
    integer :: species_number
    integer :: sbidx
    integer :: lval
    real(kind=DP), allocatable :: rad_lookup_local(:,:)
    integer :: ierr
    integer :: output_unit
    character(len=2048) :: output_format
    character(len=*), parameter :: myself = 'swri_init_sw_radial_lookup'

    ! -------------------------------------------------------------------------

    npts_x = pub_hfx_bessel_rad_nptsx

    ! jd: Find out the maximum value of x which we might need, determine the
    !     fineness of the radial grid depending on this. Using the diagonal of
    !     the cell is overkill -- technically it would be sufficient to use
    !     (exchange cutoff + r_NGWF), but this gets complicated by the fact that
    !     we calculate potentials in PPDs.
    xmax = magnitude(cell%a1 + cell%a2 + cell%a3)
    hx = xmax/(npts_x-1)

    ! jd: Split the work across cores: every core will work with points x1..x2
    npts_x_per_core = npts_x / pub_total_num_procs
    npts_x_rem = mod(npts_x,pub_total_num_procs)
    my_npts_x = npts_x_per_core
    if (pub_my_proc_id<npts_x_rem) my_npts_x = my_npts_x + 1
    x1 = pub_my_proc_id*npts_x_per_core + min(pub_my_proc_id,npts_x_rem) + 1
    x2 = (pub_my_proc_id+1)*npts_x_per_core + min(pub_my_proc_id + 1,npts_x_rem)

    ! jd: Sanity check against a situation where the last core(s) is/are left
    !     with no work
    if(x1 > npts_x .or. x2 > npts_x) then
       call utils_abort('Too many cores to sensibly divide &
            &hfx_bessel_rad_nptsx radial points in the HFx Bessel radial &
            &lookup in '//trim(myself)//' Increase the &
            &number of points or reduce the number of cores.')
    end if

    ! jd: Allocate and zero lookup array.
    allocate(swri%rad_lookup(npts_x,swri%max_num_bessels),stat=ierr)
    call utils_alloc_check(myself,'swri%rad_lookup',ierr)
    swri%rad_lookup(:,:) = 0.0_DP

    ! jd: Allocate the part of the lookup local to this core
    allocate(rad_lookup_local(my_npts_x,swri%max_num_bessels),stat=ierr)
    call utils_alloc_check(myself,'rad_lookup_local', ierr)
    rad_lookup_local = 0.0_DP

    species_number = 1 !@IDENTICAL_RADII

    ! jd: Go over all of the x's that belong to this core
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(x_rad_pt,rad,invrad,sbidx,lval) &
!$OMP SHARED(rad_lookup_local,x1,x2,hx,swri,species_number,pub_threads_max)
    do x_rad_pt = x1, x2
       rad = (x_rad_pt-1) * hx
       if(rad > SAFE_DIV_EPS) then
          ! jd: Usual case -- point not practically exactly on atom centre
          invrad = 1.0_DP / rad

          do sbidx = 1, swri%max_num_bessels
             lval = swri%sphbessels(species_number,sbidx)%lval
             if (lval == -1) exit ! jd: signals the end

             rad_lookup_local(x_rad_pt-x1+1,sbidx) = &
                  swri_sph_bess_pot_int(rad, invrad, &
                  swri%sphbessels(species_number,sbidx))
          end do

       else
          ! jd: Corner case of SWpot centre on gridpoint where we evaluate
          do sbidx = 1, swri%max_num_bessels
             lval = swri%sphbessels(species_number,sbidx)%lval
             if (lval == -1) exit ! jd: signals the end

             if(lval == 0) then
                rad_lookup_local(x_rad_pt-x1+1,sbidx) = &
                     1.0_DP / (swri%sphbessels(species_number,sbidx)%qval**2) &
                     + swri%sphbessels(species_number,sbidx)%nearpotint
             else
                rad_lookup_local(x_rad_pt-x1+1,sbidx) = 0.0_DP
             end if
          end do
       end if

    end do
!$OMP END PARALLEL DO

    ! jd: Allgather the partial results rad_lookup_local -> rad_lookup
    do sbidx = 1, swri%max_num_bessels
       call comms_allgather(swri%rad_lookup(:,sbidx),rad_lookup_local(:,sbidx),&
            my_npts_x, (/-1/), (/-1/))
    end do

    ! jd: Free the part of the lookup local to this core
    deallocate(rad_lookup_local,stat=ierr)
    call utils_dealloc_check(myself,'rad_lookup_local', ierr)

    ! jd: Output the interpolation to a file if in debug mode
    if(pub_debug_on_root) then
       output_format = '('//repeat('f20.15,',swri%max_num_bessels)//'f20.15)'
       do x_rad_pt = 1, npts_x

          rad = (x_rad_pt-1) * hx

          if(x_rad_pt == 1) then
             output_unit = utils_unit()
             open(unit = output_unit, file = trim(pub_rootname)//&
                  '.bessel_rad_interp', action = "write", err=10)
          end if

          ! jd: Output the calculated pseudopotential
          write(output_unit,output_format) rad, swri%rad_lookup(x_rad_pt,:)

          if(x_rad_pt == npts_x) then
             close(output_unit,err=20)
          end if

       end do ! radial points
    end if

    return

10  call utils_abort('Could not open a debug file to write the&
         & Bessel radial interpolation to.')
20  call utils_abort('Could not close a debug file with the&
         & Bessel radial interpolation.')

  end subroutine swri_init_sw_radial_lookup

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_read_block(swri, par, elements, matchar, atomblock, &
       atoma, atomb, error)
    !========================================================================!
    ! Reads a V/O matrix atomblock from a file.                              !
    ! The silent assumption is every proc can do I/O.                        !
    ! If 'error' is absent, aborts on error.                                 !
    ! If 'error' is present, returns 0 if succeeded, or >0 if failed.        !
    ! Tentative list of error meanings:                                      !
    ! 1 - error opening file.                                                !
    ! 2 - error reading file.                                                !
    ! 3 - error closing file.                                                !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    !   swri (in): SW_RI where all this happens.                             !
    !   par (in): Needed to get orig_atom.                                   !
    !   matchar (in):    'V' or 'O'.                                         !
    !   atomblock (out):  The atomblock to read into.                        !
    !                     If an error occured, garbage_real is returned in   !
    !                     all elements.                                      !
    !   atoma, atomb (in): Indices to the block in SFC global scheme.        !
    !   error (in, opt): See above.                                          !
    !========================================================================!

    use constants, only: garbage_real
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_swri_verbose
    use utils, only: utils_unit, utils_abort, utils_assert

    implicit none

    ! jd: Arguments
    type(SW_RI), intent(in)      :: swri
    type(PARAL_INFO), intent(in) :: par
    type(ELEMENT), intent(in)    :: elements(par%nat)
    character, intent(in)        :: matchar
    real(kind=DP), intent(out)   :: atomblock(swri%quality%max_sws_per_centre, &
        swri%quality%max_sws_per_centre)
    integer, intent(in)          :: atoma, atomb
    integer, intent(out), optional :: error

    ! jd: Local variables
    integer :: blk_unit
    character(len=1024) :: filename
    integer :: ierr

    ! --------------------------------------------------------------------

    call swri_pos_to_filename(filename, &
         elements(par%orig_atom(atoma))%centre, &
         elements(par%orig_atom(atomb))%centre, matchar, swri%swri_name)

    ! Open
    blk_unit = utils_unit()
    open(blk_unit, file=filename, form='unformatted', &
             action='read', iostat=ierr, status='old')
    if(ierr /= 0) then
       if(present(error)) then
          error = 1
          if(pub_swri_verbose) then
             write(stdout,'(a,i0,a,a)') &
                  '! Open error (', ierr,') for ', trim(filename)
          end if
          atomblock(:,:) = garbage_real
          return
       else
          call utils_abort(&
               'Error in swri_read_block: &
               &opening file "'//trim(filename)//'" failed with code ', ierr)
       end if
    end if

    ! Read
    read(blk_unit, iostat=ierr) atomblock
    if(ierr /= 0) then
       if(present(error)) then
          error = 2
          if(pub_swri_verbose) then
             write(stdout,'(a,i0,a,a)') &
                  '! Read error (', ierr,') for ', trim(filename)
          end if
          atomblock(:,:) = garbage_real
          return
       else
          call utils_abort(&
               'Error in swri_read_block: &
               &reading file "'//trim(filename)//'" failed with code ', ierr)
       end if
    end if

    ! Close
    close(blk_unit, iostat=ierr)
    if(ierr /= 0) then
       if(present(error)) then
          error = 3
          if(pub_swri_verbose) then
             write(stdout,'(a,i0,a,a)') &
                  '! Close error (', ierr,') for ', filename
          end if
          atomblock(:,:) = garbage_real
          return
       else
          call utils_assert(ierr == 0, &
               'Error in swri_read_block: &
               &closing file "'//trim(filename)//'" failed with code ', ierr)
       end if
    end if

    if(present(error)) error = 0
    if(pub_swri_verbose) then
       write(stdout,'(a,a)') &
            '! Success for ', trim(filename)
    end if

  end subroutine swri_read_block

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine swri_pos_to_filename(filename, posa, posb, matchar, swri_name)
    !========================================================================!
    ! Determines the atomblock filename from the positions of the atoms and  !
    ! the character ('V' or 'O') describing the metric.                      !
    ! The filename is built from the positions of the atoms rounded to 1 fig !
    ! after the decimal point. 'int' rounding is used. All filenames are     !
    ! prefixed with 'swri_assembly_prefix', which helps in the scenario where!
    ! many runs use the same atomblocks, located in one directory. The       !
    ! default is just pub_rootname.                                          !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    !   filename (out):  The returned filename.                              !
    !   posa, posb (in): The positions of the two atoms.                     !
    !   matchar (in): 'V' or 'O', denoting the metric used.                  !
    !   swri_name (in): The name of the swri.                                !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2016.04.                                  !
    !========================================================================!

    use geometry, only: POINT
    use rundat, only: pub_swri_assembly_prefix

    implicit none

    ! jd: Arguments
    character(len=*), intent(out) :: filename
    type(POINT), intent(in)       :: posa, posb
    character, intent(in)         :: matchar
    character(len=*), intent(in)  :: swri_name

    ! jd: Local variables
    integer :: posax_intpart, posbx_intpart
    integer :: posax_fracpart, posbx_fracpart
    integer :: posay_intpart, posby_intpart
    integer :: posay_fracpart, posby_fracpart
    integer :: posaz_intpart, posbz_intpart
    integer :: posaz_fracpart, posbz_fracpart

    ! -----------------------------------------------------------------------

    posax_intpart =  int(posa%x)
    posax_fracpart = int(abs(posa%x - real(posax_intpart,kind=DP)) * 1D5)
    posay_intpart =  int(posa%y)
    posay_fracpart = int(abs(posa%y - real(posay_intpart,kind=DP)) * 1D5)
    posaz_intpart =  int(posa%z)
    posaz_fracpart = int(abs(posa%z - real(posaz_intpart,kind=DP)) * 1D5)

    posbx_intpart =  int(posb%x)
    posbx_fracpart = int(abs(posb%x - real(posbx_intpart,kind=DP)) * 1D5)
    posby_intpart =  int(posb%y)
    posby_fracpart = int(abs(posb%y - real(posby_intpart,kind=DP)) * 1D5)
    posbz_intpart =  int(posb%z)
    posbz_fracpart = int(abs(posb%z - real(posbz_intpart,kind=DP)) * 1D5)


    write(filename,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5,a1,i5.5,a1,i5.5,a2,&
         &i5.5,a1,i5.5,a1,i5.5,a1,i5.5,a1,i5.5,a1,i5.5,a)') &
         trim(pub_swri_assembly_prefix)//'.', &
         posax_intpart, '.', posax_fracpart, '_', &
         posay_intpart, '.', posay_fracpart, '_', &
         posaz_intpart, '.', posaz_fracpart, '__', &
         posbx_intpart, '.', posbx_fracpart, '_', &
         posby_intpart, '.', posby_fracpart, '_', &
         posbz_intpart, '.', posbz_fracpart, &
         "."//trim(swri_name)//"."//matchar//'matrixblock'

  end subroutine swri_pos_to_filename

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module sw_resolution_of_identity
