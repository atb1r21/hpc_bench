! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                          Remote module                         !
!                                                                !
! This module contains subroutines for high-level communication  !
! of data between procs in a parallel, distributed-memory system,!
! (such as using MPI), where the data is cached locally to reduce!
! future comms.                                                  !
!                                                                !
! Said data currently includes:                                  !
! - NGWFs (multiple sets are now supported),                     !
! - DKN atomblocks,                                              !
! - generic SPAM3 matrix atomblocks,                             !
! - SWx expansion coefficients ('coeffs') -- in two versions,    !
!   one spin-dependent, and one spin-agnostic.                   !
!----------------------------------------------------------------!
! Written in January 2014 by Jacek Dziedzic.                     !
! Extended in February 2015 by Jacek Dziedzic to include a       !
! library for NGWF datasets.                                     !
! Extended in September 2016 by Jacek Dziedzic to work with      !
! complex NGWFs.                                                 !
! Extended in October 2016 by Jacek Dziedzic to work with multi- !
! ple NGWF sets, also _noadd variants added.                     !
! Extended in December 2016-January 2017 by Jacek Dziedzic with  !
! pack_data_real, and unpack_ngwf_to_fftbox.                     !
! Extended in June 2017 by Jacek Dziedzic with 'user_tag'.       !
! Extended in June 2017 by Jacek Dziedzic with the ability to    !
! exchange spin-dependent SW coeffs in addition to spin-agnostic !
! SW coeffs.                                                     !
! Extended in July/August 2018 by Jacek Dziedzic to fully support!
! mixed NGWF sets.                                               !
! Simplified in June 2019 by Jacek Dziedzic to disentangle DKN   !
! exchange from NGWF exchange.                                   !
!================================================================!
!                                                                !
! Basic usage:                                                   !
! ------------                                                   !
! DKN blocks are easy, obtain_dkn_block() is used to get hold of !
! a DKN atomblock. Mind that the block is returned by pointer.   !
! If the block is in local cache (hash table), the block will be !
! returned from there. If not, it is checked if the block is     !
! MPI-rank-local. If so, it will be extracted from the SPAM3     !
! structure, added to local cache and the pointer returned.      !
! Finally, if the block is not cached and is not local, a request!
! will be sent to the owner proc, and the subroutine will only   !
! return once the request is completed. The block will be added  !
! to cache, and the pointer returned.                            !
! Caveats:                                                       !
! 1) DKN block comms rely on a single, module-wide persistent    !
!    buffer that can hold one DKN block. Thus only use one OMP   !
!    thread for block comms, and use or copy out the DKN block   !
!    returned by pointer, as the underlying buffer becomes       !
!    volatile once another obtain_dkn_block() is initiated on    !
!    the same MPI rank. This can be awkward, but returning by    !
!    pointer helps with efficiency.                              !
! 2) Ensure you keep calling serve_dkn_blocks() on all MPI ranks !
!    that can own DKN blocks you're asking for. Examine the      !
!    'done' argument to ascertain if all DKN requests on all     !
!    ranks completed.                                            !
! 3) The returned pointer can point to data in a hash table.     !
!    Be aware of when pointers to hash tables get invalidated.   !
!    See the documentation of hash_table_mod for more info.      !
!                                                                !
! NGWFs are not that hard either. The general idea is the same,  !
! you call obtain_ngwf(), and all procs must serve NGWFs until   !
! 'done' is .true. (all requests satisfied). There are a few     !
! differences:                                                   !
! a) NGWFs are communicated 'packed', that is as arrays of DP    !
!    reals (even for complex NGWFs!). A packed NGWF contains not !
!    only the NGWF data on grid, but also the information        !
!    normally contained in its SPHERE too (centre, radius, etc)  !
!    and its TIGHTBOX info.                                      !
!    All you ever wanted to know about an NGWF in one place, and !
!    easy to send, receive, stash in a hash table, serialize...  !
!    There are routines to pack and unpack NGWFs.                !
! b) NGWFs can be obtained by pointer (like DKN blocks, with the !
!    same caveats), or by copy (obtain_ngwf_by_copy()). The      !
!    latter is less efficient, but free of the caveats.          !
! c) '_noadd' variants are provided for contexts, where the cache!
!    (hash table) is not writable -- e.g. when it is a member of !
!    an NGWF_REP that is intent(in) in a gradient calculation.   !
!    In this case no NGWFs are, of course, added to the cache.   !
! d) There are *separate* persistent NGWF buffers for each NGWF  !
!    set, ie. there is an array of MAX_NGWF_CACHES buffers.      !
!    Currently MAX_NGWF_CACHES is 6, for:                        !
!    1: valence, 2: conduction, 3: joint, 4: auxiliary,          !
!    5: jointPH (whatever this is), 6: vacuum state NGWFs in     !
!    polarisable embedding.                                      !
! e) When using mixed NGWFs (conduction, etc) and exchanging     !
!    them simultaneously, both sets must be provided, so that    !
!    deadlocks are avoided when serving. Each set is             !
!     distinguished with a 'cache_handle', and 'user_tag'        !
!    selects the set whose NGWFs are to be obtained.             !
!                                                                !
! Generic SPAM3 matrices can be communicated using remote_       !
! {obtain,serve}_matrix_blocks(). This works similarly to DKN    !
! exchange.                                                      !
!                                                                !
! Finally, this module allows for exchanging SW expansion coeffi-!
! cients ('coeffs'). For the types of coefficients that can be   !
! exchanged, see remote_serve_coeffs().                          !
!                                                                !
!----------------------------------------------------------------!


module remote

  use constants, only: DP, CRLF
  use hash_table, only: HT_HASH_TABLE
  use timer, only: timer_clock

  implicit none

  private

  ! jd: Public subroutines
  public :: remote_init
  public :: remote_exit

  public :: remote_ppd_list_of_atom
  public :: remote_n_ppds_of_atom

  public :: remote_pack_ngwf
  public :: remote_pack_data_real
  public :: remote_unpack_ngwf
  public :: remote_unpack_ngwf_no_sphere
  public :: remote_unpack_ngwf_to_fftbox

  public :: remote_obtain_ngwf
  public :: remote_obtain_ngwf_noadd
  public :: remote_obtain_ngwf_by_copy
  public :: remote_obtain_ngwf_by_copy_noadd
  public :: remote_serve_ngwfs
  public :: remote_obtain_unpacked_ngwfs

  public :: remote_obtain_dkn_block
  public :: remote_serve_dkn_blocks

  public :: remote_obtain_matrix_block
  public :: remote_serve_matrix_blocks

  public :: remote_obtain_coeffs
  public :: remote_serve_coeffs

  ! jd: Public state
  public :: remote_hash_table_info_unit

  ! --------------------------------------------------------------------------

  ! jd: Tags for comms
  integer, parameter, public :: NGWF_REQUEST_TAG      = 10000
  integer, parameter, public :: DKN_BLOCK_REQUEST_TAG = 20000
  integer, parameter, public :: COEFFS_REQUEST_TAG    = 30000
  integer, parameter, public :: MATRIX_BLOCK_REQUEST_TAG = 40000
  integer, parameter :: NGWF_DATA_TAG         = 100000
  integer, parameter :: DKN_BLOCK_DATA_TAG    = 200000
  integer, parameter :: COEFFS_DATA_TAG       = 300000
  integer, parameter :: MATRIX_BLOCK_DATA_TAG = 400000

  ! jd: Length (in doubles) of a packed representation of NGWF. This is
  !     determined by remote_init_packed_ngwf_size().
  integer, save, public      :: packed_ngwf_size = -1
  integer, parameter, public :: packed_ngwf_header_length = 29 ! [*]

  ! jd: Every 10 times an NGWF is retrieved from cache, take care to
  !     serve local NGWFs to remote processes, so that they don't starve
  !     once we can satisfy all local requests from cache. Increasing this value
  !     makes local NGWF accesses more efficient, but makes remote procs wait
  !     longer for requests issued towards the end.
  integer, parameter :: ngwf_serve_period = 10

  ! jd: Length (in doubles) of a packed representation of coefficients sent by
  !     remote_serve_coeffs and received by remote_obtain_coeffs. This should
  !     be at least 3 + (max_num_ngwfs_per_atom * max_num_ngwfs_per_atom *
  !     max_num_sws_in_expansion). This only dimensions the *local* buffer in
  !     send/recv routines, and the actual length of sent data is used to
  !     dimension the network comms itself, so it's OK to be generous, up to a
  !     point, as this will be stack-based.
  !     With lmax=4, qmax=25 we get max_num_sws_in_expansion=625.
  !     4ZP2P1P calculations use 33 NGWFs per atom. 625*33*33 = 680625.
  !     Using 1048576 as a generous overestimate. 8M should be fine even if it
  !     winds up on a per-thread stack.
  integer, parameter :: MAX_PACKED_COEFFS_SIZE = 1048576

  ! jd: Maximum size (in doubles) of any DKN block for a pair of atoms for
  !     (potentially) both spins. Essentially an compile-time upper bound for
  !     ngwf_basis%max_on_atom * ngwf_basis%max_on_atom * num_spins.
  !     This is only used to dimension a module-wide workspace array, so
  !     it's harmless to make it slightly bigger. What is sent, received and
  !     stored in hash tables is dimensioned with *actual* block sizes, so
  !     no worries. Asserts are in place to guard against this value being
  !     too small at runtime.
  integer, parameter :: MAX_DKN_BLOCK_SIZE = 33*33*2

  ! jd: The same for any exchanged matrix blocks (currently metric matrix)
  integer, parameter :: MAX_MATRIX_BLOCK_SIZE = 500*500

  real(kind=DP), target, save :: persistent_packed_dkn_blk(MAX_DKN_BLOCK_SIZE)
  real(kind=DP), target, save :: persistent_packed_matrix_blk(MAX_MATRIX_BLOCK_SIZE)
  real(kind=DP), target, save, allocatable :: persistent_packed_ngwf(:,:)
  ! jd: idx1: points in packed NGWF. idx2: 1:MAX_NGWF_CACHES.........^^^

  ! jd: Every 10 times a DKN block is retrieved from cache, take care to
  !     serve local DKN blocks to remote processes, so that they don't starve
  !     once we can satisfy all local requests from cache. Increasing this value
  !     makes local DKN accesses more efficient, but makes remote procs wait
  !     longer for requests issued towards the end.
  integer, parameter :: dkn_block_serve_period = 10

  ! jd: Ditto for matrix block comms
  integer, parameter :: matrix_block_serve_period = 10

  ! jd: Hash table for NGWF cache
  integer, parameter, public                :: MAX_N_NGWF_CACHES = 6
  integer, save                             :: remote_hash_table_info_unit
  character(len=256), save                  :: hash_table_info_filename
  type(HT_HASH_TABLE), save, public, target :: remote_ngwf_cache_hts(MAX_N_NGWF_CACHES)

  integer, parameter                        :: maxslots_ngwf_cache = 20000
  ! ^ This is only the number of slots, the actual HT is permitted to overfill.
  !   Overhead is tiny (2.5MB for 20000 x MAX_N_NGWF_CACHES), so we can be very
  !   generous and it does not need to be system-size dependent. Even for the
  !   largest systems (3000 atoms on 32 ranks) we still need only ~6000 slots
  !   (the HT itself is then ~1.5GB).

  character(len=32), parameter :: hash_table_info_rootname = 'hash_table_remote'

  logical, save :: initialised = .false.

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_pack_ngwf(packed_ngwf, ngwfs_on_grid, ngwf_basis, &
       local_ngwf_idx)
    !==========================================================================!
    ! Packs information about an NGWF (its SPHERE, corresponding parts of      !
    ! its FUNCTIONS, its tightbox info) into a packed_ngwf -- a 1D array of    !
    ! reals that is easier to transfer between procs and is conveniently real  !
    ! regardless of whether real or complex NGWFs are used.                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_ngwf (out):   Stores the resultant packed NGWF.                 !
    !   ngwfs_on_grid (in):  Usual FUNCTIONS representation. Data will be      !
    !                        extracted from an offset that will be determined  !
    !                        automatically.                                    !
    !   ngwf_basis (in):     The SPHERE is extracted from here, and so is the  !
    !                        number of points per PPD.                         !
    !   local_ngwf_idx (in): Proc-local index of the NGWF in the ngwf_basis.   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2013.                              !
    ! Adapted to complex NGWFs by Jacek Dziedzic in September 2016.            !
    ! Extended to include tightbox info by Jacek Dziedzic in January 2017.     !
    !==========================================================================!

    use datatypes, only: FUNCTIONS
    use function_basis, only: FUNC_BASIS
    use utils, only: utils_assert, utils_logical_to_real

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)   :: packed_ngwf(packed_ngwf_size)
    type(FUNCTIONS), intent(in)  :: ngwfs_on_grid
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    integer, intent(in)          :: local_ngwf_idx

    ! jd: Local variables
    integer :: i, j, dest, last_index, src_offset
    integer :: n_pts_tot
    integer :: n_ppds
    character(len=*), parameter :: myself = 'remote_pack_ngwf'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    call utils_assert(packed_ngwf_size /= -1, trim(myself)//&
         ': remote_init() must be called first.')

    n_ppds = ngwf_basis%spheres(local_ngwf_idx)%n_ppds_sphere

    packed_ngwf(1) = ngwf_basis%spheres(local_ngwf_idx)%centre%x  ! [*] begin
    packed_ngwf(2) = ngwf_basis%spheres(local_ngwf_idx)%centre%y
    packed_ngwf(3) = ngwf_basis%spheres(local_ngwf_idx)%centre%z
    packed_ngwf(4) = ngwf_basis%spheres(local_ngwf_idx)%radius
    packed_ngwf(5) = real(n_ppds,kind=DP)
    packed_ngwf(6) = real(ngwf_basis%n_pts)
    packed_ngwf(7) = real(ngwf_basis%spheres(local_ngwf_idx)%offset,kind=DP)
    packed_ngwf(8:10) = &
         utils_logical_to_real(ngwf_basis%spheres(local_ngwf_idx)%extended(:))
    packed_ngwf(11) = &
         utils_logical_to_real(ngwfs_on_grid%iscmplx)
    packed_ngwf(12) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%pad1,kind=DP)
    packed_ngwf(13) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%pad2,kind=DP)
    packed_ngwf(14) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%pad3,kind=DP)
    packed_ngwf(15) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_ppds1,kind=DP)
    packed_ngwf(16) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_ppds2,kind=DP)
    packed_ngwf(17) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_ppds3,kind=DP)
    packed_ngwf(18) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_ppds1,kind=DP)
    packed_ngwf(19) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_ppds2,kind=DP)
    packed_ngwf(20) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_ppds3,kind=DP)
    packed_ngwf(21) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_pts1,kind=DP)
    packed_ngwf(22) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_pts2,kind=DP)
    packed_ngwf(23) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_pts3,kind=DP)
    packed_ngwf(24) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_pts1,kind=DP)
    packed_ngwf(25) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_pts2,kind=DP)
    packed_ngwf(26) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_pts3,kind=DP)
    packed_ngwf(27) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%tight_pts1,kind=DP)
    packed_ngwf(28) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%tight_pts2,kind=DP)
    packed_ngwf(29) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%tight_pts3,kind=DP)
    ! [*] end

    dest = packed_ngwf_header_length + 1
    do j = 1, 2
       do i = 1, n_ppds
          packed_ngwf(dest) = &
               real(ngwf_basis%spheres(local_ngwf_idx)%ppd_list(j,i),kind=DP)
          dest = dest+1
       end do
    end do

    n_pts_tot = n_ppds*ngwf_basis%n_pts

    if(ngwfs_on_grid%iscmplx) then
       last_index = dest+2*n_pts_tot-1
    else
       last_index = dest+n_pts_tot-1
    end if

    call utils_assert(last_index <= packed_ngwf_size, trim(myself)//&
         ': Insufficient packed_ngwf_size. At least this much is needed: ', &
         last_index)

    src_offset = ngwf_basis%spheres(local_ngwf_idx)%offset

    if(ngwfs_on_grid%iscmplx) then
       packed_ngwf(dest:dest+n_pts_tot-1) = &
            real(ngwfs_on_grid%z(src_offset:src_offset+n_pts_tot-1))
       packed_ngwf(dest+n_pts_tot:last_index) = &
            aimag(ngwfs_on_grid%z(src_offset:src_offset+n_pts_tot-1))
       packed_ngwf(last_index+1:) = 0.0_DP ! clear the tail
    else
       packed_ngwf(dest:last_index) = &
            ngwfs_on_grid%d(src_offset:src_offset+n_pts_tot-1)
       packed_ngwf(last_index+1:) = 0.0_DP ! clear the tail
    end if

    call timer_clock(myself,2)

  end subroutine remote_pack_ngwf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_pack_data_real(packed_data, data_in_ppds, ngwf_basis, &
       local_ngwf_idx, override_value)
    !==========================================================================!
    ! Packs information about real-valued data in PPDs that are the PPDs of an !
    ! NGWF. The data taken from ngwf_basis include the SPHERE, corresponding   !
    ! parts of its FUNCTIONS, and tightbox information). Packed data is a      !
    ! 1D array of reals that is easier to transfer between procs and is        !
    ! conveniently real regardless of whether real or complex NGWFs are used.  !
    !                                                                          !
    ! This is very similar to remote_pack_ngwf(), except:                      !
    !  - the data is assumed to be real,                                       !
    !  - source data is extracted from offset 1, and not from an offset corres-!
    !    ponding to local_ngwf_idx -- thus it is usable with PPD buffers that  !
    !    span one NGWF.                                                        !
    !  - An optional override value can be used instead of the data in PPDs.   !
    !                                                                          !
    ! The resultant data can be safely unpacked using remote_unpack_ngwfs()    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_data (out):   Stores the resultant packed data.                 !
    !   data_in_pps (in):    Real representation. Data will be extracted from  !
    !                        an offset of 1.                                   !
    !   ngwf_basis (in):     The SPHERE is extracted from here, and so is the  !
    !                        number of points per PPD.                         !
    !   local_ngwf_idx (in): Proc-local index of the NGWF in the ngwf_basis.   !
    !   override_value (in, opt): If passed, the data in PPDs is ignored, and  !
    !                             this uniform value is used instead for all   !
    !                             points.                                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2016.                              !
    ! Extended with override_value by Jacek Dziedzic in January 2017.          !
    !==========================================================================!

    use function_basis, only: FUNC_BASIS
    use utils, only: utils_assert, utils_logical_to_real

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)   :: packed_data(packed_ngwf_size)
    real(kind=DP), intent(in)    :: data_in_ppds(:)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    integer, intent(in)          :: local_ngwf_idx
    real(kind=DP), intent(in), optional :: override_value

    ! jd: Local variables
    integer :: i, j, dest, last_index
    integer :: n_pts_tot
    integer :: n_ppds
    character(len=*), parameter :: myself = 'remote_pack_data_real'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    call utils_assert(packed_ngwf_size /= -1, trim(myself)//&
         ': remote_init() must be called first.')

    n_ppds = ngwf_basis%spheres(local_ngwf_idx)%n_ppds_sphere

    packed_data(1) = ngwf_basis%spheres(local_ngwf_idx)%centre%x  ! [*] begin
    packed_data(2) = ngwf_basis%spheres(local_ngwf_idx)%centre%y
    packed_data(3) = ngwf_basis%spheres(local_ngwf_idx)%centre%z
    packed_data(4) = ngwf_basis%spheres(local_ngwf_idx)%radius
    packed_data(5) = real(n_ppds,kind=DP)
    packed_data(6) = real(ngwf_basis%n_pts)
    packed_data(7) = real(ngwf_basis%spheres(local_ngwf_idx)%offset,kind=DP)
    packed_data(8:10) = &
         utils_logical_to_real(ngwf_basis%spheres(local_ngwf_idx)%extended(:))
    packed_data(11) = utils_logical_to_real(.false.) ! (meaning not complex)
    packed_data(12) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%pad1,kind=DP)
    packed_data(13) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%pad2,kind=DP)
    packed_data(14) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%pad3,kind=DP)
    packed_data(15) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_ppds1,kind=DP)
    packed_data(16) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_ppds2,kind=DP)
    packed_data(17) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_ppds3,kind=DP)
    packed_data(18) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_ppds1,kind=DP)
    packed_data(19) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_ppds2,kind=DP)
    packed_data(20) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_ppds3,kind=DP)
    packed_data(21) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_pts1,kind=DP)
    packed_data(22) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_pts2,kind=DP)
    packed_data(23) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%start_pts3,kind=DP)
    packed_data(24) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_pts1,kind=DP)
    packed_data(25) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_pts2,kind=DP)
    packed_data(26) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%finish_pts3,kind=DP)
    packed_data(27) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%tight_pts1,kind=DP)
    packed_data(28) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%tight_pts2,kind=DP)
    packed_data(29) = real(ngwf_basis%tight_boxes(local_ngwf_idx)%tight_pts3,kind=DP)
    ! [*] end
    dest = packed_ngwf_header_length + 1
    do j = 1, 2
       do i = 1, n_ppds
          packed_data(dest) = &
               real(ngwf_basis%spheres(local_ngwf_idx)%ppd_list(j,i),kind=DP)
          dest = dest+1
       end do
    end do

    n_pts_tot = n_ppds*ngwf_basis%n_pts

    last_index = dest+n_pts_tot-1

    call utils_assert(last_index <= packed_ngwf_size, trim(myself)//&
         ': Insufficient packed_ngwf_size. At least this &
         &much is needed: ', last_index)

    if(present(override_value)) then
       packed_data(dest:last_index) = override_value
    else
       packed_data(dest:last_index) = data_in_ppds(1:n_pts_tot)
    end if
    packed_data(last_index+1:) = 0.0_DP ! clear the tail

    call timer_clock(myself,2)

  end subroutine remote_pack_data_real

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_unpack_ngwf(packed_ngwf, ngwf_sphere, ngwfs_on_grid, &
      ngwf_tightbox, dest_offset)
    !==========================================================================!
    ! Unpacks a packed NGWF into an ngwf_sphere, ngwf_tightbox (optionally)    !
    ! and (optionally) ngwfs_on_grid.                                          !
    ! Data is written out to ngwfs_on_grid at an offset of 1 (since the remote !
    ! offset is meaningless at the receiving end), unless 'dest_offset' is     !
    ! specified, which can be used to override this. Whichever offset was used,!
    ! it is subsequently unpacked to the ngwf_sphere,                          !
    !                                                                          !
    ! If ngwfs_on_grid is omitted, it is not unpacked. This is useful in       !
    ! in hfx_common_ppds().                                                    !
    !                                                                          !
    ! Similarly, if ngwf_tightbox is omitted, it is not unpacked.              !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_ngwf (in): The packed NGWF to unpack.                           !
    !   ngwf_sphere (out): Some of the unpacked data goes here.                !
    !   ngwfs_on_grid (out, optional): The contents of the PPDs stored in      !
    !                                  packed_ngwf are unpacked here.          !
    !   ngwf_tightbox(out, optional): The tightbox is unpacked here.           !
    !   dest_offset (in, optional): Override for unpacking at an offset        !
    !                               that is different from 1.                  !
    !--------------------------------------------------------------------------!
    ! CAVEATS:                                                                 !
    !   - Someone else is responsible for calling sphere_deallocate() on the   !
    !     returned sphere. It must be allocated here, because we do not know   !
    !     the number of PPDs in advance.                                       !
    !   - Someone else is responsible for calling data_functions_dealloc()     !
    !     on the returned ngwfs_on_grid (if passed).                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2013.                              !
    ! Adapted to complex NGWFs by Jacek Dziedzic in September 2016.            !
    !==========================================================================!

    use basis, only: SPHERE, FUNCTION_TIGHT_BOX
    use datatypes, only: FUNCTIONS, data_functions_alloc
    use utils, only: utils_alloc_check, utils_real_to_logical

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)              :: packed_ngwf(packed_ngwf_size)
    type(SPHERE), intent(out)              :: ngwf_sphere
    type(FUNCTIONS), intent(out), optional :: ngwfs_on_grid
    type(FUNCTION_TIGHT_BOX), intent(out), optional :: ngwf_tightbox
    integer, intent(in), optional          :: dest_offset

    ! jd: Local variables
    integer :: i, j, src, ierr
    integer :: n_pts, n_pts_tot
    integer :: loc_offset
    logical :: iscmplx
    character(len=*), parameter :: myself = 'remote_unpack_ngwf'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    ngwf_sphere%centre%x = packed_ngwf(1)                      ! [*] begin
    ngwf_sphere%centre%y = packed_ngwf(2)
    ngwf_sphere%centre%z = packed_ngwf(3)
    ngwf_sphere%radius = packed_ngwf(4)
    ngwf_sphere%n_ppds_sphere = nint(packed_ngwf(5))
    n_pts = nint(packed_ngwf(6))
    if(present(dest_offset)) then
       loc_offset = dest_offset
    else
       loc_offset = 1
    end if
    ngwf_sphere%offset = loc_offset
    ngwf_sphere%extended(:) = utils_real_to_logical(packed_ngwf(8:10))
    iscmplx = utils_real_to_logical(packed_ngwf(11))
    if(present(ngwf_tightbox)) then
       ngwf_tightbox%pad1 = nint(packed_ngwf(12))
       ngwf_tightbox%pad2 = nint(packed_ngwf(13))
       ngwf_tightbox%pad3 = nint(packed_ngwf(14))
       ngwf_tightbox%start_ppds1 = nint(packed_ngwf(15))
       ngwf_tightbox%start_ppds2 = nint(packed_ngwf(16))
       ngwf_tightbox%start_ppds3 = nint(packed_ngwf(17))
       ngwf_tightbox%finish_ppds1 = nint(packed_ngwf(18))
       ngwf_tightbox%finish_ppds2 = nint(packed_ngwf(19))
       ngwf_tightbox%finish_ppds3 = nint(packed_ngwf(20))
       ngwf_tightbox%start_pts1 = nint(packed_ngwf(21))
       ngwf_tightbox%start_pts2 = nint(packed_ngwf(22))
       ngwf_tightbox%start_pts3 = nint(packed_ngwf(23))
       ngwf_tightbox%finish_pts1 = nint(packed_ngwf(24))
       ngwf_tightbox%finish_pts2 = nint(packed_ngwf(25))
       ngwf_tightbox%finish_pts3 = nint(packed_ngwf(26))
       ngwf_tightbox%tight_pts1 = nint(packed_ngwf(27))
       ngwf_tightbox%tight_pts2 = nint(packed_ngwf(28))
       ngwf_tightbox%tight_pts3 = nint(packed_ngwf(29))         ! [*] end
    end if

    src = packed_ngwf_header_length + 1

    allocate(ngwf_sphere%ppd_list(2,ngwf_sphere%n_ppds_sphere),stat=ierr)
    call utils_alloc_check('remote_unpack_ngwf','current_sphere%ppd_list',ierr)
    ! jd: must match string in basis_sphere_deallocate ^^^

    do j = 1, 2
       do i = 1, ngwf_sphere%n_ppds_sphere
          ngwf_sphere%ppd_list(j,i) = nint(packed_ngwf(src))
          src = src + 1
       end do
    end do

    n_pts_tot = ngwf_sphere%n_ppds_sphere*n_pts

    if(present(ngwfs_on_grid)) then

       call data_functions_alloc(ngwfs_on_grid,n_pts_tot,iscmplx)
       !    ^ also sets ngwfs_on_grid%iscmplx

       if(iscmplx) then
          ngwfs_on_grid%z(loc_offset:loc_offset+n_pts_tot-1) = cmplx( &
               packed_ngwf(src:src+n_pts_tot-1), &
               packed_ngwf(src+n_pts_tot:src+2*n_pts_tot-1), kind=DP)
       else
          ngwfs_on_grid%d(loc_offset:loc_offset+n_pts_tot-1) = &
               packed_ngwf(src:src+n_pts_tot-1)
       end if
    end if

    call timer_clock(myself,2)

  end subroutine remote_unpack_ngwf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_unpack_ngwf_no_sphere(packed_ngwf, ppd_list, n_ppds, &
       ngwfs_on_grid_fn, ngwfs_on_grid_raw, dest_offset)
    !==========================================================================!
    ! Unpacks a packed NGWF into an array of ppd indices and ngwfs_on_grid.    !
    ! Data is written out to ngwfs_on_grid at an offset of 1 (since the remote !
    ! offset is meaningless at the receiving end), unless 'dest_offset' is     !
    ! specified, which can be used to override this. Whichever offset was used,!
    ! it is subsequently unpacked to the ngwf_sphere,                          !
    !                                                                          !
    ! Avoids the use of SPHERE type, because of difficulties with its          !
    ! allocatable component in OMP environments (firstprivate and allocatables !
    ! do not mix).                                                             !
    !                                                                          !
    ! ngwfs_on_grid can be passed by:                                          !
    !   - a FUNCTIONS object (ngwfs_on_grid_fn), or                            !
    !   - raw array (ngwfs_on_grid_raw).                                       !
    ! Both are optional. The latter only supports real NGWFs.                  !
    !                                                                          !
    ! If ngwfs_on_grid is omitted, only the PPD indices are unpacked.          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_ngwf (in): The packed NGWF to unpack.                           !
    !   ppd_list (out): The indices of the PPDs are written out here.          !
    !   n_ppds (out): The number of PPDs written out goes here.                !
    !   ngwfs_on_grid_fn (out, opt): The contents of the PPDs stored in        !
    !                                packed_ngwf are unpacked here, unless     !
    !                                omitted. Allocated here.                  !
    !   ngwfs_on_grid_raw (out, opt): The contents of the PPDs stored in       !
    !                                 packed_ngwf are unpacked here, unless    !
    !                                 omitted. Storage is caller's             !
    !                                 responsibility.                          !
    !   dest_offset (in, opt): Optional offset for unpacking to ngwfs_on_grid. !
    !--------------------------------------------------------------------------!
    ! CAVEATS:                                                                 !
    !   - Caller needs to allocate ppd_list big enough.                        !
    !   - Someone else is responsible for calling data_functions_dealloc()     !
    !     on the returned ngwfs_on_grid_fn (if passed).                        !
    !   - Elements of ppd_list beyond n_ppds are set to garbage_int.           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2014.                              !
    ! Adapted to complex NGWFs by Jacek Dziedzic in September 2016.            !
    ! Extended with the 'raw' version by Jacek Dziedzic in July 2019.          !
    !==========================================================================!

    use constants, only: garbage_int
    use datatypes, only: FUNCTIONS, data_functions_alloc
    use utils, only: utils_real_to_logical

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)               :: packed_ngwf(packed_ngwf_size)
    integer, intent(out)                    :: ppd_list(:)
    integer, intent(out)                    :: n_ppds
    type(FUNCTIONS), intent(out), optional  :: ngwfs_on_grid_fn
    real(kind=DP), intent(out), optional    :: ngwfs_on_grid_raw(:)
    integer, intent(in), optional           :: dest_offset

    ! jd: Local variables
    integer :: src, last_index
    integer :: n_pts_tot, n_pts
    integer :: loc_offset
    logical :: iscmplx
    character(len=*), parameter :: myself = 'remote_unpack_ngwf_no_sphere'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    n_ppds = nint(packed_ngwf(5))                       ! [*] begin
    n_pts = nint(packed_ngwf(6))
    iscmplx = utils_real_to_logical(packed_ngwf(11))
    if(present(dest_offset)) then
       loc_offset = dest_offset
    else
       loc_offset = 1
    end if

    src = packed_ngwf_header_length + 1                ! [*] end

    ! jd: Unpack PPD indices
    ppd_list(1:n_ppds) = nint(packed_ngwf(src:src+n_ppds-1))
    ppd_list(n_ppds+1:) = garbage_int

    src = src + 2 * n_ppds ! jd: 2, because we must also skip packed ppd_list(2,:)

    n_pts_tot = n_ppds*n_pts
    last_index = src+n_pts_tot-1

    if(present(ngwfs_on_grid_fn)) then
       call data_functions_alloc(ngwfs_on_grid_fn,n_pts_tot,iscmplx)
       !    ^ also sets ngwfs_on_grid%iscmplx
       if(iscmplx) then
          ngwfs_on_grid_fn%z(loc_offset:loc_offset+n_pts_tot-1) = cmplx( &
               packed_ngwf(src:src+n_pts_tot-1), &
               packed_ngwf(src+n_pts_tot:src+2*n_pts_tot-1), kind=DP)
       else
          ngwfs_on_grid_fn%d(loc_offset:loc_offset+n_pts_tot-1) = &
               packed_ngwf(src:src+n_pts_tot-1)
       end if
    end if

    if(present(ngwfs_on_grid_raw)) then
       ngwfs_on_grid_raw(loc_offset:loc_offset+n_pts_tot-1) = &
            packed_ngwf(src:src+n_pts_tot-1)
    end if

    call timer_clock(myself,2)

  end subroutine remote_unpack_ngwf_no_sphere

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_obtain_unpacked_ngwfs(ngwfs_in_ppds, offsets_to_ngwfs, &
       n_atoms, atom_indices, n_ppds_on_atoms, ngwf_basis, ngwf_cache_handle, &
       cell_npts)
    !==========================================================================!
    ! Obtains NGWFs of a set of atoms from an NGWF hash table and unpacks them !
    ! to the raw PPD representation.                                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ngwfs_in_ppds (out): Results (unpacked NGWFs go here). Indexed with a  !
    !                        fused index: point-PPD-ord_A-a, where ord_A counts!
    !                        atoms in atom_indices.                            !
    !   offsets_to_ngwfs (out): Offsets to subsequent atoms A in ngwfs_in_ppds !
    !                           are returned here for caller's convenience.    !
    !   n_atoms (in): Number of atoms in atoms whose NGWFs we are unpacking.   !
    !   atom_indices (in): Global atom indices of the atoms in question.       !
    !   n_ppds_on_atoms (in): Number of PPDs for each of the atoms in question.!
    !   ngwf_basis (in): The basis in which the NGWFs reside. Needed for       !
    !                    counting NGWFs.                                       !
    !   ngwf_cache_handle (in): Cache handle to the NGWF set -- for extracting !
    !                           NGWFs.                                         !
    !   cell_npts (in): Pass cell%n_pts. Needed for dimensioning arrays.       !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Both outputs (ngwfs_in_ppds, offsets_to_ngwfs) are allocated here.     !
    !   It is the responsibility of the caller to deallocate them.             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2020.                                 !
    !==========================================================================!

    use constants, only: garbage_real
    use function_basis, only: FUNC_BASIS
    use hash_table, only: hash_table_lookup_ptr_nocount
    use utils, only: utils_assert, utils_alloc_check

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out), allocatable :: ngwfs_in_ppds(:) ! pt-PPD-ord_A-a
    integer, intent(out), allocatable       :: offsets_to_ngwfs(:) ! ord_A
    integer, intent(in)                     :: n_atoms
    integer, intent(in)                     :: atom_indices(:)
    integer, intent(in)                     :: n_ppds_on_atoms(:)
    type(FUNC_BASIS), intent(in)            :: ngwf_basis
    integer, intent(in)                     :: ngwf_cache_handle
    integer, intent(in)                     :: cell_npts

    ! jd: Local variables
    real(kind=DP), pointer :: packed_ngwf_aa(:) ! ptr to ht's internals
    integer :: ppd_indices_on_a(ngwf_basis%max_n_ppds_sphere)
    integer :: n_ppds_on_a
    integer :: ord_a
    integer :: global_a
    integer :: n_ngwfs_a
    integer :: ngwf_a
    integer :: global_aa_ngwf_idx
    integer :: packed_ngwf_aa_size
    integer :: n_ppds_total
    integer :: n_points_total
    integer :: dest_offs
    integer :: ierr
    character(len=*), parameter :: myself = 'remote_obtain_unpacked_ngwfs'

    ! -------------------------------------------------------------------------

    allocate(offsets_to_ngwfs(max(n_atoms,1)),stat=ierr)
    call utils_alloc_check(myself,'offsets_to_ngwfs',ierr)

    ! jd: Ditto for points in Aa NGWF PPDs. Populate offsets_to_ngwfs.
    n_ppds_total = 0
    do ord_a = 1, n_atoms
       global_a = atom_indices(ord_a)
       n_ngwfs_a = ngwf_basis%num_on_atom(global_a)
       offsets_to_ngwfs(ord_a) = n_ppds_total * cell_npts + 1
       n_ppds_total = n_ppds_total + n_ppds_on_atoms(ord_a) * n_ngwfs_a
    end do
    n_points_total = n_ppds_total * cell_npts

    allocate(ngwfs_in_ppds(max(n_points_total,1)),stat=ierr)
    call utils_alloc_check(myself,'ngwfs_in_ppds',ierr)
    ngwfs_in_ppds(:) = garbage_real

    ! jd: Unpack Aa NGWFs
    do ord_a = 1, n_atoms
       global_a = atom_indices(ord_a)

       n_ngwfs_a = ngwf_basis%num_on_atom(global_a)

       ! -----------------------------------------------------------------------
       ! jd: for all a on A                                                  aaa
       ! -----------------------------------------------------------------------
       dest_offs = offsets_to_ngwfs(ord_a)
       loop_ngwf_a:                                                          &
       do ngwf_a = 1, n_ngwfs_a
          ! jd: Get NGWF Aa
          global_aa_ngwf_idx = ngwf_basis%first_on_atom(global_a) + ngwf_a - 1
          call hash_table_lookup_ptr_nocount(packed_ngwf_aa, &
               packed_ngwf_aa_size, remote_ngwf_cache_hts(ngwf_cache_handle), &
               global_aa_ngwf_idx)
          call utils_assert(packed_ngwf_aa_size /= -1, myself//&
               ': NGWF Aa absent in cache.', global_aa_ngwf_idx, global_a, ngwf_a)

          call remote_unpack_ngwf_no_sphere(packed_ngwf_aa, ppd_indices_on_a, &
               n_ppds_on_a, ngwfs_on_grid_raw = ngwfs_in_ppds(dest_offs:))

          call utils_assert(n_ppds_on_a == n_ppds_on_atoms(ord_a), myself//&
               ': Bookkeeping error.', n_ppds_on_a, n_ppds_on_atoms(ord_a), &
               ord_a, global_a)

          dest_offs = dest_offs + n_ppds_on_a * cell_npts

       end do loop_ngwf_a

    end do

  end subroutine remote_obtain_unpacked_ngwfs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_init_packed_ngwf_size(ngwf_basis, iscmplx)
    !==========================================================================!
    ! Initialises the module-wide 'constant' ngwf_packed_size to an appropriate!
    ! value (maximum for any NGWF in the system).                              !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ngwf_basis (in): The NGWF basis to use.                                !
    !   iscmplx (in): Whether to use complex NGWFs.                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2014.                               !
    !==========================================================================!

    use function_basis, only: FUNC_BASIS
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    logical, intent(in) :: iscmplx

    ! jd: Local variables
    integer :: desired_packed_ngwf_size
    character(len=*), parameter :: myself = 'remote_init_packed_ngwf_size'

    ! -------------------------------------------------------------------------

    desired_packed_ngwf_size = &
          ngwf_basis%max_n_ppds_sphere*ngwf_basis%n_pts + & ! PPD data
          ngwf_basis%max_n_ppds_sphere * 2 + &              ! ppd_list(:,:)
          packed_ngwf_header_length

    if(iscmplx) then
       desired_packed_ngwf_size = desired_packed_ngwf_size + &
            ngwf_basis%max_n_ppds_sphere*ngwf_basis%n_pts
    end if

    ! jd: Future-proof this somewhat
    call utils_assert(packed_ngwf_size == desired_packed_ngwf_size .or. &
         packed_ngwf_size == -1, &
         trim(myself)//': Attempt to use two NGWF basis sets with different &
         &max_n_ppds_sphere. If this is on purpose, generalize remote_mod.')

    packed_ngwf_size = desired_packed_ngwf_size

  end subroutine remote_init_packed_ngwf_size

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_init_persistent_caches
    !==========================================================================!
    ! Initialises persistent caches of the remote module.                      !
    ! Currently these are only the NGWF caches.                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2015.                              !
    !==========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: SIZEOF_DOUBLE
    use hash_table, only: hash_table_init
    use rundat, only: pub_use_remote_ngwfs
    use utils, only: utils_unit, utils_abort

    implicit none

    ! jd: Arguments

    ! jd: Local variables
    integer :: i_cache

    ! -------------------------------------------------------------------------

    ! jd: Open a file for debugging hash-tables, if desired
    remote_hash_table_info_unit = utils_unit()
    write(hash_table_info_filename,'(a,a,i0,a)') &
         trim(hash_table_info_rootname), &
         '_proc_', pub_my_proc_id, '.log'
    open(remote_hash_table_info_unit, file=hash_table_info_filename, err=10)

    ! jd: Init remote NGWF caches
    if(pub_use_remote_ngwfs) then
       do i_cache = 1, MAX_N_NGWF_CACHES
          call hash_table_init(remote_ngwf_cache_hts(i_cache), &
               'NGWFS', 1, maxslots_ngwf_cache, maxslots_ngwf_cache, &
               remote_hash_table_info_unit)
       end do
    end if

    return

10  call utils_abort('Could not create file: '//trim(hash_table_info_filename))

  end subroutine remote_init_persistent_caches

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_init(ngwf_basis, iscmplx)
    !==========================================================================!
    ! Initialises the module:                                                  !
    ! - ascertains packed_ngwf_size,                                           !
    ! - initialises persistent caches (NGWF cache),                            !
    ! - allocates persistent_packed_ngwf(:,:).                                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ngwf_basis (in): The NGWF basis to use. If using multiple NGWF basis   !
    !                    sets, use the set with the largets max_n_ppds_sphere, !
    !                    unless they all have the same max_n_ppds_sphere,      !
    !                    then any will do.                                     !
    !   iscmplx (in): Whether to use complex NGWFs.                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2014.                                 !
    !==========================================================================!

    use function_basis, only: FUNC_BASIS
    use utils, only: utils_alloc_check, utils_assert, &
         utils_flushed_string_output

    implicit none

    ! jd: Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    logical, intent(in) :: iscmplx

    ! jd: Local variables
    integer  :: ierr

    ! -------------------------------------------------------------------------

    call utils_flushed_string_output(CRLF//&
         'Initialising remote NGWF capability... ')

    call utils_assert(.not. initialised, &
         'remote_init: Attempt to initialise more than once')

    call remote_init_packed_ngwf_size(ngwf_basis, iscmplx)

    call remote_init_persistent_caches

    if(.not. allocated(persistent_packed_ngwf)) then
       allocate(persistent_packed_ngwf(packed_ngwf_size,MAX_N_NGWF_CACHES), &
            stat=ierr)
       call utils_alloc_check('remote_init','persistent_packed_ngwf',ierr)
    end if

    call utils_flushed_string_output('done. '//CRLF)

    initialised = .true.

  end subroutine remote_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_free_persistent_caches
    !==========================================================================!
    ! Frees persistent caches in remote_mod (currently NGWF cache) and closes  !
    ! the hash table log file.                                                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2015.                              !
    ! Modified by Jacek Dziedzic in July 2018 to free the HTs.                 !
    !==========================================================================!

    use hash_table, only: hash_table_free
    use rundat, only: pub_use_remote_ngwfs
    use utils, only: utils_abort

    implicit none

    ! jd: Local variables
    integer :: i_cache

    ! -------------------------------------------------------------------------

    close(remote_hash_table_info_unit, err=10)

    ! jd: Free remote NGWF caches. This actually frees the structures, if you're
    !     looking for the code that only purges the *contents* (invalidates
    !     the cache), it's in ngwf_representation.
    if(pub_use_remote_ngwfs) then
       do i_cache = 1, MAX_N_NGWF_CACHES
          call hash_table_free(remote_ngwf_cache_hts(i_cache))
       end do
    end if

    return

10  call utils_abort('Could not close file: '//hash_table_info_filename)

  end subroutine remote_free_persistent_caches

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_exit
    !==========================================================================!
    ! Cleans up after the module                                               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in April 2014.                                 !
    !==========================================================================!

    use utils, only: utils_dealloc_check, utils_assert, &
         utils_flushed_string_output

    implicit none

    integer :: ierr

    ! -------------------------------------------------------------------------

    call utils_flushed_string_output(CRLF//&
         'Cleaning up after remote NGWF capability... ')

    call utils_assert(initialised,'remote_exit: not initialised yet')

    call remote_free_persistent_caches

    if(allocated(persistent_packed_ngwf)) then
       deallocate(persistent_packed_ngwf,stat=ierr)
       call utils_dealloc_check('remote_exit','persistent_packed_ngwf',ierr)
       packed_ngwf_size = -1
    end if

    call utils_flushed_string_output(&
         'done. '//CRLF)

    initialised = .false.

  end subroutine remote_exit

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_obtain_ngwf(packed_ngwf, packed_ngwfs_hts, who_is_done, &
       global_ngwf, ngwfs_on_grid, ngwf_basis, cache_handle, &
       user_tag, ngwfs_on_grid2, ngwf_basis2, cache_handle2, &
       overfill_strategy)
    !==========================================================================!
    ! Obtains an NGWF potentially from a remote proc. Keeps serving similar    !
    ! requests from other procs while waiting for the NGWF to arrive.          !
    ! The NGWF is returned and served in a packed form, which stores all that  !
    ! is needed to reconstruct it at the receiving end (centre, radius, list   !
    ! of PPDs, contents of PPDs, tightbox, other SPHERE properties).           !
    ! The data to send is taken from nwgf_basis%spheres and ngwfs_on_grid.     !
    ! Data is returned in 'packed_ngwf'.                                       !
    !                                                                          !
    ! Note that 'offset' is transmitted too, but it is meaningless at the      !
    ! receiving end.                                                           !
    !                                                                          !
    ! Caching occurs behind the scenes through packed_ngwfs_hts hash table.    !
    !                                                                          !
    ! The result is returned by pointer -- to a hash table proc data, if the   !
    ! NGWF is found in the cache; or to a persistent packed NGWF, populated    !
    ! from ngwfs_on_grid, in the case of a cache miss.                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_ngwf (out): The returned NGWF, as 1D array with a size of       !
    !                      packed_ngwf_size.                                   !
    !   packed_ngwfs_hts (inout): Array of NGWF caches (one element per NGWF   !
    !                             set).                                        !
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   global_ngwf (in): Index of NGWF to obtain.                             !
    !   ngwfs_on_grid (in): NGWFs to serve are extracted from here.            !
    !   ngwf_basis (in): Describes the NGWF basis.                             !
    !   cache_handle (in): Identifies the NGWF set (1..MAX_N_NGWF_CACHES).     !
    !                                                                          !
    ! Arguments needed when exchanging multiple NGWF sets (conduction, etc):   !
    ! All OPTIONAL.                                                            !
    !                                                                          !
    !   user_tag (in): Used to select the set of NGWFs from which NGWFs are to !
    !                  be obtained. Defaults to 1 if omitted.                  !
    !   ngwfs_on_grid2, ngwf_basis2, cache_handle2 -- equivalents of above args!
    !                                                 corresponding to set 2.  !
    !                                                                          !
    ! Minor optional argument controlling how to handle overfull HTs.          !
    !   overfill_strategy (in, opt): cf. hash_table_add()                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2014.                               !
    ! DKN serving added by Jacek Dziedzic in March 2015.                       !
    ! Mixed NGWF set support added by Jacek Dziedzic in July, August 2018.     !
    ! Overfill strategy support added by Jacek Dziedzic in March 2019.         !
    ! Stripped of concurrent DKN exchange complexity by Jacek Dziedzic 2019.05.!
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_irecv, comms_test, comms_wait
    use datatypes, only: FUNCTIONS
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_probe, hash_table_add, &
         hash_table_lookup_ptr
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), pointer, intent(out)        :: packed_ngwf(:)
    ! ^note: ifort 18.0.[01], due to issue #03143573, required a
    ! @workaround: declaration of intent(out) pointers as intent(inout)
    type(HT_HASH_TABLE), intent(inout), target :: packed_ngwfs_hts(:)
    integer                                    :: global_ngwf
    type(FUNCTIONS), intent(in), target        :: ngwfs_on_grid
    type(FUNC_BASIS), intent(in), target       :: ngwf_basis
    integer, intent(in)                        :: cache_handle
    logical, intent(inout)           :: who_is_done(0:pub_total_num_procs-1)

    ! jd: Optional arguments for mixed NGWF sets
    integer, intent(in), optional                  :: user_tag
    type(FUNCTIONS), intent(in), optional, target  :: ngwfs_on_grid2
    type(FUNC_BASIS), intent(in), optional, target :: ngwf_basis2
    integer, intent(in), optional                  :: cache_handle2

    character, intent(in), optional                :: overfill_strategy

    ! jd: Local variables
    type(FUNC_BASIS), pointer    :: basis
    type(FUNCTIONS), pointer     :: on_grid
    integer                      :: cached_data_size
    integer                      :: local_ngwf
    integer                      :: owner_proc
    integer                      :: send_handle, recv_handle
    logical                      :: ready
    logical                      :: done
    integer                      :: loc_tag
    integer, save                :: calls_without_serving = 0
    character(len=*), parameter  :: myself = 'remote_obtain_ngwf'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    if(present(user_tag)) then
       loc_tag = user_tag
    else
       loc_tag = 1
    end if

    ! jd: If cached, retrieve from cache. This should be faster than getting
    !     it from ngwfs_on_grid as no packing has to be performed.
    !     Every once in a while the NGWFs are served too, so that no remote
    !     process ever gets starved when locally all NGWFs can be
    !     retrieved from cache
    if(-1 /= hash_table_probe(packed_ngwfs_hts(loc_tag), global_ngwf)) then
       ! jd: Cached, retrieve packed NGWF from cache
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_lookup_ptr',1)
#endif
       call hash_table_lookup_ptr(packed_ngwf, cached_data_size, &
            packed_ngwfs_hts(loc_tag), global_ngwf)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_lookup_ptr',2)
       call timer_clock('remote_obtain_ngwf_freq_serve',1)
#endif
       calls_without_serving = calls_without_serving + 1
       ! jd: Satisfy the requests of others every once in a while
       if(pub_total_num_procs > 1 .and. &
            calls_without_serving > ngwf_serve_period) then
          call remote_serve_ngwfs(done, who_is_done, &
               ngwfs_on_grid, ngwf_basis, cache_handle, &
               ngwfs_on_grid2, ngwf_basis2, cache_handle2)
          calls_without_serving = 0
       end if
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_freq_serve',2)
#endif
       call timer_clock(myself,2)
       return
    end if

    ! ==========================================================================
    ! jd: Careful, the below never executes any more once we reach a situation
    !     where all NGWFs can be retrieved from the cache.
    ! ==========================================================================

    ! jd: Establish which of the two NGWF bases we want
    if(loc_tag == cache_handle) then
       basis => ngwf_basis
       on_grid => ngwfs_on_grid
    else
       if(present(cache_handle2)) then
          if(loc_tag == cache_handle2) then
             basis => ngwf_basis2
             on_grid => ngwfs_on_grid2
          else
             call utils_assert(.false., myself//': user_tag does not match &
                  &either of the cache handles.', &
                  loc_tag, cache_handle, cache_handle2)
          end if
       else
          call utils_assert(.false., myself//': user_tag does not match &
               &cache handle.', loc_tag, cache_handle)
       end if
    end if

    owner_proc = basis%proc_of_func(global_ngwf)

    ! jd: An NGWF I own, just pack it and return
    if(owner_proc == pub_my_proc_id) then
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_own_pack',1)
#endif
       local_ngwf = global_ngwf - basis%first_on_proc(pub_my_proc_id)+1

       call remote_pack_ngwf(persistent_packed_ngwf(:,loc_tag), on_grid, basis,&
            local_ngwf)
       packed_ngwf => persistent_packed_ngwf(:,loc_tag)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_own_pack',2)
       call timer_clock('remote_obtain_ngwf_own_cache_add',1)
#endif
       ! jd: Store in cache too, so that future accesses proceed without
       !     packing
       call hash_table_add(packed_ngwfs_hts(loc_tag), packed_ngwf, &
            packed_ngwf_size, global_ngwf, overfill_strategy = overfill_strategy)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_own_cache_add',2)
#endif
    else
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_wait',1)
#endif
       ! jd: Not cached or owned, request from remote proc
       call comms_send(owner_proc, global_ngwf, 1, &
            tag = NGWF_REQUEST_TAG + loc_tag, return_handle = send_handle, &
            add_to_stack = .false.)
       ! jd: Initiate receive
       call comms_irecv(owner_proc, persistent_packed_ngwf(:,loc_tag), &
            packed_ngwf_size, tag = NGWF_DATA_TAG + loc_tag, handle = recv_handle)
       packed_ngwf => persistent_packed_ngwf(:,loc_tag)

       ! jd: Wait for the receive to complete, satisfying the requests of
       !     others in the meantime.
       ready = .false.
       do while(.not. ready)

          call comms_test(ready,recv_handle)

          ! jd: Satisfy the requests of others
          call remote_serve_ngwfs(done, who_is_done, &
               ngwfs_on_grid, ngwf_basis, cache_handle, &
               ngwfs_on_grid2, ngwf_basis2, cache_handle2)

       end do
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_wait',2)
       call timer_clock('remote_obtain_ngwf_cache_add',1)
#endif
       ! jd: Store in cache
       call hash_table_add(packed_ngwfs_hts(loc_tag), packed_ngwf, &
            packed_ngwf_size, global_ngwf, overfill_strategy=overfill_strategy)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_cache_add',2)
#endif
       call comms_wait(send_handle)
    end if

    call timer_clock(myself,2)

  end subroutine remote_obtain_ngwf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_obtain_ngwf_noadd(packed_ngwf, packed_ngwfs_hts, &
       who_is_done, global_ngwf, ngwfs_on_grid, ngwf_basis, cache_handle, &
       user_tag, ngwfs_on_grid2, ngwf_basis2, cache_handle2)
    !==========================================================================!
    ! Obtains an NGWF potentially from a remote proc. Keeps serving similar    !
    ! requests from other procs while waiting for the NGWF to arrive.          !
    ! The NGWF is returned and served in a packed form, which stores all that  !
    ! is needed to reconstruct it at the receiving end (centre, radius, list   !
    ! of PPDs, contents of PPDs, tightbox, other SPHERE properties).           !
    ! The data to send is taken from nwgf_basis%spheres and ngwfs_on_grid.     !
    ! Data is returned in 'packed_ngwf'.                                       !
    !                                                                          !
    ! Note that 'offset' is transmitted too, but it is meaningless at the      !
    ! receiving end.                                                           !
    !                                                                          !
    ! Caching occurs behind the scenes through packed_ngwfs_hts hash table.    !
    !                                                                          !
    ! The result is returned by pointer -- to a hash table proc data, if the   !
    ! NGWF is found in the cache; or to a persistent packed NGWF, populated    !
    ! from ngwfs_on_grid, in the case of a cache miss.                         !
    !                                                                          !
    ! This version does not add NGWFs to the cache, only reads from it.        !
    ! To be used in context where the cache is read-only.                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_ngwf (out): The returned NGWF, as 1D array with a size of       !
    !                      packed_ngwf_size.                                   !
    !   packed_ngwfs_hts (in): Array of NGWF caches (one element per NGWF set).!
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   global_ngwf (in): Index of NGWF to obtain.                             !
    !   ngwfs_on_grid (in): NGWFs to serve are extracted from here.            !
    !   ngwf_basis (in): Describes the NGWF basis.                             !
    !   cache_handle (in): Identifies the NGWF set (1..MAX_N_NGWF_CACHES).     !
    !                                                                          !
    ! Arguments needed when exchanging multiple NGWF sets (conduction, etc):   !
    ! All OPTIONAL.                                                            !
    !                                                                          !
    !   user_tag (in): Used to select the set of NGWFs from which NGWFs are to !
    !                  be obtained. Defaults to 1 if omitted.                  !
    !   ngwfs_on_grid2, ngwf_basis2, cache_handle2 -- equivalents of above args!
    !                                                 corresponding to set 2.  !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2016.                               !
    ! Mixed NGWF set support added by Jacek Dziedzic in July, August 2018.     !
    ! Overfill strategy support added by Jacek Dziedzic in March 2019.         !
    ! Stripped of concurrent DKN exchange complexity by Jacek Dziedzic 2019.05.!
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_irecv, comms_test, comms_wait
    use datatypes, only: FUNCTIONS
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_probe_nocount, &
         hash_table_lookup_ptr_nocount
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), pointer, intent(out)        :: packed_ngwf(:)
    ! ^note: ifort 18.0.[01], due to issue #03143573, required a
    ! @workaround: declaration of intent(out) pointers as intent(inout)
    type(HT_HASH_TABLE), intent(in), target    :: packed_ngwfs_hts(:)
    integer                                    :: global_ngwf
    type(FUNCTIONS), intent(in), target        :: ngwfs_on_grid
    type(FUNC_BASIS), intent(in), target       :: ngwf_basis
    integer, intent(in)                        :: cache_handle
    logical, intent(inout)           :: who_is_done(0:pub_total_num_procs-1)
    integer, intent(in), optional    :: user_tag
    type(FUNCTIONS), intent(in), optional, target  :: ngwfs_on_grid2
    type(FUNC_BASIS), intent(in), optional, target :: ngwf_basis2
    integer, intent(in), optional                  :: cache_handle2

    ! jd: Local variables
    type(FUNC_BASIS), pointer    :: basis
    type(FUNCTIONS), pointer     :: on_grid
    integer                      :: cached_data_size
    integer                      :: local_ngwf
    integer                      :: owner_proc
    integer                      :: send_handle, recv_handle
    logical                      :: ready
    logical                      :: done
    integer                      :: loc_tag
    integer, save                :: calls_without_serving = 0
    character(len=*), parameter  :: myself = 'remote_obtain_ngwf_noadd'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    if(present(user_tag)) then
       loc_tag = user_tag
    else
       loc_tag = 1
    end if

    ! jd: If cached, retrieve from cache. This should be faster than getting
    !     it from ngwfs_on_grid as no packing has to be performed.
    !     Every once in a while the NGWFs are served too, so that no remote
    !     process ever gets starved when locally all NGWFs can be
    !     retrieved from cache
    if(-1 /= hash_table_probe_nocount(packed_ngwfs_hts(loc_tag), global_ngwf)) then
       ! jd: Cached, retrieve packed NGWF from cache
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_noadd_lookup_ptr',1)
#endif
       call hash_table_lookup_ptr_nocount(packed_ngwf, cached_data_size, &
            packed_ngwfs_hts(loc_tag), global_ngwf)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_noadd_lookup_ptr',2)
       call timer_clock('remote_obtain_ngwf_noadd_freq_serve',1)
#endif
       calls_without_serving = calls_without_serving + 1
       ! jd: Satisfy the requests of others every once in a while
       if(pub_total_num_procs > 1 .and. &
            calls_without_serving > ngwf_serve_period) then
          call remote_serve_ngwfs(done, who_is_done, &
               ngwfs_on_grid, ngwf_basis, cache_handle, &
               ngwfs_on_grid2, ngwf_basis2, cache_handle2)
          calls_without_serving = 0
       end if
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_noadd_freq_serve',2)
#endif
       call timer_clock(myself,2)
       return
    end if

    ! ==========================================================================
    ! jd: Careful, the below never executes any more once we reach a situation
    !     where all NGWFs can be retrieved from the cache.
    ! ==========================================================================

    ! jd: Establish which of the two NGWF bases we want
    if(loc_tag == cache_handle) then
       basis => ngwf_basis
       on_grid => ngwfs_on_grid
    else
       if(present(cache_handle2)) then
          if(loc_tag == cache_handle2) then
             basis => ngwf_basis2
             on_grid => ngwfs_on_grid2
          else
             call utils_assert(.false., myself//': user_tag does not match &
                  &either of the cache handles.', &
                  loc_tag, cache_handle, cache_handle2)
          end if
       else
          call utils_assert(.false., myself//': user_tag does not match &
               &cache handle.', loc_tag, cache_handle)
       end if
    end if

    owner_proc = basis%proc_of_func(global_ngwf)

    ! jd: An NGWF I own, just pack it and return
    if(owner_proc == pub_my_proc_id) then
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_noadd_own_pack',1)
#endif
       local_ngwf = global_ngwf - basis%first_on_proc(pub_my_proc_id)+1

       call remote_pack_ngwf(persistent_packed_ngwf(:,loc_tag), on_grid, basis,&
            local_ngwf)
       packed_ngwf => persistent_packed_ngwf(:,loc_tag)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_noadd_own_pack',2)
#endif
    else
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_noadd_wait',1)
#endif
       ! jd: Not cached or owned, request from remote proc
       call comms_send(owner_proc, global_ngwf, 1, &
            tag = NGWF_REQUEST_TAG + loc_tag, return_handle = send_handle, &
            add_to_stack = .false.)
       ! jd: Initiate receive
       call comms_irecv(owner_proc, persistent_packed_ngwf(:,loc_tag), &
            packed_ngwf_size, tag = NGWF_DATA_TAG + loc_tag, &
            handle = recv_handle)
       packed_ngwf => persistent_packed_ngwf(:,loc_tag)

       ! jd: Wait for the receive to complete, satisfying the requests of
       !     others in the meantime.
       ready = .false.
       do while(.not. ready)

          call comms_test(ready,recv_handle)

          ! jd: Satisfy the requests of others
          call remote_serve_ngwfs(done, who_is_done, &
               ngwfs_on_grid, ngwf_basis, cache_handle, &
               ngwfs_on_grid2, ngwf_basis2, cache_handle2)

       end do
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_noadd_wait',2)
#endif
       call comms_wait(send_handle)
    end if

    call timer_clock(myself,2)

  end subroutine remote_obtain_ngwf_noadd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_obtain_ngwf_by_copy(packed_ngwf, packed_ngwfs_hts, &
       who_is_done, global_ngwf, ngwfs_on_grid, ngwf_basis, cache_handle, &
       user_tag, ngwfs_on_grid2, ngwf_basis2, cache_handle2, overfill_strategy)
    !==========================================================================!
    ! Obtains an NGWF potentially from a remote proc. Keeps serving similar    !
    ! requests from other procs while waiting for the NGWF to arrive.          !
    ! The NGWF is returned and served in a packed form, which stores all that  !
    ! is needed to reconstruct it at the receiving end (centre, radius, list   !
    ! of PPDs, contents of PPDs, tightbox, other SPHERE properties).           !
    ! The data to send is taken from nwgf_basis%spheres and ngwfs_on_grid.     !
    ! Data is returned in 'packed_ngwf'.                                       !
    !                                                                          !
    ! Note that 'offset' is transmitted too, but it is meaningless at the      !
    ! receiving end.                                                           !
    !                                                                          !
    ! Caching occurs behind the scenes through packed_ngwfs_hts hash table.    !
    !                                                                          !
    ! The result is returned by copy.                                          !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_ngwf (out): The returned NGWF, as 1D array with a size of       !
    !                      packed_ngwf_size.                                   !
    !   packed_ngwfs_hts (inout): Array of NGWF caches (one element per NGWF   !
    !                             set).                                        !
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   global_ngwf (in): Index of NGWF to obtain.                             !
    !   ngwfs_on_grid (in): NGWFs to serve are extracted from here.            !
    !   ngwf_basis (in): Describes the NGWF basis.                             !
    !   cache_handle (in): Identifies the NGWF set (1..MAX_N_NGWF_CACHES).     !
    !                                                                          !
    ! Arguments needed when exchanging multiple NGWF sets (conduction, etc):   !
    ! All OPTIONAL.                                                            !
    !                                                                          !
    !   user_tag (in): Used to select the set of NGWFs from which NGWFs are to !
    !                  be obtained. Defaults to 1 if omitted.                  !
    !   ngwfs_on_grid2, ngwf_basis2, cache_handle2 -- equivalents of above args!
    !                                                 corresponding to set 2.  !
    !                                                                          !
    ! Minor optional argument controlling how to handle overfull HTs.          !
    !   overfill_strategy (in, opt): cf. hash_table_add()                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2014.                               !
    ! DKN serving added by Jacek Dziedzic in March 2015.                       !
    ! Mixed NGWF set support added by Jacek Dziedzic in July, August 2018.     !
    ! Overfill strategy support added by Jacek Dziedzic in March 2019.         !
    ! Stripped of concurrent DKN exchange complexity by Jacek Dziedzic 2019.05.!
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_irecv, comms_test, comms_wait
    use datatypes, only: FUNCTIONS
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_probe, hash_table_add, &
         hash_table_lookup
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)                 :: packed_ngwf(packed_ngwf_size)
    type(HT_HASH_TABLE), intent(inout), target :: packed_ngwfs_hts(:)
    logical, intent(inout) :: who_is_done(0:pub_total_num_procs-1)
    integer                                    :: global_ngwf
    type(FUNCTIONS), intent(in), target        :: ngwfs_on_grid
    type(FUNC_BASIS), intent(in), target       :: ngwf_basis
    integer, intent(in)                        :: cache_handle
    integer, intent(in), optional              :: user_tag
    type(FUNCTIONS), intent(in), optional, target  :: ngwfs_on_grid2
    type(FUNC_BASIS), intent(in), optional, target :: ngwf_basis2
    integer, intent(in), optional                  :: cache_handle2
    character, intent(in), optional                :: overfill_strategy

    ! jd: Local variables
    type(FUNC_BASIS), pointer    :: basis
    type(FUNCTIONS), pointer     :: on_grid
    integer                      :: cached_data_size
    integer                      :: local_ngwf
    integer                      :: owner_proc
    integer                      :: send_handle, recv_handle
    logical                      :: ready
    logical                      :: done
    integer                      :: loc_tag
    integer, save                :: calls_without_serving = 0
    character(len=*), parameter  :: myself = 'remote_obtain_ngwf_by_copy'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    if(present(user_tag)) then
       loc_tag = user_tag
    else
       loc_tag = 1
    end if

    packed_ngwf(:) = 0.0_DP

    ! jd: If cached, retrieve from cache. This should be faster than getting
    !     it from ngwfs_on_grid as no packing has to be performed.
    !     Every once in a while the NGWFs are served too, so that no remote
    !     process ever gets starved when locally all NGWFs can be
    !     retrieved from cache
    if(-1 /= hash_table_probe(packed_ngwfs_hts(loc_tag), global_ngwf)) then
       ! jd: Cached, retrieve packed NGWF from cache
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_lookup',1)
#endif
       call hash_table_lookup(packed_ngwf, cached_data_size, &
            packed_ngwfs_hts(loc_tag), global_ngwf)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_lookup',2)
       call timer_clock('remote_obtain_ngwf_by_copy_freq_serv',1)
#endif
       calls_without_serving = calls_without_serving + 1
       ! jd: Satisfy the requests of others every once in a while
       if(pub_total_num_procs > 1 .and. &
            calls_without_serving > ngwf_serve_period) then
          call remote_serve_ngwfs(done, who_is_done, &
               ngwfs_on_grid, ngwf_basis, cache_handle, &
               ngwfs_on_grid2, ngwf_basis2, cache_handle2)
          calls_without_serving = 0
       end if
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_freq_serv',2)
#endif
       call timer_clock(myself,2)
       return
    end if

    ! ==========================================================================
    ! jd: Careful, the below never executes any more once we reach a situation
    !     where all NGWFs can be retrieved from the cache.
    ! ==========================================================================

    ! jd: Establish which of the two NGWF bases we want
    if(loc_tag == cache_handle) then
       basis => ngwf_basis
       on_grid => ngwfs_on_grid
    else
       if(present(cache_handle2)) then
          if(loc_tag == cache_handle2) then
             basis => ngwf_basis2
             on_grid => ngwfs_on_grid2
          else
             call utils_assert(.false., myself//': user_tag does not match &
                  &either of the cache handles.', &
                  loc_tag, cache_handle, cache_handle2)
          end if
       else
          call utils_assert(.false., myself//': user_tag does not match &
               &cache handle.', loc_tag, cache_handle)
       end if
    end if

    owner_proc = basis%proc_of_func(global_ngwf)

    ! jd: An NGWF I own, just pack it and return
    if(owner_proc == pub_my_proc_id) then
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_own_pack',1)
#endif
       local_ngwf = global_ngwf - basis%first_on_proc(pub_my_proc_id)+1

       call remote_pack_ngwf(packed_ngwf, on_grid, basis, local_ngwf)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_own_pack',2)
       call timer_clock('remote_obtain_ngwf_by_copy_cache_add',1)
#endif
       ! jd: Store in cache too, so that future accesses proceed without
       !     packing
       call hash_table_add(packed_ngwfs_hts(loc_tag), packed_ngwf, &
            packed_ngwf_size, global_ngwf, overfill_strategy=overfill_strategy)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_cache_add',2)
#endif
    else
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_wait',1)
#endif
       ! jd: Not cached or owned, request from remote proc
       call comms_send(owner_proc, global_ngwf, 1, &
            tag = NGWF_REQUEST_TAG + loc_tag, return_handle = send_handle, &
            add_to_stack = .false.)
       ! jd: Initiate receive
       call comms_irecv(owner_proc, packed_ngwf, packed_ngwf_size, &
            tag = NGWF_DATA_TAG + loc_tag, handle = recv_handle)

       ! jd: Wait for the receive to complete, satisfying the requests of
       !     others in the meantime.
       ready = .false.
       do while(.not. ready)

          call comms_test(ready,recv_handle)

          ! jd: Satisfy the requests of others
          call remote_serve_ngwfs(done, who_is_done, &
               ngwfs_on_grid, ngwf_basis, cache_handle, &
               ngwfs_on_grid2, ngwf_basis2, cache_handle2)

       end do
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_wait',2)
       call timer_clock('remote_obtain_ngwf_by_copy_cache_add',1)
#endif
       ! jd: Store in cache
       call hash_table_add(packed_ngwfs_hts(loc_tag), packed_ngwf, &
            packed_ngwf_size, global_ngwf, overfill_strategy=overfill_strategy)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_cache_add',2)
#endif
       call comms_wait(send_handle)
    end if

    call timer_clock(myself,2)

  end subroutine remote_obtain_ngwf_by_copy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_obtain_ngwf_by_copy_noadd(packed_ngwf, packed_ngwfs_hts, &
       who_is_done, global_ngwf, ngwfs_on_grid, ngwf_basis, cache_handle, &
       user_tag, ngwfs_on_grid2, ngwf_basis2, cache_handle2)
    !==========================================================================!
    ! Obtains an NGWF potentially from a remote proc. Keeps serving similar    !
    ! requests from other procs while waiting for the NGWF to arrive.          !
    ! The NGWF is returned and served in a packed form, which stores all that  !
    ! is needed to reconstruct it at the receiving end (centre, radius, list   !
    ! of PPDs, contents of PPDs, tightbox, other SPHERE properties).           !
    ! The data to send is taken from nwgf_basis%spheres and ngwfs_on_grid.     !
    ! Data is returned in 'packed_ngwf'.                                       !
    !                                                                          !
    ! Note that 'offset' is transmitted too, but it is meaningless at the      !
    ! receiving end.                                                           !
    !                                                                          !
    ! Caching occurs behind the scenes through packed_ngwfs_hts hash table.    !
    !                                                                          !
    ! The result is returned by copy.                                          !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_ngwf (out): The returned NGWF, as 1D array with a size of       !
    !                      packed_ngwf_size.                                   !
    !   packed_ngwfs_hts (in): Array of NGWF caches (one element per NGWF set).!
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   global_ngwf (in): Index of NGWF to obtain.                             !
    !   ngwfs_on_grid (in): NGWFs to serve are extracted from here.            !
    !   ngwf_basis (in): Describes the NGWF basis.                             !
    !   cache_handle (in): Identifies the NGWF set (1..MAX_N_NGWF_CACHES).     !
    !                                                                          !
    ! Arguments needed when exchanging multiple NGWF sets (conduction, etc):   !
    ! All OPTIONAL.                                                            !
    !                                                                          !
    !   user_tag (in): Used to select the set of NGWFs from which NGWFs are to !
    !                  be obtained. Defaults to 1 if omitted.                  !
    !   ngwfs_on_grid2, ngwf_basis2, cache_handle2 -- equivalents of above args!
    !                                                 corresponding to set 2.  !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2016.                               !
    ! Mixed NGWF set support added by Jacek Dziedzic in July, August 2018.     !
    ! Overfill strategy support added by Jacek Dziedzic in March 2019.         !
    ! Stripped of concurrent DKN exchange complexity by Jacek Dziedzic 2019.05.!
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_irecv, comms_test, comms_wait
    use datatypes, only: FUNCTIONS
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_probe_nocount, &
         hash_table_lookup_nocount
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)                 :: packed_ngwf(packed_ngwf_size)
    type(HT_HASH_TABLE), intent(in), target    :: packed_ngwfs_hts(:)
    logical, intent(inout) :: who_is_done(0:pub_total_num_procs-1)
    integer                                    :: global_ngwf
    type(FUNCTIONS), intent(in), target        :: ngwfs_on_grid
    type(FUNC_BASIS), intent(in), target       :: ngwf_basis
    integer, intent(in)                        :: cache_handle
    integer, intent(in), optional    :: user_tag
    type(FUNCTIONS), intent(in), optional, target  :: ngwfs_on_grid2
    type(FUNC_BASIS), intent(in), optional, target :: ngwf_basis2
    integer, intent(in), optional                  :: cache_handle2

    ! jd: Local variables
    type(FUNC_BASIS), pointer    :: basis
    type(FUNCTIONS), pointer     :: on_grid
    integer                      :: cached_data_size
    integer                      :: local_ngwf
    integer                      :: owner_proc
    integer                      :: send_handle, recv_handle
    logical                      :: ready
    logical                      :: done
    integer                      :: loc_tag
    integer, save                :: calls_without_serving = 0
    character(len=*), parameter  :: myself = 'remote_obtain_ngwf_by_copy_noadd'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    if(present(user_tag)) then
       loc_tag = user_tag
    else
       loc_tag = 1
    end if

    packed_ngwf(:) = 0.0_DP

    ! jd: If cached, retrieve from cache. This should be faster than getting
    !     it from ngwfs_on_grid as no packing has to be performed.
    !     Every once in a while the NGWFs are served too, so that no remote
    !     process ever gets starved when locally all NGWFs can be
    !     retrieved from cache
    if(-1 /= hash_table_probe_nocount(packed_ngwfs_hts(loc_tag), global_ngwf)) then
       ! jd: Cached, retrieve packed NGWF from cache
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_noadd_lookup',1)
#endif
       call hash_table_lookup_nocount(packed_ngwf, cached_data_size, &
            packed_ngwfs_hts(loc_tag), global_ngwf)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_noadd_lookup',2)
       call timer_clock('remote_obtain_ngwf_by_copy_noadd_freq_serv',1)
#endif
       calls_without_serving = calls_without_serving + 1
       ! jd: Satisfy the requests of others every once in a while
       if(pub_total_num_procs > 1 .and. &
            calls_without_serving > ngwf_serve_period) then
          call remote_serve_ngwfs(done, who_is_done, &
               ngwfs_on_grid, ngwf_basis, cache_handle, &
               ngwfs_on_grid2, ngwf_basis2, cache_handle2)
          calls_without_serving = 0
       end if
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_noadd_freq_serv',2)
#endif
       call timer_clock(myself,2)
       return
    end if

    ! ==========================================================================
    ! jd: Careful, the below never executes any more once we reach a situation
    !     where all NGWFs can be retrieved from the cache.
    ! ==========================================================================

    ! jd: Establish which of the two NGWF bases we want
    if(loc_tag == cache_handle) then
       basis => ngwf_basis
       on_grid => ngwfs_on_grid
    else
       if(present(cache_handle2)) then
          if(loc_tag == cache_handle2) then
             basis => ngwf_basis2
             on_grid => ngwfs_on_grid2
          else
             call utils_assert(.false., myself//': user_tag does not match &
                  &either of the cache handles.', &
                  loc_tag, cache_handle, cache_handle2)
          end if
       else
          call utils_assert(.false., myself//': user_tag does not match &
               &cache handle.', loc_tag, cache_handle)
       end if
    end if

    owner_proc = basis%proc_of_func(global_ngwf)

    ! jd: An NGWF I own, just pack it and return
    if(owner_proc == pub_my_proc_id) then
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_noadd_own_pack',1)
#endif
       local_ngwf = global_ngwf - basis%first_on_proc(pub_my_proc_id)+1

       call remote_pack_ngwf(packed_ngwf, on_grid, basis, local_ngwf)
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_noadd_own_pack',2)
#endif
    else
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_noadd_wait',1)
#endif
       ! jd: Not cached or owned, request from remote proc
       call comms_send(owner_proc, global_ngwf, 1, &
            tag = NGWF_REQUEST_TAG + loc_tag, return_handle = send_handle, &
            add_to_stack = .false.)
       ! jd: Initiate receive
       call comms_irecv(owner_proc, packed_ngwf, packed_ngwf_size, &
            tag = NGWF_DATA_TAG + loc_tag, handle = recv_handle)

       ! jd: Wait for the receive to complete, satisfying the requests of
       !     others in the meantime.
       ready = .false.
       do while(.not. ready)

          call comms_test(ready,recv_handle)

          ! jd: Satisfy the requests of others
          call remote_serve_ngwfs(done, who_is_done, &
               ngwfs_on_grid, ngwf_basis, cache_handle, &
               ngwfs_on_grid2, ngwf_basis2, cache_handle2)

       end do
#ifdef REMOTE_DETAILED_TIMINGS
       call timer_clock('remote_obtain_ngwf_by_copy_noadd_wait',2)
#endif
       call comms_wait(send_handle)
    end if

    call timer_clock(myself,2)

  end subroutine remote_obtain_ngwf_by_copy_noadd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_serve_ngwfs(done, who_is_done, &
       ngwfs_on_grid, ngwf_basis, cache_handle, &
       ngwfs_on_grid2, ngwf_basis2, cache_handle2)
    !==========================================================================!
    ! Keeps probing for requests from remote procs asking for NGWFs on a grid  !
    ! and satisfies them by sending requested entries from ngwfs_on_grid.      !
    ! The NGWFs are sent in a packed form.                                     !
    !                                                                          !
    ! In the mixed NGWFs scenario (conduction, etc) *both* sets are served, or !
    ! else deadlocks would occur. Whenever exchanging more than one NGWF set,  !
    ! always pass *both* sets, and identify each with its cache_handle.        !
    ! *Both* sets will be served here.                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   done (out): Set to .true. if all procs sent a 'done' notification.     !
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   ngwfs_on_grid (in): NGWFs to serve are extracted from here.            !
    !   ngwf_basis (in): Describes the NGWF basis, spheres and tightbox info   !
    !                    is taken from here.                                   !
    !   cache_handle (in): Identifies the NGWF set (1..MAX_N_NGWF_CACHES).     !
    !                                                                          !
    ! Arguments needed when exchanging multiple NGWF sets (conduction, etc):   !
    ! All OPTIONAL.                                                            !
    !                                                                          !
    !   ngwfs_on_grid2, ngwf_basis2, cache_handle2 -- equivalents of above args!
    !                                                 corresponding to set 2.  !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2014.                               !
    ! DKN serving added by Jacek Dziedzic in March 2015.                       !
    ! Generalised to two basis sets by Jacek Dziedzic in July, August 2018.    !
    ! Stripped of concurrent DKN exchange by Jacek Dziedzic, May 2019.
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
        comms_recv, comms_probe, comms_wait
    use datatypes, only: FUNCTIONS
    use function_basis, only: FUNC_BASIS
    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    logical, intent(out)                   :: done
    type(FUNCTIONS), intent(in), target    :: ngwfs_on_grid
    logical, intent(inout)                 :: who_is_done(0:pub_total_num_procs-1)
    type(FUNC_BASIS), intent(in), target   :: ngwf_basis
    integer, intent(in)                    :: cache_handle
    type(FUNCTIONS), intent(in), optional, target  :: ngwfs_on_grid2
    type(FUNC_BASIS), intent(in), optional, target :: ngwf_basis2
    integer, intent(in), optional                  :: cache_handle2

    ! Local variables
    type(FUNC_BASIS), pointer    :: basis
    type(FUNCTIONS), pointer     :: on_grid
    logical :: who_wants_ngwfs(0:pub_total_num_procs-1, MAX_N_NGWF_CACHES)
    integer :: who_wants_which_ngwfs(0:pub_total_num_procs-1, MAX_N_NGWF_CACHES)
    integer :: ngwf_to_send, ngwf_to_send_loc
    integer :: requested_ngwf
    integer :: send_handle
    integer :: proc
    integer :: tag
    real(kind=DP), dimension(packed_ngwf_size) :: packed_ngwf_to_send
    character(len=*), parameter :: myself = 'remote_serve_ngwfs'

    ! -------------------------------------------------------------------------

    ! Do not time this subroutine. All you'll get will be timer overhead.
#ifndef MPI
    call utils_abort(myself//': Calling me is pointless in serial mode.')
#endif

    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    if(present(ngwfs_on_grid2) .neqv. present(ngwf_basis2)) then
       call utils_abort(myself//': Optional arguments for mixed NGWF bases &
            &must either all be present or all be absent [1].')
    end if
    if(present(ngwfs_on_grid2) .neqv. present(cache_handle2)) then
       call utils_abort(myself//': Optional arguments for mixed NGWF bases &
            &must either all be present or all be absent [2].')
    end if

    ! jd: Probe for requests from others
    who_wants_ngwfs(:,:) = .false.
    do proc = 0, pub_total_num_procs - 1
       do tag = 1, MAX_N_NGWF_CACHES
          call comms_probe(who_wants_ngwfs(proc,tag), proc, NGWF_REQUEST_TAG + tag)
       end do
    end do

    ! jd: Find out what NGWFs the requests refer to, by actually receiving
    !     the requests
    who_wants_which_ngwfs(:,:) = -1
    do proc = 0, pub_total_num_procs - 1
       do tag = 1, MAX_N_NGWF_CACHES
          if(who_wants_ngwfs(proc,tag)) then

             call comms_recv(proc, requested_ngwf, 1, NGWF_REQUEST_TAG + tag)
             who_wants_which_ngwfs(proc,tag) = requested_ngwf

             ! jd: Is it a 'done' notification?
             if(requested_ngwf == -1) then
                who_is_done(proc) = .true.
                who_wants_ngwfs(proc,tag) = .false.
             end if
          end if
       end do
    end do

    done = all(who_is_done)

    ! jd: Satisfy the requests from other procs
    do proc = 0, pub_total_num_procs - 1
       do tag = 1, MAX_N_NGWF_CACHES
          if(who_wants_ngwfs(proc,tag)) then
             ngwf_to_send = who_wants_which_ngwfs(proc,tag)

             ! jd: Establish which of the two NGWF bases we want
             if(tag == cache_handle) then
                basis => ngwf_basis
                on_grid => ngwfs_on_grid
             else
                if(present(cache_handle2)) then
                   if(tag == cache_handle2) then
                      basis => ngwf_basis2
                      on_grid => ngwfs_on_grid2
                   else
                      call utils_assert(.false., &
                           myself//': tag requested by remote proc does not &
                           &match either of the cache handles.', &
                           tag, cache_handle, cache_handle2)
                   end if
                else
                   call utils_assert(.false., myself//': tag requested by remote&
                        & proc does not match cache handle.', tag, cache_handle)
                end if
             end if

             ngwf_to_send_loc = &
                  ngwf_to_send - basis%first_on_proc(pub_my_proc_id) + 1

             call remote_pack_ngwf(packed_ngwf_to_send, on_grid, basis, &
                  ngwf_to_send_loc)

             call comms_send(proc, packed_ngwf_to_send, packed_ngwf_size, &
                  NGWF_DATA_TAG + tag, send_handle, add_to_stack = .false.)
             ! @optimizeme
             ! The above only initiates a send, we need to wait before we can
             ! overwrite the buffer with a new packed_ngwf_to_send. If this leads
             ! to convoy effects, consider having an array of packed_ngwfs here
             ! and testing in a loop at the end. This will be more efficient, but
             ! comes with a memory penalty.
             call comms_wait(send_handle)

          end if
       end do
    end do

  end subroutine remote_serve_ngwfs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_obtain_dkn_block(packed_dkn_block, &
       packed_dkn_blocks_ht, who_is_done, atomx, atomy, denskern, &
       max_on_atom_row, max_on_atom_col, user_tag, overfill_strategy)
    !==========================================================================!
    ! Obtains density kernel blocks potentially from a remote proc, the one    !
    ! that owns 'atomy'. Keeps serving similar requests from other procs while !
    ! waiting for the density kernel block to arrive.                          !
    ! The blocks are served from 'denskern' and returned in 'packed_dkn_block'.!
    ! The obtained blocks always have the size                                 !
    ! max_on_atom_row * max_on_atom_col * num_spins. Use reshape on recipient. !
    ! Caching occurs behind the scenes through packed_dkn_blocks_ht hash table.!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_dkn_block (out): The returned density kernel block, as 1D array !
    !                           with a suitable size.                          !
    !   packed_dkn_blocks_ht (in/out): The hash table storing received blocks, !
    !                                  and returning them when needed again.   !
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   atomx (in): Index (row) of DKN block to obtain.                        !
    !   atomy (in): Index (col) of DKN block to obtain.                        !
    !   denskern (in): The density kernel matrix from which blocks are served  !
    !                  to other procs while we wait for the block we requested.!
    !   max_on_atom_row (in): Since we will soon be supporting mixed bases     !
    !                         here, caller must supply max_on_atom from the    !
    !                         row basis here.                                  !
    !   max_on_atom_col (in): Since we will soon be supporting mixed bases     !
    !                         here, caller must supply max_on_atom from the    !
    !                         col basis here.                                  !
    !   user_tag (in, opt): In the *future* may be used to distinguish sets of !
    !                       DKN blocks if multiple sets are communicated       !
    !                       concurrently, similarly to multiple NGWF sets.     !
    !                       Currently this is not supported ('serve' would have!
    !                       to scan all sets, like for NGWFs). For now, just   !
    !                       pass consistent values between here and 'serve'.   !
    !   overfill_strategy (in, opt): cf. hash_table_add()                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2014.                               !
    ! Extended in preparation for non-square blocks (e.g. hybrid TDDFT) by     !
    ! Jacek Dziedzic in August 2018.                                           !
    ! Overfill strategy support added by Jacek Dziedzic in March 2019.         !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_irecv, comms_test, comms_wait
    use hash_table, only: HT_HASH_TABLE, hash_table_probe, hash_table_add, &
         hash_table_lookup_ptr
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_get_block, sparse_get_par
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), pointer, intent(out) :: packed_dkn_block(:)
    ! ^note: ifort 18.0.[01], due to issue #03143573, required a
    ! @workaround: declaration of intent(out) pointers as intent(inout)
    type(HT_HASH_TABLE), intent(inout), target :: packed_dkn_blocks_ht
    logical, intent(inout)      :: who_is_done(0:pub_total_num_procs-1)
    integer, intent(in)                        :: atomx, atomy
    type(SPAM3), intent(in)                    :: denskern(pub_num_spins)
    integer, intent(in), optional              :: user_tag
    integer, intent(in)                        :: max_on_atom_row
    integer, intent(in)                        :: max_on_atom_col
    character, intent(in), optional            :: overfill_strategy

    ! jd: Local variables
    real(kind=DP) :: dkn_blk(max_on_atom_row, max_on_atom_col, pub_num_spins)

    integer :: cached_data_size
    integer :: owner_proc
    integer :: send_handle1, send_handle2, recv_handle
    logical :: ready
    logical :: done
    integer :: is
    integer :: packed_dkn_block_size
    integer :: loc_tag
    integer, save :: calls_without_serving = 0
    character(len=*), parameter :: myself = 'remote_obtain_dkn_block'
    type(PARAL_INFO), pointer :: par

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    call sparse_get_par(par, denskern(1))

    if(present(user_tag)) then
       loc_tag = user_tag
    else
       loc_tag = 1
    end if

    ! jd: If cached, retrieve from cache, by pointer. This should be faster
    !     than sparse_get_block() even if I own the DKN block, as no copying
    !     is performed. Every once in a while the DKN blocks are served too,
    !     so that no remote process ever gets starved when locally all blocks
    !     can be retrieved from cache
    if(-1 /= hash_table_probe(packed_dkn_blocks_ht, atomx, atomy)) then
       ! jd: Cached, retrieve block from cache, by pointer
       call hash_table_lookup_ptr(packed_dkn_block, cached_data_size, &
            packed_dkn_blocks_ht, atomx, atomy)
       calls_without_serving = calls_without_serving + 1
       ! jd: Satisfy the requests of others every once in a while
       if(pub_total_num_procs > 1 .and. &
            calls_without_serving > dkn_block_serve_period) then
          call remote_serve_dkn_blocks(done, who_is_done, denskern, &
               max_on_atom_row, max_on_atom_col, loc_tag)
          calls_without_serving = 0
       end if
       call timer_clock(myself,2)
       return
    end if

    ! ==========================================================================
    ! jd: Careful, the below never executes any more once we reach a situation
    !     where all DKN blocks can be retrieved from the cache.
    ! ==========================================================================

    owner_proc = par%proc_of_atom(par%orig_atom(atomy))

    packed_dkn_block_size = max_on_atom_row * max_on_atom_col * pub_num_spins

    ! jd: A dkn block I own, just return it, but store in cache too, so that
    !     future accesses happen by pointer, not by copy
    if(owner_proc == pub_my_proc_id) then
       dkn_blk(:,:,:) = 0.0_DP
       do is = 1, pub_num_spins
          call sparse_get_block(dkn_blk(:,:,is),denskern(is), atomx, atomy)
       end do

       call utils_assert(packed_dkn_block_size <= MAX_DKN_BLOCK_SIZE, &
            trim(myself)//': MAX_DKN_BLOCK_SIZE constant too small.')

       persistent_packed_dkn_blk(1:packed_dkn_block_size) = &
            reshape(dkn_blk(:,:,:),(/packed_dkn_block_size/))
       packed_dkn_block => persistent_packed_dkn_blk

       ! jd: Store in cache
       call hash_table_add(packed_dkn_blocks_ht, packed_dkn_block, &
            packed_dkn_block_size, atomx, atomy, &
            overfill_strategy = overfill_strategy)

    else
       ! jd: Not cached, not owned, request from remote proc
       call comms_send(owner_proc, atomx, 1, &
            tag = DKN_BLOCK_REQUEST_TAG + loc_tag, &
            return_handle = send_handle1, add_to_stack = .false.)
       call comms_send(owner_proc, atomy, 1, &
            tag = DKN_BLOCK_REQUEST_TAG + loc_tag, &
            return_handle = send_handle2, add_to_stack = .false.)

       ! jd: Initiate receive
       call comms_irecv(owner_proc, persistent_packed_dkn_blk, &
            packed_dkn_block_size, tag = DKN_BLOCK_DATA_TAG + loc_tag, &
            handle = recv_handle)
       packed_dkn_block => persistent_packed_dkn_blk

       ! jd: Wait for the receive to complete, satisfying the requests of
       !     others in the meantime.
       ready = .false.
       do while(.not. ready)

          call comms_test(ready,recv_handle)

          ! jd: Satisfy the requests of others
          call remote_serve_dkn_blocks(done, who_is_done, denskern, &
               max_on_atom_row, max_on_atom_col, loc_tag)

       end do

       ! jd: Store in cache
       call hash_table_add(packed_dkn_blocks_ht, packed_dkn_block, &
            packed_dkn_block_size, atomx, atomy, &
            overfill_strategy = overfill_strategy)

       call comms_wait(send_handle1)
       call comms_wait(send_handle2)

    end if

    call timer_clock(myself,2)

  end subroutine remote_obtain_dkn_block

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_serve_dkn_blocks(done, who_is_done, denskern, &
       max_on_atom_row, max_on_atom_col, user_tag)
    !==========================================================================!
    ! Keeps probing for requests from remote procs asking for density kernel   !
    ! blocks and satisfies them by sending requested entries from denskern.    !
    ! The blocks always have the size                                          !
    ! max_on_atom_row * max_on_atom_col * num_spins.                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   done (out): Set to .true. if all procs sent a 'done' notification.     !
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   denskern (in): The density kernel matrix from which blocks are served. !
    !   max_on_atom_row (in): Since we support mixed bases, caller must supply !
    !                         max_on_atom from the row basis here.             !
    !   max_on_atom_col (in): Since we support mixed bases, caller must supply !
    !                         max_on_atom from the col basis here.             !
    !   user_tag (in, opt): In the *future* may be used to distinguish sets of !
    !                       DKN blocks if multiple sets are communicated       !
    !                       concurrently, similarly to multiple NGWF sets.     !
    !                       Currently this is not supported ('serve' would have!
    !                       to scan all sets, like for NGWFs). For now, just   !
    !                       pass consistent values between here and 'serve'.   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2014.                               !
    ! Extended in preparation for non-square blocks (e.g. hybrid TDDFT) by     !
    ! Jacek Dziedzic in August 2018.                                           !
    !==========================================================================!

    use comms, only: pub_total_num_procs, comms_send, comms_recv, comms_probe, &
         comms_wait
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_get_block
    use utils, only: utils_alloc_check, utils_dealloc_check
#ifndef MPI
    use utils, only: utils_abort
#endif

    implicit none

    ! Arguments
    logical, intent(out)          :: done
    logical, intent(inout)        :: who_is_done(0:pub_total_num_procs-1)
    type(SPAM3), intent(in)       :: denskern(pub_num_spins)
    integer, intent(in)           :: max_on_atom_row
    integer, intent(in)           :: max_on_atom_col
    integer, intent(in), optional :: user_tag

    ! Local variables
    logical :: who_wants_dkn_blocks(0:pub_total_num_procs-1)
    integer :: who_wants_which_dkn_blocks_x(0:pub_total_num_procs-1)
    integer :: who_wants_which_dkn_blocks_y(0:pub_total_num_procs-1)
    integer :: dkn_block_to_send_x, dkn_block_to_send_y
    integer :: requested_dkn_block_x, requested_dkn_block_y
    integer :: send_handle
    integer :: proc
    integer :: is
    integer :: ierr
    integer :: packed_dkn_block_size
    integer :: loc_tag
    real(kind=DP) :: dkn_blk(max_on_atom_row, max_on_atom_col, pub_num_spins)
    real(kind=DP), allocatable :: packed_dkn_block_to_send(:)
    character(len=*), parameter :: myself = 'remote_serve_dkn_blocks'

    ! -------------------------------------------------------------------------

    ! Do not time this subroutine. All you'll get will be timer overhead.

    if(present(user_tag)) then
       loc_tag = user_tag
    else
       loc_tag = 1
    end if

#ifndef MPI
    call utils_abort('remote_serve_dkn_blocks: Calling me is pointless in &
         &serial mode')
#endif

    packed_dkn_block_size = max_on_atom_row * max_on_atom_col * pub_num_spins

    ! jd: Probe for requests from others
    who_wants_dkn_blocks(:) = .false.
    do proc = 0, pub_total_num_procs - 1
       call comms_probe(who_wants_dkn_blocks(proc), proc, &
            DKN_BLOCK_REQUEST_TAG + loc_tag)
    end do

    allocate(packed_dkn_block_to_send(packed_dkn_block_size),stat=ierr)
    call utils_alloc_check(myself, 'packed_dkn_block_to_send', ierr)

    ! jd: Find out what DKN blocks do the requests refer to, by actually
    !     receiving the requests
    who_wants_which_dkn_blocks_x(:) = -1
    who_wants_which_dkn_blocks_y(:) = -1
    do proc = 0, pub_total_num_procs - 1
       if(who_wants_dkn_blocks(proc)) then

          call comms_recv(proc, requested_dkn_block_x, 1, &
               DKN_BLOCK_REQUEST_TAG + loc_tag)
          call comms_recv(proc, requested_dkn_block_y, 1, &
               DKN_BLOCK_REQUEST_TAG + loc_tag)

          who_wants_which_dkn_blocks_x(proc) = requested_dkn_block_x
          who_wants_which_dkn_blocks_y(proc) = requested_dkn_block_y

          ! jd: Is it a 'done' notification?
          if(requested_dkn_block_x == -1) then
             who_is_done(proc) = .true.
             who_wants_dkn_blocks(proc) = .false.
          end if
       end if
    end do

    done = all(who_is_done)

    ! jd: Satisfy the requests from other procs
    do proc = 0, pub_total_num_procs - 1
       if(who_wants_dkn_blocks(proc)) then
          dkn_block_to_send_x = who_wants_which_dkn_blocks_x(proc)
          dkn_block_to_send_y = who_wants_which_dkn_blocks_y(proc)

          dkn_blk(:,:,:) = 0.0_DP

          do is = 1, pub_num_spins
             call sparse_get_block(dkn_blk(:,:,is),denskern(is), &
                  dkn_block_to_send_x, dkn_block_to_send_y)
          end do

          packed_dkn_block_to_send = &
               reshape(dkn_blk(:,:,:),(/packed_dkn_block_size/))
          call comms_send(proc, packed_dkn_block_to_send, &
               packed_dkn_block_size, DKN_BLOCK_DATA_TAG + loc_tag, &
               send_handle, add_to_stack = .false.)
          ! @optimizeme
          ! The above only initiates a send, we need to wait before we can
          ! overwrite the buffer with a new packed_dkn_block_to_send. If this
          ! leads to convoy effects, consider having an array of
          ! packed_dkn_blocks here and testing in a loop at the end. This will
          ! be more efficient, but comes with a memory penalty.
          call comms_wait(send_handle)

       end if
    end do

    deallocate(packed_dkn_block_to_send,stat=ierr)
    call utils_dealloc_check(myself, 'packed_dkn_block_to_send', ierr)

  end subroutine remote_serve_dkn_blocks

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_obtain_matrix_block(packed_matrix_block, &
       packed_matrix_blocks_ht, who_is_done, atomx, atomy, matrix, &
       max_on_atom_row, max_on_atom_col, overfill_strategy)
    !==========================================================================!
    ! Obtains SPAM3 matrix blocks potentially from a remote proc, the one      !
    ! that owns 'atomy'. Keeps serving similar requests from other procs while !
    ! waiting for the matrix block to arrive. The blocks are served from       !
    ! 'matrix' and returned by pointer in 'packed_matrix_block'. The obtained  !
    ! blocks always have the size max_on_atom_row * max_on_atom_col.           !
    ! Use reshape on the receiving end if necessary.                           !
    ! Caching occurs behind the scenes through the packed_matrix_blocks_ht     !
    ! hash table. If 'overfill_strategy' is 'F', you can rely on the blocks    !
    ! always being in the cache.                                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   packed_matrix_block (out): The returned matrix block, as pointer to a  !
    !                              1D array with a suitable size.              !
    !   packed_matrix_blocks_ht (in/out): The hash table storing received      !
    !                                     blocks and returning them when       !
    !                                     needed later.                        !
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   atomx (in): Index (row) of the matrix block to obtain.                 !
    !   atomy (in): Index (col) of the matrix block to obtain.                 !
    !   matrix (in): The SPAM3 matrix from which blocks are served to other    !
    !                procs while we wait for the blocks we requested.          !
    !   max_on_atom_row (in): Since we support mixed bases, caller must supply !
    !                         max_on_atom from the row basis here.             !
    !   max_on_atom_col (in): Since we support mixed bases, caller must supply !
    !                         max_on_atom from the col basis here.             !
    !   overfill_strategy (in, opt): cf. hash_table_add()                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2019 using remote_obtain_dkn_block() as !
    ! template.                                                                !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_irecv, comms_test, comms_wait
    use hash_table, only: HT_HASH_TABLE, hash_table_probe, hash_table_add, &
         hash_table_lookup_ptr
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_get_block, sparse_get_par
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), pointer, intent(out) :: packed_matrix_block(:)
    type(HT_HASH_TABLE), intent(inout), target :: packed_matrix_blocks_ht
    logical, intent(inout)      :: who_is_done(0:pub_total_num_procs-1)
    integer, intent(in)                        :: atomx, atomy
    type(SPAM3), intent(in)                    :: matrix
    integer, intent(in)                        :: max_on_atom_row
    integer, intent(in)                        :: max_on_atom_col
    character, intent(in), optional            :: overfill_strategy

    ! jd: Local variables
    real(kind=DP) :: matrix_blk(max_on_atom_row, max_on_atom_col)

    integer :: cached_data_size
    integer :: owner_proc
    integer :: send_handle1, send_handle2, recv_handle
    logical :: ready
    logical :: done
    integer :: packed_block_size
    integer, save :: calls_without_serving = 0
    character(len=*), parameter :: myself = 'remote_obtain_matrix_block'
    type(PARAL_INFO), pointer :: par

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    call sparse_get_par(par, matrix)

    ! jd: If cached, retrieve from cache, by pointer. This should be faster
    !     than sparse_get_block() even if I own the matrix block, as no copying
    !     is performed. Every once in a while the matrix blocks are served too,
    !     so that no remote process ever gets starved when locally all blocks
    !     can be retrieved from cache
    if(-1 /= hash_table_probe(packed_matrix_blocks_ht, atomx, atomy)) then
       ! jd: Cached, retrieve block from cache, by pointer
       call hash_table_lookup_ptr(packed_matrix_block, cached_data_size, &
            packed_matrix_blocks_ht, atomx, atomy)
       calls_without_serving = calls_without_serving + 1
       ! jd: Satisfy the requests of others every once in a while
       if(pub_total_num_procs > 1 .and. &
            calls_without_serving > matrix_block_serve_period) then
          call remote_serve_matrix_blocks(done, who_is_done, matrix, &
               max_on_atom_row, max_on_atom_col)
          calls_without_serving = 0
       end if
       call timer_clock(myself,2)
       return
    end if

    ! ==========================================================================
    ! jd: Careful, the below never executes any more once we reach a situation
    !     where all matrix blocks can be retrieved from the cache.
    ! ==========================================================================

    owner_proc = par%proc_of_atom(par%orig_atom(atomy))

    packed_block_size = max_on_atom_row * max_on_atom_col

    ! jd: A matrix block I own, just return it, but store in cache too, so that
    !     future accesses happen by pointer, not by copy
    if(owner_proc == pub_my_proc_id) then
       matrix_blk(:,:) = 0.0_DP
       call sparse_get_block(matrix_blk(:,:), matrix, atomx, atomy)

       call utils_assert(packed_block_size <= MAX_MATRIX_BLOCK_SIZE, &
            trim(myself)//': MAX_MATRIX_BLOCK_SIZE constant too small.')

       persistent_packed_matrix_blk(1:packed_block_size) = &
            reshape(matrix_blk(:,:),(/packed_block_size/))
       packed_matrix_block => persistent_packed_matrix_blk

       ! jd: Store in cache
       call hash_table_add(packed_matrix_blocks_ht, packed_matrix_block, &
            packed_block_size, atomx, atomy, &
            overfill_strategy = overfill_strategy)

    else
       ! jd: Not cached, not owned, request from remote proc
       call comms_send(owner_proc, atomx, 1, &
            tag = MATRIX_BLOCK_REQUEST_TAG, &
            return_handle = send_handle1, add_to_stack = .false.)
       call comms_send(owner_proc, atomy, 1, &
            tag = MATRIX_BLOCK_REQUEST_TAG, &
            return_handle = send_handle2, add_to_stack = .false.)

       ! jd: Initiate receive
       call comms_irecv(owner_proc, persistent_packed_matrix_blk, &
            packed_block_size, tag = MATRIX_BLOCK_DATA_TAG, &
            handle = recv_handle)
       packed_matrix_block => persistent_packed_matrix_blk

       ! jd: Wait for the receive to complete, satisfying the requests of
       !     others in the meantime.
       ready = .false.
       do while(.not. ready)

          call comms_test(ready,recv_handle)

          ! jd: Satisfy the requests of others
          call remote_serve_matrix_blocks(done, who_is_done, matrix, &
               max_on_atom_row, max_on_atom_col)

       end do

       ! jd: Store in cache
       call hash_table_add(packed_matrix_blocks_ht, packed_matrix_block, &
            packed_block_size, atomx, atomy, &
            overfill_strategy = overfill_strategy)

       call comms_wait(send_handle1)
       call comms_wait(send_handle2)

    end if

    call timer_clock(myself,2)

  end subroutine remote_obtain_matrix_block

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_serve_matrix_blocks(done, who_is_done, matrix, &
       max_on_atom_row, max_on_atom_col)
    !==========================================================================!
    ! Keeps probing for requests from remote procs asking for SPAM3 matrix     !
    ! blocks and satisfies them by sending requested entries from 'matrix'.    !
    ! The blocks always have the size max_on_atom_row * max_on_atom_col.       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   done (out): Set to .true. if all procs sent a 'done' notification.     !
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   matrix (in): The SPAM3 matrix from which blocks are served.            !
    !   max_on_atom_row (in): Since we support mixed bases, caller must supply !
    !                         max_on_atom from the row basis here.             !
    !   max_on_atom_col (in): Since we support mixed bases, caller must supply !
    !                         max_on_atom from the col basis here.             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2019 using remote_serve_dkn_blocks() as !
    ! template.                                                                !
    !==========================================================================!

    use comms, only: pub_total_num_procs, comms_send, comms_recv, comms_probe, &
         comms_wait
    use sparse, only: SPAM3, sparse_get_block
    use utils, only: utils_alloc_check, utils_dealloc_check
#ifndef MPI
    use utils, only: utils_abort
#endif

    implicit none

    ! Arguments
    logical, intent(out)          :: done
    logical, intent(inout)        :: who_is_done(0:pub_total_num_procs-1)
    type(SPAM3), intent(in)       :: matrix
    integer, intent(in)           :: max_on_atom_row
    integer, intent(in)           :: max_on_atom_col

    ! Local variables
    logical :: who_wants_matrix_blocks(0:pub_total_num_procs-1)
    integer :: who_wants_which_matrix_blocks_x(0:pub_total_num_procs-1)
    integer :: who_wants_which_matrix_blocks_y(0:pub_total_num_procs-1)
    integer :: matrix_block_to_send_x, matrix_block_to_send_y
    integer :: requested_matrix_block_x, requested_matrix_block_y
    integer :: send_handle
    integer :: proc
    integer :: ierr
    integer :: packed_block_size
    real(kind=DP) :: matrix_blk(max_on_atom_row, max_on_atom_col)
    real(kind=DP), allocatable :: packed_matrix_block_to_send(:)
    character(len=*), parameter :: myself = 'remote_serve_matrix_blocks'

    ! -------------------------------------------------------------------------

    ! Do not time this subroutine. All you'll get will be timer overhead.

#ifndef MPI
    call utils_abort(myself//': Calling me is pointless in &
         &serial mode')
#endif

    packed_block_size = max_on_atom_row * max_on_atom_col

    ! jd: Probe for requests from others
    who_wants_matrix_blocks(:) = .false.
    do proc = 0, pub_total_num_procs - 1
       call comms_probe(who_wants_matrix_blocks(proc), proc, &
            MATRIX_BLOCK_REQUEST_TAG)
    end do

    allocate(packed_matrix_block_to_send(packed_block_size),stat=ierr)
    call utils_alloc_check(myself, 'packed_matrix_block_to_send', ierr)

    ! jd: Find out what matrix blocks do the requests refer to, by actually
    !     receiving the requests
    who_wants_which_matrix_blocks_x(:) = -1
    who_wants_which_matrix_blocks_y(:) = -1
    do proc = 0, pub_total_num_procs - 1
       if(who_wants_matrix_blocks(proc)) then

          call comms_recv(proc, requested_matrix_block_x, 1, &
               MATRIX_BLOCK_REQUEST_TAG)
          call comms_recv(proc, requested_matrix_block_y, 1, &
               MATRIX_BLOCK_REQUEST_TAG)

          who_wants_which_matrix_blocks_x(proc) = requested_matrix_block_x
          who_wants_which_matrix_blocks_y(proc) = requested_matrix_block_y

          ! jd: Is it a 'done' notification?
          if(requested_matrix_block_x == -1) then
             who_is_done(proc) = .true.
             who_wants_matrix_blocks(proc) = .false.
          end if
       end if
    end do

    done = all(who_is_done)

    ! jd: Satisfy the requests from other procs
    do proc = 0, pub_total_num_procs - 1
       if(who_wants_matrix_blocks(proc)) then
          matrix_block_to_send_x = who_wants_which_matrix_blocks_x(proc)
          matrix_block_to_send_y = who_wants_which_matrix_blocks_y(proc)

          matrix_blk(:,:) = 0.0_DP

          call sparse_get_block(matrix_blk(:,:),matrix, &
               matrix_block_to_send_x, matrix_block_to_send_y)

          packed_matrix_block_to_send = &
               reshape(matrix_blk(:,:),(/packed_block_size/))
          call comms_send(proc, packed_matrix_block_to_send, &
               packed_block_size, MATRIX_BLOCK_DATA_TAG, &
               send_handle, add_to_stack = .false.)
          ! @optimizeme
          ! The above only initiates a send, we need to wait before we can
          ! overwrite the buffer with a new packed_matrix_block_to_send. If this
          ! leads to convoy effects, consider having an array of
          ! packed_matrix_blocks here and testing in a loop at the end. This will
          ! be more efficient, but comes with a memory penalty.
          call comms_wait(send_handle)

       end if
    end do

    deallocate(packed_matrix_block_to_send,stat=ierr)
    call utils_dealloc_check(myself, 'packed_matrix_block_to_send', ierr)

  end subroutine remote_serve_matrix_blocks

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_obtain_coeffs(coeffs_ht, who_is_done, &
       global_b, global_c, coeffs_kind, swex_h, ngwf_basis, elements, swri_h, par, &
       is, ngwf_basis2, ket_basis_selector, overfill_strategy)
    !==========================================================================!
    ! Obtains already calculated expansion coefficients for a pair of NGWFs,   !
    ! BbCc, potentially from a remote proc.                                    !
    ! For a description of coefficients, see header in remote_serve_coeffs().  !
    !                                                                          !
    ! Traditionally (until mixed NGWF basis support), we could employ symmetry !
    ! as the expansion for (Bb,Cc) was the same as that for (Cc,Bb). This is   !
    ! still true in the single NGWF basis scenario. In this scenario, Bb and Cc!
    ! *must* be sorted on entry, so that Bb <= Cc (this will be checked).      !
    !                                                                          !
    ! For mixed bases, this symmetry no longer holds, and we store and communi-!
    ! cate (Bb,Cc) where Bb is not necessarily <= Cc (so, all of them, modulo  !
    ! subsystem-SWRI scenarios). Whether we are in the symmetric or the        !
    ! non-symmetric scenario is determined from SWEX_BB_CC_SYMMETRIC(:) and    !
    ! the swex handle -- essentially it's hardcoded which swexes are symmetric !
    ! and which aren't (currently all swexes except REP_SWEX_HFX_OTHER are).   !
    ! Do *not* confuse this with cache handles!                                !
    !                                                                          !
    ! Who 'owns' the coefficients? Normally it's the proc owning the NGWF Bb,  !
    ! and this is always the case when expanding the entire system (provided   !
    ! symmetry, see above). However, when expanding a subsystem, we need to    !
    ! handle the scenario where B is not in the expansion ('not in the SWRI'), !
    ! while C is. If, additionally, C<B, the expansion coefficients will be    !
    ! added on, and reside on C, so the owner of NGWF Cc.                      !
    !                                                                          !
    ! This subroutine also keeps serving similar requests from other procs,    !
    ! but *only* while waiting for the coefficients it itself asked for to     !
    ! arrive. It does not serve every N requests, like is done with NGWFs or   !
    ! DKN, because we do not mix coeff exchange with computation.              !
    ! The coefficients are served from and received into coeffs_ht.            !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   coeffs_ht (in/out): The hash table containing the coefficients (both   !
    !                       the served and received).                          !
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   global_b (in): Global index of B atom.                                 !
    !   global_c (in): Global index of C atom.                                 !
    !   coeffs_kind (in): The type of the coefficients to serve (SW_V, SW_O,   !
    !                     SW_V+SW_D, SW_O+SW_D).
    !   swex_h (in): A handle to the SWEX associated with this expansion.      !
    !                This goes into one of the keys indexing coeffs and allows !
    !                storing coeffs in a single HT disassociated from a rep    !
    !                in the presence of Bessel averaging. It is also used to   !
    !                determine if we have Bb/Cc symmetry.                      !
    !   ngwf_basis (in): Required to find out owners of NGWFs.                 !
    !   elements (in): Required, unfortunately, to find out owners of coeffs.  !
    !   swri_h (in): Required to find out owners of coeffs.                    !
    !   is (in, opt): Pass spin index here if the coeffs are spin-dependent.   !
    !                 See header in remote_serve_coeffs() for more detail.     !
    !                                                                          !
    !   *Optional* arguments for mixed NGWF bases.                             !
    !   ngwf_basis2 (in): Second NGWF basis, for finding out owners of NGWFs.  !
    !   ket_basis_selector (in): Selects the NGWF basis (1 or 2) in which      !
    !                            - Bb is indexed (ket_basis_selector(1)),      !
    !                            - Cc is indexed (ket_basis_selector(2)).      !
    !                            ket_basis_selector(3) is ignored here.        !
    !   overfill_strategy (in, opt): cf. hash_table_add()                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2014.                               !
    ! Improved by Jacek Dziedzic in July 2016 to batch coeffs over all NGWF    !
    ! pairs on an atom pair.                                                   !
    ! Extended by Jacek Dziedzic in June 2017 to support exchange of spin-     !
    ! dependent and spin-agnostic coeffs.                                      !
    ! Mixed NGWF set support added by Jacek Dziedzic in July, August 2018.     !
    ! Overfill strategy support added by Jacek Dziedzic in March 2019.         !
    !==========================================================================!

    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_send, &
         comms_irecv, comms_test
    use constants, only: SW_D, SWEX_BB_CC_SYMMETRIC
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_probe, hash_table_add
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_assert, utils_abort

    implicit none

    ! jd: Arguments
    logical, intent(inout) :: who_is_done(0:pub_total_num_procs-1)
    type(HT_HASH_TABLE), intent(inout), target     :: coeffs_ht
    integer, intent(in)                            :: global_b
    integer, intent(in)                            :: global_c
    integer, intent(in)                            :: coeffs_kind
    integer, intent(in)                            :: swex_h
    type(FUNC_BASIS), intent(in), target           :: ngwf_basis
    type(PARAL_INFO), intent(in)                   :: par
    type(ELEMENT), intent(in)                      :: elements(par%nat)
    integer, intent(in)                            :: swri_h
    integer, intent(in), optional                  :: is
    type(FUNC_BASIS), intent(in), optional, target :: ngwf_basis2
    integer, intent(in), optional                  :: ket_basis_selector(3)
    character, intent(in), optional                :: overfill_strategy

    ! jd: Local variables
    type(FUNC_BASIS), pointer :: beta_ngwf_basis
    type(FUNC_BASIS), pointer :: gamma_ngwf_basis
    integer :: owner_proc
    integer :: send_handle1, send_handle2, recv_handle
    integer :: received_global_b, received_global_c
    integer :: first_ngwf_idx_of_b, first_ngwf_idx_of_c
    integer :: first_ngwf_idx_of_received_b, first_ngwf_idx_of_received_c
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx
    integer :: ngwf_b, ngwf_c
    integer :: orig_atom_b, orig_atom_c
    integer :: buf_start, buf_end
    integer :: number_of_coeffs_in_ngwf_pair
    integer :: n_elems
    logical :: cached
    logical :: ready
    logical :: done
    real(kind=DP), dimension(MAX_PACKED_COEFFS_SIZE) :: received_packed_coeffs
    character(len=*), parameter :: myself = 'remote_obtain_coeffs'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    call utils_assert(present(ngwf_basis2) .eqv. present(ket_basis_selector), &
         myself//': Optional arguments for HFx with mixed NGWF bases must &
         &either all be present or all be absent.')

    if(SWEX_BB_CC_SYMMETRIC(swex_h)) then
       call utils_assert(global_b <= global_c, &
            trim(myself)//': Where there is Bb-Cc symmetry between NGWFs, atom &
            &indices must be sorted on entry.')
    end if

    ! jd: Select the right destination NGWF_REP and its SW_EX for holding the
    !     expansion.
    if(present(ket_basis_selector)) then
       if(ket_basis_selector(1) == 1) then
          beta_ngwf_basis => ngwf_basis
       else if(ket_basis_selector(1) == 2) then
          beta_ngwf_basis => ngwf_basis2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(1)')
       end if
       if(ket_basis_selector(2) == 1) then
          gamma_ngwf_basis => ngwf_basis
       else if(ket_basis_selector(2) == 2) then
          gamma_ngwf_basis => ngwf_basis2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(2)')
       end if
    else
       beta_ngwf_basis => ngwf_basis
       gamma_ngwf_basis => ngwf_basis
    end if

    orig_atom_b = par%orig_atom(global_b)
    orig_atom_c = par%orig_atom(global_c)

    first_ngwf_idx_of_b = beta_ngwf_basis%first_on_atom(global_b)
    first_ngwf_idx_of_c = gamma_ngwf_basis%first_on_atom(global_c)

    ! jd: C coeffs have 5 keys (Bb, Cc, kind, swex_h, spin)
    !     D coeffs have 4 keys (Bb, Cc, kind, swex_h)
    if(coeffs_ht%n_keys == 5) then
       call utils_assert(present(is), &
            myself//': Inconsistent arguments [1]')
    else if(coeffs_ht%n_keys == 4) then
       call utils_assert(.not. present(is), &
            myself//': Inconsistent arguments [2]')
    else
       call utils_abort(myself//': Unexpected number of keys in (d)coeffs_ht', &
            coeffs_ht%n_keys)
    end if

    ! jd: In usual scenarios the owner of the coeffs is the proc of atom B.
    !     Only if
    !     1) B is not part of the SWRI, but C is,
    !     and
    !     2) we are not dealing with D-coeffs,
    !     then C is the owner.
    !     Condition 2) arises because D-coeffs do not use Bb<Cc ordering.
    if(elements(orig_atom_b)%in_swri(swri_h)) then
       owner_proc = beta_ngwf_basis%proc_of_func(&
            beta_ngwf_basis%first_on_atom(global_b))
    else if(elements(orig_atom_c)%in_swri(swri_h)) then
       if(coeffs_kind < SW_D) then
          owner_proc = gamma_ngwf_basis%proc_of_func(&
               gamma_ngwf_basis%first_on_atom(global_c))
       else
          owner_proc = beta_ngwf_basis%proc_of_func(&
               beta_ngwf_basis%first_on_atom(global_b))
       end if
    else
       call utils_abort(myself//': Looking for coeffs between two atoms neither&
            &of which is inside the SWRI. Global atom indices, SWRI &
            & handle and coeffs_kind follow.', global_b, global_c, swri_h, &
            coeffs_kind)
    end if

    if(owner_proc == pub_my_proc_id) then
       ! ***********************************************************************
       ! jd: An expansion I own, nothing to do, just ensure I really have it
       ! ***********************************************************************

       ! -----------------------------------------------------------------------
       ! jd: for all b on B                                                  bbb
       ! -----------------------------------------------------------------------
       loop_ngwf_b:                                                            &
       do ngwf_b = 1, beta_ngwf_basis%num_on_atom(global_b)
          global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

          ! --------------------------------------------------------------------
          ! jd: for all c on C                                               ccc
          ! --------------------------------------------------------------------
          loop_ngwf_c:                                                         &
          do ngwf_c = 1, gamma_ngwf_basis%num_on_atom(global_c)
             global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

             ! jd: For atoms B==C, we only need NGWFs b<=c, provided there is
             !     Bb-Cc symmetry and coeffs are not D-coeffs
             if(global_b == global_c .and. &
                 global_bb_ngwf_idx > global_cc_ngwf_idx .and. &
                 SWEX_BB_CC_SYMMETRIC(swex_h) .and. coeffs_kind < SW_D) cycle

             n_elems = hash_table_probe(coeffs_ht, global_bb_ngwf_idx, &
                  global_cc_ngwf_idx, coeffs_kind, swex_h, is) ! is can be absent

             call utils_assert(n_elems /= -1, &
                  trim(myself)//': (d)coeffs_ht misbehaves. SWEX coefficients &
                  &for a pair of NGWFs were not found. This is &
                  &indicative of a bug. Sorted NGWF indices and coeffs_kind &
                  &follow. ', global_bb_ngwf_idx, global_cc_ngwf_idx, &
                  coeffs_kind, swex_h)

          end do loop_ngwf_c
       end do loop_ngwf_b

    else
       ! ***********************************************************************
       ! jd: An expansion I don't own, see if cached.
       ! ***********************************************************************

       cached = .true.
       ! -----------------------------------------------------------------------
       ! jd: for all b on B                                                  bbb
       ! -----------------------------------------------------------------------
       loop_ngwf_b2:                                                           &
       do ngwf_b = 1, beta_ngwf_basis%num_on_atom(global_b)
          global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

          ! --------------------------------------------------------------------
          ! jd: for all c on C                                               ccc
          ! --------------------------------------------------------------------
          loop_ngwf_c2:                                                        &
          do ngwf_c = 1, gamma_ngwf_basis%num_on_atom(global_c)
             global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

             ! jd: For atoms B==C, we only need NGWFs b<=c, provided there is
             !     Bb-Cc symmetry and coeffs are not D-coeffs
             if(global_b == global_c .and. &
                 global_bb_ngwf_idx > global_cc_ngwf_idx .and. &
                 SWEX_BB_CC_SYMMETRIC(swex_h) .and. coeffs_kind < SW_D) cycle

             n_elems = hash_table_probe(coeffs_ht, global_bb_ngwf_idx, &
                  global_cc_ngwf_idx, coeffs_kind, swex_h, is) ! is can be absent

             if(n_elems == -1) then
                ! jd: At least one NGWF for atom pair B-C is absent, so get
                !     all the coeffs.
                cached = .false.
                exit loop_ngwf_b2
             end if

          end do loop_ngwf_c2
       end do loop_ngwf_b2

       if(cached) then
          ! jd: All cached, good. Nothing to do.
          call timer_clock(myself,2)
          return
       else
          ! ********************************************************************
          ! jd: Not owned or cached, request from remote proc.
          ! ********************************************************************
          call comms_send(owner_proc, global_b, 1, &
               tag = COEFFS_REQUEST_TAG, return_handle = send_handle1, &
               add_to_stack = .false.)
          call comms_send(owner_proc, global_c, 1, &
               tag = COEFFS_REQUEST_TAG, return_handle = send_handle2, &
               add_to_stack = .false.)

          ! jd: Initiate receive
          call comms_irecv(owner_proc, received_packed_coeffs, &
               tag = COEFFS_DATA_TAG, handle = recv_handle)

          ! jd: Wait for the receive to complete, satisfying the requests of
          !     others in the meantime.
          ready = .false.
          do while(.not. ready)

             call comms_test(ready,recv_handle)

             ! jd: Satisfy the requests of others
             call remote_serve_coeffs(done, who_is_done, ngwf_basis, coeffs_ht,&
                  coeffs_kind, swex_h, is, ngwf_basis2, ket_basis_selector) ! is can be absent

          end do

          ! ********************************************************************
          ! jd: Unpack what was received.
          ! ********************************************************************

          received_global_b = nint(received_packed_coeffs(1))
          received_global_c = nint(received_packed_coeffs(2))

          call utils_assert(received_global_b == global_b, &
               trim(myself)//': Received a different B atom than asked for', &
               global_b, received_global_b)
          call utils_assert(received_global_c == global_c, &
               trim(myself)//': Received a different C atom than asked for', &
               global_c, received_global_c)
          first_ngwf_idx_of_received_b = &
               beta_ngwf_basis%first_on_atom(received_global_b)
          first_ngwf_idx_of_received_c = &
               gamma_ngwf_basis%first_on_atom(received_global_c)
          number_of_coeffs_in_ngwf_pair = nint(received_packed_coeffs(3))

          buf_start = 4
          buf_end = buf_start + number_of_coeffs_in_ngwf_pair-1

          ! --------------------------------------------------------------------
          ! jd: for all b on B                                               bbb
          ! --------------------------------------------------------------------
          loop_ngwf_b3:                                                        &
          do ngwf_b = 1, beta_ngwf_basis%num_on_atom(received_global_b)
             global_bb_ngwf_idx = first_ngwf_idx_of_received_b + ngwf_b - 1

             ! -----------------------------------------------------------------
             ! jd: for all c on C                                            ccc
             ! -----------------------------------------------------------------
             loop_ngwf_c3:                                                     &
             do ngwf_c = 1, gamma_ngwf_basis%num_on_atom(received_global_c)
                global_cc_ngwf_idx = first_ngwf_idx_of_received_c + ngwf_c - 1

                ! jd: For atoms B==C, we only need NGWFs b<=c, provided there is
                !     Bb-Cc symmetry and coeffs are not D-coeffs
                if(global_b == global_c .and. &
                    global_bb_ngwf_idx > global_cc_ngwf_idx .and. &
                    SWEX_BB_CC_SYMMETRIC(swex_h) .and. coeffs_kind < SW_D) cycle

                ! jd: Store in cache
                call hash_table_add(coeffs_ht, &
                     received_packed_coeffs(buf_start:buf_end), &
                     number_of_coeffs_in_ngwf_pair, global_bb_ngwf_idx, &
                     global_cc_ngwf_idx, coeffs_kind, swex_h, is, & ! is can be absent
                     overfill_strategy = overfill_strategy)

                buf_start = buf_start + number_of_coeffs_in_ngwf_pair
                buf_end = buf_end + number_of_coeffs_in_ngwf_pair

             end do loop_ngwf_c3
          end do loop_ngwf_b3

       end if
    end if

    call timer_clock(myself,2)

  end subroutine remote_obtain_coeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_serve_coeffs(done, who_is_done, ngwf_basis, coeffs_ht, &
       coeffs_kind, swex_h, is, ngwf_basis2, ket_basis_selector)
    !==========================================================================!
    ! Keeps probing for requests from remote procs asking for expansion coeffi-!
    ! cients and satisfies them by sending requested entries from coeffs_ht.   !
    ! 'coeffs_ht' can contain either of:                                       !
    ! a) spin-agnostic, but formally spin-dependent SW coefficients (in HFx    !
    !    in DMA in mode 'P') -- then 'is' is always 1, even in spin-polarised  !
    !    systems;                                                              !
    ! b) truly spin-dependent SW coefficients (in DMA in mode 'D') -- then 'is'!
    !    will refer to the real spin channel;                                  !
    ! c) truly spin-agnostic SW coefficients (D coefficients in HFx and in DMA !
    !    regardless of mode) -- then 'is' needs to be omitted from the argument!
    !    list.
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   done (out): Set to .true. if all procs sent a 'done' notification.     !
    !   who_is_done (in/out): An array of logicals for tracking which procs    !
    !                         sent a 'done' notification.                      !
    !   ngwf_basis (in): Required to find out owners of NGWFs.                 !
    !   coeffs_ht (in): The hash table from which coefficients are served.     !
    !   coeffs_kind (in): The type of the coefficients to serve (SW_V, SW_O,   !
    !                     SW_V+SW_D, SW_O+SW_D).                               !
    !   swex_h (in): Identifier for distinguishing multiple sets of            !
    !                coefficients (e.g. in Bessel averaging). Pass SWEX handle !
    !                here, eg. REP_SWEX_PROPERTIES_DMA_1. This value is used as!
    !                a key for looking up and storing coeffs in a hash table.  !
    !                It is also used to determine if we have Bb/Cc symmetry.   !
    !   is (in, opt): Spin channel for the coefficients -- only needed in      !
    !                 cases a) and b) above.                                   !
    !                                                                          !
    !   *Optional* arguments for mixed NGWF bases.                             !
    !   ngwf_basis2 (in): Second NGWF basis, for finding out owners of NGWFs.  !
    !   ket_basis_selector (in): Selects the NGWF basis in which               !
    !                            - Bb is indexed (ket_basis_selector(1)),      !
    !                            - Cc is indexed (ket_basis_selector(2)).      !
    !                            ket_basis_selector(3) is ignored here.        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2014.                               !
    ! Improved by Jacek Dziedzic in July 2016 to batch coeffs over all NGWF    !
    ! pairs on an atom pair.                                                   !
    ! Generalized by Jacek Dziedzic in June 2017 to spin-dependent coeffs.     !
    ! Mixed NGWF set support added by Jacek Dziedzic in July, August 2018.     !
    !==========================================================================!

    use comms, only: pub_total_num_procs, comms_send, comms_recv, comms_probe, &
         comms_wait
    use constants, only: SWEX_BB_CC_SYMMETRIC, SW_D
    use function_basis, only: FUNC_BASIS
    use hash_table, only: HT_HASH_TABLE, hash_table_lookup, hash_table_probe
    use utils, only: utils_assert
    use utils, only: utils_abort

    implicit none

    ! Arguments
    logical, intent(out)                           :: done
    logical, intent(inout)                         :: who_is_done(0:pub_total_num_procs-1)
    type(FUNC_BASIS), intent(in), target           :: ngwf_basis
    type(HT_HASH_TABLE), intent(inout), target     :: coeffs_ht
    integer, intent(in)                            :: coeffs_kind
    integer, intent(in)                            :: swex_h
    integer, intent(in), optional                  :: is
    type(FUNC_BASIS), intent(in), optional, target :: ngwf_basis2
    integer, intent(in), optional                  :: ket_basis_selector(3)

    ! Local variables
    logical :: who_wants_coeffs(0:pub_total_num_procs-1)
    integer :: who_wants_which_coeffs_b(0:pub_total_num_procs-1)
    integer :: who_wants_which_coeffs_c(0:pub_total_num_procs-1)
    integer :: coeffs_to_send_b, coeffs_to_send_c
    integer :: requested_coeffs_b, requested_coeffs_c
    integer :: first_ngwf_idx_of_b, first_ngwf_idx_of_c
    integer :: global_bb_ngwf_idx, global_cc_ngwf_idx
    integer :: ngwf_b, ngwf_c
    integer :: buf_start, buf_end
    integer :: number_of_coeffs_to_send
    integer :: number_of_coeffs_in_ngwf_pair
    integer :: send_handle
    integer :: proc
    real(kind=DP), dimension(MAX_PACKED_COEFFS_SIZE) :: packed_coeffs_to_send
    type(FUNC_BASIS), pointer :: beta_ngwf_basis
    type(FUNC_BASIS), pointer :: gamma_ngwf_basis
    character(len=*), parameter :: myself = 'remote_serve_coeffs'

    ! -------------------------------------------------------------------------

    ! Do not time this subroutine. All you'll get will be timer overhead.

#ifndef MPI
    call utils_abort('remote_serve_coeffs: Calling me is pointless in serial &
         &mode')
#endif
    ! jd: Ensure consistency between optional parameters for mixed NGWF bases
    call utils_assert(present(ngwf_basis2) .eqv. present(ket_basis_selector), &
         myself//': Optional arguments for HFx with mixed NGWF bases must &
         &either all be present or all be absent.')

    ! jd: Select the right destination NGWF_REP and its SW_EX for holding the
    !     expansion.
    if(present(ket_basis_selector)) then
       if(ket_basis_selector(1) == 1) then
          beta_ngwf_basis => ngwf_basis
       else if(ket_basis_selector(1) == 2) then
          beta_ngwf_basis => ngwf_basis2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(1)')
       end if
       if(ket_basis_selector(2) == 1) then
          gamma_ngwf_basis => ngwf_basis
       else if(ket_basis_selector(2) == 2) then
          gamma_ngwf_basis => ngwf_basis2
       else
          call utils_abort(myself//': Illegal ket_basis_selector(2)')
       end if
    else
       beta_ngwf_basis => ngwf_basis
       gamma_ngwf_basis => ngwf_basis
    end if

    ! jd: Probe for requests from others
    who_wants_coeffs(:) = .false.
    do proc = 0, pub_total_num_procs - 1
       call comms_probe(who_wants_coeffs(proc), proc, COEFFS_REQUEST_TAG)
    end do

    ! jd: Find out what coeffs do the requests refer to, by actually
    !     receiving the requests
    who_wants_which_coeffs_b(:) = -1
    who_wants_which_coeffs_c(:) = -1
    do proc = 0, pub_total_num_procs - 1
       if(who_wants_coeffs(proc)) then

          call comms_recv(proc, requested_coeffs_b, 1, COEFFS_REQUEST_TAG)
          call comms_recv(proc, requested_coeffs_c, 1, COEFFS_REQUEST_TAG)

          who_wants_which_coeffs_b(proc) = requested_coeffs_b
          who_wants_which_coeffs_c(proc) = requested_coeffs_c

          ! jd: Is it a 'done' notification?
          if(requested_coeffs_b == -1) then
             who_is_done(proc) = .true.
             who_wants_coeffs(proc) = .false.
          end if
       end if
    end do

    done = all(who_is_done)

    ! jd: Satisfy the requests from other procs
    do proc = 0, pub_total_num_procs - 1
       if(who_wants_coeffs(proc)) then
          coeffs_to_send_b = who_wants_which_coeffs_b(proc)
          coeffs_to_send_c = who_wants_which_coeffs_c(proc)
          number_of_coeffs_in_ngwf_pair = -1

          ! jd: Prepare a batch of coeffs for atom pair B-C. The batch
          !     contains coeffs for all NGWFs b, c.
          packed_coeffs_to_send(1) = coeffs_to_send_b
          packed_coeffs_to_send(2) = coeffs_to_send_c
          first_ngwf_idx_of_b = beta_ngwf_basis%first_on_atom(coeffs_to_send_b)
          first_ngwf_idx_of_c = gamma_ngwf_basis%first_on_atom(coeffs_to_send_c)
          number_of_coeffs_to_send = 0

          ! -----------------------------------------------------------------
          ! jd: for all b on B                                            bbb
          ! -----------------------------------------------------------------
          loop_ngwf_b:                                                      &
          do ngwf_b = 1, beta_ngwf_basis%num_on_atom(coeffs_to_send_b)
             global_bb_ngwf_idx = first_ngwf_idx_of_b + ngwf_b - 1

             ! --------------------------------------------------------------
             ! jd: for all c on C                                         ccc
             ! --------------------------------------------------------------
             loop_ngwf_c:                                                   &
             do ngwf_c = 1, gamma_ngwf_basis%num_on_atom(coeffs_to_send_c)
                global_cc_ngwf_idx = first_ngwf_idx_of_c + ngwf_c - 1

                ! jd: For atoms B==C, we only need NGWFs b<=c, provided there is
                !     Bb-Cc symmetry and coeffs are not D-coeffs
                if(coeffs_to_send_b == coeffs_to_send_c .and. &
                    global_bb_ngwf_idx > global_cc_ngwf_idx .and. &
                    SWEX_BB_CC_SYMMETRIC(swex_h) .and. coeffs_kind < SW_D) cycle

                number_of_coeffs_in_ngwf_pair = &
                     hash_table_probe(coeffs_ht, global_bb_ngwf_idx, &
                     global_cc_ngwf_idx, coeffs_kind, swex_h, is)
                call utils_assert(number_of_coeffs_in_ngwf_pair > 0, myself//&
                     ': Bad number of coeffs', number_of_coeffs_in_ngwf_pair)

                buf_start = 4+number_of_coeffs_to_send
                number_of_coeffs_to_send = &
                     number_of_coeffs_to_send + number_of_coeffs_in_ngwf_pair
                buf_end = 4+number_of_coeffs_to_send-1

                ! jd: Ensure we're not exceeding the fixed-size send buffer here
                !     and similarly for recv buffer at the recv end.
                call utils_assert(buf_end <= MAX_PACKED_COEFFS_SIZE, &
                     myself//': Insufficient MAX_PACKED_COEFFS_SIZE.', &
                     buf_end, MAX_PACKED_COEFFS_SIZE)

                call hash_table_lookup(&
                     packed_coeffs_to_send(buf_start:buf_end), &
                     number_of_coeffs_in_ngwf_pair, &
                     coeffs_ht, global_bb_ngwf_idx, global_cc_ngwf_idx, &
                     coeffs_kind, swex_h, is)

             end do loop_ngwf_c
          end do loop_ngwf_b
          if(number_of_coeffs_in_ngwf_pair == -1) then
             call utils_abort(myself//': Atom with 0 NGWFs or logic error.')
          end if
          packed_coeffs_to_send(3) = number_of_coeffs_in_ngwf_pair

          ! jd: Batch of coeffs prepared. Now send.
          call comms_send(proc, packed_coeffs_to_send, &
               number_of_coeffs_to_send+3, &
               COEFFS_DATA_TAG, send_handle, add_to_stack = .false.)
          ! @optimizeme
          ! The above only initiates a send, we need to wait before we can
          ! overwrite the buffer with a new packed_coeffs_to_send. If this leads
          ! to convoy effects, consider having an array of packed_coeffs here
          ! and testing in a loop at the end. This will be more efficient, but
          ! comes with a memory penalty.
          call comms_wait(send_handle)

       end if
    end do

  end subroutine remote_serve_coeffs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_unpack_ngwf_to_fftbox(fftbox_ngwf, packed_ngwf, fftbox, &
       cell, centre_for_target_fftbox, fftbox_start_in_cell_vec, &
       delta_mid_complement, gridpoint_offset)
    !==========================================================================!
    ! Unpacks a packed NGWF to an FFTbox. The returned FFTbox data is allocated!
    ! here and will have to be freed by the caller. This subroutine is only    !
    ! used in a highly-experimental codepath where pub_swx_dbl_grid is .true.  !
    !                                                                          !
    ! A number of optional arguments modify the behaviour of this subroutine   !
    ! in complicated ways (see below). Not for the faint of heart.             !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   fftbox_ngwf (out): This will undergo allocation, and result will be    !
    !                      returned here.                                      !
    !   packed_ngwf (in): The NGWF to be unpacked, in packed representation.   !
    !   fftbox (in): Describes the FFTbox.                                     !
    !   cell (in): The usual.                                                  !
    !                                                                          !
    ! Optional arguments: (only used if passed)                                !
    !                                                                          !
    !   centre_for_target_fftbox (in): If this POINT is specified, the returned!
    !                                  FFTbox will not be centred on the NGWF  !
    !                                  that is unpacked, but on a point        !
    !                                  *mid way* between the centre of the NGWF!
    !                                  and centre_for_target_fftbox.           !
    !                                  This is needed for constructing FFTboxes!
    !                                  around an intersection of 2 NGWFs.      !
    !                                  Due to grid granularity, the FFTbox will!
    !                                  not be centred perfectly on the mid way !
    !                                  point, but at the nearest gridpoint.    !
    !   fftbox_start_in_cell_vec (out): This will return the upper left corner !
    !                                   of the FFTbox of the unpacked NGWF.    !
    !                                   It is unaffected by centre_for_target_ !
    !                                   fftbox.                                !
    !   delta_mid_complement (out): If specified, together with centre_for_tar !
    !                               get_fftbox, it will return the offset, in  !
    !                               grid point fractions ('gp') between the    !
    !                               gridpoint corresponding to the mid way     !
    !                               point and one of the original FFTboxes.    !
    !   gridpoint_offset (in): A different way of offseting where the target   !
    !                          is unpacked to. This is measured in grid points,!
    !                          is a direct offset (no mid-way tricks), and is  !
    !                          taken into account in fftbox_start_in_cell_vec. !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in December 2016.                              !
    ! Parts of code adapted from function_ops_batch_col_start due to Nick Hine.!
    !==========================================================================!

    use datatypes, only: FFTBOX_DATA, FUNCTIONS, data_fftbox_alloc, &
         data_functions_dealloc
    use basis, only: SPHERE, FUNCTION_TIGHT_BOX, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell, basis_location_func_wrt_cell, &
         basis_copy_function_to_box, basis_sphere_deallocate
    use fft_box, only: FFTBOX_INFO
    use geometry, only: POINT, operator(*), operator(+)
    use simulation_cell, only: CELL_INFO

    implicit none

    ! jd: gp -- gridpoints, vec -- vector in real space

    ! Arguments
    type(FFTBOX_DATA), intent(out)     :: fftbox_ngwf
    real(kind=DP), intent(in)          :: packed_ngwf(packed_ngwf_size)
    type(FFTBOX_INFO), intent(in)      :: fftbox
    type(CELL_INFO), intent(in)        :: cell
    type(POINT), intent(in), optional  :: centre_for_target_fftbox
    type(POINT), intent(out), optional :: fftbox_start_in_cell_vec
    integer, intent(out), optional :: delta_mid_complement(3)
    integer, intent(in), optional  :: gridpoint_offset(3)

    ! Local variables
    type(SPHERE)             :: ngwf_sphere
    type(FUNCTIONS)          :: ngwfs_on_grid
    type(FUNCTION_TIGHT_BOX) :: ngwf_tightbox
    integer                  :: fftbox_start_in_cell_gp(3)
    integer                  :: target_fftbox_start_in_cell_gp(3)
    integer                  :: tb_start_in_cell_gp(3)
    integer                  :: tb_start_in_fftbox_gp(3)
    integer                  :: delta_mid(3)
    ! {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    real(kind=DP)            :: cbg1, cbg2, cbg3

    ! -------------------------------------------------------------------------

    ! jd: 1) --- Unpack NGWF to PPDs and sphere ---
    call remote_unpack_ngwf(packed_ngwf, ngwf_sphere, ngwfs_on_grid, &
         ngwf_tightbox, dest_offset = 1)

    ! jd: 2) --- Do what function_ops_batch_col_start() does, but     ---
    !        --- for a single sphere, potentially for a remote NGWF.  ---

    ! Centre (of unpacked ngwf sphere) wrt fftbox in terms of grid spacings
    call basis_point_wrt_box(cbg1, cbg2, cbg3, &  ! output
         ngwf_sphere%centre, fftbox%total_pt1, &
         fftbox%total_pt2, fftbox%total_pt3, cell, fftbox)

    ! Start of fftbox wrt cell in terms of grid-point number
    call basis_start_of_box_wrt_cell(fftbox_start_in_cell_gp(1), &
         fftbox_start_in_cell_gp(2), fftbox_start_in_cell_gp(3), &
         ngwf_sphere%centre, cbg1, cbg2, cbg3, cell)

    if(present(centre_for_target_fftbox)) then
       ! Centre (of target fftbox) wrt its fftbox in terms of grid spacings
       call basis_point_wrt_box(cbg1, cbg2, cbg3, &  ! output
            centre_for_target_fftbox, fftbox%total_pt1, &
            fftbox%total_pt2, fftbox%total_pt3, cell, fftbox)
       ! Start of target fftbox wrt cell in terms of grid-point number
       call basis_start_of_box_wrt_cell(target_fftbox_start_in_cell_gp(1), &
            target_fftbox_start_in_cell_gp(2), &
            target_fftbox_start_in_cell_gp(3), &
            centre_for_target_fftbox, cbg1, cbg2, cbg3, cell)
    end if

    ! jd: Tightbox of unpacked NGWF wrt start of simulation cell, in gridpoints
    call basis_location_func_wrt_cell(tb_start_in_cell_gp(1), &
         tb_start_in_cell_gp(2), tb_start_in_cell_gp(3), ngwf_tightbox, cell)

    ! jd: Tightbox of unpacked NGWF wrt its fftbox, in gridpoints
    tb_start_in_fftbox_gp(1) = tb_start_in_cell_gp(1) - fftbox_start_in_cell_gp(1) + 1
    tb_start_in_fftbox_gp(2) = tb_start_in_cell_gp(2) - fftbox_start_in_cell_gp(2) + 1
    tb_start_in_fftbox_gp(3) = tb_start_in_cell_gp(3) - fftbox_start_in_cell_gp(3) + 1

    ! jd: If unpacking to a different target FFT box, offset accordingly
    if(present(centre_for_target_fftbox)) then
       delta_mid = (target_fftbox_start_in_cell_gp - fftbox_start_in_cell_gp) / 2
       if(present(delta_mid_complement)) then
          delta_mid_complement = &
               (target_fftbox_start_in_cell_gp - fftbox_start_in_cell_gp) - delta_mid
       end if
       tb_start_in_fftbox_gp = tb_start_in_fftbox_gp - delta_mid
    end if

    if(present(gridpoint_offset)) then
       tb_start_in_fftbox_gp = tb_start_in_fftbox_gp + gridpoint_offset
    end if

    ! jd: 3) --- Allocate returned FFT box ---
    call data_fftbox_alloc(fftbox_ngwf, fftbox%total_ld1, fftbox%total_ld2, &
         fftbox%total_pt3,ngwfs_on_grid%iscmplx)

    ! jd: 4) --- Transfer the data from PPDs to FFT box ---
    call basis_copy_function_to_box(fftbox_ngwf, & ! target
         tb_start_in_fftbox_gp(1), tb_start_in_fftbox_gp(2), &
         tb_start_in_fftbox_gp(3), ngwf_tightbox, ngwfs_on_grid, ngwf_sphere, &
         cell, fftbox, 1)

    ! jd: 5) --- Deallocate NGWFs on grid and sphere
    call basis_sphere_deallocate(ngwf_sphere)
    call data_functions_dealloc(ngwfs_on_grid)

    if(present(fftbox_start_in_cell_vec)) then
       fftbox_start_in_cell_vec = &
            real(fftbox_start_in_cell_gp(1)-1 - &
            gridpoint_offset(1),kind=DP) * fftbox%d1 * fftbox%a1_unit+&
            real(fftbox_start_in_cell_gp(2)-1 - &
            gridpoint_offset(2),kind=DP) * fftbox%d2 * fftbox%a2_unit+&
            real(fftbox_start_in_cell_gp(3)-1 - &
            gridpoint_offset(3),kind=DP) * fftbox%d3 * fftbox%a3_unit
    end if

  end subroutine remote_unpack_ngwf_to_fftbox

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_ppd_list_of_atom(ppd_list, n_ppds, & ! out
       global_atom, ngwf_basis, cache_handle)            ! in
    !==========================================================================!
    ! Returns the list of PPD indices for a given atom (global index).         !
    ! The atom does not need to be local to the proc, the data is extracted    !
    ! from the NGWF cache. The cache is expected to contain the atom already.  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ppd_list (out): On exit the first n_ppds elements will contain the PPD !
    !                   indices of the atom. Remaining elements will be -1.    !
    !   n_ppds (out): On exit will contain the number of PPDs on the atom.     !
    !   ngwf_basis (in): The NGWF basis in which the atom lives.               !
    !   global_atom (in): The global atom index whose PPD list we want.        !
    !   cache_handle (in): The handle to the NGWF cache for ngwf_basis.        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2019.                                !
    !==========================================================================!

    use function_basis, only: FUNC_BASIS
    use hash_table, only: hash_table_lookup_ptr
    use utils, only: utils_assert, utils_int_to_str

    implicit none

    ! jd: Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    integer, intent(out)         :: ppd_list(ngwf_basis%max_n_ppds_sphere)
    integer, intent(out)         :: n_ppds
    integer, intent(in)          :: global_atom
    integer, intent(in)          :: cache_handle

    ! jd: Local variables
    real(kind=DP), pointer :: packed_ngwf_ptr(:) ! ptr to ht's internals
    integer :: packed_ngwf_size
    integer :: first_ngwf_idx_of_atom
    character(len=*), parameter :: myself = 'remote_ppd_list_of_atom'

    ! -------------------------------------------------------------------------

    ! NB. This subroutine is called from OMP contexts. Do not utils_trace() it.

    first_ngwf_idx_of_atom = ngwf_basis%first_on_atom(global_atom)

    call hash_table_lookup_ptr(packed_ngwf_ptr, packed_ngwf_size, &
         remote_ngwf_cache_hts(cache_handle), first_ngwf_idx_of_atom)
    call utils_assert(packed_ngwf_size /= -1, &
         myself//': NGWF '//trim(utils_int_to_str(first_ngwf_idx_of_atom))//&
         ' (1st NGWF of atom '//trim(utils_int_to_str(global_atom))//&
         ') not present in NGWF cache (handle '//&
         trim(utils_int_to_str(cache_handle))//&
         '). Insufficient NGWF cache size or a bug.')
    call remote_unpack_ngwf_no_sphere(packed_ngwf_ptr, ppd_list, n_ppds)
    ppd_list(n_ppds+1:) = -1

  end subroutine remote_ppd_list_of_atom

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine remote_n_ppds_of_atom(n_ppds, global_atom, ngwf_basis, cache_handle)            ! in
    !==========================================================================!
    ! Returns the number of PPDs for a given atom (global index).              !
    ! The atom does not need to be local to the proc, the data is extracted    !
    ! from the NGWF cache. The cache is expected to contain the atom already.  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   n_ppds (out): On exit will contain the number of PPDs on the atom.     !
    !   ngwf_basis (in): The NGWF basis in which the atom lives.               !
    !   global_atom (in): The global atom index whose PPD list we want.        !
    !   cache_handle (in): The handle to the NGWF cache for ngwf_basis.        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in March 2020.                                 !
    !==========================================================================!

    use function_basis, only: FUNC_BASIS
    use hash_table, only: hash_table_lookup_ptr
    use utils, only: utils_assert, utils_int_to_str

    implicit none

    ! jd: Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    integer, intent(out)         :: n_ppds
    integer, intent(in)          :: global_atom
    integer, intent(in)          :: cache_handle

    ! jd: Local variables
    integer                :: ppd_list(ngwf_basis%max_n_ppds_sphere)
    real(kind=DP), pointer :: packed_ngwf_ptr(:) ! ptr to ht's internals
    integer :: packed_ngwf_size
    integer :: first_ngwf_idx_of_atom
    character(len=*), parameter :: myself = 'remote_n_ppds_atom'

    ! -------------------------------------------------------------------------

    ! NB. This subroutine is called from OMP contexts. Do not utils_trace() it.

    first_ngwf_idx_of_atom = ngwf_basis%first_on_atom(global_atom)

    call hash_table_lookup_ptr(packed_ngwf_ptr, packed_ngwf_size, &
         remote_ngwf_cache_hts(cache_handle), first_ngwf_idx_of_atom)
    call utils_assert(packed_ngwf_size /= -1, &
         myself//': NGWF '//trim(utils_int_to_str(first_ngwf_idx_of_atom))//&
         ' (1st NGWF of atom '//trim(utils_int_to_str(global_atom))//&
         ') not present in NGWF cache (handle '//&
         trim(utils_int_to_str(cache_handle))//&
         '). Insufficient NGWF cache size or a bug.')
    call remote_unpack_ngwf_no_sphere(packed_ngwf_ptr, ppd_list, n_ppds)

  end subroutine remote_n_ppds_of_atom

end module remote
