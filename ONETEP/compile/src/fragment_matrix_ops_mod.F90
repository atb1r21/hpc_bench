! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                   F R A G M E N T   M A T R I X                             !
!                              M O D U L E                                    !
!=============================================================================!
!                                                                             !
! This module contains the matrix manipulation routines used by the EDA code. !
!                                                                             !
!-----------------------------------------------------------------------------!
! Originally written by Max Phipps August 2014                                !
! Modified by Max Phipps December 2015 to improve efficiency and module       !
!      structure                                                              !
!-----------------------------------------------------------------------------!

module fragment_matrix_ops

  use constants, only: DP
  use dense, only: DEM
  use fragment_data, only: pub_super2frag_ngwf_log
  use rundat, only: PUB_1K
  use sparse_embed, only: SPAM3_EMBED

  implicit none

  private

  ! subroutines:

  ! full matrix blocking routines:
  public :: fmo_freeze_mat_spam_nxn

  ! fragment block extraction routines:
  public :: fmo_copy_frag_block_spam_nxn
  public :: fmo_invert_frag_block_dens_nxn


  ! fragment extraction and
  ! mask construction routines:
  public :: fmo_get_frag_blk_dens_nxn
  public :: fmo_get_frag_blk_dens_mxn
  public :: fmo_get_frag_blk_dens_nxm
  public :: fmo_get_frag_blk_dens_mxm
  public :: fmo_put_frag_blk_dens_nxn

  public :: fmo_construct_masks
  public :: fmo_allocate_masks
  public :: fmo_destroy_masks

  ! supermolecule masks:
  type DEM_FRAG_MASKS
     type(DEM)              :: nxn
     type(DEM), allocatable :: nxm(:) ! (spin)
     type(DEM), allocatable :: mxn(:) ! (spin)
     type(DEM), allocatable :: mxm(:) ! (spin)
     type(DEM), allocatable :: nxn_r(:) ! (fragment)
                                      ! special, for row-fragmented denskern
     type(DEM), allocatable :: nxn_c(:) ! (fragment)
                                      ! special, for column fragmentations

     ! SPAM3 versions of above matrices:
     type(SPAM3_EMBED)              :: nxn_sH
     type(SPAM3_EMBED)              :: nxn_sK
     type(SPAM3_EMBED)              :: nxn_sHKS
     type(SPAM3_EMBED)              :: nxn_sSKHKS
     type(SPAM3_EMBED), allocatable :: nxn_scK(:) ! (fragment)
     type(SPAM3_EMBED), allocatable :: nxn_srK(:) ! (fragment)
     type(SPAM3_EMBED), allocatable :: nxn_scHKS(:) ! (fragment)
     type(SPAM3_EMBED), allocatable :: nxn_srHKS(:) ! (fragment)
     type(SPAM3_EMBED), allocatable :: nxn_scSKHKS(:) ! (fragment)
     type(SPAM3_EMBED), allocatable :: nxn_srSKHKS(:) ! (fragment)

  end type DEM_FRAG_MASKS

  type(DEM_FRAG_MASKS), public :: fmo_masks


  ! auxiliary methods:
  public :: fmo_load_internal_blkdiag_denskern
  public :: fmo_apply_intrafrag_partitions
  public :: fmo_construct_overlap_mats


contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fmo_allocate_masks()

    !========================================================================!
    ! This subroutine allocates the fmo_mask matrices used to extract        !
    ! the fragment matrices from a supermolecule matrix.                     !
    !------------------------------------------------------------------------!
    ! Written by Max Phipps January 2016.                                    !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    use dense, only: dense_create
    use fragment_data, only: pub_frag_data
    use rundat, only: pub_eda_nodiag, pub_frag_iatm, pub_num_spins
    use sparse_embed, only: sparse_embed_create, sparse_embed_destroy
    use utils, only: utils_alloc_check

    implicit none

    integer     :: is ! spin counter
    integer     :: it ! fragment counter
    integer     :: ierr
    type(SPAM3_EMBED) :: SK_tmp

    ! --------------------------------------------------------------------------

    ! mjsp: Allocations:
    allocate(fmo_masks%nxm(pub_num_spins),stat=ierr)
    call utils_alloc_check('fmo_allocate_masks','fmo_masks%nxm',ierr)
    allocate(fmo_masks%mxn(pub_num_spins),stat=ierr)
    call utils_alloc_check('fmo_allocate_masks','fmo_masks%mxn',ierr)
    allocate(fmo_masks%mxm(pub_num_spins),stat=ierr)
    call utils_alloc_check('fmo_allocate_masks','fmo_masks%mxm',ierr)
    if (pub_eda_nodiag .gt. 0) then
       ! Dense:
       allocate(fmo_masks%nxn_r(pub_frag_iatm(0)),stat=ierr)
       call utils_alloc_check('fmo_allocate_masks','fmo_masks%nxn_r',ierr)
       allocate(fmo_masks%nxn_c(pub_frag_iatm(0)),stat=ierr)
       call utils_alloc_check('fmo_allocate_masks','fmo_masks%nxn_c',ierr)

       ! Sparse:
       allocate(fmo_masks%nxn_srK(pub_frag_iatm(0)),stat=ierr)
       call utils_alloc_check('fmo_allocate_masks','fmo_masks%nxn_srK',ierr)
       allocate(fmo_masks%nxn_scK(pub_frag_iatm(0)),stat=ierr)
       call utils_alloc_check('fmo_allocate_masks','fmo_masks%nxn_scK',ierr)
       allocate(fmo_masks%nxn_srHKS(pub_frag_iatm(0)),stat=ierr)
       call utils_alloc_check('fmo_allocate_masks','fmo_masks%nxn_srHKS',ierr)
       allocate(fmo_masks%nxn_scHKS(pub_frag_iatm(0)),stat=ierr)
       call utils_alloc_check('fmo_allocate_masks','fmo_masks%nxn_scHKS',ierr)
       allocate(fmo_masks%nxn_srSKHKS(pub_frag_iatm(0)),stat=ierr)
       call utils_alloc_check('fmo_allocate_masks','fmo_masks%nxn_srSKHKS',ierr)
       allocate(fmo_masks%nxn_scSKHKS(pub_frag_iatm(0)),stat=ierr)
       call utils_alloc_check('fmo_allocate_masks','fmo_masks%nxn_scSKHKS',ierr)
    end if

    ! ======= nxn ======
    ! mjsp: supermolecule nxn (build from fragment masks)
    call dense_create(fmo_masks%nxn,pub_frag_data(0)%ngwf_basis(1)%num,&
                                      pub_frag_data(0)%ngwf_basis(1)%num)
    fmo_masks%nxn_sH%structure = 'H'
    fmo_masks%nxn_sK%structure = 'K'
    fmo_masks%nxn_sHKS%structure = 'HKS'
    SK_tmp%structure = 'SK'
    call sparse_embed_create(fmo_masks%nxn_sH)
    call sparse_embed_create(fmo_masks%nxn_sK)
    call sparse_embed_create(fmo_masks%nxn_sHKS)
    call sparse_embed_create(SK_tmp)
    call sparse_embed_create(fmo_masks%nxn_sSKHKS, SK_tmp, fmo_masks%nxn_sHKS)
    call sparse_embed_destroy(SK_tmp)

    ! ======= nxn_r, nxn_c ======
    ! mjsp: supermolecule nxn_r (where r indicates fragmentation along rows)
    ! mjsp: and nxn_c (where c indicates fragmentation along columns)
    ! mjsp: (build from fragment masks)
    if (pub_eda_nodiag .gt. 0) then
       do it=1,pub_frag_iatm(0) ! loop fragments

          ! Dense:
          call dense_create(fmo_masks%nxn_r(it), &
               pub_frag_data(0)%ngwf_basis(1)%num,&
               pub_frag_data(0)%ngwf_basis(1)%num)

          call dense_create(fmo_masks%nxn_c(it), &
               pub_frag_data(0)%ngwf_basis(1)%num,&
               pub_frag_data(0)%ngwf_basis(1)%num)

          ! Sparse:
          fmo_masks%nxn_srK(it)%structure = 'K'
          fmo_masks%nxn_scK(it)%structure = 'K'
          fmo_masks%nxn_srHKS(it)%structure = 'HKS'
          fmo_masks%nxn_scHKS(it)%structure = 'HKS'
          fmo_masks%nxn_srSKHKS(it)%structure = 'SKHKS'
          fmo_masks%nxn_scSKHKS(it)%structure = 'SKHKS'
          call sparse_embed_create(fmo_masks%nxn_srK(it))
          call sparse_embed_create(fmo_masks%nxn_scK(it))
          call sparse_embed_create(fmo_masks%nxn_srHKS(it))
          call sparse_embed_create(fmo_masks%nxn_scHKS(it))
          call sparse_embed_create(fmo_masks%nxn_srSKHKS(it))
          call sparse_embed_create(fmo_masks%nxn_scSKHKS(it))

       end do ! fragment
    end if

    ! ======= nxm ======
    do is=1,pub_num_spins

       ! mjsp: supermolecule nxm (build from fragment masks)
       call dense_create(fmo_masks%nxm(is), &
                         pub_frag_data(0)%ngwf_basis(1)%num,&
                         pub_frag_data(0)%rep%n_occ(is,PUB_1K))

    end do ! spin

    ! ======= mxn ======
    do is=1,pub_num_spins

       ! mjsp: supermolecule mxn (build from fragment masks)
       call dense_create(fmo_masks%mxn(is), &
                         pub_frag_data(0)%rep%n_occ(is,PUB_1K),&
                         pub_frag_data(0)%ngwf_basis(1)%num)

    end do ! spin

    ! ======= mxm ======
    do is=1,pub_num_spins

       ! mjsp: supermolecule mxm (build from fragment masks)
       call dense_create(fmo_masks%mxm(is), &
                         pub_frag_data(0)%rep%n_occ(is,PUB_1K),&
                         pub_frag_data(0)%rep%n_occ(is,PUB_1K))

    end do ! spin

  end subroutine fmo_allocate_masks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_construct_masks()

    !========================================================================!
    ! This subroutine constructs the fmo_mask object's matrices that are     !
    ! used to extract the fragment matrices from a supermolecule matrix.     !
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    ! Written by Max Phipps January 2016.                                    !
    !========================================================================!

    use constants, only: EDA_POLFRAGLOC_DEVEL, EDA_CTFRAGLOC_DEVEL
    use dense, only: dense_create, dense_destroy, dense_axpy, dense_scale, &
         dense_transpose, dense_convert
    use fragment_data, only: pub_frag_data
    use rundat, only: pub_eda_nodiag, pub_frag_iatm, pub_num_spins, &
         pub_frag_counter, pub_frag_counter2, pub_eda_mode
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(DEM)                            :: mask_buf_nxn
    type(DEM), allocatable, dimension(:) :: mask_buf_nxm, mask_buf_mxn, mask_buf_mxm
    integer :: it, it2 ! fragment counter
    integer :: is ! spin counter
    integer :: ierr

    ! Start Timer
    call timer_clock('fragment_construct_masks',1)

    ! --------------------------------------------------------------------------


    ! mjsp: Initialisations:
    it2 = 0

    ! mjsp: (re-)initialise the masks:
    call dense_scale(fmo_masks%nxn,0.0_DP)
    do is=1,pub_num_spins
       call dense_scale(fmo_masks%nxm(is),0.0_DP)
       call dense_scale(fmo_masks%mxn(is),0.0_DP)
       call dense_scale(fmo_masks%mxm(is),0.0_DP)
    end do

    ! mjsp: create buffer fragment masks used to construct
    ! the supermolecule mask

    ! mjsp: allocations:
    allocate(mask_buf_nxm(pub_num_spins),stat=ierr)
    call utils_alloc_check('fmo_construct_masks','mask_buf_nxm',ierr)
    allocate(mask_buf_mxn(pub_num_spins),stat=ierr)
    call utils_alloc_check('fmo_construct_masks','mask_buf_mxn',ierr)
    allocate(mask_buf_mxm(pub_num_spins),stat=ierr)
    call utils_alloc_check('fmo_construct_masks','mask_buf_mxm',ierr)

    ! mjsp: create matrices:
    call dense_create(mask_buf_nxn, &
                      pub_frag_data(0)%ngwf_basis(1)%num,&
                      pub_frag_data(0)%ngwf_basis(1)%num)
    do is=1,pub_num_spins
       call dense_create(mask_buf_nxm(is), &
                         pub_frag_data(0)%ngwf_basis(1)%num,&
                         pub_frag_data(0)%rep%n_occ(is,PUB_1K))
       call dense_create(mask_buf_mxn(is), &
                         pub_frag_data(0)%rep%n_occ(is,PUB_1K),&
                         pub_frag_data(0)%ngwf_basis(1)%num)
       call dense_create(mask_buf_mxm(is), &
                         pub_frag_data(0)%rep%n_occ(is,PUB_1K),&
                         pub_frag_data(0)%rep%n_occ(is,PUB_1K))
    end do

    do it=1,pub_frag_iatm(0)

       ! if fragment-specific polarisation mode and
       ! we have not arrived at the fragment we wish to polarise:
       if ((pub_eda_mode == EDA_POLFRAGLOC_DEVEL) .and. &
           (it .ne. pub_frag_counter)) then

          ! omit blocks that relate to this fragment:
          cycle

       end if

       ! if fragment pair delocalisations mode
       ! and this fragment belongs to the set of fragments we wish
       ! to combine:
       if ((pub_eda_mode == EDA_CTFRAGLOC_DEVEL) &
            .and. ((it .eq. pub_frag_counter) &
            .or. (it .eq. pub_frag_counter2))) then

          ! if we have arrived at the fragment (pub_frag_counter) that
          ! we wish to join to another fragment (pub_frag_counter2)
          ! and polarise:
          if (it .eq. pub_frag_counter) then

             ! set optional second fragment ID:
             it2 = pub_frag_counter2

          ! if we have arrived at the fragment (pub_frag_counter2)
          ! that has already been combined above with another fragment
          ! (pub_frag_counter):
          else if (it .eq. pub_frag_counter2) then

             ! omit blocks that relate to this fragment:
             cycle

          end if

       else
          ! a single fragment we wish to polarise:

          ! fragment ID2 is redundant:
          it2 = 0

       end if


       ! ======= nxn ======
       ! mjsp: clear buffer
       call dense_scale(mask_buf_nxn,0.0_DP)

       ! mjsp: fragment mask
       call fmo_get_frag_blk_dens_nxn(mask_buf_nxn,&
            frag_id=it,frag_id2=it2)

       ! mjsp: supermolecule nxn (build from fragment masks)
       call dense_axpy(fmo_masks%nxn,mask_buf_nxn,+1.0_DP)

       do is=1,pub_num_spins

          ! ======= nxm ======
          ! mjsp: clear buffer
          call dense_scale(mask_buf_nxm(is),0.0_DP)

          ! mjsp: fragment mask
          call fmo_get_frag_blk_dens_nxm(mask_buf_nxm(is),is=is,&
               frag_id=it,frag_id2=it2)

          ! mjsp: supermolecule nxm (build from fragment masks)
          call dense_axpy(fmo_masks%nxm(is),mask_buf_nxm(is),+1.0_DP)

          ! ======= mxn ======
          ! mjsp: clear buffer
          call dense_scale(mask_buf_mxn(is),0.0_DP)

          ! mjsp: fragment mask
          call fmo_get_frag_blk_dens_mxn(mask_buf_mxn(is),is=is,&
               frag_id=it,frag_id2=it2)

          ! mjsp: supermolecule mxn (build from fragment masks)
          call dense_axpy(fmo_masks%mxn(is),mask_buf_mxn(is),+1.0_DP)

          ! ======= mxm ======
          ! mjsp: clear buffer
          call dense_scale(mask_buf_mxm(is),0.0_DP)

          ! mjsp: fragment mask
          call fmo_get_frag_blk_dens_mxm(mask_buf_mxm(is),is=is,&
               frag_id=it,frag_id2=it2)

          ! mjsp: supermolecule mxm (build from fragment masks)
          call dense_axpy(fmo_masks%mxm(is),mask_buf_mxm(is),+1.0_DP)

       end do ! spin

    end do

    ! ======= nxn_r ======
    if (pub_eda_nodiag .gt. 0) then
       do it=1,pub_frag_iatm(0) ! loop fragments

          ! mjsp: fragment mask (rows-only fragmentation)
          call fmo_get_frag_blk_dens_nxn(fmo_masks%nxn_r(it),&
               frag_id=it,frag_id2=it2, fragment_rows_only=.true.)

          ! mjsp: nxn_c is simply the transpose of nxn_r
          call dense_transpose(fmo_masks%nxn_c(it), fmo_masks%nxn_r(it))

          ! mjsp: Sparse masks:
          call dense_convert(fmo_masks%nxn_scK(it), fmo_masks%nxn_c(it))
          call dense_convert(fmo_masks%nxn_srK(it), fmo_masks%nxn_r(it))
          call dense_convert(fmo_masks%nxn_scHKS(it), fmo_masks%nxn_c(it))
          call dense_convert(fmo_masks%nxn_srHKS(it), fmo_masks%nxn_r(it))
          call dense_convert(fmo_masks%nxn_scSKHKS(it), fmo_masks%nxn_c(it))
          call dense_convert(fmo_masks%nxn_srSKHKS(it), fmo_masks%nxn_r(it))

       end do

       ! mjsp: Sparse nxn mask:
       call dense_convert(fmo_masks%nxn_sH, fmo_masks%nxn)
       call dense_convert(fmo_masks%nxn_sK, fmo_masks%nxn)
       call dense_convert(fmo_masks%nxn_sHKS, fmo_masks%nxn)
       call dense_convert(fmo_masks%nxn_sSKHKS, fmo_masks%nxn)

    end if

    ! cleanup
    call dense_destroy(mask_buf_nxn)
    do is=1,pub_num_spins
       call dense_destroy(mask_buf_nxm(is))
       call dense_destroy(mask_buf_mxn(is))
       call dense_destroy(mask_buf_mxm(is))
    end do
    deallocate(mask_buf_nxm,stat=ierr)
    call utils_dealloc_check('fmo_construct_masks','mask_buf_nxm',ierr)
    deallocate(mask_buf_mxn,stat=ierr)
    call utils_dealloc_check('fmo_construct_masks','mask_buf_mxn',ierr)
    deallocate(mask_buf_mxm,stat=ierr)
    call utils_dealloc_check('fmo_construct_masks','mask_buf_mxm',ierr)

    ! Stop Timer
    call timer_clock('fragment_construct_masks',2)

  end subroutine fmo_construct_masks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_destroy_masks()

    !========================================================================!
    ! This subroutine destroys the fmo_mask matrices used to extract         !
    ! the fragment matrices from a supermolecule matrix.                     !
    !------------------------------------------------------------------------!
    ! Written by Max Phipps January 2016.                                    !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    use dense, only: dense_destroy
    use rundat, only: pub_eda_nodiag, pub_frag_iatm, pub_num_spins
    use sparse_embed, only: sparse_embed_destroy
    use utils, only: utils_dealloc_check

    implicit none

    integer :: is, it, ierr

    ! mjsp: destroy masks
    call dense_destroy(fmo_masks%nxn)
    do is=pub_num_spins,1,-1
       call dense_destroy(fmo_masks%nxm(is))
       call dense_destroy(fmo_masks%mxn(is))
       call dense_destroy(fmo_masks%mxm(is))
    end do

    ! mjsp: Deallocations
    deallocate(fmo_masks%nxm,stat=ierr)
    call utils_dealloc_check('fmo_destroy_masks','fmo_masks%nxm',ierr)
    deallocate(fmo_masks%mxn,stat=ierr)
    call utils_dealloc_check('fmo_destroy_masks','fmo_masks%mxn',ierr)
    deallocate(fmo_masks%mxm,stat=ierr)
    call utils_dealloc_check('fmo_destroy_masks','fmo_masks%mxm',ierr)

    ! mjsp: Special row-only and column-only fragmented nxn mask
    if (pub_eda_nodiag .gt. 0) then

       call sparse_embed_destroy(fmo_masks%nxn_sH)
       call sparse_embed_destroy(fmo_masks%nxn_sK)
       call sparse_embed_destroy(fmo_masks%nxn_sHKS)
       call sparse_embed_destroy(fmo_masks%nxn_sSKHKS)

       do it=1,pub_frag_iatm(0) ! loop fragments

          call dense_destroy(fmo_masks%nxn_r(it))
          call dense_destroy(fmo_masks%nxn_c(it))
          call sparse_embed_destroy(fmo_masks%nxn_scK(it))
          call sparse_embed_destroy(fmo_masks%nxn_srK(it))
          call sparse_embed_destroy(fmo_masks%nxn_scHKS(it))
          call sparse_embed_destroy(fmo_masks%nxn_srHKS(it))
          call sparse_embed_destroy(fmo_masks%nxn_scSKHKS(it))
          call sparse_embed_destroy(fmo_masks%nxn_srSKHKS(it))

       end do

       deallocate(fmo_masks%nxn_r,stat=ierr)
       call utils_dealloc_check('fmo_destroy_masks','fmo_masks%nxn_r',ierr)

       deallocate(fmo_masks%nxn_c,stat=ierr)
       call utils_dealloc_check('fmo_destroy_masks','fmo_masks%nxn_c',ierr)

       deallocate(fmo_masks%nxn_scK,stat=ierr)
       call utils_dealloc_check('fmo_destroy_masks','fmo_masks%nxn_scK',ierr)

       deallocate(fmo_masks%nxn_srK,stat=ierr)
       call utils_dealloc_check('fmo_destroy_masks','fmo_masks%nxn_srK',ierr)

       deallocate(fmo_masks%nxn_scHKS,stat=ierr)
       call utils_dealloc_check('fmo_destroy_masks','fmo_masks%nxn_scHKS',ierr)

       deallocate(fmo_masks%nxn_srHKS,stat=ierr)
       call utils_dealloc_check('fmo_destroy_masks','fmo_masks%nxn_srHKS',ierr)

       deallocate(fmo_masks%nxn_scSKHKS,stat=ierr)
       call utils_dealloc_check('fmo_destroy_masks','fmo_masks%nxn_scSKHKS',ierr)

       deallocate(fmo_masks%nxn_srSKHKS,stat=ierr)
       call utils_dealloc_check('fmo_destroy_masks','fmo_masks%nxn_srSKHKS',ierr)

    end if

  end subroutine fmo_destroy_masks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_freeze_mat_spam_nxn(mat, frag_id, frag_id2)

    !========================================================================!
    ! This subroutine freezes the off diagonal blocks of a SPAM3 format      !
    ! NGWFxNGWF indexed matrix, in order to remove fragment-fragment         !
    ! interactions in the optimisation routine.                              !
    ! If frag_id=0, then all interfragment elements will be set to zero.     !
    ! If frag_id and/or frag_id2 are /= 0, then these two fragments will     !
    !      be treated as a single fragment and the interfragment elements    !
    !      of the fragments, including this 'joined' fragment, will be set   !
    !      to zero.                                                          !
    ! NOTE: This subroutine assumes the matrix is square.                    !
    !------------------------------------------------------------------------!
    ! mat         (inout) : Matrix (NGWF x NGWF indexing) to be zeroed.      !
    ! frag_id     (in)    : Fragment specific freezing.                      !
    !                       Used calculate fragment-specific polarisation    !
    !                       contributions by zeroing all but the frag_id     !
    !                       fragment matrix block.                           !
    ! frag_id2    (in)    : ID of an additional fragment block to be         !
    !                       'joined' to the frag_id fragment.                !
    !------------------------------------------------------------------------!
    ! Written by Max Phipps April 2015.                                      !
    ! Modified to clean up unused parameters by Max Phipps January 2016.     !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!


    use datatypes, only: COEF ! agrecocmplx
    use comms, only: pub_my_proc_id
    use fragment_data, only: NGWF_INDEX_MAP_TYPE, pub_frag_data, &
         fragment_data_get_ngwf_index_map
    use rundat, only: pub_frag_iatm
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_copy, sparse_embed_scale, sparse_embed_destroy
    use sparse, only: sparse_put_element, sparse_get_element
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)    :: mat
    integer, intent(in)           :: frag_id
    integer, intent(in), optional :: frag_id2

    ! Local Variables
    type(SPAM3_EMBED) :: buf_mat
    integer :: it, ngwf_offset, ngwf_offset2
    integer :: ingwf_frag_i, ingwf_frag_j, el_i, el_j
    integer :: el_i_proc, el_j_proc
    !real(kind=DP) :: element
    ! agrecocmplx: use COEF type
    type(COEF) :: element
    integer :: frag_id2_val

    type(NGWF_INDEX_MAP_TYPE) :: ngwf_index_map

    ! Start Timer
    call timer_clock('fragment_matrix_ops_frz_spam_nxn',1)

    ! --------------------------------------------------------------------------

    ngwf_index_map = fragment_data_get_ngwf_index_map()

    frag_id2_val = 0  ! redundant if frag_id2 optional argument not used
    if (present(frag_id2)) then
       if (frag_id2 .ne. 0) frag_id2_val = frag_id2
    end if

    ! jme: calls to sparse_get/put_element will preserve this component
    ! throughout the subroutine
    element%iscmplx = mat%iscmplx

    ! copy matrix to buffer
    call sparse_embed_create(buf_mat, mat)
    call sparse_embed_copy(buf_mat, mat)

    ! clear matrix
    call sparse_embed_scale(mat, 0.0_DP)

    ! copy the fragment blocks in buffer back into the (now cleared) matrix:

    ngwf_offset = 0

    do it=1,pub_frag_iatm(0) ! loop fragments

       ! iterate NGWFs on this fragment
       do ingwf_frag_i=1,pub_frag_data(it)%ngwf_basis(1)%num ! COL

          ! Find the 'i'th index in the supermolecule matrix
          ! that this element is associated with.
          ! (fragment ingwf_frag_i -> supermolecule el_i)
          el_i = ngwf_index_map%supers_by_cum_frag( &
               ingwf_frag_i + ngwf_offset)

          ! mjsp: The NGWFs are distributed column-wise over the procs.
          ! We check if el_i is on our proc, and then loop all column
          ! rows and copy over the appropriate matrix elements from the
          ! fragment calculation.
          ! find if on this proc:
          el_i_proc = pub_frag_data(0)%ngwf_basis(1)%proc_of_func(el_i)

          if (pub_my_proc_id == el_i_proc) then

             do ingwf_frag_j=1,pub_frag_data(it)%ngwf_basis(1)%num ! ROW

                ! mjsp: Find the 'j'th index in the supermolecule matrix
                ! that this element is associated with
                ! (fragment ingwf_frag_j -> supermolecule el_j)
                el_j = ngwf_index_map%supers_by_cum_frag( &
                     ingwf_frag_j + ngwf_offset)

                ! get element from buf_mat
                call sparse_get_element(element, buf_mat%p, el_j, el_i)

                ! put element from buf_mat into mat
                call sparse_put_element(element, mat%p, &
                     el_j, & ! row
                     el_i) ! col

             end do ! NGWFs on fragment

          end if ! if on proc

       end do ! NGWFs on fragment

       ! modify NGWF offset
       ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num

    end do  ! fragments


    ! mjsp: if frag_id2 is present then we should include
    ! mjsp: the off-diagonal blocks between fragments frag_id and
    ! mjsp: frag_id2, in order to 'join' fragments frag_id and frag_id2:
    if (frag_id2_val .ne. 0) then

       ! mjsp: calculate the ngwf_offset values:
       ! frag_id:
       ngwf_offset  = 0
       do it=1,pub_frag_iatm(0) ! loop fragments
          if (frag_id .ne. it) then
             ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num
             cycle
          else
             exit
          end if
       end do
       ! frag_id2:
       ngwf_offset2 = 0
       do it=1,pub_frag_iatm(0) ! loop fragments
          if (frag_id2_val .ne. it) then
             ngwf_offset2 = ngwf_offset2 + pub_frag_data(it)%ngwf_basis(1)%num
             cycle
          else
             exit
          end if
       end do

       ! iterate NGWFs on frag_id fragment
       do ingwf_frag_j=1,pub_frag_data(frag_id)%ngwf_basis(1)%num ! COL

          ! Find the 'j'th index in the supermolecule matrix
          ! that this element is associated with.
          ! (fragment ingwf_frag_j -> supermolecule el_j)
          el_j = ngwf_index_map%supers_by_cum_frag( &
               ingwf_frag_j + ngwf_offset)

          ! iterate NGWFs on frag_id2 fragment
          do ingwf_frag_i=1,pub_frag_data(frag_id2_val)%ngwf_basis(1)%num ! ROW

             ! mjsp: Find the 'i'th index in the supermolecule matrix
             ! that this element is associated with
             ! (fragment ingwf_frag_i -> supermolecule el_i)
             el_i = ngwf_index_map%supers_by_cum_frag( &
                  ingwf_frag_i + ngwf_offset2)


             ! mjsp: The NGWFs are distributed column-wise over the procs.
             ! We check if el_j is on our proc, and then loop all column
             ! rows and copy over the appropriate matrix elements from the
             ! fragment calculation.
             ! find if on this proc:
             el_j_proc = pub_frag_data(0)%ngwf_basis(1)%proc_of_func(el_j)

             if (pub_my_proc_id == el_j_proc) then

                ! get element from buf_mat
                call sparse_get_element(element, buf_mat%p, el_i, el_j)

                ! put element from buf_mat into mat
                call sparse_put_element(element, mat%p, &
                     el_i, & ! row
                     el_j) ! col

             end if ! if on proc


             ! Upper right subblock:
             el_i_proc = pub_frag_data(0)%ngwf_basis(1)%proc_of_func(el_i)

             if (pub_my_proc_id == el_i_proc) then

                ! get element from buf_mat
                call sparse_get_element(element, buf_mat%p, el_j, el_i)

                ! put element from buf_mat into mat
                call sparse_put_element(element, mat%p, &
                     el_j, & ! row
                     el_i) ! col

             end if ! if on proc


          end do ! NGWFs on fragment

       end do ! NGWFs on fragment

    end if

    ! cleanup
    call sparse_embed_destroy(buf_mat)

    ! Stop Timer
    call timer_clock('fragment_matrix_ops_frz_spam_nxn',2)


  end subroutine fmo_freeze_mat_spam_nxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_copy_frag_block_spam_nxn(mat_dst, mat_in, frag_id, &
         frag_id2)

    !========================================================================!
    ! This subroutine copies the fragment block data from one NGWFxNGWF      !
    ! indexed matrix to another.                                             !
    ! This effectively results in a duplicated matrix with elements not      !
    ! related to the fragment (given by frag_id) set to zero.                !
    !------------------------------------------------------------------------!
    ! mat_dst   (inout)   : Matrix (NGWF x NGWF indexing) to receive.        !
    ! mat_in    (inout)   : Matrix (NGWF x NGWF indexing) to be copied.      !
    ! frag_id   (in)      : ID of the fragment block to be copied.           !
    ! frag_id2  (in)      : ID of the additional fragment block to be copied.!
    !------------------------------------------------------------------------!
    ! Written by Max Phipps April 2015.                                      !
    ! Rewritten by Max Phipps December 2015.                                 !
    ! Adapted to allow complex NGWFs by JM Escartin (summer 2016).           !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    use datatypes, only: COEF
    use comms, only: pub_my_proc_id
    use fragment_data, only: NGWF_INDEX_MAP_TYPE, pub_frag_data, &
         fragment_data_get_ngwf_index_map
    use rundat, only: pub_frag_iatm
    use sparse, only: sparse_put_element, sparse_get_element
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)    :: mat_dst
    type(SPAM3_EMBED), intent(in)       :: mat_in
    integer, intent(in)           :: frag_id
    integer, intent(in), optional :: frag_id2

    ! Local Variables
    integer :: it, ngwf_offset, ngwf_offset2
    integer :: ingwf_frag_i, ingwf_frag_j, el_i, el_j
    integer :: el_i_proc, el_j_proc
    type(COEF) :: element
    integer :: frag_id2_val, num_frags_copied, num_frags_copied_target
    character(len=*), parameter :: myself = 'fmo_copy_frag_block_spam_nxn'

    type(NGWF_INDEX_MAP_TYPE) :: ngwf_index_map

    ! Start Timer
    call timer_clock('fragment_matrix_ops_spam_nxn',1)

    ! --------------------------------------------------------------------------

    ngwf_index_map = fragment_data_get_ngwf_index_map()

    num_frags_copied = 0
    frag_id2_val = 0  ! redundant if frag_id2 optional argument not used
    num_frags_copied_target = 1
    if (present(frag_id2)) then
       if (frag_id2 .ne. 0) then
          frag_id2_val = frag_id2
          num_frags_copied_target = 2
       end if
    end if

    ! jme: initialize complex character of the copy
    call utils_assert(mat_dst%iscmplx .eqv. mat_in%iscmplx, &
         'Error in '//myself//': real and complex matrices are being mixed.')
    element%iscmplx = mat_in%iscmplx

    ! copy the fragment blocks from mat_in to mat_dst

    ngwf_offset = 0

    do it=1,pub_frag_iatm(0) ! loop fragments

       ! only copy the block(s) that relates to the fragment(s):
       if (.not.((it .eq. frag_id) .or. (it .eq. frag_id2_val))) then
          ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num
          cycle
       end if

       ! iterate NGWFs on this fragment
       do ingwf_frag_j=1,pub_frag_data(it)%ngwf_basis(1)%num ! COL

          ! Find the 'j'th index in the supermolecule matrix
          ! that this element is associated with.
          ! (fragment ingwf_frag_j -> supermolecule el_j)
          el_j = ngwf_index_map%supers_by_cum_frag( &
               ingwf_frag_j + ngwf_offset)

          ! mjsp: The NGWFs are distributed column-wise over the procs.
          ! We check if el_j is on our proc, and then loop all column
          ! rows and copy over the appropriate matrix elements from the
          ! fragment calculation.
          ! find if on this proc:
          el_j_proc = pub_frag_data(0)%ngwf_basis(1)%proc_of_func(el_j)

          if (pub_my_proc_id == el_j_proc) then

             do ingwf_frag_i=1,pub_frag_data(it)%ngwf_basis(1)%num ! ROW

                ! mjsp: Find the 'j'th index in the supermolecule matrix
                ! that this element is associated with
                ! (fragment ingwf_frag_i -> supermolecule el_i)
                el_i = ngwf_index_map%supers_by_cum_frag( &
                     ingwf_frag_i + ngwf_offset)

                ! get element from mat_in
                call sparse_get_element(element, mat_in%p, el_i, el_j)

                ! put element from mat_in into mat_dst
                call sparse_put_element(element, mat_dst%p, &
                     el_i, & ! row
                     el_j) ! col

             end do ! NGWFs on fragment

          end if ! if on proc

       end do ! NGWFs on fragment

       num_frags_copied = num_frags_copied + 1

       ! modify NGWF offset
       ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num

       if (num_frags_copied .eq. num_frags_copied_target) &
          ! exit loop as we have copied this fragment's data across if we have
          ! arrived here
          exit

    end do  ! fragments


    ! mjsp: if frag_id2 is present then we should include
    ! mjsp: the off-diagonal blocks between fragments frag_id and
    ! mjsp: frag_id2, in order to 'join' fragments frag_id and frag_id2:
    if (frag_id2_val .ne. 0) then

       ! mjsp: calculate the ngwf_offset values:
       ! frag_id:
       ngwf_offset  = 0
       do it=1,pub_frag_iatm(0) ! loop fragments
          if (frag_id .ne. it) then
             ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num
             cycle
          else
             exit
          end if
       end do
       ! frag_id2:
       ngwf_offset2 = 0
       do it=1,pub_frag_iatm(0) ! loop fragments
          if (frag_id2_val .ne. it) then
             ngwf_offset2 = ngwf_offset2 + pub_frag_data(it)%ngwf_basis(1)%num
             cycle
          else
             exit
          end if
       end do

       ! iterate NGWFs on frag_id fragment
       do ingwf_frag_j=1,pub_frag_data(frag_id)%ngwf_basis(1)%num

          ! Find the 'j'th index in the supermolecule matrix
          ! that this element is associated with.
          ! (fragment ingwf_frag_j -> supermolecule el_j)
          el_j = ngwf_index_map%supers_by_cum_frag( &
               ingwf_frag_j + ngwf_offset)

          ! iterate NGWFs on frag_id2 fragment
          do ingwf_frag_i=1,pub_frag_data(frag_id2_val)%ngwf_basis(1)%num

             ! mjsp: Find the 'i'th index in the supermolecule matrix
             ! that this element is associated with
             ! (fragment ingwf_frag_i -> supermolecule el_i)
             el_i = ngwf_index_map%supers_by_cum_frag( &
                  ingwf_frag_i + ngwf_offset2)


             ! mjsp: The NGWFs are distributed column-wise over the procs.
             ! We check if el_j is on our proc, and then loop all column
             ! rows and copy over the appropriate matrix elements from the
             ! fragment calculation.
             ! find if on this proc:
             el_j_proc = pub_frag_data(0)%ngwf_basis(1)%proc_of_func(el_j)

             if (pub_my_proc_id == el_j_proc) then

                ! get element from mat_in
                call sparse_get_element(element, mat_in%p, el_i, el_j)

                ! put element from mat_in into mat_dst
                call sparse_put_element(element, mat_dst%p, &
                     el_i, & ! row
                     el_j) ! col

             end if ! if on proc


             ! Upper right subblock:
             el_i_proc = pub_frag_data(0)%ngwf_basis(1)%proc_of_func(el_i)

             if (pub_my_proc_id == el_i_proc) then

                ! get element from mat_in
                call sparse_get_element(element, mat_in%p, el_j, el_i)

                ! put element from mat_in into mat_dst
                call sparse_put_element(element, mat_dst%p, &
                     el_j, & ! row
                     el_i) ! col

             end if ! if on proc


          end do ! NGWFs on fragment

       end do ! NGWFs on fragment

    end if

    ! Stop Timer
    call timer_clock('fragment_matrix_ops_spam_nxn',2)

  end subroutine fmo_copy_frag_block_spam_nxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_invert_frag_block_dens_nxn(mat_in)

    !========================================================================!
    ! This subroutine inverts the fragment blocks of a dense matrix passed   !
    ! in.  This is equivalent to inverting the full matrix but by iterating  !
    ! over the fragment blocks we achieve this at reduced computational      !
    ! expense.  The dense matrix returned will have dimensions               !
    ! N_ngwfs(X) x N_ngwfs(X) where X is the fragment requested via frag_id. !
    !                                                                        !
    ! NOTE: If pub_frag_counter2 .ne. 0 then the two fragments               !
    ! pub_frag_counter and pub_frag_counter2 will be treated as one          !
    ! (combined) fragment.                                                   !
    ! NOTE: The input matrix is assumed to be symmetric.                     !
    !------------------------------------------------------------------------!
    ! mat_in    (inout)   : Matrix to be fragment-block inverted             !
    !                       (N(frag)xN(frag))                                !
    !------------------------------------------------------------------------!
    ! Written by Max Phipps August 2015.                                     !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    use datatypes, only: COEF ! agrecocmplx
    use dense, only: DEM, dense_put_element, dense_get_element, dense_invert, &
         dense_create, dense_destroy
    use fragment_data, only: NGWF_INDEX_MAP_TYPE, pub_frag_data, &
         fragment_data_get_ngwf_index_map
    use rundat, only: pub_frag_iatm, pub_frag_counter, pub_frag_counter2
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(DEM), intent(inout)      :: mat_in

    ! Local Variables
    integer :: it, it2, ngwf_offset, num_frag_ngwf
    integer :: ingwf_frag_i, ingwf_super_i, ingwf_frag_j, ingwf_super_j
    !real(kind=DP) :: element
    ! agrecocmplx: change to COEF type
    type(COEF) :: element
    type(DEM) :: mat_buf1  ! buffer matrix
    logical :: iscmplx

    type(NGWF_INDEX_MAP_TYPE) :: ngwf_index_map

    ! Start Timer
    call timer_clock('fragment_matrix_invert',1)

    ! --------------------------------------------------------------------------

    ngwf_index_map = fragment_data_get_ngwf_index_map()

    ngwf_offset = 0

    ! agrecocmplx
    iscmplx = mat_in%iscmplx
    ! jme
    element%iscmplx = iscmplx

    do it=1,pub_frag_iatm(0) ! loop fragments

       ! number of NGWFs on this fragment:
       num_frag_ngwf = pub_frag_data(it)%ngwf_basis(1)%num

       ! initialisations
       it2 = 0

       ! intrafragmental (charge transfer) delocalisations:
       if (pub_frag_counter2 .ne. 0) then
          if (it .eq. pub_frag_counter2) then

             ! iterate fragment NGWF offset and MO offset
             ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num

             ! omit blocks that relate to pub_frag_counter2:
             ! (these will be grouped with pub_frag_counter fragment)
             cycle

          else if (it .eq. pub_frag_counter) then

             ! join pub_frag_counter and pub_frag_counter2 fragments:

             ! set optional second fragment ID:
             it2 = pub_frag_counter2

             ! number of NGWFs on this fragment pair:
             num_frag_ngwf = pub_frag_data(it)%ngwf_basis(1)%num + &
                             pub_frag_data(it2)%ngwf_basis(1)%num

          end if
       end if

       ! mjsp: 1. Initialisation
       ! mjsp: allocate buffer matrix
       call dense_create(mat_buf1,num_frag_ngwf,num_frag_ngwf, iscmplx=iscmplx)

       ! 2. copy elements of the supermolecule matrix to our buffer matrix:
!       call timer_clock('fragment_matrix_invert',1)
!       call fmo_get_frag_blk_dens_nxn(mat_buf1, mat_in, it, it2)
!       call timer_clock('fragment_matrix_invert',2)

       ! mjsp: NOTE: This loop has not been optimised (compare to
       ! code in nxn loops)
       ! mjsp: iterate NGWFs on this fragment
       do ingwf_frag_i=1,num_frag_ngwf

          ! Find the 'i'th index in the supermolecule matrix
          ! that this element is associated with.
          ! (fragment ingwf_frag_i -> supermolecule ingwf_super_i)
          ingwf_super_i = ngwf_index_map%supers_by_cum_frag( &
               ingwf_frag_i + ngwf_offset)

          do ingwf_frag_j=1,num_frag_ngwf

             ! mjsp: Find the 'j'th index in the supermolecule matrix
             ! that this element is associated with
             ! (fragment ingwf_frag_j -> supermolecule ingwf_super_j)
             ingwf_super_j = ngwf_index_map%supers_by_cum_frag( &
                  ingwf_frag_j + ngwf_offset)

             ! get element from mat_in
             call dense_get_element(element, mat_in, &
                       ingwf_super_j, ingwf_super_i)

             ! put element from mat_in into mat_buf1
             ! (in the basis indexing found in the fragment system)
             call dense_put_element(element, mat_buf1, &
                  ingwf_frag_j, & ! row
                  ingwf_frag_i) ! col

          end do ! NGWFs on fragment

       end do ! NGWFs on fragment

       ! mjsp: iterate NGWF offset
       ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num


       ! 3. invert buffer matrix:
       call dense_invert(mat_buf1)

       ! 4. replace bufmat elements back into mat_in:
       call timer_clock('fragment_matrix_invert',2)
       call fmo_put_frag_blk_dens_nxn(mat_in, mat_buf1, it, it2, &
                 frag_loc=.true.)
       call timer_clock('fragment_matrix_invert',1)

       ! 5. cleanup
       call dense_destroy(mat_buf1)

    end do  ! fragments

    ! Stop Timer
    call timer_clock('fragment_matrix_invert',2)

  end subroutine fmo_invert_frag_block_dens_nxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_get_frag_blk_dens_nxn(mat_dst, mat_in, frag_id, frag_id2, &
                frag_loc, fragment_rows_only)

    !========================================================================!
    ! This subroutine copies the fragment block data from one dense matrix   !
    ! to another dense matrix.                                               !
    ! if frag_loc==.true. then the dense matrix returned                     !
    !   will has dimensions N_ngwfs(X) x N_ngwfs(X) where X is the fragment  !
    !   requested via frag_i.                                                !
    ! if frag_loc==.false. then the dense matrix returned                    !
    !   will has dimensions N_ngwfs(A) x N_ngwfs(A) where A is the           !
    !   supermolecule (constructed from the fragments X)                     !
    ! @SYMMETRY: Assumes symmetric matrices for speed-up (currently ON)      !
    !------------------------------------------------------------------------!
    ! mat_dst   (inout)   : Matrix to receive.                               !
    ! mat_in    (inout)   : Matrix to be copied. If mat_in is not present    !
    !                       then identity matrix will be copied to mat_dst.  !
    ! frag_id   (in)      : ID of the fragment block to be copied.           !
    ! frag_id2  (in)      : ID of the additional fragment block to be copied.!
    ! frag_loc  (in)      : Whether to copy to the fragment pattern          !
    !                       (used for avoiding rank deficiency issues).      !
    ! fragment_rows_only (in) : Whether to fragment the rows but not the     !
    !                           columns (special fragmentation)              !
    !------------------------------------------------------------------------!
    ! Written by Max Phipps June 2015.                                       !
    ! Optimised by Max Phipps, December 2015                                 !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    ! mjsp: TODO: frag_rows_only_log
    ! mjsp: fragment_only_rows not implemented for joining two fragments
    ! mjsp: together, and therefore diagonalisation-free
    ! mjsp: fragment-specific EDA analyses are not yet available.

    use dense, only: DEM, dense_copy, dense_put_element, dense_get_element, &
         dense_create, dense_destroy, dense_transpose
    use fragment_data, only: pub_frag_data
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout)        :: mat_dst
    type(DEM), intent(in), optional :: mat_in
    integer, intent(in)             :: frag_id
    integer, intent(in), optional   :: frag_id2
    logical, intent(in), optional   :: frag_loc
    logical, intent(in), optional   :: fragment_rows_only

    ! Local Variables
    integer :: frag_id2_val
    logical :: frag_loc_log
    logical :: frag_rows_only_log
    type(DEM) :: mat_buf1, mat_buf2 !, mat_buf3  ! buffer matrices
    integer :: num_in, num_dst
    real(kind=DP) :: col(pub_frag_data(0)%ngwf_basis(1)%num)
#ifndef SCALAPACK
    real(kind=DP) :: element
    integer :: ingwf_frag_i, ingwf_super_i, ingwf_frag_j, ingwf_super_j
    logical :: copy_row_log
#endif
    character(len=*), parameter :: myself = 'fmo_get_frag_blk_dens_nxn'

    ! Start Timer
    call timer_clock('fragment_matrix_ops_blk_nxn',1)

    ! --------------------------------------------------------------------------

    ! jmecmplx
    call utils_assert(.not. mat_dst%iscmplx, &
         'Error in '//myself//': not yet ready for complex NGWFs.')
    if (present(mat_in)) then
       call utils_assert(.not. mat_in%iscmplx, &
            'Error in '//myself//': not yet ready for complex NGWFs.')
    end if

    frag_id2_val = 0  ! redundant if frag_id2 optional argument not used
    if (present(frag_id2)) then
       if (frag_id2 .ne. 0) then
          frag_id2_val = frag_id2
       end if
    end if
    frag_loc_log = .false.
    if (present(frag_loc)) frag_loc_log = frag_loc
    frag_rows_only_log = .false.
    if (present(fragment_rows_only)) frag_rows_only_log = fragment_rows_only

#ifdef SCALAPACK

    ! mjsp: Optimised ScaLAPACK nxn fragment block routine
    ! mjsp: to loop columns rather than loop
    ! mjsp: element-by-element
    !
    ! 1. Get fragment columns
    ! 2. Transpose
    ! 3. Get fragment columns (as rows)
    ! (4. Transpose ) ! @SYMMETRY: slightly faster routine that
                      ! exploits matrix symmetry

    ! mjsp: Initialise:

    ! mjsp: Construct buffer matrices
    ! mjsp: (dimensions calculated because argument matrix
    ! may not be present)
    num_in = pub_frag_data(0)%ngwf_basis(1)%num
    if (frag_loc_log) then
       ! mjsp: dimensions of output matrix
       num_dst = pub_frag_data(frag_id)%ngwf_basis(1)%num
       ! mjsp: if fragment joining
       if (frag_id2_val .gt. 0) &
          num_dst = num_dst + pub_frag_data(frag_id2_val)%ngwf_basis(1)%num

       call dense_create(mat_buf1, num_in, num_dst)
       call dense_create(mat_buf2, num_dst, num_in)
!       call dense_create(mat_buf3, num_dst, num_dst) ! @SYMMETRY_OFF
    else
       ! num_dst = pub_frag_data(0)%ngwf_basis%num
       call dense_create(mat_buf1, num_in, num_in)
       call dense_create(mat_buf2, num_in, num_in)
!       call dense_create(mat_buf3, num_in, num_in) ! @SYMMETRY_OFF
    end if

    col(:) = 1.0_DP

    ! mjsp: Copying mat_in
    if (present(mat_in)) then

       ! mjsp: extract columns for frag_id
       call internal_nxn_cols1(mat_buf1, mat_in, frag_id, frag_id2_val)

    ! mjsp: Else identity
    else

       ! mjsp: extract columns for frag_id
       call internal_nxn_cols1(dem_dst=mat_buf1,frag_iter=frag_id,&
            frag_iter2=frag_id2_val)

    end if

    ! mjsp: For special row-only fragmentation mask, then only
    ! mjsp: fragment the rows and not the columns.
    if (frag_rows_only_log) then

       ! mjsp: Copy matrix (mat_buf1 -> mat_dst)
       call dense_copy(mat_dst,mat_buf1)

    ! mjsp: For everything else, also fragment the rows
    else

       ! mjsp: Transpose matrix (mat_buf1 -> mat_buf2)
       call dense_transpose(mat_buf2,mat_buf1)


       ! mjsp: Extract the columns again
        call internal_nxn_cols2(mat_dst, mat_buf2, frag_id, &
             frag_id2_val) ! @SYMMETRY_ON
!       call internal_nxn_cols2(mat_buf3, mat_buf2, frag_id, &
!            frag_id2_val) ! @SYMMETRY_OFF


!       call dense_transpose(mat_dst,mat_buf3) ! @SYMMETRY_OFF


       ! cleanup
       call dense_destroy(mat_buf1)
       call dense_destroy(mat_buf2)
!       call dense_destroy(mat_buf3) ! @SYMMETRY_OFF

    end if

#else

    ! mjsp: Initialise
    ingwf_frag_i = 1
    ingwf_frag_j = 1
    element = 1.0_DP  ! effectively identity matrix
    copy_row_log = .true.

    ! mjsp: Only copying one fragment
    if (frag_id2_val .eq. 0) then

       ! mjsp: frag_id x frag_id block:


       ! mjsp: For special row-only fragmentation mask, then only
       ! mjsp: fragment the rows and not the columns.
       if (frag_rows_only_log) then

          ! mjsp: NOTE: symmetry not used.

          ! mjsp: iterate NGWFs on the supermolecule
          do ingwf_super_j=1,pub_frag_data(0)%ngwf_basis(1)%num ! COL

             ! mjsp: if frag_id NGWF is on supermolecule
             if (pub_super2frag_ngwf_log(ingwf_super_j,frag_id)) then

                ! mjsp: iterate all NGWFs on the supermolecule (all rows copied)
                do ingwf_super_i=1,pub_frag_data(0)%ngwf_basis(1)%num ! ROW

                   ! mjsp: if mat_in not present, then
                   ! copy identity matrix (element=1.0_DP)
                   if (present(mat_in)) &
                      ! mjsp: get element from mat_in
                      call dense_get_element(element, mat_in, &
                           ingwf_super_i, & ! row
                           ingwf_super_j) ! col

                   ! supermolecule dimensions to fragment dimensions:
                   if (frag_loc_log) then

                      ! mjsp: put element from mat_in into mat_dst
                      call dense_put_element(element, mat_dst, &
                           ingwf_frag_i, & ! row
                           ingwf_frag_j) ! col

                   ! supermolecule dimensions to supermolecule dimensions:
                   else
                      ! mjsp: put element from mat_in into mat_dst
                      call dense_put_element(element, mat_dst, &
                           ingwf_super_i, & ! row
                           ingwf_super_j) ! col
                   end if

                   ! mjsp: Iterate fragment index counter
                   ingwf_frag_i = ingwf_frag_i + 1

                end do ! i++

                ! mjsp: Iterate fragment index counter
                ingwf_frag_j = ingwf_frag_j + 1
                ingwf_frag_i = 1

             end if
          end do ! j++


       ! mjsp: For (.not. frag_rows_only_log),
       ! mjsp: also fragment the rows:
       else


          ! mjsp: iterate NGWFs on the supermolecule
          do ingwf_super_j=1,pub_frag_data(0)%ngwf_basis(1)%num ! COL

             ! mjsp: if frag_id NGWF is on supermolecule
             if (pub_super2frag_ngwf_log(ingwf_super_j,frag_id)) then

                ! mjsp: iterate NGWFs on the supermolecule
                ! mjsp: NOTE-symmetry
                do ingwf_super_i=ingwf_super_j,pub_frag_data(0)%ngwf_basis(1)%num ! ROW

                   ! mjsp: if frag_id NGWF is on supermolecule
                   if (pub_super2frag_ngwf_log(ingwf_super_i,frag_id)) then

                      ! mjsp: if mat_in not present, then
                      ! copy identity matrix (element=1.0_DP)
                      if (present(mat_in)) &
                         ! mjsp: get element from mat_in
                         call dense_get_element(element, mat_in, &
                              ingwf_super_i, & ! row
                              ingwf_super_j) ! col

                      ! supermolecule dimensions to fragment dimensions:
                      if (frag_loc_log) then

                         ! mjsp: put element from mat_in into mat_dst
                         call dense_put_element(element, mat_dst, &
                              ingwf_frag_i, & ! row
                              ingwf_frag_j) ! col
                         call dense_put_element(element, mat_dst, &
                              ingwf_frag_j, & ! row
                              ingwf_frag_i) ! col

                      ! supermolecule dimensions to supermolecule dimensions:
                      else
                         ! mjsp: put element from mat_in into mat_dst
                         call dense_put_element(element, mat_dst, &
                              ingwf_super_i, & ! row
                              ingwf_super_j) ! col
                         call dense_put_element(element, mat_dst, &
                              ingwf_super_j, & ! row
                              ingwf_super_i) ! col
                      end if

                      ! mjsp: Iterate fragment index counter
                      ingwf_frag_i = ingwf_frag_i + 1

                   end if
                end do ! i++

                ! mjsp: Iterate fragment index counter
                ingwf_frag_j = ingwf_frag_j + 1
                ingwf_frag_i = ingwf_frag_j ! symmetry

             end if
          end do ! j++

       end if ! frag_rows_only_log

    ! mjsp: If joining two fragments together:
    else

       ! mjsp: frag_id x frag_id block +
       ! mjsp: frag_id x frag_id2 block +
       ! mjsp: frag_id2 x frag_id block +
       ! mjsp: frag_id2 x frag_id2 block

       ! mjsp: iterate NGWFs on the supermolecule
       do ingwf_super_j=1,pub_frag_data(0)%ngwf_basis(1)%num ! COL

          ! mjsp: if frag_id or frag_id2 NGWF is on supermolecule
          if ((pub_super2frag_ngwf_log(ingwf_super_j,frag_id)) .or. &
              (pub_super2frag_ngwf_log(ingwf_super_j,frag_id2))) then

             ! mjsp: iterate NGWFs on the supermolecule
             ! mjsp: NOTE-symmetry
             do ingwf_super_i=ingwf_super_j,pub_frag_data(0)%ngwf_basis(1)%num ! ROW

                ! mjsp: if frag_id or frag_id2 NGWF is on supermolecule
                if ((pub_super2frag_ngwf_log(ingwf_super_i,frag_id)) .or. &
                    (pub_super2frag_ngwf_log(ingwf_super_i,frag_id2))) then

                   ! mjsp: if mat_in not present, then
                   ! copy identity matrix (element=1.0_DP)
                   if (present(mat_in)) &
                      ! mjsp: get element from mat_in
                      call dense_get_element(element, mat_in, &
                           ingwf_super_i, & ! row
                           ingwf_super_j) ! col

                   ! supermolecule dimensions to fragment dimensions:
                   if (frag_loc_log) then

                       ! mjsp: put element from mat_in into mat_dst
                       call dense_put_element(element, mat_dst, &
                            ingwf_frag_i, & ! row
                            ingwf_frag_j) ! col
                       call dense_put_element(element, mat_dst, &
                            ingwf_frag_j, & ! row
                            ingwf_frag_i) ! col

                   ! supermolecule dimensions to supermolecule dimensions:
                   else
                       call dense_put_element(element, mat_dst, &
                            ingwf_super_i, & ! row
                            ingwf_super_j) ! col
                       call dense_put_element(element, mat_dst, &
                            ingwf_super_j, & ! row
                            ingwf_super_i) ! col
                   end if

                   ! mjsp: Iterate fragment index counter
                   ingwf_frag_i = ingwf_frag_i + 1

                end if
             end do ! i++

             ! mjsp: Iterate fragment index counter
             ingwf_frag_j = ingwf_frag_j + 1
             ingwf_frag_i = ingwf_frag_j ! symmetry

          end if
       end do ! j++

    end if

#endif

    ! --------------------------------------------------------------------------

    ! Stop Timer
    call timer_clock('fragment_matrix_ops_blk_nxn',2)

#ifdef SCALAPACK

  contains

    subroutine internal_nxn_cols1(dem_dst, dem_in, frag_iter, frag_iter2)

       ! mjsp: Extracts the fragment columns from the
       ! supermolecule matrix.

       use dense, only: dense_put_col, dense_get_col

       implicit none

       ! Arguments
       type(DEM), intent(inout)         :: dem_dst
       type(DEM), intent(in), optional  :: dem_in
       integer, intent(in)              :: frag_iter
       integer, intent(in)              :: frag_iter2

       ! Local Variables
       integer :: ingwf_frag_j, ingwf_super_j
       real(kind=DP) :: col(pub_frag_data(0)%ngwf_basis(1)%num)
       integer :: frag_iter2_val

       col(:) = 1.0_DP

       ! mjsp: iterate NGWFs on the supermolecule
       ingwf_frag_j = 1
       do ingwf_super_j=1,pub_frag_data(0)%ngwf_basis(1)%num ! COL

          ! mjsp: if ingwf_super_j is on frag_iter fragment:
          if (pub_super2frag_ngwf_log(ingwf_super_j,frag_iter)) then

             ! mjsp: if dem_in not present, then
             ! copy identity matrix (col(:)=1.0_DP)
             if (present(dem_in)) &
               ! mjsp: get column at ingwf_super_j:
               call dense_get_col(col,dem_in,ingwf_super_j)

             ! mjsp: copy to dem_dst
             ! supermolecule dimensions to fragment dimensions:
             if (frag_loc_log) then
               call dense_put_col(col,dem_dst, &
                       ingwf_frag_j) ! col
             ! supermolecule dimensions to supermolecule dimensions:
             else
               call dense_put_col(col,dem_dst, &
                       ingwf_super_j) ! col
             end if

             ! mjsp: Iterate fragment index counter
             ingwf_frag_j = ingwf_frag_j + 1

          end if

          ! mjsp: If two fragments:
          if (frag_iter2 .gt. 0) then

             ! mjsp: if ingwf_super_j is on frag_iter2 fragment:
             if (pub_super2frag_ngwf_log(ingwf_super_j,frag_iter2)) then

                ! mjsp: if dem_in not present, then
                ! copy identity matrix (col(:)=1.0_DP)
                if (present(dem_in)) &
                  ! mjsp: get column at ingwf_super_j:
                  call dense_get_col(col,dem_in,ingwf_super_j)

                ! mjsp: copy to dem_dst
                ! supermolecule dimensions to fragment dimensions:
                if (frag_loc_log) then
                  call dense_put_col(col,dem_dst, &
                          ingwf_frag_j) ! col
                ! supermolecule dimensions to supermolecule dimensions:
                else
                  call dense_put_col(col,dem_dst, &
                          ingwf_super_j) ! col
                end if

                ! mjsp: Iterate fragment index counter
                ingwf_frag_j = ingwf_frag_j + 1

             end if

          end if

       end do ! j++

    end subroutine internal_nxn_cols1

    subroutine internal_nxn_cols2(dem_dst, dem_in, frag_iter, frag_iter2)

       ! mjsp: Extracts the fragment columns from the
       ! fragment buffer matrix.

       use dense, only: dense_put_col, dense_get_col

       implicit none

       ! Arguments
       type(DEM), intent(inout)  :: dem_dst
       type(DEM), intent(in)     :: dem_in
       integer, intent(in)       :: frag_iter
       integer, intent(in)       :: frag_iter2

       ! Local Variables
       integer :: ingwf_frag_j, ingwf_super_j
       real(kind=DP) :: col(pub_frag_data(0)%ngwf_basis(1)%num)

       ! mjsp: iterate NGWFs on the supermolecule
       ingwf_frag_j = 1
       do ingwf_super_j=1,pub_frag_data(0)%ngwf_basis(1)%num ! COL

          ! mjsp: if ingwf_super_j is on frag_iter fragment:
          if (pub_super2frag_ngwf_log(ingwf_super_j,frag_iter)) then

              ! supermolecule dimensions to fragment dimensions:
              if (frag_loc_log) then

                ! mjsp: get column at ingwf_super_j:
                call dense_get_col(col(1:num_dst),&
                     dem_in,ingwf_super_j)

                call dense_put_col(col(1:num_dst),&
                     dem_dst,&
                     ingwf_frag_j) ! col

              ! supermolecule dimensions to supermolecule dimensions:
              else

                ! mjsp: get column at ingwf_super_j:
                call dense_get_col(col(:),&
                     dem_in,ingwf_super_j)

                call dense_put_col(col(:),&
                     dem_dst,&
                     ingwf_super_j) ! col

              end if

             ! mjsp: Iterate fragment index counter
             ingwf_frag_j = ingwf_frag_j + 1

          end if

          ! mjsp: If two fragments:
          if (frag_iter2 .gt. 0) then

             ! mjsp: if ingwf_super_j is on frag_iter2 fragment:
             if (pub_super2frag_ngwf_log(ingwf_super_j,frag_iter2)) then

                 ! supermolecule dimensions to fragment dimensions:
                 if (frag_loc_log) then

                   ! mjsp: get column at ingwf_super_j:
                   call dense_get_col(col(1:num_dst),&
                        dem_in,ingwf_super_j)

                   call dense_put_col(col(1:num_dst),&
                        dem_dst,&
                        ingwf_frag_j) ! col

                 ! supermolecule dimensions to supermolecule dimensions:
                 else

                   ! mjsp: get column at ingwf_super_j:
                   call dense_get_col(col(:),&
                        dem_in,ingwf_super_j)

                   call dense_put_col(col(:),&
                        dem_dst,&
                        ingwf_super_j) ! col

                 end if

                ! mjsp: Iterate fragment index counter
                ingwf_frag_j = ingwf_frag_j + 1

             end if

          end if

       end do ! j++

    end subroutine internal_nxn_cols2

#endif


  end subroutine fmo_get_frag_blk_dens_nxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_get_frag_blk_dens_mxn(mat_dst, mat_in, is, frag_id, frag_id2)

    !========================================================================!
    ! This subroutine copies the fragment block data from one dense matrix   !
    ! to another dense matrix, where the dimensions/indexing is              !
    !     N_occMO(A) x N_ngwfs(A)                                            !
    ! where A is the supermolecule (i.e. elements of the supermolecule       !
    ! matrix not belonging to the fragment frag_id are set to zero)          !
    !------------------------------------------------------------------------!
    ! mat_dst   (inout)   : Matrix to receive. (M(frag)xN(frag))             !
    ! mat_in    (inout)   : Matrix to be copied. If mat_in is not present    !
    !                       then identity matrix will be copied to mat_dst.  !
    ! is        (in)      : The spin of the matrices passed.                 !
    ! frag_id   (in)      : ID of the fragment block to be copied.           !
    ! frag_id2  (in)      : ID of the additional fragment block to be copied.!
    !------------------------------------------------------------------------!
    ! Written by Max Phipps June 2015.                                       !
    ! Optimised by Max Phipps, December 2015                                 !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    use datatypes, only: COEF, data_coef_assign
    use dense, only: DEM, dense_put_element, dense_get_element
    use fragment_data, only: pub_frag_data
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout)        :: mat_dst
    type(DEM), intent(in), optional :: mat_in
    integer, intent(in)             :: is
    integer, intent(in)             :: frag_id
    integer, intent(in), optional   :: frag_id2

    ! Local Variables
    integer :: it
    type(COEF) :: element
    integer :: frag_id2_val
    integer :: imo_frag_i, imo_super_i
    integer :: ingwf_frag_j, ingwf_super_j
    integer :: mo_start, mo_end
    integer :: join_skip_mo, it_i
    character(len=*), parameter :: myself = 'fmo_get_frag_blk_dens_mxn'

    ! Start Timer
    call timer_clock('fragment_matrix_ops_blk_mxn',1)

    ! --------------------------------------------------------------------------

    ! jme: initialize complex character of the copy
    if (present(mat_in)) then
       call utils_assert(mat_dst%iscmplx .eqv. mat_in%iscmplx, &
            'Error in '//myself//': real and complex matrices are being mixed.')
    end if
    element%iscmplx = mat_dst%iscmplx

    frag_id2_val = 0  ! redundant if frag_id2 optional argument not used
    if (present(frag_id2)) then
       if (frag_id2 .ne. 0) then
          frag_id2_val = frag_id2
       end if
    end if

    ! mjsp: Initialise
    call data_coef_assign(element, 1.0_DP)   ! effectively identity matrix

    ! mjsp: Only copying one fragment
    if (frag_id2_val .eq. 0) then

       ! mjsp: frag_id x frag_id block:

       ! mjsp: Initialise
       imo_frag_i = 1
       ingwf_frag_j = 1
       mo_start = 1
       mo_end = 1

       ! mjsp: Calculate mo_start (index that this
       ! fragment's MO range begins at)
       it = 1
       do while (.not.(it .eq. frag_id))
          mo_start = mo_start + pub_frag_data(it)%rep%n_occ(is,PUB_1K)
          it = it + 1
       end do

       ! mjsp: Calculate mo_end
       mo_end = mo_start-1+pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)

       ! mjsp: iterate NGWFs on the supermolecule
       do ingwf_super_j=1,pub_frag_data(0)%ngwf_basis(1)%num ! COL

          ! mjsp: if frag_id NGWF is on supermolecule
          if (pub_super2frag_ngwf_log(ingwf_super_j,frag_id)) then

             ! mjsp: iterate occupied MOs on the fragment
             do imo_super_i=mo_start,mo_end ! ROW

                ! mjsp: if mat_in not present, then
                ! copy identity matrix (element=1.0_DP)
                if (present(mat_in)) &
                   ! mjsp: get element from mat_in
                   call dense_get_element(element, mat_in, &
                        imo_super_i, & ! row
                        ingwf_super_j) ! col

                ! mjsp: put element from mat_in into mat_dst
                call dense_put_element(element, mat_dst, &
!                     imo_frag_i, & ! row
!                     ingwf_frag_j) ! col
                     imo_super_i, & ! row
                     ingwf_super_j) ! col

!                ! mjsp: Iterate fragment index counter
!                imo_frag_i = imo_frag_i + 1

             end do

!             ! mjsp: Iterate fragment index counter
!             ingwf_frag_j = ingwf_frag_j + 1
!             imo_frag_i = 1

          end if

       end do

    ! mjsp: If joining two fragments together:
    else

       ! mjsp: frag_id x frag_id block +
       ! mjsp: frag_id x frag_id2 block +
       ! mjsp: frag_id2 x frag_id block +
       ! mjsp: frag_id2 x frag_id2 block

       ! mjsp: Initialise
       imo_frag_i = 1
       ingwf_frag_j = 1
       mo_start = 1
       mo_end = 1

       ! mjsp: Calculate mo_start (index that this
       ! fragment's MO range begins at)
       it = 1
       do while (.not.(it .eq. frag_id))
          mo_start = mo_start + pub_frag_data(it)%rep%n_occ(is,PUB_1K)
          it = it + 1
       end do

       ! mjsp: For joining fragments together, calculate mo_skip
       ! mjsp: e.g. if joining fragment 1 and 3, skip fragment 2
       join_skip_mo = 0
       do it_i=frag_id+1,frag_id2-1 ! loop fragment(s) that are between the fragment pair
          join_skip_mo = join_skip_mo + pub_frag_data(it_i)%rep%n_occ(is,PUB_1K)
       end do

       ! mjsp: Calculate mo_end
       mo_end = mo_start-1+pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)+ &
                           pub_frag_data(frag_id2)%rep%n_occ(is,PUB_1K)

       ! mjsp: iterate NGWFs on the supermolecule
       do ingwf_super_j=1,pub_frag_data(0)%ngwf_basis(1)%num ! COL

          ! mjsp: if frag_id or frag_id2 NGWF is on supermolecule
          if ((pub_super2frag_ngwf_log(ingwf_super_j,frag_id)) .or. &
              (pub_super2frag_ngwf_log(ingwf_super_j,frag_id2))) then

             ! mjsp: iterate occupied MOs on the fragments
             imo_super_i = mo_start
             do while (imo_super_i .le. (mo_end+join_skip_mo)) ! COL

                ! mjsp: if mat_in not present, then
                ! copy identity matrix (element=1.0_DP)
                if (present(mat_in)) &
                   ! mjsp: get element from mat_in
                   call dense_get_element(element, mat_in, &
                        imo_super_i, & ! row
                        ingwf_super_j) ! col

                ! mjsp: put element from mat_in into mat_dst
                call dense_put_element(element, mat_dst, &
!                     imo_frag_i, & ! row
!                     ingwf_frag_j) ! col
                     imo_super_i, & ! row
                     ingwf_super_j) ! col

!                ! mjsp: Iterate fragment index counter
!                imo_frag_i = imo_frag_i + 1

                ! mjsp: skip MOs if neccesary
                if ((imo_super_i-mo_start+1) .eq. &
                     pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)) then
                   imo_super_i = imo_super_i + 1 + join_skip_mo
                else
                   imo_super_i = imo_super_i + 1
                end if

             end do

!             ! mjsp: Iterate fragment index counter
!             ingwf_frag_j = ingwf_frag_j + 1
!             imo_frag_i = 1

          end if

       end do

    end if

    ! Stop Timer
    call timer_clock('fragment_matrix_ops_blk_mxn',2)

  end subroutine fmo_get_frag_blk_dens_mxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_get_frag_blk_dens_nxm(mat_dst, mat_in, is, frag_id, frag_id2)

    !========================================================================!
    ! This subroutine copies the fragment block data from one dense matrix   !
    ! to another dense matrix, where the dimensions/indexing is              !
    !     N_ngwfs(A) x N_occMO(A)                                            !
    ! where A is the supermolecule (i.e. elements of the supermolecule       !
    ! matrix not belonging to the fragment frag_id are set to zero)          !
    !------------------------------------------------------------------------!
    ! mat_dst   (inout)   : Matrix to receive. (N(frag)xM(frag))             !
    ! mat_in    (inout)   : Matrix to be copied. If mat_in is not present    !
    !                       then identity matrix will be copied to mat_dst.  !
    ! is        (in)      : The spin of the matrices passed.                 !
    ! frag_id   (in)      : ID of the fragment block to be copied.           !
    ! frag_id2  (in)      : ID of the additional fragment block to be copied.!
    !------------------------------------------------------------------------!
    ! Written by Max Phipps June 2015.                                       !
    ! Optimised by Max Phipps, December 2015                                 !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    use datatypes, only: COEF, data_coef_assign
    use dense, only: DEM, dense_get_element, dense_put_element
    use fragment_data, only: pub_frag_data
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout)        :: mat_dst
    type(DEM), intent(in), optional :: mat_in
    integer, intent(in)             :: is
    integer, intent(in)             :: frag_id
    integer, intent(in), optional   :: frag_id2

    ! Local Variables
    integer :: it, it_i
    type(COEF) :: element
    integer :: frag_id2_val
    integer :: join_skip_mo
    integer :: ingwf_frag_i, ingwf_super_i
    integer :: imo_frag_j, imo_super_j
    integer :: mo_start, mo_end
    character(len=*), parameter :: myself = 'fmo_get_frag_blk_dens_nxm'

    ! Start Timer
    call timer_clock('fragment_matrix_ops_blk_nxm',1)

    ! --------------------------------------------------------------------------

    ! jme: initialize complex character of the copy
    if (present(mat_in)) then
       call utils_assert(mat_dst%iscmplx .eqv. mat_in%iscmplx, &
            'Error in '//myself//': real and complex matrices are being mixed.')
    end if
    element%iscmplx = mat_dst%iscmplx

    frag_id2_val = 0  ! redundant if frag_id2 optional argument not used
    if (present(frag_id2)) then
       if (frag_id2 .ne. 0) then
          frag_id2_val = frag_id2
       end if
    end if

    ! mjsp: Initialise
    call data_coef_assign(element, 1.0_DP)   ! effectively identity matrix

    ! mjsp: Only copying one fragment
    if (frag_id2_val .eq. 0) then

       ! mjsp: frag_id x frag_id block:

       ! mjsp: Initialise
       ingwf_frag_i = 1
       imo_frag_j = 1
       mo_start = 1
       mo_end = 1

       ! mjsp: Calculate mo_start (index that this
       ! fragment's MO range begins at)
       it = 1
       do while (.not.(it .eq. frag_id))
          mo_start = mo_start + pub_frag_data(it)%rep%n_occ(is,PUB_1K)
          it = it + 1
       end do

       ! mjsp: Calculate mo_end
       mo_end = mo_start-1+pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)

       ! mjsp: iterate occupied MOs on the fragment
       do imo_super_j=mo_start,mo_end ! COL

          ! mjsp: iterate NGWFs on the supermolecule
          do ingwf_super_i=1,pub_frag_data(0)%ngwf_basis(1)%num ! ROW

             ! mjsp: if frag_id NGWF is on supermolecule
             if (pub_super2frag_ngwf_log(ingwf_super_i,frag_id)) then

                ! mjsp: if mat_in not present, then
                ! copy identity matrix (element=1.0_DP)
                if (present(mat_in)) &
                  ! mjsp: get element from mat_in
                  call dense_get_element(element, mat_in, &
                       ingwf_super_i, & ! row
                       imo_super_j) ! col

                ! mjsp: put element from mat_in into mat_dst
                call dense_put_element(element, mat_dst, &
!                     ingwf_frag_i, & ! row
!                     imo_frag_j) ! col
                     ingwf_super_i, & ! row
                     imo_super_j) ! col

!                ! mjsp: Iterate fragment index counter
!                ingwf_frag_i = ingwf_frag_i + 1

             end if

          end do

!          ! mjsp: Iterate fragment index counter
!          imo_frag_j = imo_frag_j + 1
!          ingwf_frag_i = 1

       end do

    ! mjsp: If joining two fragments together:
    else

       ! mjsp: frag_id x frag_id block +
       ! mjsp: frag_id x frag_id2 block +
       ! mjsp: frag_id2 x frag_id block +
       ! mjsp: frag_id2 x frag_id2 block

       ! mjsp: Initialise
       ingwf_frag_i = 1
       imo_frag_j = 1
       mo_start = 1
       mo_end = 1

       ! mjsp: Calculate mo_start (index that this
       ! fragment's MO range begins at)
       it = 1
       do while (.not.(it .eq. frag_id))
          mo_start = mo_start + pub_frag_data(it)%rep%n_occ(is,PUB_1K)
          it = it + 1
       end do

       ! mjsp: For joining fragments together, calculate mo_skip
       ! mjsp: e.g. if joining fragment 1 and 3, skip fragment 2
       join_skip_mo = 0
       do it_i=frag_id+1,frag_id2-1 ! loop fragment(s) that are between the fragment pair
          join_skip_mo = join_skip_mo + pub_frag_data(it_i)%rep%n_occ(is,PUB_1K)
       end do

       ! mjsp: Calculate mo_end
       mo_end = mo_start-1+pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)+ &
                           pub_frag_data(frag_id2)%rep%n_occ(is,PUB_1K)


       ! mjsp: iterate occupied MOs on the fragments
       imo_super_j = mo_start
       do while (imo_super_j .le. (mo_end+join_skip_mo)) ! COL

          ! mjsp: iterate NGWFs on the supermolecule
          do ingwf_super_i=1,pub_frag_data(0)%ngwf_basis(1)%num ! ROW

             ! mjsp: if frag_id or frag_id2 NGWF is on supermolecule
             if ((pub_super2frag_ngwf_log(ingwf_super_i,frag_id)) .or. &
                 (pub_super2frag_ngwf_log(ingwf_super_i,frag_id2))) then

                ! mjsp: if mat_in not present, then
                ! copy identity matrix (element=1.0_DP)
                if (present(mat_in)) &
                   ! mjsp: get element from mat_in
                   call dense_get_element(element, mat_in, &
                        ingwf_super_i, & ! row
                        imo_super_j) ! col

                ! mjsp: put element from mat_in into mat_dst
                call dense_put_element(element, mat_dst, &
!                     ingwf_frag_i, & ! row
!                     imo_frag_j) ! col
                     ingwf_super_i, & ! row
                     imo_super_j) ! col

!                ! mjsp: Iterate fragment index counter
!                ingwf_frag_i = ingwf_frag_i + 1

             end if

          end do

!          ! mjsp: Iterate fragment index counter
!          imo_frag_j = imo_frag_j + 1
!          ingwf_frag_i = 1

          ! mjsp: skip MOs if neccesary
          if ((imo_super_j-mo_start+1) .eq. &
               pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)) then
             imo_super_j = imo_super_j + 1 + join_skip_mo
          else
             imo_super_j = imo_super_j + 1
          end if

       end do

    end if

    ! Stop Timer
    call timer_clock('fragment_matrix_ops_blk_nxm',2)

  end subroutine fmo_get_frag_blk_dens_nxm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_get_frag_blk_dens_mxm(mat_dst, mat_in, is, frag_id, frag_id2)

    !========================================================================!
    ! This subroutine copies the fragment block data from one dense matrix   !
    ! to another dense matrix, where the dimensions/indexing is              !
    !     N_occMO(A) x N_occMO(A)                                            !
    ! where A is the supermolecule (i.e. elements of the supermolecule       !
    ! matrix not belonging to the fragment frag_id are set to zero)          !
    !------------------------------------------------------------------------!
    ! mat_dst   (inout)   : Matrix to receive. (M(frag)xM(frag))             !
    ! mat_in    (inout)   : Matrix to be copied. If mat_in is not present    !
    !                       then identity matrix will be copied to mat_dst.  !
    ! is        (in)      : The spin of the matrices passed.                 !
    ! frag_id   (in)      : ID of the fragment block to be copied.           !
    ! frag_id2  (in)      : ID of the additional fragment block to be copied.!
    !------------------------------------------------------------------------!
    ! Written by Max Phipps June 2015.                                       !
    ! Optimised by Max Phipps, December 2015                                 !
    !========================================================================!

    use datatypes, only: COEF, data_coef_assign
    use dense, only: DEM, dense_get_element, dense_put_element
    use fragment_data, only: pub_frag_data
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout)        :: mat_dst
    type(DEM), intent(in), optional :: mat_in
    integer, intent(in)             :: is
    integer, intent(in)             :: frag_id
    integer, intent(in), optional   :: frag_id2

    ! Local Variables
    integer :: it, it_i
    integer :: frag_id2_val
    integer :: join_skip_mo
    type(COEF) :: element
    integer :: imo_frag_i, imo_super_i
    integer :: imo_frag_j, imo_super_j
    integer :: mo_start, mo_end
    character(len=*), parameter :: myself = 'fmo_get_frag_blk_dens_mxm'

    ! Start Timer
    call timer_clock('fragment_matrix_ops_blk_mxm',1)

    ! --------------------------------------------------------------------------

    ! jme: initialize complex character of the copy
    if (present(mat_in)) then
       call utils_assert(mat_dst%iscmplx .eqv. mat_in%iscmplx, &
            'Error in '//myself//': real and complex matrices are being mixed.')
    end if
    element%iscmplx = mat_dst%iscmplx


    frag_id2_val = 0  ! redundant if frag_id2 optional argument not used
    if (present(frag_id2)) then
       if (frag_id2 .ne. 0) then
          frag_id2_val = frag_id2
       end if
    end if

    ! mjsp: Initialise
    call data_coef_assign(element, 1.0_DP) ! effectively identity matrix

    ! mjsp: Only copying one fragment
    if (frag_id2_val .eq. 0) then

       ! mjsp: frag_id x frag_id block:

       ! mjsp: Initialise
       imo_frag_i = 1
       imo_frag_j = 1
       mo_start = 1
       mo_end = 1

       ! mjsp: Calculate mo_start (index that this
       ! fragment's MO range begins at)
       it = 1
       do while (.not.(it .eq. frag_id))
          mo_start = mo_start + pub_frag_data(it)%rep%n_occ(is,PUB_1K)
          it = it + 1
       end do

       ! mjsp: Calculate mo_end
       mo_end = mo_start-1+pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)

       ! mjsp: iterate occupied MOs on the fragment
       do imo_super_j=mo_start,mo_end ! COL

          ! mjsp: iterate occupied MOs on the fragment
          do imo_super_i=mo_start,mo_end ! ROW

             ! mjsp: if mat_in not present, then
             ! copy identity matrix (element=1.0_DP)
             if (present(mat_in)) &
               ! mjsp: get element from mat_in
               call dense_get_element(element, mat_in, &
                    imo_super_i, & ! row
                    imo_super_j) ! col

             ! mjsp: put element from mat_in into mat_dst
             ! mjsp: NOTE-symmetry
             call dense_put_element(element, mat_dst, &
!                  imo_frag_i, & ! row
!                  imo_frag_j) ! col
                  imo_super_i, & ! row
                  imo_super_j) ! col

!             ! mjsp: Iterate fragment index counter
!             imo_frag_i = imo_frag_i + 1

          end do ! i++

!          ! mjsp: Iterate fragment index counter
!          imo_frag_j = imo_frag_j + 1
!          imo_frag_i = 1

       end do ! j++

    ! mjsp: If joining two fragments together:
    else

       ! mjsp: frag_id x frag_id block +
       ! mjsp: frag_id x frag_id2 block +
       ! mjsp: frag_id2 x frag_id block +
       ! mjsp: frag_id2 x frag_id2 block

       ! mjsp: Initialise
       imo_frag_i = 1
       imo_frag_j = 1
       mo_start = 1
       mo_end = 1

       ! mjsp: Calculate mo_start (index that this
       ! fragment's MO range begins at)
       it = 1
       do while (.not.(it .eq. frag_id))
          mo_start = mo_start + pub_frag_data(it)%rep%n_occ(is,PUB_1K)
          it = it + 1
       end do

       ! mjsp: For joining fragments together, calculate mo_skip
       ! mjsp: e.g. if joining fragment 1 and 3, skip fragment 2
       join_skip_mo = 0
       do it_i=frag_id+1,frag_id2-1 ! loop fragment(s) that are between the fragment pair
          join_skip_mo = join_skip_mo + pub_frag_data(it_i)%rep%n_occ(is,PUB_1K)
       end do

       ! mjsp: Calculate mo_end
       mo_end = mo_start-1+pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)+ &
                               pub_frag_data(frag_id2)%rep%n_occ(is,PUB_1K)

       ! mjsp: iterate occupied MOs on the fragment
       imo_super_j = mo_start
       do while (imo_super_j .le. (mo_end+join_skip_mo)) ! COL

          ! mjsp: iterate occupied MOs on the fragment
          imo_super_i = mo_start
          do while (imo_super_i .le. (mo_end+join_skip_mo)) ! ROW

             ! mjsp: if mat_in not present, then
             ! copy identity matrix (element=1.0_DP)
             if (present(mat_in)) &
               ! mjsp: get element from mat_in
               call dense_get_element(element, mat_in, &
                    imo_super_i, & ! row
                    imo_super_j) ! col

             ! mjsp: put element from mat_in into mat_dst
             call dense_put_element(element, mat_dst, &
!                  imo_frag_i, & ! row
!                  imo_frag_j) ! col
                  imo_super_i, & ! row
                  imo_super_j) ! col

!             ! mjsp: Iterate fragment index counter
!             imo_frag_i = imo_frag_i + 1

             ! mjsp: skip MOs if neccesary
             if ((imo_super_i-mo_start+1) .eq. &
                  pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)) then
                imo_super_i = imo_super_i + 1 + join_skip_mo
             else
                imo_super_i = imo_super_i + 1
             end if

          end do ! i++

!          ! mjsp: Iterate fragment index counter
!          imo_frag_j = imo_frag_j + 1
!          imo_frag_i = imo_frag_j ! symmetry

          ! mjsp: skip MOs if neccesary
          if ((imo_super_j-mo_start+1) .eq. &
               pub_frag_data(frag_id)%rep%n_occ(is,PUB_1K)) then
             imo_super_j = imo_super_j + 1 + join_skip_mo
          else
             imo_super_j = imo_super_j + 1
          end if

       end do ! j++

    end if


    ! Stop Timer
    call timer_clock('fragment_matrix_ops_blk_mxm',2)

  end subroutine fmo_get_frag_blk_dens_mxm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_put_frag_blk_dens_nxn(mat_dst, mat_in, frag_id, frag_id2, &
                  frag_loc)

    !========================================================================!
    ! This subroutine copies a dense matrix fragment block to a              !
    ! supermolecule dense matrix.                                            !
    ! where X is the fragment requested via frag_id.                         !
    ! if frag_loc==.true. then the dense matrix input                        !
    !   will has dimensions N_ngwfs(X) x N_ngwfs(X) where X is the fragment  !
    !   requested via frag_i.                                                !
    ! if frag_loc==.false. then the dense matrix input                       !
    !   will has dimensions N_ngwfs(A) x N_ngwfs(A) where A is the           !
    !   supermolecule (constructed from the fragments X)                     !
    ! NOTE: This is essentially a reversal of fmo_get_frag_blk_dens_nxn()    !
    !------------------------------------------------------------------------!
    ! mat_dst   (inout)   : Matrix to receive.                               !
    ! mat_in    (inout)   : Matrix to be copied. (N(frag)xN(frag))           !
    ! frag_id   (in)      : ID of the fragment block to be copied.           !
    ! frag_id2  (in)      : ID of the additional fragment block to be copied.!
    ! frag_loc  (in)      : Whether to copy to the fragment pattern          !
    !                       (used for avoiding rank deficiency issues).      !
    !------------------------------------------------------------------------!
    ! Written by Max Phipps June 2015.                                       !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    use datatypes, only: COEF    ! agrecocmplx
    use dense, only: DEM, dense_put_element, dense_get_element
    use fragment_data, only: pub_frag_data
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout)      :: mat_dst
    type(DEM), intent(in)         :: mat_in
    integer, intent(in)           :: frag_id
    integer, intent(in), optional :: frag_id2
    logical, intent(in), optional :: frag_loc

    ! Local Variables
    integer :: ingwf_frag_i, ingwf_frag_j
    !real(kind=DP) :: element
    ! agrecocmplx: use COEF type
    type(COEF) :: element
    integer :: frag_id2_val
    logical :: frag_loc_log
    integer :: ingwf_super_i, ingwf_super_j
    character(len=*), parameter :: myself = 'fmo_put_frag_blk_dens_nxn'

    ! Start Timer
    call timer_clock('fragment_matrix_ops_put_blk_nxn',1)

    ! --------------------------------------------------------------------------

    ! jme: initialize complex character of the copy
    call utils_assert(mat_dst%iscmplx .eqv. mat_in%iscmplx, &
         'Error in '//myself//': real and complex matrices are being mixed.')
    element%iscmplx = mat_in%iscmplx

    frag_id2_val = 0  ! redundant if frag_id2 optional argument not used
    if (present(frag_id2)) then
       if (frag_id2 .ne. 0) then
          frag_id2_val = frag_id2
       end if
    end if
    frag_loc_log = .false.
    if (present(frag_loc)) frag_loc_log = frag_loc

    ! mjsp: Initialise
    ingwf_frag_i = 1
    ingwf_frag_j = 1

    ! mjsp: Only copying one fragment
    if (frag_id2_val .eq. 0) then

       ! mjsp: frag_id x frag_id block:

       ! mjsp: iterate NGWFs on the supermolecule
       do ingwf_super_j=1,pub_frag_data(0)%ngwf_basis(1)%num ! COL

          ! mjsp: if frag_id NGWF is on supermolecule
          if (pub_super2frag_ngwf_log(ingwf_super_j,frag_id)) then

             ! mjsp: iterate NGWFs on the supermolecule
             do ingwf_super_i=1,pub_frag_data(0)%ngwf_basis(1)%num ! ROW

                ! mjsp: if frag_id NGWF is on supermolecule
                if (pub_super2frag_ngwf_log(ingwf_super_i,frag_id)) then

                   ! supermolecule dimensions to fragment dimensions:
                   if (frag_loc_log) then

                      ! mjsp: get element from mat_in
                      call dense_get_element(element, mat_in, &
                           ingwf_frag_i, & ! row
                           ingwf_frag_j) ! col

                   ! supermolecule dimensions to supermolecule dimensions:
                   else
                      ! mjsp: get element from mat_in
                      call dense_get_element(element, mat_in, &
                           ingwf_super_i, & ! row
                           ingwf_super_j) ! col
                   end if

                   ! mjsp: put element from mat_in into mat_dst
                   call dense_put_element(element, mat_dst, &
                        ingwf_super_i, & ! row
                        ingwf_super_j) ! col

                   ! mjsp: Iterate fragment index counter
                   ingwf_frag_i = ingwf_frag_i + 1

                end if
             end do ! i++

             ! mjsp: Iterate fragment index counter
             ingwf_frag_j = ingwf_frag_j + 1
             ingwf_frag_i = 1

          end if
       end do ! j++

    ! mjsp: If joining two fragments together:
    else

       ! mjsp: frag_id x frag_id block +
       ! mjsp: frag_id x frag_id2 block +
       ! mjsp: frag_id2 x frag_id block +
       ! mjsp: frag_id2 x frag_id2 block

       ! mjsp: iterate NGWFs on the supermolecule
       do ingwf_super_j=1,pub_frag_data(0)%ngwf_basis(1)%num ! COL

          ! mjsp: if frag_id or frag_id2 NGWF is on supermolecule
          if ((pub_super2frag_ngwf_log(ingwf_super_j,frag_id)) .or. &
              (pub_super2frag_ngwf_log(ingwf_super_j,frag_id2))) then

             ! mjsp: iterate NGWFs on the supermolecule
             do ingwf_super_i=1,pub_frag_data(0)%ngwf_basis(1)%num ! ROW

                ! mjsp: if frag_id or frag_id2 NGWF is on supermolecule
                if ((pub_super2frag_ngwf_log(ingwf_super_i,frag_id)) .or. &
                    (pub_super2frag_ngwf_log(ingwf_super_i,frag_id2))) then

                   ! supermolecule dimensions to fragment dimensions:
                   if (frag_loc_log) then

                      ! mjsp: get element from mat_in
                      call dense_get_element(element, mat_in, &
                           ingwf_frag_i, & ! row
                           ingwf_frag_j) ! col

                   ! supermolecule dimensions to supermolecule dimensions:
                   else
                      ! mjsp: get element from mat_in
                      call dense_get_element(element, mat_in, &
                           ingwf_super_i, & ! row
                           ingwf_super_j) ! col
                   end if

                   ! mjsp: put element from mat_in into mat_dst
                   call dense_put_element(element, mat_dst, &
                        ingwf_super_i, & ! row
                        ingwf_super_j) ! col

                   ! mjsp: Iterate fragment index counter
                   ingwf_frag_i = ingwf_frag_i + 1

                end if
             end do ! i++

             ! mjsp: Iterate fragment index counter
             ingwf_frag_j = ingwf_frag_j + 1
             ingwf_frag_i = 1

          end if
       end do ! j++

    end if

    ! --------------------------------------------------------------------------

    ! Stop Timer
    call timer_clock('fragment_matrix_ops_put_blk_nxn',2)

  end subroutine fmo_put_frag_blk_dens_nxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_load_internal_blkdiag_denskern(denskern, ngwf_basis)

    !========================================================================!
    ! This subroutine inserts the fragment density kernel elements           !
    ! into the appropriate elements of the supermolecule density kernel.     !
    !------------------------------------------------------------------------!
    ! Written by Max Phipps June 2015.                                       !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    use datatypes, only: COEF, data_coef_scale
    use comms, only: pub_my_proc_id
    use constants, only: stdout
    use dense, only: dense_get_element
    use fragment_data, only: NGWF_INDEX_MAP_TYPE, pub_frag_data, &
         fragment_data_get_ngwf_index_map
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_debug_on_root, &
         pub_num_spins, pub_frag_iatm, pub_eda_have_sp_frags
    use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_scale
    use sparse, only: sparse_put_element
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED_ARRAY), intent(inout) :: denskern
    type(FUNC_BASIS), intent(in) :: ngwf_basis(1)

    ! Local Variables
    integer :: is, it, ngwf_offset
    integer :: ingwf_frag_i, ingwf_frag_j, el_i, el_j
    integer :: el_i_proc
    type(COEF) :: element
    logical :: frag_is_sp  ! if a fragment is spin-polarised or not
    character(len=*), parameter :: myself='fmo_load_internal_blkdiag_denskern'

    type(NGWF_INDEX_MAP_TYPE) :: ngwf_index_map

    ! Start Timer
    call timer_clock('fragment_matrix_load_dkn',1)

    ! --------------------------------------------------------------------------

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering '//myself

    ngwf_index_map = fragment_data_get_ngwf_index_map()

    ! jme: initialize complex character of elements
#ifdef SCALAPACK
    element%iscmplx = pub_frag_data(1)%denskern_dens(1)%iscmplx
#else
    element%iscmplx = pub_frag_data(1)%denskern_dens_loc(1)%iscmplx
#endif

    do is=1,pub_num_spins ! loop spins

       call sparse_embed_scale(denskern%m(is,PUB_1K),0.0_DP)

       ngwf_offset = 0

       do it=1,pub_frag_iatm(0) ! loop fragments

          ! determine if this fragment is spin=2
          frag_is_sp = .false.
          if (pub_frag_data(it)%num_spins == 2) then
            frag_is_sp = .true.
          end if

          ! iterate NGWFs on this fragment
          do ingwf_frag_i=1,pub_frag_data(it)%ngwf_basis(1)%num

             ! Find the 'i'th index in the supermolecule matrix
             ! that this element is associated with.
             ! (fragment ingwf_frag_i -> supermolecule el_i)
             el_i = ngwf_index_map%supers_by_cum_frag( &
                  ingwf_frag_i + ngwf_offset)

             if(pub_debug_on_root) &
                 write(stdout,'(a,i3,a,i3,a,i3)') "DEBUG: Fragment ID #", it,&
                           " NGWF index (fragment): ", ingwf_frag_i, &
                      " NGWF index (supermolecule): ", el_i


             do ingwf_frag_j=1,pub_frag_data(it)%ngwf_basis(1)%num

                ! mjsp: Find the 'j'th index in the supermolecule matrix
                ! that this element is associated with
                ! (fragment ingwf_frag_j -> supermolecule el_j)
                el_j = ngwf_index_map%supers_by_cum_frag( &
                     ingwf_frag_j + ngwf_offset)

                ! mjsp: Deal with spin polarisation
                ! mjsp: NOTE: "fragment has spin=2 and we have a
                ! spin=1 supermolecule" case is prevented by
                ! having set the supermolecule number of spins
                ! to the maximum of the fragment number of
                ! spins.

                ! mjsp: if have any spin polarised fragments
                if (pub_eda_have_sp_frags) then

                   ! mjsp: get the element from the (spin=1) frozen density kernel:
#ifdef SCALAPACK
                   call dense_get_element(element, pub_frag_data(it)%denskern_dens(1), &
                      ingwf_frag_j,ingwf_frag_i)
#else
                   call dense_get_element(element, pub_frag_data(it)%denskern_dens_loc(1), &
                      ingwf_frag_j,ingwf_frag_i)
#endif

                   ! mjsp: if this fragment has spin=1
                   if (frag_is_sp .eqv. .false.) then

                      ! partition the density of this fragment
                      ! into its spin channel components
                      ! TODO: spin channels
                      call data_coef_scale(element, 0.5_DP)

                   end if

                ! mjsp: if no spin polarised fragments
                else

                   ! mjsp: get the element from the frozen density kernel:
#ifdef SCALAPACK
                   call dense_get_element(element, pub_frag_data(it)%denskern_dens(is), &
                      ingwf_frag_j,ingwf_frag_i)
#else
                   call dense_get_element(element, pub_frag_data(it)%denskern_dens_loc(is), &
                      ingwf_frag_j,ingwf_frag_i)
#endif

                end if

                ! mjsp: The NGWFs are distributed column-wise
                ! over the procs.  We check if el_i is on our
                ! proc, and then loop all column rows and copy
                ! over the appropriate matrix elements from the
                ! fragment calculation.  find if on this proc:
                el_i_proc = ngwf_basis(1)%proc_of_func(el_i)

                if (pub_my_proc_id == el_i_proc) then

                  ! mjsp: copy the density kernel element over to the
                  ! supermolecule
                  call sparse_put_element(element, denskern%m(is,PUB_1K)%p, &
                       el_j, & ! row
                       el_i) ! col

                end if ! if on proc


             end do ! NGWFs on fragment

          end do ! NGWFs on fragment

          ! modify NGWF offset
          ngwf_offset = ngwf_offset + pub_frag_data(it)%ngwf_basis(1)%num

       end do  ! fragments

    end do  ! spin

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving fmo_load_internal_blkdiag_denskern'

    ! Stop Timer
    call timer_clock('fragment_matrix_load_dkn',2)

  end subroutine fmo_load_internal_blkdiag_denskern

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_construct_overlap_mats(rep, ngwf_basis)

    !========================================================================!
    ! This subroutine constructs the fragment overlap and inverse overlap    !
    ! matrices.                                                              !
    !------------------------------------------------------------------------!
    ! Written by Max Phipps January 2016 (from                               !
    !      hamiltonian:hamiltonian_dens_indep_matrices)                      !
    ! Modified for embedding by Joseph Prentice, September 2018              !
    !========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, VERBOSE, EDA_CTFRAGLOC_DEVEL
    use dense, only: DEM, dense_create, dense_destroy, dense_convert
    use function_basis, only: FUNC_BASIS
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_eda_mode, pub_max_resid_hotelling, pub_output_detail, &
         pub_frag_counter, pub_frag_counter2, pub_maxit_hotelling
    use sparse_embed, only: sparse_embed_copy, sparse_embed_hotelling_init, &
         sparse_embed_hotelling_invert
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout)   :: rep
    type(FUNC_BASIS), intent(in)    :: ngwf_basis(1)

    ! Local Variables
    type(DEM) :: inv_overlap_b_dens
    real(kind=DP) :: final_max_resid

    ! mjsp: Prepare overlap and inverse overlap matrices for later
    ! mjsp: locally-projected equations (SCF-MI):

    ! make copies of the full matrices in sparse and dense form
    ! for use in the SCF-MI module:
    ! full overlap...
    call sparse_embed_copy(rep%overlap_scfmi_full,rep%overlap)
    ! full inverse overlap...
    call sparse_embed_copy(rep%inv_overlap_scfmi_full,rep%inv_overlap)

    ! For Stoll equations, zero the overlap matrix and reconstruct the
    ! inverse overlap using this zeroed overlap matrix:

    ! mjsp: Reconstruct rep%overlap by zeroing the off-diagonal
    ! (interfragmental) blocks:
    if (pub_eda_mode == EDA_CTFRAGLOC_DEVEL) then
       ! if fragment pair delocalisation stage, then treat pub_frag_counter
       ! and pub_frag_counter2 as a single fragment.
       call fmo_freeze_mat_spam_nxn(rep%overlap, pub_frag_counter, &
                pub_frag_counter2)
    else
       call fmo_freeze_mat_spam_nxn(rep%overlap,0)
    end if



    ! mjsp: Perform the overlap matrix inversion...

    ! if this is the first call to this routine...
    if (.not.rep%inv_overlap_b_init) then
       ! cks: initialise inverse overlap guess if Hotelling recursion
       ! cks: will be used to approximate it.
       ! ndmh: call sparse_mod version of hotelling's algorithm initialisation
       if (pub_maxit_hotelling > 0) call sparse_embed_hotelling_init(rep%inv_overlap, &
            rep%overlap)
       rep%inv_overlap_b_init = .true.
    end if

    ! cks : approximate inverse overlap by Hotelling's recursion
    if (pub_maxit_hotelling > 0) then

       if (pub_output_detail>=VERBOSE .and. pub_on_root) then
          write(stdout,'(/a)')'======== Calculation of Fragment NGWF S^-1 using &
               &Hotelling algorithm ============ '
          write(stdout,'(a)')'         Iteration  |  Frobenius norm and Maximum &
               &abs value    (of I-S*S_n^-1)   '
       end if

       ! ndmh: call sparse_mod version of hotelling's algorithm
       call sparse_embed_hotelling_invert(rep%inv_overlap, & !inout
            rep%overlap,show_output=(pub_output_detail>=VERBOSE), &
            max_resid_converged=pub_max_resid_hotelling, &
            num_iter=pub_maxit_hotelling,final_max_resid=final_max_resid)

       ! ndmh: test if Hotelling failed to invert the matrix
       if (final_max_resid > 0.95_DP) then

          if (pub_output_detail>=VERBOSE.and.pub_on_root) then
             write(stdout,'(15x,a)')'Resetting Hotelling algorithm'
          end if

          ! ndmh: re-initialise Hotelling algorithm from scratch
          call sparse_embed_hotelling_init(rep%inv_overlap,rep%overlap)

          ! ndmh: call sparse_mod version of hotelling's algorithm
          call sparse_embed_hotelling_invert(rep%inv_overlap, & !inout
               rep%overlap,show_output=(pub_output_detail>=VERBOSE), &
               max_resid_converged=pub_max_resid_hotelling, &
               num_iter=pub_maxit_hotelling,final_max_resid=final_max_resid)

          ! ndmh: check if it is still broken
          if (final_max_resid > 0.95_DP) then
             call utils_abort('Error in hamiltonian_dens_indep_matrices: &
                  &Inversion of overlap matrix failed')
          end if

       end if


       if (pub_output_detail>=VERBOSE.and.pub_on_root) then
          write(stdout,'(a)')'===================================&
               &============================================='
       end if

    else if (pub_maxit_hotelling == 0) then

       ! mjsp: Reconstruct inverse overlap from fragment-localised overlap
       ! matrices (this will effectively be a block diagonal matrix of the
       ! individual fragments' inverse overlaps):
       call dense_create(inv_overlap_b_dens, ngwf_basis(1)%num, ngwf_basis(1)%num)
       call dense_convert(inv_overlap_b_dens,rep%overlap)
       call fmo_invert_frag_block_dens_nxn(inv_overlap_b_dens)
       ! NOTE: fmo_invert_frag_block_dens_nxn equivalent to:
       !call dense_invert(pub_scfmi_data%inv_overlap_b_dens)
       ! but with reduced computational expense by looping over the fragment
       ! blocks.

       ! copy inverted fragment-localised overlap matrix to
       ! sparse:
       ! (i.e. construct rep%inv_overlap)
       call dense_convert(rep%inv_overlap, inv_overlap_b_dens)
       call dense_destroy(inv_overlap_b_dens)

    else
       call utils_abort('Error in hamiltonian_dens_indep_matrices: &
            &negative maxit_hotelling ',pub_maxit_hotelling)
    end if


  end subroutine fmo_construct_overlap_mats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fmo_apply_intrafrag_partitions()

    !========================================================================!
    ! This subroutine applies the modifications to the fragment data         !
    ! neccesary to partition bonds in multi-fragment molecules.              !
    !------------------------------------------------------------------------!
    ! Written by Max Phipps November 2015.                                   !
    !========================================================================!

    use constants, only: stdout
    use rundat, only: pub_debug_on_root

    implicit none

    ! Local Variables
    integer :: is, it

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Entering fmo_apply_intrafrag_partitions'

!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>
!! TODO: Fragment partitioning is under development
!! This is a placeholder
!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving fmo_apply_intrafrag_partitions'

  end subroutine fmo_apply_intrafrag_partitions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module fragment_matrix_ops


