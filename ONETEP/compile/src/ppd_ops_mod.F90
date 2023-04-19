! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!-----------------------------------------------------------------------!
!                                                                       !
! Extracted from hf_exchange by Jacek Dziedzic in January 2014.         !
!                                                                       !
!-----------------------------------------------------------------------!
!                                                                       !
! This module implements a number of operations on PPDs.                !
! Some of these operate only on PPD indices (ppd_set_intersection),     !
! while others (the majority) operate on indices and PPD contents.      !
!                                                                       !
!-----------------------------------------------------------------------!


module ppd_ops

  use constants, only: DP

  implicit none
  private

  public :: ppd_product
  public :: ppd_product_unpacked
  public :: ppd_product_one_ppd
  public :: ppd_set_union
  public :: ppd_set_union_unordered
  public :: ppd_set_union_unordered_fast
  public :: ppd_set_intersection
  public :: ppd_set_intersection_unpacked
  public :: ppd_extract_contents_of_subset
  public :: ppd_accumulate

  ! jd: This is the maximum number of PPDs that can ever be encountered in a
  !     ppd_set_union operation. If this value is insufficient, there is
  !     an assert to trap that. Don't make it too big, as MAX_PPDS_OF_INTEREST
  !     integers are used in a local stack-based array. I don't think it's
  !     worth it to generalise this to a keyword (it would no longer be a
  !     compile-time constant and would move the datastructures from the stack
  !     to the heap). We only allocate arrays of integer indices dimensioned
  !     to MAX_PPDS_OF_INTEREST. Keep it that way.
  !     It is public, because it's used in hf_exchange too.
  integer, parameter, public :: MAX_PPDS_OF_INTEREST = 50000

  ! ---------------------------------------------------------------------------

contains

  ! ---------------------------------------------------------------------------
  ! ----- P U B L I C   S U B R O U T I N E S   A N D   F U N C T I O N S -----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_product(packed_ngwf_bra, packed_ngwf_ket, cell_npts, &     ! in
       prod_on_grid, ppd_indices_in_prod, n_ppds_in_prod, &           ! out, opt
       use_1_for_bra)                                                 ! in, opt
    !==================================================================!
    ! This function calculates a pointwise product of points in the    !
    ! ppds in common between a bra in packed NGWF representation and   !
    ! a ket in packed NGWF representation. It can also be used to      !
    ! determine just the indices of the PPDs in the product.           !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   packed_ngwf_bra (input): The bra NGWF in a packed represent'n. !
    !   packed_ngwf_ket (input): The ket NGWF in a packed represent'n. !
    !   cell_npts (input): cell%n_pts.                                 !
    !   prod_on_grid (out, opt): The resultant product in ppd rep.     !
    !                            Can be omitted, if not needed by the  !
    !                            caller.                               !
    !   ppd_indices_in_prod (out, opt): The indices of the ppds in the !
    !                                   product. Can be omitted, if not!
    !                                   needed by the caller.          !
    !   n_ppds_in_prod (out, opt): The number of ppds in the product.  !
    !                              Can be omitted if not needed by the !
    !                              caller.                             !
    !   use_1_for_bra (input, opt): If passed and .true., the contents !
    !                               of the bra will be ignored and     !
    !                               instead replaced by 1.0.           !
    !------------------------------------------------------------------!
    ! Cloned from integrals::internal_brappd_ketppd and modified by    !
    ! Jacek Dziedzic in June 2012.                                     !
    ! Reworked to use packed NGWFs by Jacek Dziedzic in November 2013. !
    ! Changed not to assume max_n_ppds everywhere in February 2014.    !
    !==================================================================!

    use basis, only: SPHERE, basis_sphere_deallocate
    use datatypes, only: FUNCTIONS, data_functions_dealloc
    use remote, only: remote_unpack_ngwf
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(in)               :: packed_ngwf_bra(:)
    real(kind=DP), intent(in)               :: packed_ngwf_ket(:)
    integer, intent(in)                     :: cell_npts
    real(kind=DP), intent(out), optional    :: prod_on_grid(:)
    integer, intent(out), optional          :: ppd_indices_in_prod(:)
    integer, intent(out), optional          :: n_ppds_in_prod
    logical, intent(in), optional           :: use_1_for_bra

    ! Local Variables
    integer :: ibra, iket              ! index in bra and ket ppd lists
    integer :: bra_ppd, ket_ppd        ! ppd index of bra and ket ppds
    integer :: bra_start, ket_start    ! index of first point of ppd
    integer :: prod_start              ! jd: index for output

    type(FUNCTIONS)  :: bra_on_grid ! only 1 NGWF
    type(FUNCTIONS)  :: ket_on_grid ! only 1 NGWF
    type(SPHERE)     :: bra_sphere
    type(SPHERE)     :: ket_sphere
    integer          :: n_bra_ppds, n_ket_ppds
    logical          :: local_use_1_for_bra
    character(len=*), parameter :: myself = 'ppd_product'

    ! -------------------------------------------------------------------------

    ! NB. This subroutine is called from OMP contexts. Do not utils_trace() it.

    if(present(ppd_indices_in_prod)) then
       call utils_assert(present(n_ppds_in_prod), myself//&
            ': n_ppds_in_prod must not be omitted when ppd_indices_in_prod&
            &is present.')
    end if

    call remote_unpack_ngwf(packed_ngwf_bra, bra_sphere, bra_on_grid)
    call remote_unpack_ngwf(packed_ngwf_ket, ket_sphere, ket_on_grid)
    n_bra_ppds = bra_sphere%n_ppds_sphere
    n_ket_ppds = ket_sphere%n_ppds_sphere

    if(present(use_1_for_bra)) then
       local_use_1_for_bra = use_1_for_bra
    else
       local_use_1_for_bra = .false.
    end if

    call utils_assert(.not. (bra_on_grid%iscmplx .or. ket_on_grid%iscmplx), &
         myself//': Not ready for complex NGWFs -- changes would propagate to &
         &hf_exchange', bra_on_grid%iscmplx, ket_on_grid%iscmplx)

    ! Start at the beginning of the ket ppd list
    iket = 1
    ket_ppd = ket_sphere%ppd_list(1,iket)

    ! jd: Product contains zero ppds initially
    prod_start = 1
    if(present(n_ppds_in_prod)) n_ppds_in_prod = 0

    do ibra=1,bra_sphere%n_ppds_sphere
       ! ndmh: find ppd number and start position for ibra
       bra_ppd = bra_sphere%ppd_list(1,ibra)
       bra_start = bra_sphere%offset + (ibra-1)*cell_npts
       do
          ! ndmh: keep moving on while ket_ppd is less than bra_ppd
          if (ket_ppd < bra_ppd) then
             iket = iket + 1
             if (iket > ket_sphere%n_ppds_sphere) exit
             ket_ppd = ket_sphere%ppd_list(1,iket)
             cycle
          end if

          ! ndmh: ppd numbers match, so add product of these ppds
          if (ket_ppd == bra_ppd) then
             ket_start = ket_sphere%offset + (iket-1)*cell_npts
             if(.not. local_use_1_for_bra) then
                if(present(prod_on_grid)) then
                   prod_on_grid(prod_start:prod_start+cell_npts-1) = &
                        bra_on_grid%d(bra_start:bra_start+cell_npts-1) * &
                        ket_on_grid%d(ket_start:ket_start+cell_npts-1)
                end if
             else
                if(present(prod_on_grid)) then
                   prod_on_grid(prod_start:prod_start+cell_npts-1) = &
                        ket_on_grid%d(ket_start:ket_start+cell_npts-1)
                end if
             end if
             if(present(n_ppds_in_prod)) n_ppds_in_prod = n_ppds_in_prod + 1
             prod_start = prod_start + cell_npts
             if(present(ppd_indices_in_prod)) then
                ppd_indices_in_prod(n_ppds_in_prod) = ket_ppd
             end if
          end if

          ! ndmh: move on to next bra ppd as soon as we have found a match
          ! ndmh: or moved beyond bra_ppd
          if (ket_ppd >= bra_ppd) exit

       end do
    end do

    ! jd: Zero the tail of the array to keep overzealous compilers satisfied
    if(present(prod_on_grid)) then
       prod_on_grid(prod_start:) = 0D0
    end if

    call basis_sphere_deallocate(bra_sphere)
    call basis_sphere_deallocate(ket_sphere)
    call data_functions_dealloc(bra_on_grid)
    call data_functions_dealloc(ket_on_grid)

  end subroutine ppd_product

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_product_unpacked(&
       prod_on_grid, ppd_indices_in_prod, n_ppds_in_prod, &
       bra_on_grid, ppd_indices_in_bra, n_ppds_in_bra, &
       ket_on_grid, ppd_indices_in_ket, n_ppds_in_ket, &
       cell_npts, use_1_for_bra)
    !==================================================================!
    ! This function calculates a pointwise product of points in the    !
    ! ppds in common between a bra in PPD representation and a ket in  !
    ! PPD representation.                                              !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   prod_on_grid (output): The resultant product in ppd rep.       !
    !   ppd_indices_in_prod (output): The indices of the ppds in the   !
    !                                 product.                         !
    !   n_ppds_in_prod (output): The number of ppds in the product.    !
    !   bra_on_grid, ppd_indices_in_bra, n_ppds_in_bra: Define the bra !
    !                                                   NGWF.          !
    !   ket_on_grid, ppd_indices_in_ket, n_ppds_in_ket: Define the ket !
    !                                                   NGWF.          !
    !   cell_npts (input): cell%n_pts.                                 !
    !   use_1_for_bra (input, opt): If passed and .true., the contents !
    !                               of the bra will be ignored and     !
    !                               instead replaced by 1.0.           !
    !------------------------------------------------------------------!
    ! Caveats:                                                         !
    !   spheres are inaccessible, so offsets are assumed to be 1.      !
    !------------------------------------------------------------------!
    ! Cloned from ppd_product() by Jacek Dziedzic in February 2017.    !
    !==================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(out)    :: prod_on_grid(:)
    integer, intent(out)          :: ppd_indices_in_prod(:)
    integer, intent(out)          :: n_ppds_in_prod
    real(kind=DP), intent(in)     :: bra_on_grid(:)
    integer, intent(in)           :: ppd_indices_in_bra(:)
    integer, intent(in)           :: n_ppds_in_bra
    real(kind=DP), intent(in)     :: ket_on_grid(:)
    integer, intent(in)           :: ppd_indices_in_ket(:)
    integer, intent(in)           :: n_ppds_in_ket
    integer, intent(in)           :: cell_npts
    logical, intent(in), optional :: use_1_for_bra

    ! Local Variables
    integer :: ibra, iket              ! index in bra and ket ppd lists
    integer :: bra_ppd, ket_ppd        ! ppd index of bra and ket ppds
    integer :: bra_start, ket_start    ! index of first point of ppd
    integer :: prod_start              ! jd: index for output

    logical          :: local_use_1_for_bra
    character(len=*), parameter :: myself = 'ppd_product_unpacked'

    ! -------------------------------------------------------------------------

    if(present(use_1_for_bra)) then
       local_use_1_for_bra = use_1_for_bra
    else
       local_use_1_for_bra = .false.
    end if

    ! Start at the beginning of the ket ppd list
    iket = 1
    ket_ppd = ppd_indices_in_ket(iket)

    ! jd: Product contains zero ppds initially
    prod_start = 1
    n_ppds_in_prod = 0

    do ibra = 1, n_ppds_in_bra
       ! ndmh: find ppd number and start position for ibra
       bra_ppd = ppd_indices_in_bra(ibra)
       bra_start = 1 + (ibra-1)*cell_npts
       do
          ! ndmh: keep moving on while ket_ppd is less than bra_ppd
          if (ket_ppd < bra_ppd) then
             iket = iket + 1
             if (iket > n_ppds_in_ket) exit
             ket_ppd = ppd_indices_in_ket(iket)
             cycle
          end if

          ! ndmh: ppd numbers match, so add product of these ppds
          if (ket_ppd == bra_ppd) then
             ket_start = 1 + (iket-1)*cell_npts
             if(.not. local_use_1_for_bra) then
                prod_on_grid(prod_start:prod_start+cell_npts-1) = &
                     bra_on_grid(bra_start:bra_start+cell_npts-1) * &
                     ket_on_grid(ket_start:ket_start+cell_npts-1)
             else
                prod_on_grid(prod_start:prod_start+cell_npts-1) = &
                     ket_on_grid(ket_start:ket_start+cell_npts-1)
             end if
             n_ppds_in_prod = n_ppds_in_prod + 1
             prod_start = prod_start + cell_npts
             ppd_indices_in_prod(n_ppds_in_prod) = ket_ppd
          end if

          ! ndmh: move on to next bra ppd as soon as we have found a match
          ! ndmh: or moved beyond bra_ppd
          if (ket_ppd >= bra_ppd) exit

       end do
    end do

    ! jd: Zero the tail of the array to keep overzealous compilers satisfied
    prod_on_grid(prod_start:) = 0D0

  end subroutine ppd_product_unpacked

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_product_one_ppd(packed_ngwf_bra, packed_ngwf_ket, &
       cell_npts, cur_ppd, ppd_product, cache_offset_bra, cache_offset_ket)
    !========================================================================!
    ! This function calculates a pointwise product of points in a single PPD !
    ! shared by two NGWFs in a packed NGWF representation. This is not very  !
    ! efficient, since the PPD must be looked up in the index first.         !
    !                                                                        !
    ! However, a simple optimisation of storing the offset to the PPD        !
    ! between calls that refer to the same atom numbers (e.g. calls involving!
    ! different NGWFs a and d on the same atoms A and D) alleviates some of  !
    ! the cost.                                                              !
    !                                                                        !
    ! This subroutine makes some assumptions about the internal guts of      !
    ! packed NGWFs -- namely, how data is laid out after the header. Still,  !
    ! it uses remote's packed_ngwf_header_length to minimise the ugliness.   !
    ! All this is done in the name of efficiency.                            !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    !   packed_ngwf_bra (in): The bra NGWF in a packed representation.       !
    !   packed_ngwf_ket (in): The ket NGWF in a packed representation.       !
    !   cell_npts (in): cell%n_pts.                                          !
    !   cur_ppd (in): Index of the PPD.                                      !
    !   ppd_product (out): The product goes here. Caller is responsible for  !
    !                      providing an array that is cell_npts elements big.!
    !   cache_offset_bra (inout): } If -1 is passed, the PPD will be looked  !
    !   cache_offset_ket (inout): } up in the packed NGWFs in O(n_ppds) time.!
    !                             } On exit, the offset to the PPD will be   !
    !                             } returned. If value other than -1 is      !
    !                             } passed, it will be reused as an offset.  !
    !------------------------------------------------------------------------!

    use remote, only: packed_ngwf_header_length
    use utils, only: utils_abort, utils_int_to_str

    implicit none

    ! Arguments
    real(kind=DP), intent(in)               :: packed_ngwf_bra(:)
    real(kind=DP), intent(in)               :: packed_ngwf_ket(:)
    integer, intent(in)                     :: cell_npts
    integer, intent(in)                     :: cur_ppd
    real(kind=DP), intent(out)              :: ppd_product(cell_npts)
    integer, intent(inout)                  :: cache_offset_bra
    integer, intent(inout)                  :: cache_offset_ket

    ! jd: Local variables
    logical :: found
    integer :: i_bra_ppd, i_ket_ppd
    integer :: n_bra_ppds, n_ket_ppds
    integer :: bra_ppd, ket_ppd
    integer :: src_bra, src_ket, src
    character(len=*), parameter :: myself = 'ppd_product_one_ppd'

    ! -------------------------------------------------------------------------

    n_bra_ppds = nint(packed_ngwf_bra(5))
    n_ket_ppds = nint(packed_ngwf_ket(5))

    if(cache_offset_bra == -1) then
       ! jd: Find cur_ppd in bra NGWF
       src = packed_ngwf_header_length + 1
       found = .false.
       do i_bra_ppd = 1, n_bra_ppds
          bra_ppd = nint(packed_ngwf_bra(src))
          if(bra_ppd == cur_ppd) then
             found = .true.
             exit
          end if
          src = src + 1
       end do
       if(.not. found) then
          call utils_abort(myself//': PPD '//trim(utils_int_to_str(cur_ppd))//&
               ' not found in bra NGWF.')
       end if
    else
       i_bra_ppd = cache_offset_bra
    end if

    ! jd: Find cur_ppd in ket NGWF
    if(cache_offset_ket == -1) then
       src = packed_ngwf_header_length + 1
       found = .false.
       do i_ket_ppd = 1, n_ket_ppds
          ket_ppd = nint(packed_ngwf_ket(src))
          if(ket_ppd == cur_ppd) then
             found = .true.
             exit
          end if
          src = src + 1
       end do
       if(.not. found) then
          call utils_abort(myself//': PPD '//trim(utils_int_to_str(cur_ppd))//&
               ' not found in ket NGWF.')
       end if
    else
       i_ket_ppd = cache_offset_ket
    end if

    ! jd: 2, because we must also skip packed ppd_list(2,:)
    src_bra = packed_ngwf_header_length + 1 + 2 * n_bra_ppds + &
         (i_bra_ppd-1) * cell_npts
    src_ket = packed_ngwf_header_length + 1 + 2 * n_ket_ppds + &
         (i_ket_ppd-1) * cell_npts

    ppd_product(1:cell_npts) = &
         packed_ngwf_bra(src_bra:src_bra+cell_npts-1) * &
         packed_ngwf_ket(src_ket:src_ket+cell_npts-1)

    cache_offset_bra = i_bra_ppd
    cache_offset_ket = i_ket_ppd

  end subroutine ppd_product_one_ppd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_set_union(ppd_set_1, n_ppds_1, ppd_set_2, n_ppds_2)
    !==================================================================!
    ! Merges two sets of PPDs. On exit ppd_set_1 contains all PPDs     !
    ! present in ppd_set_1 and ppd_set_2 on entry, with no duplicates  !
    ! and sorted, while n_ppds_1 is set to the number of PPDs.         !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   ppd_set_1 (in/out): On entry, an array of PPD indices to merge.!
    !                       On exit, an array of merged PPDs.          !
    !   n_ppds_1 (in/out):  On entry, the number of PPDs in ppd_set_1. !
    !                       On exit, the no. of PPDs in new ppd_set_1. !
    !   ppd_set_2 (in):     The second array of PPD indices to merge.  !
    !   n_ppds_2 (in):      The number of PPDs in ppd_set_2.           !
    !------------------------------------------------------------------!
    ! CAVEAT:                                                          !
    !   Both input sets must be sorted, which is always the case with  !
    !   PPDs coming from single spheres in ONETEP. However, if the PPD !
    !   sets span multiple spheres, ppd_set_union_unordered must be    !
    !   used.                                                          !
    !------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in mid 2013.                           !
    !==================================================================!

    use rundat, only: pub_debug
    use utils, only: utils_assert, utils_abort

    implicit none

    ! jd: Arguments
    integer, intent(inout) :: ppd_set_1(:)
    integer, intent(inout) :: n_ppds_1
    integer, intent(in)    :: ppd_set_2(:)
    integer, intent(in)    :: n_ppds_2

    ! jd: Local variables
    integer                :: workspace(1:MAX_PPDS_OF_INTEREST)
    integer                :: pos1, pos2, pos3
    integer                :: val1, val2
    logical                :: found
    character(len=*), parameter :: myself = 'ppd_set_union'

    ! -------------------------------------------------------------------------

    pos1 = 1
    pos2 = 1
    pos3 = 1

    ! Merge sets ppd_set_out and ppd_set_in into workspace
    do while (pos1 <= n_ppds_1 .or. pos2 <= n_ppds_2)
       if(pos1 <= n_ppds_1) then
          val1 = ppd_set_1(pos1)
       else
          val1 = 999999999
       end if
       if(pos2 <= n_ppds_2) then
          val2 = ppd_set_2(pos2)
       else
          val2 = 999999999
       end if

       if(val1 < val2) then
          workspace(pos3) = ppd_set_1(pos1)
          pos1 = pos1 + 1
       end if
       if(val2 < val1) then
          workspace(pos3) = ppd_set_2(pos2)
          pos2 = pos2 + 1
       end if
       if(val1 == val2) then
          workspace(pos3) = ppd_set_1(pos1)
          pos1 = pos1 + 1
          pos2 = pos2 + 1
       end if
       pos3 = pos3 + 1
       if(pub_debug) then
          if(pos3>=MAX_PPDS_OF_INTEREST) then
             call utils_abort(myself//': MAX_PPDS_OF_INTEREST exceeded: ', pos3)
          end if
       end if
    end do

    pos3 = pos3-1

    call utils_assert(pos3<=MAX_PPDS_OF_INTEREST, &
         myself//': MAX_PPDS_OF_INTEREST exceeded: ', pos3)
    ! Return the merged set
    ppd_set_1(1:pos3) = workspace(1:pos3)

    ! Set n_ppds_1 on output
    n_ppds_1 = pos3

    if(pub_debug) then
       ! jd: Sanity-check to see if indices are sorted and positive
       do pos1 = 1, n_ppds_1
          if(pos1 < n_ppds_1) then
             if(.not. (ppd_set_1(pos1) <= ppd_set_1(pos1+1))) then
                call utils_abort(myself//': Internal error [1].')
             end if
          end if
          if(ppd_set_1(pos1) <=0 ) then
             call utils_assert(.false., &
                  myself//': Internal error [2].', ppd_set_1(pos1))
          end if
       end do
       ! jd: Also check all indices in ppd_set_2 are now in ppd_set_1
       do pos2 = 1, n_ppds_2
          val2 = ppd_set_2(pos2)
          found = .false.
          do pos1 = 1, n_ppds_1
             val1 = ppd_set_1(pos1)
             if(val1 == val2) then
                found = .true.
                exit
             end if
          end do
          if(.not. found) then
             call utils_assert(.false.,&
                  myself//': Internal error [3].', val2)
          end if
       end do
    end if

  end subroutine ppd_set_union

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_set_union_unordered(uppd_set_1, n_ppds_1, ppd_set_2, n_ppds_2)
    !==================================================================!
    ! Merges PPDs in two sets, without assuming either is ordered.     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   uppd_set_1 (in/out):On entry, an array of PPD indices to merge.!
    !                       On exit, an array of merged PPDs.          !
    !   n_ppds_1 (in/out):  On entry, the number of PPDs in ppd_set_1. !
    !                       On exit, the no. of PPDs in new ppd_set_1. !
    !   ppd_set_2 (in):     The second array of PPD indices to merge.  !
    !   n_ppds_2 (in):      The number of PPDs in ppd_set_2.           !
    !------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in mid 2013.                           !
    !==================================================================!

    implicit none

    ! jd: Arguments
    integer, intent(inout) :: uppd_set_1(:)
    integer, intent(inout) :: n_ppds_1
    integer, intent(in)    :: ppd_set_2(:)
    integer, intent(in)    :: n_ppds_2

    ! jd: Local variables
    integer                :: pos1, pos2, pos3
    integer                :: val1, val2
    logical                :: found

    ! -------------------------------------------------------------------------

    pos3 = n_ppds_1 + 1

    do pos2 = 1, n_ppds_2
       val2 = ppd_set_2(pos2)
       found = .false.
       do pos1 = 1, pos3-1
          val1 = uppd_set_1(pos1)
          if(val1 == val2) then
             found = .true.
             exit
          end if
       end do
       if(.not. found) then
          uppd_set_1(pos3) = val2
          pos3 = pos3 + 1
       end if
    end do

    n_ppds_1 = pos3-1

  end subroutine ppd_set_union_unordered

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_set_union_unordered_fast(uppd_set_1, n_ppds_1, ppd_set_2, &
       n_ppds_2, ppd_index_lookup)
    !==================================================================!
    ! Merges PPDs in two sets, without assuming either is ordered.     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   uppd_set_1 (in/out):On entry, an array of PPD indices to merge.!
    !                       On exit, an array of merged PPDs.          !
    !   n_ppds_1 (in/out):  On entry, the number of PPDs in ppd_set_1. !
    !                       On exit, the no. of PPDs in new ppd_set_1. !
    !   ppd_set_2 (in):     The second array of PPD indices to merge.  !
    !   n_ppds_2 (in):      The number of PPDs in ppd_set_2.           !
    !   ppd_index_lookup (in/out): Array of logicals, one for each PPD !
    !                              in cell. Bitmap way to remember     !
    !                              which PPDs have been added to make  !
    !                              calculating the union faster. Pass  !
    !                              all .false. on first call.          !
    !------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2019.                     !
    !==================================================================!

    implicit none

    ! jd: Arguments
    integer, intent(inout) :: uppd_set_1(:)
    integer, intent(inout) :: n_ppds_1
    integer, intent(in)    :: ppd_set_2(:)
    integer, intent(in)    :: n_ppds_2
    logical, intent(inout) :: ppd_index_lookup(:) ! 1:cell%n_ppds

    ! jd: Local variables
    integer                :: pos2, pos3
    integer                :: val2

    ! -------------------------------------------------------------------------

    pos3 = n_ppds_1 + 1

    do pos2 = 1, n_ppds_2
       val2 = ppd_set_2(pos2)
       if(ppd_index_lookup(val2)) cycle ! jd: Lookup indicates PPD already added
       ! jd: Otherwise it's a new index, no point in scanning
       uppd_set_1(pos3) = val2
       ppd_index_lookup(val2) = .true. ! jd: Update lookup
       pos3 = pos3 + 1
    end do

    n_ppds_1 = pos3-1

  end subroutine ppd_set_union_unordered_fast

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_set_intersection(ppd_indices, n_common_ppds, & ! out
       packed_ngwf_bra, packed_ngwf_ket)
    !==================================================================!
    ! Finds common PPDs of two NGWFs in a packed representation        !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   ppd_indices (output): The common PPD indices found.            !
    !   n_common_ppds (output):  The number of PPDs found.             !
    !   packed_ngwf_bra (input): The bra NGWF in packed representation.!
    !   packed_ngwf_ket (input): The ket NGWF in packed representation.!
    !------------------------------------------------------------------!
    ! Cloned from integrals::internal_brappd_ketppd and modified by    !
    ! Jacek Dziedzic in June 2012. Modified to use packed NGWFs in     !
    ! November 2013.                                                   !
    !==================================================================!

    use basis, only: SPHERE, basis_sphere_deallocate
    use remote, only: remote_unpack_ngwf

    implicit none

    ! Arguments
    integer, intent(out)         :: ppd_indices(:)
    integer, intent(out)         :: n_common_ppds
    real(kind=DP), intent(in)    :: packed_ngwf_bra(:)
    real(kind=DP), intent(in)    :: packed_ngwf_ket(:)

    ! Local Variables
    type(SPHERE) :: bra_sphere
    type(SPHERE) :: ket_sphere
    integer :: n_bra_ppds
    integer :: n_ket_ppds

    ! -------------------------------------------------------------------------

    call remote_unpack_ngwf(packed_ngwf_bra, bra_sphere)
    call remote_unpack_ngwf(packed_ngwf_ket, ket_sphere)
    n_bra_ppds = bra_sphere%n_ppds_sphere
    n_ket_ppds = ket_sphere%n_ppds_sphere

    call ppd_set_intersection_unpacked(ppd_indices, n_common_ppds, &
         bra_sphere%ppd_list(1,:), n_bra_ppds, &
         ket_sphere%ppd_list(1,:), n_ket_ppds)

    call basis_sphere_deallocate(bra_sphere)
    call basis_sphere_deallocate(ket_sphere)

  end subroutine ppd_set_intersection

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_set_intersection_unpacked(ppd_indices, n_common_ppds, & ! out
       ppd_indices_in_1, n_ppds_in_1, ppd_indices_in_2, n_ppds_in_2)
    !==================================================================!
    ! Finds common PPDs of two sets of PPDs.                           !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   ppd_indices (output): The common PPD indices found.            !
    !   n_common_ppds (output):  The number of PPDs found.             !
    !   ppd_indices_in_1: PPD indices for the first set.               !
    !   n_ppds_1: Number of PPDs in the first set.                     !
    !   ppd_indices_in_2: PPD indices for the second set.              !
    !   n_ppds_2: Number of PPDs in the second set.                    !
    !------------------------------------------------------------------!
    ! Extracted from ppd_set_intersection by Jacek Dziedzic, Feb 2017, !
    ! which used parts of itegrals::internal_brappd_ketppd() by ndmh.  !
    !==================================================================!

    implicit none

    ! Arguments
    integer, intent(out)         :: ppd_indices(:)
    integer, intent(out)         :: n_common_ppds
    integer, intent(in)          :: ppd_indices_in_1(:)
    integer, intent(in)          :: n_ppds_in_1
    integer, intent(in)          :: ppd_indices_in_2(:)
    integer, intent(in)          :: n_ppds_in_2

    ! Local Variables
    integer :: ibra, iket              ! index in bra and ket ppd lists
    integer :: bra_ppd, ket_ppd        ! ppd index of bra and ket ppds

    ! -------------------------------------------------------------------------

    ! Start at the beginning of the ket ppd list
    iket = 1
    ket_ppd = ppd_indices_in_2(iket)

    n_common_ppds = 0

    do ibra=1, n_ppds_in_1
       ! ndmh: find ppd number and start position for ibra
       bra_ppd = ppd_indices_in_1(ibra)
       do
          ! ndmh: keep moving on while ket_ppd is less than bra_ppd
          if (ket_ppd < bra_ppd) then
             iket = iket + 1
             if (iket > n_ppds_in_2) exit
             ket_ppd = ppd_indices_in_2(iket)
             cycle
          end if

          ! ndmh: ppd numbers match, so add the index of this ppd
          if (ket_ppd == bra_ppd) then
             n_common_ppds = n_common_ppds + 1
             ppd_indices(n_common_ppds) = ket_ppd
          end if

          ! ndmh: move on to next bra ppd as soon as we have found a match
          ! ndmh: or moved beyond bra_ppd
          if (ket_ppd >= bra_ppd) exit

       end do
    end do

    ! jd: Our responsibility to clear the remainder of the output array
    ppd_indices(n_common_ppds+1:) = -1

  end subroutine ppd_set_intersection_unpacked

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_extract_contents_of_subset(subset_out, &         ! out
       indices_of_subset, n_indices_of_subset, &                  ! in
       set_in, indices_of_set_in, n_indices_of_set_in, cell_npts) ! in
    !==================================================================!
    ! Extracts contents of certain PPDs from a set of PPDs.            !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !   subset_out (out): The contents of extracted PPDs.              !
    !   indices_of_subset (in): The indices of PPDs to be extracted.   !
    !   n_indices_of_subset (in): The length of the above.             !
    !   set_in (in): The contents of PPDs to extract from.             !
    !   indices_of_set_in (in): The indices of PPDs to extract from.   !
    !   n_indices_of_set_in (in): The length of the above.             !
    !   cell_npts (in): cell%n_pts.                                    !
    !------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on Friday the 13th, July 2012.         !
    !==================================================================!

    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)   :: subset_out(:)
    integer, intent(in)          :: indices_of_subset(:)
    integer, intent(in)          :: n_indices_of_subset
    real(kind=DP), intent(in)    :: set_in(:)
    integer, intent(in)          :: indices_of_set_in(:)
    integer, intent(in)          :: n_indices_of_set_in
    integer, intent(in)          :: cell_npts

    ! jd: Local variables
    integer :: i, j, k
    integer :: idx1, idx2
    logical :: found

    ! -------------------------------------------------------------------------

    subset_out(:) = -1

    j = 0
    do i = 1, n_indices_of_subset
       idx1 = indices_of_subset(i)
       ! jd: Look for idx1 in set_in
       do k = j+1, n_indices_of_set_in
          idx2 = indices_of_set_in(k)
          if(idx1 == idx2) then
             j = k
             found = .true.
             exit
          end if
       end do

       if(pub_debug) call utils_assert(found, &
            'Internal error in ppd_extract_contents_of_subset')

       subset_out((i-1)*cell_npts+1:i*cell_npts) = &
            set_in((k-1)*cell_npts+1:k*cell_npts)
    end do

  end subroutine ppd_extract_contents_of_subset

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_accumulate(&
       ppd_sum, ppd_indices_in_sum, n_ppds_in_sum, &
       ppd_term, ppd_indices_in_term, n_ppds_in_term, &
       is_ppd_sum_ordered, cell_npts, scalar_factor, sum_start_ppd, dest_shift)
    !==========================================================================!
    ! Adds the data in ppd_term * scalar_factor to ppd_sum. If scalar_factor   !
    ! is omitted, 1.0 is used by default.                                      !
    !                                                                          !
    ! If the PPDs in ppd_sum come from multiple spheres and are not ordered,   !
    ! pass .false. for is_ppd_sum_ordered.                                     !
    !                                                                          !
    ! If the optional 'sum_start_ppd' is passed, 'ppd_sum' will be scanned     !
    ! not from the first PPD, but from this (ordinal) PPD. On exit this will   !
    ! contain the index (ordinal) of the first PPD in the sum that was present !
    ! in 'ppd_term'. This can be used to avoid unnecessarily re-scanning       !
    ! 'ppd_sum' from the beginning when accumulating terms using the same PPDs,!
    ! e.g. different NGWFs of the same atom, or different spins.               !
    !                                                                          !
    ! If the optional 'dest_shift' is provided, it will be added to the offset !
    ! in the destination.                                                      !
    !                                                                          !
    ! NB: Every PPD in ppd_term must already exist in ppd_sum. This is checked.!
    !     If it doesn't, merge (ppd_set_union) the two sets in advance or      !
    !     obtain an intersection first.                                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in August 2012.                                !
    ! Extended by Jacek Dziedzic in February 2017 with 'scalar_factor'.        !
    ! Extended by Jacek Dziedzic in March 2020 with 'sum_start_ppd'.           !
    ! Extended by Jacek Dziedzic in November 2020 with 'dest_shift', switched  !
    ! to LONG.                                                                 !
    !==========================================================================!

    use constants, only: LONG
    use rundat, only: pub_debug
    use utils, only: utils_assert, utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(inout)        :: ppd_sum(:)
    integer, intent(in)                 :: ppd_indices_in_sum(:)
    integer, intent(in)                 :: n_ppds_in_sum
    real(kind=DP), intent(in)           :: ppd_term(:)
    integer, intent(in)                 :: ppd_indices_in_term(:)
    integer, intent(in)                 :: n_ppds_in_term
    logical, intent(in)                 :: is_ppd_sum_ordered
    integer, intent(in)                 :: cell_npts
    real(kind=DP), intent(in), optional :: scalar_factor
    integer, intent(inout), optional    :: sum_start_ppd
    integer(kind=LONG), intent(in), optional :: dest_shift

    ! jd: Local variables
    integer :: cur_ppd_n, cur_ppd_sum_n ! ordinal numbers of PPDs
    integer :: cur_ppd, cur_ppd_sum     ! indices of PPDs
    integer(kind=LONG) :: offset_src, offset_dest  ! offsets to data
    integer(kind=LONG) :: src_end, dest_end  ! last datum processed
    integer(kind=LONG) :: loc_dest_shift
    integer :: loc_sum_start_ppd        ! where we left off in ppd_sum
    real(kind=DP) :: factor
    logical :: found
    logical :: first_found_ppd

    ! --------------------------------------------------------------------------

    ! NB. This subroutine is called from OMP contexts. Do not utils_trace() it.

    if(pub_debug) then
      ! jd: Temporarily commenting out, these take too long.
!       call utils_sanity_check(ppd_sum,'ppd_term')
!       call utils_sanity_check(ppd_sum,'ppd_sum')
    end if

    if(present(dest_shift)) then
       loc_dest_shift = dest_shift
    else
       loc_dest_shift = 1_LONG
    end if

    factor = 1D0
    if(present(scalar_factor)) factor = scalar_factor
    if(present(sum_start_ppd)) then
       loc_sum_start_ppd = sum_start_ppd
    else
       loc_sum_start_ppd = 1
    end if
    first_found_ppd = .true.

    do cur_ppd_n = 1, n_ppds_in_term
       cur_ppd = ppd_indices_in_term(cur_ppd_n)
       offset_src = (int(cur_ppd_n,kind=LONG) - 1_LONG) * &
            int(cell_npts,kind=LONG) + 1_LONG
       found = .false.
       do cur_ppd_sum_n = loc_sum_start_ppd, n_ppds_in_sum
          cur_ppd_sum = ppd_indices_in_sum(cur_ppd_sum_n)
          if(cur_ppd_sum == cur_ppd) then
             offset_dest = loc_dest_shift + &
                  int(cur_ppd_sum_n-1,kind=LONG) * int(cell_npts,kind=LONG)
             dest_end = offset_dest + int(cell_npts-1,kind=LONG)
             src_end = offset_src + int(cell_npts-1,kind=LONG)

             if(pub_debug) then
                if(offset_dest < lbound(ppd_sum,1,kind=LONG) .or. &
                     dest_end > ubound(ppd_sum,1,kind=LONG)) then
                   call utils_abort("Destination array bound overrun", &
                        avoid_mpi_calls = .true.)
                end if
                if(offset_src < lbound(ppd_term,1,kind=LONG) .or. &
                     src_end > ubound(ppd_term,1,kind=LONG)) then
                   call utils_abort("Source array bound overrun", &
                        avoid_mpi_calls = .true.)
                end if

                if(dest_end - offset_dest /= src_end - offset_src) then
                   call utils_abort("Source/destination array size mismatch", &
                   avoid_mpi_calls = .true.)
                end if
             end if

             ppd_sum(offset_dest:dest_end) = &
                  ppd_sum(offset_dest:dest_end) + &
                  ppd_term(offset_src:src_end) * factor

             if(first_found_ppd) then
                ! jd: If this is the first PPD to be found, save sum_start_ppd
                if(present(sum_start_ppd)) then
                   sum_start_ppd = cur_ppd_sum_n
                   first_found_ppd = .false.
                end if
             end if
             found = .true.
             if(is_ppd_sum_ordered) then
                loc_sum_start_ppd = cur_ppd_sum_n + 1
             end if
             exit
          end if
       end do
       call utils_assert(found,'Internal error in ppd_accumulate.',cur_ppd)

    end do

  end subroutine ppd_accumulate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ppd_ops
