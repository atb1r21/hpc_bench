! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

module sparse_embed

!==========================================================================!
!                       Module: sparse_embed_mod.F90                       !
!==========================================================================!
! This module is designed to handle embedded SPAM3 matrices. This is done  !
! via the SPAM3_EMBED type. Procedures here generally call sparse_mod.F90  !
! for individual sub-blocks of the SPAM3_EMBED matrix.                     !
!==========================================================================!
! Written by Robert Charlton, August 2017, based largely on pre-existing   !
! subroutines within sparse_mod and sparse_array_mod.                      !
! Further contributions by Joseph Prentice.                                !
!==========================================================================!

  use constants, only: DP, stdout
  use rundat, only: pub_embed_debug, PUB_1K
  use sparse, only: SPAM3
  use sparse_array, only: SPAM3_ARRAY

  implicit none

  private

  ! Public subroutines and functions
  public :: sparse_embed_create
  public :: sparse_embed_destroy
  public :: sparse_embed_scale
  public :: sparse_embed_product
  public :: sparse_embed_axpy
  public :: sparse_embed_copy
  public :: sparse_embed_check_hermitian
  public :: sparse_embed_trace
  public :: sparse_embed_transpose
  public :: sparse_embed_conjugate
  public :: sparse_embed_num_rows
  public :: sparse_embed_num_cols
  public :: sparse_embed_create_sub
  public :: sparse_embed_extract_sub
  public :: sparse_embed_num_element
  public :: sparse_embed_extremal_eigenvalue
  public :: sparse_embed_hotelling_init
  public :: sparse_embed_hotelling_invert
  public :: sparse_embed_rms_element
  public :: sparse_embed_max_abs_element
  public :: sparse_embed_transpose_structure
  public :: sparse_embed_outer_product
  public :: sparse_embed_set_to_garbage
  public :: sparse_embed_entrywise_norm
  public :: sparse_embed_any_isnan
  public :: sparse_embed_1norm
  public :: sparse_embed_read
  public :: sparse_embed_write

  ! rc2013: conversion routines
  public :: sparse_embed_array2mat
  public :: sparse_embed_mat2array
  public :: sparse_embed_convert

  public :: sparse_embed_diagnose

  ! rc2013: embedded sparse_array utilities
  ! sparse_array_ops routines
  public :: sparse_embed_array_create
  public :: sparse_embed_array_destroy

  ! sparse_array_inquiry routines
  public :: sparse_embed_array_is_dense
  public :: sparse_embed_array_num_rows

  ! sparse_array_ops routines
  public :: sparse_embed_array_copy
  public :: sparse_embed_array_scale
  public :: sparse_embed_array_axpy
  public :: sparse_embed_array_product
  public :: sparse_embed_array_trace
  public :: sparse_embed_array_transpose
  public :: sparse_embed_array_extremal_eigenvalues
  ! agreco
  !public :: sparse_array_conjugate

  ! sparse_array_utils routines
  public :: sparse_embed_array_write
  public :: sparse_embed_array_read
  !public :: sparse_array_convert_unsegment_real
  !public :: sparse_array_convert_segment_real
  !public :: sparse_array_convert_unsegment_complex
  !public :: sparse_array_convert_segment_complex
  public :: sparse_embed_array_check

  ! rc2013: conversion routines
  public :: sparse_embed_array_create_sub
  public :: sparse_embed_array_extract_sub
  public :: sparse_embed_array_to_sparse_array
  public :: sparse_embed_extract_from_array
  public :: sparse_embed_destroy_extracted_array

  ! Sparse product routines
  interface sparse_embed_product
     module procedure sparse_embed_product_array_array
     module procedure sparse_embed_product_matrix_array
     module procedure sparse_embed_product_array_matrix
     module procedure sparse_embed_product_matrix_matrix
  end interface sparse_embed_product

  ! create routines
  interface sparse_embed_array_create
     module procedure sparse_embed_array_create_array_array
     module procedure sparse_embed_array_create_matrix_array
     module procedure sparse_embed_array_create_array_matrix
     module procedure sparse_embed_array_create_matrix_matrix
  end interface sparse_embed_array_create

  ! AXPY routines
  interface sparse_embed_array_axpy
     module procedure sparse_embed_array_axpy_real_scalar
     module procedure sparse_embed_array_axpy_real_array
     module procedure sparse_embed_array_axpy_matrix_real_scalar
     module procedure sparse_embed_array_axpy_matrix_real_array
  end interface sparse_embed_array_axpy

  interface sparse_embed_axpy
     module procedure sparse_embed_axpy_real
     module procedure sparse_embed_axpy_complex
  end interface sparse_embed_axpy

  interface sparse_embed_scale
     module procedure sparse_embed_scale_real
     module procedure sparse_embed_scale_complex
     module procedure sparse_embed_scale_complex2
  end interface sparse_embed_scale

  ! Scale and shift routines
  interface sparse_embed_array_scale
     module procedure sparse_embed_array_scale_real_scalar
     module procedure sparse_embed_array_scale_real_array
  end interface sparse_embed_array_scale

  ! Copy routines
  interface sparse_embed_array_copy
     module procedure sparse_embed_array_copy_array
     module procedure sparse_embed_array_copy_matrix
  end interface sparse_embed_array_copy

  ! Sparse product routines
  interface sparse_embed_array_product
     module procedure sparse_embed_array_product_array_array
     module procedure sparse_embed_array_product_matrix_array
     module procedure sparse_embed_array_product_array_matrix
  end interface sparse_embed_array_product

  ! Sparse trace routines
  interface sparse_embed_array_trace
     module procedure sparse_embed_array_trace_array_matrix
     module procedure sparse_embed_array_trace_array_array
  end interface sparse_embed_array_trace

  ! Sparse outer product routines
  interface sparse_embed_outer_product
     module procedure sparse_embed_outer_product_real
     module procedure sparse_embed_outer_product_complex
  end interface sparse_embed_outer_product

    ! Conversion routines
  interface sparse_embed_convert
     module procedure sparse_embed_spam3tofull_real
     module procedure sparse_embed_spam3tofull_complex
     module procedure sparse_embed_fulltospam3_real
     module procedure sparse_embed_fulltospam3_complex
  end interface sparse_embed_convert

  ! Sparse write routines
  interface sparse_embed_write
     module procedure sparse_embed_write_scalar
     module procedure sparse_embed_write_vector
     module procedure sparse_embed_write_matrix
  end interface sparse_embed_write

  ! Sparse write routines
  interface sparse_embed_read
     module procedure sparse_embed_read_scalar
     module procedure sparse_embed_read_vector
     module procedure sparse_embed_read_matrix
  end interface sparse_embed_read

  ! Type definition for the array of matrix elements.
  type, public :: SPAM3_EMBED
     integer :: mrows
     integer :: ncols
     type(SPAM3), allocatable :: m(:,:)
     ! rc2013: structure code without embedding terms
     character(len=30) :: structure
     logical :: iscmplx ! TRUE if matrix is complex, otherwise FALSE
     ! rc2013: pointer to SPAM3_EMBED%m(1,1), for elegance
     type(SPAM3), pointer :: p
  end type SPAM3_EMBED

  ! Type definition for pointer to spin components of a SPAM3_ARRAY
  type, private :: SPIN_PTR
      type(SPAM3), pointer :: p
  end type SPIN_PTR

  ! Type definition for the array of matrix elements.
  type, public :: SPAM3_EMBED_ARRAY
     integer :: mrows
     integer :: ncols
     integer :: num_spins = 0
     integer :: num_kpoints = 0
     ! rc2013: each component of m holds a SPAM3_EMBED matrix
     type(SPAM3_EMBED), allocatable :: m(:,:)
     ! rc2013: structure code without embedding terms
     !  -- same as in each SPAM3_EMBED matrix
     character(len=30) :: structure
     !! rc2013: pointer to SPAM3_EMBED_ARRAY%m(:,PUB_1K)%m(1,1), for sanity
     !type(SPIN_PTR), allocatable, dimension(:) :: spin_p
     !! rc2013: desperation: pointer to an array
     !type(SPAM3_EMBED), dimension(:), pointer :: p
  end type SPAM3_EMBED_ARRAY

contains


  !==========================================================================!
  ! This subroutine creates an embedding SPAM3 structure. The appropriate    !
  ! sub-blocks are passed to sparse_create.                                  !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !  new_spam (inout)  : The sparse matrix to be created                     !
  !  spama    (input)  : An optional sparse matrix whose structure is used   !
  !  spamb    (input)  : A second (optional) sparse matrix whose structure   !
  !                      will be used.                                       !
  !  mrows    (input)  : (Optional) No. of rows in new_spam.                 !
  !  ncols    (input)  : (Optional) No. of columns in new_spam.              !
  !  iscmplx  (input)  : A flag to indicate the data stored in this matrix   !
  !                      is complex.                                         !
  ! structure (input)  : Structure label for new_spam. Only used if spama &  !
  !                      spamb are not present.                              !
  !  trans    (input)  : Flag: should we take the transpose structure from   !
  !                      spama to build new_spam? (false by default)         !
  ! arow/bcol (input)  : Optional choice of subscript(s) for structure       !
  !                      codes of row/column matrices to be extracted.       !
  ! brow/acol (input)  : Optional parameters for extracting matrix blocks    !
  !                      from matrix products.                               !
  !(r/c)struc (input)  : Optional specifier for row/col structures. This     !
  !                      forces new_spam's structure codes to match rstruc or!
  !                      cstruc, regardless of matrix size.                  !
  !                      Obsolete (really?)!                                 !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 30/05/2017.                                  !
  ! Modified to correctly include unions of sparsity patterns for products   !
  ! in embedding systems by Joseph Prentice, July 2019                       !
  !==========================================================================!

  subroutine sparse_embed_create(new_spam, spama, spamb, &
       mrows, ncols, iscmplx, rlib, trans, arow, acol, &
       brow, bcol, rstruc, cstruc, arrname) !@adddoc

    use comms, only: pub_total_num_procs
    use constants, only: LONG
    use rundat, only: pub_ngwf_regions_ngroups, pub_debug_on_root
    use sparse, only: sparse_create, sparse_transpose_structure, &
        sparse_destroy, sparse_count_union, sparse_index_union, & ! jcap
        sparse_mat_in_library
    use utils, only: utils_alloc_check, utils_abort, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), target, intent(inout):: new_spam ! Sparse matrix to be created
    type(SPAM3_EMBED), optional, intent(in) :: spama    ! First (optional) matrix to use
    type(SPAM3_EMBED), optional, intent(in) :: spamb    ! Other (optional) matrix to use
    logical, optional, intent(in)           :: iscmplx  ! Flag for complex data (opt)
    integer, optional, intent(out)          :: rlib(:,:)  ! Optional library index
    logical, optional, intent(in)           :: trans    ! Flag for transpose structure
    integer, optional, intent(in)           :: mrows
    integer, optional, intent(in)           :: ncols
    integer, optional, intent(in)           :: arow
    integer, optional, intent(in)           :: acol
    integer, optional, intent(in)           :: brow
    integer, optional, intent(in)           :: bcol
    integer, optional, intent(in)           :: rstruc   ! Structure code to use
    integer, optional, intent(in)           :: cstruc
    character(len=*), optional, intent(in)  :: arrname

    ! Local variables
    logical :: loc_cmplx   ! Local copy of optional iscmplx argument
    integer :: ierr        ! Error flag
    integer :: isub, jsub, irow, jcol, jrow, icol
    integer :: loc_mrows, loc_ncols
    character(len=1) :: isub_str, jsub_str, ksub_str
    character(len=2) :: reg_struc
    logical :: loc_trans
    integer, allocatable :: fake_rlib(:,:)  ! To avoid segmentation faults

    ! jcap: for calculating unions of sparsity patterns for products in embedding calculations
    type(SPAM3) :: tmp_spam
    integer :: tmp_lib1, tmp_lib2, new_lib, ksub
    integer(kind=LONG) :: nze
    integer :: nzb, my_nze, my_nzb
    integer :: seg_nze(0:pub_total_num_procs-1), seg_nzb(0:pub_total_num_procs-1)
    character(len=30) :: tmp_struc
    logical :: spam_exists

    ! rc2013: check the combinations of optional arguments
    call internal_check_args

    if(present(trans)) then
       loc_trans = trans
    else
       loc_trans = .false.
    endif

    if(present(iscmplx)) then
       loc_cmplx = iscmplx
    else
       loc_cmplx = .false.
    end if
    new_spam%iscmplx = loc_cmplx

    ! jcap: default to total number of regions
    loc_mrows = pub_ngwf_regions_ngroups
    loc_ncols = pub_ngwf_regions_ngroups
    if (present(mrows)) loc_mrows=mrows
    if (present(ncols)) loc_ncols=ncols

    ! rc2013: Allocate all structures and parameters
    call sparse_embed_alloc(new_spam, loc_mrows, loc_ncols, spama, spamb, &
         fake_rlib, arow, bcol, arrname=arrname)

    ! rc2013: cycle over all subsystems
    do jsub=1,new_spam%ncols

       do isub=1,new_spam%mrows
          ! rc2013: set the matrix block iterators + structure codes
          call internal_choose_structure

          ! rc2013: structures will depend on which matrices are present
          if((.not. present(spama)) .and. (.not. present(spamb))) then
             new_spam%m(isub,jsub)%structure = &
                  trim(new_spam%structure)//reg_struc
             call sparse_create(new_spam%m(isub,jsub), &
                  iscmplx=loc_cmplx, rlib=fake_rlib(isub,jsub))

          else if(present(spama) .and. (.not. present(spamb))) then
             if(loc_trans) then
                ! rc2013: get the transpose from spama only
                call sparse_transpose_structure( &
                     new_spam%m(isub,jsub)%structure, spama%m(jcol,irow))
                call sparse_create(new_spam%m(isub,jsub), &
                     iscmplx=loc_cmplx, rlib=fake_rlib(isub,jsub))
             else
                call sparse_create(new_spam%m(isub,jsub),spama%m(irow,jcol),&
                     iscmplx=loc_cmplx, rlib=fake_rlib(isub,jsub))
             endif

          else
             ! rc2013: check if the matrix is in the library
             call sparse_mat_in_library(trim(new_spam%structure)//reg_struc, &
                  spam_exists)

             ! jcap: if we are doing a single-region calculation, this bit is easy
             ! rc2013: or if the matrix is already in the library
             if ( ((new_spam%mrows.eq.1) .and. (new_spam%ncols.eq.1)) &
                  .or. (spam_exists) ) then
                ! rc2013: pass the structure code for embedding matrix
                call sparse_create(new_spam%m(isub,jsub),spama%m(irow,icol),&
                     spamb%m(jrow,jcol), iscmplx=loc_cmplx, rlib=fake_rlib(isub,jsub), &
                     embed_struc=trim(new_spam%structure)//reg_struc)
             else
                ! jcap: If we have more than one region, this becomes
                ! more complicated, as we have to calculate the union
                ! of several sparse structures to get the correct
                ! product structure. E.g. if C=A*B, C_11 = A_11 * B_11
                ! + A_12 * B_21, so we need the union of the
                ! structures of A_11 * B_11 and A_12 * B_21

                ! jcap: First, create temporary structures from the
                ! products of sub-region matrices, and immediately
                ! destroy them
                call sparse_create(tmp_spam, spama%m(irow,1), spamb%m(1,jcol), &
                     iscmplx=loc_cmplx, rlib=tmp_lib1)
                call sparse_destroy(tmp_spam)

                do ksub=1,new_spam%mrows-1

                   call sparse_create(tmp_spam,spama%m(irow,ksub+1),spamb%m(ksub+1,jcol), &
                        iscmplx=loc_cmplx,rlib=tmp_lib2)
                   call sparse_destroy(tmp_spam)

                   ! jcap: Now, calculate the union of these structures
                   call sparse_count_union(tmp_lib1,tmp_lib2,nze,nzb,&
                        my_nze,my_nzb,seg_nze,seg_nzb)

                   ! jcap: if we haven't finished yet, use a temporary
                   ! sparse structure
                   write(ksub_str,'(i1)') ksub
                   tmp_struc='TMP'//ksub_str//trim(new_spam%structure)//reg_struc
                   ! jcap: if we are on the last cycle round, use the
                   ! correct structure
                   if (ksub.eq.new_spam%mrows-1) &
                        tmp_struc=trim(new_spam%structure)//reg_struc

                   call sparse_index_union(tmp_lib1,tmp_lib2,tmp_struc,&
                        nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb,rlib=new_lib)

                   ! jcap: reassign the library numbers
                   tmp_lib1=new_lib

                end do

                ! jcap: Finally, create the sparse matrix using the
                ! new structure
                new_spam%m(isub,jsub)%structure=trim(new_spam%structure)//reg_struc
                call sparse_create(new_spam%m(isub,jsub), iscmplx=loc_cmplx, &
                     rlib=fake_rlib(isub,jsub))

             end if

          end if
       end do
    end do

    ! rc2013: assign the SPAM3 pointer
    new_spam%p => new_spam%m(1,1)

    ! rc2013: sort out libs
    if(present(rlib)) rlib = fake_rlib
    deallocate(fake_rlib, stat=ierr)
    call utils_dealloc_check('sparse_embed_alloc','fake_rlib',ierr)

  contains

    subroutine internal_choose_structure

      implicit none

      ! rc2013: set the value of row iterator
      if(present(arow)) then
         irow=arow
      else
         irow=isub
      endif
      ! rc2013: now repeat with column data
      if(present(bcol)) then
         jcol=bcol
      else
         jcol=jsub
      endif
      ! These are only necessary for extracting sub-blocks of matrix products
      ! Actual value doesn't matter: as long as jrow=icol should be fine
      ! jcap: this does actually seem to matter, at least in systems
      ! where the regions don't overlap. jrow=isub should be correct.
      if(present(brow)) then
         jrow=brow
      else
         jrow=isub
      endif
      if(present(acol)) then
         icol=acol
      else
         icol = jrow
      endif
      write(isub_str,'(i1)') irow
      write(jsub_str,'(i1)') jcol

      ! rc2013: now set the matrix block structure codes
      if(present(rstruc)) then
         write(isub_str,'(i1)') rstruc
      else
         write(isub_str,'(i1)') irow
      endif
      if(present(cstruc)) then
         write(jsub_str,'(i1)') cstruc
      else
         write(jsub_str,'(i1)') jcol
      endif

      ! rc2013: finally write out the structure we need
      ! rc2013: if this is a 1X1 matrix forget the region code
      if(pub_ngwf_regions_ngroups == 1) then
         reg_struc = ''
      else
         reg_struc = isub_str//jsub_str
      end if


    end subroutine internal_choose_structure

    subroutine internal_check_args

      implicit none

      ! rc2013: check that the arguments which have been provided are sane
      ! These combinations will crash!
      if (present(spamb) .and. (.not. present(spama))) &
           call utils_abort('Error in sparse_embed_create: wrong combination &
           &of optional arguments (spamb requires spama).')
      if (present(spamb) .and. (present(trans))) &
           call utils_abort('Error in sparse_embed_create: wrong combination &
           &of optional arguments (trans cannot be used with spamb ).')
      if (present(trans) .and. (.not. present(spama))) &
           call utils_abort('Error in sparse_embed_create: wrong combination &
           &of optional arguments (trans requires spama ).')

      ! These will just throw warnings
      !if( present(spama) .and. (present(mrows) .or. present(ncols)) ) then
      !     write(stdout,*)('Warning in sparse_embed_create: spama specified &
      !     along with mrows/ncols. The size of spama will be ignored when &
      !     new_spam is allocated. If this is what you intended, carry on. &
      !     If you actually wanted to create a subsystem matrix, consider &
      !     specifying which row (arow) or column (bcol) you want to use.')
      !end if
      !if( present(spama) .or. present(spamb) ) then
      !   if(present(structure)) write(stdout,*) 'WARNING: a matrix structure &
      !        & has been provided as input to sparse_embed_create, but &
      !        & spama and/or spamb are also present. The structure codes of &
      !        & these matrices will not be used when creating new_spam. &
      !        & If you want to specify the size of your matrix, &
      !        & consider providing mrows and ncols instead.'
      !end if

    end subroutine internal_check_args

  end subroutine sparse_embed_create


  !============================================================================!
  ! This subroutine destroys an embedded sparse matrix structure, freeing up   !
  ! the memory allocated to internal arrays. Sub-blocks are passed on to       !
  ! sparse_destroy.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   old_spam (inout) : The sparse matrix to be destroyed.                    !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 30/05/2017.                                    !
  !============================================================================!

  subroutine sparse_embed_destroy(old_spam,arrname)

    use sparse, only: sparse_destroy
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: old_spam
    character(len=*), optional, intent(in) :: arrname

    ! Local variables
    integer :: isub, jsub, mrows, ncols, ierr

    ! rc2013: size of embedding structure
    mrows = old_spam%mrows
    ncols = old_spam%ncols

    ! Deallocate memory assigned to internal structures
    do jsub=1,ncols
       do isub=1,mrows
          call sparse_destroy(old_spam%m(isub,jsub))
       end do
    end do

    ! Deallocate pointer
    nullify(old_spam%p)

    ! Deallocate the embedding array
    deallocate(old_spam%m, stat=ierr)
    ! rc2013: make sure that tag matches the one in sparse_embed_alloc
    if(present(arrname)) then
      call utils_dealloc_check('sparse_embed_alloc', arrname, ierr)
    else
      call utils_dealloc_check('sparse_embed_alloc', 'new_spam%m', ierr)
    end if

    old_spam%mrows = 0
    old_spam%ncols = 0

  end subroutine sparse_embed_destroy


  !============================================================================!
  ! This subroutine copies an embedded sparse matrix structure into a new      !
  ! matrix. Each subblock is passed on to sparse_copy.                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest (inout) : The destination sparse matrix.                            !
  !   src  (input) : The source sparse matrix.                                 !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 08/06/2017.                                    !
  ! Added cmplx_to_real pass-through - Jolyon Aarons, 05/2020.                 !
  !============================================================================!

  subroutine sparse_embed_copy(dest,src,cmplx_to_real)

    use sparse, only: sparse_copy
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: dest
    type(SPAM3_EMBED), intent(in)    :: src
    logical, optional, intent(in)    :: cmplx_to_real

    ! Local variables
    integer :: isub, jsub, mrows, ncols

    ! rc2013: check that the matrices match
    call utils_assert(src%mrows == dest%mrows, &
         'Error in sparse_embed_copy: rows(src) =/= rows(dest).')
    call utils_assert(src%ncols == dest%ncols, &
         'Error in sparse_embed_copy: cols(src) =/= cols(dest).')

    ! rc2013: size of embedding structure
    mrows = src%mrows
    ncols = src%ncols

    ! rc2013: pass each block to sparse_copy
    if(present(cmplx_to_real)) then
       do jsub=1,ncols
          do isub=1,mrows
             call sparse_copy(dest%m(isub,jsub), src%m(isub,jsub), cmplx_to_real)
          end do
       end do
    else
       do jsub=1,ncols
          do isub=1,mrows
             call sparse_copy(dest%m(isub,jsub), src%m(isub,jsub))
          end do
       end do
    end if

  end subroutine sparse_embed_copy

  !============================================================================!
  ! This function checks if a complex sparse matrix is hermitian, by computing !
  ! the RMS of the difference mat_dagger-mat                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The sparse complex matrix to check                      !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, October 2016                                      !
  ! Adapted for SPAM3_EMBED by Robert Charlton, 12/09/2018.                    !
  !============================================================================!

  real(kind=DP) function sparse_embed_check_hermitian(mat)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in)       :: mat           ! input complex matrix

    ! Local variables
    type(SPAM3_EMBED)                   :: mat_buffer    ! mat_dagger - mat

    ! agrecocmplx: use only with complex matrices
    call utils_assert(mat%iscmplx, 'ERROR in routine &
         &sparse_embed_check_hermitian: matrix must be complex')

    ! agrecocmplx: compute transpose conjugate of mat
    call sparse_embed_create(mat_buffer, mat)
    call sparse_embed_transpose(mat_buffer, mat)

    ! agrecocmplx: compute mat_dagger - mat
    call sparse_embed_axpy(mat_buffer, mat, -1.0_DP)

    ! agrecocmplx: compute RMS of mat_dagger - mat
    sparse_embed_check_hermitian = sparse_embed_rms_element(mat_buffer)

    call sparse_embed_destroy(mat_buffer)

  end function sparse_embed_check_hermitian



  !==========================================================================!
  ! This subroutine calculates the product of two sparse matrices with       !
  ! embedding structures.                                                    !
  !     C_ij = A_ik . B_kj                                                   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   cmat  (inout) : The sparse matrix C.                                   !
  !   amat  (input) : The sparse matrix A.                                   !
  !   bmat  (input) : The sparse matrix B.                                   !
  ! allow_mix_type  : Optional parameter to pass to sparse_product.          !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 30/05/2017.                                  !
  !==========================================================================!

  subroutine sparse_embed_product_array_array(cmat,amat,bmat,allow_mix_types)

    use sparse, only: sparse_product, sparse_axpy
    use utils, only: utils_assert

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    type(SPAM3_EMBED), intent(inout) :: cmat
    type(SPAM3_EMBED), intent(in)    :: amat
    type(SPAM3_EMBED), intent(in)    :: bmat
    logical, intent(in), optional :: allow_mix_types

    ! Local Variables
    type(SPAM3_EMBED) :: tmp_mat
    integer :: nsub, mrows, ncols
    integer :: isub, jsub, ksub
    integer :: ierr

    ! rc2013: check that the array sizes match
    call utils_assert(cmat%mrows == amat%mrows, &
         'Error in sparse_embed_product: rows(C) =/= rows(A).')
    call utils_assert(cmat%ncols == bmat%ncols, &
         'Error in sparse_embed_product: cols(C) =/= cols(B).')
    call utils_assert(amat%ncols == bmat%mrows, &
         'Error in sparse_embed_product: rows(B) =/= cols(A).')

    mrows = amat%mrows
    ncols = bmat%ncols
    nsub  = amat%ncols ! A_cols = B_rows

    ! rc2013: build tmp_mat using structures of cmat
    call sparse_embed_create(tmp_mat,cmat,iscmplx=cmat%iscmplx)

    ! rc2013: dump matrix info to output if requested
    if(pub_embed_debug) then
       call sparse_embed_diagnose(amat, 'amat')
       call sparse_embed_diagnose(bmat, 'bmat')
       call sparse_embed_diagnose(cmat, 'cmat')
    endif

    ! Do the multiplication
    do jsub=1,ncols
       do isub=1,mrows
          ! rc2013: sum up all contributions from subsystem products
          do ksub=1,nsub
             call sparse_product(cmat%m(isub,jsub),amat%m(isub,ksub), &
                  bmat%m(ksub,jsub),allow_mix_types)
             call sparse_axpy(tmp_mat%m(isub,jsub),cmat%m(isub,jsub),1.0_DP)
          end do
       end do
    end do
    ! rc2013: copy tmp_mat back to cmat
    call sparse_embed_copy(cmat,tmp_mat)

    ! rc2013: deallocate temporary matrices
    call sparse_embed_destroy(tmp_mat)

  end subroutine sparse_embed_product_array_array

  !----------------------------------------------------------------------------!

  subroutine sparse_embed_product_array_matrix(cmat,amat,bmat,row, col)

    use sparse, only: sparse_product
    use utils, only: utils_assert

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    type(SPAM3_EMBED), intent(inout) :: cmat
    type(SPAM3_EMBED), intent(in)    :: amat
    type(SPAM3),       intent(in)    :: bmat
    integer, optional, intent(in)    :: row
    integer, optional, intent(in)    :: col

    ! Local variables
    integer :: isub

    if(pub_embed_debug) then
       call sparse_embed_diagnose(amat, 'amat')
       call sparse_embed_diagnose(cmat, 'cmat')
    endif

    ! Do the multiplication
    if(present(row) .and. present(col)) then
       call sparse_product(cmat%m(row,col),amat%m(row,col),bmat)
    else if(present(row)) then
       ! rc2013: check that the array sizes match
       call utils_assert(cmat%ncols == amat%ncols, &
            'Error in sparse_embed_product_array_matrix: cols(C) =/= cols(B).')
       do isub=1,cmat%ncols
          call sparse_product(cmat%m(row,isub),amat%m(row,isub),bmat)
       end do
    else if(present(col)) then
       ! rc2013: check that the array sizes match
       call utils_assert(cmat%mrows == amat%mrows, &
            'Error in sparse_embed_product_array_matrix: rows(C) =/= rows(A).')
       do isub=1,cmat%mrows
          call sparse_product(cmat%m(isub,col),amat%m(isub,col),bmat)
       end do
    else
       call utils_assert(amat%ncols == 1, &
            'Error in sparse_embed_product_array_matrix: cols(A) =/= 1.')
       call utils_assert(cmat%ncols == 1, &
            'Error in sparse_embed_product_array_matrix: cols(C) =/= 1.')
       call utils_assert(cmat%mrows == amat%mrows, &
            'Error in sparse_embed_product_array_matrix: rows(C) =/= rows(A).')
       do isub=1,cmat%mrows
          call sparse_product(cmat%m(isub,1),amat%m(isub,1),bmat)
       end do
    endif

  end subroutine sparse_embed_product_array_matrix

  !----------------------------------------------------------------------------!

  ! rc2013: why is this laid out differently from array_matrix???
  subroutine sparse_embed_product_matrix_array(cmat,amat,bmat,row,col)

    use sparse, only: sparse_product
    use utils, only: utils_assert

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    type(SPAM3_EMBED), intent(inout) :: cmat
    type(SPAM3),       intent(in)    :: amat
    type(SPAM3_EMBED), intent(in)    :: bmat
    integer, optional, intent(in)    :: row
    integer, optional, intent(in)    :: col

    ! Local variables
    integer :: isub

    if(pub_embed_debug) then
       call sparse_embed_diagnose(bmat, 'bmat')
       call sparse_embed_diagnose(cmat, 'cmat')
    endif

    ! Do the multiplication
    if(present(row) .and. present(col)) then
       call sparse_product(cmat%m(row,col),amat,bmat%m(row,col))
    else if(present(row)) then
       ! rc2013: check that the array sizes match
       call utils_assert(cmat%ncols == bmat%ncols, &
            'Error in sparse_embed_product_matrix_array: cols(C) =/= cols(B).')
       do isub=1,cmat%ncols
          call sparse_product(cmat%m(row,isub),amat,bmat%m(row,isub))
       end do
    else if(present(col)) then
       ! rc2013: check that the array sizes match
       call utils_assert(cmat%mrows == bmat%mrows, &
            'Error in sparse_embed_product_matrix_array: rows(C) =/= rows(A).')
       do isub=1,cmat%mrows
          call sparse_product(cmat%m(isub,col),amat,bmat%m(isub,col))
       end do
    else
       call utils_assert(bmat%mrows == 1, &
            'Error in sparse_embed_product_matrix_array: rows(B) =/= 1.')
       call utils_assert(cmat%mrows == 1, &
            'Error in sparse_embed_product_matrix_array: rows(C) =/= 1.')
       call utils_assert(cmat%ncols == bmat%ncols, &
            'Error in sparse_embed_product_matrix_array: cols(C) =/= cols(B).')
       do isub=1,cmat%ncols
          call sparse_product(cmat%m(1,isub),amat,bmat%m(1,isub))
       enddo
    endif

  end subroutine sparse_embed_product_matrix_array

  !----------------------------------------------------------------------------!

  subroutine sparse_embed_product_matrix_matrix(cmat,amat,bmat)

    use sparse, only: sparse_product
    use utils, only: utils_assert

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    type(SPAM3_EMBED), intent(inout) :: cmat
    type(SPAM3), intent(in) :: amat
    type(SPAM3), intent(in) :: bmat

    ! rc2013: check that the array sizes match
    call utils_assert(cmat%mrows == 1, 'Error in &
         &sparse_embed_product: rows(C) =/= 1.')
    call utils_assert(cmat%ncols == 1, 'Error in &
         &sparse_embed_product: cols(C) =/= 1.')

    call sparse_product(cmat%m(1,1), amat, bmat)

  end subroutine sparse_embed_product_matrix_matrix


  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of an         !
  ! embedding matrix. Diagonal blocks are passed to sparse_scale.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 30/05/2017.                                    !
  !============================================================================!

  subroutine sparse_embed_scale_real(mat,alpha,beta)

    use parallel_strategy, only: PARAL_INFO
    use sparse, only: sparse_scale, sparse_get_par

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)    :: mat    ! The matrix to be operated on
    real(kind=DP), intent(in)           :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local variables
    integer :: isub, jsub
    integer :: mrows, ncols
    type(PARAL_INFO), pointer :: row_par, col_par

    ! rc2013: get the embedding dimensions
    mrows = mat%mrows
    ncols = mat%ncols

    ! rc2013: cycle over all subsystem blocks
    do jsub=1,ncols
       do isub=1,mrows
          ! rc2013: parallel strategies match for diagonal blocks
          call sparse_get_par(row_par, mat%m(isub,jsub), 'R')
          call sparse_get_par(col_par, mat%m(isub,jsub), 'C')
          if ( row_par%par_index == col_par%par_index ) then
             call sparse_scale(mat%m(isub,jsub), alpha, beta)
          else
             ! rc2013: the matrix is not square, so we can't apply beta
             mat%m(isub,jsub)%dmtx = alpha * mat%m(isub,jsub)%dmtx
          end if
       end do
    end do

  end subroutine sparse_embed_scale_real

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of an         !
  ! embedding matrix. Diagonal blocks are passed to sparse_scale.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 30/05/2017.                                    !
  !============================================================================!

  subroutine sparse_embed_scale_complex(mat,alpha,beta)

    use parallel_strategy, only: PARAL_INFO
    use sparse, only: sparse_scale, sparse_get_par

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)    :: mat    ! The matrix to be operated on
    complex(kind=DP), intent(in)        :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local variables
    integer :: isub, jsub
    integer :: mrows, ncols
    type(PARAL_INFO), pointer :: row_par, col_par

    ! rc2013: get the embedding dimensions
    mrows = mat%mrows
    ncols = mat%ncols

    ! rc2013: cycle over all subsystem blocks
    do jsub=1,ncols
       do isub=1,mrows
          ! rc2013: parallel strategies match for diagonal blocks
          call sparse_get_par(row_par, mat%m(isub,jsub), 'R')
          call sparse_get_par(col_par, mat%m(isub,jsub), 'C')
          if ( row_par%par_index == col_par%par_index ) then
             call sparse_scale(mat%m(isub,jsub), alpha, beta)
          else
             ! rc2013: the matrix is not square, so we can't apply beta
             mat%m(isub,jsub)%zmtx = alpha * mat%m(isub,jsub)%zmtx
          end if
       end do
    end do

  end subroutine sparse_embed_scale_complex

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of an         !
  ! embedding matrix. Diagonal blocks are passed to sparse_scale.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional complex shift parameter                    !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 30/05/2017.                                    !
  !============================================================================!

  subroutine sparse_embed_scale_complex2(mat,alpha,beta)

    use parallel_strategy, only: PARAL_INFO
    use sparse, only: sparse_scale, sparse_get_par

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)    :: mat    ! The matrix to be operated on
    complex(kind=DP), intent(in)           :: alpha  ! The scaling parameter
    complex(kind=DP), intent(in) :: beta   ! The shift parameter

    ! Local variables
    integer :: isub, jsub
    integer :: mrows, ncols
    type(PARAL_INFO), pointer :: row_par, col_par

    ! rc2013: get the embedding dimensions
    mrows = mat%mrows
    ncols = mat%ncols

    ! rc2013: cycle over all subsystem blocks
    do jsub=1,ncols
       do isub=1,mrows
          ! rc2013: parallel strategies match for diagonal blocks
          call sparse_get_par(row_par, mat%m(isub,jsub), 'R')
          call sparse_get_par(col_par, mat%m(isub,jsub), 'C')
          if ( row_par%par_index == col_par%par_index ) then
             call sparse_scale(mat%m(isub,jsub), alpha, beta)
          else
             ! rc2013: the matrix is not square, so we can't apply beta
             mat%m(isub,jsub)%dmtx = alpha * mat%m(isub,jsub)%dmtx
          end if
       end do
    end do

  end subroutine sparse_embed_scale_complex2

  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  ! by passing the system subblocks down to sparse_axpy.                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat(mdl%nsub,mdl%nsub)  (inout) : The sparse matrix y                   !
  !   xmat(mdl%nsub,mdl%nsub)  (input) : The sparse matrix x                   !
  !   alpha                    (input) : The (real) parameter alpha            !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 20/07/17.                                      !
  !============================================================================!

  subroutine sparse_embed_axpy_real(ymat,xmat,alpha)

    use sparse, only: sparse_axpy
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: ymat   ! The sparse matrix y
    type(SPAM3_EMBED), intent(in)    :: xmat   ! The sparse matrix x
    real(kind=DP), intent(in)        :: alpha  ! The parameter alpha

    ! Local variables
    integer :: isub, jsub
    integer :: mrows, ncols

    ! rc2013: get the embedding dimensions
    mrows = ymat%mrows
    ncols = ymat%ncols

    ! rc2013: check that the array sizes match
    call utils_assert(ymat%mrows == xmat%mrows, 'Error in &
         &sparse_axpy_embed: ymat and xmat are different sizes.')
    call utils_assert(ymat%ncols == xmat%ncols, 'Error in &
         &sparse_axpy_embed: ymat and xmat are different sizes.')

    ! rc2013: loop over subsystems
    do jsub=1,ncols
       do isub=1,mrows
          call sparse_axpy(ymat%m(isub,jsub),xmat%m(isub,jsub),alpha)
       end do
    end do

  end subroutine sparse_embed_axpy_real

  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  ! by passing the system subblocks down to sparse_axpy.                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat(mdl%nsub,mdl%nsub)  (inout) : The sparse matrix y                   !
  !   xmat(mdl%nsub,mdl%nsub)  (input) : The sparse matrix x                   !
  !   alpha                    (input) : The (real) parameter alpha            !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 20/07/17.                                      !
  ! Generalised for complex alpha by Joseph Prentice, August 2018              !
  !============================================================================!

  subroutine sparse_embed_axpy_complex(ymat,xmat,alpha)

    use sparse, only: sparse_axpy
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: ymat   ! The sparse matrix y
    type(SPAM3_EMBED), intent(in)    :: xmat   ! The sparse matrix x
    complex(kind=DP), intent(in)        :: alpha  ! The parameter alpha

    ! Local variables
    integer :: isub, jsub
    integer :: mrows, ncols

    ! rc2013: get the embedding dimensions
    mrows = ymat%mrows
    ncols = ymat%ncols

    ! rc2013: check that the array sizes match
    call utils_assert(ymat%mrows == xmat%mrows, 'Error in &
         &sparse_axpy_embed: ymat and xmat are different sizes.')
    call utils_assert(ymat%ncols == xmat%ncols, 'Error in &
         &sparse_axpy_embed: ymat and xmat are different sizes.')

    ! rc2013: loop over subsystems
    do jsub=1,ncols
       do isub=1,mrows
          call sparse_axpy(ymat%m(isub,jsub),xmat%m(isub,jsub),alpha)
       end do
    end do

  end subroutine sparse_embed_axpy_complex

  !==========================================================================!
  ! This subroutine forms the transpose of an embedding matrix.              !
  !--------------------------------------------------------------------------!
  !   dest  (output) : The resulting transpose                               !
  !   src   (input)  : The matrix to be transposed                           !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 20/07/2017.                                  !
  !==========================================================================!

  subroutine sparse_embed_transpose(dest,src)

    use sparse, only: sparse_transpose
    use utils, only: utils_assert

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    type(SPAM3_EMBED), intent(inout) :: dest
    type(SPAM3_EMBED), intent(in)    :: src

    ! Local Variables
    integer :: mrows, ncols
    integer :: isub, jsub

    ! rc2013: check that the array sizes match
    call utils_assert(dest%mrows == src%mrows, 'Error in &
         &sparse_embed_transpose: src and dest are different sizes.')
    call utils_assert(dest%ncols == src%ncols, 'Error in &
         &sparse_embed_transpose: src and dest are different sizes.')

    ! rc2013: number of subsystems. Generalise for oddly shaped matrices.
    mrows = dest%mrows
    ncols = dest%ncols

    ! Do the multiplication
    do jsub=1,ncols
       do isub=1,mrows
          ! rc2013: straightforward transpose for diagonal blocks
          if (isub == jsub) then
             call sparse_transpose(dest%m(isub,jsub),src%m(isub,jsub))
          else
             ! rc2013: swap the off-diagonal blocks around from src to dest
             call sparse_transpose(dest%m(jsub,isub),src%m(isub,jsub))
          end if
       end do
    end do

  end subroutine sparse_embed_transpose



  !==========================================================================!
  ! This subroutine calculates the trace of a sparse matrix or the product of!
  ! two sparse matrices. Sub-blocks are passed down to sparse_trace.         !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   amat  (input) : The embedding sparse matrix.                           !
  !   bmat  (input) : The (optional) embedding sparse matrix B.              !
  !   opA   (input) : An optional input of some kind...                      !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 23/07/2017.                                  !
  ! Converted into subroutine by Joseph Prentice, 29/01/2022                 !
  !==========================================================================!

  subroutine sparse_embed_trace(trace,amat,bmat,opA)

    use sparse, only: sparse_trace
    use utils, only: utils_assert

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    real(kind=DP), intent(out) :: trace
    type(SPAM3_EMBED), intent(in) :: amat
    type(SPAM3_EMBED), optional, intent(in) :: bmat
    character(len=1), optional, intent(in) :: opA  ! operation to be applied to A
    ! rc2013: need to check if there's anything peculiar about opA...

    ! Local Variables
    integer :: mrows, ncols, isub, jsub
    real(kind=DP) :: tr_sum

    ! rc2013: check that the array sizes match
    if(present(bmat)) then
       call utils_assert(amat%mrows == bmat%ncols, 'Error in &
            &sparse_embed_trace: rows(A) =/= cols(B).')
       call utils_assert(amat%ncols == bmat%mrows, 'Error in &
            &sparse_embed_trace: cols(A) =/= rows(B).')
    else
       ! rc2013: make sure that A is a square matrix of SPAM3s
       call utils_assert(amat%mrows == amat%ncols, 'Error in &
            &sparse_embed_trace: A is not square.')
    end if

    ! rc2013: number of subsystems
    mrows = amat%mrows
    ncols = amat%ncols
    tr_sum = 0.0_DP

    ! Do the trace multiplication -- here rows needed in outer loop
    do isub=1,mrows
       ! rc2013: calculate the trace of each block product individually
       if(present(bmat)) then
          ! rc2013: additional loop for matrix product
          do jsub=1,ncols
             tr_sum = tr_sum + &
                  sparse_trace(amat%m(isub,jsub),bmat%m(jsub,isub), opA=opA)
          end do
       else
          tr_sum = tr_sum + sparse_trace(amat%m(isub,isub), opA=opA)
       end if
    end do

    ! rc2013: return the total trace
    trace = tr_sum

  end subroutine sparse_embed_trace


  !============================================================================!
  ! This function returns the number of non-zero elements in an embedded       !
  ! sparse matrix in a floating-point real (so as to avoid integer overflows). !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 25/07/2017.                                    !
  ! Based on sparse_num_element_mat by Nicholas Hine.                          !
  !============================================================================!

  function sparse_embed_num_element(mat) result(num_sum)

    use  sparse, only: sparse_num_element

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat   ! The embedding sparse matrix

    ! Local variables
    integer :: isub, jsub, mrows, ncols

    ! Result
    real(kind=DP) :: num_sum

    ! rc2013: number of subsystems
    mrows = mat%mrows
    ncols = mat%ncols
    num_sum = 0.0_DP

    ! rc2013: Add up the number of elements in each block individually
    do jsub=1,ncols
       do isub=1,mrows
          num_sum = num_sum + sparse_num_element(mat%m(isub,jsub))
       enddo
    enddo

  end function sparse_embed_num_element

  !============================================================================!
  ! This function returns the number of rows in a sparse matrix by consulting  !
  ! the relevant library entry.                                                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 28/07/2017, based on sparse_num_rows by        !
  ! Nicholas Hine.                                                             !
  !============================================================================!

  function sparse_embed_num_rows(mat) result(num_rows)

    use sparse, only: sparse_num_rows

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat   ! The sparse matrix

    ! Local variables
    integer :: isub

    ! Result
    integer :: num_rows

    ! Very simple... just sum the rows from block-diagonal SPAM3's
    num_rows = 0.0_DP
    do isub=1,mat%mrows
       num_rows = num_rows + sparse_num_rows(mat%m(isub,isub))
    enddo

  end function sparse_embed_num_rows

  !============================================================================!
  ! This function returns the number of columns in an embedding sparse matrix  !
  ! by consulting  the relevant library entries.                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 28/07/2017, based on sparse_num_cols by        !
  ! Nicholas Hine.                                                             !
  !============================================================================!

  function sparse_embed_num_cols(mat) result(num_cols)

    use sparse, only: sparse_num_cols

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat   ! The sparse matrix

    ! Local variables
    integer :: isub

    ! Result
    real(kind=DP) :: num_cols

    ! Very simple...
    num_cols = 0.0_DP
    do isub=1,mat%ncols
       num_cols = num_cols + sparse_num_cols(mat%m(isub,isub))
    enddo

  end function sparse_embed_num_cols

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine allocates a new array of sparse matrices and assigns a     !
  ! structure to them.                                                         !
  !----------------------------------------------------------------------------!
  ! Required arguments:                                                        !
  !   new_spam (output) : The array of sparse matrices to be created           !
  !----------------------------------------------------------------------------!
  ! Optional arguments:
  !   spama    (input)  : An array of sparse matrices whose structure is used  !
  !   spamb    (input)  : A second array of sparse matrices whose structure    !
  !                       will be used.                                        !
  !   spammata, spammatb (input) : non-array alternatives for spama and spamb  !
  !   n_spins   (input) : Force this spin dimension in the array.              !
  !   n_kpoints (input) : Force this number of kpoints in the array.           !
  !   structure (input) : Force this structure in the new array.               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 02/08/2017.                                    !
  ! Plagiarised from sparse_array_alloc. For now this only handles SPAM3s, not !
  ! SPAM3_ARRAYs.                                                              !
  !============================================================================!

  subroutine sparse_embed_alloc(new_spam, mrows, ncols, spama, spamb, &
        rlib, row, col, arrname)

     use utils, only: utils_alloc_check, utils_assert
     implicit none

     ! Arguments
     type(SPAM3_EMBED), intent(inout) :: new_spam     ! Sparse matrix array to be created
     integer, intent(in)              :: mrows        ! Optional number of rows.
     integer, intent(in)              :: ncols        ! Optional number of columns.
     type(SPAM3_EMBED), optional, intent(in) :: spama  ! First (optional) matrix array to use
     type(SPAM3_EMBED), optional, intent(in) :: spamb  ! Other (optional) matrix array to use
     !character(len=*), optional, intent(in) :: structure  ! Optional structure string.
     integer, intent(inout), allocatable :: rlib(:,:)
     integer, optional, intent(in) :: row
     integer, optional, intent(in) :: col
     character(len=*), optional, intent(in) :: arrname

     ! Local variables
     integer :: ierr        ! Error flag

     ! Set the number of rows and columns in the new matrix
     new_spam%mrows = mrows
     new_spam%ncols = ncols
     if ( (present(spama)) .and. (present(spamb)) ) then
        ! rc2013: size determined by product of A and B
        ! --> but only if we don't specify what row/col we want!
        if(.not.(present(row) .or. present(col))) &
             call utils_assert(spama%ncols == spamb%mrows, &
             'Error in sparse_embed_alloc: cols(A) =/= rows(B)')
        new_spam%mrows = spama%mrows
        new_spam%ncols = spamb%ncols
     else if ( present(spama) ) then
        new_spam%mrows = spama%mrows
        new_spam%ncols = spama%ncols
     end if

     ! rc2013: overwrite the above parameters if row/col are specified
     if(present(row)) new_spam%mrows=1
     if(present(col)) new_spam%ncols=1

     ! Set the structure
     if ( (present(spama)) .and. (present(spamb)) ) then
        ! rc2013: structure determined by product of A and B
        new_spam%structure = trim(spama%structure)//trim(spamb%structure)
     else if ( present(spama) ) then
        new_spam%structure = trim(spama%structure)
     end if
     ! Otherwise use the structure already in the matrix

     ! Now let's allocate the matrix array
     allocate(new_spam%m(new_spam%mrows, new_spam%ncols), stat=ierr)
     ! Tag differs from variable so that it matches the tag in
     ! sparse_array_destroy
     if(present(arrname)) then
       call utils_alloc_check('sparse_embed_alloc', arrname, ierr)
     else
       call utils_alloc_check('sparse_embed_alloc', 'new_spam%m', ierr)
     end if

     ! rc2013: allocate empty array to avoid segmentation faults
     allocate(rlib(new_spam%mrows,new_spam%ncols),stat=ierr)
     call utils_alloc_check('sparse_embed_alloc','fake_rlib',ierr)

  end subroutine sparse_embed_alloc

  !============================================================================!
  ! This function returns true if every element in the matrix is nonzero, or   !
  ! false otherwise.                                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array   (input)  : The array to be assessed                              !
  !----------------------------------------------------------------------------!
  ! Written by J.M. Escartin, March 2015.                                      !
  ! Adapted for embedding by Robert Charlton, 15/08/2017.                      !
  !============================================================================!

  function sparse_embed_is_dense(array) result(is_dense)

     use sparse, only: sparse_is_dense

     implicit none

     ! Argument
     type(SPAM3_EMBED), intent(in) :: array ! The array of sparse matrices

     ! Result
     logical :: is_dense
     integer :: isub,jsub

     ! rc2013: could do with a more rigorous check
     ! rc2013: HACK!
     do jsub=1,array%ncols
        do isub=1,array%mrows
           is_dense = sparse_is_dense(array%m(1,1))
           if(.not. is_dense) return
        end do
     end do

  end function sparse_embed_is_dense

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine returns the structure code which corresponds to the        !
  ! transpose of the sparse matrix passed in. For symmetric structures this    !
  ! is the same as the structure of the matrix itself, but for rectangular     !
  ! matrices it may be different.                                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be transposed                             !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, March 2011                                       !
  ! Modified for embedding by Joseph Prentice, June 2018                       !
  !============================================================================!

  subroutine sparse_embed_transpose_structure(struc_trans,mat)

    use sparse, only: sparse_transpose_structure

    implicit none

    ! Arguments
    character(len=30), intent(out) :: struc_trans
    type(SPAM3_EMBED), intent(in) :: mat   ! The sparse matrix

    ! Local variable
    character(len=30) :: struc
    integer :: len_struc

    ! Get transpose strcture
    ! jcap: bit of a hack - get transpose structure of 1,1 sub block
    ! (which will always exist and always be square), strip off the
    ! embedding suffix and use this
    if (mat%mrows == 1 .and. mat%ncols == 1) then
       ! rc2013: if this is a non-embedding run then sparse_transpose will do
       call sparse_transpose_structure(struc_trans,mat%m(1,1))
    else
       call sparse_transpose_structure(struc,mat%m(1,1))
       len_struc=len_trim(struc)-2  ! -2 to remove embedding suffix
       struc_trans = struc(1:len_struc)
    end if

  end subroutine sparse_embed_transpose_structure

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine calculates the outer product matrix uv of two vectors. The !
  ! input matrix sparsity pattern is respected in the final answer, so lots of !
  ! nonzero terms of the outer product may get skipped. It is the user's       !
  ! responsibility to ensure that a suitable sparsity pattern is passed in     !
  ! that will not result in unwanted truncation.                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !  mat (inout) : The sparse matrix                                           !
  !  uvec (in)   : The real row vector u                                       !
  !  vvec (in)   : The real column vector v                                    !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, June 2018                                      !
  !============================================================================!

  subroutine sparse_embed_outer_product_real(mat,uvec,vvec)

    use sparse, only : sparse_num_rows,sparse_num_cols,sparse_outer_product
    use utils, only : utils_alloc_check, utils_dealloc_check

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: mat   ! The sparse matrix
    real(kind=DP), intent(in) :: uvec(:),vvec(:)

    ! Local variables
    integer :: isub,jsub,past_rows,past_cols,num_rows,num_cols
    integer :: ierr
    real(kind=DP), allocatable :: utmp(:),vtmp(:)

    ! jcap: Loop over regions
    past_rows=0
    do isub=1,mat%mrows
       num_rows=int(sparse_num_rows(mat%m(isub,1)))
       ! Create temporary vector, only including section needed for
       ! this embedding region
       allocate(utmp(num_rows),stat=ierr)
       call utils_alloc_check('sparse_embed_outer_product','utmp',ierr)
       utmp=uvec((past_rows+1):(past_rows+num_rows))
       past_cols=0
       do jsub=1,mat%ncols
          num_cols=int(sparse_num_cols(mat%m(isub,jsub)))
          ! Create temporary vector, only including section needed for
          ! this embedding region
          allocate(vtmp(num_cols),stat=ierr)
          call utils_alloc_check('sparse_embed_outer_product','vtmp',ierr)
          vtmp=vvec((past_cols+1):(past_cols+num_cols))
          call sparse_outer_product(mat%m(isub,jsub),utmp,vtmp)
          deallocate(vtmp,stat=ierr)
          call utils_dealloc_check('sparse_embed_outer_product','vtmp',ierr)
          past_cols=past_cols+num_cols
       end do
       deallocate(utmp,stat=ierr)
       call utils_dealloc_check('sparse_embed_outer_product','utmp',ierr)
       past_rows=past_rows+num_rows
    end do

  end subroutine sparse_embed_outer_product_real

!..............................................................................

  subroutine sparse_embed_outer_product_complex(mat,uvec,vvec)

    use sparse, only : sparse_num_rows,sparse_num_cols,sparse_outer_product
    use utils, only : utils_alloc_check, utils_dealloc_check

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: mat   ! The sparse matrix
    complex(kind=DP), intent(in) :: uvec(:),vvec(:)

    ! Local variables
    integer :: isub,jsub,past_rows,past_cols,num_rows,num_cols
    integer :: ierr
    complex(kind=DP), allocatable :: utmp(:),vtmp(:)

    ! jcap: Loop over regions
    past_rows=0
    do isub=1,mat%mrows
       num_rows=int(sparse_num_rows(mat%m(isub,1)))
       ! Create temporary vector, only including section needed for
       ! this embedding region
       allocate(utmp(num_rows),stat=ierr)
       call utils_alloc_check('sparse_embed_outer_product','utmp',ierr)
       utmp=uvec((past_rows+1):(past_rows+num_rows))
       past_cols=0
       do jsub=1,mat%ncols
          num_cols=int(sparse_num_cols(mat%m(isub,jsub)))
          ! Create temporary vector, only including section needed for
          ! this embedding region
          allocate(vtmp(num_cols),stat=ierr)
          call utils_alloc_check('sparse_embed_outer_product','vtmp',ierr)
          vtmp=vvec((past_cols+1):(past_cols+num_cols))
          call sparse_outer_product(mat%m(isub,jsub),utmp,vtmp)
          deallocate(vtmp,stat=ierr)
          call utils_dealloc_check('sparse_embed_outer_product','vtmp',ierr)
          past_cols=past_cols+num_cols
       end do
       deallocate(utmp,stat=ierr)
       call utils_dealloc_check('sparse_embed_outer_product','utmp',ierr)
       past_rows=past_rows+num_rows
    end do

  end subroutine sparse_embed_outer_product_complex


  !==========================================================================!
  !                HOTELLING ALGORITHM FOR MATRIX INVERSION                  !
  !==========================================================================!

  !==========================================================================!
  ! This subroutine initialises the inverse of a sparse matrix so that it is !
  ! ready to be passed to sparse_hotelling_invert and inverted using the     !
  ! Hotelling Algorithm.                                                     !
  ! Based on the theory in T. Ozaki, Phys. Rev. B., vol 64, page 195110.     !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   inv_mat  (input) : The inverse to initialise                           !
  !   mat      (input) : The matrix to be inverted                           !
  !--------------------------------------------------------------------------!
  ! Originally written by Chris-Kriton Skylaris on 31/7/2003.                !
  ! Modified to support SPAM 2 by Peter Haynes on 21/7/2004.                 !
  ! Adapted for SPAM3 matrices by Nicholas Hine, June 2009.                  !
  ! Moved to sparse_mod by Nicholas Hine, September 2009.                    !
  ! Rewritten for embedding by Robert Charlton, 02/06/2017.                  !
  ! Moved to sparse_embed_mod by Robert Charlton, 07/12/2017.                !
  !==========================================================================!

  subroutine sparse_embed_hotelling_init(inv_mat,mat)

    use comms, only: pub_on_root, comms_reduce, pub_my_proc_id
    use sparse, only: sparse_first_elem_on_proc, &
        sparse_num_rows, sparse_get_col, sparse_clr_col, sparse_last_elem_on_proc, &
        sparse_element_exists
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: inv_mat
    type(SPAM3_EMBED), intent(in)    :: mat

    ! Local Variables
    real(kind=DP) :: sigma, sum_col
    real(kind=DP), allocatable, dimension(:) :: drow
    complex(kind=DP), allocatable, dimension(:) :: zrow
    integer :: nn, row, col, first_col, num_cols, last_col
    integer :: ierr
    integer :: isub, jsub, mrows, ncols

    ! agrecocmplx: either both real or both complex
    call utils_assert((mat%iscmplx .eqv. inv_mat%iscmplx), &
         'Error in sparse_embed_hotelling_init: incompatible argument types.')

    ! rc2013: get matrix dimensions
    mrows = mat%mrows
    ncols = mat%ncols

    ! rc2013: loop over all subsystems
    sigma = 0.0_DP
    do jsub=1,ncols
       ! Loop over all cols of S on this processor
       first_col = sparse_first_elem_on_proc(pub_my_proc_id,mat%m(1,jsub),'C')
       last_col = sparse_last_elem_on_proc(pub_my_proc_id,mat%m(1,jsub),'C')
       do col = first_col, last_col
          sum_col = 0.0_DP
          do isub=1,mrows
             nn = sparse_num_rows(mat%m(isub,jsub))
             if (mat%m(isub,jsub)%iscmplx) then
                allocate(zrow(nn),stat=ierr)
                call utils_alloc_check('sparse_embed_hotelling_init','zrow',ierr)
                zrow = 0.0_DP
             else
                allocate(drow(nn),stat=ierr)
                call utils_alloc_check('sparse_embed_hotelling_init','drow',ierr)
                drow = 0.0_DP
             end if

             if (mat%m(isub,jsub)%iscmplx) then

                ! Get column of S
                call sparse_get_col(zrow,mat%m(isub,jsub),col)

                ! Sum column of S
                do row=1,nn
                   ! jd: Ignore 'inaccessible' elements in dense blocks, whose values
                   !     are undefined.
                   if(sparse_element_exists(mat%m(isub,jsub),row,col)) then
                      sum_col = sum_col + abs(zrow(row))
                   end if
                end do

                ! Clear column of S
                call sparse_clr_col(zrow,mat%m(isub,jsub),col)

             else

                ! Get column of S
                call sparse_get_col(drow,mat%m(isub,jsub),col)

                ! Sum column of S
                do row=1,nn
                   sum_col = sum_col + abs(drow(row))
                end do

                ! Clear column of S
                call sparse_clr_col(drow,mat%m(isub,jsub),col)

             end if

             ! Deallocate workspace
             if (mat%m(isub,jsub)%iscmplx) then
                deallocate(zrow,stat=ierr)
                call utils_dealloc_check('sparse_embed_hotelling_init','zrow',ierr)
             else
                deallocate(drow,stat=ierr)
                call utils_dealloc_check('sparse_embed_hotelling_init','drow',ierr)
             end if

          end do
          sigma = max(sigma,sum_col)
       end do
    end do
    call comms_reduce('MAX',sigma)

    call sparse_embed_copy(inv_mat,mat)

    if (sigma > epsilon(1.0_DP)) then
       sigma = 1.0_DP / (sigma*sigma)
    else
       sigma = 0.001_DP
       if (pub_on_root) write(stdout,'(a)') &
            'WARNING in sparse_embed_hotelling_init: zero overlap matrix'
    end if

    ! cks: scale by sigma
    call sparse_embed_scale(inv_mat,sigma)

  end subroutine sparse_embed_hotelling_init

  !==========================================================================!
  ! This subroutine applies successive quadratically convergent Hotelling    !
  ! iterations to improve an approximate inverse overlap matrix, based on    !
  ! the theory in the paper by T. Ozaki, Phys. Rev. B., vol 64, page 195110. !
  ! The iterations continue (unless they exceed the num_iter input           !
  ! parameter) until convergence to machine precision is reached - but       !
  ! taking into account limitations arising from the truncation of the       !
  ! inverse overlap matrix.                                                  !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   inv_mat  (input) : The inverse to calculate                            !
  !   mat      (input) : The matrix to be inverted                           !
  !--------------------------------------------------------------------------!
  ! Originally written by Chris-Kriton Skylaris on 31/7/2003                 !
  ! Modified to support SPAM 2 by Peter Haynes on 21/7/2004.                 !
  ! Modified to test for convergence by Chris-Kriton Skylaris on 4/10/2004.  !
  ! Minor modifications for parallel SPAM 2 by Peter Haynes                  !
  ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009.               !
  ! Moved to sparse_mod by Nicholas Hine, September 2009.                    !
  ! Adapted for embedding by Robert Charlton, 29/05/17.                      !
  ! Moved to sparse_embed_mod by Robert Charlton, 07/12/2017.                !
  !==========================================================================!

  subroutine sparse_embed_hotelling_invert(inv_mat,mat,show_output, &
       max_resid_converged, num_iter, final_max_resid, ireg)

    use comms, only: pub_on_root
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    type(SPAM3_EMBED), intent(inout) :: inv_mat
    type(SPAM3_EMBED), intent(in)    :: mat
    real(kind=DP), intent(in) :: max_resid_converged
    logical, intent(in) :: show_output
    integer, intent(in) :: num_iter
    real(kind=DP), intent(out), optional :: final_max_resid
    integer, intent(in), optional :: ireg

    ! Local Variables
    integer :: iter
    logical :: quit_early  ! pdh: flag
    real(kind=DP) :: max_resid  ! maximum value of residual
    real(kind=DP) :: frob_norm
    real(kind=DP) :: previous_frob_norm
    type(SPAM3_EMBED) :: mim, inv_tmp
    integer :: ierr

    ! Start timer
    call timer_clock('sparse_embed_hotelling_invert_embed',1)

    ! agrecocmplx: either both real or both complex
    call utils_assert((mat%iscmplx .eqv. inv_mat%iscmplx), &
         'Error in sparse_embed_hotelling_invert: incompatible argument types.')

    ! Initialisations
    previous_frob_norm = huge(1.0_DP)
    quit_early = .false.

    ! Allocate workspaces
    call sparse_embed_create(inv_tmp, inv_mat)
    call sparse_embed_create(mim, mat, inv_mat, rstruc=ireg, cstruc=ireg)

    ! Start of Hotelling iteration loop
    hotel_sparse_loop: do iter=1,num_iter

       ! MIM := M*(M^-1)
       call sparse_embed_product(mim,mat,inv_mat)

       ! cks: MIM := I-M*(M^-1)
       call sparse_embed_scale(mim,-1.0_DP,1.0_DP)

       ! cks: Calculate Frobenius norm of MIM -->0
       frob_norm = sparse_embed_rms_element(mim) * &
            sqrt(sparse_embed_num_element(mim))

       ! cks: Maximum element of residual I-M*M_n^-1
       max_resid = sparse_embed_max_abs_element(mim)

       if (show_output.and.pub_on_root) then
           write(stdout,'(t12,i5,tr5,e16.8,tr3,e16.8)') iter, frob_norm, &
                max_resid
       end if

       ! cks: Test for convergence due to machine precision
       ! cks: or inverse overlap truncation
       if (frob_norm >= previous_frob_norm .or. &
            max_resid <= max_resid_converged) then
          quit_early = .true.
          exit hotel_sparse_loop
       else
          previous_frob_norm = frob_norm
       endif

       ! cks: MIM := I + MIM := 2I -M*(M^-1)
       call sparse_embed_scale(mim, 1.0_DP, 1.0_DP)

       ! (M^-1) := (M^-1)*(2*\delta_{jk}-M*(M^-1))
       call sparse_embed_copy(inv_tmp, inv_mat)
       call sparse_embed_product(inv_mat,inv_tmp,mim)

    end do hotel_sparse_loop

    if (present(final_max_resid)) final_max_resid = max_resid

    if (pub_on_root .and. .not. quit_early) write(stdout,*) &
         'WARNING: max Hotelling iterations exceeded. Last norm=',frob_norm

    ! Deallocate workspaces
    call sparse_embed_destroy(inv_tmp)
    call sparse_embed_destroy(mim)

    ! Stop timer
    call timer_clock('sparse_embed_hotelling_invert_embed',2)

  end subroutine sparse_embed_hotelling_invert


  !==========================================================================!
  ! This subroutine sets a SPAM3_EMBED matrix to garbage.                    !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   amat  (input) : The embedding sparse matrix.                           !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 13/08/2018.                                  !
  !==========================================================================!

  subroutine sparse_embed_set_to_garbage(matrix)

    use sparse, only: sparse_set_to_garbage

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: matrix

    ! Local Variables
    integer :: isub, jsub

    ! Do the trace multiplication -- here rows needed in outer loop
    do jsub=1,matrix%ncols
       do isub=1,matrix%mrows
          call sparse_set_to_garbage(matrix%m(isub,jsub))
       end do
    end do

  end subroutine sparse_embed_set_to_garbage

  !============================================================================!
  ! This function checks if any elements of a matrix are NaN.                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input) : The sparse matrix to check                               !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, 30/04/2020.                                      !
  !============================================================================!

  function sparse_embed_any_isnan(mat) result(isnan)

    use comms, only: comms_reduce
    use sparse, only: sparse_any_isnan
    implicit none

    ! Arguments
    type(SPAM3_embed), intent(in) :: mat
    logical                       :: isnan

    ! Local Variables
    logical                       :: loc_isnan
    integer                       :: isub, jsub

    loc_isnan=.false.

    ! Check if local data contains NaN
    do jsub=1,mat%ncols
       do isub=1,mat%mrows
          loc_isnan=loc_isnan.or.sparse_any_isnan(mat%m(isub,jsub))
       end do
    end do

    ! Check if any proc has found a NaN
    call comms_reduce('OR',loc_isnan)

    isnan = loc_isnan

  end function sparse_embed_any_isnan


  !============================================================================!
  ! This function calculates the entrywise norm of the data in a given sparse  !
  ! matrix.                                                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input) : The sparse matrix to check                               !
  !   ord   (input) : order to raise entries to                                !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 19/09/2018.                                    !
  ! Based on SPAM3 version by Jolyon Aarons.                                   !
  !============================================================================!

  function sparse_embed_entrywise_norm(mat,ord) result(norm)

    use comms, only: comms_reduce, comms_bcast, pub_root_proc_id

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in)     :: mat
    integer,     intent(in)     :: ord

    ! Local variables
    real(kind=DP)               :: norm
    integer                     :: isub, jsub

    norm = 0.0_dp
    do jsub=1,mat%ncols
       do isub=1,mat%mrows
          if(mat%m(isub,jsub)%iscmplx) then
             norm = norm + sum(abs(mat%m(isub,jsub)%zmtx)**ord)
          else
             norm = norm + sum(abs(mat%m(isub,jsub)%dmtx)**ord)
          end if
       end do
    end do

    ! Sum over all procs
    call comms_reduce("SUM",norm)

    norm = norm**(1.0_dp/real(ord,dp))

    ! Synchronise to ensure all procs hold same data to avoid parallel desync
    call comms_bcast(pub_root_proc_id,norm)

  end function sparse_embed_entrywise_norm


  !============================================================================!
  ! Returns the 1-(induced) norm of a matrix.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   norm (output)  : The calculated 1-norm of the matrix                     !
  !   mat (input)    : The real matrix whose 1-norm is required                !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 19/09/2018.                                    !
  ! Based on SPAM3 version by Jolyon Aarons.                                   !
  !============================================================================!

  subroutine sparse_embed_1norm(norm, mat)

    use sparse, only: sparse_1norm

    implicit none

    ! Arguments
    real(kind=dp), intent(out)     :: norm         ! The 1-norm to compute
    type(SPAM3_EMBED), intent(in)  :: mat

    ! Local variables
    real(kind=DP)               :: tmp_norm
    integer                     :: isub, jsub

    norm = 0.0_dp
    ! rc2013: untested with embedding and almost certainly wrong
    do jsub=1,mat%ncols
       do isub=1,mat%mrows
          call sparse_1norm(tmp_norm, mat%m(isub,jsub))
          norm = norm + tmp_norm
       end do
    end do

  end subroutine sparse_embed_1norm


!-!-----------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine writes a SPAM3_EMBED object to a file.                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   matname  (input) : The (optional) name of the matrix to be written       !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 19th September 2018.                           !
  !============================================================================!

  subroutine sparse_embed_write_scalar(mat, filename, matname, iunit, append)

    use comms, only: pub_on_root
    use constants, only: stdout, normal
    use rundat, only: pub_output_detail
    use sparse, only: sparse_write
    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat       ! Array to write
    character(len=*),  intent(in) :: filename  ! Filename to use
    character(len=*), optional, intent(in) :: matname   ! Name of matrix to print
    integer, optional, intent(in) :: iunit     ! I/O unit
    logical, optional, intent(in) :: append

    ! Local variables
    integer :: isub, jsub
    character(len=256) :: loc_filename      ! Filename to use
    character(len=1) :: isub_str, jsub_str

    ! jcap: although sparse_write can deal with matrices of
    ! SPAM3s, to avoid parallelisation issues I think it's best
    ! to do the regions one at a time
    do jsub = 1, mat%ncols
       do isub = 1, mat%mrows
          ! rc2013: include a region suffix if needed
          if(mat%mrows == 1 .and. mat%ncols == 1) then
             loc_filename = trim(filename)
          else
             write(isub_str,'(i1)') isub
             write(jsub_str,'(i1)') jsub
             loc_filename = trim(filename)//'_'//isub_str//jsub_str
          end if

          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) then
             write(stdout,'(/5a)',advance='no') &
                  'Writing ', matname,' to file "', trim(loc_filename),'" ...'
          end if
          call sparse_write(mat%m(isub,jsub), loc_filename, iunit, append)
          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) write(stdout,'(a)') ' done'
       end do
    end do

  end subroutine sparse_embed_write_scalar

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine writes an vector of SPAM3_EMBED object to a file.          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   matname  (input) : The name of the matrix to be written                  !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 26th October 2018.                             !
  !============================================================================!

  subroutine sparse_embed_write_vector(mat, filename, matname, iunit, append)

    use comms, only: pub_on_root
    use constants, only: stdout, normal
    use rundat, only: pub_output_detail
    use sparse, only: sparse_write, SPAM3, sparse_create, sparse_copy, &
         sparse_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat(:)    ! Array to write
    character(len=*),  intent(in) :: filename  ! Filename to use
    character(len=*), optional, intent(in) :: matname   ! Name of matrix to print
    integer, optional, intent(in) :: iunit     ! I/O unit
    logical, optional, intent(in) :: append   ! append (avoid replace)

    ! Local variables
    integer :: ii, isub, jsub, ierr
    character(len=256) :: loc_filename      ! Filename to use
    character(len=1) :: isub_str, jsub_str
    type(SPAM3), allocatable :: mat_vec(:,:)

    allocate(mat_vec(1,size(mat)),stat=ierr)
    call utils_alloc_check('sparse_write_vector','mat_vec',ierr)
    ! rc2013: assume that all matrices in the array have the same system size
    do jsub = 1, mat(1)%ncols
       do isub = 1, mat(1)%mrows
          do ii=1,size(mat)
             ! rc2013: if this is a single subsystem write as normal
             call sparse_create(mat_vec(1,ii),mat(ii)%m(isub,jsub))
             call sparse_copy(mat_vec(1,ii),mat(ii)%m(isub,jsub))
          end do
          ! rc2013: include a region suffix if needed
          if(mat(1)%mrows == 1 .and. mat(1)%ncols == 1) then
             loc_filename = trim(filename)
          else
             write(isub_str,'(i1)') isub
             write(jsub_str,'(i1)') jsub
             loc_filename = trim(filename)//'_'//isub_str//jsub_str
          end if

          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) then
             write(stdout,'(/5a)',advance='no') &
                  'Writing ', matname,' to file "', trim(loc_filename),'" ...'
          end if
          call sparse_write(mat_vec,loc_filename,iunit,append)
          do ii=1,size(mat)
             call sparse_destroy(mat_vec(1,ii))
          end do
          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) write(stdout,'(a)') ' done'
       end do
    end do
    deallocate(mat_vec,stat=ierr)
    call utils_dealloc_check('sparse_write_vector','mat_vec',ierr)

  end subroutine sparse_embed_write_vector

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine writes a 2D array of SPAM3_EMBED objects to a file.        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   matname  (input) : The name of the matrix to be written                  !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 26th October 2018.                             !
  !============================================================================!

  subroutine sparse_embed_write_matrix(mat, filename, matname, iunit, append)

    use comms, only: pub_on_root
    use constants, only: stdout, normal
    use rundat, only: pub_output_detail
    use sparse, only: sparse_write, SPAM3, sparse_create, sparse_copy, &
         sparse_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat(:,:)  ! Array to write
    character(len=*),  intent(in) :: filename  ! Filename to use
    character(len=*), optional, intent(in) :: matname   ! Name of matrix to print
    integer, optional, intent(in) :: iunit     ! I/O unit
    logical, optional, intent(in) :: append    ! append (avoid replace)

    ! Local variables
    integer :: ii, jj, isub, jsub, ierr
    character(len=256) :: loc_filename      ! Filename to use
    character(len=1) :: isub_str, jsub_str
    type(SPAM3), allocatable :: mat_arr(:,:)

    allocate(mat_arr(size(mat,dim=1),size(mat,dim=2)),stat=ierr)
    call utils_alloc_check('sparse_embed_write_matrix','mat_arr',ierr)
    ! rc2013: assume that all matrices in the array have the same system size
    do jsub = 1, mat(1,1)%ncols
       do isub = 1, mat(1,1)%mrows
          ! rc2013: loop over SPAM3_EMBED matrices
          do jj=1,size(mat,dim=2)
             do ii=1,size(mat,dim=1)
                ! rc2013: extract the component we want
                call sparse_create(mat_arr(ii,jj),mat(ii,jj)%m(isub,jsub))
                call sparse_copy(mat_arr(ii,jj),mat(ii,jj)%m(isub,jsub))
             end do
          end do
          ! rc2013: include a region suffix if needed
          if(mat(1,1)%mrows == 1 .and. mat(1,1)%ncols == 1) then
             loc_filename = trim(filename)
          else
             write(isub_str,'(i1)') isub
             write(jsub_str,'(i1)') jsub
             loc_filename = trim(filename)//'_'//isub_str//jsub_str
          end if

          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) then
             write(stdout,'(/5a)',advance='no') &
                  'Writing ', matname,' to file "', trim(loc_filename),'" ...'
          end if
          call sparse_write(mat_arr,loc_filename,iunit,append)
          do jj=1,size(mat,dim=2)
             do ii=1,size(mat,dim=1)
                call sparse_destroy(mat_arr(ii,jj))
             end do
          end do
          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) write(stdout,'(a)') ' done'
       end do
    end do
    deallocate(mat_arr,stat=ierr)
    call utils_dealloc_check('sparse_embed_write_matrix','mat_arr',ierr)

  end subroutine sparse_embed_write_matrix

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine reads a SPAM3_EMBED object from a file.                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (inout) : The sparse matrix to be read                          !
  !   filename (input) : The filename to use                                   !
  !   matname  (input) : The name of the matrix to be written                  !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 19th September 2018.                           !
  !============================================================================!

  subroutine sparse_embed_read_scalar(mat, filename, matname, iunit, get_unit)

    use comms, only: pub_on_root
    use constants, only: stdout, normal
    use rundat, only: pub_output_detail
    use sparse, only: sparse_read

    implicit none

    type(SPAM3_EMBED), intent(inout) :: mat       ! Array to write
    character(len=*),  intent(in)    :: filename  ! Filename to use
    character(len=*), optional, intent(in) :: matname   ! Name of matrix to print
    integer, optional, intent(inout) :: iunit     ! I/O unit
    logical, optional, intent(inout) :: get_unit

    ! Local variables
    integer :: isub, jsub
    character(len=256) :: loc_filename      ! Filename to use
    character(len=1) :: isub_str, jsub_str

    ! jcap: although sparse_write can deal with matrices of
    ! SPAM3s, to avoid parallelisation issues I think it's
    ! best to do the regions one at a time
    ! jcap: get unit when sparse_read is first called, and
    ! then use this unit
    do jsub=1,mat%ncols
       do isub=1,mat%mrows
          ! rc2013: include a region suffix if needed
          if(mat%mrows == 1 .and. mat%ncols == 1) then
             loc_filename = adjustl(trim(filename))
          else
             write(isub_str,'(i1)') isub
             write(jsub_str,'(i1)') jsub
             loc_filename = trim(filename)//'_'//isub_str//jsub_str
          end if
          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) then
             write(stdout,'(/5a)',advance='no') &
                  'Reading ', matname,' from file "', trim(loc_filename),'" ...'
          end if
          call sparse_read(mat%m(isub,jsub),loc_filename,&
               unit=iunit,get_unit=get_unit)
          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) write(stdout,'(a)') ' done'
       end do
    end do

  end subroutine sparse_embed_read_scalar

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine writes a vector of SPAM3_EMBED object to a file.           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   matname  (input) : The name of the matrix to be written                  !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 26th October 2018.                             !
  !============================================================================!

  subroutine sparse_embed_read_vector(mat, filename, matname, iunit, get_unit)

    use comms, only: pub_on_root
    use constants, only: stdout, normal
    use rundat, only: pub_output_detail
    use sparse, only: sparse_read, SPAM3, sparse_create, sparse_copy, &
         sparse_destroy
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)  :: mat(:)    ! Array to write
    character(len=*),  intent(in)     :: filename  ! Filename to use
    character(len=*), optional, intent(in) :: matname   ! Name of matrix to print
    integer, optional, intent(inout)  :: iunit     ! I/O unit
    logical, optional, intent(inout) :: get_unit

    ! Local variables
    integer :: ii, isub, jsub, ierr
    character(len=256) :: loc_filename      ! Filename to use
    character(len=1) :: isub_str, jsub_str
    logical :: append
    type(SPAM3), allocatable :: mat_copy(:,:)

    allocate(mat_copy(1,size(mat)),stat=ierr)
    call utils_alloc_check('sparse_embed_read_vector','mat_copy',ierr)
    ! rc2013: assume that all matrices in the array have the same system size
    do jsub = 1, mat(1)%ncols
       do isub = 1, mat(1)%mrows
          do ii=1,size(mat)
             ! rc2013: if this is a single subsystem write as normal
             call sparse_create(mat_copy(1,ii),mat(ii)%m(isub,jsub))
          end do
          ! rc2013: include a region suffix if needed
          if(mat(1)%mrows == 1 .and. mat(1)%ncols == 1) then
             loc_filename = adjustl(trim(filename))
          else
             write(isub_str,'(i1)') isub
             write(jsub_str,'(i1)') jsub
             loc_filename = trim(filename)//'_'//isub_str//jsub_str
          end if
          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) then
             write(stdout,'(/5a)',advance='no') &
                  'Reading ', matname,' from file "', trim(loc_filename),'" ...'
          end if
          call sparse_read(mat_copy,loc_filename,iunit,get_unit)
          do ii=1,size(mat)
             call sparse_copy(mat(ii)%m(isub,jsub),mat_copy(1,ii))
             call sparse_destroy(mat_copy(1,ii))
          end do
          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) write(stdout,'(a)') ' done'
       end do
    end do
    deallocate(mat_copy,stat=ierr)
    call utils_dealloc_check('sparse_embed_read_vector','mat_copy',ierr)

  end subroutine sparse_embed_read_vector

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine writes a 2D array of SPAM3_EMBED objects to a file.        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   matname  (input) : The name of the matrix to be written                  !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  ! read_diagonal_blocks (input) : read in diagonal blocks only, leaving off-  !
  ! diagonals zero.                                                            !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 26th October 2018.                             !
  !============================================================================!

  subroutine sparse_embed_read_matrix(mat, filename, matname, iunit, get_unit, &
       read_diagonal_blocks)

    use comms, only: pub_on_root
    use constants, only: stdout, normal
    use rundat, only: pub_output_detail
    use sparse, only: sparse_read, SPAM3, sparse_create, sparse_copy, &
         sparse_destroy, sparse_scale
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: mat(:,:)    ! Array to write
    character(len=*),  intent(in)    :: filename  ! Filename to use
    character(len=*), optional, intent(in) :: matname   ! Name of matrix to print
    integer, optional, intent(inout) :: iunit     ! I/O unit
    logical, optional, intent(inout) :: get_unit
    logical, optional, intent(in)    :: read_diagonal_blocks

    ! Local variables
    integer :: ii, jj, isub, jsub, ierr
    character(len=256) :: loc_filename      ! Filename to use
    character(len=1) :: isub_str, jsub_str
    type(SPAM3), allocatable :: mat_copy(:,:)
    logical :: loc_read_diagonal_blocks

    loc_read_diagonal_blocks = .false.
    if(present(read_diagonal_blocks)) &
         loc_read_diagonal_blocks = read_diagonal_blocks

    allocate(mat_copy(size(mat,dim=1),size(mat,dim=2)),stat=ierr)
    call utils_alloc_check('sparse_embed_read_matrix','mat_copy',ierr)
    ! rc2013: assume that all matrices in the array have the same system size
    do jsub = 1, mat(1,1)%ncols
       do isub = 1, mat(1,1)%mrows
          if(loc_read_diagonal_blocks .and. isub .ne. jsub) then
             do jj=1,size(mat,dim=2)
                do ii=1,size(mat,dim=1)
                   call sparse_scale(mat(ii,jj)%m(isub,jsub),0.0_DP)
                end do
             end do
             cycle
          end if
          ! rc2013: loop over SPAM3_EMBED matrices
          do jj=1,size(mat,dim=2)
             do ii=1,size(mat,dim=1)
                ! rc2013: extract the component we want
                call sparse_create(mat_copy(ii,jj),mat(ii,jj)%m(isub,jsub))
             end do
          end do
          ! rc2013: include a region suffix if needed
          if(mat(1,1)%mrows == 1 .and. mat(1,1)%ncols == 1) then
             loc_filename = adjustl(trim(filename))
          else
             write(isub_str,'(i1)') isub
             write(jsub_str,'(i1)') jsub
             loc_filename = trim(filename)//'_'//isub_str//jsub_str
          end if
          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) then
             write(stdout,'(/5a)',advance='no') &
                  'Reading ', matname,' from file "', trim(loc_filename),'" ...'
          end if
          call sparse_read(mat_copy,loc_filename,iunit,get_unit)
          do jj=1,size(mat,dim=2)
             do ii=1,size(mat,dim=1)
                call sparse_copy(mat(ii,jj)%m(isub,jsub),mat_copy(ii,jj))
                call sparse_destroy(mat_copy(ii,jj))
             end do
          end do
          if (present(matname) .and. pub_on_root .and. &
               (pub_output_detail >= NORMAL)) write(stdout,'(a)') ' done'
       end do
    end do
    deallocate(mat_copy,stat=ierr)
    call utils_dealloc_check('sparse_embed_read_matrix','mat_copy',ierr)

  end subroutine sparse_embed_read_matrix


!-!-----------------------------------------------------------------------------


!------------------------------------------------------------------------------

  !============================================================================!
  !                   CONVERSION ROUTINES FOR SPARSE ARRAYS                    !
  !============================================================================!

!------------------------------------------------------------------------------

  !==========================================================================!
  ! This subroutine takes an embedded SPAM3_ARRAY and transfers it to a 3D   !
  ! SPAM3 embedding matrix with an extra index for spin.                     !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   outmat   (inout) : The embedding sparse matrix.                        !
  !   inarray  (input) : The embedding sparse array.                         !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 25/07/2017.                                  !
  ! This has been designed specifically for ngwf_gradient_mod so probably    !
  ! needs work before it can be used elsewhere.                              !
  ! This is most likely just a placeholder until sparse_arrays are fully     !
  ! incorporated.
  !==========================================================================!

  subroutine sparse_embed_array2mat(outmat, inarray)

    use rundat, only: PUB_1K
    use sparse, only: sparse_copy

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)    :: outmat(:)
    type(SPAM3_EMBED_ARRAY), intent(in) :: inarray

    ! Local Variables
    integer :: num_spins, isub, jsub, is

    ! rc2013: number of subsystems
    num_spins = size(outmat,dim=1)

    outmat%mrows = inarray%mrows
    outmat%ncols = inarray%ncols

    ! rc2013: copy the elements of inarray across to outmat
    do is=1,inarray%num_spins
        do jsub=1,inarray%ncols
            do isub=1,inarray%mrows
                call sparse_copy(outmat(is)%m(isub,jsub), &
                    inarray%m(is,PUB_1K)%m(isub,jsub))
            enddo
        enddo
    enddo

  end subroutine sparse_embed_array2mat

!------------------------------------------------------------------------------

  !==========================================================================!
  ! This subroutine builds an embedding matrix from the elements of another  !
  ! SPAM3 matrix array.                                                      !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   dest  (inout) : The embedded sparse matrix to be built.                !
  !   src   (input) : The embedded sparse matrix to be copied.               !
  !   row   (input) : Index to extract row vector.                           !
  !   col   (input) : Index to extract column vector.                        !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 21/07/2017.                                  !
  !==========================================================================!

  subroutine sparse_embed_extract_sub(dest,src,row,col)

    use sparse, only: sparse_copy
    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: dest
    type(SPAM3_EMBED), intent(in)    :: src
    integer, optional, intent(in)    :: row ! Which row we want to extract
    integer, optional, intent(in)    :: col ! Which column we want to extract

    ! Local Variables
    integer :: mrows, ncols
    integer :: isub

    ! rc2013: pass info to dest based on inputs
    if(present(row) .and. present(col)) then
        ! We only want one 'element'
        mrows = 1
        ncols = 1
        call utils_assert(dest%mrows == mrows, &
            'Error in sparse_embed_extract_sub: rows(dest) =/= 1.')
        call utils_assert(dest%ncols == ncols, &
            'Error in sparse_embed_extract_sub: cols(dest) =/= 1.')
        call sparse_copy(dest%m(1,1),src%m(row,col))

    else if(.not.present(row)) then
        ! Extract column 'vector'
        mrows = src%mrows
        ncols = 1
        call utils_assert(dest%mrows == mrows, &
            'Error in sparse_embed_extract_sub: rows(dest) =/= mrows.')
        call utils_assert(dest%ncols == ncols, &
            'Error in sparse_embed_extract_sub: cols(dest) =/= 1.')
        ! Extract the column vector specified by col
        do isub = 1,mrows
            call sparse_copy(dest%m(isub,1),src%m(isub,col))
        end do

    else if(.not.present(col)) then
        ! Extract row 'vector'
        mrows = 1
        ncols = src%ncols
        call utils_assert(dest%mrows == mrows, &
            'Error in sparse_embed_extract_sub: rows(dest) =/= 1.')
        call utils_assert(dest%ncols == ncols, &
            'Error in sparse_embed_extract_sub: cols(dest) =/= ncols.')
        ! Extract the row vector specified by row
        do isub = 1,ncols
            call sparse_copy(dest%m(1,isub),src%m(row,isub))
        end do
    else
        ! rc2013: nothing to be done here, pass the buck
        ! Doesn't neceessarilly need to abort but do so for testing
        call utils_abort('Error in sparse_embed_extract_sub called but &
            &no row/column information has been provided.')
        call sparse_embed_create(dest,src)
    end if

  end subroutine sparse_embed_extract_sub

!------------------------------------------------------------------------------

  !==========================================================================!
  ! This subroutine takes an SPAM3_EMBED structure and transfers it to a     !
  ! SPAM3_ARRAY.                                                             !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   outmat   (inout) : The embedding sparse matrix.                        !
  !   inarray  (input) : The embedding sparse array.                         !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 13/08/2017.                                  !
  ! Inverse of sparse_embed_array2mat.                                       !
  !==========================================================================!

  subroutine sparse_embed_mat2array(outarray, inmat)

    use rundat, only: PUB_1K
    use sparse, only: sparse_copy

    implicit none

    ! Arguments
    type(SPAM3_EMBED_ARRAY), intent(inout) :: outarray
    type(SPAM3_EMBED), intent(in)    :: inmat(:)

    ! Local Variables
    integer :: mrows, ncols, num_spins, isub, jsub, is

    ! rc2013: number of subsystems
    num_spins = size(inmat,dim=1)
    mrows = inmat(1)%mrows
    ncols = inmat(1)%ncols

    ! rc2013: copy the elements of inarray across to outmat
    do is=1,num_spins
        do jsub=1,ncols
            do isub=1,mrows
                call sparse_copy(outarray%m(is,PUB_1K)%m(isub,jsub), &
                    inmat(is)%m(isub,jsub))
            enddo
        enddo
    enddo

  end subroutine sparse_embed_mat2array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine obtains the maximum eigenvalues of an array of sparse      !
  ! matrices by an iterative conjugate gradients procedure.                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array (input)  : The array of matrices whose eigenvalues are desired     !
  !   met (input)    : The metric to use                                       !
  !   eval (output)  : The eigenvalue estimate                                 !
  !   tol (input)    : Optional tolerance for eigenvalue estimate              !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 15/08/2017.                                    !
  ! Slightly modified to fix calls with optional args by Joseph Prentice,      !
  ! June 2018                                                                  !
  !============================================================================!

  subroutine sparse_embed_extremal_eigenvalue(array, met, eval, tol, min_val, &
          evec, allow_warnings)

    use datatypes, only: FUNCTIONS, data_functions_alloc, data_functions_copy, &
         data_functions_dealloc
    use sparse, only: sparse_extremal_eigenvalue, &
         sparse_extremal_eigenvalue_array, sparse_num_cols

    implicit none

    ! Argument
    type(SPAM3_EMBED), intent(in)       :: array     ! The array of sparse matrices
    type(SPAM3_EMBED), intent(in)       :: met       ! The metric to be used
    real(kind=DP), intent(out)          :: eval      ! Eigenvalue
    real(kind=DP), optional, intent(in) :: tol       ! Optional tolerance
    logical, optional, intent(in)       :: min_val
    ! agrecocmplx
    type(FUNCTIONS), optional, intent(inout) :: evec
    logical, intent(in), optional :: allow_warnings

    ! Local variables
    integer :: isub, num_cols
    real(kind=DP) :: local_tol
    logical :: loc_warnings,local_min_val
    ! jcap: local copy of evec split into parts for each region
    type(FUNCTIONS) :: loc_evec(array%ncols)

    ! jcap: create local copies of optional arguments, if they exist,
    ! with default values
    if (present(tol)) then
       local_tol = tol
    else
       local_tol = 0.01_DP
    end if
    if (present(allow_warnings)) then
       loc_warnings = allow_warnings
    else
       loc_warnings = .false.
    end if
    if (present(min_val)) then
       local_min_val = min_val
    else
       local_min_val = .false.
    endif

    ! jcap: allocate loc_evec
    if (present(evec)) then
       do isub=1,array%ncols
          call data_functions_alloc(loc_evec(isub), &
               int(sparse_num_cols(array%m(1,isub))),evec%iscmplx)
       end do
    end if

    ! Compute each of the eigenvalues -- very simple.
    ! rc2013: simply call the relevant sparse routine. Eventually this
    ! routine should replace sparse_extremal_eigenvalue
    if (present(evec)) then
       if(array%mrows == 1 .and. array%ncols == 1) then
          call sparse_extremal_eigenvalue(array%m(1,1), met%m(1,1), eval, tol=local_tol, &
               min_val=local_min_val, evec=loc_evec(1), allow_warnings=loc_warnings)
       else
          call sparse_extremal_eigenvalue_array(array%m, met%m, eval, tol=local_tol, &
               min_val=local_min_val, evec=loc_evec, allow_warnings=loc_warnings)
       end if
    else
       if(array%mrows == 1 .and. array%ncols == 1) then
          call sparse_extremal_eigenvalue(array%m(1,1), met%m(1,1), eval, tol=local_tol, &
               min_val=local_min_val, allow_warnings=loc_warnings)
       else
          call sparse_extremal_eigenvalue_array(array%m, met%m, eval, tol=local_tol, &
               min_val=local_min_val, allow_warnings=loc_warnings)
       end if
    end if

    ! jcap: copy loc_evec back into evec
    if (present(evec)) then
       num_cols=0
       do isub=1,array%ncols
          call data_functions_copy(evec,loc_evec(isub),starty=num_cols+1,&
               length=int(sparse_num_cols(array%m(1,isub))))
          num_cols=num_cols+int(sparse_num_cols(array%m(1,isub)))
          call data_functions_dealloc(loc_evec(isub))
       end do
    end if

  end subroutine sparse_embed_extremal_eigenvalue

!------------------------------------------------------------------------------


  !============================================================================!
  ! This surborutine prints out the blocking information for an embedding      !
  ! matrix. Useful as diagnostic tool when combining matrices.                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat (inout) : The embedding sparse matrix.                               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 20/08/2017.                                    !
  !============================================================================!

  subroutine sparse_embed_diagnose(mat,mat_name)

    use sparse, only: sparse_diagnose

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in)  :: mat   ! The sparse matrix
    character(len=*), intent(in)   :: mat_name

    ! Local variables
    integer :: isub, jsub, mrows, ncols
    integer :: lib

    ! Write out embedding matrix size
    write(stdout,'(a40,a20)') 'DEBUG: Writing matrix information for ', mat_name
    mrows = mat%mrows
    ncols = mat%ncols
    write(stdout,'(a15,i3,a10,i3)') '     mrows = ', mrows, ', ncols = ', ncols

    do jsub=1,ncols
        do isub=1,mrows
            write(stdout,'(a7,i3,a9,i3,a12,a30)') 'isub = ', isub, ', jsub = ', jsub
            call sparse_diagnose(mat%m(isub,jsub))
        end do
    end do

  end subroutine sparse_embed_diagnose


  !============================================================================!
  !                         EMBEDDING SPARSE ARRAYS                            !
  !             REWRITTEN TO USE SPAM3_EMBED_ARRAY STRUCTURES.                 !
  !============================================================================!

  !============================================================================!
  ! This subroutine creates an embedding SPAM3 array structure. The appropriate!
  ! sub-blocks are passed to sparse_create.                                    !
  !----------------------------------------------------------------------------!
  ! Required arguments:                                                        !
  !   new_spam (output) : The array of sparse matrices to be created           !
  !----------------------------------------------------------------------------!
  ! Optional arguments:                                                        !
  !   spama    (input)  : An array of sparse matrices whose structure is used  !
  !   spamb    (input)  : A second array of sparse matrices whose structure    !
  !                       will be used.                                        !
  !   iscmplx  (input)  : A flag to indicate the data stored in the matrices   !
  !                       of this array is complex.                            !
  !   n_kpoints (input) : Force this number of kpoints in the array.           !
  !   n_spins   (input) : Force this spin dimension in the array.              !
  !   structure (input) : Force this structure in the new array.               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 18/07/2017. Based on sparse_create_embed.      !
  !============================================================================!

  subroutine sparse_embed_array_create_array_array(new_spam,spama,spamb, &
       iscmplx, rlib, n_spins, n_kpoints, structure, mrows, ncols, trans, &
       arow, acol, brow, bcol, rstruc, cstruc)

     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), target,  intent(inout) :: new_spam ! Sparse matrix to be created
     type(SPAM3_EMBED_ARRAY), intent(in)    :: spama    ! First (optional) matrix to use
     type(SPAM3_EMBED_ARRAY), intent(in)    :: spamb    ! Other (optional) matrix to use
     logical, optional, intent(in)     :: iscmplx       ! Flag for complex data (opt)
     integer, optional, intent(out)    :: rlib(:,:)          ! Optional library index
     integer, optional, intent(in) :: n_spins           ! Optional spin dimension.
     integer, optional, intent(in) :: n_kpoints         ! Optional number of kpoints.
     character(len=*), optional, intent(in)   :: structure        ! Structure code
     integer, optional, intent(in) :: mrows
     integer, optional, intent(in) :: ncols
     integer, optional, intent(in) :: arow
     integer, optional, intent(in) :: acol
     integer, optional, intent(in) :: brow
     integer, optional, intent(in) :: bcol
     logical, optional, intent(in) :: trans
     integer, optional, intent(in) :: rstruc
     integer, optional, intent(in) :: cstruc

     ! Local variables
     integer :: is, ik

     call sparse_embed_array_alloc(new_spam, spama=spama, spamb=spamb, &
          n_spins=n_spins, n_kpoints=n_kpoints, structure=structure)

     ! Fill the array using sparse_embed_create.
     do ik = 1, new_spam%num_kpoints
        do is = 1, new_spam%num_spins
            call sparse_embed_create(new_spam%m(is,ik), spama%m(1,1), &
                 spamb%m(1,1), mrows, ncols, iscmplx, rlib, &
                 trans, arow, acol, brow, bcol, rstruc, cstruc)
        end do
     end do

     ! rc2013: now create an array of pointers to new_spam%m(:,PUB_1K)
     ! This will reduce the need to change places where SPAM3_ARRAY is not
     ! fully integrated (I hope)
     !new_spam%p => new_spam%m(:,PUB_1K)
     !do is = 1, new_spam%num_spins
     !   new_spam%spin_p(is)%p => new_spam%m(is,PUB_1K)%m(1,1)
     !end do

  end subroutine sparse_embed_array_create_array_array


  !============================================================================!
  ! This subroutine creates an embedding SPAM3 array structure. The appropriate!
  ! sub-blocks are passed to sparse_create.                                    !
  !----------------------------------------------------------------------------!
  ! Required arguments:                                                        !
  !   new_spam (output) : The array of sparse matrices to be created           !
  !----------------------------------------------------------------------------!
  ! Optional arguments:                                                        !
  !   spama    (input)  : An array of sparse matrices whose structure is used  !
  !   spamb    (input)  : A second array of sparse matrices whose structure    !
  !                       will be used.                                        !
  !   iscmplx  (input)  : A flag to indicate the data stored in the matrices   !
  !                       of this array is complex.                            !
  !   n_kpoints (input) : Force this number of kpoints in the array.           !
  !   n_spins   (input) : Force this spin dimension in the array.              !
  !   structure (input) : Force this structure in the new array.               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 18/07/2017. Based on sparse_create_embed.      !
  !============================================================================!

  subroutine sparse_embed_array_create_array_matrix(new_spam,spama,spammatb, &
       iscmplx, rlib, n_spins, n_kpoints, structure, mrows, ncols, trans, &
       arow, acol, brow, bcol, rstruc, cstruc)

     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), target, intent(inout)  :: new_spam ! Sparse matrix to be created
     type(SPAM3_EMBED_ARRAY), intent(in)     :: spama    ! First (optional) matrix to use
     type(SPAM3_EMBED), intent(in), optional :: spammatb ! Other (optional) matrix to use
     logical, optional, intent(in)     :: iscmplx        ! Flag for complex data (opt)
     integer, optional, intent(out)    :: rlib(:,:)           ! Optional library index
     integer, optional, intent(in) :: n_spins          ! Optional spin dimension.
     integer, optional, intent(in) :: n_kpoints        ! Optional number of kpoints.
     character(len=*), optional, intent(in)   :: structure        ! Structure code
     integer, optional, intent(in) :: mrows
     integer, optional, intent(in) :: ncols
     integer, optional, intent(in) :: arow
     integer, optional, intent(in) :: acol
     integer, optional, intent(in) :: brow
     integer, optional, intent(in) :: bcol
     logical, optional, intent(in) :: trans
     integer, optional, intent(in) :: rstruc
     integer, optional, intent(in) :: cstruc

     ! Local variables
     integer :: is, ik

     if (present(structure)) new_spam%structure = structure
     call sparse_embed_array_alloc(new_spam, spama=spama, spammatb=spammatb, &
          n_spins=n_spins, n_kpoints=n_kpoints, structure=structure)

     ! Fill the array using sparse_embed_create.
     do ik = 1, new_spam%num_kpoints
        do is = 1, new_spam%num_spins
            call sparse_embed_create(new_spam%m(is,ik), spama%m(1,1), &
                 spammatb, mrows, ncols, iscmplx, rlib, &
                 trans, arow, acol, brow, bcol, rstruc, cstruc)
        end do
     end do

     ! rc2013: now create an array of pointers to new_spam%m(:,PUB_1K)
     ! This will reduce the need to change places where SPAM3_ARRAY is not
     ! fully integrated
     !do is = 1, new_spam%num_spins
     !   new_spam%spin_p(is)%p => new_spam%m(is,PUB_1K)%m(1,1)
     !end do

  end subroutine sparse_embed_array_create_array_matrix


  !============================================================================!
  ! This subroutine creates an embedding SPAM3 array structure. The appropriate!
  ! sub-blocks are passed to sparse_create.                                    !
  !----------------------------------------------------------------------------!
  ! Required arguments:                                                        !
  !   new_spam (output) : The array of sparse matrices to be created           !
  !----------------------------------------------------------------------------!
  ! Optional arguments:                                                        !
  !   spama    (input)  : An array of sparse matrices whose structure is used  !
  !   spamb    (input)  : A second array of sparse matrices whose structure    !
  !                       will be used.                                        !
  !   iscmplx  (input)  : A flag to indicate the data stored in the matrices   !
  !                       of this array is complex.                            !
  !   n_kpoints (input) : Force this number of kpoints in the array.           !
  !   n_spins   (input) : Force this spin dimension in the array.              !
  !   structure (input) : Force this structure in the new array.               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 18/07/2017. Based on sparse_create_embed.      !
  !============================================================================!

  subroutine sparse_embed_array_create_matrix_array(new_spam,spammata,spamb, &
       iscmplx, rlib, n_spins, n_kpoints, structure, mrows, ncols, trans, &
       arow, acol, brow, bcol, rstruc, cstruc)

     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), target, intent(inout) :: new_spam ! Sparse matrix to be created
     type(SPAM3_EMBED), intent(in)          :: spammata ! First matrix to use
     type(SPAM3_EMBED_ARRAY), intent(in)    :: spamb    ! Other matrix to use
     logical, optional, intent(in)     :: iscmplx       ! Flag for complex data (opt)
     integer, optional, intent(out)    :: rlib(:,:)          ! Optional library index
     integer, optional, intent(in) :: n_spins           ! Optional spin dimension.
     integer, optional, intent(in) :: n_kpoints         ! Optional number of kpoints.
     character(len=*), optional, intent(in)   :: structure        ! Structure code
     integer, optional, intent(in) :: mrows
     integer, optional, intent(in) :: ncols
     logical, optional, intent(in) :: trans       ! Flag for complex data (opt)
     integer, optional, intent(in) :: arow
     integer, optional, intent(in) :: acol
     integer, optional, intent(in) :: brow
     integer, optional, intent(in) :: bcol
     integer, optional, intent(in) :: rstruc
     integer, optional, intent(in) :: cstruc

     ! Local variables
     integer :: is, ik

     if (present(structure)) new_spam%structure = structure
     call sparse_embed_array_alloc(new_spam, spammata=spammata, spamb=spamb, &
         n_spins=n_spins, n_kpoints=n_kpoints, &
         structure=structure)

     ! Fill the array using sparse_embed_create
     do ik = 1, new_spam%num_kpoints
        do is = 1, new_spam%num_spins
            call sparse_embed_create(new_spam%m(is,ik), spammata, spamb%m(1,1), &
                 mrows, ncols, iscmplx, rlib, &
                 trans, arow, acol, brow, bcol, rstruc, cstruc)
        end do
     end do

     ! rc2013: now create an array of pointers to new_spam%m(:,PUB_1K)
     ! This will reduce the need to change places where SPAM3_ARRAY is not
     ! fully integrated
     !do is = 1, new_spam%num_spins
     !   new_spam%spin_p(is)%p => new_spam%m(is,PUB_1K)%m(1,1)
     !end do

  end subroutine sparse_embed_array_create_matrix_array


  !============================================================================!
  ! This subroutine creates an embedding SPAM3 array structure. The appropriate!
  ! sub-blocks are passed to sparse_create.                                    !
  !----------------------------------------------------------------------------!
  ! Required arguments:                                                        !
  !   new_spam (output) : The array of sparse matrices to be created           !
  !----------------------------------------------------------------------------!
  ! Optional arguments:                                                        !
  !   spama    (input)  : An array of sparse matrices whose structure is used  !
  !   spamb    (input)  : A second array of sparse matrices whose structure    !
  !                       will be used.                                        !
  !   iscmplx  (input)  : A flag to indicate the data stored in the matrices   !
  !                       of this array is complex.                            !
  !   n_kpoints (input) : Force this number of kpoints in the array.           !
  !   n_spins   (input) : Force this spin dimension in the array.              !
  !   structure (input) : Force this structure in the new array.               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 18/07/2017. Based on sparse_create_embed.      !
  !============================================================================!

  subroutine sparse_embed_array_create_matrix_matrix(new_spam,spammata,spammatb, &
       iscmplx, rlib, n_spins, n_kpoints, structure, mrows, ncols, trans, &
       arow, acol, brow, bcol, rstruc, cstruc)

     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), target, intent(inout)  :: new_spam ! Sparse matrix to be created
     type(SPAM3_EMBED), intent(in), optional :: spammata ! First (optional) matrix to use
     type(SPAM3_EMBED), intent(in), optional :: spammatb ! Other (optional) matrix to use
     logical, optional, intent(in)     :: iscmplx        ! Flag for complex data (opt)
     integer, optional, intent(out)    :: rlib(:,:)           ! Optional library index
     integer, optional, intent(in) :: n_spins          ! Optional spin dimension.
     integer, optional, intent(in) :: n_kpoints        ! Optional number of kpoints.
     character(len=*), optional, intent(in)   :: structure        ! Structure code
     integer, optional, intent(in) :: mrows
     integer, optional, intent(in) :: ncols
     logical, optional, intent(in) :: trans       ! Flag for complex data (opt)
     integer, optional, intent(in) :: arow
     integer, optional, intent(in) :: acol
     integer, optional, intent(in) :: brow
     integer, optional, intent(in) :: bcol
     integer, optional, intent(in) :: rstruc
     integer, optional, intent(in) :: cstruc

     ! Local variables
     integer :: is, ik

     call sparse_embed_array_alloc(new_spam, spammata=spammata, spammatb=spammatb, &
          n_spins=n_spins, n_kpoints=n_kpoints, structure=structure, &
          mrows=mrows, ncols=ncols)

     ! rc2013: cycle over all subsystems
     do ik = 1, new_spam%num_kpoints
        do is = 1, new_spam%num_spins
           call sparse_embed_create(new_spam%m(is,ik), spammata, spammatb, &
                mrows, ncols, iscmplx, rlib, &
                trans, arow, acol, brow, bcol, rstruc, cstruc)
        end do
     end do

     ! rc2013: now create an array of pointers to new_spam%m(:,PUB_1K)
     ! This will reduce the need to change places where SPAM3_ARRAY is not
     ! fully integrated
     !do is = 1, new_spam%num_spins
     !   new_spam%spin_p(is)%p => new_spam%m(is,PUB_1K)%m(1,1)
     !end do

  end subroutine sparse_embed_array_create_matrix_matrix

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine allocates a new array of sparse matrices and assigns a     !
  ! structure to them.                                                         !
  !----------------------------------------------------------------------------!
  ! Required arguments:                                                        !
  !   new_spam (output) : The array of sparse matrices to be created           !
  !----------------------------------------------------------------------------!
  ! Optional arguments:
  !   spama    (input)  : An array of sparse matrices whose structure is used  !
  !   spamb    (input)  : A second array of sparse matrices whose structure    !
  !                       will be used.                                        !
  !   spammata, spammatb (input) : non-array alternatives for spama and spamb  !
  !   n_spins   (input) : Force this spin dimension in the array.              !
  !   n_kpoints (input) : Force this number of kpoints in the array.           !
  !   structure (input) : Force this structure in the new array.               !
  !============================================================================!

  subroutine sparse_embed_array_alloc(new_spam, spama, spamb, spammata, spammatb, &
       n_spins, n_kpoints, structure, mrows, ncols)

    use rundat, only: pub_ngwf_regions_ngroups
    use utils, only: utils_abort, utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED_ARRAY), target, intent(inout) :: new_spam      ! Sparse matrix array to be created
    type(SPAM3_EMBED_ARRAY), optional, intent(in) :: spama  ! First (optional) matrix array to use
    type(SPAM3_EMBED_ARRAY), optional, intent(in) :: spamb  ! Other (optional) matrix array to use
    type(SPAM3_EMBED), optional, intent(in) :: spammata     ! Optional matrix to use (instead of spama)
    type(SPAM3_EMBED), optional, intent(in) :: spammatb     ! Optional matrix to use (instead of spamb)
    integer, optional, intent(in) :: n_spins          ! Optional spin dimension.
    integer, optional, intent(in) :: n_kpoints        ! Optional number of kpoints.
    character(len=*), optional, intent(in) :: structure  ! Optional structure string.
    integer, optional, intent(in) :: mrows          ! Optional spin dimension.
    integer, optional, intent(in) :: ncols        ! Optional number of kpoints.

    ! Local variables
    integer :: ierr        ! Error flag
    integer :: loc_mrows, loc_ncols

    ! jcap: default to total number of regions
    loc_mrows = pub_ngwf_regions_ngroups
    loc_ncols = pub_ngwf_regions_ngroups
    if (present(mrows)) loc_mrows=mrows
    if (present(ncols)) loc_ncols=ncols

    ! Check combination of optional sparse array arguments.
    if ( (present(spamb) .or. present(spammatb)) .and. &
         (.not. ( present(spama) .or. present(spammata) ) ) ) then
       call utils_abort('Error in sparse_embed_array_alloc: &
            &spamb or spammatb present but spama and spammata missing.')
    end if
    if ( present(spama) .and. present(spammata) ) then
       call utils_abort('Error in sparse_embed_array_alloc: &
            &you cannot use spammata and spama simultaneously.')
    end if
    if ( present(spamb) .and. present(spammatb) ) then
       call utils_abort('Error in sparse_embed_array_alloc: &
            &you cannot use spammatb and spamb simultaneously.')
    end if

    ! Check consistency of optional sparse array arguments.
    if (present(spama)) then
       if ( sparse_embed_array_check(spama) /= 1 ) then
          call utils_abort('Error in sparse_embed_array_alloc: &
               &spama fails the sparse_array check.')
       end if
    end if
    if (present(spamb)) then
       if ( sparse_embed_array_check(spamb) /= 1 ) then
          call utils_abort('Error in sparse_embed_array_alloc: &
               &spamb fails the sparse_array check.')
       end if
    end if

    ! Check that the new sparse array is really empty.
    ierr = sparse_embed_array_check(new_spam)
    call utils_assert(ierr == 0, 'Error in sparse_embed_array_alloc: &
         &new_spam not empty.  sparse_embed_array_check returned: ', ierr)

    ! Set the number of kpoints of the new object.
    if ( present(n_kpoints) ) then
       new_spam%num_kpoints = n_kpoints
    else  if ( (present(spama)) .and. (present(spamb)) ) then
       if ( spama%num_kpoints == spamb%num_kpoints ) then
          new_spam%num_kpoints = spama%num_kpoints
       else if ( ( spama%num_kpoints == 1 ) .or. (spamb%num_kpoints == 1 ) ) then
          new_spam%num_kpoints = spama%num_kpoints * spamb%num_kpoints
       else
          call utils_abort('Error in sparse_embed_array_alloc: &
               &mismatch in number of kpoints of spama and spamb')
       end if
    else if ( present(spama) ) then
       new_spam%num_kpoints = spama%num_kpoints
    else if ( present(spamb) ) then
       new_spam%num_kpoints = spamb%num_kpoints
    else
       call utils_abort('Error in sparse_embed_array_alloc: &
            &undefined number of kpoints to allocate')
    end if

    ! Set the spin dimension of the new object.
    if ( present(n_spins) ) then
       new_spam%num_spins = n_spins
    else  if ( (present(spama)) .and. (present(spamb)) ) then
       if ( spama%num_spins == spamb%num_spins ) then
          new_spam%num_spins = spama%num_spins
       else if ( ( spama%num_spins == 1 ) .or. (spamb%num_spins == 1 ) ) then
          new_spam%num_spins = spama%num_spins * spamb%num_spins
       else
          call utils_abort('Error in sparse_embed_array_alloc: &
               &mismatch in number of spins of spama and spamb')
       end if
    else if ( present(spama) ) then
       new_spam%num_spins = spama%num_spins
    else if ( present(spamb) ) then
       new_spam%num_spins = spamb%num_spins
    else
       call utils_abort('Error in sparse_embed_array_alloc: &
            &undefined number of spins to allocate')
    end if

    ! rc2013: set the number of sub-block rows and columns
    new_spam%mrows = loc_mrows
    new_spam%ncols = loc_ncols
    if(present(spama) .and. present(spamb)) then
       call utils_assert(spama%ncols == spamb%mrows, &
            'Error in sparse_embed_array_alloc: cols(A) =/= rows(B)')
       new_spam%mrows = spama%mrows
       new_spam%ncols = spamb%ncols
       ! rc2013: allocate "macro" structure too
       new_spam%structure = trim(spama%structure)//trim(spamb%structure)
    else if(present(spama)) then
       new_spam%mrows = spama%mrows
       new_spam%ncols = spama%ncols
       new_spam%structure = trim(spama%structure)
    else if(present(spammata) .and. present(spammatb)) then
       call utils_assert(spammata%ncols == spammatb%mrows, &
            'Error in sparse_embed_array_alloc: cols(A) =/= rows(B)')
       new_spam%mrows = spammata%mrows
       new_spam%ncols = spammatb%ncols
       new_spam%structure = trim(spammata%structure)//trim(spammata%structure)
    else if(present(spammata)) then
       new_spam%mrows = spammata%mrows
       new_spam%ncols = spammata%ncols
    end if

    ! rc2013: allocate "macro" structure
    if(present(spama) .and. present(spamb)) then
       new_spam%structure = trim(spama%structure)//trim(spamb%structure)
    else if(present(spama) .and. present(spammatb)) then
       new_spam%structure = trim(spama%structure)//trim(spammatb%structure)
    else if(present(spama)) then
       new_spam%structure = trim(spama%structure)
    else if(present(spammata) .and. present(spammatb)) then
       new_spam%structure = trim(spammata%structure)//trim(spammata%structure)
    else if(present(spammata)) then
       new_spam%structure = trim(spammata%structure)
    else if (present(structure)) then
       new_spam%structure = trim(structure)
    else
       call utils_abort('Error in sparse_embed_array_alloc: &
            &invalid structure input.')
    end if


    ! Allocate the array element.
    allocate(new_spam%m(new_spam%num_spins, new_spam%num_kpoints), stat=ierr)
    ! Tag differs from variable so that it matches the tag in
    ! sparse_embed_array_destroy
    call utils_alloc_check('sparse_embed_array_alloc', 'array%m', ierr)

    ! Initialise structures
    if (present(structure)) new_spam%m(:,:)%structure = structure

  end subroutine sparse_embed_array_alloc

  !============================================================================!
  ! This subroutine destroys an array of SPAM3_EMBED matrices, freeing up      !
  ! the memory allocated to internal arrays.                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array (inout) : The array of sparse matrices to be destroyed.            !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 18/07/2017.                                    !
  !============================================================================!

  subroutine sparse_embed_array_destroy(array)

     use utils, only: utils_dealloc_check
     implicit none

     ! Argument
     type(SPAM3_EMBED_ARRAY), intent(inout) :: array

     ! Local variables
     integer :: is, ik
     integer :: ierr

     ! Deallocate memory assigned to internal structures
     do ik = 1, array%num_kpoints
        do is = 1, array%num_spins
            ! Destroy each of the SPAM3_EMBED objects.
            call sparse_embed_destroy(array%m(is,ik))
        end do
     end do

     !do is = 1, array%num_spins
     !   nullify(array%spin_p(is)%p)
     !end do

     ! Deallocate the array itself.
     deallocate(array%m, stat=ierr)
     array%num_spins = 0
     array%num_kpoints = 0
     call utils_dealloc_check('sparse_embed_array_destroy', 'array%m', ierr)

  end subroutine sparse_embed_array_destroy

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues               !
  ! of the matrices in an array of sparse matrices (scalar parameters).        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array  (inout) : The array of sparse matrices to be rescaled and shifted !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional real shift parameter                       !
  !============================================================================!

  subroutine sparse_embed_array_scale_real_scalar(array, alpha, beta)

     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout)    :: array  ! The array to be operated on
     real(kind=DP), intent(in)           :: alpha  ! The scaling parameter
     real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

     ! Local variables
     integer :: ispin, kpoint

     ! Check consistency of sparse array argument.
     if ( sparse_embed_array_check(array) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_scale_real_scalar: &
             &array fails the sparse_array check.')
     end if

     ! Scale each of the matrices.
     do kpoint = 1, array%num_kpoints
        do ispin = 1, array%num_spins
           call sparse_embed_scale(array%m(ispin,kpoint), alpha, beta)
        end do
     end do

  end subroutine sparse_embed_array_scale_real_scalar


  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues               !
  ! of the matrices in an array of sparse matrices (parameters form arrays).   !
  ! The shapes of the array, alpha, and beta (if present) must conform.        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array  (inout) : The array of sparse matrices to be rescaled and shifted !
  !   alpha  (input) : The real array of scaling parameters.                   !
  !   beta   (input) : The optional real array of shift parameters.            !
  !============================================================================!

  subroutine sparse_embed_array_scale_real_array(array, alpha, beta)

     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout)                    :: array  ! The array to be operated on
     real(kind=DP), dimension(:,:), intent(in)           :: alpha  ! The scaling parameters
     real(kind=DP), dimension(:,:), optional, intent(in) :: beta   ! The shift parameters

     ! Local variables
     integer :: is, ik

     ! Check shapes
     if ( present(beta) ) then
        if ( utils_compare_vectors(shape(alpha), shape(beta)) /= 0 ) then
           call utils_abort('Error in sparse_embed_array_scale_real_array: &
                &shapes of alpha and beta do not match.')
        end if
     end if
     if ( utils_compare_vectors(shape(array%m), shape(alpha)) /= 0 ) then
        call utils_abort('Error in sparse_embed_array_scale_real_array: &
             &shapes of array%m and alpha do not match.')
     end if

     ! Scale each of the matrices.
     do ik = 1, array%num_kpoints
        do is = 1, array%num_spins
           if( present(beta) ) then
              call sparse_embed_scale( array%m(is,ik), alpha(is,ik), beta(is,ik) )
           else
              call sparse_embed_scale( array%m(is,ik), alpha(is,ik) )
           end if
        end do
     end do

  end subroutine sparse_embed_array_scale_real_array

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs an axpy operation on the arrays of matrices:      !
  !   y := y + alpha * x                                                       !
  ! The shapes of ymat and xmat must conform.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The array of sparse matrices y                           !
  !   xmat  (input) : The array of sparse matrices x                           !
  !   alpha (input) : The real scalar parameter alpha                          !
  !============================================================================!

  subroutine sparse_embed_array_axpy_real_scalar(ymat, xmat, alpha)

     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: ymat    ! The array of sparse matrices y
     type(SPAM3_EMBED_ARRAY), intent(in) :: xmat       ! The array of sparse matrices x
     real(kind=DP), intent(in) :: alpha          ! The parameter alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_embed_array_check(xmat) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_axpy_real_array: &
             &xmat fails the sparse_array check.')
     else if ( sparse_embed_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_axpy_real_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Check that shapes conform.
     if ( utils_compare_vectors(shape(xmat%m), shape(ymat%m)) /= 0 ) then
        call utils_abort('Error in sparse_embed_array_axpy_real_array: &
             &shapes of xmat%m and ymat%m do not match.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_embed_axpy_real(ymat%m(ispin,kpoint), xmat%m(ispin,kpoint), alpha)
        end do
     end do

  end subroutine sparse_embed_array_axpy_real_scalar

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs an axpy operation on the arrays of matrices:      !
  !   y := y + alpha * x                                                       !
  ! The shapes of all three arrays must conform.                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The array of sparse matrices y                           !
  !   xmat  (input) : The array of sparse matrices x                           !
  !   alpha (input) : The array of real parameters alpha                       !
  !============================================================================!

  subroutine sparse_embed_array_axpy_real_array(ymat, xmat, alpha)

     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: ymat           ! Array of sparse matrices y
     type(SPAM3_EMBED_ARRAY), intent(in) :: xmat              ! Array of sparse matrices x
     real(kind=DP), dimension(:,:), intent(in) :: alpha ! Array of parameters alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_embed_array_check(xmat) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_axpy_real_array: &
             &xmat fails the sparse_array check.')
     else if ( sparse_embed_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_axpy_real_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Check that shapes conform.
     if ( utils_compare_vectors(shape(xmat%m), shape(ymat%m)) /= 0 ) then
        call utils_abort('Error in sparse_embed_array_axpy_real_array: &
             &shapes of xmat%m and ymat%m do not match.')
     end if
     if ( utils_compare_vectors(shape(alpha), shape(xmat%m)) /= 0 ) then
        call utils_abort('Error in sparse_embed_array_axpy_real_array: &
             &shapes of alpha and xmat%m do not match.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_embed_axpy_real(ymat%m(ispin,kpoint), xmat%m(ispin,kpoint), alpha(ispin,kpoint))
        end do
     end do

  end subroutine sparse_embed_array_axpy_real_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs an axpy operation on the arrays of matrices:      !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The array of sparse matrices y                           !
  !   xmat  (input) : A sparse matrix x                                        !
  !   alpha (input) : The real scalar parameter alpha                          !
  !============================================================================!

  subroutine sparse_embed_array_axpy_matrix_real_scalar(ymat, xmat, alpha)

     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: ymat    ! The array of sparse matrices y
     type(SPAM3_EMBED), intent(in) :: xmat             ! The sparse matrix x
     real(kind=DP), intent(in) :: alpha          ! The parameter alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_embed_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_axpy_matrix_real_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_embed_axpy_real(ymat%m(ispin,kpoint), xmat, alpha)
        end do
     end do

  end subroutine sparse_embed_array_axpy_matrix_real_scalar

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs an axpy operation on the arrays of matrices:      !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The array of sparse matrices y                           !
  !   xmat  (input) : The sparse matrix x                                      !
  !   alpha (input) : The array of real parameters alpha                       !
  !============================================================================!

  subroutine sparse_embed_array_axpy_matrix_real_array(ymat, xmat, alpha)

     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: ymat           ! Array of sparse matrices y
     type(SPAM3_EMBED), intent(in) :: xmat                    ! Sparse matrix x
     real(kind=DP), dimension(:,:), intent(in) :: alpha ! Array of parameters alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_embed_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_axpy_matrix_real_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Check that shapes conform.
     if ( utils_compare_vectors(shape(alpha), shape(ymat%m)) /= 0 ) then
        call utils_abort('Error in sparse_embed_array_axpy_matrix_real_array: &
             &shapes of alpha and ymat%m do not match.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_embed_axpy_real(ymat%m(ispin,kpoint), xmat, alpha(ispin,kpoint))
        end do
     end do

  end subroutine sparse_embed_array_axpy_matrix_real_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues               !
  ! of the matrices in an array of sparse matrices (scalar parameters).        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array  (inout) : The array of sparse matrices to be rescaled and shifted !
  !   alpha  (input) : The complex scaling parameter                           !
  !   beta   (input) : The optional real shift parameter                       !
  !============================================================================!

  !subroutine sparse_embed_array_scale_complex_scalar(array, alpha, beta)

  !   use sparse, only: sparse_scale
  !   use utils, only: utils_abort
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout)    :: array  ! The array to be operated on
  !   complex(kind=DP), intent(in)        :: alpha  ! The scaling parameter
  !   real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

  !   ! Local variables
  !   integer :: ispin, kpoint

  !   ! Check consistency of sparse array argument.
  !   if ( sparse_embed_array_check(array) /= 1 ) then
  !      call utils_abort('Error in sparse_embed_array_scale_complex_scalar: &
  !           &array fails the sparse_array check.')
  !   end if

  !   ! Scale each of the matrices.
  !   do kpoint = 1, array%num_kpoints
  !      do ispin = 1, array%num_spins
  !         call sparse_embed_scale(array%m(ispin,kpoint), alpha, beta)
  !      end do
  !   end do

  !end subroutine sparse_embed_array_scale_complex_scalar

!-!-----------------------------------------------------------------------------

  !!============================================================================!
  !! This subroutine scales and optionally shifts the eigenvalues               !
  !! of the matrices in an array of sparse matrices (parameters form arrays).   !
  !! The shapes of the array, alpha, and beta (if present) must conform.        !
  !!----------------------------------------------------------------------------!
  !! Arguments:                                                                 !
  !!   array  (inout) : The array of sparse matrices to be rescaled and shifted !
  !!   alpha  (input) : The complex array of scaling parameters.                !
  !!   beta   (input) : The optional real array of scaling parameters.          !
  !!============================================================================!

  !subroutine sparse_embed_array_scale_complex_array(array, alpha, beta)

  !   use sparse, only: sparse_scale
  !   use utils, only: utils_abort, utils_compare_vectors
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout)                    :: array  ! The array to be operated on
  !   complex(kind=DP), dimension(:,:), intent(in)        :: alpha  ! The scaling parameters
  !   real(kind=DP), dimension(:,:), optional, intent(in) :: beta   ! The shift parameters

  !   ! Local variables
  !   integer :: ispin, kpoint

  !   ! Check shapes
  !   if ( present(beta) ) then
  !      if ( utils_compare_vectors(shape(alpha), shape(beta)) /= 0 ) then
  !         call utils_abort('Error in sparse_embed_array_scale_complex_array: &
  !              &shapes of alpha and beta do not match.')
  !      end if
  !   end if
  !   if ( utils_compare_vectors(shape(array%m), shape(alpha)) /= 0 ) then
  !      call utils_abort('Error in sparse_embed_array_scale_complex_array: &
  !           &shapes of array%m and alpha do not match.')
  !   end if

  !   ! Scale each of the matrices.
  !   do kpoint = 1, array%num_kpoints
  !      do ispin = 1, array%num_spins
  !         call sparse_embed_scale( array%m(ispin,kpoint), alpha(ispin,kpoint), beta(ispin,kpoint))
  !      end do
  !   end do

  !end subroutine sparse_embed_array_scale_complex_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs an axpy operation on the arrays of matrices:      !
  !   y := y + alpha * x                                                       !
  ! The shapes of ymat and xmat must conform.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The array of sparse matrices y                           !
  !   xmat  (input) : The array of sparse matrices x                           !
  !   alpha (input) : The complex parameter alpha                              !
  !============================================================================!

  !subroutine sparse_embed_array_axpy_complex_scalar(ymat, xmat, alpha)

  !   use sparse, only: sparse_axpy
  !   use utils, only: utils_abort, utils_compare_vectors
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout) :: ymat    ! The array of sparse matrices y
  !   type(SPAM3_EMBED_ARRAY), intent(in) :: xmat       ! The array of sparse matrices x
  !   complex(kind=DP), intent(in) :: alpha       ! The parameter alpha

  !   ! Local variables
  !   integer :: ispin, kpoint

  !   ! Check arrays of sparse matrices.
  !   if ( sparse_embed_array_check(xmat) /= 1 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_complex_scalar: &
  !           &xmat fails the sparse_array check.')
  !   else if ( sparse_embed_array_check(ymat) /= 1 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_complex_scalar: &
  !           &ymat fails the sparse_array check.')
  !   end if

  !   ! Check that shapes conform.
  !   if ( utils_compare_vectors(shape(xmat%m), shape(ymat%m)) /= 0 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_complex_scalar: &
  !           &shapes of xmat%m and ymat%m do not match.')
  !   end if

  !   ! Axpy each of the matrices
  !   do kpoint = 1, ymat%num_kpoints
  !      do ispin = 1, ymat%num_spins
  !         call sparse_embed_axpy(ymat%m(ispin,kpoint), xmat%m(ispin,kpoint), alpha)
  !      end do
  !   end do

  !end subroutine sparse_embed_array_axpy_complex_scalar

!-!-----------------------------------------------------------------------------

  !!============================================================================!
  !! This subroutine performs an axpy operation on the arrays of matrices:      !
  !!   y := y + alpha * x                                                       !
  !! The shapes of all three arrays must conform.                               !
  !!----------------------------------------------------------------------------!
  !! Arguments:                                                                 !
  !!   ymat  (inout) : The array of sparse matrices y                           !
  !!   xmat  (input) : The array of sparse matrices x                           !
  !!   alpha (input) : The array of complex parameters alpha                    !
  !!============================================================================!

  !subroutine sparse_embed_array_axpy_complex_array(ymat, xmat, alpha)

  !   use sparse, only: sparse_axpy
  !   use utils, only: utils_abort, utils_compare_vectors
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout) :: ymat              ! Array of sparse matrices y
  !   type(SPAM3_EMBED_ARRAY), intent(in) :: xmat                 ! Array of sparse matrices x
  !   complex(kind=DP), dimension(:,:), intent(in) :: alpha ! Array of parameters alpha

  !   ! Local variables
  !   integer :: ispin, kpoint

  !   ! Check arrays of sparse matrices.
  !   if ( sparse_embed_array_check(xmat) /= 1 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_complex_array: &
  !           &xmat fails the sparse_array check.')
  !   else if ( sparse_embed_array_check(ymat) /= 1 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_complex_array: &
  !           &ymat fails the sparse_array check.')
  !   end if

  !   ! Check that shapes conform.
  !   if ( utils_compare_vectors(shape(xmat%m), shape(ymat%m)) /= 0 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_complex_array: &
  !           &shapes of xmat%m and ymat%m do not match.')
  !   end if
  !   if ( utils_compare_vectors(shape(alpha), shape(xmat%m)) /= 0 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_complex_array: &
  !           &shapes of alpha and xmat%m do not match.')
  !   end if

  !   ! Axpy each of the matrices
  !   do kpoint = 1, ymat%num_kpoints
  !      do ispin = 1, ymat%num_spins
  !         call sparse_embed_axpy(ymat%m(ispin,kpoint), xmat%m(ispin,kpoint), alpha(ispin,kpoint))
  !      end do
  !   end do

  !end subroutine sparse_embed_array_axpy_complex_array

!-!-----------------------------------------------------------------------------

  !!============================================================================!
  !! This subroutine performs an axpy operation on the arrays of matrices:      !
  !!   y := y + alpha * x                                                       !
  !!----------------------------------------------------------------------------!
  !! Arguments:                                                                 !
  !!   ymat  (inout) : The array of sparse matrices y                           !
  !!   xmat  (input) : A sparse matrix x                                        !
  !!   alpha (input) : The complex scalar parameter alpha                       !
  !!============================================================================!

  !subroutine sparse_embed_array_axpy_matrix_complex_scalar(ymat, xmat, alpha)

  !   use sparse, only: sparse_axpy
  !   use utils, only: utils_abort, utils_compare_vectors
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout) :: ymat    ! The array of sparse matrices y
  !   type(SPAM3_EMBED), intent(in) :: xmat             ! The sparse matrix x
  !   complex(kind=DP), intent(in) :: alpha       ! The parameter alpha

  !   ! Local variables
  !   integer :: ispin, kpoint

  !   ! Check arrays of sparse matrices.
  !   if ( sparse_embed_array_check(ymat) /= 1 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_matrix_complex_array: &
  !           &ymat fails the sparse_array check.')
  !   end if

  !   ! Axpy each of the matrices
  !   do kpoint = 1, ymat%num_kpoints
  !      do ispin = 1, ymat%num_spins
  !         call sparse_embed_axpy(ymat%m(ispin,kpoint), xmat, alpha)
  !      end do
  !   end do

  !end subroutine sparse_embed_array_axpy_matrix_complex_scalar

!-!-----------------------------------------------------------------------------

  !!============================================================================!
  !! This subroutine performs an axpy operation on the arrays of matrices:      !
  !!   y := y + alpha * x                                                       !
  !!----------------------------------------------------------------------------!
  !! Arguments:                                                                 !
  !!   ymat  (inout) : The array of sparse matrices y                           !
  !!   xmat  (input) : The sparse matrix x                                      !
  !!   alpha (input) : The array of complex parameters alpha                    !
  !!============================================================================!

  !subroutine sparse_embed_array_axpy_matrix_complex_array(ymat, xmat, alpha)

  !   use sparse, only: sparse_axpy
  !   use utils, only: utils_abort, utils_compare_vectors
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout) :: ymat           ! Array of sparse matrices y
  !   type(SPAM3_EMBED), intent(in) :: xmat                    ! Sparse matrix x
  !   complex(kind=DP), dimension(:,:), intent(in) :: alpha ! Array of parameters alpha

  !   ! Local variables
  !   integer :: ispin, kpoint

  !   ! Check arrays of sparse matrices.
  !   if ( sparse_embed_array_check(ymat) /= 1 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_matrix_complex_array: &
  !           &ymat fails the sparse_array check.')
  !   end if

  !   ! Check that shapes conform.
  !   if ( utils_compare_vectors(shape(alpha), shape(ymat%m)) /= 0 ) then
  !      call utils_abort('Error in sparse_embed_array_axpy_matrix_complex_array: &
  !           &shapes of alpha and ymat%m do not match.')
  !   end if

  !   ! Axpy each of the matrices
  !   do kpoint = 1, ymat%num_kpoints
  !      do ispin = 1, ymat%num_spins
  !         call sparse_embed_axpy(ymat%m(ispin,kpoint), xmat, alpha(ispin,kpoint))
  !      end do
  !   end do

  !end subroutine sparse_embed_array_axpy_matrix_complex_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine copies the data from one array of sparse matrices          !
  ! to another.                                                                !
  ! The copy is performed if either the sizes match or the mismatch is due to  !
  ! one or multiple sizes being 1 in the source, in which case the             !
  ! source is copied to all corresponding elements in the destination array.   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest (inout) : The destination array of sparse matrices                  !
  !   src  (input) : The source array of sparse matrices                       !
  !============================================================================!

  subroutine sparse_embed_array_copy_array(dest, src)

     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: dest   ! The destination array of sparse matrices
     type(SPAM3_EMBED_ARRAY), intent(in) :: src       ! The source array of sparse matrices

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_embed_array_check(src) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_copy_array: &
             &src fails the sparse_array check.')
     else if ( sparse_embed_array_check(dest) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_copy_array: &
             &dest fails the sparse_array check.')
     end if

     ! Perform the copy, guessing how to match matrices if needed and possible.
     if ( ( dest%num_kpoints == src%num_kpoints ) .and. &
          ( dest%num_spins == src%num_spins ) ) then
        ! Both arrays have the same shape.
        do kpoint = 1, src%num_kpoints
           do ispin = 1, src%num_spins
              call sparse_embed_copy( dest%m(ispin, kpoint) , src%m(ispin, kpoint) )
           end do
        end do
     else if ( ( src%num_kpoints == 1 ) .and. ( src%num_spins == 1 ) ) then
        ! The source array is 1x1.
        do kpoint = 1, dest%num_kpoints
           do ispin = 1, dest%num_spins
              call sparse_embed_copy( dest%m(ispin, kpoint), src%m(1,1) )
           end do
        end do
     else if ( ( src%num_kpoints == 1) .and. ( dest%num_spins == src%num_spins ) ) then
        ! The source array only has one kpoint
        do kpoint = 1, dest%num_kpoints
           do ispin = 1, src%num_spins
              call sparse_embed_copy( dest%m(ispin, kpoint), src%m(ispin, 1) )
           end do
        end do
     else if ( ( src%num_kpoints == dest%num_spins) .and. ( src%num_spins == 1) ) then
        ! The source array only has one spin
        do kpoint = 1, src%num_kpoints
           do ispin = 1, dest%num_spins
              call sparse_embed_copy( dest%m(ispin, kpoint), src%m(1, kpoint) )
           end do
        end do
     else if ( src%num_kpoints /= dest%num_kpoints ) then
        ! We don't know how to copy one array onto the other.
        call utils_abort('Error in sparse_embed_array_copy_array: &
             &cannot guess how to match kpoints dimensions')
     else if ( src%num_spins /= dest%num_spins ) then
        ! We don't know how to copy one array onto the other.
        call utils_abort('Error in sparse_embed_array_copy_array: &
             &cannot guess how to match spin dimensions')
     else
        ! We should never arrive here.
        call utils_abort('Undetermined error in sparse_embed_array_copy_array.')
     end if

  end subroutine sparse_embed_array_copy_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine copies the data from a sparse matrix to all matrices in an !
  ! array of sparse matrices.                                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest (inout) : The destination array of sparse matrices                  !
  !   src  (input) : The source sparse matrix                                  !
  !============================================================================!

  subroutine sparse_embed_array_copy_matrix(dest, src)

     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: dest   ! The destination array of sparse matrices
     type(SPAM3_EMBED), intent(in) :: src             ! The source sparse matrix

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_embed_array_check(dest) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_copy_matrix: &
             &dest fails the sparse_array check.')
     end if

     ! Perform the copy.
     do kpoint = 1, dest%num_kpoints
        do ispin = 1, dest%num_spins
           call sparse_embed_copy( dest%m(ispin, kpoint), src )
        end do
     end do

  end subroutine sparse_embed_array_copy_matrix


!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs a matmul operation on the matrices in the arrays. !
  !   C := A . B                                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   arrayC  (inout) : The array of sparse matrices C                         !
  !   arrayA  (input) : The array of sparse matrices A                         !
  !   arrayB  (input) : The array of sparse matrices B                         !
  !============================================================================!

  subroutine sparse_embed_array_product_array_array(arrayC, arrayA, arrayB)

     use utils, only: utils_abort

     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: arrayC ! The array of sparse matrices C
     type(SPAM3_EMBED_ARRAY), intent(in)    :: arrayA ! The array of sparse matrices A
     type(SPAM3_EMBED_ARRAY), intent(in)    :: arrayB ! The array of sparse matrices B

     ! Local variables
     integer :: ispin, kpoint                     ! Iterators
     logical :: loc_iscmplx ! Local copy of optional iscmplx argument
     integer :: ierr        ! Error flag
     integer :: iter
     type(SPAM3_EMBED_ARRAY) :: tmp_array

     ! Check arrays of sparse matrices.
     if ( sparse_embed_array_check(arrayA) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_product_array_array: &
             &arrayA fails the sparse_array check.')
     else if ( sparse_embed_array_check(arrayB) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_product_array_array: &
             &arrayB fails the sparse_array check.')
     else if ( sparse_embed_array_check(arrayC) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_product_array_array: &
             &arrayC fails the sparse_array check.')
     end if

     ! Check compatibility of numbers of kpoints.
     if ( .not. ( ( (arrayA%num_kpoints == arrayB%num_kpoints) .and.  &
                    (arrayA%num_kpoints == arrayC%num_kpoints) ) .or. &
                  ( ( arrayA%num_kpoints * arrayB%num_kpoints == arrayC%num_kpoints ) .and.  &
                    ( (arrayA%num_kpoints == 1) .or. (arrayB%num_kpoints == 1) ) ) ) ) then
        call utils_abort('Error in sparse_embed_array_product_array_array: &
             &cannot guess how to match kpoints.')
     end if

     ! Check compatibility of spin dimensions.
     if ( .not. ( ( (arrayA%num_spins == arrayB%num_spins) .and.  &
                    (arrayA%num_spins == arrayC%num_spins) ) .or. &
                  ( ( arrayA%num_spins * arrayB%num_spins == arrayC%num_spins ) ) ) ) then
        call utils_abort('Error in sparse_embed_array_(product_array_array: &
             &cannot guess how to match spin dimensions.')
     end if

     ! Multiply the matrices
     do kpoint = 1, arrayC%num_kpoints
        do ispin = 1, arrayC%num_spins
           call sparse_embed_product(arrayC%m(ispin, kpoint), &
                arrayA%m(min(ispin,arrayA%num_spins), min(kpoint,arrayA%num_kpoints)), &
                arrayB%m(min(ispin,arrayB%num_spins), min(kpoint,arrayB%num_kpoints)))
        end do
     end do

  end subroutine sparse_embed_array_product_array_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs a matmul operation on the matrices in the arrays. !
  !   C := A . B                                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   arrayC  (inout) : The array of sparse matrices C                         !
  !   amat    (input) : A sparse matrix A                                      !
  !   arrayB  (input) : The array of sparse matrices B                         !
  !============================================================================!

  subroutine sparse_embed_array_product_matrix_array(arrayC, amat, arrayB)

     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: arrayC ! The array of sparse matrices C
     type(SPAM3_EMBED), intent(in)          :: amat   ! The sparse matrix A
     type(SPAM3_EMBED_ARRAY), intent(in)    :: arrayB ! The array of sparse matrices B

     ! Local variables
     integer :: ispin, kpoint                     ! Iterators

     ! Check arrays of sparse matrices.
     if ( sparse_embed_array_check(arrayB) /= 1 ) then
        call utils_abort('Error in sparse_array_product_matrix_array: &
             &arrayB fails the sparse_array check.')
     else if ( sparse_embed_array_check(arrayC) /= 1 ) then
        call utils_abort('Error in sparse_array_product_matrix_array: &
             &arrayC fails the sparse_array check.')
     end if

     ! Check compatibility of shapes of arrays.
     !if ( utils_compare_vectors(shape(arrayB%m), shape(arrayC%m)) /= 0 ) then
     !   call utils_abort('Error in sparse_array_product_matrix_array: &
     !        &dimensions of arrayB and arrayC do not match.')
     !end if

     ! Multiply the matrices
     do kpoint = 1, arrayC%num_kpoints
        do ispin = 1, arrayC%num_spins
           call sparse_embed_product(arrayC%m(ispin,kpoint), amat, arrayB%m(ispin,kpoint) )
        end do
     end do

  end subroutine sparse_embed_array_product_matrix_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs a matmul operation on the matrices in the arrays. !
  !   C := A . B                                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   arrayC  (inout) : The array of sparse matrices C                         !
  !   arrayA  (input) : The array of sparse matrices A                         !
  !   bmat    (input) : A sparse matrix B                                      !
  !============================================================================!

  subroutine sparse_embed_array_product_array_matrix(arrayC, arrayA, bmat)

     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: arrayC ! The array of sparse matrices C
     type(SPAM3_EMBED_ARRAY), intent(in)    :: arrayA ! The array of sparse matrices A
     type(SPAM3_EMBED), intent(in)          :: bmat   ! The sparse matrix B

     ! Local variables
     integer :: ispin, kpoint                     ! Iterators

     ! Check arrays of sparse matrices.
     if ( sparse_embed_array_check(arrayA) /= 1 ) then
        call utils_abort('Error in sparse_array_product_array_matrix: &
             &arrayA fails the sparse_array check.')
     else if ( sparse_embed_array_check(arrayC) /= 1 ) then
        call utils_abort('Error in sparse_array_product_array_matrix: &
             &arrayC fails the sparse_array check.')
     end if

     ! Check compatibility of shapes of arrays.
     !if ( utils_compare_vectors(shape(arrayA%m), shape(arrayC%m)) /= 0 ) then
     !   call utils_abort('Error in sparse_array_product_array_matrix: &
     !        &dimensions of arrayA and arrayC do not match.')
     !end if

     ! Multiply the matrices
     do kpoint = 1, arrayC%num_kpoints
        do ispin = 1, arrayC%num_spins
           call sparse_embed_product(arrayC%m(ispin,kpoint), arrayA%m(ispin,kpoint), bmat)
        end do
     end do

  end subroutine sparse_embed_array_product_array_matrix

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs a trace on a sparse matrix or the product of two  !
  ! sparse matrices (the whole product is not calculated).                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (input) : The array of sparse matrices A                           !
  !   bmat  (input) : The compulsory array of sparse matrices B                !
  !                                                                            !
  ! Modified by Andrea Greco to allow operation on amat, May 2016.             !
  ! Currently supports only conjugacy of amat as operation to be performed.    !
  ! Converted to subroutine by Joseph Prentice, 29/01/2022.                    !
  !============================================================================!


  subroutine sparse_embed_array_trace_array_array(trace, amat, bmat, opA)

     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(in)    :: amat  ! The array of sparse matrices A
     type(SPAM3_EMBED_ARRAY), intent(in)    :: bmat  ! The array of sparse matrices B
     character(len=1), optional, intent(in) :: opA   ! Operation to be applied to A
     ! Result
     real(kind=DP), dimension(amat%num_spins,bmat%num_kpoints), intent(out) :: trace ! The result

     ! Local variables
     integer :: is, ik                ! Iterators
     integer :: ierr                  ! Error flag
     ! agrecocmplx
     character(len=1) :: loc_opA

     ! Check consistency of sparse arrays.
     if ( sparse_embed_array_check(amat) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_trace: &
             &amat fails the sparse_array check.')
     end if
     if ( sparse_embed_array_check(bmat) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_trace: &
             &bmat fails the sparse_array check.')
     end if

     ! Check that the dimensions of the sparse arrays match.
     if ( utils_compare_vectors(shape(amat%m), shape(bmat%m)) /= 0 ) then
        call utils_abort('Error in sparse_embed_array_trace: &
             &dimensions of amat and bmat do not match.')
     end if

     ! agrecocmplx
     if (present(opA)) then
        loc_opA = opA
     else
        loc_opA = 'N'
     end if

     ! Compute each of the traces.
     ! agrecocmplx: apply operation to amat if required
     do ik = 1, amat%num_kpoints
        do is = 1, amat%num_spins
           call sparse_embed_trace(trace(is, ik), amat%m(is,ik), bmat%m(is,ik), &
                opA=loc_opA)
        end do
     end do

  end subroutine sparse_embed_array_trace_array_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs a trace on a sparse matrix or the product of two  !
  ! sparse matrices (the whole product is not calculated).                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (input) : The array of sparse matrices A                           !
  !   bmat  (input) : The optional sparse matrix B                             !
  !                                                                            !
  ! Modified by Andrea Greco to allow operation on A, May 2016.                !
  ! Currently supports only conjugacy as operation to be performed.            !
  ! Converted to subroutine by Joseph Prentice, 29/01/2022.                    !
  !============================================================================!


  subroutine sparse_embed_array_trace_array_matrix(trace, amat, bmat, opA)

     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(in)          :: amat  ! The array of sparse matrices A
     type(SPAM3_EMBED), optional, intent(in)      :: bmat  ! The sparse matrix B
     ! agrecocmplx
     character(len=1), optional, intent(in) :: opA   ! Operation to be performed on A
     ! Result
     real(kind=DP), dimension(amat%num_spins,amat%num_kpoints), intent(out) :: trace  ! The result

     ! Local variables
     integer :: is, ik                ! Iterators
     integer :: ierr                  ! Error flag
     ! agrecocmplx
     character(len=1) :: loc_opA

     ! Check consistency of sparse arrays.
     if ( sparse_embed_array_check(amat) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_trace: &
             &amat fails the sparse_array check.')
     end if

     ! agrecocmplx
     if (present(opA)) then
        loc_opA = opA
     else
        loc_opA = 'N'
     end if

     ! Compute each of the traces.
     ! agrecocmplx: apply operation to amat if specified
     do ik = 1, amat%num_kpoints
        do is = 1, amat%num_spins
           call sparse_embed_trace(trace(is, ik), amat%m(is,ik), bmat, opA=loc_opA)
        end do
     end do

  end subroutine sparse_embed_array_trace_array_matrix

!------------------------------------------------------------------------------


  !============================================================================!
  ! This subroutine takes an array of sparse matrices and returns the list of  !
  ! their transposed matrices.  The order of the matrices in the array is      !
  ! left unchanged.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest  (inout)  : The resulting array of transposed matrices.             !
  !   src   (input)  : The array of matrices to be transposed.                 !
  !============================================================================!

  subroutine sparse_embed_array_transpose(dest, src)

     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(inout) :: dest    ! Resulting array of transposed matrices
     type(SPAM3_EMBED_ARRAY), intent(in)    :: src     ! Array of matrices to be transposed

     ! Local variables
     integer :: ispin, kpoint

     ! Check consistency of sparse arrays.
     if ( sparse_embed_array_check(src) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_transpose: &
             &src fails the sparse_array check.')
     else if ( sparse_embed_array_check(dest) /= 1 ) then
        call utils_abort('Error in sparse_embed_array_transpose: &
             &dest fails the sparse_array check.')
     end if

     ! Check that the dimensions match.
     if ( utils_compare_vectors(shape(src%m), shape(dest%m)) /= 0 ) then
        call utils_abort('Error in sparse_embed_array_transpose: &
             &mismatch of array dimensions.')
     end if

     ! Transpose the array element by element.
     do kpoint = 1, src%num_kpoints
        do ispin = 1, src%num_spins
           call sparse_embed_transpose(dest%m(ispin, kpoint), src%m(ispin, kpoint))
        end do
     end do

  end subroutine sparse_embed_array_transpose

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine runs the sparse expasion over an array of sparse matrices. !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat  (inout)  : The array of matrices to be expanded                     !
  !   pattern (in)  : pattern argument to be passed to sparse_expand           !
  !   sym  (in)     : Whether the matrices are symmetric.
  !============================================================================!

  !subroutine sparse_embed_array_expand(mat, pattern, sym)

  !   use sparse, only: sparse_expand
  !   use utils, only: utils_abort
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout) :: mat    ! The matrix
  !   integer,intent(in) :: pattern
  !   logical,optional,intent(in) :: sym

  !   ! Local variables
  !   integer :: ispin, kpoint

  !   ! Check consistency of the sparse array.
  !   if ( sparse_embed_array_check(mat) /= 1 ) then
  !      call utils_abort('Error in sparse_embed_array_expand: &
  !           &mat fails the sparse_array check.')
  !   end if

  !   ! Expand each of the matrices of the array.
  !   do kpoint = 1, mat%num_kpoints
  !      do ispin = 1, mat%num_spins
  !         call sparse_embed_expand(mat%m(ispin, kpoint), pattern, sym)
  !      end do
  !   end do

  !end subroutine sparse_embed_array_expand

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine obtains the maximum eigenvalues of an array of sparse      !
  ! matrices by an iterative conjugate gradients procedure.                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array (input)  : The array of matrices whose eigenvalues are desired     !
  !   met (input)    : The metric to use                                       !
  !   eval (output)  : The eigenvalue estimate                                 !
  !   tol (input)    : Optional tolerance for eigenvalue estimate              !
  !----------------------------------------------------------------------------!
  ! Written by J.M. Escartin, March 2015.
  !============================================================================!

  subroutine sparse_embed_array_extremal_eigenvalues(array, met, eval, tol)

     use utils, only: utils_abort
     implicit none

     ! Argument
     type(SPAM3_EMBED_ARRAY), intent(in)       :: array     ! The array of sparse matrices
     type(SPAM3_EMBED), intent(in)             :: met       ! The metric to be used
     real(kind=DP), intent(out)          :: eval(:,:) ! Array of eigenvalues
     real(kind=DP), optional, intent(in) :: tol       ! Optional tolerance

     ! Local variables
     integer :: ispin, kpoint

     ! Check size of array of eigenvalues
     if ( ( array%num_kpoints /= size(eval, dim=2) ) .or. &
          ( array%num_spins   /= size(eval, dim=1) ) ) then
        call utils_abort('Error in sparse_embed_array_extremal_eigenvalues: &
             &wrong sizes in eval array.')
     end if

     ! Compute each of the eigenvalues.
     do kpoint = 1, array%num_kpoints
        do ispin = 1, array%num_spins
           call sparse_embed_extremal_eigenvalue(array%m(ispin,kpoint), met, eval(ispin,kpoint), tol)
        end do
     end do

  end subroutine sparse_embed_array_extremal_eigenvalues

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs a conjugacy operation on the arrays of matrices:  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (inout) : The array of sparse matrices a                           !
  !   bmat  (inout) : if present, stores the conjugate of array a              !
  !                                                                            !
  ! Written by Andrea Greco in April 2016.                                     !
  !============================================================================!

  !subroutine sparse_embed_array_conjugate(amat, bmat)

  !   use sparse, only: sparse_conjugate
  !   use utils, only: utils_abort
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout)           :: amat    ! The array of sparse matrices a
  !   type(SPAM3_EMBED_ARRAY), optional, intent(inout) :: bmat    ! Optional array to store conjugate

  !   ! Local variables
  !   integer :: ispin, kpoint

  !   ! Check arrays of sparse matrices.
  !   if ( sparse_embed_array_check(amat) /= 1 ) then
  !      call utils_abort('Error in sparse_embed_array_conjugate: &
  !           &amat fails the sparse_array check.')
  !   end if

  !   ! the conjugate of amat is stored in bmat
  !   if (present(bmat)) then
  !      if ( sparse_embed_array_check(bmat) /= 1 ) then
  !         call utils_abort('Error in sparse_embed_array_conjugate: &
  !                 &bmat fails the sparse_array check.')
  !      end if

  !      ! check dimensions are compatible
  !      if (amat%num_kpoints /= bmat%num_kpoints) then
  !         call utils_abort('Error in sparse_embed_array_conjugate: &
  !                 &kpoints dimensions do not match')
  !      end if

  !      if (amat%num_spins /= bmat%num_spins) then
  !         call utils_abort('Error in sparse_embed_array_conjugate: &
  !                 &spins dimensions do not match')
  !      end if

  !      do kpoint = 1, amat%num_kpoints
  !         do ispin = 1, amat%num_spins
  !            call sparse_embed_conjugate(amat%m(ispin,kpoint),bmat%m(ispin,kpoint))
  !         end do
  !      end do

  !   ! Conjugate each of the input matrices directly
  !   else
  !      do kpoint = 1, amat%num_kpoints
  !         do ispin = 1, amat%num_spins
  !            call sparse_embed_conjugate(amat%m(ispin,kpoint))
  !         end do
  !      end do
  !   end if

  !end subroutine sparse_embed_array_conjugate


!-!-----------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine writes a SPAM3_EMBED_ARRAY object to a file.                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 19th September 2018.                           !
  ! Note that this routine is specifically intended to work with 1 subsystem;  !
  ! this is really just a wrapper for the regular sparse_write functionality.  !
  !============================================================================!

  subroutine sparse_embed_array_write(mat, filename, iunit)

    use sparse, only: sparse_copy, sparse_write
    use sparse_array, only: sparse_array_create, sparse_array_destroy

    implicit none

    ! Arguments
    type(SPAM3_EMBED_ARRAY), intent(in) :: mat       ! Array to write
    character(len=*),  intent(in) :: filename  ! Filename to use
    integer, optional, intent(in) :: iunit     ! I/O unit

    ! Local variables
    integer :: ik, is
    type(SPAM3_ARRAY) :: mat_array
    character(len=30) :: tmp_struc

    ! rc2013: work around for embedding matrix structures
    tmp_struc = trim(mat%structure)
    call sparse_array_create(mat_array, n_spins=mat%num_spins, &
         n_kpoints=mat%num_kpoints, structure=tmp_struc)
    do is = 1, mat%num_spins
       do ik = 1, mat%num_kpoints
          ! rc2013: copy the mat info across to a SPAM3_ARRAY structure
          call sparse_copy(mat_array%m(is,ik), mat%m(is,ik)%m(1,1))
       end do
    end do

    ! rc2013: now let's write the matrix to file
    call sparse_write(mat_array%m, filename, iunit)
    call sparse_array_destroy(mat_array)

  end subroutine sparse_embed_array_write

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine reads a SPAM3_EMBED_ARRAY object from a file.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (inout) : The sparse matrix to be read                          !
  !   filename (input) : The filename to use                                   !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 19th September 2018.                           !
  ! Note that this routine is specifically intended to work with 1 subsystem;  !
  ! this is really just a wrapper for the regular sparse_read functionality.   !
  !============================================================================!

  subroutine sparse_embed_array_read(mat, filename, iunit)

    use sparse, only: sparse_copy, sparse_read
    use sparse_array, only: sparse_array_create, sparse_array_destroy

    implicit none

    type(SPAM3_EMBED_ARRAY), intent(inout), target :: mat       ! Array to write
    character(len=*),  intent(in)    :: filename  ! Filename to use
    integer, optional, intent(inout)    :: iunit     ! I/O unit

    ! Local variables
    integer :: is, ik
    logical :: get_unit
    type(SPAM3_ARRAY) :: mat_array
    character(len=30) :: tmp_struc

    ! rc2013: work around for embedding matrix structures
    tmp_struc = trim(mat%structure)
    call sparse_array_create(mat_array, n_spins=mat%num_spins, &
         n_kpoints=mat%num_kpoints, structure=tmp_struc)

    ! rc2013: now let's write the matrix to file
    call sparse_read(mat_array%m, filename, iunit)

    ! rc2013: copy the mat info across to the SPAM3_EMBED_ARRAY structure
    do is = 1, mat%num_spins
       do ik = 1, mat%num_kpoints
          call sparse_copy(mat%m(is,ik)%m(1,1), mat_array%m(is,ik))
       end do
    end do
    call sparse_array_destroy(mat_array)

  end subroutine sparse_embed_array_read

!-!-----------------------------------------------------------------------------

  !!============================================================================!
  !! Wrapper for sparse_convert_unsegment_real                                  !
  !!----------------------------------------------------------------------------!
  !! Arguments:                                                                 !
  !!   mat      (input)   : The sparse array to be converted                    !
  !!   idxinfo  (output)  : The block-index in a non-segmented form             !
  !!   dmtx     (output)  : The actual non-zero matrix elements                 !
  !!============================================================================!

  !subroutine sparse_embed_array_convert_unsegment_real(mat, idxlen, dmtxlen, &
  !     idxinfo, dmtx)

  !   use sparse, only: sparse_convert_unsegment_real
  !   use utils, only: utils_assert
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(in)          :: mat
  !   integer, intent(out)                   :: idxlen
  !   integer, intent(out)                   :: dmtxlen
  !   integer, intent(inout), optional       :: idxinfo(:,:)
  !   real(kind=DP), intent(inout), optional :: dmtx(:,:)

  !   ! Local variables
  !   integer :: ik

  !   ! Check arguments
  !   call utils_assert(present(idxinfo).eqv.present(dmtx), &
  !        'Error in sparse_embed_array_convert_unsegment_real: &
  !        &only all/none optional arguments allowed.')

  !   if (present(dmtx)) then
  !      do ik = 1, mat%num_kpoints
  !         call sparse_embed_convert_unsegment_real(mat%m(:,ik), idxlen, &
  !              dmtxlen, idxinfo(:,ik), dmtx(:,ik))
  !      end do
  !   else
  !      call sparse_embed_convert_unsegment_real(mat%m(:,1), idxlen, dmtxlen)
  !   end if

  !end subroutine sparse_embed_array_convert_unsegment_real

!-!-----------------------------------------------------------------------------

  !!============================================================================!
  !! Wrapper for sparse_convert_segment_real                                    !
  !!----------------------------------------------------------------------------!
  !! Arguments:                                                                 !
  !!   mat      (output)  : The sparse array to be converted                    !
  !!   idxinfo  (input)   : The block-index in a non-segmented form             !
  !!   dmtx     (input)   : The actual non-zero matrix elements                 !
  !!============================================================================!

  !subroutine sparse_embed_array_convert_segment_real(mat, idxinfo, dmtx)

  !   use sparse, only: sparse_convert_segment_real
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout) :: mat        ! Matrix to fill
  !   integer, intent(in)        :: idxinfo(:,:)
  !   real(kind=DP), intent(in)  :: dmtx(:,:)

  !   ! Local variables
  !   integer :: ik

  !   do ik = 1, mat%num_kpoints
  !      call sparse_embed_convert_segment_real(mat%m(:,ik), idxinfo(:,ik), dmtx(:,ik))
  !   end do

  !end subroutine sparse_embed_array_convert_segment_real

!-!-----------------------------------------------------------------------------

  !!============================================================================!
  !! Wrapper for sparse_convert_unsegment_complex                               !
  !!----------------------------------------------------------------------------!
  !! Arguments:                                                                 !
  !!   mat      (input)   : The sparse array to be converted                    !
  !!   idxinfo  (output)  : The block-index in a non-segmented form             !
  !!   zmtx     (output)  : The actual non-zero matrix elements                 !
  !!============================================================================!

  !subroutine sparse_embed_array_convert_unsegment_complex(mat, idxlen, zmtxlen, &
  !     idxinfo, zmtx)

  !   use sparse, only: sparse_convert_unsegment_complex
  !   use utils, only: utils_assert
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(in)             :: mat
  !   integer, intent(out)                      :: idxlen
  !   integer, intent(out)                      :: zmtxlen
  !   integer, intent(inout), optional          :: idxinfo(:,:)
  !   complex(kind=DP), intent(inout), optional :: zmtx(:,:)

  !   ! Local variables
  !   integer :: ik

  !   ! Check arguments
  !   call utils_assert(present(idxinfo).eqv.present(zmtx), &
  !        'Error in sparse_embed_array_convert_unsegment_complex: &
  !        &only all/none optional arguments allowed.')

  !   if (present(zmtx)) then
  !      do ik = 1, mat%num_kpoints
  !         call sparse_embed_convert_unsegment_complex(mat%m(:,ik), idxlen, &
  !              zmtxlen, idxinfo(:,ik), zmtx(:,ik))
  !      end do
  !   else
  !      call sparse_embed_convert_unsegment_complex(mat%m(:,1), idxlen, zmtxlen)
  !   end if

  !end subroutine sparse_embed_array_convert_unsegment_complex

!-!-----------------------------------------------------------------------------

  !!============================================================================!
  !! Wrapper for sparse_convert_segment_complex                                 !
  !!----------------------------------------------------------------------------!
  !! Arguments:                                                                 !
  !!   mat      (output)  : The sparse array to be converted                    !
  !!   idxinfo  (input)   : The block-index in a non-segmented form             !
  !!   zmtx     (input)   : The actual non-zero matrix elements                 !
  !!============================================================================!

  !subroutine sparse_embed_array_convert_segment_complex(mat, idxinfo, zmtx)

  !   use sparse, only: sparse_convert_segment_complex
  !   implicit none

  !   ! Arguments
  !   type(SPAM3_EMBED_ARRAY), intent(inout) :: mat        ! Matrix to fill
  !   integer, intent(in)              :: idxinfo(:,:)
  !   complex(kind=DP), intent(in)     :: zmtx(:,:)

  !   ! Local variables
  !   integer :: ik

  !   do ik = 1, mat%num_kpoints
  !      call sparse_embed_convert_segment_complex(mat%m(:,ik), idxinfo(:,ik), &
  !           zmtx(:,ik))
  !   end do

  !end subroutine sparse_embed_array_convert_segment_complex

!------------------------------------------------------------------------------

  !============================================================================!
  ! This function checks the components of an array object and their           !
  ! consistency.                                                               !
  !----------------------------------------------------------------------------!
  ! Returned values are:                                                       !
  !   0     empty object                                                       !
  !   1     filled, consistent object                                          !
  !  -1,-2  otherwise (not properly allocated, inconsistent, etc.)             !
  !============================================================================!

  pure integer function sparse_embed_array_check(array) result(check)

     implicit none

     ! Arguments
     type(SPAM3_EMBED_ARRAY), intent(in) :: array

     if ( (array%num_kpoints == 0) .and. (array%num_spins == 0) .and.    &
          (.not.(allocated(array%m))) ) then
        check = 0
     else if ( (array%num_kpoints > 0) .and. (array%num_spins > 0) .and. &
          (allocated(array%m)) ) then
        if ( ( array%num_kpoints == size(array%m, dim=2) ) .and.         &
             ( array%num_spins   == size(array%m, dim=1) ) ) then
           check = 1
        else
           check = -1
        end if
     else
        check = -2
     end if

  end function sparse_embed_array_check

!------------------------------------------------------------------------------

  !============================================================================!
  ! This function returns true if every element in the matrix is nonzero, or   !
  ! false otherwise.                                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array   (input)  : The array to be assessed                              !
  !----------------------------------------------------------------------------!
  ! Written by J.M. Escartin, March 2015.                                      !
  ! Adapted for embedding by Robert Charlton, 15/08/2017.                      !
  !============================================================================!

  function sparse_embed_array_is_dense(array) result(is_dense)

     implicit none

     ! Argument
     type(SPAM3_EMBED_ARRAY), intent(in) :: array ! The array of sparse matrices

     ! Result
     logical :: is_dense

     is_dense = sparse_embed_is_dense(array%m(1,1))

  end function sparse_embed_array_is_dense

  ! EMBED_FIX!
  !============================================================================!
  ! This function returns the number of rows in an array of sparse matrices    !
  ! by consulting the relevant library entry.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array   (input)  : The array to be assessed                              !
  !----------------------------------------------------------------------------!
  ! Written by J.M. Escartin, March 2015.                                      !
  !============================================================================!

  function sparse_embed_array_num_rows(array) result(nrows)

     implicit none

     ! Argument
     type(SPAM3_EMBED_ARRAY), intent(in) :: array   ! The array of sparse matrices

     ! Result
!     real(DP), dimension(:,:), allocatable :: nrows
     real(kind=DP) :: nrows

     ! Local variables
!     integer :: ispin, kpoint
!     integer :: ierr

!     ! Allocate result.
!     allocate(nrows(array%num_spins, array%num_kpoints), stat=ierr)
!     call utils_alloc_check('sparse_array_num_rows', 'num_rows', ierr)

!     ! Compute each of the numbers of rows.
!     do kpoint = 1, array%num_kpoints
!        do ispin = 1, array%num_spins
!           nrows(ispin, kpoint) = sparse_num_rows(array%m(ispin,kpoint))
!        end do
!     end do
     ! Just grab the number of rows from the first matrix in the array.
     nrows = sparse_embed_num_rows(array%m(1,1))

  end function sparse_embed_array_num_rows

  !============================================================================!
  ! This function returns the number of non-zero elements in an embedded       !
  ! sparse matrix in a floating-point real (so as to avoid integer overflows). !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 25/07/2017.                                    !
  ! Based on sparse_num_element_mat by Nicholas Hine.                          !
  !============================================================================!

  function sparse_embed_array_num_element(mat) result(num_sum)

    use sparse, only: sparse_num_element

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat   ! The embedding sparse matrix

    ! Local variables
    integer :: isub, jsub, mrows, ncols

    ! Result
    real(kind=DP) :: num_sum

    ! rc2013: number of subsystems
    mrows = mat%mrows
    ncols = mat%ncols
    num_sum = 0.0_DP

    ! Very simple...
    do jsub=1,ncols
       do isub=1,mrows
          num_sum = num_sum + sparse_num_element(mat%m(isub,jsub))
       enddo
    enddo

  end function sparse_embed_array_num_element


!------------------------------------------------------------------------------

  !==========================================================================!
  ! This subroutine builds an embedding matrix from the elements of another  !
  ! SPAM3 matrix array.                                                      !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   dest  (inout) : The embedded sparse matrix to be built.                !
  !   src   (input) : The embedded sparse matrix to be copied.               !
  !   row   (input) : Index to extract row vector.                           !
  !   col   (input) : Index to extract column vector.                        !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 04/09/2017.                                  !
  !==========================================================================!

  subroutine sparse_embed_array_extract_sub(dest,src,row,col)

    use utils, only: utils_abort

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    type(SPAM3_EMBED_ARRAY), intent(inout) :: dest
    type(SPAM3_EMBED_ARRAY), intent(in)    :: src
    integer, optional, intent(in)          :: row ! Which row we want to extract
    integer, optional, intent(in)          :: col ! Which column we want to extract

    ! Local Variables
    integer :: is, ik

    ! rc2013: for now only do this if spin/k-points match
    if ( ( dest%num_kpoints == src%num_kpoints ) .and. &
       ( dest%num_spins == src%num_spins ) ) then
       ! rc2013: copy each k/spin block
       do ik = 1, dest%num_kpoints
          do is = 1, dest%num_spins
             call sparse_embed_extract_sub(dest%m(is,ik),src%m(is,ik), row, col)
          end do
       end do
    else
       call utils_abort('Error in sparse_embed_array_extract_sub: &
            &src and dest do not match.')
    endif

  end subroutine sparse_embed_array_extract_sub


!------------------------------------------------------------------------------

  !==========================================================================!
  ! This subroutine builds an embedding matrix from the elements of another  !
  ! SPAM3 matrix array.                                                      !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   dest  (inout) : The embedded sparse matrix to be built.                !
  !   src   (input) : The embedded sparse matrix to be copied.               !
  !   row   (input) : Index to extract row vector.                           !
  !   col   (input) : Index to extract column vector.                        !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 26/09/2017.                                  !
  !==========================================================================!

  subroutine sparse_embed_array_create_sub(dest,src,row,col,diag)

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    type(SPAM3_EMBED_ARRAY), intent(inout) :: dest
    type(SPAM3_EMBED_ARRAY), intent(in)    :: src
    integer, optional, intent(in)          :: row ! Which row we want to extract
    integer, optional, intent(in)          :: col ! Which column we want to extract
    logical, optional, intent(in)          :: diag

    ! Local Variables
    integer :: is, ik

    ! rc2013: allocate SPAM3_EMBED_ARRAY
    call sparse_embed_array_alloc(dest, spama=src)

    ! rc2013: for now only do this if spin/k-points match
    do ik = 1, dest%num_kpoints
       do is = 1, dest%num_spins
          call sparse_embed_create_sub(dest%m(is,ik),src%m(is,ik), row, col, diag)
       end do
    end do

  end subroutine sparse_embed_array_create_sub
!------------------------------------------------------------------------------

  !==========================================================================!
  ! This subroutine creates an embedding matrix from the elements of another !
  ! SPAM3 matrix array.                                                      !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   dest  (inout) : The embedded sparse matrix to be built.                !
  !   src   (input) : The embedded sparse matrix to be copied.               !
  !   row   (input) : Index to extract row vector.                           !
  !   col   (input) : Index to extract column vector.                        !
  !  diag   (input) : Flag for building the subsytem matrix dest out of the  !
  !                     diagonal components of src.                          !
  !--------------------------------------------------------------------------!
  ! Written by Robert Charlton, 26/09/2017.                                  !
  ! Should probably merge with regular sparse_embed_create.                  !
  !==========================================================================!

  subroutine sparse_embed_create_sub(dest,src,row,col,diag)

    use sparse, only: sparse_create
    use utils, only: utils_abort

    implicit none

    ! Arguments
    ! rc2013: embedding blocks included
    type(SPAM3_EMBED), intent(inout) :: dest
    type(SPAM3_EMBED), intent(in)    :: src
    integer, optional, intent(in)    :: row ! Which row we want to extract
    integer, optional, intent(in)    :: col ! Which column we want to extract
    logical, optional, intent(in)    :: diag

    ! Local Variables
    integer :: mrows, ncols
    integer :: isub
    integer, allocatable :: fake_rlib(:,:)
    logical :: loc_diag

    ! rc2013: set local version of diagonal flag
    if(present(diag)) then
       loc_diag = diag
    else
       loc_diag = .false.
    endif

    ! rc2013: pass info to dest based on inputs
    if(present(row) .and. present(col)) then
       ! We only want one 'element'
       call sparse_embed_alloc(dest, 1, 1, src, rlib=fake_rlib)
       call sparse_create(dest%m(1,1),src%m(row,col))

    else if(.not.present(row)) then
       mrows = src%mrows
       ncols = 1
       ! Extract column 'vector'
       call sparse_embed_alloc(dest, mrows, ncols, src, &
            rlib=fake_rlib)
       ! Extract the column vector specified by col
       do isub = 1,mrows
          ! rc2013: real test -- do we get the structures right?
          if(loc_diag) then
               call sparse_create(dest%m(isub,1),src%m(isub,isub),rlib=fake_rlib(isub,1))
          else
               call sparse_create(dest%m(isub,1),src%m(isub,col),rlib=fake_rlib(isub,1))
          endif
       end do

    else if(.not.present(col)) then
       ! Extract row 'vector'
       mrows = 1
       ncols = src%ncols
       call sparse_embed_alloc(dest, mrows, ncols, src, rlib=fake_rlib)
       ! Extract the row vector specified by row
       do isub = 1,ncols
          if(loc_diag) then
             call sparse_create(dest%m(1,isub),src%m(isub,isub),rlib=fake_rlib(1,isub))
          else
             call sparse_create(dest%m(1,isub),src%m(row,isub),rlib=fake_rlib(1,isub))
          endif
       end do
    else
       ! rc2013: nothing to be done here, pass the buck
       call utils_abort('Error in sparse_embed_extract_sub: called but &
            &no row/column information has been provided.')
       call sparse_embed_create(dest,src)
    end if

  end subroutine sparse_embed_create_sub

!------------------------------------------------------------------------------

  !==========================================================================!
  ! This function calculates the RMS element value of a embedding sparse     !
  ! matrix                                                                   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   mat  (input) : The matrix to be assessed                               !
  !--------------------------------------------------------------------------!
  ! Written by Joseph Prentice, May 2018.                                    !
  !==========================================================================!

  real(kind=DP) function sparse_embed_rms_element(mat)

    use sparse, only: sparse_rms_element,sparse_num_element

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat  ! The sparse embed matrix

    ! Local variables
    integer :: isub,jsub  ! Region loop counters
    real(kind=DP) :: rms_sum,el_sum,num_el  ! Summation and temporary variables

    rms_sum=0.d0
    el_sum=0.d0
    ! Loop over regions
    do jsub=1,mat%ncols
       do isub=1,mat%mrows

          ! Find the number of non-zero elements in that region
          num_el=sparse_num_element(mat%m(isub,jsub))
          ! Sum the total number of non-zero elements
          el_sum=el_sum+num_el
          ! Sum the total mod squared of each region's elements
          rms_sum=rms_sum+num_el*sparse_rms_element(mat%m(isub,jsub))**2

       end do
    end do

    ! Calculate full RMS value
    sparse_embed_rms_element = sqrt(rms_sum/el_sum)

  end function sparse_embed_rms_element

!------------------------------------------------------------------------------

  !==========================================================================!
  ! This function calculates the maximum absolute element value of a sparse  !
  ! embedded matrix                                                          !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   mat  (input) : The matrix to be assessed                               !
  !--------------------------------------------------------------------------!
  ! Written by Joseph Prentice, May 2018.                                    !
  !==========================================================================!

  real(kind=DP) function sparse_embed_max_abs_element(mat)

    use sparse, only: sparse_max_abs_element

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat  ! The sparse embed matrix

    ! Local variables
    integer :: isub,jsub  ! Region loop counters
    real(kind=DP) :: max_el,current_max_el    ! Temporary variables

    max_el=0.d0
    ! Loop over regions
    do jsub=1,mat%ncols
       do isub=1,mat%mrows

          ! Find the maximum element in that region
          current_max_el=sparse_max_abs_element(mat%m(isub,jsub))
          ! Store it if it is larger than the current largest element
          if (current_max_el.gt.max_el) max_el=current_max_el

       end do
    end do

    ! Return value
    sparse_embed_max_abs_element = max_el

  end function sparse_embed_max_abs_element

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3_EMBED into a     !
  ! dense full square real array                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block embedded sparse matrix                !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, August 2018.                                   !
  !============================================================================!

  subroutine sparse_embed_spam3tofull_real(dest,src)

    use sparse, only: sparse_num_rows, sparse_num_cols, sparse_convert
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: dest(:,:)     ! The full square matrix
    type(SPAM3_EMBED), intent(in) :: src        ! The block sparse matrix

    ! Local variables
    integer :: prev_rows, prev_cols
    integer :: num_rows, num_cols
    integer :: isub,jsub,ierr
    real(kind=DP), allocatable :: temp_mat(:,:)

    prev_rows=0
    ! jcap: loop over row blocks
    do isub=1,src%mrows
       num_rows=sparse_num_rows(src%m(isub,1))

       prev_cols=0
       ! jcap: loop over column blocks
       do jsub=1,src%ncols
          num_cols=sparse_num_cols(src%m(isub,jsub))

          allocate(temp_mat(num_rows,num_cols),stat=ierr)
          call utils_alloc_check('sparse_embed_spam3tofull_real','temp_mat',ierr)

          call sparse_convert(temp_mat,src%m(isub,jsub))
          dest((prev_rows+1):(prev_rows+num_rows),(prev_cols+1):(prev_rows+num_cols))=&
               temp_mat

          deallocate(temp_mat,stat=ierr)
          call utils_dealloc_check('sparse_embed_spam3tofull_real','temp_mat',ierr)

          ! jcap: keep track of number of cols already filled
          prev_cols=prev_cols+num_cols
       end do
       ! jcap: keep track of number of rows already filled
       prev_rows=prev_rows+num_rows
    end do

  end subroutine sparse_embed_spam3tofull_real

  !============================================================================!
  ! This subroutine converts a dense full square real array into a block       !
  ! sparse matrix in a SPAM3_EMBED                                             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, August 2018.                                   !
  !============================================================================!

  subroutine sparse_embed_fulltospam3_real(dest,src)

    use sparse, only: sparse_num_rows, sparse_num_cols, sparse_convert
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: dest ! The block sparse matrix
    real(kind=DP), intent(in) :: src(:,:)    ! The full square matrix

    ! Local variables
    integer :: prev_rows, prev_cols
    integer :: num_rows, num_cols
    integer :: isub,jsub,ierr
    real(kind=DP), allocatable :: temp_mat(:,:)

    prev_rows=0
    ! jcap: loop over row blocks
    do isub=1,dest%mrows
       num_rows=sparse_num_rows(dest%m(isub,1))

       prev_cols=0
       ! jcap: loop over column blocks
       do jsub=1,dest%ncols
          num_cols=sparse_num_cols(dest%m(isub,jsub))

          allocate(temp_mat(num_rows,num_cols),stat=ierr)
          call utils_alloc_check('sparse_embed_fulltospam3_real','temp_mat',ierr)

          temp_mat=src((prev_rows+1):(prev_rows+num_rows),&
               (prev_cols+1):(prev_rows+num_cols))
          call sparse_convert(dest%m(isub,jsub),temp_mat)

          deallocate(temp_mat,stat=ierr)
          call utils_dealloc_check('sparse_embed_fulltospam3_real','temp_mat',ierr)

          ! jcap: keep track of number of cols already filled
          prev_cols=prev_cols+num_cols
       end do
       ! jcap: keep track of number of rows already filled
       prev_rows=prev_rows+num_rows
    end do


  end subroutine sparse_embed_fulltospam3_real

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3_EMBED into a     !
  ! dense full square real array                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block embedded sparse matrix                !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, August 2018.                                   !
  !============================================================================!

  subroutine sparse_embed_spam3tofull_complex(dest,src)

    use sparse, only: sparse_num_rows, sparse_num_cols, sparse_convert
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: dest(:,:)     ! The full square matrix
    type(SPAM3_EMBED), intent(in) :: src        ! The block sparse matrix

    ! Local variables
    integer :: prev_rows, prev_cols
    integer :: num_rows, num_cols
    integer :: isub,jsub,ierr
    complex(kind=DP), allocatable :: temp_mat(:,:)

    prev_rows=0
    ! jcap: loop over row blocks
    do isub=1,src%mrows
       num_rows=sparse_num_rows(src%m(isub,1))

       prev_cols=0
       ! jcap: loop over column blocks
       do jsub=1,src%ncols
          num_cols=sparse_num_cols(src%m(isub,jsub))

          allocate(temp_mat(num_rows,num_cols),stat=ierr)
          call utils_alloc_check('sparse_embed_spam3tofull_complex','temp_mat',ierr)

          call sparse_convert(temp_mat,src%m(isub,jsub))
          dest((prev_rows+1):(prev_rows+num_rows),(prev_cols+1):(prev_rows+num_cols))=&
               temp_mat

          deallocate(temp_mat,stat=ierr)
          call utils_dealloc_check('sparse_embed_spam3tofull_complex','temp_mat',ierr)

          ! jcap: keep track of number of cols already filled
          prev_cols=prev_cols+num_cols
       end do
       ! jcap: keep track of number of rows already filled
       prev_rows=prev_rows+num_rows
    end do

  end subroutine sparse_embed_spam3tofull_complex

  !============================================================================!
  ! This subroutine converts a dense full square real array into a block       !
  ! sparse matrix in a SPAM3_EMBED                                             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, August 2018.                                   !
  !============================================================================!

  subroutine sparse_embed_fulltospam3_complex(dest,src)

    use sparse, only: sparse_num_rows, sparse_num_cols, sparse_convert
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: dest ! The block sparse matrix
    complex(kind=DP), intent(in) :: src(:,:)    ! The full square matrix

    ! Local variables
    integer :: prev_rows, prev_cols
    integer :: num_rows, num_cols
    integer :: isub,jsub,ierr
    complex(kind=DP), allocatable :: temp_mat(:,:)

    prev_rows=0
    ! jcap: loop over row blocks
    do isub=1,dest%mrows
       num_rows=sparse_num_rows(dest%m(isub,1))

       prev_cols=0
       ! jcap: loop over column blocks
       do jsub=1,dest%ncols
          num_cols=sparse_num_cols(dest%m(isub,jsub))

          allocate(temp_mat(num_rows,num_cols),stat=ierr)
          call utils_alloc_check('sparse_embed_fulltospam3_complex','temp_mat',ierr)

          temp_mat=src((prev_rows+1):(prev_rows+num_rows),&
               (prev_cols+1):(prev_rows+num_cols))
          call sparse_convert(dest%m(isub,jsub),temp_mat)

          deallocate(temp_mat,stat=ierr)
          call utils_dealloc_check('sparse_embed_fulltospam3_complex','temp_mat',ierr)

          ! jcap: keep track of number of cols already filled
          prev_cols=prev_cols+num_cols
       end do
       ! jcap: keep track of number of rows already filled
       prev_rows=prev_rows+num_rows
    end do


  end subroutine sparse_embed_fulltospam3_complex

!----------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine converts a SPAM3_EMBED_ARRAY into a SPAM3_ARRAY, given a   !
  ! particular region to extract                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination SPAM3_ARRAY                            !
  !   src    (input)  : The source SPAM3_EMBED_ARRAY                           !
  !   ireg   (input)  : The region to extract from the SPAM3_EMBED_ARRAY       !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, September 2018.                                !
  !============================================================================!

  subroutine sparse_embed_array_to_sparse_array(dest, src, ireg)

    use rundat, only: pub_num_spins, pub_num_kpoints
    use sparse, only: sparse_create, sparse_copy
    use sparse_array, only: SPAM3_ARRAY
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(SPAM3_ARRAY), intent(inout)    :: dest
    type(SPAM3_EMBED_ARRAY), intent(in) :: src
    integer, intent(in)                 :: ireg

    ! Local variables
    integer :: is,ik,ierr

    ! jcap: set number of kpoints and spins
    dest%num_spins=pub_num_spins
    dest%num_kpoints=pub_num_kpoints

    ! jcap: allocate SPAM3_ARRAY
    allocate(dest%m(dest%num_spins,dest%num_kpoints),stat=ierr)
    call utils_alloc_check('sparse_embed_array_to_sparse_array','dest%m',ierr)

    do is=1,dest%num_spins
       do ik=1,dest%num_kpoints
          call sparse_create(dest%m(is,ik),src%m(is,ik)%m(ireg,ireg))
          call sparse_copy(dest%m(is,ik),src%m(is,ik)%m(ireg,ireg))
       end do
    end do

  end subroutine sparse_embed_array_to_sparse_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sparse_embed_conjugate(amat, bmat)

    use sparse, only: sparse_conjugate
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout)  :: amat  ! The matrix to be operated on
    type(SPAM3_EMBED), optional, intent(inout)  :: bmat  ! Conjugate stored here if present

    ! Local variables
    integer :: isub,jsub

    ! This only makes sense if matrices are complex...
    call utils_assert(amat%iscmplx, &
         'Error in sparse_embed_conjugate: amat must be complex.')

    ! store the conjugate of amat in bmat
    if (present(bmat)) then
       call utils_assert(bmat%iscmplx, &
            'Error in sparse_embed_conjugate: bmat must be complex.')
       call utils_assert(amat%mrows==bmat%mrows, &
            'Error in sparse_embed_conjugate: amat and bmat must &
            &have the same number of embedding rows')
       call utils_assert(amat%ncols==bmat%ncols, &
            'Error in sparse_embed_conjugate: amat and bmat must &
            &have the same number of embedding columns')
       ! jcap: loop over regions
       do isub=1,amat%mrows
          do jsub=1,amat%ncols
             call sparse_conjugate(amat%m(isub,jsub),bmat%m(isub,jsub))
          end do
       end do
    ! modify amat directly
    else
       ! jcap: loop over regions
       do isub=1,amat%mrows
          do jsub=1,amat%ncols
             call sparse_conjugate(amat%m(isub,jsub))
          end do
       end do
    end if

  end subroutine sparse_embed_conjugate

!----------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine converts an array of SPAM3_EMBEDs into an array of SPAM3s, !
  ! given a particular section to extract                                      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination SPAM3s                                 !
  !   src    (input)  : The source SPAM3_EMBEDs                                !
  !   isub,jsub(input): The region to extract                                  !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, October 2018.                                  !
  !============================================================================!

  subroutine sparse_embed_extract_from_array(dst,src,isub,jsub)

    use sparse, only: sparse_create, sparse_copy
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout)    :: dst(:)
    type(SPAM3_EMBED), intent(in) :: src(:)
    integer, intent(in), optional :: isub,jsub

    ! Local variables
    integer :: size_src, size_dst, i, loc_isub, loc_jsub

    size_src=size(src)
    size_dst=size(dst)

    call utils_assert(size_src.eq.size_dst,'Error in sparse_embed_extract_from_array: '&
         &'sizes of arrays do not match!')

    loc_isub=1
    loc_jsub=1
    if (present(isub)) loc_isub=isub
    if (present(jsub)) loc_jsub=jsub

    do i=1,size_src
       call sparse_create(dst(i),src(i)%m(loc_isub,loc_jsub))
       call sparse_copy(dst(i),src(i)%m(loc_isub,loc_jsub))
    end do

  end subroutine sparse_embed_extract_from_array

!----------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine destroys an array of SPAM3s created by                     !
  ! sparse_embed_extract_from_array, after optionally copying them across      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination SPAM3_EMBEDs                           !
  !   src    (input)  : The source SPAM3s                                      !
  !   copy   (input)  : Whether to copy the data across                        !
  !   isub,jsub(input): The region to extract                                  !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, October 2018.                                  !
  !============================================================================!

  subroutine sparse_embed_destroy_extracted_array(src,dst,copy,isub,jsub)

    use sparse, only: sparse_destroy, sparse_copy
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout)       :: src(:)
    type(SPAM3_EMBED), intent(inout), optional :: dst(:)
    logical, intent(in), optional    :: copy
    integer, intent(in), optional    :: isub,jsub

    ! Local variables
    integer :: size_src, size_dst, i, loc_isub, loc_jsub
    logical :: loc_copy

    if (present(copy)) &
         call utils_assert(present(dst),&
         'Error in sparse_embed_destroy_extracted_array: incorrect combination of arguments')

    size_src=size(src)
    if (present(dst)) then
       size_dst=size(dst)
       call utils_assert(size_src.eq.size_dst,'Error in sparse_embed_extract_from_array: '&
            &'sizes of arrays do not match!')
    end if

    loc_copy=.false.
    if (present(copy)) loc_copy=copy

    if (loc_copy) then
       loc_isub=1
       loc_jsub=1
       if (present(isub)) loc_isub=isub
       if (present(jsub)) loc_jsub=jsub
    end if

    do i=1,size_src
       if (loc_copy) call sparse_copy(dst(i)%m(loc_isub,loc_jsub),src(i))
       call sparse_destroy(src(i))
    end do

  end subroutine sparse_embed_destroy_extracted_array


  !============================================================================!
  ! This function returns the index of the largest region in a SPAM3_EMBED     !
  ! structure.
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Robert Charlton, 16/03/2019.                                    !
  !============================================================================!

  function sparse_embed_largest_region(mat) result(region)

    use  sparse, only: sparse_num_rows, sparse_num_cols

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(in) :: mat   ! The embedding sparse matrix

    ! Local variables
    integer :: isub, mrows
    real(kind=DP) :: num, max_num

    ! Result
    integer :: region

    ! rc2013: only checks diagonal blocks
    mrows = mat%mrows
    max_num = 0.0_DP

    ! jcap: For single region calculations, the largest region is the
    ! only region
    if (mrows.eq.1) then
       region=1
    else
       ! rc2013: Add up the number of elements in each block individually
       do isub=1,mrows
       ! jcap: make this real to avoid integer overflow
          num = real(sparse_num_rows(mat%m(isub,isub)),DP) * &
               real(sparse_num_cols(mat%m(isub,isub)),DP)
          if(num .gt. max_num) then
             max_num = num
             region = isub
          end if
       end do
    end if

  end function sparse_embed_largest_region

end module sparse_embed
