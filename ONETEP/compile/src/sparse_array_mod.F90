! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

module sparse_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The module in this file was originally created by Jose M. Escartin in     !
!   Fall 2014, based on the sparse module.                                    !
!                                                                             !
!   TCM Group, Cavendish Laboratory, University of Cambridge                  !
!   19 JJ Thomson Av, Cambridge CB3 0HE, UK                                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use constants, only: DP

  use sparse, only: SPAM3
  implicit none

  private

  ! Public subroutines and functions

  ! sparse_array_base routines
  public :: sparse_array_create
  public :: sparse_array_destroy

  ! sparse_array_inquiry routines
  public :: sparse_array_is_dense
  public :: sparse_array_num_rows

  ! sparse_array_ops routines
  public :: sparse_array_copy
  public :: sparse_array_scale
  public :: sparse_array_axpy
  public :: sparse_array_product
  public :: sparse_array_trace
  public :: sparse_array_transpose
  public :: sparse_array_expand
  public :: sparse_array_extremal_eigenvalues
  ! agreco
  public :: sparse_array_conjugate

  ! sparse_array_utils routines
  public :: sparse_array_write
  public :: sparse_array_read
  public :: sparse_array_convert_unsegment_real
  public :: sparse_array_convert_segment_real
  public :: sparse_array_convert_unsegment_complex
  public :: sparse_array_convert_segment_complex
  public :: sparse_array_check

  ! create routines
  interface sparse_array_create
     module procedure sparse_array_create_array_array
     module procedure sparse_array_create_matrix_array
     module procedure sparse_array_create_array_matrix
     module procedure sparse_array_create_matrix_matrix
  end interface sparse_array_create

  ! AXPY routines
  interface sparse_array_axpy
     module procedure sparse_array_axpy_real_scalar
     module procedure sparse_array_axpy_real_array
     module procedure sparse_array_axpy_complex_scalar
     module procedure sparse_array_axpy_complex_array
     module procedure sparse_array_axpy_matrix_real_scalar
     module procedure sparse_array_axpy_matrix_real_array
     module procedure sparse_array_axpy_matrix_complex_scalar
     module procedure sparse_array_axpy_matrix_complex_array
  end interface sparse_array_axpy

  ! Scale and shift routines
  interface sparse_array_scale
     module procedure sparse_array_scale_real_scalar
     module procedure sparse_array_scale_real_array
     module procedure sparse_array_scale_complex_scalar
     module procedure sparse_array_scale_complex_array
  end interface sparse_array_scale

  ! Copy routines
  interface sparse_array_copy
     module procedure sparse_array_copy_array
     module procedure sparse_array_copy_matrix
  end interface sparse_array_copy

  ! Sparse product routines
  interface sparse_array_product
     module procedure sparse_array_product_array_array
     module procedure sparse_array_product_matrix_array
     module procedure sparse_array_product_array_matrix
  end interface sparse_array_product

  ! Sparse trace routines
  interface sparse_array_trace
     module procedure sparse_array_trace_array_array
     module procedure sparse_array_trace_array_matrix
  end interface

  ! Type definition for the array of matrix elements.
  type, public :: SPAM3_ARRAY
     integer :: num_spins = 0
     integer :: num_kpoints = 0
     type(SPAM3), allocatable :: m(:,:)
  end type SPAM3_ARRAY

contains


!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine creates a new array of sparse matrices by calling          !
  ! sparse_create for each of its elements.                                    !
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
  !============================================================================!

  subroutine sparse_array_create_array_array(new_spam, spama, spamb, &
       iscmplx, rlib, n_spins, n_kpoints, structure)

     use sparse, only: sparse_create
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(inout) :: new_spam      ! Sparse matrix array to be created
     type(SPAM3_ARRAY), intent(in) :: spama            ! First matrix array to use
     type(SPAM3_ARRAY), intent(in) :: spamb            ! Second matrix array to use
     logical, optional, intent(in) :: iscmplx          ! Flag for complex data (opt)
     integer, optional, intent(out) :: rlib            ! Optional library index
     integer, optional, intent(in) :: n_spins          ! Optional spin dimension.
     integer, optional, intent(in) :: n_kpoints        ! Optional number of kpoints.
     character(len=*), optional, intent(in) :: structure  ! Optional structure string.

     ! Local variables
     integer :: is, ik

     call sparse_array_alloc(new_spam, spama=spama, spamb=spamb, &
          n_spins=n_spins, n_kpoints=n_kpoints, structure=structure)

     ! Fill the array using sparse_create.
     do ik = 1, new_spam%num_kpoints
       do is = 1, new_spam%num_spins
          call sparse_create(new_spam%m(is,ik), spama%m(1,1), spamb%m(1,1), &
               iscmplx, rlib)
       end do
     end do

  end subroutine sparse_array_create_array_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine creates a new array of sparse matrices by calling          !
  ! sparse_create for each of its elements.                                    !
  !----------------------------------------------------------------------------!
  ! Required arguments:                                                        !
  !   new_spam (output) : The array of sparse matrices to be created           !
  !----------------------------------------------------------------------------!
  ! Optional arguments:                                                        !
  !   spama    (input)  : An array of sparse matrices whose structure is used  !
  !   spammatb (input)  : An optional matrix whose structure will be used.     !
  !   iscmplx  (input)  : A flag to indicate the data stored in the matrices   !
  !                       of this array is complex.                            !
  !   n_kpoints (input) : Force this number of kpoints in the array.           !
  !   n_spins   (input) : Force this spin dimension in the array.              !
  !   structure (input) : Force this structure in the new array.               !
  !============================================================================!

  subroutine sparse_array_create_array_matrix(new_spam, spama, spammatb, &
       iscmplx, rlib, n_spins, n_kpoints, structure)

     use sparse, only: sparse_create
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(inout) :: new_spam      ! Sparse matrix array to be created
     type(SPAM3_ARRAY), intent(in) :: spama            ! First matrix array to use
     type(SPAM3), intent(in), optional :: spammatb     ! (Second) matrix to use.
     logical, optional, intent(in) :: iscmplx          ! Flag for complex data (opt)
     integer, optional, intent(out) :: rlib            ! Optional library index
     integer, optional, intent(in) :: n_spins          ! Optional spin dimension.
     integer, optional, intent(in) :: n_kpoints        ! Optional number of kpoints.
     character(len=*), optional, intent(in) :: structure  ! Optional structure string.

     ! Local variables
     integer :: is, ik

     call sparse_array_alloc(new_spam, spama=spama, spammatb=spammatb, &
          n_spins=n_spins, n_kpoints=n_kpoints, structure=structure)

     ! Fill the array using sparse_create.
     do ik = 1, new_spam%num_kpoints
       do is = 1, new_spam%num_spins
          call sparse_create(new_spam%m(is,ik), spama%m(1,1), spammatb, &
               iscmplx, rlib)
       end do
     end do

  end subroutine sparse_array_create_array_matrix

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine creates a new array of sparse matrices by calling          !
  ! sparse_create for each of its elements.                                    !
  !----------------------------------------------------------------------------!
  ! Required arguments:                                                        !
  !   new_spam (output) : The array of sparse matrices to be created           !
  !----------------------------------------------------------------------------!
  ! Optional arguments:                                                        !
  !   spammata (input)  : A sparse matrix whose structure is used              !
  !   spamb    (input)  : An array of sparse matrices whose structure is used. !
  !   iscmplx  (input)  : A flag to indicate the data stored in the matrices   !
  !                       of this array is complex.                            !
  !   n_kpoints (input) : Force this number of kpoints in the array.           !
  !   n_spins   (input) : Force this spin dimension in the array.              !
  !   structure (input) : Force this structure in the new array.               !
  !============================================================================!

  subroutine sparse_array_create_matrix_array(new_spam, spammata, spamb, &
       iscmplx, rlib, n_spins, n_kpoints, structure)

     use sparse, only: sparse_create
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(inout) :: new_spam      ! Sparse matrix array to be created
     type(SPAM3), intent(in) :: spammata               ! (First) matrix to use
     type(SPAM3_ARRAY), intent(in) :: spamb            ! (Second) matrix array to use.
     logical, optional, intent(in) :: iscmplx          ! Flag for complex data (opt)
     integer, optional, intent(out) :: rlib            ! Optional library index
     integer, optional, intent(in) :: n_spins          ! Optional spin dimension.
     integer, optional, intent(in) :: n_kpoints        ! Optional number of kpoints.
     character(len=*), optional, intent(in) :: structure  ! Optional structure string.

     ! Local variables
     integer :: is, ik

     call sparse_array_alloc(new_spam, spammata=spammata, spamb=spamb, &
          n_spins=n_spins, n_kpoints=n_kpoints, structure=structure)

     ! Fill the array using sparse_create.
     do ik = 1, new_spam%num_kpoints
       do is = 1, new_spam%num_spins
          call sparse_create(new_spam%m(is,ik), spammata, spamb%m(1,1), &
               iscmplx, rlib)
       end do
     end do

  end subroutine sparse_array_create_matrix_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine creates a new array of sparse matrices by calling          !
  ! sparse_create for each of its elements.                                    !
  !----------------------------------------------------------------------------!
  ! Required arguments:                                                        !
  !   new_spam (output) : The array of sparse matrices to be created           !
  !----------------------------------------------------------------------------!
  ! Optional arguments:                                                        !
  !   spammata (input)  : An optional sparse matrix whose structure            !
  !                       will be used                                         !
  !   spammatb (input)  : An optional sparse matrix whose structure            !
  !                       will be used                                         !
  !   iscmplx  (input)  : A flag to indicate the data stored in the matrices   !
  !                       of this array is complex.                            !
  !   n_kpoints (input) : Number of kpoints in the array.                      !
  !   n_spins   (input) : Spin dimension of the array.                         !
  !   structure (input) : Force this structure in the new array.               !
  !============================================================================!

  subroutine sparse_array_create_matrix_matrix(new_spam, spammata, spammatb, &
       iscmplx, rlib, n_spins, n_kpoints, structure)

     use sparse, only: sparse_create
     use utils, only: utils_assert
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(inout) :: new_spam      ! Sparse matrix array to be created
     type(SPAM3), intent(in), optional :: spammata     ! First matrix to use
     type(SPAM3), intent(in), optional :: spammatb     ! Second matrix to use.
     logical, optional, intent(in) :: iscmplx          ! Flag for complex data (opt)
     integer, optional, intent(out) :: rlib            ! Optional library index
     integer, intent(in) :: n_spins                    ! Spin dimension.
     integer, intent(in) :: n_kpoints                  ! Number of kpoints.
     character(len=*), optional, intent(in) :: structure  ! Optional structure string.

     ! Local variables
     integer :: is, ik

     call sparse_array_alloc(new_spam, spammata=spammata, spammatb=spammatb, &
          n_spins=n_spins, n_kpoints=n_kpoints, structure=structure)

     ! Fill the array using sparse_create.
     do ik = 1, new_spam%num_kpoints
       do is = 1, new_spam%num_spins
          call sparse_create(new_spam%m(is,ik), spammata, spammatb, &
               iscmplx, rlib)
       end do
     end do

  end subroutine sparse_array_create_matrix_matrix

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

  subroutine sparse_array_alloc(new_spam, spama, spamb, spammata, spammatb, &
       n_spins, n_kpoints, structure)

     use utils, only: utils_abort, utils_alloc_check, utils_assert
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(inout) :: new_spam      ! Sparse matrix array to be created
     type(SPAM3_ARRAY), optional, intent(in) :: spama  ! First (optional) matrix array to use
     type(SPAM3_ARRAY), optional, intent(in) :: spamb  ! Other (optional) matrix array to use
     type(SPAM3), optional, intent(in) :: spammata     ! Optional matrix to use (instead of spama)
     type(SPAM3), optional, intent(in) :: spammatb     ! Optional matrix to use (instead of spamb)
     integer, optional, intent(in) :: n_spins          ! Optional spin dimension.
     integer, optional, intent(in) :: n_kpoints        ! Optional number of kpoints.
     character(len=*), optional, intent(in) :: structure  ! Optional structure string.

     ! Local variables
     integer :: ierr        ! Error flag

     ! Check combination of optional sparse array arguments.
     if ( (present(spamb) .or. present(spammatb)) .and. &
          (.not. ( present(spama) .or. present(spammata) ) ) ) then
        call utils_abort('Error in sparse_array_alloc: &
             &spamb or spammatb present but spama and spammata missing.')
     end if
     if ( present(spama) .and. present(spammata) ) then
        call utils_abort('Error in sparse_array_alloc: &
             &you cannot use spammata and spama simultaneously.')
     end if
     if ( present(spamb) .and. present(spammatb) ) then
        call utils_abort('Error in sparse_array_alloc: &
             &you cannot use spammatb and spamb simultaneously.')
     end if

     ! Check consistency of optional sparse array arguments.
     if (present(spama)) then
        if ( sparse_array_check(spama) /= 1 ) then
           call utils_abort('Error in sparse_array_alloc: &
                &spama fails the sparse_array check.')
        end if
     end if
     if (present(spamb)) then
        if ( sparse_array_check(spamb) /= 1 ) then
           call utils_abort('Error in sparse_array_alloc: &
                &spamb fails the sparse_array check.')
        end if
     end if

     ! Check that the new sparse array is really empty.
     ierr = sparse_array_check(new_spam)
     call utils_assert(ierr == 0, 'Error in sparse_array_alloc: &
          &new_spam not empty.  sparse_array_check returned: ', ierr)

     ! Set the number of kpoints of the new object.
     if ( present(n_kpoints) ) then
        new_spam%num_kpoints = n_kpoints
     else  if ( (present(spama)) .and. (present(spamb)) ) then
        if ( spama%num_kpoints == spamb%num_kpoints ) then
           new_spam%num_kpoints = spama%num_kpoints
        else if ( ( spama%num_kpoints == 1 ) .or. (spamb%num_kpoints == 1 ) ) then
           new_spam%num_kpoints = spama%num_kpoints * spamb%num_kpoints
        else
           call utils_abort('Error in sparse_array_alloc: &
                &mismatch in number of kpoints of spama and spamb')
        end if
     else if ( present(spama) ) then
        new_spam%num_kpoints = spama%num_kpoints
     else if ( present(spamb) ) then
        new_spam%num_kpoints = spamb%num_kpoints
     else
        call utils_abort('Error in sparse_array_alloc: &
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
           call utils_abort('Error in sparse_array_alloc: &
                &mismatch in number of spins of spama and spamb')
        end if
     else if ( present(spama) ) then
        new_spam%num_spins = spama%num_spins
     else if ( present(spamb) ) then
        new_spam%num_spins = spamb%num_spins
     else
        call utils_abort('Error in sparse_array_alloc: &
             &undefined number of spins to allocate')
     end if

     ! Allocate the array element.
     allocate(new_spam%m(new_spam%num_spins, new_spam%num_kpoints), stat=ierr)
     ! Tag differs from variable so that it matches the tag in
     ! sparse_array_destroy
     call utils_alloc_check('sparse_array_alloc', 'array%m', ierr)

     ! Initialise structures
     if (present(structure)) new_spam%m(:,:)%structure = structure

  end subroutine sparse_array_alloc

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine destroys an array of sparse matrices, freeing up           !
  ! the memory allocated to internal arrays.                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array (inout) : The array of sparse matrices to be destroyed.            !
  !============================================================================!

  subroutine sparse_array_destroy(array)

     use sparse, only: sparse_destroy
     use utils, only: utils_dealloc_check
     implicit none

     ! Argument
     type(SPAM3_array), intent(inout) :: array

     ! Local variables
     integer :: is, ik
     integer :: ierr

     ! Destroy each of the SPAM3 objects.
     do ik = 1, array%num_kpoints
        do is = 1, array%num_spins
           call sparse_destroy(array%m(is,ik))
        end do
     end do

     ! Deallocate the array itself.
     deallocate(array%m, stat=ierr)
     call utils_dealloc_check('sparse_array_destroy', 'array%m', ierr)

     ! Reset the sizes.
     array%num_spins = 0
     array%num_kpoints = 0

  end subroutine sparse_array_destroy

!------------------------------------------------------------------------------

  !============================================================================!
  ! This function returns true if every element in the matrix is nonzero, or   !
  ! false otherwise.                                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array   (input)  : The array to be assessed                              !
  !----------------------------------------------------------------------------!
  ! Written by J.M. Escartin, March 2015.                                      !
  !============================================================================!

  function sparse_array_is_dense(array) result(is_dense)

     use sparse, only: sparse_is_dense
     implicit none

     ! Argument
     type(SPAM3_ARRAY), intent(in) :: array ! The array of sparse matrices

     ! Result
     logical :: is_dense

     is_dense = sparse_is_dense(array%m(1,1))

  end function sparse_array_is_dense

!------------------------------------------------------------------------------

  !============================================================================!
  ! This function returns the number of rows in an array of sparse matrices    !
  ! by consulting the relevant library entry.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array   (input)  : The array to be assessed                              !
  !----------------------------------------------------------------------------!
  ! Written by J.M. Escartin, March 2015.                                      !
  !============================================================================!

  function sparse_array_num_rows(array) result(nrows)

     use sparse, only: sparse_num_rows
     implicit none

     ! Argument
     type(SPAM3_ARRAY), intent(in) :: array   ! The array of sparse matrices

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
     nrows = sparse_num_rows(array%m(1,1))

  end function sparse_array_num_rows

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

  subroutine sparse_array_scale_real_scalar(array, alpha, beta)

     use sparse, only: sparse_scale
     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout)    :: array  ! The array to be operated on
     real(kind=DP), intent(in)           :: alpha  ! The scaling parameter
     real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

     ! Local variables
     integer :: ispin, kpoint

     ! Check consistency of sparse array argument.
     if ( sparse_array_check(array) /= 1 ) then
        call utils_abort('Error in sparse_array_scale_real_scalar: &
             &array fails the sparse_array check.')
     end if

     ! Scale each of the matrices.
     do kpoint = 1, array%num_kpoints
        do ispin = 1, array%num_spins
           call sparse_scale(array%m(ispin,kpoint), alpha, beta)
        end do
     end do

  end subroutine sparse_array_scale_real_scalar


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

  subroutine sparse_array_scale_real_array(array, alpha, beta)

     use sparse, only: sparse_scale
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout)                    :: array  ! The array to be operated on
     real(kind=DP), dimension(:,:), intent(in)           :: alpha  ! The scaling parameters
     real(kind=DP), dimension(:,:), optional, intent(in) :: beta   ! The shift parameters

     ! Local variables
     integer :: is, ik

     ! Check shapes
     if ( present(beta) ) then
        if ( utils_compare_vectors(shape(alpha), shape(beta)) /= 0 ) then
           call utils_abort('Error in sparse_array_scale_real_array: &
                &shapes of alpha and beta do not match.')
        end if
     end if
     if ( utils_compare_vectors(shape(array%m), shape(alpha)) /= 0 ) then
        call utils_abort('Error in sparse_array_scale_real_array: &
             &shapes of array%m and alpha do not match.')
     end if

     ! Scale each of the matrices.
     do ik = 1, array%num_kpoints
        do is = 1, array%num_spins
           if( present(beta) ) then
              call sparse_scale( array%m(is,ik), alpha(is,ik), beta(is,ik) )
           else
              call sparse_scale( array%m(is,ik), alpha(is,ik) )
           end if
        end do
     end do

  end subroutine sparse_array_scale_real_array

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

  subroutine sparse_array_axpy_real_scalar(ymat, xmat, alpha)

     use sparse, only: sparse_axpy
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: ymat    ! The array of sparse matrices y
     type(SPAM3_array), intent(in) :: xmat       ! The array of sparse matrices x
     real(kind=DP), intent(in) :: alpha          ! The parameter alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(xmat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_real_array: &
             &xmat fails the sparse_array check.')
     else if ( sparse_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_real_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Check that shapes conform.
     if ( utils_compare_vectors(shape(xmat%m), shape(ymat%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_axpy_real_array: &
             &shapes of xmat%m and ymat%m do not match.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_axpy(ymat%m(ispin,kpoint), xmat%m(ispin,kpoint), alpha)
        end do
     end do

  end subroutine sparse_array_axpy_real_scalar

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

  subroutine sparse_array_axpy_real_array(ymat, xmat, alpha)

     use sparse, only: sparse_axpy
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: ymat           ! Array of sparse matrices y
     type(SPAM3_array), intent(in) :: xmat              ! Array of sparse matrices x
     real(kind=DP), dimension(:,:), intent(in) :: alpha ! Array of parameters alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(xmat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_real_array: &
             &xmat fails the sparse_array check.')
     else if ( sparse_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_real_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Check that shapes conform.
     if ( utils_compare_vectors(shape(xmat%m), shape(ymat%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_axpy_real_array: &
             &shapes of xmat%m and ymat%m do not match.')
     end if
     if ( utils_compare_vectors(shape(alpha), shape(xmat%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_axpy_real_array: &
             &shapes of alpha and xmat%m do not match.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_axpy(ymat%m(ispin,kpoint), xmat%m(ispin,kpoint), alpha(ispin,kpoint))
        end do
     end do

  end subroutine sparse_array_axpy_real_array

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

  subroutine sparse_array_axpy_matrix_real_scalar(ymat, xmat, alpha)

     use sparse, only: sparse_axpy
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: ymat    ! The array of sparse matrices y
     type(SPAM3), intent(in) :: xmat             ! The sparse matrix x
     real(kind=DP), intent(in) :: alpha          ! The parameter alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_matrix_real_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_axpy(ymat%m(ispin,kpoint), xmat, alpha)
        end do
     end do

  end subroutine sparse_array_axpy_matrix_real_scalar

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

  subroutine sparse_array_axpy_matrix_real_array(ymat, xmat, alpha)

     use sparse, only: sparse_axpy
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: ymat           ! Array of sparse matrices y
     type(SPAM3), intent(in) :: xmat                    ! Sparse matrix x
     real(kind=DP), dimension(:,:), intent(in) :: alpha ! Array of parameters alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_matrix_real_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Check that shapes conform.
     if ( utils_compare_vectors(shape(alpha), shape(ymat%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_axpy_matrix_real_array: &
             &shapes of alpha and ymat%m do not match.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_axpy(ymat%m(ispin,kpoint), xmat, alpha(ispin,kpoint))
        end do
     end do

  end subroutine sparse_array_axpy_matrix_real_array

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

  subroutine sparse_array_scale_complex_scalar(array, alpha, beta)

     use sparse, only: sparse_scale
     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout)    :: array  ! The array to be operated on
     complex(kind=DP), intent(in)        :: alpha  ! The scaling parameter
     real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

     ! Local variables
     integer :: ispin, kpoint

     ! Check consistency of sparse array argument.
     if ( sparse_array_check(array) /= 1 ) then
        call utils_abort('Error in sparse_array_scale_complex_scalar: &
             &array fails the sparse_array check.')
     end if

     ! Scale each of the matrices.
     do kpoint = 1, array%num_kpoints
        do ispin = 1, array%num_spins
           call sparse_scale(array%m(ispin,kpoint), alpha, beta)
        end do
     end do

  end subroutine sparse_array_scale_complex_scalar

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues               !
  ! of the matrices in an array of sparse matrices (parameters form arrays).   !
  ! The shapes of the array, alpha, and beta (if present) must conform.        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   array  (inout) : The array of sparse matrices to be rescaled and shifted !
  !   alpha  (input) : The complex array of scaling parameters.                !
  !   beta   (input) : The optional real array of scaling parameters.          !
  !============================================================================!

  subroutine sparse_array_scale_complex_array(array, alpha, beta)

     use sparse, only: sparse_scale
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout)                    :: array  ! The array to be operated on
     complex(kind=DP), dimension(:,:), intent(in)        :: alpha  ! The scaling parameters
     real(kind=DP), dimension(:,:), optional, intent(in) :: beta   ! The shift parameters

     ! Local variables
     integer :: ispin, kpoint

     ! Check shapes
     if ( present(beta) ) then
        if ( utils_compare_vectors(shape(alpha), shape(beta)) /= 0 ) then
           call utils_abort('Error in sparse_array_scale_complex_array: &
                &shapes of alpha and beta do not match.')
        end if
     end if
     if ( utils_compare_vectors(shape(array%m), shape(alpha)) /= 0 ) then
        call utils_abort('Error in sparse_array_scale_complex_array: &
             &shapes of array%m and alpha do not match.')
     end if

     ! Scale each of the matrices.
     do kpoint = 1, array%num_kpoints
        do ispin = 1, array%num_spins
           call sparse_scale( array%m(ispin,kpoint), alpha(ispin,kpoint), beta(ispin,kpoint))
        end do
     end do

  end subroutine sparse_array_scale_complex_array

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

  subroutine sparse_array_axpy_complex_scalar(ymat, xmat, alpha)

     use sparse, only: sparse_axpy
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: ymat    ! The array of sparse matrices y
     type(SPAM3_array), intent(in) :: xmat       ! The array of sparse matrices x
     complex(kind=DP), intent(in) :: alpha       ! The parameter alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(xmat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_complex_scalar: &
             &xmat fails the sparse_array check.')
     else if ( sparse_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_complex_scalar: &
             &ymat fails the sparse_array check.')
     end if

     ! Check that shapes conform.
     if ( utils_compare_vectors(shape(xmat%m), shape(ymat%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_axpy_complex_scalar: &
             &shapes of xmat%m and ymat%m do not match.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_axpy(ymat%m(ispin,kpoint), xmat%m(ispin,kpoint), alpha)
        end do
     end do

  end subroutine sparse_array_axpy_complex_scalar

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs an axpy operation on the arrays of matrices:      !
  !   y := y + alpha * x                                                       !
  ! The shapes of all three arrays must conform.                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The array of sparse matrices y                           !
  !   xmat  (input) : The array of sparse matrices x                           !
  !   alpha (input) : The array of complex parameters alpha                    !
  !============================================================================!

  subroutine sparse_array_axpy_complex_array(ymat, xmat, alpha)

     use sparse, only: sparse_axpy
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: ymat              ! Array of sparse matrices y
     type(SPAM3_array), intent(in) :: xmat                 ! Array of sparse matrices x
     complex(kind=DP), dimension(:,:), intent(in) :: alpha ! Array of parameters alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(xmat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_complex_array: &
             &xmat fails the sparse_array check.')
     else if ( sparse_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_complex_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Check that shapes conform.
     if ( utils_compare_vectors(shape(xmat%m), shape(ymat%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_axpy_complex_array: &
             &shapes of xmat%m and ymat%m do not match.')
     end if
     if ( utils_compare_vectors(shape(alpha), shape(xmat%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_axpy_complex_array: &
             &shapes of alpha and xmat%m do not match.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_axpy(ymat%m(ispin,kpoint), xmat%m(ispin,kpoint), alpha(ispin,kpoint))
        end do
     end do

  end subroutine sparse_array_axpy_complex_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs an axpy operation on the arrays of matrices:      !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The array of sparse matrices y                           !
  !   xmat  (input) : A sparse matrix x                                        !
  !   alpha (input) : The complex scalar parameter alpha                       !
  !============================================================================!

  subroutine sparse_array_axpy_matrix_complex_scalar(ymat, xmat, alpha)

     use sparse, only: sparse_axpy
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: ymat    ! The array of sparse matrices y
     type(SPAM3), intent(in) :: xmat             ! The sparse matrix x
     complex(kind=DP), intent(in) :: alpha       ! The parameter alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_matrix_complex_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_axpy(ymat%m(ispin,kpoint), xmat, alpha)
        end do
     end do

  end subroutine sparse_array_axpy_matrix_complex_scalar

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine performs an axpy operation on the arrays of matrices:      !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The array of sparse matrices y                           !
  !   xmat  (input) : The sparse matrix x                                      !
  !   alpha (input) : The array of complex parameters alpha                    !
  !============================================================================!

  subroutine sparse_array_axpy_matrix_complex_array(ymat, xmat, alpha)

     use sparse, only: sparse_axpy
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: ymat           ! Array of sparse matrices y
     type(SPAM3), intent(in) :: xmat                    ! Sparse matrix x
     complex(kind=DP), dimension(:,:), intent(in) :: alpha ! Array of parameters alpha

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(ymat) /= 1 ) then
        call utils_abort('Error in sparse_array_axpy_matrix_complex_array: &
             &ymat fails the sparse_array check.')
     end if

     ! Check that shapes conform.
     if ( utils_compare_vectors(shape(alpha), shape(ymat%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_axpy_matrix_complex_array: &
             &shapes of alpha and ymat%m do not match.')
     end if

     ! Axpy each of the matrices
     do kpoint = 1, ymat%num_kpoints
        do ispin = 1, ymat%num_spins
           call sparse_axpy(ymat%m(ispin,kpoint), xmat, alpha(ispin,kpoint))
        end do
     end do

  end subroutine sparse_array_axpy_matrix_complex_array

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

  subroutine sparse_array_copy_array(dest, src)

     use sparse, only: sparse_copy
     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: dest   ! The destination array of sparse matrices
     type(SPAM3_array), intent(in) :: src       ! The source array of sparse matrices

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(src) /= 1 ) then
        call utils_abort('Error in sparse_array_copy_array: &
             &src fails the sparse_array check.')
     else if ( sparse_array_check(dest) /= 1 ) then
        call utils_abort('Error in sparse_array_copy_array: &
             &dest fails the sparse_array check.')
     end if

     ! Perform the copy, guessing how to match matrices if needed and possible.
     if ( ( dest%num_kpoints == src%num_kpoints ) .and. &
          ( dest%num_spins == src%num_spins ) ) then
        ! Both arrays have the same shape.
        do kpoint = 1, src%num_kpoints
           do ispin = 1, src%num_spins
              call sparse_copy( dest%m(ispin, kpoint) , src%m(ispin, kpoint) )
           end do
        end do
     else if ( ( src%num_kpoints == 1 ) .and. ( src%num_spins == 1 ) ) then
        ! The source array is 1x1.
        do kpoint = 1, dest%num_kpoints
           do ispin = 1, dest%num_spins
              call sparse_copy( dest%m(ispin, kpoint), src%m(1,1) )
           end do
        end do
     else if ( ( src%num_kpoints == 1) .and. ( dest%num_spins == src%num_spins ) ) then
        ! The source array only has one kpoint
        do kpoint = 1, dest%num_kpoints
           do ispin = 1, src%num_spins
              call sparse_copy( dest%m(ispin, kpoint), src%m(ispin, 1) )
           end do
        end do
     else if ( ( src%num_kpoints == dest%num_spins) .and. ( src%num_spins == 1) ) then
        ! The source array only has one spin
        do kpoint = 1, src%num_kpoints
           do ispin = 1, dest%num_spins
              call sparse_copy( dest%m(ispin, kpoint), src%m(1, kpoint) )
           end do
        end do
     else if ( src%num_kpoints /= dest%num_kpoints ) then
        ! We don't know how to copy one array onto the other.
        call utils_abort('Error in sparse_array_copy_array: &
             &cannot guess how to match kpoints dimensions')
     else if ( src%num_spins /= dest%num_spins ) then
        ! We don't know how to copy one array onto the other.
        call utils_abort('Error in sparse_array_copy_array: &
             &cannot guess how to match spin dimensions')
     else
        ! We should never arrive here.
        call utils_abort('Undetermined error in sparse_array_copy_array.')
     end if

  end subroutine sparse_array_copy_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine copies the data from a sparse matrix to all matrices in an !
  ! array of sparse matrices.                                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest (inout) : The destination array of sparse matrices                  !
  !   src  (input) : The source sparse matrix                                  !
  !============================================================================!

  subroutine sparse_array_copy_matrix(dest, src)

     use sparse, only: sparse_copy
     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: dest   ! The destination array of sparse matrices
     type(SPAM3), intent(in) :: src             ! The source sparse matrix

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(dest) /= 1 ) then
        call utils_abort('Error in sparse_array_copy_matrix: &
             &dest fails the sparse_array check.')
     end if

     ! Perform the copy.
     do kpoint = 1, dest%num_kpoints
        do ispin = 1, dest%num_spins
           call sparse_copy( dest%m(ispin, kpoint), src )
        end do
     end do

  end subroutine sparse_array_copy_matrix

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

  subroutine sparse_array_product_array_array(arrayC, arrayA, arrayB, &
       allow_mix_types)

     use sparse, only: sparse_product
     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: arrayC   ! The array of sparse matrices C
     type(SPAM3_array), intent(in) :: arrayA      ! The array of sparse matrices A
     type(SPAM3_array), intent(in) :: arrayB      ! The array of sparse matrices B
     logical, intent(in), optional :: allow_mix_types

     ! Local variables
     integer :: ispin, kpoint                     ! Iterators

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(arrayA) /= 1 ) then
        call utils_abort('Error in sparse_array_product_array_array: &
             &arrayA fails the sparse_array check.')
     else if ( sparse_array_check(arrayB) /= 1 ) then
        call utils_abort('Error in sparse_array_product_array_array: &
             &arrayB fails the sparse_array check.')
     else if ( sparse_array_check(arrayC) /= 1 ) then
        call utils_abort('Error in sparse_array_product_array_array: &
             &arrayC fails the sparse_array check.')
     end if

     ! Check compatibility of numbers of kpoints.
     if ( .not. ( ( (arrayA%num_kpoints == arrayB%num_kpoints) .and.  &
                    (arrayA%num_kpoints == arrayC%num_kpoints) ) .or. &
                  ( ( arrayA%num_kpoints * arrayB%num_kpoints == arrayC%num_kpoints ) .and.  &
                    ( (arrayA%num_kpoints == 1) .or. (arrayB%num_kpoints == 1) ) ) ) ) then
        call utils_abort('Error in sparse_array_product_array_array: &
             &cannot guess how to match kpoints.')
     end if

     ! Check compatibility of spin dimensions.
     if ( .not. ( ( (arrayA%num_spins == arrayB%num_spins) .and.  &
                    (arrayA%num_spins == arrayC%num_spins) ) .or. &
                  ( ( arrayA%num_spins * arrayB%num_spins == arrayC%num_spins ) ) ) ) then
        call utils_abort('Error in sparse_array_product_array_array: &
             &cannot guess how to match spin dimensions.')
     end if

     ! Multiply the matrices
     do kpoint = 1, arrayC%num_kpoints
        do ispin = 1, arrayC%num_spins
           call sparse_product(arrayC%m(ispin, kpoint), &
                arrayA%m(min(ispin,arrayA%num_spins), &
                         min(kpoint,arrayA%num_kpoints)), &
                arrayB%m(min(ispin,arrayB%num_spins), &
                         min(kpoint,arrayB%num_kpoints)), &
                allow_mix_types)
        end do
     end do

  end subroutine sparse_array_product_array_array

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

  subroutine sparse_array_product_matrix_array(arrayC, amat, arrayB, &
       allow_mix_types)

     use sparse, only: sparse_product
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: arrayC   ! The array of sparse matrices C
     type(SPAM3), intent(in) :: amat              ! The sparse matrix A
     type(SPAM3_array), intent(in) :: arrayB      ! The array of sparse matrices B
     logical, intent(in), optional :: allow_mix_types

     ! Local variables
     integer :: ispin, kpoint                     ! Iterators

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(arrayB) /= 1 ) then
        call utils_abort('Error in sparse_array_product_matrix_array: &
             &arrayB fails the sparse_array check.')
     else if ( sparse_array_check(arrayC) /= 1 ) then
        call utils_abort('Error in sparse_array_product_matrix_array: &
             &arrayC fails the sparse_array check.')
     end if

     ! Check compatibility of shapes of arrays.
     if ( utils_compare_vectors(shape(arrayB%m), shape(arrayC%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_product_matrix_array: &
             &dimensions of arrayB and arrayC do not match.')
     end if

     ! Multiply the matrices
     do kpoint = 1, arrayC%num_kpoints
        do ispin = 1, arrayC%num_spins
           call sparse_product(arrayC%m(ispin,kpoint), &
                amat, arrayB%m(ispin,kpoint), allow_mix_types)
        end do
     end do

  end subroutine sparse_array_product_matrix_array

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

  subroutine sparse_array_product_array_matrix(arrayC, arrayA, bmat, &
       allow_mix_types)

     use sparse, only: sparse_product
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: arrayC   ! The array of sparse matrices C
     type(SPAM3_array), intent(in) :: arrayA      ! The array of sparse matrices A
     type(SPAM3), intent(in) :: bmat              ! The sparse matrix B
     logical, intent(in), optional :: allow_mix_types

     ! Local variables
     integer :: ispin, kpoint                     ! Iterators

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(arrayA) /= 1 ) then
        call utils_abort('Error in sparse_array_product_array_matrix: &
             &arrayA fails the sparse_array check.')
     else if ( sparse_array_check(arrayC) /= 1 ) then
        call utils_abort('Error in sparse_array_product_array_matrix: &
             &arrayC fails the sparse_array check.')
     end if

     ! Check compatibility of shapes of arrays.
     if ( utils_compare_vectors(shape(arrayA%m), shape(arrayC%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_product_array_matrix: &
             &dimensions of arrayA and arrayC do not match.')
     end if

     ! Multiply the matrices
     do kpoint = 1, arrayC%num_kpoints
        do ispin = 1, arrayC%num_spins
           call sparse_product(arrayC%m(ispin,kpoint), &
                arrayA%m(ispin,kpoint), bmat, allow_mix_types)
        end do
     end do

  end subroutine sparse_array_product_array_matrix

!------------------------------------------------------------------------------

  !============================================================================!
  ! This function performs a trace on a sparse matrix or the product of two    !
  ! sparse matrices (the whole product is not calculated).                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (input) : The array of sparse matrices A                           !
  !   bmat  (input) : The compulsory array of sparse matrices B                !
  !                                                                            !
  ! Modified by Andrea Greco to allow operation on amat, May 2016.             !
  ! Currently supports only conjugacy of amat as operation to be performed.    !
  !============================================================================!


  function sparse_array_trace_array_array(amat, bmat, opA) result(trace)

     use sparse, only: sparse_trace
     use utils, only: utils_abort, utils_alloc_check, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(in)          :: amat  ! The array of sparse matrices A
     type(SPAM3_array), intent(in)          :: bmat  ! The array of sparse matrices B
     character(len=1), optional, intent(in) :: opA   ! Operation to be applied to A

     ! Result
     real(DP), dimension(:,:), allocatable :: trace       ! The result

     ! Local variables
     integer :: is, ik                ! Iterators
     integer :: ierr                  ! Error flag
     ! agrecocmplx
     character(len=1) :: loc_opA

     ! Check consistency of sparse arrays.
     if ( sparse_array_check(amat) /= 1 ) then
        call utils_abort('Error in sparse_array_trace: &
             &amat fails the sparse_array check.')
     end if
     if ( sparse_array_check(bmat) /= 1 ) then
        call utils_abort('Error in sparse_array_trace: &
             &bmat fails the sparse_array check.')
     end if

     ! Check that the dimensions of the sparse arrays match.
     if ( utils_compare_vectors(shape(amat%m), shape(bmat%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_trace: &
             &dimensions of amat and bmat do not match.')
     end if

     ! Allocate result.
     allocate(trace(amat%num_spins, amat%num_kpoints), stat=ierr)
     !call utils_alloc_check('sparse_array_trace', 'trace', ierr)

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
           trace(is, ik) = sparse_trace(amat%m(is,ik), bmat%m(is,ik), opA=loc_opA)
        end do
     end do

  end function sparse_array_trace_array_array

!------------------------------------------------------------------------------

  !============================================================================!
  ! This function performs a trace on a sparse matrix or the product of two    !
  ! sparse matrices (the whole product is not calculated).                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (input) : The array of sparse matrices A                           !
  !   bmat  (input) : The optional sparse matrix B                             !
  !                                                                            !
  ! Modified by Andrea Greco to allow operation on A, May 2016.                !
  ! Currently supports only conjugacy as operation to be performed.            !
  !============================================================================!


  function sparse_array_trace_array_matrix(amat, bmat, opA) result(trace)

     use sparse, only: sparse_trace
     use utils, only: utils_abort, utils_alloc_check, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(in)          :: amat  ! The array of sparse matrices A
     type(SPAM3), optional, intent(in)      :: bmat  ! The sparse matrix B
     ! agrecocmplx
     character(len=1), optional, intent(in) :: opA   ! Operation to be performed on A

     ! Result
     real(DP), dimension(:,:), allocatable :: trace       ! The result

     ! Local variables
     integer :: is, ik                ! Iterators
     integer :: ierr                  ! Error flag
     ! agrecocmplx
     character(len=1) :: loc_opA

     ! Check consistency of sparse arrays.
     if ( sparse_array_check(amat) /= 1 ) then
        call utils_abort('Error in sparse_array_trace: &
             &amat fails the sparse_array check.')
     end if

     ! Allocate result.
     allocate(trace(amat%num_spins, amat%num_kpoints), stat=ierr)
     !call utils_alloc_check('sparse_array_trace', 'trace', ierr)

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
           trace(is, ik) = sparse_trace(amat%m(is,ik), bmat, opA=loc_opA)
        end do
     end do

  end function sparse_array_trace_array_matrix

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

  subroutine sparse_array_transpose(dest, src)

     use sparse, only: sparse_transpose
     use utils, only: utils_abort, utils_compare_vectors
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: dest    ! Resulting array of transposed matrices
     type(SPAM3_array), intent(in)    :: src     ! Array of matrices to be transposed

     ! Local variables
     integer :: ispin, kpoint

     ! Check consistency of sparse arrays.
     if ( sparse_array_check(src) /= 1 ) then
        call utils_abort('Error in sparse_array_transpose: &
             &src fails the sparse_array check.')
     else if ( sparse_array_check(dest) /= 1 ) then
        call utils_abort('Error in sparse_array_transpose: &
             &dest fails the sparse_array check.')
     end if

     ! Check that the dimensions match.
     if ( utils_compare_vectors(shape(src%m), shape(dest%m)) /= 0 ) then
        call utils_abort('Error in sparse_array_transpose: &
             &mismatch of array dimensions.')
     end if

     ! Transpose the array element by element.
     do kpoint = 1, src%num_kpoints
        do ispin = 1, src%num_spins
           call sparse_transpose(dest%m(ispin, kpoint), src%m(ispin, kpoint))
        end do
     end do

  end subroutine sparse_array_transpose

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine runs the sparse expasion over an array of sparse matrices. !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat  (inout)  : The array of matrices to be expanded                     !
  !   pattern (in)  : pattern argument to be passed to sparse_expand           !
  !   sym  (in)     : Whether the matrices are symmetric.
  !============================================================================!

  subroutine sparse_array_expand(mat, pattern, sym)

     use sparse, only: sparse_expand
     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout) :: mat    ! The matrix
     integer,intent(in) :: pattern
     logical,optional,intent(in) :: sym

     ! Local variables
     integer :: ispin, kpoint

     ! Check consistency of the sparse array.
     if ( sparse_array_check(mat) /= 1 ) then
        call utils_abort('Error in sparse_array_expand: &
             &mat fails the sparse_array check.')
     end if

     ! Expand each of the matrices of the array.
     do kpoint = 1, mat%num_kpoints
        do ispin = 1, mat%num_spins
           call sparse_expand(mat%m(ispin, kpoint), pattern, sym)
        end do
     end do

  end subroutine sparse_array_expand

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

  subroutine sparse_array_extremal_eigenvalues(array, met, eval, tol)

     use sparse, only: sparse_extremal_eigenvalue
     use utils, only: utils_abort
     implicit none

     ! Argument
     type(SPAM3_ARRAY), intent(in)       :: array     ! The array of sparse matrices
     type(SPAM3), intent(in)             :: met       ! The metric to be used
     real(kind=DP), intent(out)          :: eval(:,:) ! Array of eigenvalues
     real(kind=DP), optional, intent(in) :: tol       ! Optional tolerance

     ! Local variables
     integer :: ispin, kpoint

     ! Check size of array of eigenvalues
     if ( ( array%num_kpoints /= size(eval, dim=2) ) .or. &
          ( array%num_spins   /= size(eval, dim=1) ) ) then
        call utils_abort('Error in sparse_array_extremal_eigenvalues: &
             &wrong sizes in eval array.')
     end if

     ! Compute each of the eigenvalues.
     do kpoint = 1, array%num_kpoints
        do ispin = 1, array%num_spins
           call sparse_extremal_eigenvalue(array%m(ispin,kpoint), met, eval(ispin,kpoint), tol)
        end do
     end do

  end subroutine sparse_array_extremal_eigenvalues

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

  subroutine sparse_array_conjugate(amat, bmat)

     use sparse, only: sparse_conjugate
     use utils, only: utils_abort
     implicit none

     ! Arguments
     type(SPAM3_array), intent(inout)           :: amat    ! The array of sparse matrices a
     type(SPAM3_array), optional, intent(inout) :: bmat    ! Optional array to store conjugate

     ! Local variables
     integer :: ispin, kpoint

     ! Check arrays of sparse matrices.
     if ( sparse_array_check(amat) /= 1 ) then
        call utils_abort('Error in sparse_array_conjugate: &
             &amat fails the sparse_array check.')
     end if

     ! the conjugate of amat is stored in bmat
     if (present(bmat)) then
        if ( sparse_array_check(bmat) /= 1 ) then
           call utils_abort('Error in sparse_array_conjugate: &
                   &bmat fails the sparse_array check.')
        end if

        ! check dimensions are compatible
        if (amat%num_kpoints /= bmat%num_kpoints) then
           call utils_abort('Error in sparse_array_conjugate: &
                   &kpoints dimensions do not match')
        end if

        if (amat%num_spins /= bmat%num_spins) then
           call utils_abort('Error in sparse_array_conjugate: &
                   &spins dimensions do not match')
        end if

        do kpoint = 1, amat%num_kpoints
           do ispin = 1, amat%num_spins
              call sparse_conjugate(amat%m(ispin,kpoint),bmat%m(ispin,kpoint))
           end do
        end do

     ! Conjugate each of the input matrices directly
     else
        do kpoint = 1, amat%num_kpoints
           do ispin = 1, amat%num_spins
              call sparse_conjugate(amat%m(ispin,kpoint))
           end do
        end do
     end if

  end subroutine sparse_array_conjugate


!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine writes a SPAM3_ARRAY object to a file.                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  !============================================================================!

  subroutine sparse_array_write(mat, filename, iunit)

     use sparse, only: sparse_write
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(in) :: mat       ! Array to write
     character(len=*),  intent(in) :: filename  ! Filename to use
     integer, optional, intent(in) :: iunit     ! I/O unit

     call sparse_write(mat%m, filename, iunit)

  end subroutine sparse_array_write

!------------------------------------------------------------------------------

  !============================================================================!
  ! This subroutine reads a SPAM3_ARRAY object from a file.                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (inout) : The sparse matrix to be read                          !
  !   filename (input) : The filename to use                                   !
  !   iunit    (input) : The (optional) unit number to use if file is open     !
  !============================================================================!

  subroutine sparse_array_read(mat, filename, iunit)

     use sparse, only: sparse_read
     implicit none

     type(SPAM3_ARRAY), intent(inout), target :: mat       ! Array to write
     character(len=*),  intent(in)            :: filename  ! Filename to use
     integer, optional, intent(inout)         :: iunit     ! I/O unit

     call sparse_read(mat%m, filename, iunit)

  end subroutine sparse_array_read

!------------------------------------------------------------------------------

  !============================================================================!
  ! Wrapper for sparse_convert_unsegment_real                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input)   : The sparse array to be converted                    !
  !   idxinfo  (output)  : The block-index in a non-segmented form             !
  !   dmtx     (output)  : The actual non-zero matrix elements                 !
  !============================================================================!

  subroutine sparse_array_convert_unsegment_real(mat, idxlen, dmtxlen, &
       idxinfo, dmtx)

     use sparse, only: sparse_convert_unsegment_real
     use utils, only: utils_assert
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(in)          :: mat
     integer, intent(out)                   :: idxlen
     integer, intent(out)                   :: dmtxlen
     integer, intent(inout), optional       :: idxinfo(:,:)
     real(kind=DP), intent(inout), optional :: dmtx(:,:)

     ! Local variables
     integer :: ik

     ! Check arguments
     call utils_assert(present(idxinfo).eqv.present(dmtx), &
          'Error in sparse_array_convert_unsegment_real: &
          &only all/none optional arguments allowed.')

     if (present(dmtx)) then
        do ik = 1, mat%num_kpoints
           call sparse_convert_unsegment_real(mat%m(:,ik), idxlen, &
                dmtxlen, idxinfo(:,ik), dmtx(:,ik))
        end do
     else
        call sparse_convert_unsegment_real(mat%m(:,1), idxlen, dmtxlen)
     end if

  end subroutine sparse_array_convert_unsegment_real

!------------------------------------------------------------------------------

  !============================================================================!
  ! Wrapper for sparse_convert_segment_real                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (output)  : The sparse array to be converted                    !
  !   idxinfo  (input)   : The block-index in a non-segmented form             !
  !   dmtx     (input)   : The actual non-zero matrix elements                 !
  !============================================================================!

  subroutine sparse_array_convert_segment_real(mat, idxinfo, dmtx)

     use sparse, only: sparse_convert_segment_real
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(inout) :: mat        ! Matrix to fill
     integer, intent(in)        :: idxinfo(:,:)
     real(kind=DP), intent(in)  :: dmtx(:,:)

     ! Local variables
     integer :: ik

     do ik = 1, mat%num_kpoints
        call sparse_convert_segment_real(mat%m(:,ik), idxinfo(:,ik), dmtx(:,ik))
     end do

  end subroutine sparse_array_convert_segment_real

!------------------------------------------------------------------------------

  !============================================================================!
  ! Wrapper for sparse_convert_unsegment_complex                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input)   : The sparse array to be converted                    !
  !   idxinfo  (output)  : The block-index in a non-segmented form             !
  !   zmtx     (output)  : The actual non-zero matrix elements                 !
  !============================================================================!

  subroutine sparse_array_convert_unsegment_complex(mat, idxlen, zmtxlen, &
       idxinfo, zmtx)

     use sparse, only: sparse_convert_unsegment_complex
     use utils, only: utils_assert
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(in)             :: mat
     integer, intent(out)                      :: idxlen
     integer, intent(out)                      :: zmtxlen
     integer, intent(inout), optional          :: idxinfo(:,:)
     complex(kind=DP), intent(inout), optional :: zmtx(:,:)

     ! Local variables
     integer :: ik

     ! Check arguments
     call utils_assert(present(idxinfo).eqv.present(zmtx), &
          'Error in sparse_array_convert_unsegment_complex: &
          &only all/none optional arguments allowed.')

     if (present(zmtx)) then
        do ik = 1, mat%num_kpoints
           call sparse_convert_unsegment_complex(mat%m(:,ik), idxlen, &
                zmtxlen, idxinfo(:,ik), zmtx(:,ik))
        end do
     else
        call sparse_convert_unsegment_complex(mat%m(:,1), idxlen, zmtxlen)
     end if

  end subroutine sparse_array_convert_unsegment_complex

!------------------------------------------------------------------------------

  !============================================================================!
  ! Wrapper for sparse_convert_segment_complex                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (output)  : The sparse array to be converted                    !
  !   idxinfo  (input)   : The block-index in a non-segmented form             !
  !   zmtx     (input)   : The actual non-zero matrix elements                 !
  !============================================================================!

  subroutine sparse_array_convert_segment_complex(mat, idxinfo, zmtx)

     use sparse, only: sparse_convert_segment_complex
     implicit none

     ! Arguments
     type(SPAM3_ARRAY), intent(inout) :: mat        ! Matrix to fill
     integer, intent(in)              :: idxinfo(:,:)
     complex(kind=DP), intent(in)     :: zmtx(:,:)

     ! Local variables
     integer :: ik

     do ik = 1, mat%num_kpoints
        call sparse_convert_segment_complex(mat%m(:,ik), idxinfo(:,ik), &
             zmtx(:,ik))
     end do

  end subroutine sparse_array_convert_segment_complex

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

  pure integer function sparse_array_check(array) result(check)

     implicit none

     ! Arguments
     type(SPAM3_array), intent(in) :: array

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

  end function sparse_array_check

end module sparse_array
