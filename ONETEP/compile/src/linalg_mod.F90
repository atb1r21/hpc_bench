! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   This module contains routines for simple serial linear algebra
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Nicholas Hine, Robert Bell, David O'Regan,
!   Quintin Hill and Jose M. Escartin.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Name changed from wrappers to linalg by Valerio Vitale, 23/07/2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module linalg

  implicit none

  private

  public :: linalg_dot
  public :: linalg_ddot
  public :: linalg_zdotc

  public :: linalg_dsygv_lt
  public :: linalg_dsygv_lt_2
  public :: linalg_dcopy
  public :: linalg_invert_sym_matrix
  public :: linalg_invert_sym_cmatrix
  public :: linalg_dscal
  public :: linalg_daxpy
  public :: linalg_dgemm ! jd: actually not used anywhere
  public :: linalg_dgemm_serial_square ! jd
  public :: linalg_dgemm_serial        ! jd
  public :: linalg_zgemm_serial        ! jd
  public :: linalg_1d_fft
  public :: linalg_dsyev_lt ! lpl
  public :: linalg_zheev_lt
  public :: linalg_zhegv_lt
  public :: linalg_dgesv
  public :: linalg_dgelss

  public :: linalg_hdiag_serial
  public :: linalg_diag_serial
  public :: linalg_mat_mul_serial
  public :: linalg_sym_invert_serial
  public :: linalg_invert_serial

  interface linalg_dot
     module procedure linalg_ddot
     module procedure linalg_zdotc
  end interface linalg_dot

  interface linalg_hdiag_serial         ! hermitian matrix
     module procedure linalg_dsyev_lt   ! real, normal
     module procedure linalg_dsygv_lt   ! real, generalised
     module procedure linalg_zheev_lt   ! complex, normal
     module procedure linalg_zhegv_lt   ! complex, generalised
  end interface

  interface linalg_diag_serial          ! general matrix
!     module procedure linalg_dgeev_lt   ! real, normal
!     module procedure linalg_dggev_lt   ! real, generalised
      module procedure linalg_zgeev   ! complex, normal
!     module procedure linalg_zggev_lt   ! complex, generalised
  end interface

  interface linalg_mat_mul_serial
     ! matrix / matrix
     module procedure linalg_mat_mul_mm_dd ! real, real
     module procedure linalg_mat_mul_mm_zz ! complex, complex
  end interface

  interface linalg_invert_serial             ! general matrix
     module procedure linalg_invert_matrix   ! real
     module procedure linalg_invert_cmatrix  ! complex
  end interface

  interface linalg_sym_invert_serial             ! symmetric matrix
     module procedure linalg_invert_sym_matrix   ! real
     module procedure linalg_invert_sym_cmatrix  ! complex
  end interface



#ifdef MKL_FFTW3
#define FFTW3
#endif
#ifdef FFTW3_NO_OMP
#ifndef FFTW3
#define FFTW3
#endif
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dgemm(cc, & ! output
       aa, bb, num)               ! input

    !=========================================================================!
    ! This subroutine returns in cc the product C = A x B of square matrices  !
    ! A and B. The multiplication operation is parallelised. All processors   !
    ! hold copies of the whole matrices A and B and also C (on return from    !
    ! this subroutine).                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  cc (out)     : matrix product (matrix C)                               !
    !  aa  (input)  : A square matrix                                         !
    !  bb  (input)  : B square matrix                                         !
    !  num (input)  : dimension, all matrices are num x num square matrices   !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !  1) All matrices are square but not necessarily symmetric               !
    !                                                                         !
    !  2) The dummy array cc is overwritten during the calculation so it must !
    !     under no circustances point to the same memory as aa or bb          !
    !-------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 2/3/2004                            !
    !=========================================================================!

    use comms, only: comms_bcast, pub_my_proc_id, pub_total_num_procs
    use constants, only: DP
    use timer, only: timer_clock
    implicit none

    integer, intent(in) :: num
    real(kind=DP), intent(out):: cc(num, num) ! cks: not the same array as aa or bb!
    real(kind=DP), intent(in) :: aa(num, num)
    real(kind=DP), intent(in) :: bb(num, num)

    ! cks: <<local variables>>
    integer :: col_bb_batch_size
    integer :: col_bb_start
    integer :: col_bb_end
    integer :: local_num

    integer :: proc_count
    integer :: col_proc_start
    integer :: col_proc_end

    ! BLAS subroutine
    external :: dgemm

    call timer_clock('linalg_dgemm',1)

    ! cks: parallelisation of the dgemm - each proc
    ! cks: does a subset of the columns of bb
    col_bb_batch_size =num/pub_total_num_procs
    col_bb_start      =1 +pub_my_proc_id*col_bb_batch_size
    col_bb_end        =col_bb_start +col_bb_batch_size -1
    if (pub_my_proc_id .eq. (pub_total_num_procs-1) ) col_bb_end= num
    local_num         =col_bb_end -col_bb_start +1

    call timer_clock('linalg_dgemm_call_dgemm',1)

    ! cks: calculate and return aa x bb
    call dgemm('n', 'n', num, local_num, num, 1.0_DP, aa, num, &
         bb(:, col_bb_start: col_bb_end), num, 0.0_DP, &
         cc(:, col_bb_start: col_bb_end), num )

    call timer_clock('linalg_dgemm_call_dgemm',2)


    call timer_clock('linalg_dgemm_bcasts', 1)

    ! cks: every proc broadcasts the columns of cc it has calculated
    ! cks: to the other procs
    do proc_count =0, pub_total_num_procs -1

       col_proc_start =1 +proc_count*col_bb_batch_size
       col_proc_end   =col_proc_start +col_bb_batch_size -1
       if (proc_count .eq. (pub_total_num_procs-1) ) col_proc_end= num
       local_num =col_proc_end -col_proc_start +1

       call comms_bcast(proc_count, cc(:, &
            col_proc_start : col_proc_end), num*local_num)

    end do

    call timer_clock('linalg_dgemm_bcasts', 2)

    call timer_clock('linalg_dgemm',2)

  end subroutine linalg_dgemm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dgemm_serial_square(cc, & ! output
       aa, bb, num)                           ! input

    !=========================================================================!
    ! This subroutine returns in cc the product C = A x B of square matrices  !
    ! A and B. The multiplication operation is peformed in serial. This is    !
    ! convenient when one wants to multiply different matrices on different   !
    ! computational procs, e.g. atomblocks.                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  cc (out)     : matrix product (matrix C)                               !
    !  aa  (input)  : A square matrix                                         !
    !  bb  (input)  : B square matrix                                         !
    !  num (input)  : dimension, all matrices are num x num square matrices   !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !  1) All matrices are square but not necessarily symmetric.              !
    !                                                                         !
    !  2) The dummy array cc is overwritten during the calculation so it must !
    !     under no circustances point to the same memory as aa or bb.         !
    !-------------------------------------------------------------------------!
    ! Cloned by Jacek Dziedzic in May 2012 from wrappers_dgemm written by     !
    ! Chris-Kriton Skylaris on 2/3/2004 and stripped of parallelism.          !
    !=========================================================================!

    use constants, only: DP
    use timer, only: timer_clock
    implicit none

    integer, intent(in) :: num
    real(kind=DP), intent(out):: cc(num, num) ! cks: not the same array as aa or bb!
    real(kind=DP), intent(in) :: aa(num, num)
    real(kind=DP), intent(in) :: bb(num, num)

    ! BLAS subroutine
    external :: dgemm

    ! ------------------------------------------------------------------------
    call timer_clock('linalg_dgemm_serial_square',1)

    ! cks: calculate and return aa x bb
    call dgemm('n', 'n', num, num, num, 1.0_DP, aa, num, &
         bb, num, 0.0_DP, cc, num )

    call timer_clock('linalg_dgemm_serial_square',2)

  end subroutine linalg_dgemm_serial_square

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dgemm_serial(cc, & ! output
       aa, bb, m, n, k, lda, ldb, ldc) ! input

    !=========================================================================!
    ! This subroutine returns in cc the product C = A x B of matrices A and B.!
    ! The multiplication operation is peformed in serial.                     !
    ! Note that 'C' is *added to*, not overwritten.                           !
    ! The subroutine is not timed. It is expected that callers time it outside!
    ! of the call, only if needed.                                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  cc  (inout)     : matrix product (matrix C).                           !
    !  aa  (input)     : A matrix.                                            !
    !  bb  (input)     : B matrix.                                            !
    !  m, n, k (input) : matrix dinensions.                                   !
    !  Matrix dimensions -- follow dgemm(), and so:                           !
    !    aa(lda, k). Leading m-by-k matrix contains A on entry.               !
    !    bb(ldb, n). Leading k-by-n matrix contains B on entry.               !
    !    cc(ldc, n). Leading m-by-n matrix contains C on entry and on exit.   !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !  1) The dummy array cc is overwritten during the calculation so it must !
    !     under no circustances point to the same memory as aa or bb.         !
    !-------------------------------------------------------------------------!
    ! Cloned and generalized by Jacek Dziedzic in March 2021 from             !
    ! linalg_dgemm_serial, which then became linalg_dgemm_serial_square.      !
    !=========================================================================!

    use constants, only: DP
    implicit none

    integer, intent(in) :: m, n, k, lda, ldb, ldc
    real(kind=DP), intent(inout) :: cc(ldc, n)
    real(kind=DP), intent(in)    :: aa(lda, k)
    real(kind=DP), intent(in)    :: bb(ldb, n)

    ! BLAS subroutine
    external :: dgemm

    ! ------------------------------------------------------------------------

    ! cks: calculate and return aa x bb
    call dgemm('n', 'n', m, n, k, 1.0_DP, aa, lda, &
         bb, ldb, 1.0_DP, cc, ldc)

  end subroutine linalg_dgemm_serial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_zgemm_serial(cc, & ! output
       aa, bb, m, n, k, lda, ldb, ldc) ! input

    !=========================================================================!
    ! This subroutine returns in cc the product C = A x B of matrices A and B.!
    ! The multiplication operation is peformed in serial.                     !
    ! Note that 'C' is *added to*, not overwritten.                           !
    ! The subroutine is not timed. It is expected that callers time it outside!
    ! of the call, only if needed.                                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  cc  (inout)     : matrix product (matrix C).                           !
    !  aa  (input)     : A matrix.                                            !
    !  bb  (input)     : B matrix.                                            !
    !  m, n, k (input) : matrix dinensions.                                   !
    !  Matrix dimensions -- follow dgemm(), and so:                           !
    !    aa(lda, k). Leading m-by-k matrix contains A on entry.               !
    !    bb(ldb, n). Leading k-by-n matrix contains B on entry.               !
    !    cc(ldc, n). Leading m-by-n matrix contains C on entry and on exit.   !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !  1) The dummy array cc is overwritten during the calculation so it must !
    !     under no circustances point to the same memory as aa or bb.         !
    !-------------------------------------------------------------------------!
    ! Cloned by Jacek Dziedzic in March 2021 from linalg_dgemm_serial.        !
    !=========================================================================!

    use constants, only: DP
    implicit none

    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(kind=DP), intent(inout) :: cc(ldc, n)
    complex(kind=DP), intent(in)    :: aa(lda, k)
    complex(kind=DP), intent(in)    :: bb(ldb, n)

    ! BLAS subroutine
    external :: zgemm

    ! ------------------------------------------------------------------------

    ! cks: calculate and return aa x bb
    call zgemm('n', 'n', m, n, k, (1.0_DP, 0.0_DP), aa, lda, &
         bb, ldb, (1.0_DP, 0.0_DP), cc, ldc)

  end subroutine linalg_zgemm_serial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function linalg_ddot(n,x,incx,y,incy)

    !===========================================!
    ! ddot wrapper.                             !
    !-------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000. !
    !===========================================!

    use constants, only: DP
    implicit none

    real(kind=DP) :: linalg_ddot
    integer, intent(in) :: n, incx, incy
    real(kind=DP), intent(in), dimension(:) :: x, y

    real(kind=DP), external :: ddot

    linalg_ddot=ddot(n,x,incx,y,incy)

  end function linalg_ddot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function linalg_zdotc(n,x,incx,y,incy)

    !===========================================!
    ! zdotc wrapper.                            !
    !-------------------------------------------!
    ! Based on ddot wrapper by C.-K. Skylaris   !
    ! Written by J.M. Escartin in June 2015.    !
    !===========================================!

    use constants, only: DP
    implicit none

    complex(kind=DP) :: linalg_zdotc
    integer, intent(in) :: n, incx, incy
    complex(kind=DP), intent(in), dimension(:) :: x, y

    complex(kind=DP), external :: zdotc

    linalg_zdotc=zdotc(n,x,incx,y,incy)

  end function linalg_zdotc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dsygv_lt(eigenvectors,eigenvalues, &
       overlap_square,num,ut)

    !=======================================================!
    ! This subroutine solves the FC=SCe type of generalised !
    ! eigenvalue problem and obtains both eigenvalues and   !
    ! eigenvectors. F is stored as a lower triangle in a    !
    ! full square matrix. The eigenvalues are returned in   !
    ! ascending order.                                      !
    !-------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.             !
    !=======================================================!

    use constants, only: DP
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: num
    real(kind=DP), intent(inout) :: eigenvectors(:,:)
    real(kind=DP), intent(inout) :: overlap_square(:,:)
    real(kind=DP), intent(out) :: eigenvalues(:)
    logical, optional :: ut

    ! Local Variables
    real(kind=DP), dimension(:), allocatable :: work
    character :: lut
    integer :: info
    integer :: lda,ldb
    integer :: ierr ! error flag

    ! LAPACK subroutine
    external :: dsygv

    ! ndmh: check arguments
    if (size(eigenvalues)<num) then
       call utils_abort('Error in linalg_dsygv_lt: eigenvalues array is &
            &too small')
    end if
    if ((size(overlap_square,1)<num).or.(size(overlap_square,2)<num)) then
       call utils_abort('Error in linalg_dsygv_lt: overlap matrix is &
            &too small')
    end if
    if ((size(eigenvectors,1)<num).or.(size(eigenvectors,2)<num)) then
       call utils_abort('Error in linalg_dsygv_lt: eigenvectors matrix is &
            &too small')
    end if
    lut = 'L'
    if (present(ut)) then
       if (ut) lut = 'U'
    end if

    ! ndmh: find leading dimension of arrays
    lda = size(eigenvectors,1)
    ldb = size(overlap_square,1)

    allocate(work(3*num), stat=ierr)
    call utils_alloc_check('linalg_dsygv_lt','work',ierr)

    ! cks: we want to solve the FC=SCe type of generalised eigenvalue
    !      problem and obtain both eigenvalues and eigenvectors.
    !      F is stored as a lower triangle in a full square matrix.
    !      The eigenvalues are returned in ascending order.
#ifndef ESSL
    call dsygv(1,'V',lut,num,eigenvectors,lda,overlap_square,ldb,&
         eigenvalues,work,3*num,info)
#else
    info = 0
    call dsygv(1,eigenvectors,lda,overlap_square,ldb,&
         eigenvalues,eigenvectors,lda,num,work,3*num)
#endif

    call utils_assert(info == 0, 'Error in linalg_dsygv_lt: &
         &DSYGV returned info = ', info)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('linalg_dsygv_lt','work',ierr)

  end subroutine linalg_dsygv_lt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dsygv_lt_2(eigenvectors, eigenvalues, &
       overlap_square, num)

    !=======================================================!
    ! This subroutine solves the KSC=Cn type of generalised !
    ! eigenvalue problem and obtains both eigenvalues and   !
    ! eigenvectors. F is stored as a lower triangle in a    !
    ! full square matrix. The eigenvalues are returned in   !
    ! ascending order.                                      !
    !-------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.             !
    !=======================================================!

    use constants, only: DP
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: num
    real(kind=DP), intent(inout) :: eigenvectors(:,:)
    real(kind=DP), intent(inout) :: overlap_square(:,:)
    real(kind=DP), intent(out) :: eigenvalues(:)

    ! Local Variables
    real(kind=DP), dimension(:), allocatable :: work
    integer :: info
    integer :: lda,ldb
    integer :: ierr ! error flag

    ! LAPACK subroutine
    external :: dsygv

    call timer_clock('linalg_dsygv_lt_2', 1)

    ! ndmh: check arguments
    if (size(eigenvectors)<num) then
       call utils_abort('Error in linalg_dsygv_lt_2: eigenvectors array is &
            &too small')
    end if
    if ((size(overlap_square,1)<num).or.(size(overlap_square,2)<num)) then
       call utils_abort('Error in linalg_dsygv_lt_2: overlap matrix is &
            &too small')
    end if
    if ((size(eigenvectors,1)<num).or.(size(eigenvectors,2)<num)) then
       call utils_abort('Error in linalg_dsygv_lt_2: eigenvectors matrix is &
            &too small')
    end if

    ! ndmh: find leading dimension of arrays
    lda = size(eigenvectors,1)
    ldb = size(overlap_square,1)

    allocate(work(3*num),stat=ierr)
    call utils_alloc_check('linalg_dsygv_lt_2','work',ierr)

    ! cks: we want to solve the KSC=Cn type of generalised eigenvalue
    !      problem and obtain both eigenvalues and eigenvectors.
    !      K is stored as a lower triangle in a full square matrix.
    !      The eigenvalues are returned in ascending order.
    call dsygv(2,'V','L',num,eigenvectors,lda,overlap_square,ldb,&
         eigenvalues,work,3*num,info)

    call utils_assert(info == 0, 'Error in linalg_dsygv_lt_2: &
         &DSYGV returned info = ', info)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('linalg_dsygv_lt_2','work',ierr)

    call timer_clock('linalg_dsygv_lt_2', 2)

  end subroutine linalg_dsygv_lt_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_invert_sym_matrix(matrix,num,cond)

    !========================================================!
    ! This subroutine returns the inverse of a full square   !
    ! symmatric matrix.                                      !
    !--------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 17/4/2001          !
    ! Extended by Jacek Dziedzic on 10/05/2012 to optionally !
    ! return a condition number.                             !
    !========================================================!

    use constants, only: DP
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! jd: Arguments
    integer, intent(in) :: num
    real(kind=DP), intent(inout) :: matrix(:,:) ! inverse on exit
    real(kind=DP), intent(out), optional :: cond ! jd: Returned condition number

    ! cks: internal declarations
    integer :: work_length, info, row, col, lda
    integer, allocatable, dimension(:) :: ipiv
    real(kind=DP), allocatable, dimension(:) :: work_array
    integer :: ierr ! error flag
    integer, allocatable, dimension(:) :: iwork_array ! jd: For cond only
    real(kind=DP) :: anorm ! jd: 1-norm of input matrix
    real(kind=DP) :: rcond ! jd: Inverse of the condition number

    ! LAPACK function and subroutines
    real(kind=DP), external :: dlange ! jd: LAPACK auxiliary routine for 1-norms
    external :: dsytrf, dsycon, dsytri

    ! ------------------------------------------------------------------------

    call timer_clock('linalg_invert_sym_matrix',1)

    work_length=3*num
    lda = size(matrix,1)

    call utils_assert(lda >= num, 'Error in linalg_invert_sym_matrix(): &
         &invalid matrix sizes.')

    allocate(ipiv(num),stat=ierr)
    call utils_alloc_check('linalg_invert_sym_matrix','ipiv',ierr)
    allocate(work_array(work_length),stat=ierr)
    call utils_alloc_check('linalg_invert_sym_matrix','work_array',ierr)

    ! jd: If condition number is desired, calculate the 1-norm a priori,
    !     first allocating the extra workspace needed by DSYCON
    if(present(cond)) then
       allocate(iwork_array(num),stat=ierr)
       call utils_alloc_check('linalg_invert_sym_matrix','iwork_array',ierr)
       anorm = dlange('1', num, num, matrix, lda, work_array)
    end if

    ! cks: compute the factorization of a real symmetric matrix A using
    !      the Bunch-Kaufman diagonal pivoting method
    call dsytrf('L', num, matrix, lda, ipiv, work_array, work_length, info)

    call utils_assert(info == 0, 'Error in linalg_invert_sym_matrix(): &
         &DSYTRF returned info = ', info)

    ! jd: Calculate the inverse condition number, if requested
    if(present(cond)) then
       call dsycon('L', num, matrix, lda, ipiv, anorm, rcond, work_array, &
            iwork_array, info)

       call utils_assert(info == 0, 'Error in linalg_invert_sym_matrix(): &
            &DSYCON returned info = ', info)

       cond = 1.0_DP/rcond
    end if

    ! cks: compute the inverse of a real symmetric indefinite matrix A using
    ! cks: the factorization A = U*D*U**T or A = L*D*L**T computed by DSYTRF
    call dsytri('L', num, matrix, lda, ipiv, work_array, info )

    call utils_assert(info == 0, 'Error in linalg_invert_sym_matrix(): &
         &DSYTRI returned info = ', info)

    ! cks: fill in the upper triangular part of the symmetric matrix
    do row=1,num
       do col=1,row-1
          matrix(col,row)=matrix(row,col)
       enddo
    enddo

    if(present(cond)) then
       deallocate(iwork_array,stat=ierr)
       call utils_dealloc_check('linalg_invert_sym_matrix','iwork_array',ierr)
    end if
    deallocate(work_array,stat=ierr)
    call utils_dealloc_check('linalg_invert_sym_matrix','work_array',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('linalg_invert_sym_matrix','ipiv',ierr)

    call timer_clock('linalg_invert_sym_matrix',2)

  end subroutine linalg_invert_sym_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine linalg_invert_sym_cmatrix(cmatrix,num)

    !=======================================================!
    ! This subroutine returns the inverse of a full square  !
    ! complex indefinite symmetric matrix.                  !
    !-------------------------------------------------------!
    ! Written by David O'Regan in April 2009 based on       !
    ! linalg_invert_sym_matrix by Chris-Kriton Skylaris   !
    !=======================================================!

    use constants, only: DP
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    integer, intent(in) :: num
    complex(kind=DP), intent(inout) :: cmatrix(num,num) ! inverse on exit

    ! cks: internal declarations
    integer :: work_length, info, row, col
    integer, allocatable, dimension(:) :: ipiv
    complex(kind=DP), allocatable, dimension(:) :: work_array
    integer :: ierr ! error flag

    ! LAPACK subroutines
    external :: zsytrf, zsytri

    !call timer_clock('linalg_invert_sym_cmatrix',1)

    work_length=3*num

    allocate(ipiv(num),stat=ierr)
    call utils_alloc_check('linalg_invert_sym_cmatrix','ipiv',ierr)
    allocate(work_array(work_length),stat=ierr)
    call utils_alloc_check('linalg_invert_sym_cmatrix','work_array',ierr)

    ! cks: compute the factorization of a complex symmetric matrix A using
    !      the Bunch-Kaufman diagonal pivoting method
    call zsytrf('L', num, cmatrix, num, ipiv, work_array, work_length, info)

    call utils_assert(info == 0, 'Error in linalg_invert_sym_cmatrix(): &
         &ZSYTRF returned info = ', info)

    ! cks: compute the inverse of a complex symmetric indefinite matrix A using
    ! cks: the factorization A = U*D*U**T or A = L*D*L**T computed by CSYTRF
    call zsytri('L', num, cmatrix, num, ipiv, work_array, info )

    call utils_assert(info == 0, 'Error in linalg_invert_sym_cmatrix(): &
         &ZSYTRI returned info = ', info)

    ! cks: fill in the upper triangular part of the complex symmetric matrix
    do row=1,num
       do col=1,row-1
          cmatrix(col,row)=cmatrix(row,col)
       enddo
    enddo

    deallocate(work_array,stat=ierr)
    call utils_dealloc_check('linalg_invert_sym_cmatrix','work_array',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('linalg_invert_sym_cmatrix','ipiv',ierr)

    !call timer_clock('linalg_invert_sym_cmatrix',2)

  end subroutine linalg_invert_sym_cmatrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg_invert_matrix(matrix,info)

    !=======================================================!
    ! This subroutine returns the inverse of a full square  !
    ! real matrix.                                          !
    ! The matrix is overwritten on exit.                    !
    !-------------------------------------------------------!
    ! Written by Robert Bell, May 2014 based on             !
    ! linalg_invert_sym_matrix by Chris-Kriton Skylaris     !
    !=======================================================!

    use constants, only: DP
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    real(kind=DP), intent(inout) :: matrix(:,:) ! inverse on exit
    integer, intent(out) :: info

    ! cks: internal declarations
    integer :: num
    integer :: work_length
    integer, allocatable, dimension(:) :: ipiv
    real(kind=DP), allocatable, dimension(:) :: work_array
    integer :: ierr ! error flag

    ! LAPACK subroutines
    external :: dgetrf, dgetri

    num = size(matrix,dim=1)
    call utils_assert(num == size(matrix,dim=2),'Error in linalg_invert_&
         &matrix: matrix must be square',num,size(matrix,dim=2))

    work_length=3*num

    allocate(ipiv(num),stat=ierr)
    call utils_alloc_check('linalg_invert_matrix','ipiv',ierr)
    allocate(work_array(work_length),stat=ierr)
    call utils_alloc_check('linalg_invert_matrix','work_array',ierr)

    ! rab207: compute the LU factorization of a complex symmetric matrix A
    call dgetrf(num, num, matrix, num, ipiv, info)
    if (info /= 0) goto 100 ! rab207: return info rather than erroring

    ! rab207: compute the inverse of a complex symmetric indefinite matrix A
    call dgetri( num, matrix, num, ipiv, work_array, work_length, info )

100 continue ! ensure deallocation is done

    deallocate(work_array,stat=ierr)
    call utils_dealloc_check('linalg_invert_matrix','work_array',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('linalg_invert_matrix','ipiv',ierr)

  end subroutine linalg_invert_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_invert_cmatrix(cmatrix,info)

    !=======================================================!
    ! This subroutine returns the inverse of a full square  !
    ! complex matrix.                                       !
    ! The matrix is overwritten on exit.                    !
    !-------------------------------------------------------!
    ! Written by Robert Bell, May 2014 based on             !
    ! linalg_invert_sym_matrix by Chris-Kriton Skylaris   !
    !=======================================================!

    use constants, only: DP
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    complex(kind=DP), intent(inout) :: cmatrix(:,:) ! inverse on exit
    integer, intent(out) :: info

    ! cks: internal declarations
    integer :: num
    integer :: work_length
    integer, allocatable, dimension(:) :: ipiv
    complex(kind=DP), allocatable, dimension(:) :: work_array
    integer :: ierr ! error flag

    ! LAPACK subroutines
    external :: zgetrf, zgetri

    num = size(cmatrix,dim=1)
    call utils_assert(num == size(cmatrix,dim=2),'Error in linalg_invert_&
         &cmatrix: cmatrix must be square',num,size(cmatrix,dim=2))

    work_length=3*num

    allocate(ipiv(num),stat=ierr)
    call utils_alloc_check('linalg_invert_cmatrix','ipiv',ierr)
    allocate(work_array(work_length),stat=ierr)
    call utils_alloc_check('linalg_invert_cmatrix','work_array',ierr)

    ! rab207: compute the LU factorization of a complex symmetric matrix A
    call zgetrf(num, num, cmatrix, num, ipiv, info)
    if (info /= 0) goto 100 ! rab207: return info rather than erroring

    ! rab207: compute the inverse of a complex symmetric indefinite matrix A
    call zgetri( num, cmatrix, num, ipiv, work_array, work_length, info )

100 continue ! ensure deallocation is done

    deallocate(work_array,stat=ierr)
    call utils_dealloc_check('linalg_invert_cmatrix','work_array',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('linalg_invert_cmatrix','ipiv',ierr)

  end subroutine linalg_invert_cmatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dcopy(nn,xx,incx,yy,incy)

    !================================================!
    ! dcopy wrapper.                                 !
    !------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/7/2001. !
    !================================================!

    ! cks: BLAS vector copy. yy<--xx
    use constants, only: DP
    implicit none

    integer :: nn, incx, incy
    real(kind=DP), intent(in) :: xx(nn)
    real(kind=DP), intent(out) :: yy(nn)

    ! BLAS subroutine
    external :: dcopy

    call dcopy(nn,xx,incx,yy,incy)

  end subroutine linalg_dcopy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dscal(nn, aa, xx, incx )

    !================================================!
    ! dscal wrapper.                                 !
    !------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 27/7/2001. !
    !================================================!

    ! cks: wrapper for BLAS scale of vector  by scalar
    use constants, only: DP
    implicit none

    integer, intent(in) :: nn, incx
    real(kind=DP), intent(in) :: aa
    real(kind=DP), intent(out) :: xx(nn)

    ! BLAS subroutine
    external :: dscal

    call dscal(nn, aa, xx, incx)

  end subroutine linalg_dscal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_daxpy(nn,aa,xx,incx,yy,incy)

    !================================================!
    ! daxpy wrapper.                                 !
    !------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 22/8/2001. !
    !================================================!

    use constants, only: DP
    implicit none

    integer, intent(in) :: nn, incx, incy
    real(kind=DP), intent(in) :: aa
    real(kind=DP), intent(in) :: xx(nn)
    real(kind=DP), intent(inout) :: yy(nn)

    ! BLAS subroutine
    external :: daxpy

    call daxpy(nn,aa,xx,incx,yy,incy)

  end subroutine linalg_daxpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_1d_fft(dir, isize, in_array, out_array)

    !================================================!
    ! Wrapper for 1-dimensional FFT                  !
    !------------------------------------------------!
    ! Written by David O'Regan in March 2009         !
    ! based on fftbench by A. Mostofi                !
    ! Modified by Jolyon Aarons in May 2017 to use   !
    ! temporary buffers because the FFTW3 version    !
    ! was giving strange results (due to the intent  !
    ! out nature of the input vector in the FFTW3    !
    ! plan creation routine).
    !================================================!

    use constants, only: DP, stdout, LONG, cmplx_0
    use rundat, only: pub_debug_on_root
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

#ifdef FFTW
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = 1
  integer, parameter :: FFTW_REAL_TO_COMPLEX = -1
  integer, parameter :: FFTW_COMPLEX_TO_REAL = 1
  integer, parameter :: FFTW_ESTIMATE = 0
  integer, parameter :: FFTW_MEASURE = 1
  integer, parameter :: FFTW_OUT_OF_PLACE = 0
  integer, parameter :: FFTW_IN_PLACE = 8
  integer, parameter :: FFTW_USE_WISDOM = 16
  integer(kind=LONG) :: the_plan
#endif

#ifdef FFTW3
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = 1
  !integer, parameter :: FFTW_MEASURE = 0
  integer, parameter :: FFTW_ESTIMATE = 64
  integer(kind=LONG) :: the_plan
#endif

#ifdef ACML
  integer, parameter :: ACML_MODE_PLAN = 100
  integer, parameter :: ACML_MODE_FORWARD = 1
  integer, parameter :: ACML_MODE_BACKWARD = -1
  complex(kind=DP), allocatable   :: ACML_COMM(:)
  integer            :: ACML_INFO
#endif

  character, intent(in) :: dir
  integer, intent(in) :: isize
  complex(kind=DP), intent(in) :: in_array(isize)
  complex(kind=DP), intent(inout) :: out_array(isize)
  integer :: zwork_len
  complex(kind=DP), allocatable :: zwork(:)
  complex(kind=DP), allocatable :: intmp(:), outtmp(:)
  integer :: ierr

    call timer_clock('linalg_1d_fft', 1)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering linalg_1d_fft'

    allocate(intmp(isize),stat=ierr)
    call utils_alloc_check('linalg_1d_fft','intmp',ierr)
    intmp=cmplx_0
    allocate(outtmp(isize),stat=ierr)
    call utils_alloc_check('linalg_1d_fft','outtmp',ierr)
    outtmp=cmplx_0


!    ! Copy original into out_array
!    out_array = in_array
    ! Copy original into temporary array
    intmp=in_array


    ! Initialise FFT routine
    call fourier_1d_init(dir,isize)

    ! Do forwards or backwards FFT
    call fourier_1d_apply(dir,isize,intmp, outtmp)

    out_array=outtmp

    ! Finalise FFT routine
    call fourier_1d_exit

    ! Workspace deallocation

   if(allocated(zwork)) then
       deallocate(zwork,stat=ierr)
       call utils_dealloc_check('fourier_1d_fft','zwork',ierr)
    end if
#ifdef ACML
    deallocate(ACML_COMM,stat=ierr)
    call utils_dealloc_check('fourier_1d_fft','ACML_COMM',ierr)
#endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving linalg_1d_fft'

    call timer_clock('linalg_1d_fft', 2)


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

contains

  subroutine fourier_1d_init(dir,n)

    use utils, only: utils_alloc_check, utils_assert

    implicit none

    character, intent(in) :: dir
    integer, intent(in) :: n

    ! Local variables

    integer :: ierr

#ifdef FFTW
    external :: fftw_f77_create_plan
#endif
#ifdef FFTW3
    external :: dfftw_plan_dft_1d
#endif
#ifdef ACML
    external :: zfft1d
#endif

    zwork_len = n

       allocate(zwork(zwork_len),stat=ierr)
       call utils_alloc_check('fourier_1d_fft','zwork',ierr)

#ifdef ACML
       allocate(ACML_COMM(3*zwork_len+100),stat=ierr)
       call utils_alloc_check('fourier_1d_fft','ACML_COMM',ierr)
#endif
    ! Platform-dependent initialisation

#ifdef FFTW
    if (dir == 'B' .or. dir == 'b') then
    call fftw_f77_create_plan(the_plan,n, &
         FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
    else
    call fftw_f77_create_plan(the_plan,n, &
         FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
    endif
#endif

#ifdef FFTW3
    if (dir == 'B' .or. dir == 'b') then
    call dfftw_plan_dft_1d(the_plan,n,intmp,outtmp,FFTW_BACKWARD,FFTW_ESTIMATE)
    else
    call dfftw_plan_dft_1d(the_plan,n,intmp,outtmp,FFTW_FORWARD,FFTW_ESTIMATE)
    endif
#endif

#ifdef ACML
    call zfft1d(ACML_MODE_PLAN,n,zwork,ACML_COMM,ACML_INFO)
    call utils_assert(ACML_INFO == 0, 'Error in linalg_1d_fft &
          &(linalg_mod.F90): zfft1d failed with code ',ACML_INFO)
#endif

  end subroutine fourier_1d_init

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_1d_apply(dir,n,inarray,outarray)

    use constants, only: DP
    use utils, only: utils_assert

    implicit none

    character, intent(in) :: dir
    integer, intent(in) :: n
    complex(kind=DP), intent(inout) :: inarray(n)
    complex(kind=DP), intent(inout) :: outarray(n)

    real(kind=DP) :: scale
#ifdef ACML
    ! jd: Quick fix: ACML_MODE_PLAN was not writable.
    integer :: local_acml_mode_plan
#endif

    ! Platform-dependent call

#ifdef FFTW
    external :: fftw_f77_one
    outarray=inarray
    call fftw_f77_one(the_plan,outarray,zwork)
    if (dir == 'B' .or. dir == 'b') then
       scale = 1.0_dp / n
       outarray = outarray * scale
    end if
#endif

#ifdef FFTW3
    external :: dfftw_execute_dft

    call dfftw_execute_dft(the_plan,inarray,outarray)
    if (dir == 'B' .or. dir == 'b') then
       scale = 1.0_DP / n
       outarray = outarray * scale
    end if
#endif

#ifdef ACML
    external :: zfft1d
    outarray=inarray
    if (dir == 'B' .or. dir == 'b') then
       local_acml_mode_plan = 1
       scale = 1.0_DP / n
       call zfft1d(LOCAL_ACML_MODE_PLAN,n,outarray,ACML_COMM,ACML_INFO)
       call utils_assert(ACML_INFO == 0, 'Error in linalg_1d_fft &
            &(linalg_mod.F90): zfft1d failed with code ',ACML_INFO)
    else
       local_acml_mode_plan = -1
       call zfft1d(LOCAL_ACML_MODE_PLAN,n,outarray,ACML_COMM,ACML_INFO)
       call utils_assert(ACML_INFO == 0, 'Error in linalg_1d_fft &
            &(linalg_mod.F90): zfft1d failed with code ',ACML_INFO)
    end if
#endif

  end subroutine fourier_1d_apply

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_1d_exit

    use utils, only: utils_dealloc_check

    implicit none

    ! Platform-dependent finalisation

#ifdef FFTW
    external :: fftw_f77_destroy_plan

    call fftw_f77_destroy_plan(the_plan)
#endif

#ifdef FFTW3
    external :: dfftw_destroy_plan

    call dfftw_destroy_plan(the_plan)
#endif


  end subroutine fourier_1d_exit

  end subroutine linalg_1d_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dsyev_lt(eigenvectors,eigenvalues,num)
  ! lpl: Solves normal symmetric eigenvalue equation Av=av

    use constants, only: DP, stdout
    use rundat, only: pub_debug_on_root
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    implicit none

    integer, intent(in) :: num
    real(kind=DP), intent(inout) :: eigenvectors(num,num)
    real(kind=DP), intent(out) :: eigenvalues(num)

    real(kind=DP), dimension(:), allocatable :: work
    integer :: info
    integer :: ierr ! error flag

    ! LAPACK subroutine
    external :: dsyev

    call timer_clock('linalg_dsyev_lt', 1)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering linalg_dsyev_lt'

    allocate(work(5*num), stat=ierr)
    call utils_alloc_check('linalg_dsyev_lt','work',ierr)

    call dsyev('V','U',num,eigenvectors,num,eigenvalues,work,5*num,info)

    call utils_assert(info == 0, 'Error in linalg_dsyev_lt: &
         &DSYEV returned info = ', info)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('linalg_dsyev_lt','work',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving linalg_dsyev_lt'

    call timer_clock('linalg_dsyev_lt', 2)

  end subroutine linalg_dsyev_lt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine linalg_zheev_lt(eigenvectors,eigenvalues,num)
  ! rab207: Solves complex hermitian eigenvalue equation Av=av
  ! based on code by Louis Lee
  ! Rewritten by JM Escartin, 07/10/2016.

    use constants, only: DP, stdout
    use rundat, only: pub_debug_on_root
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_safe_nint
    implicit none

    integer, intent(in) :: num
    complex(kind=DP), intent(inout) :: eigenvectors(num,num)
    real(kind=DP), intent(out) :: eigenvalues(num)

    real(kind=DP), dimension(:), allocatable :: rwork
    complex(kind=DP), dimension(:), allocatable :: work
    complex(kind=DP) :: zlwork(1)
    integer :: lwork
    integer :: info
    integer :: ierr ! error flag

    ! LAPACK subroutine
    external :: zheev

    call timer_clock('linalg_zheev_lt', 1)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering linalg_zheev_lt'

    !!! jme: fixed?
    !!! @FIXME rab207: segfaults in zgeev with size(rwork) = lwork
    !!! Possibly a bug in MKL?

    ! Allocate second workspace
    allocate(rwork(max(1,3*num-2)), stat=ierr)
    call utils_alloc_check('linalg_zheev_lt','rwork',ierr)

    ! Find optimal size of first workspace.
    call zheev('V', 'U', num, eigenvectors, num, eigenvalues, zlwork, -1, &
         rwork, info)

    call utils_assert(info == 0, 'Error in linalg_zheev_lt: &
         &when querying for optimal LWORK, ZHEEV returned info = ', info)
    lwork = utils_safe_nint(real(zlwork(1)))

    ! Allocate first workspace
    allocate(work(lwork), stat=ierr)
    call utils_alloc_check('linalg_zheev_lt','work',ierr)

    ! LAPACK call
    call zheev('V', 'U', num, eigenvectors, num, eigenvalues, work, lwork, &
         rwork, info)

    call utils_assert(info == 0, 'Error in linalg_zheev_lt: &
         &ZHEEV returned info = ', info)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('linalg_zheev_lt','work',ierr)
    deallocate(rwork,stat=ierr)
    call utils_dealloc_check('linalg_zheev_lt','rwork',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving linalg_zheev_lt'

    call timer_clock('linalg_zheev_lt', 2)

  end subroutine linalg_zheev_lt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg_zgeev(eigenvectors,eigenvalues,num)
  ! rab207: Solves complex general matrix eigenvalue equation Av=av
  ! based on code by Louis Lee
  ! computes only the RIGHT eigenvectors
  ! Rewritten by JM Escartin, 07/10/2016.

    use constants, only: DP, stdout
    use rundat, only: pub_debug_on_root
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_safe_nint
    implicit none

    integer, intent(in) :: num
    complex(kind=DP), intent(inout) :: eigenvectors(num,num)
    complex(kind=DP), intent(out)   :: eigenvalues(num)

    real(kind=DP), dimension(:), allocatable :: rwork
    complex(kind=DP), allocatable :: work(:)
    complex(kind=DP), allocatable :: A(:,:)
    complex(kind=DP) :: zlwork(1)
    integer :: lwork
    integer :: info
    integer :: ierr ! error flag
    complex(kind=DP) :: vl_dummy(1)  ! left eigenvectors not referenced

    ! LAPACK subroutine
    external :: zgeev

    call timer_clock('linalg_zgeev', 1)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering linalg_zgeev'

    !!! jme: fixed???
    !!! @FIXME rab207: segfaults in zgeev with size(rwork) = lwork
    !!! Possibly a bug in MKL?

    ! Matrix to be diagonalised
    allocate(A(num,num),stat=ierr)
    call utils_alloc_check('linalg_zgeev','A',ierr)
    A(:,:) = eigenvectors(:,:)

    ! Allocate second workspace
    allocate(rwork(2*num), stat=ierr)
    call utils_alloc_check('linalg_zgeev','rwork',ierr)

    ! Find optimal size of first workspace.
    call zgeev('N','V', num, A, num, eigenvalues, vl_dummy, 1, eigenvectors, &
         num, zlwork, -1, rwork, info)

    call utils_assert(info == 0, 'Error in linalg_zgeev: &
         &when querying for optimal LWORK, ZGEEV returned info = ', info)
    lwork = utils_safe_nint(real(zlwork(1)))

    ! Allocate first workspace
    allocate(work(lwork), stat=ierr)
    call utils_alloc_check('linalg_zgeev','work',ierr)

    ! LAPACK call
    call zgeev('N','V',num, A, num, eigenvalues, vl_dummy, 1, eigenvectors, &
         num, work, lwork, rwork, info)

    call utils_assert(info == 0, 'Error in linalg_zgeev: &
         &ZGEEV returned info = ', info)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('linalg_zgeev','work',ierr)
    deallocate(rwork,stat=ierr)
    call utils_dealloc_check('linalg_zgeev','rwork',ierr)
    deallocate(A,stat=ierr)
    call utils_dealloc_check('linalg_zgeev','A',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving linalg_zgeev'

    call timer_clock('linalg_zgeev', 2)

  end subroutine linalg_zgeev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine linalg_zhegv_lt(eigenvectors,eigenvalues,metric,num,itype)
  ! rab207: Solves complex Hermitian generalised eigenvalue equation Av=aBv
  ! based on code by Louis Lee
  ! Rewritten by JM Escartin, 07/10/2016.

    use constants, only: DP, stdout
    use rundat, only: pub_debug_on_root
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_safe_nint
    implicit none

    integer, intent(in) :: num
    complex(kind=DP), intent(inout) :: eigenvectors(num,num)
    complex(kind=DP), intent(inout) :: metric(num,num)
    real(kind=DP), intent(out) :: eigenvalues(num)
    integer, intent(in), optional :: itype

    complex(kind=DP), dimension(:), allocatable :: work
    real(kind=DP), dimension(:), allocatable :: rwork
    complex(kind=DP) :: zlwork(1)
    integer :: info, itype_loc, lwork
    integer :: ierr ! error flag

    ! LAPACK subroutine
    external :: zhegv

    call timer_clock('linalg_zhegv_lt', 1)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering linalg_zhegv_lt'

    itype_loc = 1
    if (present(itype)) itype_loc = itype

    !!!! jme: fixed?
    !!!! @FIXME rab207: segfaults in zgeev with size(rwork) = lwork
    !!!! Possibly a bug in MKL?
    ! Allocate second workspace.
    allocate(rwork(2*num),stat=ierr)
    call utils_alloc_check('linalg_zhegv_lt','rwork',ierr)

    ! Find optimal size of first workspace.
    call zhegv(itype_loc, 'V', 'U', num, eigenvectors, num, metric, num, &
         eigenvalues, zlwork, -1, rwork, info)

    call utils_assert(info == 0, 'Error in linalg_zhegv_lt: &
         &when querying for optimal LWORK, ZHEGV returned info = ', info)
    lwork = utils_safe_nint(real(zlwork(1)))

    ! Allocate first workspace.
    allocate(work(lwork), stat=ierr)
    call utils_alloc_check('linalg_zhegv_lt','work',ierr)

    ! LAPACK call
    call zhegv(itype_loc, 'V', 'U', num, eigenvectors, num, metric, num, &
         eigenvalues, work, lwork, rwork, info)

    call utils_assert(info == 0, 'Error in linalg_zhegv_lt: &
         &ZHEGV returned info = ', info)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('linalg_zhegv_lt','work',ierr)
    deallocate(rwork,stat=ierr)
    call utils_dealloc_check('linalg_zhegv_lt','rwork',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving linalg_zhegv_lt'

    call timer_clock('linalg_zhegv_lt', 2)

  end subroutine linalg_zhegv_lt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_mat_mul_mm_params(shapeC,shapeA,shapeB,opA,opB, & ! input
                                        LDA,LDB,LDC,M,N,K)              ! output

    !=========================================================================!
    ! Determines the parameters for BLAS *GEMM calls                          !
    ! Private subroutine                                                      !
    !-------------------------------------------------------------------------!
    ! Written by Robert Bell, 4 April 2014                                    !
    !=========================================================================!

     use utils, only: utils_assert

     implicit none

     ! arguments
     integer, intent(in)          :: shapeC(2), shapeA(2), shapeB(2)
     character(len=1), intent(in) :: opA, opB
     integer, intent(out)         :: LDA, LDB, LDC, M, N, K

     ! internal
     integer :: shape_opA(2), shape_opB(2)

     LDA=shapeA(1)
     LDB=shapeB(1)
     LDC=shapeC(1)

     ! swap axes if op = 'T' or 'C'
     if (opA == 'T' .or. opA == 'C') then
        shape_opA(:) = (/shapeA(2),shapeA(1)/)
     else
        shape_opA(:) = shapeA(:)
     endif
     if (opB == 'T' .or. opB == 'C') then
        shape_opB(:) = (/shapeB(2),shapeB(1)/)
     else
        shape_opB(:) = shapeB(:)
     endif

     M = shape_opA(1)
     K = shape_opA(2)
     N = shape_opB(2)

     ! check parameters
     call utils_assert(shape_opB(1)==K, 'Error in linalg_mat_mul_mm_params: &
         &rows op(A) /= cols op(B)',shape_opA(2),shape_opB(1))
     call utils_assert(shapeC(1)==M, 'Error in linalg_mat_mul_mm_params: &
         &rows op(A) /= rows C',shape_opA(1),shapeC(1))
     call utils_assert(shapeC(2)==N, 'Error in linalg_mat_mul_mm_params: &
         &cols op(B) /= cols C',shape_opB(2),shapeC(2))

  end subroutine linalg_mat_mul_mm_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_mat_mul_mm_dd(C,A,B,alpha,beta,opA,opB)

    !=========================================================================!
    ! Compute the matrix/matrix multiply:                                     !
    !                                                                         !
    !  C := alpha*A*B + beta*C                                                !
    !                                                                         !
    ! Calls standard BLAS code DGEMM                                          !
    ! Version for real A, real B matrices                                     !
    !-------------------------------------------------------------------------!
    ! Written by Robert Bell, 4 April 2014                                    !
    !=========================================================================!

    use constants, only: DP

    implicit none

    ! arguments
    real(kind=DP), intent(inout) :: C(:,:)
    real(kind=DP), intent(in)    :: A(:,:)
    real(kind=DP), intent(in)    :: B(:,:)
    real(kind=DP), intent(in), optional    :: alpha, beta
    character(len=1), intent(in), optional :: opA, opB

    ! internal
    integer          :: M, N, K, LDA, LDB, LDC
    character(len=1) :: opA_loc, opB_loc
    real(kind=DP)    :: alpha_loc, beta_loc
    integer          :: shape_A(2), shape_B(2), shape_C(2)

    ! BLAS subroutine
    external :: dgemm

    opA_loc = 'N'
    opB_loc = 'N'
    if (present(opA)) opA_loc = opA
    if (present(opB)) opB_loc = opB

    alpha_loc = 1.0_DP
    if (present(alpha)) alpha_loc = alpha
    beta_loc = 0.0_DP
    if (present(beta)) beta_loc = beta

    shape_A = shape(A)
    shape_B = shape(B)
    shape_C = shape(C)

    call linalg_mat_mul_mm_params(shape_C,shape_A,shape_B,opA_loc,opB_loc, &
                                 LDA,LDB,LDC,M,N,K)
    call dgemm(opA_loc,opB_loc,M,N,K,alpha_loc,A,LDA,B,LDB,beta_loc,C,LDC)

  end subroutine linalg_mat_mul_mm_dd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_mat_mul_mm_zz(C,A,B,alpha,beta,opA,opB)

    !=========================================================================!
    ! Compute the matrix/matrix multiply:                                     !
    !                                                                         !
    !  C := alpha*A*B + beta*C                                                !
    !                                                                         !
    ! Calls standard BLAS code ZGEMM                                          !
    ! Version for complex A, complex B matrices                               !
    !-------------------------------------------------------------------------!
    ! Written by Robert Bell, 4 April 2014                                    !
    !=========================================================================!

    use constants, only: DP, cmplx_0, cmplx_1

    implicit none

    ! arguments
    complex(kind=DP), intent(inout) :: C(:,:)
    complex(kind=DP), intent(in)    :: A(:,:)
    complex(kind=DP), intent(in)    :: B(:,:)
    complex(kind=DP), intent(in), optional :: alpha, beta
    character(len=1), intent(in), optional :: opA, opB

    ! internal
    integer          :: M, N, K, LDA, LDB, LDC
    character(len=1) :: opA_loc, opB_loc
    complex(kind=DP) :: alpha_loc, beta_loc
    integer          :: shape_A(2), shape_B(2), shape_C(2)

    ! BLAS subroutine
    external :: zgemm

    opA_loc = 'N'
    opB_loc = 'N'
    if (present(opA)) opA_loc = opA
    if (present(opB)) opB_loc = opB

    alpha_loc = cmplx_1
    if (present(alpha)) alpha_loc = alpha
    beta_loc = cmplx_0
    if (present(beta)) beta_loc = beta

    shape_A = shape(A)
    shape_B = shape(B)
    shape_C = shape(C)

    call linalg_mat_mul_mm_params(shape_C,shape_A,shape_B,opA_loc,opB_loc, &
                                 LDA,LDB,LDC,M,N,K)
    call zgemm(opA_loc,opB_loc,M,N,K,alpha_loc,A,LDA,B,LDB,beta_loc,C,LDC)

  end subroutine linalg_mat_mul_mm_zz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dgesv(xx, & ! output
       aa, bb, num, nrhs)         ! input

    !=========================================================================!
    ! This subroutine solves a linear system AX=B where A is a general matrix.!
    ! B and X are N by NRHS matrices. NRHS is the number of right hand sides. !
    ! Uses LAPACK.                                                            !
    !-------------------------------------------------------------------------!
    ! Written by Quintin Hill on 21/12/2007  for dposv                        !
    ! Modified by Quintin Hill on 30/01/2007 for dgesv                        !
    !=========================================================================!

    use constants, only: DP
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    implicit none

    integer, intent(in) :: num
    integer, intent(in) :: nrhs
    real(kind=DP), intent(out):: xx(num,nrhs) !X
    real(kind=DP), intent(in) :: aa(num, num) !A
    real(kind=DP), intent(in) :: bb(num,nrhs) !B
    integer, allocatable :: pivots(:)
    real(kind=DP), allocatable :: aainternal(:,:)
    integer :: ierr
    integer :: info

    ! LAPACK subroutine
    external :: dgesv

#ifdef NO_OMP_IN_DGESV
!$  integer :: omp_nthreads_save
!$  logical :: omp_dynamic_save
!$  integer, external :: omp_get_max_threads
!$  logical, external :: omp_get_dynamic
#endif

    allocate(aainternal(num,num),stat=ierr)
    call utils_alloc_check('linalg_dgesv','aainternal',ierr)
    allocate(pivots(num), stat=ierr)
    call utils_alloc_check('linalg_dgesv','pivots',ierr)

    aainternal = aa
    xx = bb

    call timer_clock('linalg_dgesv',1)

#ifdef NO_OMP_IN_DGESV
! On Archer2 libsci, only with gfortran, dgesv() perfroms poorly when
! threads are used. Reported to EPCC and HPE as Query Q1767358.
! The workaround is to disable threading within libsci for this call,
! and to restore it later.
!$  omp_nthreads_save = omp_get_max_threads()
!$  omp_dynamic_save = omp_get_dynamic()
!$  call omp_set_dynamic(.false.)
!$  call omp_set_num_threads(1)
#endif

    call dgesv(num, nrhs, aainternal, num, pivots, xx, num, info)

    call utils_assert(info == 0, 'Error in linalg_dgesv(): &
         &DGESV returned info = ', info)

#ifdef NO_OMP_IN_DGESV
!$  call omp_set_num_threads(omp_nthreads_save)
!$  call omp_set_dynamic(omp_dynamic_save)
#endif

    call timer_clock('linalg_dgesv',2)

    deallocate(aainternal, stat=ierr)
    call utils_dealloc_check('linalg_dgesv','aainternal',ierr)
    deallocate(pivots, stat=ierr)
    call utils_dealloc_check('linalg_dgesv','pivots',ierr)


  end subroutine linalg_dgesv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linalg_dgelss(xx, & ! output
       aa, bb, num, nrhs)         ! input

    !=========================================================================!
    ! This subroutine solves a linear system AX=B where A is a general matrix.!
    ! B and X are N by NRHS matrices. NRHS is the number of right hand sides. !
    ! Uses LAPACK.                                                            !
    !-------------------------------------------------------------------------!
    ! Written by Quintin Hill on 27/01/2009  for dgelss                       !
    !=========================================================================!

    use constants, only: DP
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    implicit none

    integer, intent(in) :: num
    integer, intent(in) :: nrhs
    real(kind=DP), intent(out):: xx(num,nrhs) !X
    real(kind=DP), intent(in) :: aa(num, num) !A
    real(kind=DP), intent(in) :: bb(num,nrhs) !B
    real(kind=DP), allocatable :: ss(:)
    real(kind=DP), allocatable :: aainternal(:,:)
    real(kind=DP), allocatable :: work(:)
    real(kind=DP) :: rcond
    integer :: rank
    integer :: ierr
    integer :: info
    integer :: lwork

    ! LAPACK subroutine
    external :: dgelss

    lwork = 3*num + max(2*num,nrhs)

    allocate(aainternal(num,num),stat=ierr)
    call utils_alloc_check('linalg_dgelss','aainternal',ierr)
    allocate(ss(num), stat=ierr)
    call utils_alloc_check('linalg_dgelss','ss',ierr)
    allocate(work(lwork), stat=ierr)
    call utils_alloc_check('linalg_dgelss','work',ierr)

    rcond = -1.0_DP
    aainternal = aa
    xx = bb

    call timer_clock('linalg_dgelss',1)

    call dgelss(num, num, nrhs, aainternal, num, xx, num, ss, rcond, rank, &
         work, lwork, info)

    call utils_assert(info == 0, 'Error in linalg_dgelss: &
         &DGELSS returned info = ', info)
    call timer_clock('linalg_dgelss',2)

    deallocate(work, stat=ierr)
    call utils_dealloc_check('linalg_dgelss','work',ierr)
    deallocate(ss, stat=ierr)
    call utils_dealloc_check('linalg_dgelss','ss',ierr)
    deallocate(aainternal, stat=ierr)
    call utils_dealloc_check('linalg_dgelss','aainternal',ierr)

  end subroutine linalg_dgelss

end module linalg
