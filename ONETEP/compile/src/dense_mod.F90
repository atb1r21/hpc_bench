! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Dense Matrix Algebra Module
!
!   The subroutines in this file were written by
!
!   Nicholas D.M. Hine
!
!   Subsequent additions by Jolyon Aarons, Jacek Dziedzic, Alvaro
!   Ruiz-Serrano, Rob Bell, David Turban, Jian-Hao Li, Jose M Escartin,
!   and Andrea Greco.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dense

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  ! Type definition for the structure of a dense matrix
  type DEM

    ! Variables required for both LAPACK and ScaLAPACK matrices
    logical :: iscmplx
    integer :: nrows, mcols
    real(kind=DP),allocatable :: dmtx(:,:)
    complex(kind=DP),allocatable :: zmtx(:,:)

#ifdef SCALAPACK
    ! Variables only required for the BLACS distributed matrices used
    ! by ScaLAPACK
    integer :: blacs_desc(50)=0
    integer :: blacs_nb
    integer :: blacs_ld
    integer :: blacs_ncol
#endif

  end type DEM

  ! Public Subroutines and types
  public :: DEM
  public :: dense_init
  public :: dense_exit
  public :: dense_create
  public :: dense_destroy
  public :: dense_product
  public :: dense_vector_product
  public :: dense_eigensolve
  public :: dense_svdsolve
  public :: dense_invert
  public :: dense_convert
  public :: dense_axpy
  public :: dense_rank1_update
  public :: dense_scale
  public :: dense_norm
  public :: dense_norm_vec
  public :: dense_dot_product
  public :: dense_trace
  public :: dense_transpose
  public :: dense_write
  public :: dense_copy
  public :: dense_get_element
  public :: dense_put_element
  public :: dense_get_col
  public :: dense_put_col
  public :: dense_determinant
  public :: dense_mat_sqrt
  public :: dense_vec_axpy
  ! agrecocmplx
  public :: dense_safe_axpy
  public :: dense_max_imag_part
  public :: dense_take_real_part

  public :: dense_redist

  ! lpl: Currently only used in npa_mod
  public :: dense_normal_eigensolve
  public :: dense_point2point_element

  ! Create routines
  interface dense_create
     module procedure dense_create_new
     module procedure dense_create_copy
  end interface

  ! AXPY routines
  interface dense_axpy
     module procedure dense_axpy_real
     module procedure dense_axpy_complex
  end interface

  ! agrecocmplx
  ! safe AXPY routines
  interface dense_safe_axpy
     module procedure dense_safe_axpy_real
     module procedure dense_safe_axpy_complex
  end interface

  ! Rank one update routines
  interface dense_rank1_update
     module procedure dense_rank1_update_real
     module procedure dense_rank1_update_complex
  end interface

  ! Scale and shift routines
  interface dense_scale
     module procedure dense_scale_real
     module procedure dense_scale_complex
     module procedure dense_scale_complex2
  end interface

  ! Conversion routines
  interface dense_convert
     module procedure dense_convert_densetosparse
     module procedure dense_convert_sparsetodense
     module procedure dense_convert_densetosparseembed
     module procedure dense_convert_sparseembedtodense
     module procedure dense_convert_sparseembedtosparse
     module procedure dense_convert_sparsetosparseembed
  end interface

  ! Element operation routines
  interface dense_get_element
     module procedure dense_get_element_real
     module procedure dense_get_element_complex
     module procedure dense_get_element_coef
  end interface
  interface dense_put_element
     module procedure dense_put_element_real
     module procedure dense_put_element_complex
     module procedure dense_put_element_coef
  end interface
  interface dense_get_col
     module procedure dense_get_col_real
     module procedure dense_get_col_complex
  end interface
  interface dense_put_col
     module procedure dense_put_col_real
     module procedure dense_put_col_complex
  end interface

  ! lpl: dense_point2point
  interface dense_point2point_element
     module procedure dense_point2point_element_real
     module procedure dense_point2point_element_complex
  end interface

  ! Determinant routines
  interface dense_determinant
     module procedure dense_determinant_real
     module procedure dense_determinant_complex
  end interface

#ifdef SCALAPACK
  ! Module variables describing BLACS context
  integer :: blacs_nprow, blacs_npcol
  integer :: blacs_myrow, blacs_mycol
  integer :: blacs_iam, blacs_nprocs
  integer :: blacs_context, blacs_root_context, blacs_global_context
  integer,dimension(:),allocatable :: im_contexts
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine initialises the variables in the dense matrix algebra      !
  ! module.                                                                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_init

#ifdef SCALAPACK
    use comms, only: pub_total_num_procs,pub_root_proc_id,pub_image_comm,&
                pub_world_comm, comms_bcast,pub_on_root,pub_imroots_comm,&
                comms_barrier
    use rundat, only: pub_num_images
#endif

    use constants, only: stdout
    use image_comms, only: pub_my_image, pub_roots_in_world, pub_image_size,&
                pub_my_rank_in_imroots
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check,&
                utils_flush

    implicit none

#ifdef SCALAPACK

    integer :: i, j, k, ierr, procnum
    integer :: t_row, t_col, t_context
    integer, dimension(1,1) :: map
    integer, allocatable, dimension(:,:) :: usermap

    ! BLACS subroutines
    external :: blacs_setup, blacs_get, blacs_pinfo, &
         blacs_gridinit, blacs_gridinfo, blacs_gridmap

    ! Choose a suitable value for nprow
    blacs_npcol = pub_total_num_procs
    do blacs_nprow=int(sqrt(real(pub_total_num_procs,kind=DP))),1,-1
       if (int(pub_total_num_procs/blacs_nprow)*blacs_nprow== &
            pub_total_num_procs) then
          blacs_npcol = pub_total_num_procs/blacs_nprow
          exit
       end if
    end do

    ! Setup BLACS grid
    call blacs_pinfo(blacs_iam, blacs_nprocs)
    if (blacs_nprocs < 1) then
       call blacs_setup(blacs_iam, blacs_nprow*blacs_npcol)
    end if

    ! Get the blacs context corresponding to the system default, i.e.
    ! mpi_comm_world
    call blacs_get(-1, 0, blacs_global_context)

    if(pub_num_images > 1)then
        allocate(im_contexts(pub_num_images))
        ! If we have more than one ONETEP image we need to map them onto BLACS
        ! grids and also receive the new BLACS context that corresponds to them
        if(pub_image_size /= blacs_nprow*blacs_npcol) call utils_abort(&
         "Error in dense_init: BLACS context size does not match image size",&
         pub_my_image, pub_image_size, blacs_nprow, blacs_npcol)

        do k=1,pub_num_images
           ! kkbd: TODO Redo this so we use image info to calculate these
           ! quantities without broadcasts...
           t_row = blacs_nprow
           t_col = blacs_npcol

           if(pub_on_root)then
              call comms_barrier(pub_imroots_comm)
              call comms_bcast(k-1,t_row,length=1,comm=pub_imroots_comm)
              call comms_barrier(pub_imroots_comm)
              call comms_bcast(k-1,t_col,length=1,comm=pub_imroots_comm)
           end if
           call comms_barrier(pub_image_comm)
           call comms_bcast(0,t_row,length=1,comm=pub_image_comm)
           call comms_bcast(0,t_col,length=1,comm=pub_image_comm)

           allocate(usermap(t_row, t_col),stat=ierr)
           call utils_alloc_check('dense_init','usermap',ierr)

           procnum = pub_roots_in_world(k)
           do j=1,t_col
              do i=1,t_row
                 usermap(i,j) = procnum
                 procnum = procnum + 1
              end do
           end do

           t_context = blacs_global_context

           call blacs_gridmap(t_context, usermap, blacs_nprow, &
             blacs_nprow, blacs_npcol)

           im_contexts(k) = t_context

           if(pub_my_image == k-1)then
              blacs_context = t_context
           end if

           deallocate(usermap,stat=ierr)
           call utils_dealloc_check('dense_init','usermap',ierr)
        end do
    else
        ! kkbd: Do normal dense_init stuff.
        !       Don't allocate im_contexts becuase EDA runs come through here
        !       more than once... Really they shouldn't do that, but it's tied
        !       to energy_and_force_init_cell. A potentially better solution is
        !       to exit dense_init at the start if a dense_init_done flag is
        !       .true., where it's set to true in the body of the routine and
        !       set to .false. on dense_exit. Then uncomment the other lines
        !       in this block.
        !allocate(im_contexts(1))
        blacs_context = blacs_global_context
        call blacs_gridinit(blacs_context, 'R', blacs_nprow, blacs_npcol)
        !im_contexts(1) = blacs_context
    end if

    call blacs_gridinfo(blacs_context, blacs_nprow, blacs_npcol, &
     blacs_myrow, blacs_mycol)

    call blacs_get(-1, 0, blacs_root_context)
    map(1,1)=pub_root_proc_id
    call blacs_gridmap(blacs_root_context, map, 1, 1, 1)

#endif

  end subroutine dense_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine cleans up variables associated with the dense matrix       !
  ! algebra module.                                                            !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_exit

    use utils, only: utils_dealloc_check
    use rundat, only: pub_num_images

    implicit none

#ifdef SCALAPACK
    integer :: nprow, npcol, myrow, mycol, ierr

    ! BLACS subroutines
    external :: blacs_gridinfo, blacs_gridexit

    ! Remove BLACS contexts
    ! kkbd TODO exit image contexts
    if(blacs_myrow /= -1)then
      call blacs_gridexit(blacs_context)
    end if

    call blacs_gridinfo(blacs_root_context, nprow, npcol, myrow, mycol)
    if(nprow/=-1) then
       call blacs_gridexit(blacs_root_context)
    end if

    if (pub_num_images > 1) then
       deallocate(im_contexts,stat=ierr)
       call utils_dealloc_check("dense_exit","im_contexts",ierr)
    end if

#endif

  end subroutine dense_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine allocates storage for a dense matrix, and, if applicable,  !
  ! sets up a BLACS descriptor for it.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The dense matrix to create                            !
  !   nrows    (input) : The number of rows in the matrix                      !
  !   mcols    (input) : The number of columns in the matrix                   !
  !   iscmplx  (input) : Whether the matrix is to be complex                   !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_create_new(mat,nrows,mcols,iscmplx,local_to_root)

    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: mat
    integer, intent(in) :: nrows
    integer, intent(in) :: mcols
    logical, intent(in), optional :: iscmplx
    logical, intent(in), optional :: local_to_root

    ! Local Variables
    logical :: loc_iscmplx ! Local copy of optional iscmplx argument
    logical :: loc_locroot
    integer :: ierr        ! Error flag
#ifdef SCALAPACK
    external :: blacs_gridinfo, descset, descinit ! BLACS & ScaLAPACK subroutines
    integer, external :: numroc
    integer :: nprow, npcol, myrow, mycol
#endif

    ! Set default for optional argument
    loc_iscmplx = .false.
    loc_locroot=.false.
    if (present(iscmplx)) loc_iscmplx = iscmplx
    if(present(local_to_root)) loc_locroot=local_to_root
    mat%iscmplx = loc_iscmplx
    mat%nrows = nrows
    mat%mcols = mcols

    ! Check arguments
    call utils_assert(nrows>=0, 'Internal error in dense_create_new: nrows < 0')
    call utils_assert(mcols>=0, 'Internal error in dense_create_new: mcols < 0')

#ifdef SCALAPACK

    ! Choose a suitable blocking factor


    ! Initialise ld and ncol

    if(loc_locroot) then
       mat%blacs_nb = 1
    else
       !       mat%blacs_nb = min(mat%nrows/blacs_nprow,mat%mcols/blacs_npcol,1)
       mat%blacs_nb = max(min(mat%nrows/blacs_nprow,mat%mcols/blacs_npcol),1)
    end if

    ! Initialise ld and ncol
    if(loc_locroot) then
       mat%blacs_ld = numroc(mat%nrows,mat%blacs_nb,blacs_myrow,0,1)
       mat%blacs_ncol = numroc(mat%mcols,mat%blacs_nb,blacs_mycol,0,1)
    else
       mat%blacs_ld = max(numroc(mat%nrows,mat%blacs_nb,blacs_myrow,0,blacs_nprow),1)
       !    mat%blacs_ld = numroc(mat%nrows,mat%blacs_nb,blacs_myrow,0,blacs_nprow)
       mat%blacs_ncol = numroc(mat%mcols,mat%blacs_nb,blacs_mycol,0,blacs_npcol)
    end if

    ! Create matrix description
    if(loc_locroot) then
       call blacs_gridinfo(blacs_root_context, nprow, npcol, myrow, mycol)

       if(nprow/=-1) then
          call descinit(mat%blacs_desc, mat%nrows,mat%mcols, mat%blacs_nb, mat%blacs_nb, &
               & 0, 0, blacs_root_context, mat%blacs_ld, ierr)
       else
          call descset(mat%blacs_desc, mat%nrows,mat%mcols,mat%blacs_nb, mat%blacs_nb, &
               & 0, 0, -1, size(mat%dmtx))
          ierr=0
       end if
    else

       call descinit(mat%blacs_desc, mat%nrows, mat%mcols, mat%blacs_nb, mat%blacs_nb, &
            0, 0, blacs_context, mat%blacs_ld, ierr)
    end if

    ! Check return value ierr
    call utils_assert(ierr == 0, 'Error in dense_create_new: &
            &BLACS descinit call returned ierr=',ierr)

    ! Allocate arrays
    if (mat%iscmplx) then
       allocate(mat%zmtx(mat%blacs_ld,max(mat%blacs_ncol,1)),stat=ierr)
       call utils_alloc_check('dense_create_new','mat%zmtx',ierr)
    else
       allocate(mat%dmtx(mat%blacs_ld,max(mat%blacs_ncol,1)),stat=ierr)
       call utils_alloc_check('dense_create_new','mat%dmtx',ierr)
    end if

#else

    ! Allocate arrays
    if (mat%iscmplx) then
       allocate(mat%zmtx(mat%nrows,mat%mcols),stat=ierr)
       call utils_alloc_check('dense_create_new','mat%zmtx',ierr)
    else
       allocate(mat%dmtx(mat%nrows,mat%mcols),stat=ierr)
       call utils_alloc_check('dense_create_new','mat%dmtx',ierr)
    end if

#endif

    ! Set contents to zero
    if (mat%iscmplx) then
       mat%zmtx = cmplx(0.0_DP,0.0_DP,kind=DP)
    else
       mat%dmtx = 0.0_DP
    end if

  end subroutine dense_create_new


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine allocates storage for a dense matrix, and, if applicable,  !
  ! sets up a BLACS descriptor for it.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The dense matrix to create                            !
  !   nrows    (input) : The number of rows in the matrix                      !
  !   mcols    (input) : The number of columns in the matrix                   !
  !   iscmplx  (input) : Whether the matrix is to be complex                   !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_create_copy(mat,origmat,iscmplx)

    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: mat
    type(DEM), intent(in   ) :: origmat
    logical, intent(in), optional :: iscmplx

    ! Local Variables
    logical :: loc_iscmplx ! Local copy of optional iscmplx argument
    integer :: ierr        ! Error flag
#ifdef SCALAPACK
    integer, external :: numroc
    external :: descinit ! ScaLAPACK subroutine
#endif

    ! Set default for optional argument
    loc_iscmplx = origmat%iscmplx
    if (present(iscmplx)) loc_iscmplx = iscmplx
    mat%iscmplx = loc_iscmplx
    mat%nrows = origmat%nrows
    mat%mcols = origmat%mcols

    ! Check arguments
    call utils_assert(mat%nrows>=0, &
         'Internal error in dense_create_new: nrows < 0')
    call utils_assert(mat%mcols>=0, &
         'Internal error in dense_create_new: mcols < 0')

#ifdef SCALAPACK

    ! Choose a suitable blocking factor
    mat%blacs_nb = max(min(mat%nrows/blacs_nprow,mat%mcols/blacs_npcol),1)

    ! Initialise ld and ncol
    mat%blacs_ld = max(numroc(mat%nrows,mat%blacs_nb,blacs_myrow,0,blacs_nprow),1)
    mat%blacs_ncol = numroc(mat%mcols,mat%blacs_nb,blacs_mycol,0,blacs_npcol)

    ! Create matrix description
    call descinit(mat%blacs_desc, mat%nrows, mat%mcols, mat%blacs_nb, mat%blacs_nb, &
         0, 0, blacs_context, mat%blacs_ld, ierr)
    ! Check return value ierr
    call utils_assert(ierr == 0, 'Error in dense_create_new: &
            &BLACS descinit call returned ierr=',ierr)

    ! Allocate arrays
    if (mat%iscmplx) then
       allocate(mat%zmtx(mat%blacs_ld,max(mat%blacs_ncol,1)),stat=ierr)
       call utils_alloc_check('dense_create_new','mat%zmtx',ierr)
    else
       allocate(mat%dmtx(mat%blacs_ld,max(mat%blacs_ncol,1)),stat=ierr)
       call utils_alloc_check('dense_create_new','mat%dmtx',ierr)
    end if

#else

    ! Allocate arrays
    if (mat%iscmplx) then
       allocate(mat%zmtx(mat%nrows,mat%mcols),stat=ierr)
       call utils_alloc_check('dense_create_new','mat%zmtx',ierr)
    else
       allocate(mat%dmtx(mat%nrows,mat%mcols),stat=ierr)
       call utils_alloc_check('dense_create_new','mat%dmtx',ierr)
    end if

#endif

    ! Set contents to zero
    if (mat%iscmplx) then
       mat%zmtx = cmplx(0.0_DP,0.0_DP,kind=DP)
    else
       mat%dmtx = 0.0_DP
    end if

  end subroutine dense_create_copy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !============================================================================!
  ! This subroutine deallocates storage for a dense matrix.                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat     (input) : The dense matrix to deallocate.                        !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_destroy(mat)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: mat

    ! Local Variables
    integer :: ierr        ! Error flag

#ifdef SCALAPACK
    ! Deallocate arrays
    if (mat%iscmplx) then
       deallocate(mat%zmtx,stat=ierr)
       call utils_dealloc_check('dense_destroy','mat%zmtx',ierr)
    else
       deallocate(mat%dmtx,stat=ierr)
       call utils_dealloc_check('dense_destroy','mat%dmtx',ierr)
    end if
#else
    ! Deallocate arrays
    if (mat%iscmplx) then
       deallocate(mat%zmtx,stat=ierr)
       call utils_dealloc_check('dense_destroy','mat%zmtx',ierr)
    else
       deallocate(mat%dmtx,stat=ierr)
       call utils_dealloc_check('dense_destroy','mat%dmtx',ierr)
    end if
#endif

  end subroutine dense_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine calculates the product of two dense matrices A and B and   !
  ! stores it in the dense matrix C.                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   cmat      (inout) : The dense matrix C                                   !
  !   amat      (in)    : The dense matrix A                                   !
  !   bmat      (in)    : The dense matrix B                                   !
  !   opA, opB  (in,optional) : operations to apply to A/B = 'N'/'T'/'C'       !
  !                             (N:none ; T:transpose ; C:Hermitian conjugate) !
  !   first_k, last_k (in,optional) : Range of k-values to sum over in product !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23/02/2010, based partially on wrappers_dgemm,   !
  ! which was written by Chris-Kriton Skylaris on 2/3/2004.                    !
  ! Modified by Robert Bell to take opA instead of transpose_amat to allow     !
  ! both Hermitian conjugate and transpose operations. 28/04/2014              !
  !============================================================================!

  subroutine dense_product(cmat,amat,bmat,opA,opB,first_k,last_k)

    use comms, only: comms_bcast, pub_my_proc_id, pub_total_num_procs
    use constants, only: cmplx_0, cmplx_1
    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: cmat
    type(DEM), intent(in) :: amat
    type(DEM), intent(in) :: bmat
    character(len=1), optional, intent(in) :: opA
    character(len=1), optional, intent(in) :: opB
    integer, optional, intent(in) :: first_k, last_k

    ! Local Variables
    integer :: m,n,k ! sizes of the matrices (see DGEMM manpages)
#ifndef SCALAPACK
    integer :: batch_size, local_start, local_end, local_n
    integer :: ldb
    integer :: proc
#endif
    integer :: ia,ja,ib,jb
    character :: opA_loc,opB_loc

    ! (P)BLAS subroutines
#ifdef SCALAPACK
    external :: PZGEMM, PDGEMM
#else
    external :: ZGEMM, DGEMM
#endif

    ! Deal with optional argument for transposing matrices
    opA_loc = 'n'
    if (present(opA)) then
       opA_loc=opA
    end if
    opB_loc = 'n'
    if (present(opB)) then
       opB_loc = opB
    end if

    ! Check matrices all real or all complex
    if ((bmat%iscmplx.neqv.amat%iscmplx) .or. &
         (cmat%iscmplx.neqv.amat%iscmplx)) then
       call utils_abort('Error in dense_product: &
            &matrices must be all real or all complex')
    end if

    ! Check Matrix Sizes
    m = cmat%nrows
    n = cmat%mcols

    if (opA_loc=='n'.or.opA_loc=='N') then
       ! ab: number of rows of A should be equal to number of rows of C
       call utils_assert(amat%nrows == m, 'Error in dense_product: &
               &amat%nrows /= cmat%nrows')
       k = amat%mcols
    else
       ! ab: number of cols of A should be equal to number of rows of C
       call utils_assert(amat%mcols == m, 'Error in dense_product: &
               &amat%mcols /= cmat%nrows')
       k = amat%nrows
    end if
    if (opB_loc=='n'.or.opB_loc=='N') then
       call utils_assert(bmat%mcols == n, 'Error in dense_product: &
               &bmat%mcols /= cmat%mcols')
       call utils_assert(bmat%nrows == k, 'Error in dense_product: &
               &bmat%nrows /= amat%mcols')
    else ! opB_loc=='t' .or. opB_loc=='c'
       call utils_assert(bmat%nrows == n, 'Error in dense_product: &
               &bmat%nrows /= cmat%mcols')
       call utils_assert(bmat%mcols == k, 'Error in dense_product: &
               &bmat%mcols /= amat%mcols')
    end if

    ! Deal with optional arguments for only multiplying certain row ranges
    if (present(first_k).and.present(last_k)) then
       ! Check range is valid
       if ((first_k < 0).or.(last_k > k).or.(last_k < first_k)) then
          call utils_abort('Error in dense_product: &
               &invalid k range specification')
       end if
       ! Start from first_k'th column
       ia = 1
       ja = first_k
       if (opB_loc=='n'.or.opB_loc=='N') then
          ib = first_k
          jb = 1
       else
          ib = 1
          jb = first_k
       end if
       ! Set range of k
       k = last_k - first_k + 1
    else
       ! Use whole range
       ia = 1
       ja = 1
       ib = 1
       jb = 1
    end if

    ! Perform matrix product with  BLAS    [LAPACK]  (d/z)gemm
    !                          or PBLAS [ScaLAPACK] p(d/z)gemm
#ifdef SCALAPACK

    if (amat%iscmplx) then
       call PZGEMM(opA_loc, opB_loc, m, n, k, cmplx_1, &
            amat%zmtx(1,1), ia, ja, amat%blacs_desc, &
            bmat%zmtx(1,1), ib, jb, bmat%blacs_desc, cmplx_0, &
            cmat%zmtx(1,1), 1, 1, cmat%blacs_desc)
    else
       call PDGEMM(opA_loc, opB_loc, m, n, k, 1.0_DP, &
            amat%dmtx(1,1), ia, ja, amat%blacs_desc, &
            bmat%dmtx(1,1), ib, jb, bmat%blacs_desc, 0.0_DP, &
            cmat%dmtx(1,1), 1, 1, cmat%blacs_desc)
    end if

#else

    ! Semi-parallelised DGEMM call - procs work out their columns of result
    batch_size = n/pub_total_num_procs
    local_start = pub_my_proc_id*batch_size + 1
    local_end = local_start + batch_size - 1
    if (pub_my_proc_id == pub_total_num_procs-1) local_end = n
    local_n = local_end - local_start + 1

    ! Sort out what local cols/rows of B to select, based on opB
    if (opB_loc=='n'.or.opB_loc=='N') then
       jb = local_start
       ldb = bmat%nrows
    else
       ib = local_start
       ldb = n
    end if

    ! Perform the LAPACK call
    if (amat%iscmplx) then
       call ZGEMM(opA_loc, opB_loc, m, local_n, k, cmplx_1, amat%zmtx(ia,ja), m, &
            bmat%zmtx(ib,jb), ldb, cmplx_0, &
            cmat%zmtx(1,local_start), m)
    else
       call DGEMM(opA_loc, opB_loc, m, local_n, k, 1.0_DP, amat%dmtx(ia,ja), m, &
            bmat%dmtx(ib,jb), ldb, 0.0_DP, &
            cmat%dmtx(1,local_start), m)
    end if

    ! Broadcast results from each proc to fill whole matrix of C
    do proc=0,pub_total_num_procs-1
       local_start = proc*batch_size + 1
       local_end = local_start + batch_size - 1
       if (proc == pub_total_num_procs-1) local_end = n
       local_n = local_end - local_start + 1

       if (amat%iscmplx) then
          call comms_bcast(proc,cmat%zmtx(1,local_start),m*local_n)
       else
          call comms_bcast(proc,cmat%dmtx(1,local_start),m*local_n)
       end if
    end do

#endif

  end subroutine dense_product

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine calculates the product of two dense matrices A and B and   !
  ! stores it in the dense matrix C.                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   cmat      (inout) : The dense matrix C                                   !
  !   amat      (in)    : The dense matrix A                                   !
  !   bmat      (in)    : The dense matrix B                                   !
  !   transpose_bmat (in,optional)  : Whether to transpose the B matrix        !
  !   first_k, last_k (in,optional) : Range of k-values to sum over in product !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23/02/2010, based partially on wrappers_dgemm,   !
  ! which was written by Chris-Kriton Skylaris on 2/3/2004.                    !
  !============================================================================!

  subroutine dense_vector_product(cvec,amat,bvec)

    use comms, only: comms_bcast
    use constants, only: cmplx_0, cmplx_1
    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: cvec
    type(DEM), intent(in) :: amat
    type(DEM), intent(in) :: bvec

    ! Local Variables
    integer :: m,n ! sizes of the matrices (see DGEMV manpages)
    integer :: ia,ja,ib,jb

    ! (P)BLAS subroutines
#ifdef SCALAPACK
    external :: pzgemv, pdgemv
#else
    external :: zgemv, dgemv
#endif

    ! Check matrices all real or all complex
    if ((bvec%iscmplx.neqv.amat%iscmplx) .or. &
         (cvec%iscmplx.neqv.amat%iscmplx)) then
       call utils_abort('Error in dense_vector_product: &
            &matrices and vectors must be all real or all complex')
    end if

    ! Check Matrix Sizes
    m = amat%nrows
    n = amat%mcols
!
!    call utils_assert(amat%nrows == m, 'Error in dense_product: &
!            &amat%nrows /= cmat%nrows')
!    if (transb=='n') then
!       call utils_assert(bmat%mcols == n, 'Error in dense_product: &
!               &bmat%mcols /= cmat%mcols')
!       call utils_assert(bmat%nrows == k, 'Error in dense_product: &
!               &bmat%nrows /= amat%mcols')
!    else ! transb=='t'
!       call utils_assert(bmat%nrows == n, 'Error in dense_product: &
!               &bmat%nrows /= cmat%mcols')
!       call utils_assert(bmat%mcols == k, 'Error in dense_product: &
!               &bmat%mcols /= amat%mcols')
!    end if
!

    ! Use whole range
    ia = 1
    ja = 1
    ib = 1
    jb = 1

    ! Perform matrix product with LAPACK (d/z)gemm or ScaLAPACK p(d/z)gemm
#ifdef SCALAPACK
    if (amat%iscmplx) then
       call pzgemv('N',m,n,cmplx_1,amat%zmtx,1,1,amat%blacs_desc,bvec%zmtx,1,1,bvec%blacs_desc,1, &
            & cmplx_0,cvec%zmtx,1,1,cvec%blacs_desc,1)
    else
       call pdgemv('N',m,n,1.0_dp,amat%dmtx,1,1,amat%blacs_desc,bvec%dmtx,1,1,bvec%blacs_desc,1, &
            & 0.0_dp,cvec%dmtx,1,1,cvec%blacs_desc,1)
    end if
#else
    if (amat%iscmplx) then
       call zgemv('N',m,n,cmplx_1,amat%zmtx(1,1),size(amat%zmtx,1),bvec%zmtx(1,1),1,cmplx_0,cvec%zmtx(1,1),1)
    else
       call dgemv('N',m,n,1.0_dp,amat%dmtx(1,1),size(amat%dmtx,1),bvec%dmtx(1,1),1,0.0_dp,cvec%dmtx(1,1),1)
    end if
#endif
  end subroutine dense_vector_product


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine finds the solution to the generalised eigenproblem defined !
  ! by two dense matrices amat and bmat, and returns the real eigenvalues and  !
  ! a dense matrix storing the eigenvectors.                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   eigenvals (out)   : Real array storing the eigenvalues of amat.          !
  !   amat (input)      : Dense Symmetric or Hermitian matrix.                 !
  !                     : WARNING: values changed on output!                   !
  !   bmat (input)      : Dense Symmetric or Hermitian matrix.                 !
  !                     : WARNING: values changed on output!                   !
  !   ibtype (input)    : Which type of eigenproblem to solve (see LAPACK help)!
  !   eigenvecs (out)   : Dense matrix storing eigenvectors defining solution. !
  !   arg_orfac (in)    :                                                      !
  !   arg_abstol (in)   :                                                      !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  ! Fixed to broadcast resulting eigenvectors from root proc by Nicholas Hine  !
  ! on 10 March 2010, to fix problem of different eigenvectors on different    !
  ! procs.                                                                     !
  ! Added warning that a and b matrices are overwritten, 25 September 2013     !
  !============================================================================!

  subroutine dense_eigensolve(n_eig,eigenvals,amat,bmat,ibtype,eigenvecs,arg_orfac,arg_abstol)

    use comms, only: comms_bcast, pub_on_root, pub_root_proc_id
    use constants, only: stdout
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
#ifdef SCALAPACK
#ifdef ELPA2
    use comms, only: pub_world_comm
    use constants, only: cmplx_1
    use elpa2, only: solve_evp_complex_2stage, solve_evp_real_2stage
    use elpa1, only: get_elpa_row_col_comms, invert_trm_real, cholesky_real, &
         & invert_trm_complex, cholesky_complex
#endif
#endif

    implicit none

    ! Arguments
    integer, intent(in) :: n_eig
    real(kind=DP), intent(out) :: eigenvals(:)
    type(DEM), intent(inout) :: eigenvecs
    type(DEM), intent(inout) :: amat
    type(DEM), intent(inout) :: bmat
    integer, intent(in) :: ibtype
    real(kind=DP), optional, intent(in) :: arg_orfac ! ars
    real(kind=DP), optional, intent(in) :: arg_abstol

    ! Local Variables
#ifdef SCALAPACK
    integer :: nzfound
    logical :: orthogonalisation_memory
    real(kind=DP), external :: pdlamch
#ifdef NO_OMP_IN_SCALAPACK_DIAG
!$  integer :: omp_nthreads_save
!$  logical :: omp_dynamic_save
!$  integer, external :: omp_get_max_threads
!$  logical, external :: omp_get_dynamic
#endif
#endif
    real(kind=DP) :: abstol, orfac
    integer :: nfound
    integer :: ierr, info
    character :: jobz, range

    ! Work arrays and sizes
    integer :: lwork
#ifdef SCALAPACK
    integer :: liwork, lrwork
#endif
    integer, allocatable :: iwork(:), ifail(:)
    real(kind=DP), allocatable :: work(:)
    complex(kind=DP), allocatable :: zwork(:)
#ifdef SCALAPACK
    integer, allocatable :: iclustr(:)
    real(kind=DP), allocatable :: gap(:)
#endif
#ifdef ELPA2
    integer :: mpi_comm_rows, mpi_comm_cols
    type(DEM) :: tmpmat
    type(DEM) :: tmpmat2
    real(kind=dp), dimension(:,:), allocatable :: tmparr
#endif

    ! (Sca)LAPACK subroutines
#ifdef SCALAPACK
#ifdef ELPA2
#else
    external :: PZHEGVX, PDSYGVX
#endif
#else
    external :: ZHEGVX, DSYGVX
#endif

    ! Start timer
    call timer_clock('dense_eigensolve',1)

    ! Check arguments
    if (n_eig<=0) then
       jobz = 'N'
       range = 'A'
    else
       jobz = 'V'
       if ((n_eig>0).and.(n_eig<amat%nrows)) then
          range = 'I'
       else
          range = 'A'
       end if
    end if

#ifdef SCALAPACK
#ifdef ELPA2
    ! jd: utils/checkout does not support ifndef, so using a dummy ifdef/else construct
#else
    ! ndmh: call ScaLAPACK routine to solve the generalised eigenproblem
    ! ars: set orfac
    if (present(arg_orfac)) then
       if (arg_orfac.le.tiny(1.0_DP)) then
          orfac = 0.0_DP
          orthogonalisation_memory = .false.
       else
          orfac = arg_orfac
          orthogonalisation_memory = .true.
       end if
    else
       orfac = 1e-4_DP
       orthogonalisation_memory = .true.
    end if
#endif
#endif

    ! Allocate work arrays
    call internal_allocate_work

#ifdef SCALAPACK

#if defined(NO_OMP_IN_SCALAPACK_DIAG) && !defined(ELPA2)
! On Archer2 libsci, both with gfortran and Cray, returns extremely noisy
! results when ScaLAPACK is used together with OpenMP. This renders EDFT
! useless on Archer2 without a workaround. For more details see
! https://bitbucket.org/onetep/onetep/issues/1889/on-archer2-scalapack-with-omp-gives-wrong
! The workaround is to disable threading within ScaLAPACK for this call,
! and to restore it later.
! We cannot do this with an #ifdef inside an ELPA2 #ifdef or #else because this
! confuses the 'checkout' script.

!$  omp_nthreads_save = omp_get_max_threads()
!$  omp_dynamic_save = omp_get_dynamic()
!$  call omp_set_dynamic(.false.)
!$  call omp_set_num_threads(1)

#endif

#ifdef ELPA2
    ! Get elpa communicators
    call get_elpa_row_col_comms(pub_world_comm, blacs_myrow, blacs_mycol, &
         & mpi_comm_rows, mpi_comm_cols)
    ! Make temporary arrays
    call dense_create_copy(tmpmat,bmat)
    call dense_create_copy(tmpmat2,amat)
    call dense_copy(tmpmat,bmat)
    call dense_copy(tmpmat2,amat)
    if (amat%iscmplx) then
       ! Firstly convert to standard eigenproblem
       ! Cholesky decomposition of B
       call pzpotrf('U', tmpmat%nrows, tmpmat%zmtx(1,1), 1, 1, tmpmat%blacs_desc, info)
!       call cholesky_complex(tmpmat%nrows, tmpmat%zmtx(1,1), tmpmat%nrows, tmpmat%blacs_nb, mpi_comm_rows, mpi_comm_cols)
       ! Invert triangular U
       call pztrtri('U', 'N', tmpmat%nrows, tmpmat%zmtx(1,1), 1, 1,  tmpmat%blacs_desc, info)
!       call invert_trm_complex(tmpmat%nrows, tmpmat%zmtx(1,1), tmpmat%nrows, tmpmat%blacs_nb, mpi_comm_rows, mpi_comm_cols)
       ! tmpmat2 -> U^{-T} A U^{-1}
!       call dense_product(tmpmat2,tmpmat,amat,.true.,.false.)
!       call dense_product(tmpmat3,tmpmat2,tmpmat,.false.,.false.)
       call pztrmm('L','U','C','N',tmpmat%nrows,tmpmat%nrows,cmplx_1,tmpmat%zmtx(1,1),1,1,tmpmat%blacs_desc,&
            & tmpmat2%zmtx(1,1), 1, 1, tmpmat2%blacs_desc)
       call pztrmm('R','U','N','N',tmpmat%nrows,tmpmat%nrows,cmplx_1,tmpmat%zmtx(1,1),1,1,tmpmat%blacs_desc,&
            & tmpmat2%zmtx(1,1), 1, 1, tmpmat2%blacs_desc)

       ! Solve eigenproblem
!       call solve_evp_complex_2stage(tmpmat3%nrows, n_eig, tmpmat3%zmtx, tmpmat%nrows, eigenvals, &
!            &eigenvecs%zmtx(1,1), tmpmat3%nrows, tmpmat3%blacs_nb, mpi_comm_rows, mpi_comm_cols, pub_world_comm)
       call solve_evp_complex_2stage(tmpmat2%nrows, n_eig, tmpmat2%zmtx, tmpmat2%blacs_ld, eigenvals, &
            &eigenvecs%zmtx(1,1), eigenvecs%blacs_ld, tmpmat2%blacs_nb, mpi_comm_rows, mpi_comm_cols, pub_world_comm)
       ! Recover eigenvectors
       call pztrmm('L','U','N','N',tmpmat%nrows,tmpmat%nrows,cmplx_1,tmpmat%zmtx(1,1),1,1,tmpmat%blacs_desc,&
            & eigenvecs%zmtx(1,1), 1, 1, eigenvecs%blacs_desc)
    else
       ! Firstly convert to standard eigenproblem
       ! Cholesky decomposition of B
       call pdpotrf('U', tmpmat%nrows, tmpmat%dmtx(1,1), 1, 1, tmpmat%blacs_desc, info)
!       call cholesky_real(tmpmat%nrows, tmpmat%dmtx(1,1), tmpmat%nrows, tmpmat%blacs_nb, mpi_comm_rows, mpi_comm_cols)
       ! Invert triangular U
       call pdtrtri('U', 'N', tmpmat%nrows, tmpmat%dmtx(1,1), 1, 1,  tmpmat%blacs_desc, info)
!       call invert_trm_real(tmpmat%nrows, tmpmat%dmtx(1,1), tmpmat%nrows, tmpmat%blacs_nb, mpi_comm_rows, mpi_comm_cols)
       ! tmpmat2 -> U^{-T} A U^{-1}
!       call dense_product(tmpmat2,tmpmat,amat,.true.,.false.)
!       call dense_product(tmpmat3,tmpmat2,tmpmat,.false.,.false.)
       call pdtrmm('L','U','T','N',tmpmat%nrows,tmpmat%nrows,1.0_dp,tmpmat%dmtx(1,1),1,1,tmpmat%blacs_desc,&
            & tmpmat2%dmtx(1,1), 1, 1, tmpmat2%blacs_desc)
       call pdtrmm('R','U','N','N',tmpmat%nrows,tmpmat%nrows,1.0_dp,tmpmat%dmtx(1,1),1,1,tmpmat%blacs_desc,&
            & tmpmat2%dmtx(1,1), 1, 1, tmpmat2%blacs_desc)

       ! Solve eigenproblem
!       call solve_evp_complex_2stage(tmpmat3%nrows, n_eig, tmpmat3%zmtx, tmpmat%nrows, eigenvals, &
!            &eigenvecs%zmtx(1,1), tmpmat3%nrows, tmpmat3%blacs_nb, mpi_comm_rows, mpi_comm_cols, pub_world_comm)
       call solve_evp_real_2stage(tmpmat2%nrows, n_eig, tmpmat2%dmtx, tmpmat2%blacs_ld, eigenvals, &
            &eigenvecs%dmtx(1,1), eigenvecs%blacs_ld, tmpmat2%blacs_nb, mpi_comm_rows, mpi_comm_cols, pub_world_comm)
       ! Recover eigenvectors
       call pdtrmm('L','U','N','N',tmpmat%nrows,tmpmat%nrows,1.0_dp,tmpmat%dmtx(1,1),1,1,tmpmat%blacs_desc,&
            & eigenvecs%dmtx(1,1), 1, 1, eigenvecs%blacs_desc)
    end if
    call dense_destroy(tmpmat)
    call dense_destroy(tmpmat2)
    call mpi_comm_free(mpi_comm_rows,info)
    call mpi_comm_free(mpi_comm_cols,info)
#else
    ! ars: set abstol
    if (present(arg_abstol)) then
       if (arg_abstol.le.tiny(1.0_DP)) then
          abstol = 2.0_DP*PDLAMCH(blacs_context,'S')
       else
          abstol = arg_abstol
       end if
    else
       abstol = 1e-9_DP
    end if

    if (amat%iscmplx) then
       call PZHEGVX(ibtype, jobz, range, 'L', amat%nrows, &
            amat%zmtx(1,1), 1, 1, amat%blacs_desc, &
            bmat%zmtx(1,1), 1, 1, bmat%blacs_desc, &
            0.0_DP, 0.0_DP, 1, n_eig, abstol, nfound, nzfound, &
            eigenvals, orfac, eigenvecs%zmtx(1,1), 1, 1, eigenvecs%blacs_desc, &
            zwork, lwork, work, lrwork, iwork, liwork, ifail, iclustr, gap, &
            info)
    else
       call PDSYGVX(ibtype, jobz, range, 'L', amat%nrows, &
            amat%dmtx(1,1), 1, 1, amat%blacs_desc, &
            bmat%dmtx(1,1), 1, 1, bmat%blacs_desc, &
            0.0_DP, 0.0_DP, 1, n_eig, abstol, nfound, nzfound, &
            eigenvals, orfac, eigenvecs%dmtx(1,1), 1, 1, eigenvecs%blacs_desc, &
            work, lwork, iwork, liwork, ifail, iclustr, gap, info)
    end if


#endif
#else
    if (present(arg_orfac).and.(present(arg_abstol))) then
       orfac = arg_orfac
       abstol = arg_abstol
    end if
    ! ndmh: call LAPACK routine to solve the generalised eigenproblem
    if (amat%iscmplx) then
       call ZHEGVX(ibtype,jobz,range,'L',amat%nrows, &
            amat%zmtx(1,1),amat%nrows,bmat%zmtx(1,1),bmat%nrows,&
            0.0_DP, 0.0_DP, 1, n_eig, -1.0_DP, nfound, eigenvals, &
            eigenvecs%zmtx(1,1),amat%nrows,zwork,lwork,work,iwork,ifail,info)

       ! ndmh: degenerate eigenvectors may have come out in different linear
       ! ndmh: combinations on different procs, so synchronise to eigenvectors
       ! ndmh: found on root proc
       call comms_bcast(pub_root_proc_id,eigenvecs%zmtx, &
            eigenvecs%nrows*eigenvecs%mcols)

       ! ars: broadcast eigenvalues as well
       call comms_bcast(pub_root_proc_id, eigenvals, n_eig)

    else
       call DSYGVX(ibtype,jobz,range,'L',amat%nrows, &
            amat%dmtx(1,1),amat%nrows,bmat%dmtx(1,1),bmat%nrows, &
            0.0_DP, 0.0_DP, 1, n_eig, -1.0_DP, nfound, eigenvals, &
            eigenvecs%dmtx(1,1),amat%nrows,work,lwork,iwork,ifail,info)

       ! ndmh: degenerate eigenvectors may have come out in different linear
       ! ndmh: combinations on different procs, so synchronise to eigenvectors
       ! ndmh: found on root proc
       call comms_bcast(pub_root_proc_id,eigenvecs%dmtx, &
            eigenvecs%nrows*eigenvecs%mcols)

       ! ars: broadcast eigenvalues as well
       call comms_bcast(pub_root_proc_id, eigenvals, n_eig)


    end if

#if defined(NO_OMP_IN_SCALAPACK_DIAG) && !defined(ELPA2)
!$  call omp_set_num_threads(omp_nthreads_save)
!$  call omp_set_dynamic(omp_dynamic_save)
#endif

#endif

    ! ndmh: check for errors
    if (info/=0) then
       if (amat%iscmplx) then
          if (pub_on_root) then
             write(stdout,'(a,i5)') '(P)ZHEGVX in subroutine &
                  &dense_eigensolve returned info=',info
             write(stdout,*) 'ifail=',ifail(1:min(info,size(ifail)))
#ifdef SCALAPACK
             write(stdout,*) 'iclustr=', iclustr(1:min(info,size(iclustr)))
             write(stdout,*) 'gap=', gap(1:min(info,size(gap)))
#endif
          end if
       else
          if (pub_on_root) then
             write(stdout,'(a,i5)') '(P)DSYGVX in subroutine &
                  &dense_eigensolve returned info=',info
             write(stdout,*) 'ifail=',ifail(1:min(info,size(ifail)))
#ifdef SCALAPACK
             write(stdout,*) 'iclustr=', iclustr(1:min(info,size(iclustr)))
             write(stdout,*) 'gap=', gap(1:min(info,size(gap)))
#endif
          end if
       end if
       ! ndmh: continue even if some eigenvectors did not converge
       if ((info<1).or.(info>amat%nrows)) then
          call utils_abort('Error in dense_eigensolve(). &
               &The error message should be in your output file.')
       end if
    end if

    ! Deallocate work arrays
    call internal_deallocate_work

    ! Stop timer
    call timer_clock('dense_eigensolve',2)

    return

contains

    ! Allocate diagonalisation work arrays
    subroutine internal_allocate_work

      implicit none

#ifdef SCALAPACK

      ! Local Variables
      integer :: n, nn, nb, np0, mq0, nq0, anb, sqnpc, nps
      integer :: nsytrd_lwopt, nsygst_lwopt
      integer :: nhetrd_lwopt, nhegst_lwopt
      integer :: liclustr, lgap
      integer, external :: numroc, pjlaenv, iceil

      ! Calculate optimal sizes of work arrays
      n = amat%nrows
      nb = amat%blacs_nb
      nn = max(amat%nrows,nb,2)
      np0 = numroc(nn,nb,0,0,blacs_nprow)
      nq0 = numroc(nn,nb,0,0,blacs_npcol)
      mq0 = numroc(nn,nb,0,0,blacs_npcol)
      if (amat%iscmplx) then
         anb = pjlaenv(blacs_context,3,'PZHETTRD','L',0,0,0,0)
      else
         anb = pjlaenv(blacs_context,3,'PDSYTTRD','L',0,0,0,0)
      end if
      sqnpc = int(sqrt(real(blacs_nprow*blacs_npcol,kind=DP)))
      nps = max(numroc(n,1,0,0,sqnpc),2*anb)
      nsytrd_lwopt = n + 2*(anb+1)*(4*nps+2) + (nps+3)*nps
      nsygst_lwopt = 2*np0*nb + nq0*nb + nb*nb

      if (.not.amat%iscmplx) then ! real version

         ! LWORK >= 5 * N + MAX( 5*NN, NP0 * MQ0 + 2 * NB * NB ) +
         !   ICEIL( NEIG, NPROW*NPCOL)*NN
         lwork = 5*n + max(5*n, np0*nq0 + 2*nb*nb) &
              + iceil(n, blacs_nprow*blacs_npcol)*nn
         if (orthogonalisation_memory) &
              lwork = max(lwork, 5*n + nsytrd_lwopt,nsygst_lwopt)

      else ! complex version

         !  LWORK >= N + ( NP0 + MQ0 + NB ) * NB
         !    NQ0 = NUMROC( NN, NB, 0, 0, NPCOL ).
         !    MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
         !For optimal performance, greater workspace is needed, i.e.
         !  LWORK >= MAX( LWORK, N + NHETRD_LWOPT, NHEGST_LWOPT )
         !Where LWORK is as defined above, and
         !  NHETRD_LWORK = 2*( ANB+1 )*( 4*NPS+2 ) + ( NPS + 1 ) * NPS
         !  NHEGST_LWOPT =  2*NP0*NB + NQ0*NB + NB*NB
         lwork = n + ( np0 + mq0 + nb ) * nb
         nhetrd_lwopt = 2*( anb+1 )*( 4*nps+2 ) + ( nps + 1 ) * nps
         nhegst_lwopt =  2*np0*nb + nq0*nb + nb*nb
         lwork = max(lwork,n + nhetrd_lwopt, nhegst_lwopt)

         lrwork = 4*n + max( 5*nn, np0 * mq0 ) + iceil( nn, blacs_nprow*blacs_npcol)*nn

      end if

      liwork = 6*max(n,blacs_nprow*blacs_npcol+1,4)
      liclustr = 2*blacs_nprow*blacs_npcol
      lgap = blacs_nprow*blacs_npcol

      ! Allocate work arrays
      if (amat%iscmplx) then
         allocate(zwork(lwork),stat=ierr)
         call utils_alloc_check('dense_eigensolve','zwork',ierr)
         allocate(work(lrwork),stat=ierr)
         call utils_alloc_check('dense_eigensolve','work',ierr)
      else
         allocate(work(lwork),stat=ierr)
         call utils_alloc_check('dense_eigensolve','work',ierr)
      end if
      allocate(iwork(liwork),stat=ierr)
      call utils_alloc_check('dense_eigensolve','iwork',ierr)
      allocate(iclustr(liclustr),stat=ierr)
      call utils_alloc_check('dense_eigensolve','iclustr',ierr)
      allocate(gap(liclustr),stat=ierr)
      call utils_alloc_check('dense_eigensolve','gap',ierr)
      allocate(ifail(n),stat=ierr)
      call utils_alloc_check('dense_eigensolve','ifail',ierr)

#else

      ! Allocate work arrays
      lwork = 8*amat%nrows
      if (amat%iscmplx) then
         allocate(zwork(lwork),stat=ierr)
         call utils_alloc_check('dense_eigensolve','zwork',ierr)
      end if
      allocate(work(lwork),stat=ierr)
      call utils_alloc_check('dense_eigensolve','work',ierr)
      allocate(iwork(5*amat%nrows),stat=ierr)
      call utils_alloc_check('dense_eigensolve','iwork',ierr)
      allocate(ifail(amat%nrows),stat=ierr)
      call utils_alloc_check('dense_eigensolve','ifail',ierr)


#endif

    end subroutine internal_allocate_work

    ! Deallocate diagonalisation work arrays
    subroutine internal_deallocate_work

      implicit none

#ifdef SCALAPACK

      deallocate(ifail,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','ifail',ierr)
      deallocate(gap,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','gap',ierr)
      deallocate(iclustr,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','iclustr',ierr)
      deallocate(iwork,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','iwork',ierr)
      deallocate(work,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','work',ierr)
      if (amat%iscmplx) then
         deallocate(zwork,stat=ierr)
         call utils_dealloc_check('dense_eigensolve','zwork',ierr)
      end if

#else

      deallocate(ifail,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','ifail',ierr)
      deallocate(iwork,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','iwork',ierr)
      deallocate(work,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','work',ierr)
      if (amat%iscmplx) then
         deallocate(zwork,stat=ierr)
         call utils_dealloc_check('dense_eigensolve','zwork',ierr)
      end if

#endif

    end subroutine internal_deallocate_work

  end subroutine dense_eigensolve


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine finds the solution to the singular value decomposition     !
  ! (SVD) problem:   A = U * S * Vt  (S is a diagonal matrix)                  !
  ! If the dimension of matrix A is denoted as (M, N), the routine returns the !
  ! U of dimension (M, min[M,N]), V of dimension (min[M,N], N) and the diagonal!
  ! part (singular values) of S is returned as an array of dimension min[M,N]. !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   singulars (out)   : Real array storing singular values.                  !
  !   amat (input)      : Dense matrix to be decomposed.                       !
  !                     : WARNING: values changed on output!                   !
  !   leftvecs (out)    : Dense matrix storing left eigenvectors (column-wise).!
  !   rightvecs (out)   : Dense matrix storing right eigenvectors (row-wise).  !
  !----------------------------------------------------------------------------!
  ! Written by Jian-Hao Li, October 2014.                                      !
  ! Largely rewritten by JM Escartin, 3/10/2016.                               !
  !============================================================================!
  subroutine dense_svdsolve(singulars, amat, leftvecs, rightvecs) !jhl52

    use comms, only: pub_on_root
    use constants, only: stdout
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_safe_nint

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: singulars(:)
    type(DEM), intent(inout) :: amat
    type(DEM), optional, intent(inout), target :: leftvecs
    type(DEM), optional, intent(inout), target :: rightvecs

    ! Local Variables
    integer :: ierr, info
    character :: jobu, jobvt
    type(DEM), target :: U_loc, VT_loc ! local buffers (won't be used)
    type(DEM), pointer :: U, VT        ! pointers to local or optional buffers

    ! Workspaces
    real(kind=DP), allocatable :: work(:), rwork(:)
    complex(kind=DP), allocatable :: zwork(:)

    ! Workspace size
    real(kind=DP) :: lwork(1)
    complex(kind=DP) :: zlwork(1)
    integer :: lwork_int

    ! (Sca)LAPACK subroutines
#ifdef SCALAPACK
    external :: PZGESVD, PDGESVD
#else
    external :: ZGESVD, DGESVD

    ! More local variables
    integer :: lda, ldu, ldvt
#endif


    ! Start timer
    call timer_clock('dense_svdsolve',1)

#ifdef SCALAPACK
    ! Singular vectors required?
    if (present(leftvecs)) then
       jobu  = 'V'
       U => leftvecs
    else
       jobu = 'N'
       call dense_create(U_loc, 1, 1, iscmplx=amat%iscmplx)
       U => U_loc
    end if
    if (present(rightvecs)) then
       jobvt  = 'V'
       VT => rightvecs
    else
       jobvt = 'N'
       call dense_create(VT_loc, 1, 1, iscmplx=amat%iscmplx)
       VT => VT_loc
    end if

    ! Run SVD decomposition
    if (amat%iscmplx) then
       ! Allocate second workspace
       allocate(rwork(1+4*min(amat%nrows, amat%mcols)), stat=ierr)
       call utils_alloc_check('dense_svdsolve', 'rwork', ierr)

       ! Inquiry for size of first workspace
       call PZGESVD(jobu, jobvt, amat%nrows, amat%mcols, amat%zmtx(1,1), 1, 1, &
                    amat%blacs_desc, singulars, U%zmtx(1,1), 1, 1, &
                    U%blacs_desc, VT%zmtx(1,1), 1, 1, &
                    VT%blacs_desc, zlwork, -1, rwork, info)

       call utils_assert(info==0, 'Error in dense_svdsolve: &
            &when querying for optimal LWORK, PZGESVD returned info = ', info)
       lwork_int = utils_safe_nint(real(zlwork(1)))

       ! Allocate first workspace
       allocate(zwork(lwork_int), stat=ierr)
       call utils_alloc_check('dense_svdsolve', 'zwork', ierr)

       ! Singular value decomposition
       call PZGESVD(jobu, jobvt, amat%nrows, amat%mcols, amat%zmtx(1,1), 1, 1, &
                    amat%blacs_desc, singulars, U%zmtx(1,1), 1, 1, &
                    U%blacs_desc, VT%zmtx(1,1), 1, 1, &
                    VT%blacs_desc, zwork, lwork_int, rwork, info)

       call utils_assert(info==0, 'Error in dense_svdsolve &
            &(the output file should contain the SCALAPACK error message): &
            &PZGESVD returned info = ', info)

    else
       ! Inquiry for size of workspace
       call PDGESVD(jobu, jobvt, amat%nrows, amat%mcols, amat%dmtx(1,1), 1, 1, &
                    amat%blacs_desc, singulars, U%dmtx(1,1), 1, 1, &
                    U%blacs_desc, VT%dmtx(1,1), 1, 1, &
                    VT%blacs_desc, lwork, -1, info)

       call utils_assert(info==0, 'Error in dense_svdsolve: &
            &when querying for optimal LWORK, PDGESVD returned info = ', info)
       lwork_int = utils_safe_nint(lwork(1))

       ! Allocate workspace
       allocate(work(lwork_int), stat=ierr)
       call utils_alloc_check('dense_svdsolve', 'work', ierr)

       ! Singular value decomposition
       call PDGESVD(jobu, jobvt, amat%nrows, amat%mcols, amat%dmtx(1,1), 1, 1, &
                    amat%blacs_desc, singulars, U%dmtx(1,1), 1, 1, &
                    U%blacs_desc, VT%dmtx(1,1), 1, 1, &
                    VT%blacs_desc, work, lwork_int, info)

       call utils_assert(info==0, 'Error in dense_svdsolve &
            &(the output file should contain the SCALAPACK error message): &
            &PDGESVD returned info = ', info)

    end if !amat%iscmplx

#else
    ! Singular vectors required?
    if (present(leftvecs)) then
       jobu  = 'S'
       U => leftvecs
    else
       jobu = 'N'
       call dense_create(U_loc, 1, 1, iscmplx=amat%iscmplx)
       U => U_loc
    end if
    if (present(rightvecs)) then
       jobvt  = 'S'
       VT => rightvecs
    else
       jobvt = 'N'
       call dense_create(VT_loc, 1, 1, iscmplx=amat%iscmplx)
       VT => VT_loc
    end if

    ! Set up leading dimensions of A, U, and Vt
    lda = amat%nrows
    ldu = U%nrows
    ldvt = VT%nrows

    ! Run SVD decomposition
    if (amat%iscmplx) then
       ! Allocate second workspace
       allocate(rwork(5*min(amat%nrows, amat%mcols)), stat=ierr)
       call utils_alloc_check('dense_svdsolve', 'rwork', ierr)

       ! Inquiry for size of first workspace
       call ZGESVD(jobu, jobvt, amat%nrows, amat%mcols, amat%zmtx(1,1), &
                    lda, singulars, U%zmtx(1,1), ldu, &
                    VT%zmtx(1,1), ldvt, zlwork, -1, rwork, info)

       call utils_assert(info==0, 'Error in dense_svdsolve: &
            &when querying for optimal LWORK, ZGESVD returned info = ', info)
       lwork_int = utils_safe_nint(real(zlwork(1)))

       ! Allocate first workspace
       allocate(zwork(lwork_int), stat=ierr)
       call utils_alloc_check('dense_svdsolve', 'zwork', ierr)

       ! Singular value decomposition
       call ZGESVD(jobu, jobvt, amat%nrows, amat%mcols, amat%zmtx(1,1), &
                    lda, singulars, U%zmtx(1,1), ldu, &
                    VT%zmtx(1,1), ldvt, zwork, lwork_int, rwork, info)

       call utils_assert(info==0, 'Error in dense_svdsolve &
            &(the output file should contain the LAPACK error message): &
            &ZGESVD returned info = ', info)

    else
       ! Inquiry for size of workspace
       call DGESVD(jobu, jobvt, amat%nrows, amat%mcols, amat%dmtx(1,1), &
                    lda, singulars, U%dmtx(1,1), ldu, &
                    VT%dmtx(1,1), ldvt, lwork, -1, info)

       call utils_assert(info==0, 'Error in dense_svdsolve: &
            &when querying for optimal LWORK, DGESVD returned info = ', info)
       lwork_int = utils_safe_nint(lwork(1))

       ! Allocate workspace
       allocate(work(lwork_int), stat=ierr)
       call utils_alloc_check('dense_svdsolve', 'work', ierr)

       ! Singular value decomposition
       call DGESVD(jobu, jobvt, amat%nrows, amat%mcols, amat%dmtx(1,1), &
                    lda, singulars, U%dmtx(1,1), ldu, &
                    VT%dmtx(1,1), ldvt, work, lwork_int, info)

       call utils_assert(info==0, 'Error in dense_svdsolve &
            &(the output file should contain the LAPACK error message): &
            &DGESVD returned info = ', info)

    end if !amat%iscmplx

    nullify(U)
    nullify(VT)
    if (.not.present(leftvecs)) call dense_destroy(U_loc)
    if (.not.present(rightvecs)) call dense_destroy(VT_loc)
#endif

    ! Deallocate work spaces
    if (amat%iscmplx) then
       deallocate(rwork, stat=ierr)
       call utils_dealloc_check('dense_svdsolve', 'rwork', ierr)
       deallocate(zwork, stat=ierr)
       call utils_dealloc_check('dense_svdsolve', 'work', ierr)
    else
       deallocate(work, stat=ierr)
       call utils_dealloc_check('dense_svdsolve', 'work', ierr)
    end if

    ! Stop timer
    call timer_clock('dense_svdsolve',2)

  end subroutine dense_svdsolve


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine finds the inverse of a general dense matrix.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat (input)      : Dense matrix.                                        !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 1 March 2010.                                    !
  ! April 2014 rab207: Isn't  (P)*GETRF for any matrix (not symmetric)?        !
  !============================================================================!

  subroutine dense_invert(mat,error)

    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_safe_nint

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: mat
    integer, optional, intent(out) :: error

    ! Local Variables
    integer, allocatable :: ipiv(:)
    real(kind=DP), allocatable :: work(:)
    complex(kind=DP), allocatable :: zwork(:)
#ifdef SCALAPACK
    integer, allocatable :: iwork(:)
    integer :: liwork
#endif
    integer :: lwork
    integer :: ierr
    integer :: info

    ! (Sca)LAPACK subroutines
#ifdef SCALAPACK
    external :: PZGETRI, PDGETRI, PZGETRF, PDGETRF
#else
    external :: ZGETRI, DGETRI, ZGETRF, DGETRF
#endif

    if(present(error)) error=0

    ! Check Arguments
    call utils_assert(mat%nrows == mat%mcols,'ERROR in dense_invert: &
            &Matrix is not square')

    ! Allocate work arrays
#ifdef SCALAPACK
    !---------------------------------
    ! rab207: perform workspace query
    lwork = -1
    liwork = -1
    if (mat%iscmplx) then
       allocate(zwork(1),stat=ierr)
       call utils_alloc_check('dense_invert','zwork',ierr)
    else
       allocate(work(1),stat=ierr)
       call utils_alloc_check('dense_invert','work',ierr)
    endif
    allocate(iwork(1),stat=ierr)
    call utils_alloc_check('dense_invert','iwork',ierr)

    allocate(ipiv(1),stat=ierr)
    call utils_alloc_check('dense_invert','ipiv',ierr)

    if (mat%iscmplx) then
       call PZGETRI(mat%nrows, mat%zmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv, zwork, lwork, iwork, liwork, info)
       lwork = utils_safe_nint(real(zwork(1)))
    else
       call PDGETRI(mat%nrows, mat%dmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv, work, lwork, iwork, liwork, info)
       lwork = utils_safe_nint(work(1))
    end if
    liwork = iwork(1)

    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('dense_invert','ipiv',ierr)

    deallocate(iwork,stat=ierr)
    call utils_dealloc_check('dense_invert','iwork',ierr)
    if (mat%iscmplx) then
       deallocate(zwork,stat=ierr)
       call utils_dealloc_check('dense_invert','zwork',ierr)
    else
       deallocate(work,stat=ierr)
       call utils_dealloc_check('dense_invert','work',ierr)
    endif
    ! rab207: end workspace query
    !---------------------------------

    ! rab207: sizes now determined by workspace query, 10/06/14
    allocate(ipiv(mat%nrows),stat=ierr)
    call utils_alloc_check('dense_invert','ipiv',ierr)
    allocate(iwork(liwork),stat=ierr)
    call utils_alloc_check('dense_invert','iwork',ierr)
    if(mat%iscmplx) then
       allocate(zwork(lwork),stat=ierr)
       call utils_alloc_check('dense_invert','zwork',ierr)
    else
       allocate(work(lwork),stat=ierr)
       call utils_alloc_check('dense_invert','work',ierr)
    endif
#else
    lwork = 3*mat%nrows
    allocate(ipiv(mat%nrows),stat=ierr)
    call utils_alloc_check('dense_invert','ipiv',ierr)
    if (mat%iscmplx) then
       allocate(zwork(lwork),stat=ierr)
       call utils_alloc_check('dense_invert','zwork',ierr)
    else
       allocate(work(lwork),stat=ierr)
       call utils_alloc_check('dense_invert','work',ierr)
    endif
#endif

    ! cks: compute the factorization of a real symmetric matrix A using
    !      the Bunch-Kaufman diagonal pivoting method
#ifdef SCALAPACK
    if (mat%iscmplx) then
       call PZGETRF(mat%nrows, mat%nrows, mat%zmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv, info)
    else
       call PDGETRF(mat%nrows, mat%nrows, mat%dmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv, info)
    end if
#else
    if (mat%iscmplx) then
       call ZGETRF(mat%nrows, mat%nrows, mat%zmtx(1,1), mat%nrows, ipiv, info)
    else
       call DGETRF(mat%nrows, mat%nrows, mat%dmtx(1,1), mat%nrows, ipiv, info)
    end if
#endif

    if(present(error)) then
       if(info/=0) then
       error=info
          ! Deallocate work arrays
#ifdef SCALAPACK
          if (mat%iscmplx) then
             deallocate(zwork,stat=ierr)
             call utils_dealloc_check('dense_invert','zwork',ierr)
          else
             deallocate(work,stat=ierr)
             call utils_dealloc_check('dense_invert','work',ierr)
          endif
          deallocate(iwork,stat=ierr)
          call utils_dealloc_check('dense_invert','iwork',ierr)
          deallocate(ipiv,stat=ierr)
          call utils_dealloc_check('dense_invert','ipiv',ierr)
#else
          if (mat%iscmplx) then
             deallocate(zwork,stat=ierr)
             call utils_dealloc_check('dense_invert','zwork',ierr)
          else
             deallocate(work,stat=ierr)
             call utils_dealloc_check('dense_invert','work',ierr)
          endif
          deallocate(ipiv,stat=ierr)
          call utils_dealloc_check('dense_invert','ipiv',ierr)
#endif
       return
       end if
    else
       call utils_assert(info==0, 'ERROR in dense_invert: &
            &(P)DGETRF returned info=',info)
    end if

    ! cks: compute the inverse of a real symmetric indefinite matrix A using
    ! cks: the factorization A = U*D*U**T or A = L*D*L**T computed by DSYTRF
#ifdef SCALAPACK
    if (mat%iscmplx) then
       call PZGETRI(mat%nrows, mat%zmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv,zwork, lwork, iwork, liwork, info)
    else
       call PDGETRI(mat%nrows, mat%dmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv, work, lwork, iwork, liwork, info)
    end if
#else
    if (mat%iscmplx) then
       call ZGETRI(mat%nrows, mat%zmtx(1,1), mat%nrows, ipiv, zwork,lwork, info)
    else
       call DGETRI(mat%nrows, mat%dmtx(1,1), mat%nrows, ipiv, work, lwork, info)
    end if
#endif

    if(present(error)) then
       error=info
!       return
    else
       call utils_assert(info==0, 'ERROR in dense_invert: &
            &(P)DGETRI returned info=',info)
    end if

    ! Deallocate work arrays
#ifdef SCALAPACK
    if (mat%iscmplx) then
       deallocate(zwork,stat=ierr)
       call utils_dealloc_check('dense_invert','zwork',ierr)
    else
       deallocate(work,stat=ierr)
       call utils_dealloc_check('dense_invert','work',ierr)
    endif
    deallocate(iwork,stat=ierr)
    call utils_dealloc_check('dense_invert','iwork',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('dense_invert','ipiv',ierr)
#else
    if (mat%iscmplx) then
       deallocate(zwork,stat=ierr)
       call utils_dealloc_check('dense_invert','zwork',ierr)
    else
       deallocate(work,stat=ierr)
       call utils_dealloc_check('dense_invert','work',ierr)
    endif
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('dense_invert','ipiv',ierr)
#endif

  end subroutine dense_invert


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine computes the determinant of a general real matrix.         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat (input)      : Dense Symmetric or Hermitian matrix.                 !
  !   det  (output)     : Real determinant                                     !
  !----------------------------------------------------------------------------!
  ! Written by Simon Dubois, 16 October 2013.                                  !
  !============================================================================!

  subroutine dense_determinant_real(amat,det)

    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(in)        :: amat
    real(kind=DP), intent(inout) :: det

    ! Local Variables
    type(DEM) :: buffer
    integer, allocatable :: ipiv(:)
    real(kind=DP) :: inc
    integer   :: im
    integer   :: ierr
    integer   :: info

    ! (Sca)LAPACK subroutines
#ifdef SCALAPACK
    external :: PDGETRF
#else
    external :: DGETRF
#endif

    ! Check Arguments
    call utils_assert(amat%nrows == amat%mcols,'ERROR in dense_determinant_real: &
            &Matrix is not square')

    ! Copy original matrix in buffer
    call dense_create(buffer,amat)
    call dense_copy(buffer,amat)

    ! Allocate work arrays
    allocate(ipiv(buffer%nrows),stat=ierr)
    call utils_alloc_check('dense_determinant_real','ipiv',ierr)

    ! cks: compute the factorization of a real symmetric matrix A using
    !      the Bunch-Kaufman diagonal pivoting method
#ifdef SCALAPACK
    call PDGETRF(buffer%nrows, buffer%nrows, buffer%dmtx(1,1), 1, 1, buffer%blacs_desc, &
         ipiv, info)
#else
    call DGETRF(buffer%nrows, buffer%nrows, buffer%dmtx(1,1), buffer%nrows, ipiv, info)
#endif

    call utils_assert(info==0, 'ERROR in dense_determinant_real: &
            &(P)DGETRF returned info=',info)

    ! Compute determinant
    det = 1.0_DP
    do im = 1, buffer%mcols
       call dense_get_element(inc,buffer,im,im)
       if (ipiv(im)/=im) then
          det = -det*inc
       else
          det = det*inc
       endif
    enddo

    ! Deallocate work arrays
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('dense_determinant_real','ipiv',ierr)

    ! Destroy buffer
    call dense_destroy(buffer)

  end subroutine dense_determinant_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine computes the determinant of a general complex matrix.      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat (input)      : Dense Symmetric or Hermitian matrix.                 !
  !   det  (output)     : Complex determinant                                  !
  !----------------------------------------------------------------------------!
  ! Written by Simon Dubois, 16 October 2013.                                  !
  !============================================================================!

  subroutine dense_determinant_complex(amat,det)

    use constants, only: cmplx_1
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(in) :: amat
    complex(kind=DP), intent(inout) :: det

    ! Local Variables
    type(DEM) :: buffer
    integer, allocatable :: ipiv(:)
    complex(kind=DP) :: inc
    integer   :: im
    integer   :: ierr
    integer   :: info

    ! (Sca)LAPACK subroutines
#ifdef SCALAPACK
    external :: PZGETRF
#else
    external :: ZGETRF
#endif

    ! Check Arguments
    call utils_assert(amat%nrows == amat%mcols, &
         'ERROR in dense_determinant_complex: Matrix is not square')

    ! Copy original matrix in buffer
    call dense_create(buffer,amat)
    call dense_copy(buffer,amat)

    ! Allocate work arrays
    allocate(ipiv(buffer%nrows),stat=ierr)
    call utils_alloc_check('dense_determinant_complex','ipiv',ierr)

    ! cks: compute the factorization of a real symmetric matrix A using
    !      the Bunch-Kaufman diagonal pivoting method
#ifdef SCALAPACK
    call PZGETRF(buffer%nrows, buffer%nrows, buffer%zmtx(1,1), 1, 1, buffer%blacs_desc, &
         ipiv, info)
#else
    call ZGETRF(buffer%nrows, buffer%nrows, buffer%zmtx(1,1), buffer%nrows, ipiv, info)
#endif

    call utils_assert(info==0, 'ERROR in dense_determinant_complex: &
            &(P)ZGETRF returned info=',info)

    ! Compute determinant
    det = cmplx_1
    do im = 1, buffer%mcols
       call dense_get_element(inc,buffer,im,im)
       if (ipiv(im)/=im) then
          det = -det*inc
       else
          det = det*inc
       endif
    enddo

    ! Deallocate work arrays
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('dense_determinant_complex','ipiv',ierr)

    ! Destroy buffer
    call dense_destroy(buffer)

  end subroutine dense_determinant_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dense_mat_sqrt(mat,buf1,buf2)
  !============================================================================!
  ! This subroutine computes the square-root of a positive semi-definite matrix!
  ! by rotating into and back out of the diagonal representation.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (inout) : The matrix, will be over-written on return            !
  !   buf1     (inout) : (optional) dense workspace                            !
  !   buf2     (inout) : (optional) dense workspace                            !
  !                                                                            !
  ! Workspace matrices are optional, and if not provided internal workspace    !
  ! will be used.                                                              !
  !----------------------------------------------------------------------------!
  ! Written by Robert Bell, 21/07/2014                                         !
  !============================================================================!

      use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

      implicit none

      ! arguments
      type(DEM), intent(inout) :: mat
      type(DEM), intent(inout), optional, target :: buf1, buf2

      ! internal
      integer :: ii, norb, ierr
      type(DEM), target :: ibuf1, ibuf2 ! internal buffers, if buf1, buf2 not defined
      type(DEM), pointer :: work1, work2 ! pointers to the workspace
      real(kind=DP), allocatable    :: eigvals(:)
      real(kind=DP), allocatable    :: rcol(:)
      complex(kind=DP), allocatable :: zcol(:)
      real(kind=DP), parameter :: tolerance = 0.0001_DP

      norb = mat%nrows
      allocate(eigvals(norb),stat=ierr)
      call utils_alloc_check('dense_mat_sqrt','eigvals',ierr)
      if (mat%iscmplx) then
         allocate(zcol(norb),stat=ierr)
         call utils_alloc_check('dense_mat_sqrt','zcol',ierr)
      else
         allocate(rcol(norb),stat=ierr)
         call utils_alloc_check('dense_mat_sqrt','rcol',ierr)
      endif

      ! allocate workspace if none provided
      if (present(buf1) .and. present(buf2)) then
         work1 => buf1
         work2 => buf2
      else
         call dense_create(ibuf1,norb,norb,iscmplx=mat%iscmplx)
         call dense_create(ibuf2,norb,norb,iscmplx=mat%iscmplx)
         work1 => ibuf1
         work2 => ibuf2
      endif

      ! compute eigenvectors
      call dense_copy(work1,mat)
      call dense_normal_eigensolve(norb,eigvals,mat,work1)

      call utils_assert(all(eigvals > -tolerance), &
            'Error in dense_mat_sqrt: non positive definite matrix',minval(eigvals))
      where (eigvals < 0.0_DP) eigvals = 0.0_dP

      eigvals = sqrt(eigvals)

      ! multiply eigenvectors by sqrt(eigvals)
      if (mat%iscmplx) then
         do ii=1,norb
            call dense_get_col(zcol,work1,ii)
            zcol = eigvals(ii) * zcol
            call dense_put_col(zcol,work2,ii)
         enddo
      else
         do ii=1,norb
            call dense_get_col(rcol,work1,ii)
            rcol = eigvals(ii) * rcol
            call dense_put_col(rcol,work2,ii)
         enddo
      endif

      ! mat = zbuf1^H * zbuf2
      if (mat%iscmplx) then
         call dense_product(mat,work1,work2,opB='C')
      else
         call dense_product(mat,work1,work2,opB='T')
      endif

      ! destroy workspace
      deallocate(eigvals,stat=ierr)
      call utils_dealloc_check('dense_mat_sqrt','eigvals',ierr)

      if (mat%iscmplx) then
         deallocate(zcol,stat=ierr)
         call utils_dealloc_check('dense_mat_sqrt','zcol',ierr)
      else
         deallocate(rcol,stat=ierr)
         call utils_dealloc_check('dense_mat_sqrt','rcol',ierr)
      endif
      if (.not. (present(buf1) .and. present(buf2)) ) then
         call dense_destroy(ibuf1)
         call dense_destroy(ibuf2)
      endif

      nullify(work1)
      nullify(work2)

      return

   end subroutine dense_mat_sqrt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine converts a SPAM3 matrix to a dense matrix.                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dmat      (input) : The destination dense matrix                         !
  !   smat      (input) : The source sparse matrix                             !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_convert_sparsetodense(dmat_out,smat_in)

#ifdef SCALAPACK
    use sparse, only: SPAM3,sparse_spam3toblacs
#else
    use sparse, only: SPAM3,sparse_convert
#endif
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(DEM),intent(inout) :: dmat_out
    type(SPAM3),intent(in) :: smat_in

    ! Check compatability
    if (smat_in%iscmplx.and.(.not.dmat_out%iscmplx)) then
       call utils_abort('Error in &
            &dense_convert_densetosparse: source matrix is complex but dest &
            &matrix is real')
    end if

    ! Fill dense matrix values from SPAM3 matrix
#ifdef SCALAPACK
    if (dmat_out%iscmplx) then
       call sparse_spam3toblacs(dmat_out%zmtx,smat_in,dmat_out%blacs_ld, &
            dmat_out%blacs_ncol,dmat_out%blacs_desc)
    else
       call sparse_spam3toblacs(dmat_out%dmtx,smat_in,dmat_out%blacs_ld, &
            dmat_out%blacs_ncol,dmat_out%blacs_desc)
    end if
#else
    if (dmat_out%iscmplx) then
       call sparse_convert(dmat_out%zmtx,smat_in)
    else
       call sparse_convert(dmat_out%dmtx,smat_in)
    end if
#endif

  end subroutine dense_convert_sparsetodense


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine converts a dense matrix to a SPAM3 matrix.                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dmat      (input) : The destination sparse matrix                        !
  !   smat      (input) : The source dense matrix                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_convert_densetosparse(smat_out,dmat_in)

#ifdef SCALAPACK
    use sparse, only: SPAM3,sparse_blacstospam3
#else
    use sparse, only: SPAM3,sparse_convert
#endif
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: smat_out
    type(DEM),intent(in) :: dmat_in

    ! Check compatability
    if (smat_out%iscmplx.neqv.dmat_in%iscmplx) then
       call utils_abort('Error in &
            &dense_convert_densetosparse: matrices must be both real or both &
            &complex')
    end if
    ! Fill SPAM3 matrix values from dense matrix
#ifdef SCALAPACK
    if (dmat_in%iscmplx) then
       call sparse_blacstospam3(smat_out,dmat_in%zmtx,dmat_in%blacs_ld, &
            dmat_in%blacs_ncol,dmat_in%blacs_desc)
    else
       call sparse_blacstospam3(smat_out,dmat_in%dmtx,dmat_in%blacs_ld, &
            dmat_in%blacs_ncol,dmat_in%blacs_desc)
    end if
#else
    if (dmat_in%iscmplx) then
       call sparse_convert(smat_out,dmat_in%zmtx)
    else
       call sparse_convert(smat_out,dmat_in%dmtx)
    end if
#endif

  end subroutine dense_convert_densetosparse


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine copies one dense matrix to another.                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest     (output) : The destination dense matrix                         !
  !   src       (input) : The source dense matrix                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  ! Added argument to check complex2real conversion by Andrea Greco, Oct 2016. !
  !============================================================================!

  subroutine dense_copy(dest,src,cmplx_to_real)

    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: dest      ! The destination dense matrix
    type(DEM), intent(in) :: src          ! The source dense matrix
    ! agrecocmplx
    logical, optional, intent(in) :: cmplx_to_real ! check complex2real

    ! Local variables
    logical :: loc_cmplx_to_real         ! agrecocmplx

    loc_cmplx_to_real = .false.
    if (present(cmplx_to_real)) loc_cmplx_to_real = cmplx_to_real

    ! agrecocmplx: protection  against unsafe conversion from complex to real
    call utils_assert(.not.(src%iscmplx).or.(dest%iscmplx).or.loc_cmplx_to_real, &
         'ERROR: unsafe conversion of complex matrix to real. Call routine &
         &dense_copy with cmplx_to_real=.true. to force the conversion.')

    ! Check Arguments
    if ((dest%nrows /= src%nrows).or.(dest%mcols /= src%mcols)) then
       call utils_abort('Error in dense_copy: &
            &src and dest matrix sizes do not match')
    end if

    ! Trivial axpy of whole data arrays
    if (src%iscmplx) then
       if (dest%iscmplx) then
          dest%zmtx = src%zmtx
       else
          dest%dmtx(:,:) = real(src%zmtx(:,:),kind=DP)
       end if
    else
       if (dest%iscmplx) then
          dest%zmtx(:,:) = cmplx(src%dmtx(:,:),0.0_DP,kind=DP)
       else
          dest%dmtx = src%dmtx
       end if
    end if

  end subroutine dense_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This function calculates the largest imaginary part of a dense matrix      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The dense complex matrix to check                       !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, October 2016                                      !
  !============================================================================!

  real(kind=DP) function dense_max_imag_part(mat)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(in)         :: mat      ! input complex matrix

    call utils_assert(mat%iscmplx, 'ERROR in routine dense_max_imag_part: &
         &matrix must be complex')

    dense_max_imag_part = maxval(abs(aimag(mat%zmtx(:,:))))

  end function dense_max_imag_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs a safe conversion of a complex matrix to a real   !
  ! one. The imaginary part of the complex source matrix is checked against    !
  ! a specified threshold.                                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest  (output) : The resulting dense real matrix                         !
  !   src   (input)  : The original dense complex matrix                       !
  !   threshold      : max value allowed for the imaginary part of src         !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, September 2016                                    !
  !============================================================================!

  subroutine dense_take_real_part(dest,src,threshold)

    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(in)         :: src           ! input complex matrix
    type(DEM), intent(inout)      :: dest          ! output real matrix
    real(kind=DP), intent(in)     :: threshold     ! threshold to accept input matrix

    ! Local variables
    real(kind=DP)                 :: max_imag      ! global max imaginary part
    character(len=100)            :: message       ! error message

    call utils_assert(src%iscmplx, 'ERROR in routine dense_take_real_part: &
         &source matrix must be complex')
    call utils_assert(.not.(dest%iscmplx), 'ERROR in routine &
         &dense_take_real_part: destination matrix must be real')

    ! calculates largest imaginary part of source matrix
    max_imag = dense_max_imag_part(src)

    ! agrecocmplx: abort if imaginary part of complex source matrix
    ! is too large
    if (max_imag>threshold) then
       write(message,'(a,f15.8)') 'Error in routine dense_take_real_part: &
          &largest imaginary part of source is ', max_imag
       call utils_abort(message)
    ! agrecocmplx: call dense_copy with complx_to_real=.true. since
    ! it is safe to throw away the imaginary part
    else
       call dense_copy(dest,src,cmplx_to_real=.true.)
    end if

  end subroutine dense_take_real_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs a safe axpy operation involving a complex matrix  !
  ! ymat and a real one xmat. The imaginary part of xmat is checked against a  !
  ! specified threshold. The alpha parameter is real.                          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (output) : The resulting dense real matrix                         !
  !   xmat   (input)  : The original dense complex matrix                      !
  !   alpha  (input)  : y = y + alpha * x                                      !
  !   threshold      : max value allowed for the imaginary part of src         !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, October 2016                                      !
  !============================================================================!

  subroutine dense_safe_axpy_real(ymat,xmat,alpha,threshold)

    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(in)         :: xmat          ! input complex matrix
    type(DEM), intent(inout)      :: ymat          ! output real matrix
    real(kind=DP), intent(in)     :: alpha
    real(kind=DP), intent(in)     :: threshold     ! threshold to accept input matrix

    ! Local variables
    real(kind=DP)                 :: max_imag      ! global max imaginary part
    character(len=100)            :: message       ! error message

    call utils_assert(xmat%iscmplx, 'ERROR in routine dense_safe_axpy: &
         &x matrix must be complex')
    call utils_assert(.not.(ymat%iscmplx), 'ERROR in routine &
         &dense_safe_axpy: y matrix must be real')

    ! calculates largest imaginary part of source matrix
    max_imag = dense_max_imag_part(xmat)

    ! agrecocmplx: abort if imaginary part of complex source matrix
    ! is too large
    if (max_imag>threshold) then
       write(message,'(a,f15.8)') 'Error in routine dense_safe_axpy: &
          &largest imaginary part of x is ', max_imag
       call utils_abort(message)
    ! agrecocmplx: call dense_axpy with complx_to_real=.true. since
    ! it is safe to throw away the imaginary part
    else
       call dense_axpy_real(ymat,xmat,alpha,cmplx_to_real=.true.)
    end if

  end subroutine dense_safe_axpy_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs a safe axpy operation involving a complex matrix  !
  ! ymat and a real one xmat. The imaginary part of xmat is checked against a  !
  ! specified threshold. The alpha parameter is complex.                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (output) : The resulting dense real matrix                         !
  !   xmat   (input)  : The original dense complex matrix                      !
  !   alpha  (input)  : y = y + alpha * x                                      !
  !   threshold      : max value allowed for the imaginary part of src         !
  !----------------------------------------------------------------------------!
  ! Written by Andrea Greco, October 2016                                      !
  !============================================================================!

  subroutine dense_safe_axpy_complex(ymat,xmat,alpha,threshold)

    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(in)         :: xmat          ! input complex matrix
    type(DEM), intent(inout)      :: ymat          ! output real matrix
    complex(kind=DP), intent(in)  :: alpha
    real(kind=DP), intent(in)     :: threshold     ! threshold to accept input matrix

    ! Local variables
    real(kind=DP)                 :: max_imag      ! global max imaginary part
    character(len=100)            :: message       ! error message

    call utils_assert(xmat%iscmplx, 'ERROR in routine dense_safe_axpy: &
         &x matrix must be complex')
    call utils_assert(.not.(ymat%iscmplx), 'ERROR in routine &
         &dense_safe_axpy: y matrix must be real')

    ! calculates largest imaginary part of source matrix
    ! agrecocmplx: should we check the imaginary part of only xmat,
    ! or alpha*xmat? since alpha is complex, is it possible that
    ! xmat has non-zero imag part but alpha*xmat does?
    ! currently in sparse_axpy, alpha is multiplied by the real part
    ! of xmat, instead of taking the real part of alpha*xmat
    max_imag = dense_max_imag_part(xmat)

    ! agrecocmplx: abort if imaginary part of complex source matrix
    ! is too large
    if (max_imag>threshold) then
       write(message,'(a,f15.8)') 'Error in routine dense_safe_axpy: &
          &largest imaginary part of x is ', max_imag
       call utils_abort(message)
    ! agrecocmplx: call dense_axpy with complx_to_real=.true. since
    ! it is safe to throw away the imaginary part
    else
       call dense_axpy_complex(ymat,xmat,alpha,cmplx_to_real=.true.)
    end if

  end subroutine dense_safe_axpy_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The dense matrix y                                       !
  !   xmat  (input) : The dense matrix x                                       !
  !   alpha (input) : The real parameter alpha                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  ! Modified by Andrea Greco to allow safe conversion from complex xmat to     !
  ! real ymat, October 2016.                                                   !
  !============================================================================!

  subroutine dense_axpy_real(ymat,xmat,alpha,cmplx_to_real)

    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: ymat      ! The dense matrix y
    type(DEM), intent(in) :: xmat         ! The dense matrix x
    real(kind=DP), intent(in) :: alpha    ! The parameter alpha
    ! agrecocmplx
    logical, optional, intent(in) :: cmplx_to_real

    ! Local variables
    logical :: loc_cmplx_to_real ! agrecocmplx

    ! Check Arguments
    if ((ymat%nrows /= xmat%nrows).or.(ymat%mcols /= xmat%mcols)) then
       call utils_abort('Error in dense_axpy_real: &
            &xmat and ymat sizes do not match')
    end if

    ! agrecocmplx
    loc_cmplx_to_real = .false.
    if (present(cmplx_to_real)) loc_cmplx_to_real = cmplx_to_real

    ! agrecocmplx: protection  against unsafe conversion from complex to real
    call utils_assert(.not.(xmat%iscmplx).or.(ymat%iscmplx).or.loc_cmplx_to_real, &
         'ERROR: unsafe conversion of complex matrix to real. Call routine &
         &dense_axpy_real with cmplx_to_real=.true. to force the conversion.')

    ! Trivial axpy of whole data arrays
    if (xmat%iscmplx) then
       if (ymat%iscmplx) then
          ymat%zmtx = ymat%zmtx + alpha * xmat%zmtx
       else
          ymat%dmtx(:,:) = ymat%dmtx(:,:) + &
               alpha * real(xmat%zmtx(:,:),kind=DP)
       end if
    else
       if (ymat%iscmplx) then
          ymat%zmtx(:,:) = ymat%zmtx(:,:) + &
               cmplx(alpha * xmat%dmtx(:,:), 0.0_DP, kind=DP)
       else
          ymat%dmtx = ymat%dmtx + alpha * xmat%dmtx
       end if
    end if

  end subroutine dense_axpy_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The dense matrix y                                       !
  !   xmat  (input) : The dense matrix x                                       !
  !   alpha (input) : The complex parameter alpha                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  ! Modified by Andrea Greco to protect safe conversion from complex xmat      !
  ! to real ymat, October 2016.                                                !
  !============================================================================!

  subroutine dense_axpy_complex(ymat,xmat,alpha,cmplx_to_real)

    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: ymat      ! The dense matrix y
    type(DEM), intent(in) :: xmat         ! The dense matrix x
    complex(kind=DP), intent(in) :: alpha ! The parameter alpha
    ! agrecocmplx
    logical, optional, intent(in) :: cmplx_to_real

    ! Local arguments

    logical :: loc_cmplx_to_real         ! agrecocmplx
    ! Check Arguments
    if ((ymat%nrows /= xmat%nrows).or.(ymat%mcols /= xmat%mcols)) then
       call utils_abort('Error in dense_axpy_complex: &
            &xmat and ymat sizes do not match')
    end if

    ! agrecocmplx
    loc_cmplx_to_real = .false.
    if (present(cmplx_to_real)) loc_cmplx_to_real = cmplx_to_real

    ! agrecocmplx: protection against unsafe conversion from complex to real
    call utils_assert( ((aimag(alpha)==0.0_DP) .and. (.not.(xmat%iscmplx))) &
         .or. (ymat%iscmplx) .or. loc_cmplx_to_real, &
         'ERROR: unsafe conversion of complex matrix to real. Call routine &
         &dense_axpy_complex with cmplx_to_real=.true. to force the conversion.')

    ! Trivial axpy of data arrays
    if (xmat%iscmplx) then
       if (ymat%iscmplx) then
          ymat%zmtx = ymat%zmtx + alpha * xmat%zmtx
       else
          ymat%dmtx(:,:) = ymat%dmtx(:,:) + &
               real(alpha * xmat%zmtx(:,:),kind=DP)
       end if
    else
       if (ymat%iscmplx) then
          ymat%zmtx(:,:) = ymat%zmtx(:,:) + &
               alpha * cmplx(xmat%dmtx(:,:),0.0_DP,kind=DP)
       else
          ymat%dmtx = ymat%dmtx + real(alpha,kind=DP) * xmat%dmtx
       end if
    end if

  end subroutine dense_axpy_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs an rank 1 update operation on the matrices:       !
  !   A := A + alpha * x * y^T                                                 !
  ! selecting the columns ix of x and iy of y as the vectors for the update    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (inout) : The dense matrix A                                       !
  !   xmat  (inout) : The dense matrix x                                       !
  !   ymat  (input) : The dense matrix y                                       !
  !   ix    (input) : The row of matrix x to use as the updating row           !
  !   iy    (input) : The row of matrix y to use as the updating row           !
  !   alpha (input) : The real parameter alpha                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 6 February 2012.                                 !
  !============================================================================!

  subroutine dense_rank1_update_real(amat,xmat,ymat,ix,iy,alpha,skip_size_check)

    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: amat      ! The dense matrix y
    type(DEM), intent(in) :: xmat         ! The dense matrix x
    type(DEM), intent(in) :: ymat         ! The dense matrix y
    integer, intent(in) :: ix,iy          ! Updating rows
    real(kind=DP), intent(in) :: alpha    ! The parameter alpha
    logical, intent(in), optional :: skip_size_check

    ! Local Variables
    complex(kind=DP) :: alpha_cmplx
    logical :: loc_skip_size_check

    ! (P)BLAS subroutines
#ifdef SCALAPACK
    external :: PZGERU, PDGER
#else
    external :: ZGERU, DGER
#endif


    ! Check Arguments
    loc_skip_size_check = .false.
    if (present(skip_size_check)) loc_skip_size_check = skip_size_check
    if ((ymat%nrows /= xmat%nrows).or.(ymat%mcols /= xmat%mcols)) then
       if (.not.skip_size_check) then
          call utils_abort('Error in dense_rank1_update_real: &
               &xmat and ymat sizes do not match')
       end if
    end if
    if ((amat%iscmplx.neqv.xmat%iscmplx) .or. &
         (amat%iscmplx.neqv.ymat%iscmplx)) then
       call utils_abort('Error in dense_rank1_update_real: &
            &amat, xmat and ymat must be all real or all complex')
    end if
    if ((ix<0).or.(ix>xmat%nrows)) then
       call utils_abort('Error in dense_rank1_update_real: ix is out of bounds')
    end if
    if ((iy<0).or.(iy>ymat%nrows)) then
       call utils_abort('Error in dense_rank1_update_real: iy is out of bounds')
    end if

#ifdef SCALAPACK
    if (amat%iscmplx) then
       alpha_cmplx = cmplx(alpha,0.0_DP,kind=DP)
       call PZGERU(amat%nrows, amat%mcols, alpha_cmplx, &
            xmat%zmtx(1,1), 1,ix,xmat%blacs_desc, 1, &
            ymat%zmtx(1,1), 1,iy,ymat%blacs_desc, 1, &
            amat%zmtx(1,1), 1, 1,amat%blacs_desc)
    else
       call PDGER(amat%nrows, amat%mcols, alpha, &
            xmat%dmtx(1,1), 1,ix,xmat%blacs_desc, 1, &
            ymat%dmtx(1,1), 1,iy,ymat%blacs_desc, 1, &
            amat%dmtx(1,1), 1, 1,amat%blacs_desc)
    end if
#else
    if (amat%iscmplx) then
       alpha_cmplx = cmplx(alpha,0.0_DP,kind=DP)
       call ZGERU(amat%nrows, amat%mcols, alpha_cmplx, &
            xmat%zmtx(1,ix), 1, ymat%zmtx(1,iy), 1, &
            amat%zmtx(1,1), amat%nrows)
    else
       call DGER(amat%nrows, amat%mcols, alpha, &
            xmat%dmtx(1,ix), 1, ymat%dmtx(1,iy), 1, &
            amat%dmtx(1,1), amat%nrows)
    end if
#endif

  end subroutine dense_rank1_update_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs an rank 1 update operation on the matrices:       !
  !   A := A + alpha * x * y^T                                                 !
  ! selecting the columns ix of x and iy of y as the vectors for the update    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (inout) : The dense matrix A                                       !
  !   xmat  (inout) : The dense matrix x                                       !
  !   ymat  (input) : The dense matrix y                                       !
  !   ix    (input) : The row of matrix x to use as the updating row           !
  !   iy    (input) : The row of matrix y to use as the updating row           !
  !   alpha (input) : The complex parameter alpha                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 6 February 2012.                                 !
  !============================================================================!

  subroutine dense_rank1_update_complex(amat,xmat,ymat,ix,iy,alpha)

    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: amat      ! The dense matrix y
    type(DEM), intent(in) :: xmat         ! The dense matrix x
    type(DEM), intent(in) :: ymat         ! The dense matrix y
    integer, intent(in) :: ix,iy          ! Updating rows
    complex(kind=DP), intent(in) :: alpha ! The parameter alpha

    ! (P)BLAS subroutines
#ifdef SCALAPACK
    external :: PZGERU
#else
    external :: ZGERU
#endif

    ! Check Arguments
    if ((ymat%nrows /= xmat%nrows).or.(ymat%mcols /= xmat%mcols)) then
       call utils_abort('Error in dense_rank1_update_real: &
            &xmat and ymat sizes do not match')
    end if
    if ((amat%iscmplx.neqv.xmat%iscmplx) .or. &
         (amat%iscmplx.neqv.ymat%iscmplx)) then
       call utils_abort('Error in dense_rank1_update_real: &
            &amat, xmat and ymat must be all real or all complex')
    end if
    if (.not.amat%iscmplx) then
       call utils_abort('Error in dense_rank1_update_complex: &
            &cannot update real amat with complex alpha')
    end if
    if ((ix<0).or.(ix>xmat%nrows)) then
       call utils_abort('Error in dense_rank1_update_real: ix is out of bounds')
    end if
    if ((iy<0).or.(iy>ymat%nrows)) then
       call utils_abort('Error in dense_rank1_update_real: iy is out of bounds')
    end if

#ifdef SCALAPACK
    call PZGERU(amat%nrows, amat%mcols, alpha, &
         xmat%zmtx(1,1), 1,ix,xmat%blacs_desc, 1, &
         ymat%zmtx(1,1), 1,iy,ymat%blacs_desc, 1, &
         amat%zmtx(1,1), 1, 1,amat%blacs_desc)
#else
    call ZGERU(amat%nrows, amat%mcols, alpha, &
         xmat%zmtx(1,iy), 1, ymat%zmtx(1,iy), 1, &
         amat%zmtx(1,1), amat%nrows)
#endif

  end subroutine dense_rank1_update_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine writes a dense matrix to file for debug purposes           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be written                                !
  !----------------------------------------------------------------------------!
  ! Written by Simon Dubois, Oct. 2013                                         !
  ! Modified by David Turban to enable writing to stdout, July 2015            !
  !============================================================================!

  subroutine dense_write(mat,filename,to_stdout)

    use constants, only: stdout
    use utils, only: utils_abort, utils_assert, utils_unit

    implicit none

    ! Arguments
    type(DEM), intent(in)          :: mat
    character(len=*), intent(in)   :: filename
    logical, optional, intent(in)  :: to_stdout

    ! Local Variables
    complex(kind=DP) :: zm
    real(kind=DP) :: dm
    integer :: m, n
    integer :: iel, jel
    integer :: io_unit1, io_unit2
    character(len=80) :: filemat
    logical :: loc_to_stdout

    if (present(to_stdout)) then
       loc_to_stdout = to_stdout
    else
       loc_to_stdout = .false.
    end if

    write(filemat,'(a,a)') filename,'.mat'

    io_unit1 = utils_unit()
    open(unit=io_unit1, form="formatted", file=filename, action="write")
    io_unit2 = utils_unit()
    open(unit=io_unit2, form="formatted", file=filemat, action="write")

    n = mat%nrows
    m = mat%mcols
    write(io_unit1,'(a,l1)') "iscmplx      : ", mat%iscmplx
    write(io_unit1,'(a,2(i4.4,1x))') "nrows, mcols : ", n, m

    if (mat%iscmplx) then
       do iel=1,n
          do jel=1,m
             call dense_get_element(zm,mat,iel,jel)
             write(io_unit1,'(2(i4.4,1x),2f12.8)') iel, jel, zm
             write(io_unit2,'(2f12.8,2x)',advance='no') zm
             if(loc_to_stdout) write(stdout,'(2f12.8,2x)',advance='no') zm
          enddo
          write(io_unit2,*)
          if(loc_to_stdout) write(stdout,*)
       enddo
    else
       do iel=1,n
          do jel=1,m
             call dense_get_element(dm,mat,iel,jel)
             write(io_unit1,'(2(i4.4,1x),f12.8)') iel, jel, dm
             write(io_unit2,'(f12.8,2x)',advance='no') dm
             if(loc_to_stdout) write(stdout,'(2f12.8,2x)',advance='no') dm
          enddo
          write(io_unit2,*)
          if(loc_to_stdout) write(stdout,*)
       enddo
    endif

    close(io_unit1)
    close(io_unit2)

  end subroutine dense_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This function computes the trace of a matrix (or product) in dense format. !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   A      (input)  : The matrix on which to perform the trace.              !
  !   B      (input)  : The optional matrix to perform tr(A*B).                !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Nov. 2014                                        !
  !============================================================================!

  function dense_trace(A,B) result(trace)
    use comms, only: comms_reduce, comms_bcast
#ifdef SCALAPACK
    use comms, only: pub_on_root
    !use mpi, only: mpi_wtime
#endif
    use utils, only: utils_abort, utils_assert
    implicit none
    type(DEM),           intent(inout) :: A
    type(DEM), optional, intent(inout) :: B
    real(kind=dp)                      :: trace

    integer       :: N, M, i
#ifdef SCALAPACK
    integer       :: row,col,indx,indy
    real(kind=dp) :: tmp
    type(DEM)     :: tmp_dem
    integer       :: proc
    logical,save  :: profile=.false.
    logical,save  :: element_wise=.true.
    real(kind=dp) :: bef_time, el_time, den_time
    real(kind=dp), external :: pdlatra
    integer,       external :: blacs_pnum
    external :: pddot, infog2l ! PBLAS & ScaLAPACK subroutines
#endif

    call utils_assert(.not.A%iscmplx, 'Error in dense_trace : complex &
         & matrices are not currently supported. (A)')

    N=A%nrows
    if(present(B)) then
       M=B%mcols
       call utils_assert(.not.B%iscmplx, 'Error in dense_trace : complex &
            & matrices are not currently supported. (B)')
    else
       M=A%mcols
    end if
    call utils_assert(M==N, 'Error in dense_trace : the number of rows and &
         &columns of the matrix to be traced must be the same.')

    if(present(B)) then
#ifdef SCALAPACK
       if(profile.or.element_wise) then
          !if(profile) bef_time=mpi_wtime()
          trace = 0.0_dp
          do i=1,N
             tmp=0.0_dp
             call pddot(N, tmp, A%dmtx, i, 1, A%blacs_desc, N, B%dmtx, 1, i, B%blacs_desc, 1)
             call infog2l(i,1,A%blacs_desc,blacs_nprow,blacs_npcol, &
                  & blacs_myrow,blacs_mycol,indx,indy,row,col)
             proc = blacs_pnum(blacs_context,row,col)
             call comms_bcast(proc,tmp)

             trace = trace + tmp
          end do
          !if(profile) el_time=mpi_wtime()-bef_time
       end if
       if(profile.or..not.element_wise) then
          !if(profile) bef_time=mpi_wtime()
          call dense_create(tmp_dem,A)
          call dense_product(tmp_dem,A,B)
          tmp = pdlatra(N,tmp_dem%dmtx,1,1,tmp_dem%blacs_desc)
          call dense_destroy(tmp_dem)
          if(profile) then
          !   den_time=mpi_wtime()-bef_time
          else
             trace=tmp
          end if
       end if

       if(profile) then
          profile=.false.
          !element_wise=(den_time>el_time)
          !if(pub_on_root) write(*,*) "trace profile : element_wise = ",element_wise,", vals = ", tmp,trace
       end if
#else
       trace = 0.0_dp
       do i=1,N
          trace=trace+dot_product(A%dmtx(i,:),B%dmtx(:,i))
       end do
#endif

    else
#ifdef SCALAPACK
       trace = pdlatra(N,A%dmtx,1,1,A%blacs_desc)
#else
       trace = 0.0_dp
       do i=1,N
          trace=trace+A%dmtx(i,i)
       end do
#endif
    end if

  end function dense_trace

! subroutine dense_reduce(scope,A,reduced,proc,ws_query)
!   implicit none
!   character(len=1),                intent(in)  :: op
!   type(DEM),                       intent(in)  :: A
!   real(kind=dp),     dimension(:), intent(out) :: reduced
!   logical, optional,               intent(in)  :: distrib
!   integer, optional,               intent(out) :: ws_query ! Workspace query
!
!   integer :: m,n
!   integer :: lproc
!
!   m = A%nrows
!   n = A%mcols
!   lproc=0
!   if(present(distrib)) lproc=-1
!
!   call utils_assert(op=="c".or.op=="C".or.op=="r".or."R".or. &
!        & op=="a".or.op=="A".or."g".or."G",'Error in dense_reduce :&
!        & op must equal A, C, G or R')
!   if(op=="c") op="C"
!   if(op=="r") op="R"
!   if(op=="A") op="A"
!   if(op=="g") op="G"
!
!   call blacs_gridinfo(blacs_context, nprow, npcol, myrow, mycol)
!   iroffa = mod(0, A%blacs_nb) ! offsets 0 because we're starting at first element.
!   icoffa = mod(0, A%blacs_nb) ! 0
!
!   if(op=="C") loc_lwork = numroc(n+icoffa, A%blacs_nb, mycol, iacol, npcol)
!   if(op=="R") loc_lwork = numroc(m+iroffa, A%blacs_nb, myrow, iarow, nprow)
!   if(op=="G") loc_lwork = min(numroc(n+icoffa, A%blacs_nb, mycol, iacol, npcol),&
!        &numroc(m+iroffa, A%blacs_nb, myrow, iarow, nprow))
!   if(op=="A") loc_lwork = 1
!
!   call infog2l(1,1,A%blacs_desc,nprow,npcol, &
!        & myrow,mycol,indx,indy,row,col)
!
!   lda = A%blacs_ld
!
!
!   iarow = indxg2p(1, A%blacs_nb, myrow, row, nprow)
!   iacol = indxg2p(1, A%blacs_nb, mycol, col, npcol)
!   mp = numroc(m+iroffa, A%blacs_nb, myrow, iarow, nprow)
!   nq = numroc(n+icoffa, A%blacs_nb, mycol, iacol, npcol)
!   if(myrow.eq.iarow) mp = mp-iroffa
!   if(mycol.eq.iacol) nq = nq-icoffa
!
!   if(op=="C") then
!      ! The following 2-way nested loop is completely local and hence
!      ! unrollable / OpenMP-able for anyone who might care enough to
!      ! do the optimisation.
!      if(nq>0) then
!         do j = indy,indy+nq-1
!            loc_sum=0.0_dp
!            if(mp>0) then
!               do i = indx,mp+indx-1
!                  loc_sum = loc_sum + abs(A%dmtx(i,j))
!               end do
!            end if
!            loc_sums(j-indy+1) = loc_sum
!         end do
!      end if
!
!      ! Find sum of global matrix columns and store on row 0 of
!      ! process grid
!
!      call dgsum2d(blacs_context, 'C', ' ', 1, nq, loc_sums, 1, &
!           & 0, mycol)
!
!   elseif(op=="R".or.op=="G") then
!      ! The following 2-way nested loop is completely local and hence
!      ! unrollable / OpenMP-able for anyone who might care enough to
!      ! do the optimisation.
!      if(mp>0) then
!         do i = indx,indx+mp-1
!            loc_sum=0.0_dp
!            if(nq>0) then
!               do j = indy,nq+indy-1
!                  loc_sum = loc_sum + abs(A%dmtx(i,j))
!               end do
!            end if
!            loc_sums(i-indx+1) = loc_sum
!         end do
!      end if
!
!      ! Find sum of global matrix columns and store on row 0 of
!      ! process grid
!
!      call dgsum2d(blacs_context, 'R', ' ', mp, 1, loc_sums, max(1,mp), &
!           & myrow, 0)
!
!   end if
!
!
!   if(op=="G") then
!
!      jn = min(A%blacs_nb, n)
!      loci = indx
!      locj = indy
!      !
!      !     handle first diagonal block separately
!      !
!      if(myrow==row .and. mycol==col ) then
!         do ll = 0, jn-1
!            loc_diag(loci+ll) = A(loci+ll,locj+ll)
!         end do
!      end if
!      if(myrow.eq.row) loci = loci + jn
!      if(mycol.eq.col) locj = locj + jn
!      row = mod(row+1, nprow)
!      col = mod(col+1, npcol)
!
!      !
!      !     loop over the remaining block of columns
!      !
!      do j = jn+1, n, A%blacs_nb
!         jb = min(1+n-j, A%blacs_nb)
!         !
!         if(myrow==row .and. mycol==col) then
!            do ll = 0, jb-1
!               loc_diag(loci+ll) = A(loci+ll,locj+ll)
!            end do
!         end if
!         if(myrow==row) loci = loci + jb
!         if(mycol==col) locj = locj + jb
!         row = mod(row+1, nprow)
!         col = mod(col+1, npcol)
!      end do
!
!   end if
!
! end subroutine dense_reduce

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This function computes a norm of a matrix.                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   op     (input)  : The particular norm (1=1-norm, I=inf-norm, F=Frob-norm)!
  !   A      (input)  : The matrix                                             !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Nov. 2014                                        !
  !============================================================================!

  function dense_norm(op,A) result(norm)
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check
    implicit none
    character(len=1), intent(in)    :: op
    type(DEM),        intent(inout) :: A
    real(kind=dp)                   :: norm
    real(kind=dp)                   :: tmp
    integer                         :: m,n,i
#ifdef SCALAPACK
    real(kind=dp),   allocatable    :: work(:)
    integer                         :: lwork
    integer                         :: nprow,npcol,myrow,mycol,iroffa,icoffa,indx,indy,row,col, iarow,iacol
    integer :: ierr
#endif
    character(len=1)                :: lop
#ifdef SCALAPACK
    real(kind=dp),    external      :: pdlange, pzlange
    integer,          external      :: indxg2p, numroc
    external :: blacs_gridinfo, infog2l ! BLACS & ScaLAPACK subroutines
#endif

    m = A%nrows
    n = A%mcols

    call utils_assert(op=="I".or.op=="1".or.op=="F".or.op=="2".or.op=="i".or.op=="f",'Error in dense_norm :&
         & op must equal I, 1, F or 2')
    lop=op
    if(lop=="f") lop="F"
    if(lop=="i") lop="I"
    if(lop=="2") lop="F"

#ifdef SCALAPACK
    if(lop=="I".or.lop=="1") then
       call blacs_gridinfo(blacs_context, nprow, npcol, myrow, mycol)
       iroffa = mod(0, A%blacs_nb)
       icoffa = mod(0, A%blacs_nb)
       call infog2l(1,1,A%blacs_desc,nprow,npcol, &
            & myrow,mycol,indx,indy,row,col)
       iarow = indxg2p(1, A%blacs_nb, myrow, row, nprow)
       iacol = indxg2p(1, A%blacs_nb, mycol, col, npcol)
       if(lop=="I") lwork = numroc(m+iroffa, A%blacs_nb, myrow, iarow, nprow)
       if(lop=="1") lwork = numroc(n+icoffa, A%blacs_nb, mycol, iacol, npcol)
    else
       lwork=1
    end if
    allocate(work(lwork),stat=ierr)
    call utils_alloc_check('dense_norm','work',ierr)
    if(A%iscmplx) then
       norm=pzlange(lop,m,n,A%zmtx,1,1,A%blacs_desc,work)
    else
       norm=pdlange(lop,m,n,A%dmtx,1,1,A%blacs_desc,work)
    end if

    deallocate(work,stat=ierr)
    call utils_dealloc_check('dense_norm','work',ierr)
#else
    norm=0.0_dp
    select case(lop)
    case("F")
       norm=sqrt(dense_trace(A,A))
    case("1")
       if(A%iscmplx) then
          do i=1,n
             tmp=sum(abs(A%zmtx(:,i)))
             if(tmp>norm) norm=tmp
          end do
       else
          do i=1,n
             tmp=sum(abs(A%dmtx(:,i)))
             if(tmp>norm) norm=tmp
          end do
       end if
    case("I")
       if(A%iscmplx) then
          do i=1,n
             tmp=sum(abs(A%zmtx(i,:)))
             if(tmp>norm) norm=tmp
          end do
       else
          do i=1,n
             tmp=sum(abs(A%dmtx(i,:)))
             if(tmp>norm) norm=tmp
          end do
       end if
    end select
#endif

  end function dense_norm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This function computes a norm of a vector.                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   v      (input)  : The vector                                             !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Nov. 2014                                        !
  !============================================================================!

  function dense_norm_vec(v) result(norm)
    use comms, only: comms_bcast, comms_barrier
    use utils, only: utils_assert
    implicit none
    type(DEM),        intent(in) :: v
    real(kind=dp)                :: norm
    integer                      :: m
#ifdef SCALAPACK
    integer                      :: nprow,npcol,myrow,mycol,iroffa,icoffa,indx,indy,row,col, iarow,iacol
    external   :: pdznrm2, pdnrm2
!    integer,          external   :: blacs_pnum
#endif

    m = v%nrows
    call utils_assert (v%mcols==1, 'Error in dense_norm_vec: &
         &the input dense object is not a vector.')

#ifdef SCALAPACK
    call comms_barrier()

    if(v%iscmplx) then
       call pdznrm2(m, norm, v%zmtx, 1, 1, v%blacs_desc, 1)
    else
       call pdnrm2(m, norm, v%dmtx, 1, 1, v%blacs_desc, 1)
    end if
!    call blacs_gridinfo(blacs_context, nprow, npcol, myrow, mycol)
!    call infog2l(1,1,v%blacs_desc,blacs_nprow,blacs_npcol, &
!         & blacs_myrow,blacs_mycol,indx,indy,row,col)
!    proc = blacs_pnum(blacs_context,row,col)
!    call comms_bcast(proc,norm)
    call comms_barrier()

#else
    ! Shouldn't need to reduce anything if we don't have distributed vectors... every core should have the
    ! same data.
    if(v%iscmplx) then
       norm=sqrt(sum(real(conjg(v%zmtx)*v%zmtx)))
    else
       norm=sqrt(sum(v%dmtx**2))
    end if
#endif

  end function dense_norm_vec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This function computes a dot product of two vectors.                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   v      (input)  : A vector                                               !
  !   u      (input)  : A vector                                               !
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Nov. 2014                                        !
  !============================================================================!

  function dense_dot_product(v,u) result(dot)
    use comms, only: comms_bcast, comms_barrier, comms_reduce, pub_root_proc_id
    use utils, only: utils_assert
    implicit none
    type(DEM),        intent(in) :: v, u
    real(kind=dp)                :: dot
    integer                      :: m,n
#ifdef SCALAPACK
    integer                      :: nprow,npcol,myrow,mycol,iroffa,icoffa,indx,indy,row,col, iarow,iacol
    real(kind=DP) :: dot_aux(1,1)
#endif
    real(kind=dp),    external   :: ddot,zdot
#ifdef SCALAPACK
    integer,          external   :: blacs_pnum
    external :: dgsum2d, blacs_pcoord ! BLACS subroutines
#endif

    m = v%nrows
    n = v%mcols ! 1


    call utils_assert(.not.(v%iscmplx.or.u%iscmplx),'Error in dense_dot_product :&
         & complex vectors not yet supported!')

#ifdef SCALAPACK
    call comms_barrier()


    ! This operation was very slow on Iridis4 and Polaris, having tried many different
    ! options. Finally settled on this as the best of a bad bunch. Seems fine on Archer.

    if(v%iscmplx) then
!       if(allocated(v%zmtx).and.allocated(u%zmtx)) then
!          if(size(v%zmtx)>0 .and. size(u%zmtx)>0) then
!             !          dot=dot_product(v%dmtx(:,1),u%dmtx(:,1))
!             dot=zdot(size(v%zmtx,1),v%zmtx(1,1),1,u%zmtx(1,1),1)
!          else
!             dot=0.0_dp
!          end if
!       else
!          dot=0.0_dp
!       end if
    else
       if(allocated(v%dmtx).and.allocated(u%dmtx)) then
          if(size(v%dmtx)>0 .and. size(u%dmtx)>0) then
             !          dot=dot_product(v%dmtx(:,1),u%dmtx(:,1))
             dot=ddot(size(v%dmtx,1),v%dmtx(1,1),1,u%dmtx(1,1),1)
          else
             dot=0.0_dp
          end if
       else
          dot=0.0_dp
       end if
    end if

!    call blacs_gridinfo(blacs_context, nprow, npcol, myrow, mycol)
!    call infog2l(1,1,v%blacs_desc,blacs_nprow,blacs_npcol, &
!         & blacs_myrow,blacs_mycol,indx,indy,row,col)
!    proc = blacs_pnum(blacs_context,row,col)
    call blacs_pcoord(blacs_context, pub_root_proc_id, row,col)

    call dgsum2d(blacs_context, 'C', ' ', 1, 1, dot_aux, 1, row, col)
    dot = dot_aux(1,1)

!    tmp=dot
!    dot=0.0_dp
!
 !  call pddot(m, dot, v%dmtx, 1, 1, v%blacs_desc, &
 !       & 1, u%dmtx, 1, 1, u%blacs_desc, 1)
!   call comms_reduce('SUM',dot,root=pub_root_proc_id)

!
!    if(pub_on_root) write(*,*) "dense_dot_product : ", row,col,dot


!    call blacs_gridinfo(blacs_context, nprow, npcol, myrow, mycol)
!    call infog2l(1,1,v%blacs_desc,blacs_nprow,blacs_npcol, &
!         & blacs_myrow,blacs_mycol,indx,indy,row,col)
!    proc = blacs_pnum(blacs_context,row,col)
!    call comms_bcast(proc,dot)

    call comms_barrier()
#else
    if(v%iscmplx) then
!       dot=zdot(size(v%zmtx),v%zmtx(1,1),1,u%zmtx(1,1),1)
!       dot=dot_product(v%zmtx,u%zmtx)
    else
       dot=ddot(size(v%dmtx),v%dmtx(1,1),1,u%dmtx(1,1),1)
!       dot=dot_product(v%dmtx,u%dmtx)
    end if
#endif

  end function dense_dot_product

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine computes u = u + a*v where a is a scalar.                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   v      (input)   : A vector                                              !
  !   u      (output)  : A vector                                              !
  !----------------------------------------------------------------------------!
  ! Caveats : a is real, whether u and v are or not. I guess this should be    !
  !           converted into an optional real and an optional complex argument.!
  !----------------------------------------------------------------------------!
  ! Written by Jolyon Aarons, Nov. 2014                                        !
  !============================================================================!

  subroutine dense_vec_axpy(u,a,v)
    use utils, only: utils_assert
    implicit none
    type(DEM),        intent(in) :: v
    type(DEM),        intent(inout) :: u
    real(kind=dp),    intent(in) :: a
    integer                      :: m,n
#ifdef SCALAPACK
    integer                      :: nprow,npcol,myrow,mycol,iroffa,icoffa,indx,indy,row,col, iarow,iacol
    integer,          external   :: blacs_pnum
    external :: pzaxpy, pdaxpy !! ScaLAPACK subroutines
#endif

    call utils_assert(v%iscmplx .eqv. u%iscmplx,'Error in dense_vec_axpy :&
         & u and v must both be or not be complex!')

    m = v%nrows
    n = v%mcols ! 1

#ifdef SCALAPACK
    if(u%iscmplx) then
       call pzaxpy(m, cmplx(a,0.0_dp,kind=DP), v%zmtx, 1, 1, v%blacs_desc, 1, u%zmtx, 1, 1, u%blacs_desc, 1)
    else
       call pdaxpy(m, a, v%dmtx, 1, 1, v%blacs_desc, 1, u%dmtx, 1, 1, u%blacs_desc, 1)
    end if
#else
    if(u%iscmplx) then
       u%zmtx=u%zmtx + cmplx(a,0.0_dp,kind=DP)*v%zmtx
    else
       u%dmtx=u%dmtx + a*v%dmtx
    end if
#endif

  end subroutine dense_vec_axpy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine transposes a matrix in dense format                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat   (input)  : The matrix to be transposed                            !
  !   bmat   (output) : The transposed matrix                                  !
  !----------------------------------------------------------------------------!
  ! Written by Simon Dubois, Oct. 2013                                         !
  ! Modifed by Jolyon Aarons, Nov. 2014 to use pdtran.                         !
  !============================================================================!

  subroutine dense_transpose(bmat,amat)
!    use constants, only: cmplx_0, cmplx_1
    use utils, only: utils_abort, utils_assert

    implicit none

    ! Arguments
    type(DEM), intent(in) :: amat
    type(DEM), intent(inout) :: bmat

#ifdef SCALAPACK
    ! PBLAS subroutine
    external :: pdtran
#endif

    ! Local Variables
    integer :: m, n
    integer :: iel, jel
    real(kind=DP) :: dm
    complex(kind=DP) :: zm

    ! Check matrices all real or all complex
    if (bmat%iscmplx.neqv.amat%iscmplx) then
       call utils_abort('Error in dense_transpose: &
            &both matrices must be real or complex')
    end if

    ! Check matrices sizes
    m = amat%nrows
    n = amat%mcols
    call utils_assert(bmat%nrows == n, 'Error in dense_transpose: &
            &amat%mcols /= bmat%nrows')
    call utils_assert(bmat%mcols == m, 'Error in dense_transpose: &
            &amat%nrows /= bmat%mcols')

 ! jd: @fixme
 !     Temporarily using a naive transpose, will revert to using pdtran
 !     once Jolyon identifies the pdos bug.

    ! Transpose amat
    if (amat%iscmplx) then
       do iel=1,m
          do jel=1,n
             call dense_get_element(zm,amat,iel,jel)
             call dense_put_element(zm,bmat,jel,iel)
          enddo
       enddo
       ! Note this is a transpose, not a conjugate transpose (pztranc)
       ! in order to maintain compatibility.
!       call pztranu(m,n,cmplx_1, amat%zmtx,1,1,amat%blacs_desc,cmplx_0,bmat%zmtx,1,1,bmat%blacs_desc)

    else

#ifndef SCALAPACK
       do iel=1,m
          do jel=1,n
             call dense_get_element(dm,amat,iel,jel)
             call dense_put_element(dm,bmat,jel,iel)
          enddo
       enddo
#else
       call pdtran(n,m,1.0_dp, amat%dmtx,1,1,amat%blacs_desc,0.0_dp,bmat%dmtx,1,1,bmat%blacs_desc)
#endif

    endif

  end subroutine dense_transpose

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, March 2010 based on SPAM3 version.               !
  !============================================================================!

  subroutine dense_scale_real(mat,alpha,beta)

    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(inout)            :: mat    ! The matrix to be operated on
    real(kind=DP), intent(in)           :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local Variables
    integer :: ielem

#ifdef SCALAPACK
    external :: infog2l ! ScaLAPACK subroutine
    integer :: nprow, npcol, myrow, mycol, iia, jja, iarow, iacol
#endif

    ! Rescale the elements if required
    if (alpha /= 1.0_DP) then
       if (mat%iscmplx) then
          mat%zmtx = alpha * mat%zmtx
       else
          mat%dmtx = alpha * mat%dmtx
       end if
    end if

    ! Shift the eigenvalues if required
    if (present(beta)) then

       ! Can only shift eigenvalues of square matrix
       if (mat%mcols/=mat%nrows) then
          call utils_abort('Error in dense_scale_real: &
               &cannot shift eigenvalues of non-square matrices')
       end if

       ! Add identity matrix scaled by beta to all elements
! ja531 -> Added scalapack version for when beta is present.
#ifdef SCALAPACK
       if (mat%iscmplx) then
          do ielem=1,mat%nrows
             call infog2l(ielem, ielem, mat%blacs_desc, blacs_nprow, blacs_npcol, blacs_myrow, blacs_mycol, &
                  & iia, jja, iarow, iacol)
             if( blacs_myrow.eq.iarow .and. blacs_mycol.eq.iacol ) then
                mat%zmtx(iia,jja) = mat%zmtx(iia,jja) + cmplx(beta,0.0_dp,kind=dp)
             end if
          end do
       else
          do ielem=1,mat%nrows
             call infog2l(ielem, ielem, mat%blacs_desc, blacs_nprow, blacs_npcol, blacs_myrow, blacs_mycol, &
                  & iia, jja, iarow, iacol)
!                   write(*,*) pub_my_proc_id, "HERE",myrow,iarow,mycol,iacol
             if( blacs_myrow.eq.iarow .and. blacs_mycol.eq.iacol ) then
!                write(*,*) "==>", pub_my_proc_id, ielem
                mat%dmtx(iia,jja) = mat%dmtx(iia,jja) + beta
             end if
          end do
       end if
#else
       if (mat%iscmplx) then
          do ielem=1,mat%nrows
             mat%zmtx(ielem,ielem) = mat%zmtx(ielem,ielem) + &
                  cmplx(beta,0.0_DP,kind=DP)
          end do
       else
          do ielem=1,mat%nrows
             mat%dmtx(ielem,ielem) = mat%dmtx(ielem,ielem) + beta
          end do
       end if
#endif
    end if

  end subroutine dense_scale_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, March 2010 based on SPAM3 version.               !
  !============================================================================!

  subroutine dense_scale_complex(mat,alpha,beta)

    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(inout)            :: mat    ! The matrix to be operated on
    complex(kind=DP), intent(in)        :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local Variables
    integer :: ielem

    ! This only makes sense if mat is complex...
    if (.not.mat%iscmplx) then
       call utils_abort('Error in sparse_scale_complex: mat must be complex.')
    end if

    ! Rescale the elements if required
    if (alpha /= (1.0_DP,0.0_DP)) then
       if (mat%iscmplx) then
          mat%zmtx = alpha * mat%zmtx
       end if
    end if

    ! Shift the eigenvalues if required
    if (present(beta)) then

       ! Can only shift eigenvalues of square matrix
       if (mat%mcols/=mat%nrows) then
          call utils_abort('Error in dense_scale_real: &
               &cannot shift eigenvalues of non-square matrices')
       end if

       ! Add identity matrix scaled by beta to all elements
       if (mat%iscmplx) then
          do ielem=1,mat%nrows
             mat%zmtx(ielem,ielem) = mat%zmtx(ielem,ielem) + &
                  cmplx(beta,0.0_DP,kind=DP)
          end do
       end if
    end if

  end subroutine dense_scale_complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The complex shift parameter                             !
  !----------------------------------------------------------------------------!
  ! Modified from dense_scale_complex by Jolyon Aarons, 2017.                  !
  ! Written by Nicholas Hine, March 2010 based on SPAM3 version.               !
  !============================================================================!

  subroutine dense_scale_complex2(mat,alpha,beta)
    use constants, only: cmplx_1
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(DEM), intent(inout)            :: mat    ! The matrix to be operated on
    complex(kind=DP), intent(in)        :: alpha  ! The scaling parameter
    complex(kind=DP), intent(in) :: beta   ! The shift parameter

    ! Local Variables
    integer :: ielem

    ! This only makes sense if mat is complex...
    if (.not.mat%iscmplx) then
       call utils_abort('Error in sparse_scale_complex: mat must be complex.')
    end if

    ! Rescale the elements if required
    if (alpha /= cmplx_1) then
       if (mat%iscmplx) then
          mat%zmtx = alpha * mat%zmtx
       end if
    end if

    ! Shift the eigenvalues if required
!    if (present(beta)) then

       ! Can only shift eigenvalues of square matrix
       if (mat%mcols/=mat%nrows) then
          call utils_abort('Error in dense_scale_real: &
               &cannot shift eigenvalues of non-square matrices')
       end if

       ! Add identity matrix scaled by beta to all elements
       if (mat%iscmplx) then
          do ielem=1,mat%nrows
             mat%zmtx(ielem,ielem) = mat%zmtx(ielem,ielem) + &
                  beta
          end do
       end if
!    end if

  end subroutine dense_scale_complex2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine retrieves an element from a given dense matrix, where the  !
  ! element can be local to any processor.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el  (output)   : The element of the real matrix retrieved                !
  !   mat  (input)   : The real dense matrix                                   !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_get_element_real(el,mat,jrow,icol)

    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: el          ! The element to insert
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

#ifdef SCALAPACK
    ! ScaLAPACK subroutine
    external :: pdelget
#endif

    ! Check matrix is real
    if (mat%iscmplx) then
       call utils_abort('Error in dense_get_element_real: mat must be real.')
    end if

    ! Get element directly from matrix, or call pdelget
#ifdef SCALAPACK
    call pdelget('A',' ',el,mat%dmtx,jrow,icol,mat%blacs_desc)
#else
    el = mat%dmtx(jrow,icol)
#endif

  end subroutine dense_get_element_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine inserts an element into a given dense matrix, where the    !
  ! element can be local to any processor.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el   (input)   : The element to insert into the real matrix              !
  !   mat  (inout)   : The real dense matrix                                   !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_put_element_real(el,mat,jrow,icol)

    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: el           ! The element to insert
    type(DEM), intent(inout) :: mat           ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

#ifdef SCALAPACK
    ! ScaLAPACK subroutine
    external :: pdelset
#endif

    ! Check matrix is real
    if (mat%iscmplx) then
       call utils_abort('Error in dense_put_element_real: mat must be real.')
    end if

    ! Put element directly in matrix, or call pdelset
#ifdef SCALAPACK
    call pdelset(mat%dmtx,jrow,icol,mat%blacs_desc,el)
#else
    mat%dmtx(jrow,icol) = el
#endif

  end subroutine dense_put_element_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine retrieves an element from a given dense matrix, where the  !
  ! element can be local to any processor.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el  (output)   : The element of the complex matrix retrieved             !
  !   mat  (input)   : The complex dense matrix                                !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_get_element_complex(el,mat,jrow,icol)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: el       ! The element to insert
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

#ifdef SCALAPACK
    ! ScaLAPACK subroutine
    external :: pzelget
#endif

    ! Check matrix is complex
    call utils_assert(mat%iscmplx,'Error in dense_get_element_complex: &
            &mat must be complex.')

    ! Get element directly from matrix, or call pdelget
#ifdef SCALAPACK
    call pzelget('A',' ',el,mat%zmtx,jrow,icol,mat%blacs_desc)
#else
    el = mat%zmtx(jrow,icol)
#endif

  end subroutine dense_get_element_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine inserts an element into a given dense matrix, where the    !
  ! element can be local to any processor.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el   (input)   : The element to insert into the complex matrix           !
  !   mat  (inout)   : The complex dense matrix                                !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_put_element_complex(el,mat,jrow,icol)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: el        ! The element to insert
    type(DEM), intent(inout) :: mat           ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

#ifdef SCALAPACK
    ! ScaLAPACK subroutine
    external :: pzelset
#endif

    ! Check matrix is complex
    call utils_assert(mat%iscmplx,'Error in dense_put_element_complex: &
            &mat must be complex.')

    ! Put element directly in matrix, or call pzelset
#ifdef SCALAPACK
    call pzelset(mat%zmtx,jrow,icol,mat%blacs_desc,el)
#else
    mat%zmtx(jrow,icol) = el
#endif

  end subroutine dense_put_element_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! Retrieval selector for real/complex COEFs.                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element of the complex matrix retrieved             !
  !   mat  (input)   : The complex dense matrix                                !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by JM Escartin, August 2016.                                       !
  !============================================================================!

  subroutine dense_get_element_coef(el,mat,jrow,icol)

    use datatypes, only: COEF

    implicit none

    ! Arguments
    type(COEF), intent(inout) :: el           ! The element to insert
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! jme: No need to check that element and matrix are of the same type,
    ! since the check will be done inside dense_get_element_real/complex.

    if (el%iscmplx) then
       call dense_get_element_complex(el%z, mat, jrow, icol)
    else
       call dense_get_element_real(el%d, mat, jrow, icol)
    end if

  end subroutine dense_get_element_coef


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! Insertion selector for real/complex COEFs.                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el   (input)   : The element to insert into the complex matrix           !
  !   mat  (inout)   : The complex dense matrix                                !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by JM Escartin, August 2016.                                       !
  !============================================================================!

  subroutine dense_put_element_coef(el, mat, jrow, icol)

    use datatypes, only: COEF

    implicit none

    ! Arguments
    type(COEF), intent(in) :: el              ! The element to insert
    type(DEM), intent(inout) :: mat           ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! jme: No need to check that element and matrix are of the same type,
    ! since the check will be done inside dense_put_element_real/complex.

    if (el%iscmplx) then
       call dense_put_element_complex(el%z, mat, jrow, icol)
    else
       call dense_put_element_real(el%d, mat, jrow, icol)
    end if

  end subroutine dense_put_element_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine retrieves a column from a given dense matrix.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   col (output)   : The column of the real matrix retrieved                 !
  !   mat  (input)   : The real dense matrix                                   !
  !   icol (input)   : The column index of the column required                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_get_col_real(col,mat,icol)

#ifdef SCALAPACK
    use comms, only: comms_reduce
#endif
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: col(:)      ! The column to insert
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: jrow
#ifdef SCALAPACK
    complex(kind=DP) :: tmp
    external :: pzelget, pdelget ! ScaLAPACK subroutines
#endif

    ! Check size of col
    call utils_assert(size(col)==mat%nrows,'Error in dense_get_col_real: &
            &incorrect size of input column.')
    col(:) = 0.0_DP

    ! Get column directly from matrix, or call pdelget
    if (mat%iscmplx) then
       do jrow=1,mat%nrows
#ifdef SCALAPACK
          call pzelget(' ',' ',tmp,mat%zmtx,jrow,icol,mat%blacs_desc)
          col(jrow) = real(tmp,kind=DP)
#else
          col(jrow) = real(mat%zmtx(jrow,icol),kind=DP)
#endif
       end do
#ifdef SCALAPACK
       call comms_reduce('SUM',col)
#endif
    else
       do jrow=1,mat%nrows
#ifdef SCALAPACK
          call pdelget(' ',' ',col(jrow),mat%dmtx,jrow,icol,mat%blacs_desc)
#else
          col(jrow) = mat%dmtx(jrow,icol)
#endif
       end do
#ifdef SCALAPACK
       call comms_reduce('SUM',col)
#endif
    end if

  end subroutine dense_get_col_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine inserts a column into a given dense matrix.                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   col (output)   : The column to insert into the real matrix               !
  !   mat  (input)   : The real dense matrix                                   !
  !   icol (input)   : The column index of the column required                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_put_col_real(col,mat,icol)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: col(:)       ! The column to insert
    type(DEM), intent(inout) :: mat           ! The matrix
    integer, intent(in) :: icol               ! The index of the column

#ifdef SCALAPACK
    ! ScaLAPACK subroutines
    external :: pzelset, pdelset
#endif

    ! Local variables
    integer :: jrow

    ! Check size of col
    call utils_assert(size(col)==mat%nrows,'Error in dense_put_col_real: &
            &incorrect size of input column.')

    ! Put column directly in matrix, or call pdelset
    if (mat%iscmplx) then
       do jrow=1,mat%nrows
#ifdef SCALAPACK
          call pzelset(mat%zmtx,jrow,icol,mat%blacs_desc, &
               cmplx(col(jrow),0.0_DP,kind=DP))
#else
          mat%zmtx(jrow,icol) = cmplx(col(jrow),0.0_DP,kind=DP)
#endif
       end do
    else
       do jrow=1,mat%nrows
#ifdef SCALAPACK
          call pdelset(mat%dmtx,jrow,icol,mat%blacs_desc,col(jrow))
#else
          mat%dmtx(jrow,icol) = col(jrow)
#endif
       end do
    end if

  end subroutine dense_put_col_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine retrieves an column from a given dense matrix, where the   !
  ! column can be local to any processor.                                      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   col (output)   : The column of the complex matrix retrieved              !
  !   mat  (input)   : The complex dense matrix                                !
  !   icol (input)   : The column index of the column required                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_get_col_complex(col,mat,icol)

#ifdef SCALAPACK
    use comms, only: comms_reduce
#endif
    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: col(:)   ! The column to insert
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: icol               ! The index of the column

#ifdef SCALAPACK
    ! ScaLAPACK subroutine
    external :: pzelget
#endif

    ! Local variables
    integer :: jrow

    ! Check matrix is complex
    call utils_assert(mat%iscmplx,'Error in dense_get_col_complex: &
            &mat must be complex.')
    ! Check size of col
    call utils_assert(size(col)==mat%nrows,'Error in dense_get_col_complex: &
            &incorrect size of input column.')

    col(:) = cmplx(0.0_DP,0.0_DP,kind=DP)

    ! Get column directly from matrix, or call pdelget
    do jrow=1,mat%nrows
#ifdef SCALAPACK
       call pzelget(' ',' ',col(jrow),mat%zmtx,jrow,icol,mat%blacs_desc)
#else
       col(jrow) = mat%zmtx(jrow,icol)
#endif
    end do
#ifdef SCALAPACK
    call comms_reduce('SUM',col)
#endif

  end subroutine dense_get_col_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !============================================================================!
  ! This subroutine inserts an column into a given dense matrix.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   col  (input)   : The column to insert into the complex matrix            !
  !   mat  (inout)   : The complex dense matrix                                !
  !   jrow (input)   : The column index of the row required                    !
  !   icol (input)   : The column index of the column required                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 January 2011.                                 !
  !============================================================================!

  subroutine dense_put_col_complex(col,mat,icol)

    use utils, only: utils_assert

    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: col(:)    ! The column to insert
    type(DEM), intent(inout) :: mat           ! The matrix
    integer, intent(in) :: icol               ! The index of the column

#ifdef SCALAPACK
    ! ScaLAPACK subroutines
    external :: pzelset
#endif

    ! Local variables
    integer :: jrow

    ! Check matrix is complex
    call utils_assert(mat%iscmplx,'Error in dense_put_col_complex: &
            &mat must be complex.')
    ! Check size of col
    call utils_assert(size(col)==mat%nrows,'Error in dense_put_col_complex: &
            &incorrect size of input column.')

    ! Put column directly in matrix, or call pzelset
    do jrow=1,mat%nrows
#ifdef SCALAPACK
       call pzelset(mat%zmtx,jrow,icol,mat%blacs_desc,col(jrow))
#else
       mat%zmtx(jrow,icol) = col(jrow)
#endif
    end do

  end subroutine dense_put_col_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !============================================================================!
  ! This subroutine finds the solution to the normal AV=VL eigenproblem        !
  ! of amat = A and, and returns the real eigenvalues and a dense matrix       !
  ! storing the eigenvectors.                                                  !
  ! amat is overwritten by the eigensolvers                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   eigenvals (out)   : Real array storing the eigenvalues of amat.          !
  !   amat (input)      : Dense Symmetric or Hermitian matrix.                 !
  !   eigenvecs (out)   : Dense matrix storing eigenvectors defining solution. !
  !   n_eig (out)       : # eigenvalues required.                              !
  !   start_i (in,opt)  : Starting index of nxn submatrix                      !
  !   end_i (in,opt)    : Last index of nxn submatrix                          !
  !   rtype (in,opt)    : 0=(P)???EVX 1=(P)???EV                               !
  !----------------------------------------------------------------------------!
  ! Written by Louis Lee on April 4th 2012, based on code from the routine     !
  ! dense_eigensolve by Nicholas Hine.                                         !
  ! Works when start_i=1 and end_i <= amat%nrows but otherwise untested        !
  ! Complex solver untested. Currently only used in npa_mod and etrans_mod     !
  !============================================================================!

  subroutine dense_normal_eigensolve(n_eig,eigenvals,amat,eigenvecs, &
       start_i,end_i,rtype, arg_orfac, arg_abstol)

    use comms, only: comms_bcast, pub_on_root, pub_root_proc_id
    use constants, only: stdout
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
#ifdef SCALAPACK
#ifdef ELPA2
    use comms, only: pub_world_comm
    use elpa2, only: solve_evp_complex_2stage, solve_evp_real_2stage
    use elpa1, only: get_elpa_row_col_comms
#endif
#endif

    implicit none

    ! Arguments
    integer, intent(in) :: n_eig
    real(kind=DP), intent(out) :: eigenvals(:)
    type(DEM), intent(inout) :: eigenvecs
    type(DEM), intent(inout) :: amat
    integer, optional, intent(in) :: rtype

    ! lpl: Start/End submatrix index
    integer, optional, intent(in) :: start_i, end_i

    ! ars: ORFAC and ABSTOL values for ScaLAPACK
    real(kind=DP), optional, intent(in) :: arg_orfac ! ars
    real(kind=DP), optional, intent(in) :: arg_abstol ! ars

    ! Local Variables
#ifdef SCALAPACK
    integer :: nzfound
    real(kind=DP), external :: pdlamch
#endif
    real(kind=DP) :: abstol, orfac
    integer :: nfound
    integer :: ierr, info
    character :: jobz, range

    integer :: submat_dim, start_idx, end_idx, routine_type

    ! Work arrays and sizes
    integer :: lwork, liwork
#ifdef SCALAPACK
    integer :: lrwork
#endif
    integer, allocatable :: iwork(:), ifail(:)
    real(kind=DP), allocatable :: work(:)
    complex(kind=DP), allocatable :: zwork(:)
#ifdef SCALAPACK
    integer, allocatable :: iclustr(:)
    real(kind=DP), allocatable :: gap(:)
#ifdef ELPA2
    integer :: mpi_comm_rows
    integer :: mpi_comm_cols
#endif
#endif

    ! (Sca)LAPACK subroutines
#ifdef SCALAPACK
#ifdef ELPA2
#else
    external :: PZHEEVX, PZHEEV, PDSYEVX, PDSYEV
#endif
#else
    external :: ZHEEVX, ZHEEV, DSYEVX, DSYEV
#endif

    ! Start timer
    call timer_clock('dense_normal_eigensolve',1)

    if(pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: Entering dense_normal_eigensolve'

    ! rab207: set abstol here to avoid warning in memory inquiry
    abstol = 1.0e-9_DP

    ! lpl: Submatrix arguments check
    if(.not. present(start_i)) then
       start_idx = 1
    else
       start_idx = start_i
    end if
    if(.not. present(end_i)) then
       end_idx = amat%nrows
    else
       end_idx = end_i
    end if

    if (start_idx < 0) then
       call utils_abort(' ERROR: start_idx < 0 in dense_normal_eigensolve')
    else if (end_idx > amat%nrows) then
       call utils_abort(' ERROR: end_idx>amat%nrows in dense_normal_eigensolve')
    end if

    submat_dim = end_idx - start_idx + 1
    if (n_eig > 0 .and. submat_dim > n_eig) then
       call utils_abort(' ERROR: submat_dim > n_eig in dense_normal_eigensolve')
    else if (submat_dim < 0) then
       call utils_abort(' ERROR: submat_dim < 0 in dense_normal_eigensolve')
    end if

    ! Check arguments
    if (n_eig<=0) then
       jobz = 'N'
       range = 'A'
    else
       jobz = 'V'
       if ((n_eig>0).and.(n_eig<submat_dim)) then
          range = 'I'
       else
          range = 'A'
       end if
    end if

    if (present(rtype)) then
       routine_type = rtype
       if(routine_type < 0 .or. routine_type > 3) then
          call utils_abort(' ERROR: Invalid routine_type &
               &in dense_normal_eigensolve')
       end if
    else
       routine_type = 0
    end if

    ! lpl: Memory inquiry
    if (routine_type == 2 .or. routine_type == 3) then
       if(size(eigenvals) < 2) then
          call utils_abort(' ERROR: size(eigenvals) < 2 in &
               &dense_normal_eigensolve  during inquiry run')
       end if
       call internal_allocate_work
       call internal_deallocate_work

       ! 1st index stores real size
       eigenvals(lbound(eigenvals)) = real(lwork,kind=DP)
       ! 2nd index stores int size
       eigenvals(lbound(eigenvals) + 1) = real(liwork,kind=DP)
    ! lpl: Else normal run
    else
       ! Allocate work arrays
       call internal_allocate_work

#ifdef SCALAPACK
#ifdef ELPA2
    ! Get elpa communicators
    call get_elpa_row_col_comms(pub_world_comm, blacs_myrow, blacs_mycol, &
         & mpi_comm_rows, mpi_comm_cols)
    if (amat%iscmplx) then
       ! Solve eigenproblem
       call solve_evp_complex_2stage(submat_dim, n_eig, amat%zmtx, amat%blacs_ld, eigenvals, &
            &eigenvecs%zmtx(1,1), eigenvecs%blacs_ld, amat%blacs_nb, mpi_comm_rows, mpi_comm_cols, pub_world_comm)
    else
       ! Solve eigenproblem
       call solve_evp_real_2stage(submat_dim, n_eig, amat%dmtx, amat%blacs_ld, eigenvals, &
            &eigenvecs%dmtx(1,1), eigenvecs%blacs_ld, amat%blacs_nb, mpi_comm_rows, mpi_comm_cols, pub_world_comm)
    end if
    ! Clean up
    call mpi_comm_free(mpi_comm_rows,info)
    call mpi_comm_free(mpi_comm_cols,info)
#else
       ! lpl: call ScaLAPACK routine to solve the normal eigenproblem
       ! lpl: Bisection (-X) version copied from dense_eigensolve

       ! ars: set orfac
       if (present(arg_orfac)) then
          if (arg_orfac.le.tiny(1.0_DP)) then
             orfac = 0.0_DP
          else
             orfac = arg_orfac
          end if
       else
          orfac = 1e-4_DP
       end if
       ! ars: set abstol
       if (present(arg_abstol)) then
          if (arg_abstol.le.tiny(1.0_DP)) then
             abstol = 2.0_DP*PDLAMCH(blacs_context,'S')
          else
             abstol = arg_abstol
          end if
       else
          abstol = 1e-9_DP
       end if

       ! lpl: Works correctly when start_idx = 1. Don't know what will
       !      happen otherwise, presumably it will work too
       if (start_idx /= 1) then
          if (pub_on_root) write(stdout,'(a)') &
               ' WARNING: start_idx /= 1. Not sure what &
               &will happen to PDSYEV(X) - use at your own risk!'
       end if

       if (amat%iscmplx) then
          if (routine_type == 0) then

             if(pub_debug_on_root) write(stdout,'(a)') ' DEBUG: PZHEEVX'

             call PZHEEVX(jobz, range, 'L', submat_dim, amat%zmtx(1,1), &
                  start_idx, start_idx, amat%blacs_desc, &
                  0.0_DP, 0.0_DP, 1, n_eig, abstol, nfound, nzfound, &
                  eigenvals(start_idx:end_idx), orfac, eigenvecs%zmtx(1,1), &
                  start_idx, start_idx, eigenvecs%blacs_desc, &
                  zwork, lwork, work, lrwork, iwork, liwork, ifail, &
                  iclustr, gap, info)
          else

             if(pub_debug_on_root) write(stdout,'(a)') ' DEBUG: PZHEEV'

             call PZHEEV(jobz, 'L', submat_dim, amat%zmtx(1,1), &
                  start_idx, start_idx, amat%blacs_desc, &
                  eigenvals(start_idx:end_idx), eigenvecs%zmtx(1,1), &
                  start_idx, start_idx, eigenvecs%blacs_desc, &
                  zwork, lwork, work, lrwork, info)
          end if
       else
          if(routine_type == 0) then

             if(pub_debug_on_root) write(stdout,'(a)') ' DEBUG: PDSYEVX'

             call PDSYEVX(jobz , range, 'L', submat_dim, amat%dmtx(1,1), &
                  start_idx, start_idx, amat%blacs_desc, &
                  0.0_DP, 0.0_DP, 1, n_eig, abstol, nfound, nzfound, &
                  eigenvals(start_idx:end_idx), orfac, eigenvecs%dmtx(1,1), &
                  start_idx, start_idx, eigenvecs%blacs_desc, &
                  work, lwork, iwork, liwork, ifail, iclustr, gap, info)
          else

             if(pub_debug_on_root) write(stdout,'(a)') ' DEBUG: PDSYEV'

             call PDSYEV(jobz, 'L', submat_dim, amat%dmtx(1,1), &
                  start_idx, start_idx, amat%blacs_desc, &
                  eigenvals(start_idx:end_idx), eigenvecs%dmtx(1,1), &
                  start_idx, start_idx, eigenvecs%blacs_desc, &
                  work, lwork, info)
          end if
       end if
#endif
! lpl: Else if not ScaLAPACK
#else
       if (present(arg_orfac).and.present(arg_abstol)) then
           orfac = arg_orfac
           abstol = arg_abstol
       end if
       ! ndmh: call LAPACK routine to solve the generalised eigenproblem
       if (amat%iscmplx) then
          if (routine_type == 0) then

             if(pub_debug_on_root) write(stdout,'(a)') ' DEBUG: DSYEVX'

             call ZHEEVX(jobz,range,'L',submat_dim, &
                  amat%zmtx(start_idx:end_idx,start_idx:end_idx), &
                  submat_dim, 0.0_DP, 0.0_DP, 1, n_eig, -1.0_DP, nfound, &
                  eigenvals(start_idx:end_idx), &
                  eigenvecs%zmtx(start_idx:end_idx,start_idx:end_idx), &
                  submat_dim, zwork, lwork, work, iwork, ifail, info)
          else

             if(pub_debug_on_root) write(stdout,'(a)') ' DEBUG: DSYEV'

             call ZHEEV(jobz, 'L', submat_dim, &
                  amat%zmtx(start_idx:end_idx,start_idx:end_idx), submat_dim, &
                  eigenvals(start_idx:end_idx), zwork, lwork, work, info)
          end if

          ! ndmh: degenerate eigenvectors may have come out in different linear
          ! ndmh: combinations on different procs, so synchronise to
          ! ndmh: eigenvectors found on root proc
          call comms_bcast(pub_root_proc_id,eigenvecs%zmtx, &
               eigenvecs%nrows*eigenvecs%mcols)
       else
          if (routine_type == 0) then
             call DSYEVX(jobz,range,'L',submat_dim, &
                  amat%dmtx(start_idx:end_idx,start_idx:end_idx), &
                  submat_dim, 0.0_DP, 0.0_DP, 1, n_eig, &
                  -1.0_DP, nfound, eigenvals(start_idx:end_idx), &
                  eigenvecs%dmtx(start_idx:end_idx,start_idx:end_idx), &
                  submat_dim,work,lwork,iwork,ifail,info)
          else
             call DSYEV(jobz, 'L', submat_dim, &
                  amat%dmtx(start_idx:end_idx,start_idx:end_idx), submat_dim, &
                  eigenvals(start_idx:end_idx), work, lwork, info)
             call dense_copy(eigenvecs,amat)
          end if

          ! ndmh: degenerate eigenvectors may have come out in different linear
          ! ndmh: combinations on different procs, so synchronise to
          ! ndmh:  eigenvectors found on root proc
          call comms_bcast(pub_root_proc_id,eigenvecs%dmtx, &
               eigenvecs%nrows*eigenvecs%mcols)
       end if
! lpl: END #ifdef SCALAPACK
#endif

       ! ndmh: check for errors
       if (info/=0) then
          if (amat%iscmplx) then
             if (pub_on_root) then
                write(stdout,'(a,i5)') '(P)ZHEEVX in subroutine &
                     &dense_normal_eigensolve returned info=',info
                write(stdout,*) 'ifail=',ifail(1:info)
             end if
          else
             if (pub_on_root) then
                write(stdout,'(a,i5)') '(P)DSYEVX in subroutine &
                     &dense_normal_eigensolve returned info=',info
                write(stdout,*) 'ifail=',ifail(1:info)
             end if
          end if
          ! ndmh: continue even if some eigenvectors did not converge
          if ((info<1).or.(info>amat%nrows)) then
             call utils_abort('Error in dense_eigensolve(). &
                  &The error message should be in your output file.')
          end if
       end if

       ! Deallocate work arrays
       call internal_deallocate_work
    end if ! END if (routine_type == 1 .or. routine_type == 3)

    if(pub_debug_on_root) write(stdout,'(a)') &
         ' DEBUG: Leaving dense_normal_eigensolve'

    ! Stop timer
    call timer_clock('dense_normal_eigensolve',2)

    return

contains

    ! Allocate diagonalisation work arrays
    subroutine internal_allocate_work

#ifdef SCALAPACK
      use utils, only: utils_safe_nint
      implicit none

      ! ScaLAPACK subroutines
      external :: PZHEEVX, PDSYEVX, PZHEEV, PDSYEV

      ! Local Variables
      integer :: liclustr

      ! lpl: For mem inquiry
      lwork = 0
      lrwork = 0
      liwork = 0
      liclustr = 0

      ! Calculate optimal sizes of work arrays
      if(routine_type == 0 .or. routine_type == 2) then

         liclustr = 2*blacs_nprow*blacs_npcol

         if(routine_type == 0) then
            allocate(iclustr(liclustr),stat=ierr)
            call utils_alloc_check('dense_normal_eigensolve','iclustr',ierr)
            allocate(gap(liclustr),stat=ierr)
            call utils_alloc_check('dense_normal_eigensolve','gap',ierr)
            allocate(ifail(amat%nrows),stat=ierr)
            call utils_alloc_check('dense_normal_eigensolve','ifail',ierr)
         end if

         allocate(work(1),stat=ierr)
         call utils_alloc_check('dense_normal_eigensolve','work',ierr)
         allocate(iwork(1),stat=ierr)
         call utils_alloc_check('dense_normal_eigensolve','iwork',ierr)

         if (amat%iscmplx) then
            allocate(zwork(1),stat=ierr)
            call utils_alloc_check('dense_normal_eigensolve','zwork',ierr)

            call PZHEEVX(jobz, range, 'L', submat_dim, amat%zmtx(1,1), &
                 start_idx, start_idx, amat%blacs_desc, &
                 0.0_DP, 0.0_DP, 1, n_eig, abstol, nfound, nzfound, &
                 eigenvals(start_idx:end_idx), orfac, eigenvecs%zmtx(1,1), &
                 start_idx, start_idx, eigenvecs%blacs_desc, &
                 zwork, -1, work, -1, iwork, -1, ifail, iclustr, gap, info)
            lwork = utils_safe_nint(real(zwork(1)))
            lrwork = utils_safe_nint(work(1))
            liwork = iwork(1)
            deallocate(zwork,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','zwork',ierr)
            deallocate(iwork,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','iwork',ierr)
            deallocate(work,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','work',ierr)
         else
            call PDSYEVX(jobz , range, 'L', submat_dim, amat%dmtx(1,1), &
                 start_idx, start_idx, amat%blacs_desc, &
                 0.0_DP, 0.0_DP, 1, n_eig, abstol, nfound, nzfound, &
                 eigenvals(start_idx:end_idx), orfac, eigenvecs%dmtx(1,1), &
                 start_idx, start_idx, eigenvecs%blacs_desc, &
                 work, -1, iwork, -1, ifail, iclustr, gap, info)
            lwork = utils_safe_nint(work(1))
            liwork = iwork(1)
            deallocate(iwork,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','iwork',ierr)
            deallocate(work,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','work',ierr)
         end if

         ! lpl: Allocate arrays
         if (routine_type == 0) then
            if (amat%iscmplx) then
               allocate(work(lrwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','work',ierr)
               allocate(iwork(liwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','iwork',ierr)
               allocate(zwork(lwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','zwork',ierr)
            else
               allocate(work(lwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','work',ierr)
               allocate(iwork(liwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','iwork',ierr)
            end if
         else
            lwork = lwork + lrwork + liclustr
            liwork = liwork + liclustr + amat%nrows
         end if

      else if (routine_type == 1 .or. routine_type == 3) then

         allocate(work(1),stat=ierr)
         call utils_alloc_check('dense_normal_eigensolve','work',ierr)

         if (amat%iscmplx) then
            allocate(zwork(1),stat=ierr)
            call utils_alloc_check('dense_normal_eigensolve','zwork',ierr)

            call PZHEEV(jobz, 'L', submat_dim, amat%zmtx(1,1), &
                 start_idx, start_idx, amat%blacs_desc, &
                 eigenvals(start_idx:end_idx), eigenvecs%zmtx(1,1), &
                 start_idx, start_idx, eigenvecs%blacs_desc, &
                 zwork, -1, work, -1, info)
            lwork = utils_safe_nint(real(zwork(1)))
            lrwork = utils_safe_nint(work(1))
            deallocate(zwork,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','zwork',ierr)
            deallocate(work,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','work',ierr)
         else
            call PDSYEV(jobz, 'L', submat_dim, amat%dmtx(1,1), &
                 start_idx, start_idx, amat%blacs_desc, &
                 eigenvals(start_idx:end_idx),eigenvecs%dmtx(1,1), &
                 start_idx, start_idx, eigenvecs%blacs_desc, work, -1, info)
            lwork = utils_safe_nint(work(1))
            deallocate(work,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','work',ierr)
         end if

         ! lpl: Allocate arrays
         if (routine_type == 1) then
            if (amat%iscmplx) then
               allocate(work(lrwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','work',ierr)
               allocate(zwork(lwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','zwork',ierr)
            else
               allocate(work(lwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','work',ierr)
            end if
         else
            lwork = lwork + lrwork
         end if
      else
         call utils_abort('Invalid routine_type in dense_normal_eigensolve')
      end if

#else
      implicit none

      lwork = 0
      liwork = 0

      ! Allocate work arrays
      if (routine_type == 0 .or. routine_type == 2) then
         lwork = 8*amat%nrows
         liwork = 5*amat%nrows
         if(routine_type == 0) then
            if (amat%iscmplx) then
               allocate(zwork(lwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','zwork',ierr)
            end if
            allocate(work(lwork),stat=ierr)
            call utils_alloc_check('dense_normal_eigensolve','work',ierr)
            allocate(iwork(liwork),stat=ierr)
            call utils_alloc_check('dense_normal_eigensolve','iwork',ierr)
            allocate(ifail(amat%nrows),stat=ierr)
            call utils_alloc_check('dense_normal_eigensolve','ifail',ierr)
         else
            liwork = liwork + amat%nrows
         end if
      else if(routine_type == 1 .or. routine_type == 3) then
         lwork = 5*submat_dim
         if(routine_type == 1) then
            allocate(work(lwork),stat=ierr)
            call utils_alloc_check('dense_normal_eigensolve','work',ierr)
            if (amat%iscmplx) then
               allocate(zwork(lwork),stat=ierr)
               call utils_alloc_check('dense_normal_eigensolve','zwork',ierr)
            end if
         else
            if (amat%iscmplx) lwork = 2*lwork
         end if
      else
         call utils_abort('Invalid routine_type in dense_normal_eigensolve')
      end if

#endif

    end subroutine internal_allocate_work

    ! Deallocate diagonalisation work arrays
    subroutine internal_deallocate_work

      implicit none

#ifdef SCALAPACK

      if(routine_type == 0) then
         if (amat%iscmplx) then
            deallocate(zwork,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','zwork',ierr)
         end if
         deallocate(iwork,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','iwork',ierr)
         deallocate(work,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','work',ierr)
         deallocate(ifail,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','ifail',ierr)
         deallocate(gap,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','gap',ierr)
         deallocate(iclustr,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','iclustr',ierr)

      else if(routine_type == 1) then
         if (amat%iscmplx) then
            deallocate(zwork,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','zwork',ierr)
         end if
         deallocate(work,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','work',ierr)
      end if

#else

      if(routine_type == 0) then
         deallocate(ifail,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','ifail',ierr)
         deallocate(iwork,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','iwork',ierr)
         deallocate(work,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','work',ierr)
         if (amat%iscmplx) then
            deallocate(zwork,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','zwork',ierr)
         end if
      else if(routine_type == 1) then
         deallocate(work,stat=ierr)
         call utils_dealloc_check('dense_normal_eigensolve','work',ierr)
         if (amat%iscmplx) then
            deallocate(zwork,stat=ierr)
            call utils_dealloc_check('dense_normal_eigensolve','zwork',ierr)
         end if
      end if

#endif

    end subroutine internal_deallocate_work

  end subroutine dense_normal_eigensolve


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !============================================================================!
  ! This subroutine retrieves an element from a given dense matrix, where the  !
  ! element can be local to any processor, and sends it ONLY to recv_proc      !
  ! Only relevant for ScaLAPACK - otherwise it's simply dense_get_element      !
  ! Could be merged with dense_get_element with optional recv_proc argument?   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element of the real matrix retrieved                !
  !   mat  (input)   : The real dense matrix                                   !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !   recv_proc (input) : Proc to receive element                              !
  !----------------------------------------------------------------------------!
  ! lpl: Plagiarized from dense_get_element 26-04-2012                         !
  !============================================================================!

  subroutine dense_point2point_element_real(el,mat,jrow,icol,recv_proc)

    use comms, only: comms_barrier, pub_my_proc_id
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: el          ! The element to obtain
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column
    integer, intent(in) :: recv_proc          ! The proc to receive el

#ifdef SCALAPACK
    ! lpl: Local variables
    integer :: lrindx, lcindx, rsrc, csrc, rdest, cdest
    real(kind=DP) :: el_aux(1,1)

    ! BLACS & ScaLAPACK subroutines
    external :: blacs_pcoord, infog2l, dgesd2d, dgerv2d
#endif

    ! Check matrix is real
    call utils_assert(.not. mat%iscmplx, &
          'Error in dense_point2point_element_real: mat must be real.')

    ! Get element directly from matrix, or call pdelget
#ifdef SCALAPACK
    ! lpl: Get grid coordinate for proc containing el
    call infog2l(jrow,icol,mat%blacs_desc,blacs_nprow,blacs_npcol, &
         blacs_myrow,blacs_mycol,lrindx,lcindx,rsrc,csrc)
    ! lpl: Get grid coordinate for proc receiving el
    call blacs_pcoord(blacs_context,recv_proc,rdest,cdest)
    ! lpl: Send el to recv_proc
    if((rsrc == rdest .and. csrc == cdest) &
         .and. (pub_my_proc_id == recv_proc)) then
       el = mat%dmtx(lrindx,lcindx)
    else
       if(blacs_myrow == rsrc .and. blacs_mycol == csrc) then
          call dgesd2d(blacs_context,1,1, &
               mat%dmtx(lrindx,lcindx),1,rdest,cdest)
       else if(blacs_myrow == rdest .and. blacs_mycol == cdest) then
          call dgerv2d(blacs_context,1,1,el_aux,1,rsrc,csrc)
          el = el_aux(1,1)
       end if
    end if
    call comms_barrier
#else
    if ((recv_proc==pub_my_proc_id).and.(icol==jrow)) el = 0.0_DP
    el = mat%dmtx(jrow,icol)
#endif

  end subroutine dense_point2point_element_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !============================================================================!
  ! This subroutine retrieves an element from a given dense matrix, where the  !
  ! element can be local to any processor, and sends it ONLY to recv_proc      !
  ! Only relevant for ScaLAPACK - otherwise it's simply dense_get_element      !
  ! Could be merged with dense_get_element with optional recv_proc argument?   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element of the real matrix retrieved                !
  !   mat  (input)   : The real dense matrix                                   !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !   recv_proc (input) : Proc to receive element                              !
  !----------------------------------------------------------------------------!
  ! lpl: Plagiarized from dense_get_element 26-04-2012                         !
  !============================================================================!

  subroutine dense_point2point_element_complex(el,mat,jrow,icol,recv_proc)

    use comms, only: comms_barrier, pub_my_proc_id
    use constants, only: cmplx_0
    use utils, only: utils_assert
    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: el          ! The element to obtain
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column
    integer, intent(in) :: recv_proc          ! The proc to receive el

#ifdef SCALAPACK
    ! lpl: Local variables
    integer :: lrindx, lcindx, rsrc, csrc, rdest, cdest
    complex(kind=DP) :: el_aux(1,1)

    ! BLACS & ScaLAPACK subroutines
    external :: blacs_pcoord, infog2l, zgesd2d, zgerv2d
#endif

    ! Check matrix is real
    call utils_assert(mat%iscmplx, &
         'Error in dense_point2point_element_complex: mat must be complex.')

    ! Get element directly from matrix, or call pdelget
#ifdef SCALAPACK
    ! lpl: Get grid coordinate for proc containing el
    call infog2l(jrow,icol,mat%blacs_desc,blacs_nprow,blacs_npcol, &
         blacs_myrow,blacs_mycol,lrindx,lcindx,rsrc,csrc)
    ! lpl: Get grid coordinate for proc receiving el
    call blacs_pcoord(blacs_context,recv_proc,rdest,cdest)
    ! lpl: Send el to recv_proc
    if((rsrc == rdest .and. csrc == cdest) &
         .and. (pub_my_proc_id == recv_proc)) then
       el = mat%zmtx(lrindx,lcindx)
    else
       if(blacs_myrow == rsrc .and. blacs_mycol == csrc) then
          call zgesd2d(blacs_context,1,1, &
               mat%zmtx(lrindx,lcindx),1,rdest,cdest)
       else if(blacs_myrow == rdest .and. blacs_mycol == cdest) then
          call zgerv2d(blacs_context,1,1,el_aux,1,rsrc,csrc)
          el = el_aux(1,1)
       end if
    end if
    call comms_barrier
#else
    if ((recv_proc==pub_my_proc_id).and.(icol==jrow)) el = cmplx_0
    el = mat%zmtx(jrow,icol)
#endif

  end subroutine dense_point2point_element_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine copies a matrix, but possibly across blacs contexts &      !
  ! different distributions.                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat_src  (input)                                                         !
  !   mat_dest (output)                                                        !
  !----------------------------------------------------------------------------!
  !   Written by Jolyon Aarons 03/08/14.                                       !
  !============================================================================!

  subroutine dense_redist(mat_src, mat_dest)
    use comms, only: comms_barrier, comms_bcast
#ifndef SCALAPACK
    use comms, only: pub_root_proc_id
#endif
    use utils, only: utils_assert
    implicit none
    type(DEM), intent(in)    :: mat_src
    type(DEM), intent(inout) :: mat_dest

#ifdef SCALAPACK
    ! ScaLAPACK subroutines
    external :: pdgemr2d, pzgemr2d
#endif

    if(.not.mat_src%iscmplx) then
       call utils_assert(.not.mat_dest%iscmplx, &
            'Error in dense_redist: input matrices must be both real/complex.')
    else
       call utils_assert(mat_dest%iscmplx, &
            'Error in dense_redist: input matrices must be both real/complex.')
    end if

#ifdef SCALAPACK
    if(.not.mat_src%iscmplx) then
       call pdgemr2d(mat_src%nrows,mat_src%mcols,mat_src%dmtx(1,1),1,1,mat_src%blacs_desc,mat_dest%dmtx(1,1),1,1, &
            & mat_dest%blacs_desc,blacs_context)
    else
       call pzgemr2d(mat_src%nrows,mat_src%mcols,mat_src%zmtx(1,1),1,1,mat_src%blacs_desc,mat_dest%zmtx(1,1),1,1, &
            & mat_dest%blacs_desc,blacs_context)
    end if
#else
    if(.not.mat_src%iscmplx) then
       mat_dest%dmtx(:,:) = mat_src%dmtx(:,:)
       call comms_bcast(pub_root_proc_id, mat_dest%dmtx(:,:))
    else
       mat_dest%zmtx(:,:) = mat_src%zmtx(:,:)
       call comms_bcast(pub_root_proc_id, mat_dest%zmtx(:,:))
    end if
#endif

  end subroutine dense_redist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine converts a SPAM3_EMBED matrix to a dense matrix.           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dmat      (input) : The destination dense matrix                         !
  !   smat      (input) : The source embedded sparse matrix                    !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, 29 May 2018.                                   !
  ! Modifiied to use sparse_embed routines by Robert Charlton, 21/08/2018.     !
  !============================================================================!

  subroutine dense_convert_sparseembedtodense(dmat_out,smat_in)

    use sparse_embed, only: SPAM3_EMBED,sparse_embed_num_cols,sparse_embed_num_rows
#ifdef SCALAPACK
    use sparse, only: sparse_spam3toblacs,sparse_num_cols,sparse_num_rows
#else
    use sparse, only: sparse_convert,sparse_num_cols,sparse_num_rows
#endif
    use utils, only: utils_abort,utils_alloc_check,utils_dealloc_check, &
         utils_assert
    use constants, only: DP

    implicit none

    ! Arguments
    type(DEM),intent(inout) :: dmat_out
    type(SPAM3_EMBED),intent(in) :: smat_in

    type(DEM) :: dense_buffer
    integer :: past_ncol,past_mrow,sub_col,sub_row
    integer :: col,row,ierr
#ifdef SCALAPACK
    integer :: blacs_ld,blacs_ncol,blacs_desc(50)
#endif
    real(DP),allocatable :: dtmp(:,:)
    complex(DP),allocatable :: ztmp(:,:)
    logical :: loc_iscmplx

    loc_iscmplx = smat_in%iscmplx

    if (smat_in%iscmplx.and.(.not.dmat_out%iscmplx)) &
         call utils_abort('Error in &
         &dense_convert_sparseembedtodense: source matrix is complex but &
         &dest matrix is real')

#ifdef SCALAPACK
    ! rc2013: EMBED_FIX!
    call utils_assert(smat_in%mrows == 1 .and. smat_in%ncols == 1, 'Error in &
         &dense_convert_sparseembedtodense: embedding diagonalisation not set &
         &up yet with ScaLAPACK. Consider running again without ScaLAPACK or &
         &placing all atoms in one line in %block species_ngwf_regions.')
    call dense_convert_sparsetodense(dmat_out, smat_in%p)
#else
    ! rc2013: get number of rows and columns directly from sparse_embed routines
    call dense_create(dense_buffer,int(sparse_embed_num_rows(smat_in)), &
         int(sparse_embed_num_cols(smat_in)),iscmplx=loc_iscmplx)

    ! Loop over regions
    past_ncol=0
    do col=1,smat_in%ncols
       past_mrow=0
       do row=1,smat_in%mrows

          sub_col=sparse_num_cols(smat_in%m(row,col))
          sub_row=sparse_num_rows(smat_in%m(row,col))
          if (loc_iscmplx) then
             allocate(ztmp(sub_row,sub_col),stat=ierr)
             call utils_alloc_check('dense_convert_sparseembedtodense',&
                  'ztmp',ierr)
          else
             allocate(dtmp(sub_row,sub_col),stat=ierr)
             call utils_alloc_check('dense_convert_sparseembedtodense',&
                  'dtmp',ierr)
          end if


          ! Fill dense matrix values from SPAM3 matrix
          if (loc_iscmplx) then
             call sparse_convert(ztmp,smat_in%m(row,col))
          else
             call sparse_convert(dtmp,smat_in%m(row,col))
          end if

          if (loc_iscmplx) then
             dense_buffer%zmtx(past_mrow+1:past_mrow+sub_row,&
                  past_ncol+1:past_ncol+sub_col)=ztmp
             deallocate(ztmp,stat=ierr)
             call utils_dealloc_check('dense_convert_sparseembedtodense',&
                  'ztmp',ierr)
          else
             dense_buffer%dmtx(past_mrow+1:past_mrow+sub_row,&
                  past_ncol+1:past_ncol+sub_col)=dtmp
             deallocate(dtmp,stat=ierr)
             call utils_dealloc_check('dense_convert_sparseembedtodense',&
                  'dtmp',ierr)
          end if

          past_mrow=past_mrow+sparse_num_rows(smat_in%m(row,col))
       end do
       past_ncol=past_ncol+sparse_num_cols(smat_in%m(1,col))
    end do

    call dense_copy(dmat_out,dense_buffer)

    call dense_destroy(dense_buffer)
#endif

  end subroutine dense_convert_sparseembedtodense

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine converts a dense matrix to a SPAM3_EMBED matrix.           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dmat      (input) : The destination embedded sparse matrix               !
  !   smat      (input) : The source dense matrix                              !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, 29 May 2018.                                   !
  !============================================================================!

  subroutine dense_convert_densetosparseembed(smat_out,dmat_in)

    use sparse_embed, only: SPAM3_EMBED
#ifdef SCALAPACK
    use sparse, only: SPAM3,sparse_blacstospam3,sparse_copy,sparse_create,&
         sparse_destroy,sparse_num_cols,sparse_num_rows
#else
    use sparse, only: SPAM3,sparse_convert,sparse_copy,sparse_create,&
         sparse_destroy,sparse_num_cols,sparse_num_rows
#endif
    use utils, only: utils_abort,utils_alloc_check,utils_dealloc_check,utils_assert
    use constants, only: DP

    implicit none

    ! Arguments
    type(SPAM3_EMBED),intent(inout) :: smat_out
    type(DEM),intent(in) :: dmat_in

    type(SPAM3) :: sparse_buffer
    integer :: past_ncol(smat_out%ncols+1),past_mrow(smat_out%mrows+1)
    integer :: col,row,ierr
    real(DP),allocatable :: dtmp(:,:)
    complex(DP),allocatable :: ztmp(:,:)
    logical :: loc_iscmplx
#ifdef SCALAPACK
    integer :: blacs_ld,blacs_ncol,blacs_desc(50)
#endif

#ifdef SCALAPACK
    ! rc2013: EMBED_FIX!
    call utils_assert(smat_out%mrows == 1 .and. smat_out%ncols == 1, 'Error in &
         &dense_convert_sparseembedtodense: embedding diagonalisation not set &
         &up yet with ScaLAPACK. Consider running again without ScaLAPACK or &
         &placing all atoms in one line in %block species_ngwf_regions.')
    call dense_convert_densetosparse(smat_out%p, dmat_in)
#else
    ! Find where each embedding block starts
    past_ncol(1)=0
    do col=1,smat_out%ncols
       past_ncol(col+1)=past_ncol(col)+sparse_num_cols(smat_out%m(1,col))
    end do
    past_mrow(1)=0
    do row=1,smat_out%mrows
       past_mrow(row+1)=past_mrow(row)+sparse_num_rows(smat_out%m(row,1))
    end do

    loc_iscmplx=dmat_in%iscmplx

    ! Create buffer sparse matrix
    do col=1,smat_out%ncols
       do row=1,smat_out%mrows
          call sparse_create(sparse_buffer,smat_out%m(row,col))

          if (loc_iscmplx) then
             allocate(ztmp(past_mrow(row+1)-past_mrow(row),&
                  past_ncol(col+1)-past_ncol(col)),stat=ierr)
             call utils_alloc_check('dense_convert_densetosparseembed',&
                  'ztmp',ierr)
             ztmp=dmat_in%zmtx((past_mrow(row)+1):past_mrow(row+1),&
                  (past_ncol(col)+1):past_ncol(col+1))
             call sparse_convert(sparse_buffer,ztmp)
             deallocate(ztmp,stat=ierr)
             call utils_dealloc_check('dense_convert_densetosparseembed',&
                  'ztmp',ierr)
          else
             allocate(dtmp(past_mrow(row+1)-past_mrow(row),&
                  past_ncol(col+1)-past_ncol(col)),stat=ierr)
             call utils_alloc_check('dense_convert_densetosparseembed',&
                  'dtmp',ierr)
             dtmp=dmat_in%dmtx((past_mrow(row)+1):past_mrow(row+1),&
                  (past_ncol(col)+1):past_ncol(col+1))
             call sparse_convert(sparse_buffer,dtmp)
             deallocate(dtmp,stat=ierr)
             call utils_dealloc_check('dense_convert_densetosparseembed',&
                  'dtmp',ierr)
          end if
          ! Copy sparse buffer into output sparse embed matrix and destroy buffer
          call sparse_copy(smat_out%m(row,col),sparse_buffer)
          call sparse_destroy(sparse_buffer)
       end do
    end do
#endif

  end subroutine dense_convert_densetosparseembed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine converts a SPAM3_EMBED matrix to a SPAM3 matrix.           !
  ! I realise that dense_mod is an odd place for this, but it does the         !
  ! conversion via a DEM matrix, so it needs to be in here...                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   smat_out     (output): The destination sparse matrix                     !
  !   smat_in      (input) : The source embedded sparse matrix                 !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, 29 May 2018.                                   !
  !============================================================================!

  subroutine dense_convert_sparseembedtosparse(smat_out,smat_in)

    use sparse_embed, only: SPAM3_EMBED
    use sparse, only: SPAM3, sparse_num_cols, sparse_num_rows

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: smat_out
    type(SPAM3_EMBED),intent(in) :: smat_in

    type(DEM) :: dense_buffer
    integer :: ncols,mrows
    integer :: col,row
    logical :: loc_iscmplx

    ! Find total size of dense matrix
    ncols=0
    mrows=0
    loc_iscmplx=.false.
    do col=1,smat_in%ncols
       do row=1,smat_in%mrows
          ncols=ncols+sparse_num_cols(smat_in%m(row,col))
          mrows=mrows+sparse_num_rows(smat_in%m(row,col))
          if (smat_in%m(row,col)%iscmplx) loc_iscmplx=.true.
       end do
    end do

    call dense_create(dense_buffer,mrows,ncols,iscmplx=loc_iscmplx)

    call dense_convert_sparseembedtodense(dense_buffer,smat_in)

    call dense_convert_densetosparse(smat_out,dense_buffer)

    call dense_destroy(dense_buffer)

  end subroutine dense_convert_sparseembedtosparse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine converts a SPAM3 matrix to a SPAM3_EMBED matrix.           !
  ! I realise that dense_mod is an odd place for this, but it does the         !
  ! conversion via a DEM matrix, so it needs to be in here...                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   smat_out     (output): The destination sparse embed matrix               !
  !   smat_in      (input) : The source sparse matrix                          !
  !----------------------------------------------------------------------------!
  ! Written by Joseph Prentice, 29 May 2018.                                   !
  !============================================================================!

  subroutine dense_convert_sparsetosparseembed(smat_out,smat_in)

    use sparse_embed, only: SPAM3_EMBED
    use sparse, only: SPAM3, sparse_num_cols, sparse_num_rows

    implicit none

    ! Arguments
    type(SPAM3_EMBED),intent(inout) :: smat_out
    type(SPAM3),intent(in) :: smat_in

    type(DEM) :: dense_buffer
    integer :: ncols,mrows
    logical :: loc_iscmplx

    ! Find total size of dense matrix
    ncols=int(sparse_num_cols(smat_in))
    mrows=int(sparse_num_rows(smat_in))
    loc_iscmplx=smat_in%iscmplx

    call dense_create(dense_buffer,mrows,ncols,iscmplx=loc_iscmplx)

    call dense_convert_sparsetodense(dense_buffer,smat_in)

    call dense_convert_densetosparseembed(smat_out,dense_buffer)

    call dense_destroy(dense_buffer)

  end subroutine dense_convert_sparsetosparseembed

end module dense
