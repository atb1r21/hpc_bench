! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Peter D. Haynes and Nicholas Hine
!
!      TCM Group, Cavendish laboratory
!      Madingley Road
!      Cambridge CB3 0HE
!      UK
!
!   Subsequent additions and modifications by: James C Womack,
!   Andrea Greco, and Jose M Escartin.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module kinetic

  implicit none

  private

  public :: kinetic_apply_on_box
  public :: kinetic_apply_to_func_batch
  public :: kinetic_apply_to_funcs_on_grid
  public :: kinetic_grad_on_box
  public :: kinetic_grad_on_cart_box ! JCW
  public :: kinetic_grad_to_func_batch
  public :: kinetic_grad_to_funcs_on_grid ! JCW
#ifdef GPU_PGI_KINETIC
  public :: kinetic_gpu_app2_func_batch
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine kinetic_apply_to_func_batch(kinetic_fftbox_batch, &  ! output
       funcs_on_grid, fbasis, cell, fftbox, ket_start_in_box, &   ! input
       batch_size, local_start, local_end)                        ! input

    !======================================================================!
    ! This subroutine applies the kinetic energy operator to a batch       !
    ! of functions (eg NGWFs) in fftboxes.                                 !
    !======================================================================!
    ! Arguments:                                                           !
    ! kinetic_fftbox_batch (output): Batch of fftboxes with kinetic-NGWFs  !
    ! funcs_on_grid (input) : functions on this proc in ppd representation !
    ! fbasis (input)        : Function basis type for this function set    !
    ! batch_size (input)    : Number of fftboxes in the batch              !
    ! local_start(input)    : First function of the batch                  !
    ! local_end  (input)    : Last function of the batch                   !
    !======================================================================!
    ! Written by Chris-Kriton Skylaris on 20/11/2003.                      !
    ! Modified to use new basis routines by Nicholas Hine, May 2008        !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type.  !
    ! Modified to not use workspace_mod by Nicholas Hine, November 2009.   !
    ! OpenMP parallelised by Nicholas Hine in May 2013.                    !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.  !
    ! Modified by Jose M Escartin on 28/08/2017 to loop over complex       !
    ! functions individually.                                              !
    !======================================================================!

    use basis, only: basis_copy_function_to_box
    use constants, only: DP
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
        data_fftbox_alloc, data_fftbox_dealloc, data_fftbox_copy
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box, fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_threads_fftbox
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    integer, intent(in) :: batch_size, local_start, local_end
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(in) :: funcs_on_grid
    type(FFTBOX_DATA), intent(inout) :: kinetic_fftbox_batch(batch_size)
    integer, intent(in) :: ket_start_in_box(3,batch_size)

    ! Local Variables
    integer :: row1, row2
    integer :: iii
    integer :: ierr
    integer :: batch_count
    type(FFTBOX_DATA) :: work_box1, work_box2
    complex(kind=DP), allocatable :: zwork_box(:,:,:)

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(zwork_box,iii,batch_count, ierr,row1,row2) &
!$OMP FIRSTPRIVATE(work_box1, work_box2) &
!$OMP SHARED(ket_start_in_box, fftbox, local_start, &
!$OMP      cell, local_end, fbasis, funcs_on_grid, kinetic_fftbox_batch, &
!$OMP      pub_threads_fftbox,pub_threads_num_fftboxes)

    ! ndmh: allocate workspace
    call data_fftbox_alloc(work_box1, fftbox%total_ld1, fftbox%total_ld2, &
         fftbox%total_pt3, iscmplx=funcs_on_grid%iscmplx)

    ! agrecocmplx: if using complex NGWFs, fourier transform 1 at a time,
    ! otherwise transform 2 at a time

    if (funcs_on_grid%iscmplx) then

!$OMP DO
       ! jme: if functions are complex, we address only one NGWF per
       ! iteration rather than the 2 we used to, so that
       ! (a) the balance between threads is better (before, a thread
       ! could still have two to go while another was done) and
       ! (b) the memory requirement is smaller
       do iii =local_start, local_end, 1
          row1 =iii

          batch_count = iii-local_start+1

          ! ndmh: copy first function to FFTbox
          call basis_copy_function_to_box(work_box1, &
               ket_start_in_box(1,batch_count), &
               ket_start_in_box(2,batch_count), &
               ket_start_in_box(3,batch_count), &
               fbasis%tight_boxes(row1), funcs_on_grid, &
               fbasis%spheres(row1), cell, fftbox)

           ! agrecocmplx: FFT row1 to reciprocal space
           call fourier_apply_box('Coarse', 'Forward', work_box1%z, &
                omp=pub_threads_fftbox)

           ! agrecocmplx: Apply kinetic energy operator
           call kinetic_apply_on_box(work_box1%z,fftbox)

           ! agrecocmplx: FFT back to real space
           call fourier_apply_box('Coarse', 'Backward', work_box1%z, &
                omp=pub_threads_fftbox)

           ! jme: with Intel Fortran v17.0.4, direct copy seems to require
           ! a much higher OMP_STACKSIZE than call to data_fftbox_copy
!           kinetic_fftbox_batch(batch_count)%z(:,:,:) = work_box1%z
           call data_fftbox_copy(kinetic_fftbox_batch(batch_count),work_box1)

       end do
!$OMP END DO

    else ! funcs_on_grid are real

       ! allocate additional workspace required for real functions
       allocate(zwork_box(fftbox%total_ld1,fftbox%total_ld2, &
            fftbox%total_pt3),stat=ierr)
       call utils_alloc_check('kinetic_apply_to_func_batch','zwork_box',ierr)
       call data_fftbox_alloc(work_box2, fftbox%total_ld1, fftbox%total_ld2, &
            fftbox%total_pt3, iscmplx=.false.)

!$OMP DO
       do iii =local_start, local_end, 2
          row1 = iii ; row2 = iii +1

          batch_count = iii-local_start+1

          ! ndmh: copy first function to FFTbox
          call basis_copy_function_to_box(work_box1, &
               ket_start_in_box(1,batch_count), &
               ket_start_in_box(2,batch_count), &
               ket_start_in_box(3,batch_count), &
               fbasis%tight_boxes(row1), funcs_on_grid, &
               fbasis%spheres(row1), cell, fftbox)

          if (row2 <= local_end) then

             batch_count = iii-local_start+2

             ! ndmh: copy second function to FFTbox
             call basis_copy_function_to_box(work_box2, &
                  ket_start_in_box(1,batch_count), &
                  ket_start_in_box(2,batch_count), &
                  ket_start_in_box(3,batch_count), &
                  fbasis%tight_boxes(row2), funcs_on_grid, &
                  fbasis%spheres(row2), cell, fftbox)

          else
             ! agrecocmplx
             call data_set_to_zero(work_box2)
          end if

          ! FFT 'row' functions to reciprocal space:
          ! aam: use new fourier routines
          call fourier_apply_box_pair('Coarse', 'Forward', work_box1%d, &
               work_box2%d, zwork_box, omp=pub_threads_fftbox)

          ! Apply kinetic energy operator
          call kinetic_apply_on_box(zwork_box,fftbox)

          ! FFT back to real space:
          call fourier_apply_box('Coarse', 'Backward', zwork_box, &
               omp=pub_threads_fftbox)

          ! cks: put functions in batch of fftboxes after the kinetic
          ! cks: operator has been applied on them
          batch_count = iii - local_start + 1
          kinetic_fftbox_batch(batch_count)%d(:,:,:) = real(zwork_box,kind=DP)
          if (row2 <= local_end) then
             batch_count = iii - local_start + 2
             kinetic_fftbox_batch(batch_count)%d(:,:,:) = aimag(zwork_box)
          end if

       end do
!$OMP END DO

       ! deallocate workspaces
       call data_fftbox_dealloc(work_box2)
       deallocate(zwork_box,stat=ierr)
       call utils_dealloc_check('kinetic_apply_to_func_batch','zwork_box',ierr)

    end if

    ! ndmh: deallocate workspace
    call data_fftbox_dealloc(work_box1)

!$OMP END PARALLEL

  end subroutine kinetic_apply_to_func_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine kinetic_apply_to_funcs_on_grid(kin_funcs_on_grid, &  ! output
       funcs_on_grid, fbasis, cell, fftbox)                       ! input

    !======================================================================!
    ! This subroutine applies the kinetic energy operator to a batch       !
    ! of functions (eg NGWFs) in fftboxes.                                 !
    !======================================================================!
    ! Arguments:                                                           !
    ! kin_funcs_on_grid (output): kinetic energy operator applied to funcs !
    ! funcs_on_grid (input) : functions on this proc in ppd representation !
    ! fbasis (input)        : Function basis type for this function set    !
    ! batch_size (input)    : Number of fftboxes in the batch              !
    ! local_start(input)    : First function of the batch                  !
    ! local_end  (input)    : Last function of the batch                   !
    !======================================================================!
    ! Written by Nicholas Hine, March 2011.                                !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.  !
    !======================================================================!

    ! agrecocmplx
    use basis, only: basis_ket_start_wrt_fftbox, basis_copy_function_to_box, &
         basis_extract_function_from_box
    use constants, only: DP
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_set_to_zero, &
         data_fftbox_alloc, data_fftbox_dealloc
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box, fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(inout) :: kin_funcs_on_grid
    type(FUNCTIONS), intent(in) :: funcs_on_grid

    ! Local Variables
    integer :: row1, row2
    integer :: iii
    integer :: ierr
    integer :: row_start1, row_start2, row_start3
    type(FFTBOX_DATA) :: dwork_box1, dwork_box2
    complex(kind=DP), allocatable :: zwork_box(:,:,:)


    ! ndmh: allocate workspace
    call data_fftbox_alloc(dwork_box1, fftbox%total_ld1, &
        fftbox%total_ld2, fftbox%total_pt3, iscmplx=funcs_on_grid%iscmplx)
    call data_fftbox_alloc(dwork_box2, fftbox%total_ld1, &
        fftbox%total_ld2, fftbox%total_pt3, iscmplx=funcs_on_grid%iscmplx)
    allocate(zwork_box(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_apply_to_funcs_on_grid','zwork_box',ierr)

    do iii=1,fbasis%proc_num, 2
       row1 = iii
       row2 = iii + 1

       ! cks: determine where 'row' tightbox begins wrt fftbox
       call basis_ket_start_wrt_fftbox(row_start1, row_start2, row_start3, &
            fftbox%total_pt1, fftbox%total_pt2, fftbox%total_pt3, &
            fftbox)

       ! ndmh: copy first function to FFTbox
       ! agrecocmplx: call modified routine in basis_new
       call basis_copy_function_to_box(dwork_box1, &
            row_start1, row_start2, row_start3,&
            fbasis%tight_boxes(row1), funcs_on_grid, &
            fbasis%spheres(row1), cell, fftbox)

       if (row2 <= fbasis%proc_num) then

          ! cks: determine where 'row' tightbox begins wrt fftbox
          call basis_ket_start_wrt_fftbox(row_start1, row_start2, row_start3, &
               fftbox%total_pt1, fftbox%total_pt2, fftbox%total_pt3, &
               fftbox)

          ! ndmh: copy second function to FFTbox
          ! agrecocmplx: call modified routine in basis_new
          call basis_copy_function_to_box(dwork_box2, &
               row_start1, row_start2, row_start3,&
               fbasis%tight_boxes(row2), funcs_on_grid, &
               fbasis%spheres(row2), cell, fftbox)

       else
          ! agrecocmplx
          call data_set_to_zero(dwork_box2)
       end if

       ! agrecocmplx: if using complex NGWFs, fourier
       ! transform only one at a time
       if (.not.funcs_on_grid%iscmplx) then
           ! FFT 'row' functions to reciprocal space:
           ! aam: use new fourier routines
           call fourier_apply_box_pair('Coarse', 'Forward', dwork_box1%d, &
                dwork_box2%d, zwork_box)

           ! Apply kinetic energy operator
           call kinetic_apply_on_box(zwork_box,fftbox)

           ! FFT back to real space:
           call fourier_apply_box_pair('Coarse', 'Backward', dwork_box1%d, &
                dwork_box2%d, zwork_box)
       else
           ! agrecocmplx: first function
           call fourier_apply_box('Coarse', 'Forward', dwork_box1%z)

           call kinetic_apply_on_box(dwork_box1%z,fftbox)

           call fourier_apply_box('Coarse', 'Backward', dwork_box1%z)

           ! agrecocmplx: second function only if we have it
           if (row2 <= fbasis%proc_num) then
               call fourier_apply_box('Coarse', 'Forward', dwork_box2%z)

               call kinetic_apply_on_box(dwork_box2%z,fftbox)

               call fourier_apply_box('Coarse', 'Backward', dwork_box2%z)
           end if

       end if

       ! ndmh: extract first function from FFTbox
       ! agrecocmplx: call modified routine in basis_new
       call basis_extract_function_from_box(kin_funcs_on_grid, &
            dwork_box1, fbasis%spheres(row1), fbasis%tight_boxes(row1), &
            row_start1, row_start2, row_start3, fbasis%spheres(row1)%offset, &
            cell, fftbox)

       if (row2 <= fbasis%proc_num) then

          ! ndmh: extract second function from FFTbox
          ! agrecocmplx: call modified routine in basis_new
          call basis_extract_function_from_box(kin_funcs_on_grid, &
               dwork_box2, fbasis%spheres(row2), fbasis%tight_boxes(row2), &
               row_start1, row_start2, row_start3, fbasis%spheres(row2)%offset, &
               cell, fftbox)

       else
          ! agrecocmplx
          call data_set_to_zero(dwork_box2)
       end if

    end do

    ! ndmh: deallocate workspace
    deallocate(zwork_box,stat=ierr)
    call utils_dealloc_check('kinetic_apply_to_funcs_on_grid','zwork_box',ierr)
    !deallocate(dwork_box2,stat=ierr)
    !call utils_dealloc_check('kinetic_apply_to_funcs_on_grid','dwork_box2',ierr)
    !deallocate(dwork_box1,stat=ierr)
    !call utils_dealloc_check('kinetic_apply_to_funcs_on_grid','dwork_box1',ierr)
    ! agrecocmplx: deallocate using appropriate routine
    call data_fftbox_dealloc(dwork_box2)
    call data_fftbox_dealloc(dwork_box1)

  end subroutine kinetic_apply_to_funcs_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine kinetic_apply_on_box(data_complex,fftbox)

    !=========================================================!
    ! This subroutine applies the kinetic energy operator in  !
    ! reciprocal space to a function in an FFTbox.            !
    !---------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris in 2000.    !
    ! Rewritten by Peter D. Haynes                            !
    !=========================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
!!$ use rundat, only: pub_threads_per_fftbox

    implicit none

    ! Argument
    complex(kind=DP), intent(inout) :: data_complex(:,:,:)
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local variable
    integer :: i1,i2,i3

    ! -------------------------------------------------------------------------

    ! apply the kinetic energy operator
!!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!!$OMP PRIVATE(i1, i2, i3) &
!!$OMP SHARED(data_complex, fftbox, pub_threads_per_fftbox)
!!$OMP DO
    do i3=1,fftbox%total_pt3
       do i2=1,fftbox%total_pt2
          do i1=1,fftbox%total_pt1
             data_complex(i1,i2,i3) = data_complex(i1,i2,i3) * &
                  fftbox%recip_grid(5,i1,i2,i3)
          end do
       end do
    end do
!!$OMP END PARALLEL

  end subroutine kinetic_apply_on_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine kinetic_grad_to_func_batch(grad_fftbox_batch, &      ! output
       funcs_on_grid, fbasis, cell, fftbox, ket_start_in_box, &   ! input
       batch_size, local_start, local_end)                        ! input

    !=========================================================================!
    ! This subroutine places a batch of functions (eg NGWFs) in FFTboxes and  !
    ! applies the grad operator to them, by transforming to reciprocal space, !
    ! multiplying by iG and transforming back.                                !
    !=========================================================================!
    ! Arguments:                                                              !
    ! grad_fftbox_batch (output): Batch of fftboxes                           !
    ! funcs_on_grid (input) : functions on this proc in ppd representation    !
    ! fbasis (input)        : Function basis type this function set           !
    ! batch_size (input)    : Number of fftboxes in the batch                 !
    ! local_start(input)    : First function of the batch                     !
    ! local_end  (input)    : Last function of the batch                      !
    !=========================================================================!
    ! Written by Peter Haynes on 16/7/2005.                                   !
    ! Based on kinetic_apply_to_ngwf_batch by Chris-Kriton Skylaris on        !
    !    20/11/2003.                                                          !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type.     !
    ! OpenMP parallelised by Nicholas Hine in May 2013.                       !
    ! Data parallel OpenMP by Karl Wilkinson in July 2013.                    !
    ! Modified by Andrea Greco on 14/06/2015 to allow use of complex NGWFs.   !
    !=========================================================================!

    ! agrecocmplx
    use basis, only: basis_copy_function_to_box
    use constants, only: DP
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, data_fftbox_alloc, &
         data_fftbox_dealloc, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box, fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_threads_fftbox
!$  use rundat, only: pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    integer, intent(in) :: batch_size, local_start, local_end
    type(FUNC_BASIS), intent(in) :: fbasis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(in) :: funcs_on_grid
    type(FFTBOX_DATA), intent(inout) :: grad_fftbox_batch(batch_size, 3)
    integer, intent(in) :: ket_start_in_box(3,batch_size)

    ! Local Variables
    integer :: row1, row2
    integer :: iii
    integer :: batch_count
    integer :: ierr
    integer :: dim
    type(FFTBOX_DATA) :: dwork_box1, dwork_box2
    complex(kind=DP), dimension(:,:,:), allocatable :: zwork_box
    complex(kind=DP), dimension(:,:,:,:), allocatable :: zwork_grad_box
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    loc_cmplx = funcs_on_grid%iscmplx

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(zwork_box,zwork_grad_box, iii,batch_count,ierr,row1,row2) &
!$OMP FIRSTPRIVATE(dwork_box1,dwork_box2) &
!$OMP SHARED(ket_start_in_box, fftbox, local_start, &
!$OMP      cell, local_end, fbasis, funcs_on_grid, grad_fftbox_batch, &
!$OMP      pub_threads_fftbox, pub_threads_num_fftboxes, loc_cmplx)

    ! ndmh: allocate workspace
    !allocate(dwork_box1(fftbox%total_ld1,fftbox%total_ld2, &
    !     fftbox%total_pt3),stat=ierr)
    !call utils_alloc_check('kinetic_grad_to_func_batch','dwork_box1',ierr)
    !allocate(dwork_box2(fftbox%total_ld1,fftbox%total_ld2, &
    !     fftbox%total_pt3),stat=ierr)
    !call utils_alloc_check('kinetic_grad_to_func_batch','dwork_box2',ierr)
    ! agrecocmplx: allocate using appropriate routine
    call data_fftbox_alloc(dwork_box1, fftbox%total_ld1, &
         fftbox%total_ld2, fftbox%total_pt3, iscmplx=loc_cmplx)
    call data_fftbox_alloc(dwork_box2, fftbox%total_ld1, &
         fftbox%total_ld2, fftbox%total_pt3, iscmplx=loc_cmplx)
    allocate(zwork_box(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_grad_to_func_batch','zwork_box',ierr)
    allocate(zwork_grad_box(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3,3),stat=ierr)
    call utils_alloc_check('kinetic_grad_to_func_batch','zwork_grad_box',ierr)

    zwork_box = (0.0_DP,0.0_DP)

!$OMP DO
    do iii=local_start, local_end, 2
       row1=iii ; row2=iii+1

       batch_count = iii-local_start+1

       ! ndmh: use new basis_copy routine
       ! agrecocmplx: use modified version in basis_new
       call basis_copy_function_to_box(dwork_box1, &
            ket_start_in_box(1,batch_count), &
            ket_start_in_box(2,batch_count), &
            ket_start_in_box(3,batch_count), &
            fbasis%tight_boxes(row1), funcs_on_grid, &
            fbasis%spheres(row1),cell, fftbox)

       if (row2 <= local_end) then

          batch_count = iii-local_start+2

          ! ndmh: use new basis_copy routine
          ! agrecocmplx: use modified version in basis_new
          call basis_copy_function_to_box(dwork_box2, &
               ket_start_in_box(1,batch_count), &
               ket_start_in_box(2,batch_count), &
               ket_start_in_box(3,batch_count), &
               fbasis%tight_boxes(row2), funcs_on_grid, &
               fbasis%spheres(row2), cell, fftbox)
       else
          ! agrecocmplx
          call data_set_to_zero(dwork_box2)
       end if

       ! FFT 'row' functions to reciprocal space:
       ! aam: use new fourier routines
       ! agrecocmplx: if using complex functions, perform one FFT at a time
       if (.not.loc_cmplx) then
           call fourier_apply_box_pair('Coarse', 'Forward', dwork_box1%d, &
                dwork_box2%d, zwork_box, omp=pub_threads_fftbox)

           call kinetic_grad_on_box(zwork_box,zwork_grad_box,fftbox)

           do dim=1,3
              ! FFT back to real space:
              call fourier_apply_box('Coarse','Backward',zwork_grad_box(:,:,:,dim), &
                                     omp=pub_threads_fftbox)
           end do

           ! cks: put functions in batch of fftboxes after the grad
           ! cks: operator has been applied on them
           batch_count = iii - local_start + 1
           ! agrecocmplx
           do dim=1,3
               grad_fftbox_batch(batch_count,dim)%d(:,:,:) = &
                    real(zwork_grad_box(:,:,:,dim),kind=DP)
           end do
           if (row2 <= local_end) then
              batch_count = iii - local_start + 2
              ! agrecocmplx
              do dim=1,3
                  grad_fftbox_batch(batch_count,dim)%d(:,:,:) = &
                      aimag(zwork_grad_box(:,:,:,dim))
              end do
           end if

       ! agrecocmplx: treat complex case separately
       else
           call fourier_apply_box('Coarse', 'Forward', dwork_box1%z, &
               omp=pub_threads_fftbox)

           call kinetic_grad_on_box(dwork_box1%z, zwork_grad_box,fftbox)

           do dim=1,3
               ! FFT back to real space:
               call fourier_apply_box('Coarse','Backward',zwork_grad_box(:,:,:,dim), &
                                       omp=pub_threads_fftbox)
           end do

           ! cks: put functions in batch of fftboxes after the grad
           ! cks: operator has been applied on them
           batch_count = iii - local_start + 1
           ! agrecocmplx
           do dim=1,3
               grad_fftbox_batch(batch_count,dim)%z(:,:,:) = &
                   zwork_grad_box(:,:,:,dim)
           end do
           ! agrecocmplx: second function only if we have it
           if (row2 <= local_end) then
               call fourier_apply_box('Coarse', 'Forward', dwork_box2%z, &
                      omp=pub_threads_fftbox)

               call kinetic_grad_on_box(dwork_box2%z, zwork_grad_box,fftbox)

               do dim=1,3
                   ! FFT back to real space:
                   call fourier_apply_box('Coarse','Backward',zwork_grad_box(:,:,:,dim), &
                                           omp=pub_threads_fftbox)

               end do

               batch_count = iii - local_start + 2

               do dim=1,3
                   grad_fftbox_batch(batch_count,dim)%z(:,:,:) = &
                       zwork_grad_box(:,:,:,dim)
               end do

           end if

       end if

    end do
!$OMP END DO

    ! ndmh: deallocate workspace
    deallocate(zwork_grad_box,stat=ierr)
    call utils_dealloc_check('kinetic_grad_to_func_batch','zwork_grad_box',ierr)
    deallocate(zwork_box,stat=ierr)
    call utils_dealloc_check('kinetic_grad_to_func_batch','zwork_box',ierr)
    !deallocate(dwork_box2,stat=ierr)
    !call utils_dealloc_check('kinetic_grad_to_func_batch','dwork_box2',ierr)
    !deallocate(dwork_box1,stat=ierr)
    !call utils_dealloc_check('kinetic_grad_to_func_batch','dwork_box1',ierr)
    ! agrecocmplx: deallocate using appropriate routines
    call data_fftbox_dealloc(dwork_box1)
    call data_fftbox_dealloc(dwork_box2)

!$OMP END PARALLEL

  end subroutine kinetic_grad_to_func_batch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine kinetic_grad_to_funcs_on_grid(grad_funcs_on_grid, &  ! output
       funcs_on_grid, fbasis, cell, fftbox)                       ! input
    !======================================================================!
    ! This subroutine applies the grad energy operator to NGWFs in PPD     !
    ! representation and returns the gradients in PPD represenation.       !
    ! The gradient is applied in an FFTbox, but subsequently "shaved" to   !
    ! the NGWF localization sphere when placed back in PPD representation. !
    !======================================================================!
    ! Arguments:                                                           !
    ! TODO Add details for missing arguments
    ! grad_funcs_on_grid (output): grad operator applied to funcs          !
    ! funcs_on_grid (input) : functions on this proc in ppd representation !
    ! fbasis (input)        : Function basis type for this function set    !
    !======================================================================!
    ! Modified version of kinetic_apply_to_funcs_on_grid, created by       !
    ! J. C. Womack, Jan 2016.                                              !
    !======================================================================!
    ! kinetic_apply_to_funcs_on_grid:                                      !
    ! Written by Nicholas Hine, March 2011.                                !
    ! Modified by Andrea Greco in May 2015 to allow use of complex NGWFs.  !
    !======================================================================!

    ! agrecocmplx
    use basis, only: basis_ket_start_wrt_fftbox, basis_copy_function_to_box, &
         basis_extract_function_from_box
    use constants, only: DP
    use datatypes, only: FUNCTIONS, FFTBOX_DATA, &
         data_fftbox_alloc, data_fftbox_dealloc, data_set_to_zero
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box, fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(FUNCTIONS), intent(inout) :: grad_funcs_on_grid(3)
    type(FUNCTIONS), intent(in) :: funcs_on_grid

    ! Local Variables
    integer :: row1, row2
    integer :: iii
    integer :: dim
    integer :: ierr
    integer :: row_start1, row_start2, row_start3
    ! agrecocmplx: use new FFTBOX_DATA type
    type(FFTBOX_DATA) :: dwork_box1(3), dwork_box2(3)
    complex(kind=DP), allocatable :: zwork_box(:,:,:)
    complex(kind=DP), allocatable :: zwork_grad_box(:,:,:,:)


    ! ndmh: allocate workspace
    ! agrecocmplx: allocate using appropriate routine
    do dim = 1, 3
       call data_fftbox_alloc(dwork_box1(dim), fftbox%total_ld1, &
            fftbox%total_ld2, fftbox%total_pt3, iscmplx=funcs_on_grid%iscmplx)
       call data_fftbox_alloc(dwork_box2(dim), fftbox%total_ld1, &
            fftbox%total_ld2, fftbox%total_pt3, iscmplx=funcs_on_grid%iscmplx)
    end do
    allocate(zwork_box(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_grad_to_funcs_on_grid','zwork_box',ierr)
    allocate(zwork_grad_box(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3,3),stat=ierr)
    call utils_alloc_check('kinetic_grad_to_funcs_on_grid','zwork_grad_box',&
         ierr)

    do iii=1,fbasis%proc_num, 2
       row1 = iii
       row2 = iii + 1

       ! cks: determine where 'row' tightbox begins wrt fftbox
       call basis_ket_start_wrt_fftbox(row_start1, row_start2, row_start3, &
            fftbox%total_pt1, fftbox%total_pt2, fftbox%total_pt3, &
            fftbox)

       ! ndmh: copy first function to FFTbox
       ! agrecocmplx: call modified routine in basis_new
       call basis_copy_function_to_box(dwork_box1(1), &
            row_start1, row_start2, row_start3,&
            fbasis%tight_boxes(row1), funcs_on_grid, &
            fbasis%spheres(row1), cell, fftbox)

       if (row2 <= fbasis%proc_num) then

          ! cks: determine where 'row' tightbox begins wrt fftbox
          call basis_ket_start_wrt_fftbox(row_start1, row_start2, row_start3, &
               fftbox%total_pt1, fftbox%total_pt2, fftbox%total_pt3, &
               fftbox)

          ! ndmh: copy second function to FFTbox
          ! agrecocmplx: call modified routine in basis_new
          call basis_copy_function_to_box(dwork_box2(1), &
               row_start1, row_start2, row_start3,&
               fbasis%tight_boxes(row2), funcs_on_grid, &
               fbasis%spheres(row2), cell, fftbox)

       else
          ! agrecocmplx
          call data_set_to_zero(dwork_box2(1))
       end if

       ! agrecocmplx: if using complex NGWFs, fourier
       ! transform only one at a time
       if (.not.funcs_on_grid%iscmplx) then
           ! FFT 'row' functions to reciprocal space:
           ! aam: use new fourier routines
           call fourier_apply_box_pair('Coarse', 'Forward', dwork_box1(1)%d, &
                dwork_box2(1)%d, zwork_box)

           ! Apply kinetic energy operator
           call kinetic_grad_on_box(zwork_box,zwork_grad_box,fftbox)

           do dim = 1, 3
              ! FFT back to real space:
              call fourier_apply_box_pair('Coarse', 'Backward', &
                   dwork_box1(dim)%d, dwork_box2(dim)%d, &
                   zwork_grad_box(:,:,:,dim))
           end do
       else
           ! agrecocmplx: first function
           call fourier_apply_box('Coarse', 'Forward', dwork_box1(1)%z)

           call kinetic_grad_on_box(dwork_box1(1)%z,zwork_grad_box,fftbox)

           do dim = 1, 3
              call fourier_apply_box('Coarse', 'Backward', &
                   dwork_box1(dim)%z,zwork_grad_box(:,:,:,dim))
           end do

           ! agrecocmplx: second function only if we have it
           if (row2 <= fbasis%proc_num) then
               call fourier_apply_box('Coarse', 'Forward', dwork_box2(1)%z)

               call kinetic_grad_on_box(dwork_box2(1)%z,zwork_grad_box,fftbox)

               call fourier_apply_box('Coarse', 'Backward', &
                    dwork_box2(dim)%z,zwork_grad_box(:,:,:,dim))
           end if

       end if

       ! ndmh: extract first function from FFTbox
       ! agrecocmplx: call modified routine in basis_new
       do dim = 1, 3
          call basis_extract_function_from_box(grad_funcs_on_grid(dim), &
               dwork_box1(dim), fbasis%spheres(row1), fbasis%tight_boxes(row1),&
               row_start1, row_start2, row_start3, fbasis%spheres(row1)%offset,&
               cell, fftbox)

          if (row2 <= fbasis%proc_num) then

             ! ndmh: extract second function from FFTbox
             ! agrecocmplx: call modified routine in basis_new
             call basis_extract_function_from_box(grad_funcs_on_grid(dim), &
                  dwork_box2(dim), fbasis%spheres(row2), &
                  fbasis%tight_boxes(row2), &
                  row_start1, row_start2, row_start3, &
                  fbasis%spheres(row2)%offset, cell, fftbox)

          else
             ! agrecocmplx
             call data_set_to_zero(dwork_box2(dim))
          end if
       end do

    end do

    ! ndmh: deallocate workspace
    deallocate(zwork_grad_box)
    call utils_alloc_check('kinetic_grad_to_funcs_on_grid','zwork_grad_box',&
         ierr)
    deallocate(zwork_box,stat=ierr)
    call utils_dealloc_check('kinetic_grad_to_funcs_on_grid','zwork_box',ierr)
    do dim = 1, 3
       call data_fftbox_dealloc(dwork_box2(dim))
       call data_fftbox_dealloc(dwork_box1(dim))
    end do

  end subroutine kinetic_grad_to_funcs_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine kinetic_grad_on_cart_box(data_complex_in,data_complex_out,fftbox)

    !=========================================================!
    ! This subroutine applies the Cartesian components of the !
    ! grad operator in reciprocal space to the corresponding  !
    ! Cartesian components of FFTboxes.                       !
    ! This is useful as an intermediate step in taking the    !
    ! divergence of a Cartesian vector field in reciprocal    !
    ! space, where the gradient of each of the x, y and z     !
    ! components in reciprocal space is to be transformed     !
    ! back to real space, where they can then be summed.      !
    !---------------------------------------------------------!
    ! Modified version of kinetic_grad_on_box, created by     !
    ! James C. Womack, 2015                                   !
    !---------------------------------------------------------!
    ! Written by Peter Haynes on 16/7/2005.                   !
    ! Based on kinetic_apply_on_box by Chris-Kriton Skylaris  !
    ! in 2000.                                                !
    !=========================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_threads_per_fftbox
    use utils, only: utils_use_var

    implicit none

    ! Argument
    complex(kind=DP), intent(in) :: data_complex_in(:,:,:,:)
    complex(kind=DP), intent(out) :: data_complex_out(:,:,:,:)
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local variable
    integer :: i1,i2,i3,dim

    ! -------------------------------------------------------------------------

    call utils_use_var(pub_threads_per_fftbox)

    ! TODO: Consider the number of threads per fftbox --- better to rearrange
    !       the loops so loop over dim is outside?
    ! apply the grad operator
!!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!!$OMP PRIVATE(i1, i2, i3) &
!!$OMP SHARED(data_complex_in, data_complex_out, fftbox, &
!!$OMP      pub_threads_per_fftbox)
!!$OMP DO
    do i3=1,fftbox%total_pt3
       do i2=1,fftbox%total_pt2
          do i1=1,fftbox%total_pt1
             do dim=1,3
                data_complex_out(i1,i2,i3,dim) = &
                     data_complex_in(i1,i2,i3,dim) * &
                     cmplx(0.0_DP,fftbox%recip_grid(dim,i1,i2,i3),kind=DP)
             end do
          end do
       end do
    end do
!!$OMP END DO
!!$OMP END PARALLEL

  end subroutine kinetic_grad_on_cart_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine kinetic_grad_on_box(data_complex_in,data_complex_out,fftbox)

    !=========================================================!
    ! This subroutine applies the grad operator in            !
    ! reciprocal space to a function in an FFTbox.            !
    !---------------------------------------------------------!
    ! Written by Peter Haynes on 16/7/2005.                   !
    ! Based on kinetic_apply_on_box by Chris-Kriton Skylaris  !
    ! in 2000.                                                !
    !=========================================================!

    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use rundat, only: pub_threads_per_fftbox
    use utils, only: utils_use_var

    implicit none

    ! Argument
    complex(kind=DP), intent(in) :: data_complex_in(:,:,:)
    complex(kind=DP), intent(out) :: data_complex_out(:,:,:,:)
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local variable
    integer :: i1,i2,i3,dim

    ! -------------------------------------------------------------------------

    call utils_use_var(pub_threads_per_fftbox)

    ! apply the grad operator
!!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) DEFAULT(none) &
!!$OMP PRIVATE(i1, i2, i3) &
!!$OMP SHARED(data_complex_in, data_complex_out, fftbox, &
!!$OMP      pub_threads_per_fftbox)
!!$OMP DO
    do i3=1,fftbox%total_pt3
       do i2=1,fftbox%total_pt2
          do i1=1,fftbox%total_pt1
             do dim=1,3
                data_complex_out(i1,i2,i3,dim) = data_complex_in(i1,i2,i3) * &
                     cmplx(0.0_DP,fftbox%recip_grid(dim,i1,i2,i3),kind=DP)
             end do
          end do
       end do
    end do
!!$OMP END DO
!!$OMP END PARALLEL

  end subroutine kinetic_grad_on_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef GPU_PGI_KINETIC
  subroutine kinetic_gpu_app2_func_batch(kinetic_fftbox_batch, &  ! output
       funcs_on_grid, fbasis, ket_start_in_box, batch_size, &     ! input
       local_start, local_end)                                    ! input

    !======================================================================!
    ! This subroutine applies the kinetic energy operator to a batch       !
    ! of functions (eg NGWFs) in fftboxes.                                 !
    !======================================================================!
    ! Arguments:                                                           !
    ! kinetic_fftbox_batch (output): Batch of fftboxes with kinetic-NGWFs  !
    ! funcs_on_grid (input) : functions on this proc in ppd representation !
    ! fbasis (input)        : Function basis type for this function set    !
    ! batch_size (input)    : Number of fftboxes in the batch              !
    ! local_start(input)    : First function of the batch                  !
    ! local_end  (input)    : Last function of the batch                   !
    !======================================================================!
    ! Written by Chris-Kriton Skylaris on 20/11/2003.                      !
    ! Modified to use new basis routines by Nicholas Hine, May 2008        !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type.  !
    ! Modified to not use workspace_mod by Nicholas Hine, November 2009.   !
    ! Modified to run on GPUs by Karl Wilkinson, December 2012.            !
    !======================================================================!

    use basis, only: basis_copy_function_to_box
    use constants, only: DP, stdout
    use fourier, only: fourier_apply_box_pair, fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_debug_on_root
    use simulation_cell, only: cell, fftbox
    use utils, only: utils_alloc_check, utils_dealloc_check
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor !! External dependency
    use fourier_gpu_wrapper_mod

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    integer, intent(in) :: batch_size, local_start, local_end
    real(kind=DP), intent(in)  :: &
         funcs_on_grid(fbasis%size_on_grid)
    real(kind=DP), intent(out) :: &
         kinetic_fftbox_batch(fftbox%total_ld1, fftbox%total_ld2, &
         fftbox%total_pt3, batch_size)
    integer, intent(in) :: ket_start_in_box(3,batch_size)

    ! Local Variables
    integer :: row1, row2
    integer :: iii
    integer :: ierr
    integer :: batch_count
    real(kind=DP), allocatable :: dwork_box1(:,:,:), dwork_box2(:,:,:)

    integer :: i1, i2, i3
    real(kind=DP) :: scalefac
    ! kaw: Check for potential re-use and potentially declare/allocate elsewhere.
    !      certainly move outside the batch loop.
    complex(kind=DP), dimension (:,:,:), allocatable :: temp
    integer :: istat

    if (pub_debug_on_root) write(stdout,*) 'DEBUG: kinetic_gpu_app2_func_batch'

    scalefac = 1.0_DP / (fftbox%total_ld1*fftbox%total_ld2*fftbox%total_pt3)

    ! kaw: This should be done at initiailisation time.
    fftbox_recip_grid_gpu = fftbox%recip_grid(5,:,:,:)

    ! kaw: This should be done outside the loop over batches
    allocate(temp(fftbox%total_ld1, fftbox%total_ld2, fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_gpu_app2_func_batch','temp',ierr)

    ! ndmh: allocate workspace
    allocate(dwork_box1(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_gpu_app2_func_batch','dwork_box1',ierr)
    allocate(dwork_box2(fftbox%total_ld1,fftbox%total_ld2, &
         fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_gpu_app2_func_batch','dwork_box2',ierr)

    batch_count =0
    do iii =local_start, local_end, 2
       row1 =iii ; row2 =iii +1

       batch_count = iii-local_start+1

       ! ndmh: copy first function to FFTbox
       call basis_copy_function_to_box(dwork_box1, fftbox%total_ld1, &
            fftbox%total_ld2, fftbox%total_pt3, &
            ket_start_in_box(1,batch_count), &
            ket_start_in_box(2,batch_count), &
            ket_start_in_box(3,batch_count), &
            fbasis%tight_boxes(row1), funcs_on_grid, &
            fbasis%spheres(row1), cell, fftbox)

       if (row2 <= local_end) then

          batch_count = iii-local_start+2

          ! ndmh: copy second function to FFTbox
          call basis_copy_function_to_box(dwork_box2, fftbox%total_ld1, &
               fftbox%total_ld2, fftbox%total_pt3, &
               ket_start_in_box(1,batch_count), &
               ket_start_in_box(2,batch_count), &
               ket_start_in_box(3,batch_count), &
               fbasis%tight_boxes(row2), funcs_on_grid, &
               fbasis%spheres(row2), cell, fftbox)

       else
          dwork_box2 = 0.0_DP
       end if

       ! FFT 'row' functions to reciprocal space:
       ! aam: use new fourier routines
       temp(:,:,:) = cmplx(dwork_box1(:,:,:),dwork_box2(:,:,:),kind=DP)
       coarse_work_gpu = temp

       call cufftExec(cufftplan_coarse,planType,coarse_work_gpu,coarse_work_gpu,CUFFT_FORWARD)

       ! apply the kinetic energy operator
!$acc data region
!$acc region
       do i3=1,fftbox%total_pt3
          do i2=1,fftbox%total_ld2
             do i1=1,fftbox%total_ld1
                coarse_work_gpu(i1,i2,i3) = coarse_work_gpu(i1,i2,i3) * &
                     scalefac*fftbox_recip_grid_gpu(i1,i2,i3)
             end do
          end do
       end do
!$acc end region
!$acc end data region

       ! FFT back to real space:
       call cufftExec(cufftplan_coarse,planType,coarse_work_gpu,coarse_work_gpu,CUFFT_INVERSE)

       ! cks: put functions in batch of fftboxes after the kinetic
       ! cks: operator has been applied on them
       batch_count = batch_count + 1

       temp = coarse_work_gpu
       kinetic_fftbox_batch(:,:,:,batch_count) = real(temp(:,:,:),kind=DP)
       if (row2 <= local_end) then
          batch_count = batch_count + 1
          kinetic_fftbox_batch(:,:,:,batch_count) = aimag(temp(:,:,:))
       end if
    end do

    ! ndmh: deallocate workspace
    deallocate(dwork_box2,stat=ierr)
    call utils_dealloc_check('kinetic_gpu_app2_func_batch','dwork_box2',ierr)
    deallocate(dwork_box1,stat=ierr)
    call utils_dealloc_check('kinetic_gpu_app2_func_batch','dwork_box1',ierr)
    deallocate(temp,stat=ierr)
    call utils_dealloc_check('kinetic_gpu_app2_func_batch','temp',ierr)



  end subroutine kinetic_gpu_app2_func_batch
#else
  subroutine dummy_gpu
  use fourier_gpu_wrapper_mod
  !kaw: This is a dummy use statement to make the check_dependencies script happy.
  end subroutine dummy_gpu
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module kinetic

