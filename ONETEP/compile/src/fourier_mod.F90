! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                                                                             !
! This module contains data structures and routines to perform FFTs within    !
! the ONETEP code                                                             !
!                                                                             !
! Original version by Peter Haynes, 7/7/03                                    !
! Parallel version by Peter Haynes, 21/6/04                                   !
!                                                                             !
! Modified by Nicholas Hine, July 2008 to improve cache performance, reduce   !
! memory usage, and remove unnecessary memory-copying. Also removed           !
! overloaded interface to fourier_apply to make calls explicit.               !
!                                                                             !
! Subsequent modifications by Nicholas Hine and Karl Wilkinson.               !
!=============================================================================!


module fourier

  use constants, only: DP, LONG
  use rundat, only: pub_debug_on_root

  implicit none

  private

  public :: fourier_init_threads
  public :: fourier_exit_threads
  public :: fourier_init_cell
  public :: fourier_exit_cell
  public :: fourier_apply_cell_forward
  public :: fourier_apply_cell_backward
  public :: fourier_filter_cell
  public :: fourier_interpolate_cell

  public :: fourier_init_fftbox
  public :: fourier_init_box
  public :: fourier_apply_box
  public :: fourier_apply_box_pair
  public :: fourier_exit_fftbox
  public :: fourier_exit_box
  public :: fourier_filter
  public :: fourier_interpolate
  public :: fourier_interpolate_product

#ifdef GPU_PGI
#ifdef GPU_SP_TEST
  public :: fourier_gpu_interp_prod_SP
#ifdef GPU_SP_TEST_POT
  public :: fourier_gpu_interpolate_SP
  public :: fourier_gpu_filter_SP
#endif
#endif
  public :: fourier_gpu_interp_prod
  public :: fourier_gpu_interpolate
  public :: fourier_gpu_filter
  public :: fourier_gpu_interpolate_grad
  public :: fourier_gpu_filter_grad
  public :: fourier_gpu_init
  public :: fourier_gpu_exit
#endif

#ifdef MKL_FFTW3
#define FFTW3
#endif
#ifdef FFTW3_NO_OMP
#ifndef FFTW3
#define FFTW3
#endif
#endif
#ifdef FFTW3
#ifndef FFTW3_NO_OMP
#ifndef ParallelFFT
#define ParallelFFT
#endif
#endif
#endif

  ! *********************************
  ! ***  P r i v a t e   d a t a  ***
  ! *********************************

  interface fourier_write
     module procedure fourier_write_real_2
     module procedure fourier_write_real_3
     module procedure fourier_write_complex_2
     module procedure fourier_write_complex_3
  end interface

  interface fourier_filter
     module procedure fourier_filter_1complex
     module procedure fourier_filter_2real
  end interface

  interface fourier_interpolate
     module procedure fourier_interpolate_1complex
     module procedure fourier_interpolate_2real
  end interface

  type serial_fft3d_info
     integer :: n1,n2,n3  ! The dimensions of the FFT grid
     integer :: ld1,ld2   ! The dimensions of the arrays holding the grid
     integer :: iwork_len ! Length of integer workspace
     integer :: dwork_len ! Length of real workspace
     integer :: zwork_len ! Length of complex workspace
     integer, allocatable, dimension(:) :: iwork ! integer workspace
     real(kind=DP), allocatable, dimension(:) :: dwork ! real workspace
     complex(kind=DP), allocatable, dimension(:) :: zwork ! complex workspace
#ifdef FFTW
     integer(kind=LONG) :: forward_plan, backward_plan
#endif
#ifdef FFTW3
     integer(kind=LONG) :: forward_plan, backward_plan
     integer(kind=LONG) :: n1_f_plan, n1_b_plan
     integer(kind=LONG) :: n2_f_plan, n2_b_plan
     integer(kind=LONG) :: n3_f_plan, n3_b_plan
     integer(kind=LONG) :: n3_f_plan_1, n3_b_plan_1
     integer(kind=LONG) :: n3_f_plan_2, n3_b_plan_2
     integer(kind=LONG) :: n3_f_plan_3, n3_b_plan_3

     ! kaw: Plans for data parallel FFT box FFTs
     integer(kind=LONG) :: forward_omp_plan, backward_omp_plan
     integer(kind=LONG) :: n1_f_omp_plan, n1_b_omp_plan
     integer(kind=LONG) :: n2_f_omp_plan, n2_b_omp_plan
     integer(kind=LONG) :: n3_f_omp_plan, n3_b_omp_plan
     integer(kind=LONG) :: n3_f_omp_plan_1, n3_b_omp_plan_1
     integer(kind=LONG) :: n3_f_omp_plan_2, n3_b_omp_plan_2
     integer(kind=LONG) :: n3_f_omp_plan_3, n3_b_omp_plan_3

#endif
#ifdef ACML
     integer :: tbl_size
     complex(kind=DP), allocatable, dimension(:,:) :: tbl
#endif
  end type serial_fft3d_info

  type parallel_fft3d_info
     integer :: n1,n2,n3     ! The dimensions of the whole FFT grid
     integer :: ld1,ld2,ld3  ! The dimensions of the arrays holding the grids
     integer :: iwork_len    ! Length of integer workspaces
     integer :: dwork_len    ! Length of real workspaces
     integer :: zwork_len    ! Length of complex workspaces
     integer :: num12slabs   ! Number of 12-slabs for real space
     integer :: max12slabs   ! Maximum number of 12-slabs
     integer :: num23slabs   ! Number of 23-slabs for real space
     integer :: max23slabs   ! Maximum number of 23-slabs
     integer, allocatable, dimension(:) :: idx12slab ! index of 12-slabs
     integer, allocatable, dimension(:) :: idx23slab ! index of 23-slabs
     complex(kind=DP), allocatable, dimension(:,:) :: buf12slab ! buffer for 12-slab
     complex(kind=DP), allocatable, dimension(:,:) :: buf23slab ! buffer for 23-slab
     complex(kind=DP), allocatable, dimension(:,:,:) :: fsendbuf ! forward send buf
     complex(kind=DP), allocatable, dimension(:,:,:) :: frecvbuf ! forward send buf
     complex(kind=DP), allocatable, dimension(:,:,:) :: bsendbuf ! backw'd recv buf
     complex(kind=DP), allocatable, dimension(:,:,:) :: brecvbuf ! forward recv buf
     integer, allocatable, dimension(:,:) :: iwork ! integer workspaces
     real(kind=DP), allocatable, dimension(:,:) :: dwork ! real workspaces
     complex(kind=DP), allocatable, dimension(:,:) :: zwork ! complex workspaces
#ifdef FFTW
     integer(kind=LONG) :: forward_plan(2), backward_plan(2)
     real(kind=DP), allocatable, dimension(:,:) :: rfftwbuf ! buffer for rFFTw
#endif
#ifdef FFTW3
     integer(kind=LONG) :: forward_plan(2), backward_plan(2)
#ifdef ParallelFFT
     integer(kind=LONG) :: forward_plan_omp(2), backward_plan_omp(2)
#endif
#endif
#ifdef ACML
     integer :: tbl_size
     real(kind=DP), allocatable, dimension(:) :: rbuf ! buffer for zd/dz fft
     complex(kind=DP), allocatable, dimension(:,:) :: tbl
#endif
  end type parallel_fft3d_info

  ! qoh: Information arrays for FFT boxes
  type(serial_fft3d_info) :: internal_box_coarse_info
  type(serial_fft3d_info) :: internal_box_dbl_info

  integer, parameter :: MAX_BOX_INFOS = 3
  type(serial_fft3d_info) :: internal_box_info(MAX_BOX_INFOS)
  integer :: box_info_index = 0

  ! qoh: Information arrays for cell
  integer, parameter :: MAX_CELL_INFOS = 6
  type(parallel_fft3d_info) :: internal_cell_info(MAX_CELL_INFOS)
  integer :: cell_info_index = 0


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
#endif

#ifdef FFTW3
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = 1
  integer, parameter :: FFTW_MEASURE = 0
  integer, parameter :: FFTW_ESTIMATE = 64
#endif

#ifdef ACML
  integer, parameter :: ACML_MODE_PLAN = 100
  integer, parameter :: ACML_MODE_FORWARD = 1
  integer, parameter :: ACML_MODE_BACKWARD = -1
#endif

#ifdef VENDOR
#ifdef ALPHA
  integer, parameter :: dxml_structure_size = 256
  integer, external :: zfft_init_3d,zfft_apply_3d,zfft_exit_3d
  integer, external :: zfft_init_2d,zfft_apply_2d,zfft_exit_2d
  integer, external :: dfft_init,dfft_apply,dfft_exit
#endif
#endif

contains

  ! ***************************************
  ! ***  P u b l i c   r o u t i n e s  ***
  ! ***************************************

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_init_threads

    use constants, only: stdout
    use rundat, only: pub_threads_fftbox

    implicit none

#ifdef FFTW3
#ifndef FFTW3_NO_OMP
    external :: dfftw_init_threads

    ! Local Variables
    integer :: ierr

    ! kaw: This needs to be called before ANY other FFTW routines or the code
    ! will hang.
    if (pub_threads_fftbox) then
        call dfftw_init_threads(ierr)

    if (pub_debug_on_root) write(stdout,*) 'dfftw_init_threads exited with &
        & status:',ierr

    end if
#endif
#endif

  end subroutine fourier_init_threads

  subroutine fourier_exit_threads

    use rundat, only: pub_threads_fftbox

    implicit none

    ! kaw: This needs to be called after ANY other FFTW routines or the code
    ! will hang.
#ifdef FFTW3
#ifndef FFTW3_NO_OMP
    external :: dfftw_cleanup_threads

    if (pub_threads_fftbox) then
        call dfftw_cleanup_threads
    end if
#endif
#endif

  end subroutine fourier_exit_threads

  subroutine fourier_init_fftbox(fftbox)

    use fft_box, only: FFTBOX_INFO
#ifdef BASIC_GPU_FFT
    ! kaw: Added to give access to the rungpu switch. Remove eventually as the
    !      switch should be replaced by compile time flags.
    use fourier_gpu_wrapper_mod
#endif

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(inout) :: fftbox

    ! Local variables

    ! Coarse-Grid FFTbox
    call internal_serial_init(fftbox%total_pt1,fftbox%total_pt2, &
         fftbox%total_pt3,fftbox%total_ld1,fftbox%total_ld2, &
         internal_box_coarse_info)

    ! Double-Grid FFTbox
    call internal_serial_init(fftbox%total_pt1_dbl, &
         fftbox%total_pt2_dbl,fftbox%total_pt3_dbl, &
         fftbox%total_ld1_dbl,fftbox%total_ld2_dbl, &
         internal_box_dbl_info)

#ifdef GPU_PGI
       ! kaw: Set up the FFT plans and arrays for the GPU routines.
       call fourier_gpu_init
#endif
#ifdef BASIC_GPU_FFT
    allocate(coarse_work_gpu(fftbox%total_ld1,fftbox%total_ld2, &
                             fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('fourier_init_fftbox','coarse_work_gpu',ierr)
    call cufftPlan3D(cufftplan_coarse, fftbox%total_pt3, &
                     fftbox%total_ld2,fftbox%total_ld1, &
                     planType)

    allocate(fine_work_gpu(fftbox%total_ld1_dbl,fftbox%total_ld2_dbl, &
                           fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('fourier_init_fftbox','fine_work_gpu',ierr)
    call cufftPlan3D(cufftplan_fine, fftbox%total_pt3_dbl, &
                     fftbox%total_ld2_dbl,fftbox%total_ld1_dbl, &
                     planType)
#endif

  end subroutine fourier_init_fftbox

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_init_box(box)

    use fft_box, only: FFTBOX_INFO
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(inout) :: box

    ! Local variables

    if (box%fft_index/=0) then
       call utils_abort('Error in fourier_init_box: box already initialised')
    end if

    box_info_index = box_info_index + 1

    if (box_info_index>MAX_BOX_INFOS) then
       call utils_abort('Error in fourier_init_box: too many boxes allocated')
    end if

    box%fft_index = box_info_index

    call internal_serial_init(box%total_pt1,box%total_pt2,box%total_pt3, &
         box%total_pt1,box%total_pt2,internal_box_info(box%fft_index))

  end subroutine fourier_init_box

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_init_cell(grid)

    use cell_grid, only: GRID_INFO
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid

    if (grid%fft_index/=0) then
       call utils_abort('Error in fourier_init: grid already initialised')
    end if

    cell_info_index = cell_info_index + 1

    if (cell_info_index>MAX_CELL_INFOS) then
       call utils_abort('Error in fourier_init: too many grids allocated')
    end if

    grid%fft_index = cell_info_index

    call internal_parallel_init(grid,internal_cell_info(grid%fft_index))

  end subroutine fourier_init_cell

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_apply_cell_forward(rspc,gspc,grid)
    !=========================================================================!
    ! Performs a Fast Fourier Transform from real to reciprocal space.        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rspc (input):  Slab-distributed array that contains the real space    !
    !                  data to be transformed.                                !
    !   gspc (output): Slab-distributed array that will contain the result    !
    !                  of the transform in reciprocal space.                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 7/7/03                                         !
    ! Parallel version by Peter Haynes, 21/6/04                               !
    ! Improved for const-correctness and documented, Jacek Dziedzic, 05/2010  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only : comms_free
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: rspc(grid%ld1,grid%ld2,grid%max_slabs12)
    complex(kind=DP), intent(out) :: gspc(grid%ld3,grid%ld2,grid%max_slabs23)

    ! -----------------------------------------------------------------------

    call timer_clock('fourier_apply_cell_forward',1)

    ! Make a license heartbeat call, if ACCELRYS defined
    call make_license_heartbeat

    ! Do the hard work
    call internal_fft3d_cell_forward(internal_cell_info(grid%fft_index),rspc, &
         gspc)

    call timer_clock('fourier_apply_cell_forward',2, &
         calculate_flops(internal_cell_info(grid%fft_index)))

    ! Free up comms resources
    call comms_free

  end subroutine fourier_apply_cell_forward

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_apply_cell_backward(rspc,gspc,grid)
    !=========================================================================!
    ! Performs a Fast Fourier Transform from reciprocal to real space.        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rspc (output): Slab-distributed array that contains the reciprocal    !
    !                  space data to be transformed.                          !
    !   gspc (input):  Slab-distributed array that will contain the result    !
    !                  of the transform in real space.                        !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 7/7/03                                         !
    ! Parallel version by Peter Haynes, 21/6/04                               !
    ! Improved for const-correctness and documented, Jacek Dziedzic, 05/2010  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only : comms_free
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out) :: rspc(grid%ld1,grid%ld2,grid%max_slabs12)
    complex(kind=DP), intent(in) :: gspc(grid%ld3,grid%ld2,grid%max_slabs23)

    ! -----------------------------------------------------------------------

    call timer_clock('fourier_apply_cell_backward',1)

    ! Make a license heartbeat call, if ACCELRYS defined
    call make_license_heartbeat

    ! Do the hard work
    call internal_fft3d_cell_backward(internal_cell_info(grid%fft_index), &
         rspc,gspc)

    call timer_clock('fourier_apply_cell_backward',2, &
         calculate_flops(internal_cell_info(grid%fft_index)))

    ! Free up comms resources
    call comms_free

  end subroutine fourier_apply_cell_backward

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_interpolate_cell(rspc1,rspc2,grid1,grid2,apply_nyquist)
    !=========================================================================!
    ! Performs a Fast Fourier Transform Interpolation from a specified grid   !
    ! to a larger grid.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rspc1 (input):  Slab-distributed array that contains the real space   !
    !                   data on grid 1 to be interpolated to grid 2.          !
    !   rspc2 (output): Slab-distributed array that will contain the result   !
    !                   of the filtering in real space on grid 2.             !
    !   grid1 (input):  GRID_INFO describing grid 1 (the smaller grid).       !
    !   grid2 (input):  GRID_INFO describing grid 2 (the larger grid).        !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 29/06/10                                      !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only : comms_free, comms_recv, comms_send, pub_my_proc_id
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid1
    type(GRID_INFO), intent(in) :: grid2
    real(kind=DP), intent(in) :: rspc1(grid1%ld1,grid1%ld2,grid1%max_slabs12)
    real(kind=DP), intent(out) :: rspc2(grid2%ld1,grid2%ld2,grid2%max_slabs12)
    logical, intent(in), optional :: apply_nyquist

    ! Local Variables
    integer :: ierr
    integer :: islab23,jslab23,kslab23
    integer :: n3,n2,m3,m2,e3,e2,l3,l2,k3,k2
    integer :: j1
    integer :: sendproc, recvproc
    complex(kind=DP), allocatable :: gspc1(:,:,:)
    complex(kind=DP), allocatable :: gspc2(:,:,:)
    complex(kind=DP), allocatable :: gbuf1(:,:)
    logical :: loc_apply_nyquist

    ! -----------------------------------------------------------------------

    call timer_clock('fourier_interpolate_cell',1)

    ! Optional arguments
    loc_apply_nyquist = .false.
    if (present(apply_nyquist)) loc_apply_nyquist = apply_nyquist

    ! Allocate workspaces
    allocate(gspc1(grid1%ld3,grid1%ld2,grid1%max_slabs23),stat=ierr)
    call utils_alloc_check('fourier_interpolate_cell','gspc1',ierr)
    allocate(gspc2(grid2%ld3,grid2%ld2,grid2%max_slabs23),stat=ierr)
    call utils_alloc_check('fourier_interpolate_cell','gspc2',ierr)
    allocate(gbuf1(grid1%ld3,grid1%ld2),stat=ierr)
    call utils_alloc_check('fourier_interpolate_cell','gbuf1',ierr)

    ! Make a license heartbeat call, if ACCELRYS defined
    call make_license_heartbeat

    ! Transform from real to reciprocal space on grid 1
    call internal_fft3d_cell_forward(internal_cell_info(grid1%fft_index), &
         rspc1,gspc1)

    ! Nyquist filter (grid is always going to be even)
    if (grid1%num_slabs23 > 0) then
       gspc1(grid1%n3/2+1,:,:) = (0.0_DP,0.0_DP)
       gspc1(:,grid1%n2/2+1,:) = (0.0_DP,0.0_DP)
    end if
    ! Nyquist filter for last slab
    if (pub_my_proc_id==grid1%proc_slab23(grid1%n1/2+1)) &
         gspc1(:,:,grid1%num_slabs23) = (0.0_DP,0.0_DP)

    ! Scale grid1 recip data
    gspc1 = gspc1 * (real(grid2%n1,kind=DP) * real(grid2%n2,kind=DP) &
         * real(grid2%n3,kind=DP)) / (real(grid1%n1,kind=DP) &
         * real(grid1%n2,kind=DP) * real(grid1%n3,kind=DP))

    ! Find limits in reciprocal i3,i2 directions on grids 1 and 2
    n3 = grid1%n3;         n2 = grid1%n2
    e3 = grid2%n3;         e2 = grid2%n2
    m3 = (n3+1)/2;         m2 = (n2+1)/2
    l3 = n3/2+1+mod(n3,2); l2 = n2/2+1+mod(n2,2)
    k3 = e3-n3+l3;         k2 = e2-n2+l2
    gspc2(:,:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)

    ! Loop over all slabs of grid2 (the larger grid)
    do islab23=1,grid2%max_slabs23

       ! Find if there are any slabs needing to be sent to other procs
       ! on this step
       do jslab23=1,grid1%num_slabs23
          ! Find the global index of this slab, the proc holding this slab
          ! in grid2, and the local index of this slab on that proc
          j1 = jslab23 + grid1%first_slab23(pub_my_proc_id) - 1
          sendproc = grid2%proc_slab23(j1)
          kslab23 = j1 - grid2%first_slab23(sendproc) + 1
          ! If the local index on that proc is the same as the step we
          ! are currently on, then send the slab
          if (kslab23==islab23) then
             if (sendproc/=pub_my_proc_id) then
                call comms_send(sendproc,gspc1(:,:,jslab23),tag=j1)
             end if
          end if
       end do

       ! Find index of this slab in global array
       j1 = islab23 + grid2%first_slab23(pub_my_proc_id) - 1
       ! As long as there is a slab left on this proc
       if (islab23<=grid2%num_slabs23) then
          ! And as long as we are still within the slabs of grid1
          if (j1<=(grid1%n1/2+1)) then
             ! Find which proc holds the slab required on this proc on this step
             recvproc = grid1%proc_slab23(j1)

             ! Receive slab, or copy if local
             if (recvproc/=pub_my_proc_id) then
                call comms_recv(recvproc,gbuf1(:,:),tag=j1)
             else
                jslab23 = j1 - grid1%first_slab23(pub_my_proc_id) + 1
                gbuf1(:,:) = gspc1(:,:,jslab23)
             end if

             gspc2(1:m3,1:m2,islab23)   = gbuf1(1:m3,1:m2)
             gspc2(k3:e3,1:m2,islab23)  = gbuf1(l3:n3,1:m2)
             gspc2(1:m3,k2:e2,islab23)  = gbuf1(1:m3,l2:n2)
             gspc2(k3:e3,k2:e2,islab23) = gbuf1(l3:n3,l2:n2)
          end if
       end if

    end do

    ! Transform from reciprocal to real space on grid 2
    call internal_fft3d_cell_backward(internal_cell_info(grid2%fft_index), &
         rspc2,gspc2)

    call timer_clock('fourier_interpolate_cell',2, &
         calculate_flops(internal_cell_info(grid1%fft_index)) &
         + calculate_flops(internal_cell_info(grid2%fft_index)))

    ! Free up comms resources
    call comms_free

    ! Dellocate workspaces
    deallocate(gbuf1,stat=ierr)
    call utils_dealloc_check('fourier_interpolate_cell','gbuf1',ierr)
    deallocate(gspc2,stat=ierr)
    call utils_dealloc_check('fourier_interpolate_cell','gspc2',ierr)
    deallocate(gspc1,stat=ierr)
    call utils_dealloc_check('fourier_interpolate_cell','gspc1',ierr)

  end subroutine fourier_interpolate_cell

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_filter_cell(rspc1,rspc2,grid1,grid2,apply_nyquist)
    !=========================================================================!
    ! Performs a Fast Fourier Transform Filtering from a specified grid to a  !
    ! smaller grid.                                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rspc1 (input):  Slab-distributed array that contains the real space   !
    !                   data on grid 1 to be filtered to grid 2.              !
    !   rspc2 (output): Slab-distributed array that will contain the result   !
    !                   of the filtering in real space on grid 2.             !
    !   grid1 (input):  GRID_INFO describing grid 1 (the larger grid).        !
    !   grid2 (input):  GRID_INFO describing grid 2 (the smaller grid).       !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 29/06/10                                      !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only : comms_free, comms_recv, comms_send, pub_my_proc_id
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid1
    type(GRID_INFO), intent(in) :: grid2
    real(kind=DP), intent(in) :: rspc1(grid1%ld1,grid1%ld2,grid1%max_slabs12)
    real(kind=DP), intent(out) :: rspc2(grid2%ld1,grid2%ld2,grid2%max_slabs12)
    logical, intent(in), optional :: apply_nyquist

    ! Local Variables
    integer :: ierr
    integer :: islab23,jslab23,kslab23
    integer :: n3,n2,m3,m2,e3,e2,l3,l2,k3,k2
    integer :: j1
    integer :: sendproc, recvproc
    complex(kind=DP), allocatable :: gspc1(:,:,:)
    complex(kind=DP), allocatable :: gspc2(:,:,:)
    complex(kind=DP), allocatable :: gbuf1(:,:)
    logical :: loc_apply_nyquist

    ! -----------------------------------------------------------------------

    call timer_clock('fourier_filter_cell',1)

    ! Optional arguments
    loc_apply_nyquist = .false.
    if (present(apply_nyquist)) loc_apply_nyquist = apply_nyquist

    ! Allocate workspaces
    allocate(gspc1(grid1%ld3,grid1%ld2,grid1%max_slabs23),stat=ierr)
    call utils_alloc_check('fourier_filter_cell','gspc1',ierr)
    allocate(gspc2(grid2%ld3,grid2%ld2,grid2%max_slabs23),stat=ierr)
    call utils_alloc_check('fourier_filter_cell','gspc2',ierr)
    allocate(gbuf1(grid1%ld3,grid1%ld2),stat=ierr)
    call utils_alloc_check('fourier_filter_cell','gbuf1',ierr)

    ! Make a license heartbeat call, if ACCELRYS defined
    call make_license_heartbeat

    ! Transform from real to reciprocal space on grid 1 (the larger grid)
    call internal_fft3d_cell_forward(internal_cell_info(grid1%fft_index), &
         rspc1,gspc1)

    ! Scale grid1 recip data
    gspc1 = gspc1 * (real(grid2%n1,kind=DP) * real(grid2%n2,kind=DP) &
         * real(grid2%n3,kind=DP)) / (real(grid1%n1,kind=DP) * &
         real(grid1%n2,kind=DP) * real(grid1%n3,kind=DP))

    ! Find limits in reciprocal i3,i2 directions on grids 1 and 2
    n3 = grid1%n3;         n2 = grid1%n2
    e3 = grid2%n3;         e2 = grid2%n2
    m3 = (e3+1)/2;         m2 = (e2+1)/2
    k3 = e3/2+1+mod(e3,2); k2 = e2/2+1+mod(e2,2)
    l3 = n3-e3+k3;         l2 = n2-e2+k2
    gspc2(:,:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)

    ! Loop over all slabs of grid2 (the smaller grid)
    do islab23=1,grid2%max_slabs23

       ! Find if there are any slabs needing to be sent to other procs
       ! on this step
       do jslab23=1,grid1%num_slabs23
          ! Find the global index of this slab, the proc holding this slab
          ! in grid2, and the local index of this slab on that proc
          j1 = jslab23 + grid1%first_slab23(pub_my_proc_id) - 1
          if (j1>(grid2%n1/2+1)) cycle
          sendproc = grid2%proc_slab23(j1)
          kslab23 = j1 - grid2%first_slab23(sendproc) + 1
          ! If the local index on that proc is the same as the step we
          ! are currently on, then send the slab
          if (kslab23==islab23) then
             if (sendproc/=pub_my_proc_id) then
                call comms_send(sendproc,gspc1(1,1,jslab23),grid1%ld3*grid1%ld2,tag=j1)
             end if
          end if
       end do

       ! Find which proc holds the slab required for this proc on this step
       j1 = islab23 + grid2%first_slab23(pub_my_proc_id) - 1
       if (islab23<=grid2%num_slabs23) then
          recvproc = grid1%proc_slab23(j1)

          jslab23 = j1 - grid1%first_slab23(recvproc) + 1
          if (recvproc/=pub_my_proc_id) then
             call comms_recv(recvproc,gbuf1(1,1),grid1%ld3*grid1%ld2,tag=j1)
          else
             gbuf1(:,:) = gspc1(:,:,jslab23)
          end if

          gspc2(1:m3,1:m2,islab23)   = gbuf1(1:m3,1:m2)
          gspc2(k3:e3,1:m2,islab23)  = gbuf1(l3:n3,1:m2)
          gspc2(1:m3,k2:e2,islab23)  = gbuf1(1:m3,l2:n2)
          gspc2(k3:e3,k2:e2,islab23) = gbuf1(l3:n3,l2:n2)

       end if

    end do

    ! Nyquist filter (grid is always going to be even)
    if (grid2%num_slabs23 > 0) then
       gspc2(grid2%n3/2+1,:,:) = (0.0_DP,0.0_DP)
       gspc2(:,grid2%n2/2+1,:) = (0.0_DP,0.0_DP)
    end if
    ! Nyquist filter for last slab
    if (pub_my_proc_id==grid2%proc_slab23(grid2%n1/2+1)) &
         gspc2(:,:,grid2%num_slabs23) = (0.0_DP,0.0_DP)

    ! Transform from reciprocal to real space on grid 2
    call internal_fft3d_cell_backward(internal_cell_info(grid2%fft_index), &
         rspc2,gspc2)

    call timer_clock('fourier_filter_cell',2, &
         calculate_flops(internal_cell_info(grid1%fft_index)) &
         + calculate_flops(internal_cell_info(grid2%fft_index)))

    ! Free up comms resources
    call comms_free

    ! Dellocate workspaces
    deallocate(gbuf1,stat=ierr)
    call utils_dealloc_check('fourier_filter_cell','gbuf1',ierr)
    deallocate(gspc2,stat=ierr)
    call utils_dealloc_check('fourier_filter_cell','gspc2',ierr)
    deallocate(gspc1,stat=ierr)
    call utils_dealloc_check('fourier_filter_cell','gspc1',ierr)

  end subroutine fourier_filter_cell

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_apply_box(grid,dir,rspc,gspc,scale_in,box,padbox,omp)

    use constants, only: DP
    use timer, only: timer_clock
    use fft_box, only: FFTBOX_INFO
    use utils, only: utils_assert, utils_abort
#ifdef ITC_TRACE
    use vt
#endif
#ifdef BASIC_GPU_FFT
    use fourier_gpu_wrapper_mod
#endif

    implicit none

    character, intent(in) :: grid
    character, intent(in) :: dir
    complex(kind=DP), dimension(:,:,:), intent(inout) :: rspc
    complex(kind=DP), dimension(:,:,:), intent(inout), optional :: gspc
    logical, optional, intent(in) :: scale_in
    ! Optional FFTBOX_INFO (for now)`
    type(FFTBOX_INFO), optional, intent(in) :: box
    ! Switch for padded box
    logical, optional, intent(in) :: padbox

    ! Switch for parallel transforms
    logical, optional, intent(in) :: omp

    ! Local variables
    integer :: npts ! Total number of points in FFT grid
    real(kind=DP) :: flops
    logical :: scale
    logical :: open_mp

#ifdef BASIC_GPU_FFT
    real(kind=DP) :: scalefac
#endif
    ! Check arguments

    ! Optional execution of data parallel transforms
    if (present(omp)) then
       open_mp = omp
    else
       open_mp = .false.
    end if

    ! Optional supression of scaling for backwards transform
    if (present(scale_in)) then
       scale = scale_in
    else
       scale = .true.
    end if

    if (grid /= 'C' .and. grid /= 'c' .and. grid /= 'F' .and. grid /= 'f') then
       call utils_abort('Error in fourier_apply_box: unknown grid flag "'//&
            trim(grid)//'"')
    end if

    if (dir /= 'F' .and. dir /= 'f' .and. dir /= 'B' .and. dir /= 'b') then
       call utils_abort('Error in fourier_apply_box: unknown direction flag "'//&
            trim(dir)//'"')
    end if

    if (present(gspc)) then

       call utils_assert(size(rspc,1) == size(gspc,1), 'Error in fourier_apply_&
            &box: array size mismatch in dimension 1')
       call utils_assert(size(rspc,2) == size(gspc,2), 'Error in fourier_apply_&
            &box: array size mismatch in dimension 2')
       call utils_assert(size(rspc,3) == size(gspc,3), 'Error in fourier_apply_&
            &box: array size mismatch in dimension 3')

    end if

    ! Start timer
    call timer_clock('fourier_apply_box',1)

#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_apply_box, vt_err)
#endif

    ! Copy input data into output array for in-place transform

    if (present(gspc)) then
       if (dir == 'F' .or. dir == 'f') then
          gspc(:,:,:) = rspc
       else
          rspc(:,:,:) = gspc
       end if
    end if

    ! Select appropriate routine
    if (present(box)) then
       if (present(gspc) .and. (dir == 'F' .or. dir == 'f')) then
          call internal_fft3d_box(internal_box_info(box%fft_index), &
               dir,gspc,scale,omp=open_mp)
       else
          call internal_fft3d_box(internal_box_info(box%fft_index), &
               dir,rspc,scale,omp=open_mp)
       end if
       npts = internal_box_info(box%fft_index)%n1 * &
              internal_box_info(box%fft_index)%n2 * &
              internal_box_info(box%fft_index)%n3
    else
       ! qoh: FFT box
       if (grid == 'C' .or. grid == 'c') then
          if (present(gspc) .and. (dir == 'F' .or. dir == 'f')) then
             call internal_fft3d_box(internal_box_coarse_info,dir,gspc,scale, &
                  padbox,omp=open_mp)
          else
#ifdef BASIC_GPU_FFT
             coarse_work_gpu = rspc     !kaw: This is the cuda Fortran syntax
                                        ! which will soon be common to the
                                        ! accelerator approach.
             if ((dir.eq.'f').or.(dir.eq.'F')) then
                call cufftExec(cufftplan_coarse,planType,coarse_work_gpu, &
                               coarse_work_gpu,CUFFT_FORWARD)
                rspc(:,:,:) = coarse_work_gpu
             else
                call cufftExec(cufftplan_coarse,planType,coarse_work_gpu, &
                               coarse_work_gpu,CUFFT_INVERSE)
                rspc(:,:,:) = coarse_work_gpu
                scalefac = 1.0_DP / (internal_box_coarse_info%n1* &
                           internal_box_coarse_info%n2*internal_box_coarse_info%n3)
                if (scale) then
                   rspc(:,:,:) = rspc * scalefac
                end if
             end if
#else
             call internal_fft3d_box(internal_box_coarse_info,dir,rspc,scale, &
                  padbox,omp=open_mp)
#endif
          end if
          npts = internal_box_coarse_info%n1 * internal_box_coarse_info%n2 * &
               internal_box_coarse_info%n3
       else
          if (present(gspc) .and. (dir == 'F' .or. dir == 'f')) then
             call internal_fft3d_box(internal_box_dbl_info,dir,gspc,scale, &
                  padbox,omp=open_mp)
          else
#ifdef BASIC_GPU_FFT
             fine_work_gpu = rspc
             if ((dir.eq.'f').or.(dir.eq.'F')) then
                call cufftExec(cufftplan_fine,planType,fine_work_gpu, &
                               fine_work_gpu,CUFFT_FORWARD)
                rspc(:,:,:) = fine_work_gpu
             else
                call cufftExec(cufftplan_fine,planType,fine_work_gpu, &
                               fine_work_gpu,CUFFT_INVERSE)
                rspc(:,:,:) = fine_work_gpu
                scalefac = 1.0_DP / (internal_box_dbl_info%n1* &
                           internal_box_dbl_info%n2*internal_box_dbl_info%n3)
                if (scale) then
                   rspc(:,:,:) = rspc * scalefac
                end if
             end if
#else
             call internal_fft3d_box(internal_box_dbl_info,dir,rspc,scale, &
                  padbox,omp=open_mp)
#endif
          end if
          npts = internal_box_dbl_info%n1 * internal_box_dbl_info%n2 * &
               internal_box_dbl_info%n3
       end if
    end if


    ! Stop timer
    flops = 5.0_DP * npts * log(real(npts,kind=DP)) / log(2.0_DP)
    call timer_clock('fourier_apply_box',2,flops)

#ifdef ITC_TRACE
    call VTEND(vt_fourier_apply_box, vt_err)
#endif

  end subroutine fourier_apply_box

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_apply_box_pair(grid,dir,rspc1,rspc2,gspc,omp)

    use constants, only: DP
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    character, intent(in) :: grid
    character, intent(in) :: dir
    real(kind=DP), dimension(:,:,:), intent(inout) :: rspc1, rspc2
    complex(kind=DP), dimension(:,:,:), intent(inout) :: gspc

    ! Switch for parallel transforms
    logical, optional, intent(in) :: omp

    ! Local variables

    integer :: npts ! Total number of points in FFT grid
    real(kind=DP) :: flops
    logical :: open_mp

    ! Check arguments

    if (grid /= 'C' .and. grid /= 'c' .and. grid /= 'F' .and. grid /= 'f') then
       call utils_abort('Error in fourier_apply_box_pair: unknown grid flag "'//&
            trim(grid)//'"')
    end if

    if (dir /= 'F' .and. dir /= 'f' .and. dir /= 'B' .and. dir /= 'b') then
       call utils_abort('Error in fourier_apply_box_pair: &
            &unknown direction flag "'//trim(dir)//'"')
    end if

    call utils_assert(size(rspc1,1) == size(gspc,1), 'Error in fourier_apply_&
         &box: array size mismatch in dimension 1')
    call utils_assert(size(rspc2,1) == size(gspc,1), 'Error in fourier_apply_&
         &box: array size mismatch in dimension 1')
    call utils_assert(size(rspc1,2) == size(gspc,2), 'Error in fourier_apply_&
         &box: array size mismatch in dimension 2')
    call utils_assert(size(rspc2,2) == size(gspc,2), 'Error in fourier_apply_&
         &box: array size mismatch in dimension 2')
    call utils_assert(size(rspc1,3) == size(gspc,3), 'Error in fourier_apply_&
         &box: array size mismatch in dimension 3')
    call utils_assert(size(rspc2,3) == size(gspc,3), 'Error in fourier_apply_&
         &box: array size mismatch in dimension 3')

    ! Optional execution of data parallel transforms
    if (present(omp)) then
       open_mp = omp
    else
       open_mp = .false.
    end if

    ! Start timer
    call timer_clock('fourier_apply_box_pair',1)
#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_apply_box_pair, vt_err)
#endif

    ! For forward routine, copy data into gspc array
    if (dir == 'F' .or. dir == 'f') gspc(:,:,:) = cmplx(rspc1,rspc2,kind=DP)

    ! Select appropriate routine

    if (grid == 'C' .or. grid == 'c') then
       call internal_fft3d_box(internal_box_coarse_info,dir,gspc,.true.,omp=open_mp)
       npts = internal_box_coarse_info%n1 * internal_box_coarse_info%n2 * &
            internal_box_coarse_info%n3
    else
       call internal_fft3d_box(internal_box_dbl_info,dir,gspc,.true.,omp=open_mp)
       npts = internal_box_dbl_info%n1 * internal_box_dbl_info%n2 * &
            internal_box_dbl_info%n3
    end if

    ! For backward routine, copy data out of gspc array
    if (dir == 'B' .or. dir == 'b') then
       rspc1(:,:,:) = real(gspc,kind=DP)
       rspc2(:,:,:) = aimag(gspc)
    end if

    ! Stop timer
    flops = 5.0_DP * npts * log(real(npts,kind=DP)) / log(2.0_DP)
    call timer_clock('fourier_apply_box_pair',2,flops)

#ifdef ITC_TRACE
    call VTEND(vt_fourier_apply_box_pair, vt_err)
#endif

  end subroutine fourier_apply_box_pair

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_exit_fftbox

#ifdef GPU_PGI
    use fourier_gpu_wrapper_mod
#endif
    implicit none

    ! Double FFTbox grid
    call internal_serial_exit(internal_box_dbl_info)

    ! Coarse FFTbox grid
    call internal_serial_exit(internal_box_coarse_info)
#ifdef GPU_PGI
       ! kaw: Destroy the cufft plans for fourier_interpolate_product and
       !      deallocate the gpu work arrays.
       call fourier_gpu_exit
#endif

  end subroutine fourier_exit_fftbox

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_exit_box(box)

    use fft_box, only: FFTBOX_INFO

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(inout) :: box

    call internal_serial_exit(internal_box_info(box%fft_index))
    box%fft_index = 0
    box_info_index = box_info_index - 1

  end subroutine fourier_exit_box

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_exit_cell

    implicit none

    ! Local variables
    integer :: igrid ! whole-cell grid counter

    ! Whole-simulation-cell grids
    do igrid=cell_info_index,1,-1
       call internal_parallel_exit(internal_cell_info(igrid))
    end do
    cell_info_index = 0

  end subroutine fourier_exit_cell

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_interpolate_2real(coarse_work,fine_work,in1,in2,out1,out2)

    use constants, only: DP
    use rundat, only: pub_threads_fftbox
!$  use rundat, only: pub_threads_per_fftbox
    use utils, only: utils_assert
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout), dimension(:,:,:) :: coarse_work
    complex(kind=DP), intent(inout), dimension(:,:,:) :: fine_work
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out1(:,:,:)
    real(kind=DP), intent(out) :: out2(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3,it
    real(kind=DP) :: scale

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_interpolate, vt_err)
#endif

    ! Check array sizes
    call utils_assert(.not. &
         ((size(in1,1) < n1 .or. size(in1,2) < n2 .or. size(in1,3) < n3)), &
         'Error in fourier_interpolate: invalid dimensions for array in1')
    call utils_assert(.not. &
         ((size(in2,1) < n1 .or. size(in2,2) < n2 .or. size(in2,3) < n3)), &
         'Error in fourier_interpolate: invalid dimensions for array in2')

    call utils_assert(.not. (size(out1,1) < 2*n1 .or. size(out1,2) < 2*n2 .or. &
         size(out1,3) < 2*n3), 'Error in fourier_interpolate: &
            &invalid dimensions for array out1')
    call utils_assert(.not. (size(out2,1) < 2*n1 .or. size(out2,2) < 2*n2 .or. &
         size(out2,3) < 2*n3), 'Error in fourier_interpolate: &
            &invalid dimensions for array out2')

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scale = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)



!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(it)
!$OMP DO
    do it = 1, n3
       coarse_work(1:n1,1:n2,it) = scale * cmplx(in1(1:n1,1:n2,it),in2(1:n1,1:n2,it),kind=DP)
    end do
!$OMP END DO
!$OMP END PARALLEL


    ! Fourier transform to reciprocal space on coarse grid
    call fourier_apply_box('Coarse','Forward',coarse_work, &
         scale_in=.false., omp=pub_threads_fftbox)


    ! Zero all other components
!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(it)
!$OMP DO
    do it = 1,2*n3
       fine_work(1:2*n1,1:2*n2,it) = 0.0_DP
    end do
!$OMP END DO


    ! Copy Fourier components to fine reciprocal space grid
!$OMP DO
    do it = 1, m3
       fine_work(1:m1,1:m2,it) = coarse_work(1:m1,1:m2,it)
       fine_work(k1:2*n1,1:m2,it) = coarse_work(l1:n1,1:m2,it)
       fine_work(1:m1,k2:2*n2,it) = coarse_work(1:m1,l2:n2,it)
       fine_work(k1:2*n1,k2:2*n2,it) = coarse_work(l1:n1,l2:n2,it)
    end do
!$OMP DO
    do it = l3, n3
       fine_work(1:m1,1:m2,it+n3) = coarse_work(1:m1,1:m2,it)
       fine_work(k1:2*n1,1:m2,it+n3) = coarse_work(l1:n1,1:m2,it)
       fine_work(1:m1,k2:2*n2,it+n3) = coarse_work(1:m1,l2:n2,it)
       fine_work(k1:2*n1,k2:2*n2,it+n3) = coarse_work(l1:n1,l2:n2,it)
    end do
!$OMP END PARALLEL


    ! Fourier transform to real space on fine grid
    call fourier_apply_box('Fine','Backward',fine_work,scale_in=.false., &
         padbox=.true., omp=pub_threads_fftbox)


    ! Unpack data from real and imaginary parts of complex array
    ! ndmh: a loop saves on cache misses when repeating for aimag
!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(it)
!$OMP DO
    do it = 1,2*n3
       out1(1:2*n1,1:2*n2,it) = real(fine_work(1:2*n1,1:2*n2,it),DP)
       out2(1:2*n1,1:2*n2,it) = aimag(fine_work(1:2*n1,1:2*n2,it))
    end do
!$OMP END DO
!$OMP END PARALLEL


#ifdef ITC_TRACE
    call VTEND(vt_fourier_interpolate, vt_err)
#endif

  end subroutine fourier_interpolate_2real

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_interpolate_1complex(coarse_work,fine_work,in)

    use constants, only: DP, stdout
    use rundat, only: pub_threads_fftbox, pub_threads_per_fftbox
    use utils, only: utils_assert
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout), dimension(:,:,:) :: coarse_work
    complex(kind=DP), intent(inout), dimension(:,:,:) :: fine_work
    complex(kind=DP), intent(in) :: in(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3,it
    real(kind=DP) :: scale

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_interpolate, vt_err)
#endif

    ! Check array sizes
    call utils_assert(.not. &
         ((size(in,1) < n1 .or. size(in,2) < n2 .or. size(in,3) < n3)), &
         'Error in fourier_interpolate: invalid dimensions for array in')

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scale = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)

    coarse_work = scale * in

    ! Fourier transform to reciprocal space on coarse grid
    call fourier_apply_box('Coarse','Forward',coarse_work, &
         scale_in=.false., omp=pub_threads_fftbox)

    ! Zero all other components
!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(it)
!$OMP DO
    do it = 1,2*n3
       fine_work(1:2*n1,1:2*n2,it) = 0.0_DP
    end do
!$OMP END DO


    ! Copy Fourier components to fine reciprocal space grid
!$OMP DO
    do it = 1, m3
       fine_work(1:m1,1:m2,it) = coarse_work(1:m1,1:m2,it)
       fine_work(k1:2*n1,1:m2,it) = coarse_work(l1:n1,1:m2,it)
       fine_work(1:m1,k2:2*n2,it) = coarse_work(1:m1,l2:n2,it)
       fine_work(k1:2*n1,k2:2*n2,it) = coarse_work(l1:n1,l2:n2,it)
    end do
!$OMP END DO
!$OMP DO
    do it = l3, n3
       fine_work(1:m1,1:m2,it+n3) = coarse_work(1:m1,1:m2,it)
       fine_work(k1:2*n1,1:m2,it+n3) = coarse_work(l1:n1,1:m2,it)
       fine_work(1:m1,k2:2*n2,it+n3) = coarse_work(1:m1,l2:n2,it)
       fine_work(k1:2*n1,k2:2*n2,it+n3) = coarse_work(l1:n1,l2:n2,it)
    end do
!$OMP END DO
!$OMP END PARALLEL


    ! Fourier transform to real space on fine grid
    call fourier_apply_box('Fine','Backward',fine_work,scale_in=.false., &
         padbox=.true., omp=pub_threads_fftbox)

#ifdef ITC_TRACE
    call VTEND(vt_fourier_interpolate, vt_err)
#endif

  end subroutine fourier_interpolate_1complex

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine fourier_interpolate_product(coarse_work,fine_work,in1,in2,out)

    use constants, only: DP
    use rundat, only: pub_threads_fftbox
!$  use rundat, only: pub_threads_per_fftbox
    use utils, only: utils_assert, utils_abort
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout), dimension(:,:,:) :: coarse_work
    complex(kind=DP), intent(inout), dimension(:,:,:) :: fine_work
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3
    real(kind=DP) :: scale


#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_interpolate_product, vt_err)
#endif

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3
    ! Coarse Grid

    ! Check array sizes
    call utils_assert(&
         size(in1,1) >= n1 .and. size(in1,2) >= n2 .and. size(in1,3) >= n3, &
        'Error in fourier_interpolate_product: invalid dimensions for array in1')
    call utils_assert(&
         size(in2,1) >= n1 .and. size(in2,2) >= n2 .and. size(in2,3) >= n3, &
        'Error in fourier_interpolate_product: invalid dimensions for array in2')

    if (size(out,1) < 2*n1 .or. size(out,2) < 2*n2 .or. &
         size(out,3) < 2*n3) then
       call utils_abort('Error in fourier_interpolate_product: &
            &invalid dimensions for array out')
    end if

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scale = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)

!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(j3)
!$OMP DO
    do j3 = 1, n3
       coarse_work(1:n1,1:n2,j3) = &
            scale * cmplx(in1(1:n1,1:n2,j3),in2(1:n1,1:n2,j3),kind=DP)
    end do
!$OMP END DO
!$OMP END PARALLEL

    call fourier_apply_box('Coarse', 'Forward', coarse_work, &
         scale_in=.false., omp=pub_threads_fftbox)

!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(j3)
!$OMP DO
    do j3 = 1,2*n3
       fine_work(1:2*n1,1:2*n2,j3) = 0.0_DP
    end do
!$OMP END DO

    ! Copy Fourier components to fine reciprocal space grid
!$OMP DO
    do j3 = 1, m3
       fine_work(1:m1,1:m2,j3) = coarse_work(1:m1,1:m2,j3)
       fine_work(k1:2*n1,1:m2,j3) = coarse_work(l1:n1,1:m2,j3)
       fine_work(1:m1,k2:2*n2,j3) = coarse_work(1:m1,l2:n2,j3)
       fine_work(k1:2*n1,k2:2*n2,j3) = coarse_work(l1:n1,l2:n2,j3)
    end do
!$OMP END DO
!$OMP DO
    do j3 = l3, n3
       fine_work(1:m1,1:m2,j3+n3) = coarse_work(1:m1,1:m2,j3)
       fine_work(k1:2*n1,1:m2,j3+n3) = coarse_work(l1:n1,1:m2,j3)
       fine_work(1:m1,k2:2*n2,j3+n3) = coarse_work(1:m1,l2:n2,j3)
       fine_work(k1:2*n1,k2:2*n2,j3+n3) = coarse_work(l1:n1,l2:n2,j3)
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! Fourier transform to real space on fine grid
    call fourier_apply_box('Fine','Backward',fine_work, &
         scale_in=.false., padbox=.true., omp=pub_threads_fftbox)

    ! ndmh: unpack data from real and imaginary parts of complex array
    ! ndmh: and multiply them together to find product
!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(j3)
!$OMP DO
    do j3 = 1,2*n3
       out(1:n1*2,1:n2*2,j3) = real(fine_work(1:n1*2,1:n2*2,j3),DP) &
            * aimag(fine_work(1:n1*2,1:n2*2,j3))
    end do
!$OMP END PARALLEL

#ifdef ITC_TRACE
    call VTEND(vt_fourier_interpolate_product, vt_err)
#endif

  end subroutine fourier_interpolate_product


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  subroutine fourier_filter_2real(coarse_work,fine_work,in1,in2,out1,out2)

    use constants, only: DP
    use rundat, only: pub_threads_fftbox
!$  use rundat, only: pub_threads_per_fftbox
    use utils, only: utils_assert
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout), dimension(:,:,:) :: coarse_work
    complex(kind=DP), intent(inout), dimension(:,:,:) :: fine_work
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out1(:,:,:)
    real(kind=DP), intent(out) :: out2(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j3
    real(kind=DP) :: scale

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_filter, vt_err)
#endif

    ! Check array sizes
    call utils_assert(&
         size(in1,1) >= 2*n1 .and. size(in1,2) >= 2*n2 .and. size(in1,3) >= 2*n3, &
        'Error in fourier_filter: invalid dimensions for array in1')
    call utils_assert(&
         size(in2,1) >= 2*n1 .and. size(in2,2) >= 2*n2 .and. size(in2,3) >= 2*n3, &
        'Error in fourier_filter: invalid dimensions for array in2')
    call utils_assert(&
         size(out1,1) >= n1 .and. size(out1,2) >= n2 .and. size(out1,3) >= n3, &
        'Error in fourier_filter: invalid dimensions for array out1')
    call utils_assert(&
         size(out2,1) >= n1 .and. size(out2,2) >= n2 .and. size(out2,3) >= n3, &
        'Error in fourier_filter: invalid dimensions for array out2')

    m1 = (n1+1)/2 ;     m2 = (n2+1)/2 ;         m3 = (n3+1)/2
    l1 = n1/2+2 ;       l2 = n2/2+2 ;           l3 = n3/2+2
    k1 = n1 + l1 ;      k2 = n2 + l2 ;          k3 = n3 + l3

    ! Pack data into real and imaginary parts of complex array
!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(j3)
!$OMP DO
    do j3 = 1, 2*n3
       fine_work(1:2*n1,1:2*n2,j3) = cmplx(in1(1:2*n1,1:2*n2,j3), &
            in2(1:2*n1,1:2*n2,j3),dp)
    end do
!$OMP END DO
!$OMP END PARALLEL


    ! Fourier transform to reciprocal space on fine grid
    call fourier_apply_box('Fine','Forward',fine_work,scale_in=.false., &
         padbox=.true., omp=pub_threads_fftbox)


    ! Copy Fourier components to coarse reciprocal space grid
!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(j3)
!$OMP DO
    do j3 = 1, m3
       coarse_work(1:m1,1:m2,j3) = fine_work(1:m1,1:m2,j3)
       coarse_work(l1:n1,1:m2,j3) = fine_work(k1:2*n1,1:m2,j3)
       coarse_work(1:m1,l2:n2,j3) = fine_work(1:m1,k2:2*n2,j3)
       coarse_work(l1:n1,l2:n2,j3) = fine_work(k1:2*n1,k2:2*n2,j3)
    end do
!$OMP DO
    do j3 = l3, n3
       coarse_work(1:m1,1:m2,j3) = fine_work(1:m1,1:m2,j3+n3)
       coarse_work(l1:n1,1:m2,j3) = fine_work(k1:2*n1,1:m2,j3+n3)
       coarse_work(1:m1,l2:n2,j3) = fine_work(1:m1,k2:2*n2,j3+n3)
       coarse_work(l1:n1,l2:n2,j3) = fine_work(k1:2*n1,k2:2*n2,j3+n3)
    end do
!$OMP END PARALLEL

    ! Fourier transform to real space on coarse grid
    call fourier_apply_box('Coarse','Backward',coarse_work, &
         scale_in=.false., omp=pub_threads_fftbox)


    ! Normalise
    scale = 0.125_DP / (internal_box_coarse_info%ld1*&
         internal_box_coarse_info%ld2*internal_box_coarse_info%n3)

!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(j3)
!$OMP DO
    do j3 = 1, n3
       coarse_work(1:n1,1:n2,j3) = scale * coarse_work(1:n1,1:n2,j3)
    end do

    ! Unpack data from real and imaginary parts of complex array
    ! ndmh: a loop saves on cache misses when repeating for aimag
!$OMP DO
    do j3 = 1,n3
       out1(1:n1,1:n2,j3) =  real(coarse_work(1:n1,1:n2,j3),DP)
       out2(1:n1,1:n2,j3) = aimag(coarse_work(1:n1,1:n2,j3))
    end do
!$OMP END PARALLEL


#ifdef ITC_TRACE
    call VTEND(vt_fourier_filter, vt_err)
#endif


  end subroutine fourier_filter_2real


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  subroutine fourier_filter_1complex(coarse_work,fine_work,in,out)

    use constants, only: DP, stdout
    use rundat, only: pub_threads_fftbox, pub_threads_per_fftbox
    use utils, only: utils_assert
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout), dimension(:,:,:) :: coarse_work
    complex(kind=DP), intent(inout), dimension(:,:,:) :: fine_work
    complex(kind=DP), intent(in) :: in(:,:,:)
    complex(kind=DP), intent(out) :: out(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j3
    real(kind=DP) :: scale

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_filter, vt_err)
#endif

    ! Check array sizes
    call utils_assert(&
         size(in,1) >= 2*n1 .and. size(in,2) >= 2*n2 .and. size(in,3) >= 2*n3, &
        'Error in fourier_filter: invalid dimensions for array in')
    call utils_assert(&
         size(out,1) >= n1 .and. size(out,2) >= n2 .and. size(out,3) >= n3, &
        'Error in fourier_filter: invalid dimensions for array out')

    m1 = (n1+1)/2 ;     m2 = (n2+1)/2 ;         m3 = (n3+1)/2
    l1 = n1/2+2 ;       l2 = n2/2+2 ;           l3 = n3/2+2
    k1 = n1 + l1 ;      k2 = n2 + l2 ;          k3 = n3 + l3

    fine_work = in

    ! Fourier transform to reciprocal space on fine grid
    call fourier_apply_box('Fine','Forward',fine_work,scale_in=.false., &
         padbox=.true., omp=pub_threads_fftbox)

    ! Copy Fourier components to coarse reciprocal space grid
!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(j3)
!$OMP DO
    do j3 = 1, m3
       coarse_work(1:m1,1:m2,j3) = fine_work(1:m1,1:m2,j3)
       coarse_work(l1:n1,1:m2,j3) = fine_work(k1:2*n1,1:m2,j3)
       coarse_work(1:m1,l2:n2,j3) = fine_work(1:m1,k2:2*n2,j3)
       coarse_work(l1:n1,l2:n2,j3) = fine_work(k1:2*n1,k2:2*n2,j3)
    end do
!$OMP DO
    do j3 = l3, n3
       coarse_work(1:m1,1:m2,j3) = fine_work(1:m1,1:m2,j3+n3)
       coarse_work(l1:n1,1:m2,j3) = fine_work(k1:2*n1,1:m2,j3+n3)
       coarse_work(1:m1,l2:n2,j3) = fine_work(1:m1,k2:2*n2,j3+n3)
       coarse_work(l1:n1,l2:n2,j3) = fine_work(k1:2*n1,k2:2*n2,j3+n3)
    end do
!$OMP END PARALLEL

    ! Fourier transform to real space on coarse grid
    call fourier_apply_box('Coarse','Backward',coarse_work, &
         scale_in=.false., omp=pub_threads_fftbox)

    ! Normalise
    scale = 0.125_DP / (internal_box_coarse_info%ld1*&
         internal_box_coarse_info%ld2*internal_box_coarse_info%n3)

!$OMP PARALLEL NUM_THREADS(pub_threads_per_fftbox) &
!$OMP PRIVATE(j3)
!$OMP DO
    do j3 = 1, n3
       out(1:n1,1:n2,j3) = scale * coarse_work(1:n1,1:n2,j3)
    end do
!$OMP END DO
!$OMP END PARALLEL


#ifdef ITC_TRACE
    call VTEND(vt_fourier_filter, vt_err)
#endif


  end subroutine fourier_filter_1complex


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#if 0
  subroutine fourier_interpolate(in1,in2,out1,out2)

    use constants, only: DP
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out1(:,:,:)
    real(kind=DP), intent(out) :: out2(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3
    real(kind=DP) :: scale

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_interpolate, vt_err)
#endif

    ! Check array sizes
    call utils_assert(&
         size(in1,1) >= n1 .and. size(in1,2) >= n2 .and. size(in1,3) >= n3),&
        'Error in fourier_interpolate: invalid dimensions for array in1')
    call utils_assert(&
         size(in2,1) >= n1 .and. size(in2,2) >= n2 .and. size(in2,3) >= n3),&
        'Error in fourier_interpolate: invalid dimensions for array in2')
    call utils_assert(size(out1,1) >= 2*n1 .and. size(out1,2) >= 2*n2 .and. &
         size(out1,3) >= 2*n3,'Error in fourier_interpolate: &
            &invalid dimensions for array out1')
    call utils_assert(size(out2,1) >= 2*n1 .and. size(out2,2) >= 2*n2 .and. &
         size(out2,3) >= 2*n3,'Error in fourier_interpolate: &
            &invalid dimensions for array out2')

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scale = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)
    coarse_work(1:n1,1:n2,1:n3) = scale * cmplx(in1(1:n1,1:n2,1:n3), &
         in2(1:n1,1:n2,1:n3),kind=DP)

    ! Fourier transform to reciprocal space on coarse grid

    call fourier_apply_box('Coarse','Forward',coarse_work,scale_in=.false.)

    ! Copy Fourier components to fine reciprocal space grid

    fine_work(1:m1,1:m2,1:m3) = coarse_work(1:m1,1:m2,1:m3)
    fine_work(k1:2*n1,1:m2,1:m3) = coarse_work(l1:n1,1:m2,1:m3)
    fine_work(1:m1,k2:2*n2,1:m3) = coarse_work(1:m1,l2:n2,1:m3)
    fine_work(k1:2*n1,k2:2*n2,1:m3) = coarse_work(l1:n1,l2:n2,1:m3)
    fine_work(1:m1,1:m2,k3:2*n3) = coarse_work(1:m1,1:m2,l3:n3)
    fine_work(k1:2*n1,1:m2,k3:2*n3) = coarse_work(l1:n1,1:m2,l3:n3)
    fine_work(1:m1,k2:2*n2,k3:2*n3) = coarse_work(1:m1,l2:n2,l3:n3)
    fine_work(k1:2*n1,k2:2*n2,k3:2*n3) = coarse_work(l1:n1,l2:n2,l3:n3)

    ! Zero all other components

    fine_work(:,:,l3:j3) = 0.0_DP
    fine_work(:,l2:j2,1:l3-1) = 0.0_DP
    fine_work(:,l2:j2,j3+1:) = 0.0_DP
    fine_work(l1:j1,1:l2-1,1:l3-1) = 0.0_DP
    fine_work(l1:j1,j2+1:,1:l3-1) = 0.0_DP
    fine_work(l1:j1,1:l2-1,j3+1:) = 0.0_DP
    fine_work(l1:j1,j2+1:,j3+1:) = 0.0_DP

    ! Fourier transform to real space on fine grid

    call fourier_apply_box('Fine','Backward',fine_work,scale_in=.false., &
         padbox=.true.)

    ! Unpack data from real and imaginary parts of complex array
    ! ndmh: a loop saves on cache misses when repeating for aimag
    do j3 = 1,2*n3
       do j2 = 1,2*n2
          out1(1:2*n1,j2,j3) = real(fine_work(1:2*n1,j2,j3),DP)
          out2(1:2*n1,j2,j3) = aimag(fine_work(1:2*n1,j2,j3))
       end do
    end do

#ifdef ITC_TRACE
    call VTEND(vt_fourier_interpolate, vt_err)
#endif

  end subroutine fourier_interpolate

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_interpolate_product(in1,in2,out)

    use constants, only: DP
    use utils, only: utils_assert
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3
    real(kind=DP) :: scale

#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_interpolate_product, vt_err)
#endif

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

    ! Check array sizes
    call utils_assert(&
         size(in1,1) >= n1 .and. size(in1,2) >= n2 .and. size(in1,3) >= n3),&
        'Error in fourier_interpolate_product: invalid dimensions for array in1')
    call utils_assert(&
         size(in2,1) >= n2 .and. size(in2,2) >= n2 .and. size(in2,3) >= n3),&
        'Error in fourier_interpolate_product: invalid dimensions for array in2')
    call utils_assert(size(out,1) >= 2*n1 .and. size(out,2) >= 2*n2 .and. &
         size(out1,3) >= 2*n3,'Error in fourier_interpolate: &
            &invalid dimensions for array out')

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scale = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)
    coarse_work(1:n1,1:n2,1:n3) = scale * cmplx(in1(1:n1,1:n2,1:n3), &
         in2(1:n1,1:n2,1:n3),kind=DP)

    ! Fourier transform to reciprocal space on coarse grid

    call fourier_apply_box('Coarse','Forward',coarse_work,scale_in=.false.)

    ! Copy Fourier components to fine reciprocal space grid

    fine_work(1:m1,1:m2,1:m3) = coarse_work(1:m1,1:m2,1:m3)
    fine_work(k1:2*n1,1:m2,1:m3) = coarse_work(l1:n1,1:m2,1:m3)
    fine_work(1:m1,k2:2*n2,1:m3) = coarse_work(1:m1,l2:n2,1:m3)
    fine_work(k1:2*n1,k2:2*n2,1:m3) = coarse_work(l1:n1,l2:n2,1:m3)
    fine_work(1:m1,1:m2,k3:2*n3) = coarse_work(1:m1,1:m2,l3:n3)
    fine_work(k1:2*n1,1:m2,k3:2*n3) = coarse_work(l1:n1,1:m2,l3:n3)
    fine_work(1:m1,k2:2*n2,k3:2*n3) = coarse_work(1:m1,l2:n2,l3:n3)
    fine_work(k1:2*n1,k2:2*n2,k3:2*n3) = coarse_work(l1:n1,l2:n2,l3:n3)

    ! Zero all other components

    fine_work(:,:,l3:j3) = 0.0_DP
    fine_work(:,l2:j2,1:l3-1) = 0.0_DP
    fine_work(:,l2:j2,j3+1:) = 0.0_DP
    fine_work(l1:j1,1:l2-1,1:l3-1) = 0.0_DP
    fine_work(l1:j1,j2+1:,1:l3-1) = 0.0_DP
    fine_work(l1:j1,1:l2-1,j3+1:) = 0.0_DP
    fine_work(l1:j1,j2+1:,j3+1:) = 0.0_DP

    ! Fourier transform to real space on fine grid
    call fourier_apply_box('Fine','Backward',fine_work,scale_in=.false., &
         padbox=.true.)

    ! ndmh: unpack data from real and imaginary parts of complex array
    ! ndmh: and multiply them together to find product
    do j3 = 1,2*n3
       do j2 = 1,2*n2
          out(1:2*n1,j2,j3) = real(fine_work(1:2*n1,j2,j3),DP) * &
               aimag(fine_work(1:2*n1,j2,j3))
       end do
    end do

#ifdef ITC_TRACE
    call VTEND(vt_fourier_interpolate_product, vt_err)
#endif

  end subroutine fourier_interpolate_product

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_filter(in1,in2,out1,out2)

    use constants, only: DP
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out1(:,:,:)
    real(kind=DP), intent(out) :: out2(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3
    real(kind=DP) :: scale

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

#ifdef ITC_TRACE
    call VTBEGIN(vt_fourier_filter, vt_err)
#endif

    ! Check array sizes
    call utils_assert(&
         size(in1,1) >= 2*n1 .and. size(in1,2) >= 2*n2 .and. size(in1,3) >= 2*n3),&
        'Error in fourier_filter: invalid dimensions for array in1')
    call utils_assert(&
         size(in2,1) >= 2*n1 .and. size(in2,2) >= 2*n2 .and. size(in2,3) >= 2*n3),&
        'Error in fourier_filter: invalid dimensions for array in2')
    call utils_assert(size(out1,1) >= n1 .and. size(out1,2) >= n2 .and. &
         size(out1,3) >= n3,'Error in fourier_filter: &
            &invalid dimensions for array out1')
    call utils_assert(size(out2,1) >= n1 .and. size(out2,2) >= n2 .and. &
         size(out2,3) >= n3,'Error in fourier_filter: &
            &invalid dimensions for array out2')

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1 + l1 ; k2 = n2 + l2 ; k3 = n3 + l3

    ! Pack data into real and imaginary parts of complex array

    fine_work(1:2*n1,1:2*n2,1:2*n3) = cmplx(in1(1:2*n1,1:2*n2,1:2*n3), &
         in2(1:2*n1,1:2*n2,1:2*n3),dp)

    ! Fourier transform to reciprocal space on fine grid

    call fourier_apply_box('Fine','Forward',fine_work,scale_in=.false., &
         padbox=.true.)

    ! Copy Fourier components to coarse reciprocal space grid

    coarse_work(1:m1,1:m2,1:m3) = fine_work(1:m1,1:m2,1:m3)
    coarse_work(l1:n1,1:m2,1:m3) = fine_work(k1:2*n1,1:m2,1:m3)
    coarse_work(1:m1,l2:n2,1:m3) = fine_work(1:m1,k2:2*n2,1:m3)
    coarse_work(l1:n1,l2:n2,1:m3) = fine_work(k1:2*n1,k2:2*n2,1:m3)
    coarse_work(1:m1,1:m2,l3:n3) = fine_work(1:m1,1:m2,k3:2*n3)
    coarse_work(l1:n1,1:m2,l3:n3) = fine_work(k1:2*n1,1:m2,k3:2*n3)
    coarse_work(1:m1,l2:n2,l3:n3) = fine_work(1:m1,k2:2*n2,k3:2*n3)
    coarse_work(l1:n1,l2:n2,l3:n3) = fine_work(k1:2*n1,k2:2*n2,k3:2*n3)

    ! Fourier transform to real space on coarse grid

    call fourier_apply_box('Coarse','Backward',coarse_work,scale_in=.false.)

    ! Normalise
    scale = 0.125_DP / (internal_box_coarse_info%ld1*&
         internal_box_coarse_info%ld2*internal_box_coarse_info%n3)
    coarse_work(1:n1,1:n2,1:n3) = scale * coarse_work(1:n1,1:n2,1:n3)

    ! Unpack data from real and imaginary parts of complex array
    ! ndmh: a loop saves on cache misses when repeating for aimag
    do k3 = 1,n3
       do k2 = 1,n2
          do k1 = 1,n1
             out1(k1,k2,k3) = real(coarse_work(k1,k2,k3),dp)
             out2(k1,k2,k3) = aimag(coarse_work(k1,k2,k3))
          end do
       end do
    end do

    !out1(1:n1,1:n2,1:n3) = real(coarse_work(1:n1,1:n2,1:n3),dp)
    !out2(1:n1,1:n2,1:n3) = aimag(coarse_work(1:n1,1:n2,1:n3))

#ifdef ITC_TRACE
    call VTEND(vt_fourier_filter, vt_err)
#endif


  end subroutine fourier_filter
#endif

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ! *****************************************
  ! ***  P r i v a t e   r o u t i n e s  ***
  ! *****************************************

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_serial_init(n1,n2,n3,ld1,ld2,fft3d_info)

#if defined(SUN) || defined(ACML)
    use constants, only: DP
#endif
#ifdef MKL_FFTW3
    use constants, only: I4B
#endif
    use rundat, only: pub_threads_fftbox, pub_threads_per_fftbox, &
         pub_threads_num_fftboxes
    use utils, only : utils_alloc_check, utils_assert, utils_use_var

    implicit none

    integer, intent(in) :: n1,n2,n3 ! Dimensions of FFT grid
    integer, intent(in) :: ld1,ld2  ! Array dimensions containing grid
    type(serial_fft3d_info), intent(out) :: fft3d_info

    ! External subroutines
#ifdef FFTW
    external ::  fftw3d_f77_create_plan
#endif
#ifdef FFTW3
    external :: dfftw_plan_many_dft
#ifndef FFTW3_NO_OMP
    external :: dfftw_plan_with_nthreads
#endif
#endif
#ifdef ACML
    external :: zfft3dy
#endif
#ifdef VENDOR
#ifdef SUN
    external :: zfftz3
#endif
#endif

    ! Local variables
    integer :: ierr                 ! Error flag
#ifdef ACML
    integer :: tid
#endif
#ifdef SUN
    complex(kind=DP) :: zdum
#endif
!$  integer, external :: omp_get_thread_num

#ifdef MKL_FFTW3
    ! ndmh: code to ensure that when running MKL FFTW wrappers, the
    ! ndmh: plans are designed to be shared between multiple threads
#ifndef F2008
    !dec$ attributes align : 8 :: fftw3_mkl
    ! jd: ^ the above is not F2008-compliant.
    ! https://community.intel.com/t5/Intel-oneAPI-Math-Kernel-Library/
    ! /Calls-using-FFTW-wrapper-don-t-seem-to-be-thread-safe/td-p/945735
    ! claims that this is only needed to shut up the warning and that
    ! everything should work regardless.
#endif
    common/fftw3_mkl/ignore(4),number_of_user_threads,ignore2(7)
    integer(kind=I4B) :: ignore, number_of_user_threads, ignore2
    bind (c) :: /fftw3_mkl/

    number_of_user_threads = pub_threads_num_fftboxes
#else
    call utils_use_var(pub_threads_num_fftboxes)
#endif

    ! Simple checks on validity of input arguments

    call utils_assert(n1 >=1, 'Error in internal_serial_init &
            &(fourier_mod.F90): n1 < 1')
    call utils_assert(n2 >=1, 'Error in internal_serial_init &
            &(fourier_mod.F90): n2 < 1')
    call utils_assert(n3 >=1, 'Error in internal_serial_init &
            &(fourier_mod.F90): n3 < 1')

    call utils_assert(ld1 >= n1, 'Error in internal_serial_init &
            &(fourier_mod.F90): ld1 < n1')
    call utils_assert(ld2 >= n2, 'Error in internal_serial_init &
            &(fourier_mod.F90): ld2 < n2')

    ! Platform-dependent checking

#ifdef FFTW
    call utils_assert(ld1 == n1,'Error in internal_serial_init &
            &(fourier_mod.F90): ld1 must equal n1 for FFTw')
    call utils_assert(ld2 == n2,'Error in internal_serial_init &
            &(fourier_mod.F90): ld2 must equal n2 for FFTw')
#endif

    ! Simple initialisation

    ! Store basic information in fft3d_info structure

    fft3d_info%n1 = n1
    fft3d_info%n2 = n2
    fft3d_info%n3 = n3
    fft3d_info%ld1 = ld1
    fft3d_info%ld2 = ld2
    fft3d_info%iwork_len = 0
    fft3d_info%dwork_len = 0
    fft3d_info%zwork_len = 0

    ! Platform-dependent workspace allocation

#ifdef FFTW3
    fft3d_info%iwork_len = 12
    fft3d_info%zwork_len = ld1*ld2*n3
#endif

#ifdef ACML
    fft3d_info%zwork_len = ld1*ld2*n3
    fft3d_info%tbl_size = ld1*ld2*n3 + 4*(ld1+ld2+n3) + 400
    allocate(fft3d_info%tbl(fft3d_info%tbl_size, &
         0:pub_threads_num_fftboxes-1),stat=ierr)
    call utils_alloc_check('internal_serial_init (fourier_mod.F90)', &
         'fft3d_info%tbl',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
    fft3d_info%iwork_len = dxml_structure_size
#endif
#ifdef SUN
    fft3d_info%iwork_len = 3*128
    fft3d_info%dwork_len = 2*(n1+n2+n3)
    fft3d_info%zwork_len = max(n1,n2,n3) + 16*n3
#endif
#endif

    if (fft3d_info%iwork_len > 0) then
       allocate(fft3d_info%iwork(fft3d_info%iwork_len),stat=ierr)
       call utils_alloc_check('internal_serial_init (fourier_mod.F90)', &
            'fft3d_info%iwork',ierr)
    end if

    if (fft3d_info%dwork_len > 0) then
       allocate(fft3d_info%dwork(fft3d_info%dwork_len),stat=ierr)
       call utils_alloc_check('internal_serial_init (fourier_mod.F90)', &
            'fft3d_info%dwork',ierr)
    end if

    if (fft3d_info%zwork_len > 0) then
       allocate(fft3d_info%zwork(fft3d_info%zwork_len),stat=ierr)
       call utils_alloc_check('internal_serial_init (fourier_mod.F90)', &
            'fft3d_info%zwork',ierr)
    end if

    ! Platform-dependent initialisation

#ifdef FFTW
    call fftw3d_f77_create_plan(fft3d_info%forward_plan,n1,n2,n3, &
         FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
    call fftw3d_f77_create_plan(fft3d_info%backward_plan,n1,n2,n3, &
         FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
#endif

#ifdef FFTW3
#ifndef FFTW3_NO_OMP
    call dfftw_plan_with_nthreads(1)
#endif

    fft3d_info%iwork(1) = n1
    fft3d_info%iwork(2) = n2
    fft3d_info%iwork(3) = n3
    fft3d_info%iwork(4) = ld1
    fft3d_info%iwork(5) = ld2
    fft3d_info%iwork(6) = n3
    fft3d_info%iwork(7) = n2*n3
    fft3d_info%iwork(8) = n3*n1
    fft3d_info%iwork(9) = n1*n2
    fft3d_info%iwork(10) = n2*n3
    fft3d_info%iwork(11) = n3*n1
    fft3d_info%iwork(12) = n1*n2
    ! 3D FFTs
    call dfftw_plan_many_dft(fft3d_info%forward_plan,3, &
         fft3d_info%iwork(1:3),1, &
         fft3d_info%zwork,fft3d_info%iwork(4:6),1,0, &
         fft3d_info%zwork,fft3d_info%iwork(4:6),1,0, &
         FFTW_FORWARD,FFTW_MEASURE)
    call dfftw_plan_many_dft(fft3d_info%backward_plan,3, &
         fft3d_info%iwork(1:3),1, &
         fft3d_info%zwork,fft3d_info%iwork(4:6),1,0, &
         fft3d_info%zwork,fft3d_info%iwork(4:6),1,0, &
         FFTW_BACKWARD,FFTW_MEASURE)

    ! 1D FFTs

    ! n2*n3 transforms along '1'-dir, stride 1, dist n1
    call dfftw_plan_many_dft(fft3d_info%n1_f_plan, 1, & ! plan, rank
         fft3d_info%iwork(1), &   ! n
         fft3d_info%iwork(7), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(4), &   ! inembed
         1, &                     ! istride
         fft3d_info%iwork(4), &   ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(4), &   ! onembed
         1, &                     ! ostride
         fft3d_info%iwork(4), &   ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n1_b_plan, 1, & ! plan, rank
         fft3d_info%iwork(1), &   ! n
         fft3d_info%iwork(7), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(4), &   ! inembed
         1, &                     ! istride
         fft3d_info%iwork(4), &   ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(4), &   ! onembed
         1, &                     ! ostride
         fft3d_info%iwork(4), &   ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! n3 transforms along '2'-dir, stride n1, dist n1*n2
    call dfftw_plan_many_dft(fft3d_info%n2_f_plan, 1, & ! plan, rank
         fft3d_info%iwork(2), &   ! n
         fft3d_info%iwork(3), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(5), &   ! inembed
         fft3d_info%iwork(4), &   ! istride
         fft3d_info%iwork(12), &  ! dist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(5), &   ! onembed
         fft3d_info%iwork(4), &   ! ostride
         fft3d_info%iwork(12), &  ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n2_b_plan, 1, & ! plan, rank
         fft3d_info%iwork(2), &   ! n
         fft3d_info%iwork(3), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(5), &   ! inembed
         fft3d_info%iwork(4), &   ! istride
         fft3d_info%iwork(12), &  ! dist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(5), &   ! onembed
         fft3d_info%iwork(4), &   ! ostride
         fft3d_info%iwork(12), &  ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! n1*n2 transforms along '3'-dir, stride n1*n2, dist 1
    call dfftw_plan_many_dft(fft3d_info%n3_f_plan, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(9), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n3_b_plan, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(9), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! 1 transform along '3'-dir, stride n1*n2, dist 1
    call dfftw_plan_many_dft(fft3d_info%n3_f_plan_1, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         1, &                     ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n3_b_plan_1, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         1, &                     ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! n1/4+1 transforms along '3'-dir, stride n1*n2, dist 1
    call dfftw_plan_many_dft(fft3d_info%n3_f_plan_2, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(1)/4+1, & ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n3_b_plan_2, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(1)/4+1, & ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! n1/4 transforms along '3'-dir, stride n1*n2, dist 1
    call dfftw_plan_many_dft(fft3d_info%n3_f_plan_3, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(1)/4, & ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n3_b_plan_3, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(1)/4, & ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

#ifdef ParallelFFT

       if (pub_threads_fftbox) then
          ! kaw: Create parallel plans for data parallel OpenMP.

          ! kaw: Initialise threading for plans
#ifndef FFTW3_NO_OMP
          call dfftw_plan_with_nthreads(pub_threads_per_fftbox)
#endif

          ! 3D FFTs
          call dfftw_plan_many_dft(fft3d_info%forward_omp_plan,3, &
               fft3d_info%iwork(1:3),1, &
               fft3d_info%zwork,fft3d_info%iwork(4:6),1,0, &
               fft3d_info%zwork,fft3d_info%iwork(4:6),1,0, &
               FFTW_FORWARD,FFTW_MEASURE)
          call dfftw_plan_many_dft(fft3d_info%backward_omp_plan,3, &
               fft3d_info%iwork(1:3),1, &
               fft3d_info%zwork,fft3d_info%iwork(4:6),1,0, &
               fft3d_info%zwork,fft3d_info%iwork(4:6),1,0, &
               FFTW_BACKWARD,FFTW_MEASURE)



          ! 1D FFTs
          ! n2*n3 transforms along '1'-dir, stride 1, dist n1
          call dfftw_plan_many_dft(fft3d_info%n1_f_omp_plan, 1, & ! plan, rank
               fft3d_info%iwork(1), &   ! n
               fft3d_info%iwork(7), &   ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(4), &   ! inembed
               1, &                     ! istride
               fft3d_info%iwork(4), &   ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(4), &   ! onembed
               1, &                     ! ostride
               fft3d_info%iwork(4), &   ! odist
               FFTW_FORWARD, &          ! sign
               FFTW_MEASURE)            ! flags
          call dfftw_plan_many_dft(fft3d_info%n1_b_omp_plan, 1, & ! plan, rank
               fft3d_info%iwork(1), &   ! n
               fft3d_info%iwork(7), &   ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(4), &   ! inembed
               1, &                     ! istride
               fft3d_info%iwork(4), &   ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(4), &   ! onembed
               1, &                     ! ostride
               fft3d_info%iwork(4), &   ! odist
               FFTW_BACKWARD, &         ! sign
               FFTW_MEASURE)            ! flags

          ! n3 transforms along '2'-dir, stride n1, dist n1*n2
          call dfftw_plan_many_dft(fft3d_info%n2_f_omp_plan, 1, & ! plan, rank
               fft3d_info%iwork(2), &   ! n
               fft3d_info%iwork(3), &   ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(5), &   ! inembed
               fft3d_info%iwork(4), &   ! istride
               fft3d_info%iwork(12), &  ! dist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(5), &   ! onembed
               fft3d_info%iwork(4), &   ! ostride
               fft3d_info%iwork(12), &  ! odist
               FFTW_FORWARD, &          ! sign
               FFTW_MEASURE)            ! flags
          call dfftw_plan_many_dft(fft3d_info%n2_b_omp_plan, 1, & ! plan, rank
               fft3d_info%iwork(2), &   ! n
               fft3d_info%iwork(3), &   ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(5), &   ! inembed
               fft3d_info%iwork(4), &   ! istride
               fft3d_info%iwork(12), &  ! dist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(5), &   ! onembed
               fft3d_info%iwork(4), &   ! ostride
               fft3d_info%iwork(12), &  ! odist
               FFTW_BACKWARD, &         ! sign
               FFTW_MEASURE)            ! flags

          ! n1*n2 transforms along '3'-dir, stride n1*n2, dist 1
          call dfftw_plan_many_dft(fft3d_info%n3_f_omp_plan, 1, & ! plan, rank
               fft3d_info%iwork(3), &   ! n
               fft3d_info%iwork(9), &   ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(6), &   ! inembed
               fft3d_info%iwork(9), &   ! istride
               1, &                     ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(6), &   ! onembed
               fft3d_info%iwork(9), &   ! ostride
               1, &                     ! odist
               FFTW_FORWARD, &          ! sign
               FFTW_MEASURE)            ! flags
          call dfftw_plan_many_dft(fft3d_info%n3_b_omp_plan, 1, & ! plan, rank
               fft3d_info%iwork(3), &   ! n
               fft3d_info%iwork(9), &   ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(6), &   ! inembed
               fft3d_info%iwork(9), &   ! istride
               1, &                     ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(6), &   ! onembed
               fft3d_info%iwork(9), &   ! ostride
               1, &                     ! odist
               FFTW_BACKWARD, &         ! sign
               FFTW_MEASURE)            ! flags

          ! 1 transform along '3'-dir, stride n1*n2, dist 1
          call dfftw_plan_many_dft(fft3d_info%n3_f_omp_plan_1, 1, & ! plan, rank
               fft3d_info%iwork(3), &   ! n
               1, &                     ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(6), &   ! inembed
               fft3d_info%iwork(9), &   ! istride
               1, &                     ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(6), &   ! onembed
               fft3d_info%iwork(9), &   ! ostride
               1, &                     ! odist
               FFTW_FORWARD, &          ! sign
               FFTW_MEASURE)            ! flags
          call dfftw_plan_many_dft(fft3d_info%n3_b_omp_plan_1, 1, & ! plan, rank
               fft3d_info%iwork(3), &   ! n
               1, &                     ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(6), &   ! inembed
               fft3d_info%iwork(9), &   ! istride
               1, &                     ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(6), &   ! onembed
               fft3d_info%iwork(9), &   ! ostride
               1, &                     ! odist
               FFTW_BACKWARD, &         ! sign
               FFTW_MEASURE)            ! flags

          ! n1/4+1 transforms along '3'-dir, stride n1*n2, dist 1
          call dfftw_plan_many_dft(fft3d_info%n3_f_omp_plan_2, 1, & ! plan, rank
               fft3d_info%iwork(3), &   ! n
               fft3d_info%iwork(1)/4+1, & ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(6), &   ! inembed
               fft3d_info%iwork(9), &   ! istride
               1, &                     ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(6), &   ! onembed
               fft3d_info%iwork(9), &   ! ostride
               1, &                     ! odist
               FFTW_FORWARD, &          ! sign
               FFTW_MEASURE)            ! flags
          call dfftw_plan_many_dft(fft3d_info%n3_b_omp_plan_2, 1, & ! plan, rank
               fft3d_info%iwork(3), &   ! n
               fft3d_info%iwork(1)/4+1, & ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(6), &   ! inembed
               fft3d_info%iwork(9), &   ! istride
               1, &                     ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(6), &   ! onembed
               fft3d_info%iwork(9), &   ! ostride
               1, &                     ! odist
               FFTW_BACKWARD, &         ! sign
               FFTW_MEASURE)            ! flags

          ! n1/4 transforms along '3'-dir, stride n1*n2, dist 1
          call dfftw_plan_many_dft(fft3d_info%n3_f_omp_plan_3, 1, & ! plan, rank
               fft3d_info%iwork(3), &   ! n
               fft3d_info%iwork(1)/4, & ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(6), &   ! inembed
               fft3d_info%iwork(9), &   ! istride
               1, &                     ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(6), &   ! onembed
               fft3d_info%iwork(9), &   ! ostride
               1, &                     ! odist
               FFTW_FORWARD, &          ! sign
               FFTW_MEASURE)            ! flags
          call dfftw_plan_many_dft(fft3d_info%n3_b_omp_plan_3, 1, & ! plan, rank
               fft3d_info%iwork(3), &   ! n
               fft3d_info%iwork(1)/4, & ! howmany
               fft3d_info%zwork, &      ! in
               fft3d_info%iwork(6), &   ! inembed
               fft3d_info%iwork(9), &   ! istride
               1, &                     ! idist
               fft3d_info%zwork, &      ! out
               fft3d_info%iwork(6), &   ! onembed
               fft3d_info%iwork(9), &   ! ostride
               1, &                     ! odist
               FFTW_BACKWARD, &         ! sign
               FFTW_MEASURE)            ! flags

          ! kaw: Reset to create any subsequent plans with just a single thread.
#ifndef FFTW3_NO_OMP
          call dfftw_plan_with_nthreads(1)
#endif
       end if
#endif
#endif

#ifdef ACML
!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(NONE) &
!$OMP PRIVATE(tid,ierr) SHARED(fft3d_info,ld1,ld2,n1,n2,n3)
    tid = 0
!$  tid = omp_get_thread_num()
    call zfft3dy(ACML_MODE_PLAN,1.0_DP,.true.,n1,n2,n3, &
         fft3d_info%zwork(1),1,ld1,ld1*ld2,fft3d_info%zwork(1),1,ld1,ld1*ld2, &
         fft3d_info%tbl(1,tid),fft3d_info%tbl_size,ierr)
    call utils_assert(ierr == 0, 'Error in internal_serial_init &
            &(fourier_mod.F90): zfft3dy failed with code ',ierr)
!$OMP END PARALLEL
#endif

#ifdef VENDOR
#ifdef ALPHA
    ierr = zfft_init_3d(n1,n2,n3,fft3d_info%iwork(1),.true.)
    call utils_assert(ierr == 0, 'Error in internal_serial_init &
            &(fourier_mod.F90): zfft_init_3d failed with code ',ierr)
#endif
#ifdef SUN
    call zfftz3(0,n1,n2,n3,1.0_DP,zdum,ld1,ld2,zdum,ld1,ld2, &
         fft3d_info%dwork,fft3d_info%iwork,fft3d_info%zwork, &
         2*fft3d_info%zwork_len,ierr)
    call utils_assert(ierr == 0, 'Error in internal_serial_init &
            &(fourier_mod.F90): zfftz3 failed with code ',ierr)
#endif
#endif

  end subroutine internal_serial_init

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_parallel_init(grid,fft3d_info)

    use cell_grid, only: GRID_INFO
    use comms, only: pub_total_num_procs
#if defined (SUN) || defined (ACML)
    use constants, only: DP
#endif
#ifdef MKL_FFTW3
    use constants, only: I4B
#endif
    use rundat, only: pub_threads_cellfft, pub_threads_per_cellfft
    use utils, only : utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(parallel_fft3d_info), intent(out) :: fft3d_info

    ! External subroutines
#ifdef FFTW
    external :: rfftw_f77_create_plan, fftw2d_f77_create_plan
#endif
#ifdef FFTW3
    external :: dfftw_plan_many_dft, dfftw_plan_many_dft_r2c, &
         dfftw_plan_many_dft_c2r
#ifndef FFTW3_NO_OMP
    external :: dfftw_plan_with_nthreads
#endif
#endif
#ifdef ACML
    external :: dzfft, zfft2dx, zdfft
#endif
#ifdef VENDOR
#ifdef SUN
    external :: dfftzm, zfftz2, zfftdm
#endif
#endif

    ! Local variables
    integer :: n1,n2,n3    ! Dimensions of FFT grid
    integer :: ld1,ld2,ld3 ! Array dimensions containing grids
    integer :: ierr        ! Error flag
#ifdef SUN
    real(kind=DP) :: ddum
    complex(kind=DP) :: zdum
#endif

#ifdef MKL_FFTW3
    ! ndmh: code to ensure that when running MKL FFTW wrappers, the
    ! ndmh: plans are designed to be shared between multiple threads
#ifndef F2008
    !dec$ attributes align : 8 :: fftw3_mkl
    ! jd: ^ the above is not F2008-compliant.
    ! https://community.intel.com/t5/Intel-oneAPI-Math-Kernel-Library/
    ! /Calls-using-FFTW-wrapper-don-t-seem-to-be-thread-safe/td-p/945735
    ! claims that this is only needed to shut up the warning and that
    ! everything should work regardless.
#endif
    common/fftw3_mkl/ignore(4),number_of_user_threads,ignore2(7)
    integer(kind=I4B) :: ignore, number_of_user_threads, ignore2
    bind (c) :: /fftw3_mkl/

    number_of_user_threads = pub_threads_per_cellfft
#endif

    n1 = grid%n1
    n2 = grid%n2
    n3 = grid%n3
    ld1 = grid%ld1
    ld2 = grid%ld2
    ld3 = grid%ld3

    ! Simple checks on validity of input arguments
    call utils_assert(n1 >= 1, 'Error in internal_parallel_init &
            &(fourier_mod.F90): n1 < 1')
    call utils_assert(n2 >= 1, 'Error in internal_parallel_init &
            &(fourier_mod.F90): n2 < 1')
    call utils_assert(n2 >= 1, 'Error in internal_parallel_init &
            &(fourier_mod.F90): n3 < 1')

    call utils_assert(ld1 >= n1, 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld1 < n1')
    call utils_assert(ld2 >= n2, 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld2 < n2')
    call utils_assert(ld3 >= n3, 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld3 < n3')

    ! Check parallel strategy module has distributed fine cell FFTs

    call utils_assert(grid%distributed, 'Error in internal_parallel_init &
            &(fourier_mod.F90): no parallel strategy')

    ! Platform-dependent checking

#ifdef FFTW
    call utils_assert(ld1 == n1+2, 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld1 must equal n1+2 for FFTw')
    call utils_assert(ld2 == n2, 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld2 must equal n2 for FFTw')
    call utils_assert(ld3 == n3, 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld3 must equal n3 for FFTw')
#endif

    ! Simple initialisation

    ! Store basic information in fft3d_info structure

    fft3d_info%n1 = n1
    fft3d_info%n2 = n2
    fft3d_info%n3 = n3
    fft3d_info%ld1 = ld1
    fft3d_info%ld2 = ld2
    fft3d_info%ld3 = ld3
    fft3d_info%iwork_len = 0
    fft3d_info%dwork_len = 0
    fft3d_info%zwork_len = 0

    ! Store parallel strategy information in fft3d_info structure

    fft3d_info%num12slabs = grid%num_my_slabs12
    fft3d_info%max12slabs = grid%max_slabs12

    allocate(fft3d_info%idx12slab(0:pub_total_num_procs),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%idx12slab',ierr)

    fft3d_info%idx12slab(0:pub_total_num_procs-1) = &
         grid%first_slab12(0:pub_total_num_procs-1)
    fft3d_info%idx12slab(pub_total_num_procs) = n3 + 1

    fft3d_info%num23slabs = grid%num_slabs23
    fft3d_info%max23slabs = grid%max_slabs23

    allocate(fft3d_info%idx23slab(0:pub_total_num_procs),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%idx23slab',ierr)

    fft3d_info%idx23slab(0:pub_total_num_procs-1) = &
         grid%first_slab23(0:pub_total_num_procs-1)
    fft3d_info%idx23slab(pub_total_num_procs) = ld1/2 + 1

    ! General workspace allocation

    allocate(fft3d_info%buf12slab(fft3d_info%ld1/2,fft3d_info%n2),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%buf12slab',ierr)

    allocate(fft3d_info%buf23slab(fft3d_info%ld3,fft3d_info%n2),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%buf23slab',ierr)

    allocate(fft3d_info%fsendbuf(fft3d_info%max12slabs,fft3d_info%n2, &
         fft3d_info%ld1/2),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%fsendbuf',ierr)

    allocate(fft3d_info%bsendbuf(fft3d_info%max23slabs,fft3d_info%n2, &
         fft3d_info%n3),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%bsendbuf',ierr)

    allocate(fft3d_info%frecvbuf(fft3d_info%max12slabs,fft3d_info%n2, &
         fft3d_info%num23slabs),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%frecvbuf',ierr)

    allocate(fft3d_info%brecvbuf(fft3d_info%max23slabs,fft3d_info%n2, &
         fft3d_info%num12slabs),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%brecvbuf',ierr)

    ! Platform-dependent workspace allocation

#ifdef FFTW
    allocate(fft3d_info%rfftwbuf(fft3d_info%ld1,fft3d_info%n2),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%rfftwbuf',ierr)
#endif

#ifdef FFTW3
    fft3d_info%iwork_len = 7
    fft3d_info%dwork_len = ld1*n2
    fft3d_info%zwork_len = max(n2*ld1/2,ld2*ld3)
#endif

#ifdef ACML
    fft3d_info%zwork_len = max(n2*ld1/2,ld2*ld3)
    fft3d_info%dwork_len = max(n2*ld1/2,ld2*ld3)
    fft3d_info%tbl_size = ld1*ld2*n3 + 4*(ld1+ld2+n3) + 300
    allocate(fft3d_info%tbl(fft3d_info%tbl_size,3),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%tbl',ierr)
    allocate(fft3d_info%rbuf(ld1),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%rbuf',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
    fft3d_info%iwork_len = dxml_structure_size
#endif
#ifdef SUN
    fft3d_info%iwork_len = 2*128
    fft3d_info%dwork_len = 2*max(n1,n2+n3)
    fft3d_info%zwork_len = 2*max(n1,n2,n3)
#endif
#endif

    if (fft3d_info%iwork_len > 0) then
       allocate(fft3d_info%iwork(fft3d_info%iwork_len,3),stat=ierr)
       call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
            'fft3d_info%iwork',ierr)
    end if

    if (fft3d_info%dwork_len > 0) then
       allocate(fft3d_info%dwork(fft3d_info%dwork_len,3),stat=ierr)
       call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
            'fft3d_info%dwork',ierr)
    end if

    if (fft3d_info%zwork_len > 0) then
       allocate(fft3d_info%zwork(fft3d_info%zwork_len,3),stat=ierr)
       call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
            'fft3d_info%zwork',ierr)
    end if

    ! Platform-dependent initialisation

#ifdef FFTW
    call rfftw_f77_create_plan(fft3d_info%forward_plan(1),n1, &
         FFTW_REAL_TO_COMPLEX,FFTW_MEASURE+FFTW_USE_WISDOM)
    call fftw2d_f77_create_plan(fft3d_info%forward_plan(2),n3,n2, &
         FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
    call fftw2d_f77_create_plan(fft3d_info%backward_plan(2),n3,n2, &
         FFTW_BACKWARD,FFTW_MEASURE+FFTW_USE_WISDOM)
    call rfftw_f77_create_plan(fft3d_info%backward_plan(1),n1, &
         FFTW_COMPLEX_TO_REAL,FFTW_MEASURE+FFTW_USE_WISDOM)
#endif

#ifdef FFTW3
#ifndef FFTW3_NO_OMP
    call dfftw_plan_with_nthreads(1)
#endif
    fft3d_info%iwork(1,1) = n1
    fft3d_info%iwork(2,1) = n3
    fft3d_info%iwork(3,1) = n2
    fft3d_info%iwork(4,1) = ld1
    fft3d_info%iwork(5,1) = ld3
    fft3d_info%iwork(6,1) = ld2
    fft3d_info%iwork(7,1) = ld1/2
    call dfftw_plan_many_dft_r2c(fft3d_info%forward_plan(1),1, &
         fft3d_info%iwork(1,1),n2,fft3d_info%dwork,fft3d_info%iwork(4,1),1, &
         ld1,fft3d_info%zwork,fft3d_info%iwork(7,1),1,ld1/2,FFTW_MEASURE)
    call dfftw_plan_many_dft(fft3d_info%forward_plan(2),2, &
         fft3d_info%iwork(2:3,1),1, &
         fft3d_info%zwork,fft3d_info%iwork(5:6,1),1,0, &
         fft3d_info%zwork,fft3d_info%iwork(5:6,1),1,0, &
         FFTW_FORWARD,FFTW_MEASURE)
    call dfftw_plan_many_dft(fft3d_info%backward_plan(2),2, &
         fft3d_info%iwork(2:3,1),1, &
         fft3d_info%zwork(1,1),fft3d_info%iwork(5:6,1),1,0, &
         fft3d_info%zwork(1,2),fft3d_info%iwork(5:6,1),1,0, &
         FFTW_BACKWARD,FFTW_MEASURE)
    call dfftw_plan_many_dft_c2r(fft3d_info%backward_plan(1),1, &
         fft3d_info%iwork(1,1),n2,fft3d_info%dwork,fft3d_info%iwork(7,1),1, &
         ld1/2,fft3d_info%dwork,fft3d_info%iwork(4,1),1,ld1,FFTW_MEASURE)

#ifdef ParallelFFT
       if (pub_threads_cellfft) then
          ! kaw: Create parallel plans for data parallel OpenMP.

          ! kaw: Initialise threading for plans
#ifndef FFTW3_NO_OMP
          call dfftw_plan_with_nthreads(pub_threads_per_cellfft)
#endif
          call dfftw_plan_many_dft_r2c(fft3d_info%forward_plan_omp(1),1, &
               fft3d_info%iwork(1,1),n2,fft3d_info%dwork,fft3d_info%iwork(4,1),1, &
               ld1,fft3d_info%zwork,fft3d_info%iwork(7,1),1,ld1/2,FFTW_MEASURE)
          call dfftw_plan_many_dft(fft3d_info%forward_plan_omp(2),2, &
               fft3d_info%iwork(2:3,1),1, &
               fft3d_info%zwork,fft3d_info%iwork(5:6,1),1,0, &
               fft3d_info%zwork,fft3d_info%iwork(5:6,1),1,0, &
               FFTW_FORWARD,FFTW_MEASURE)
          call dfftw_plan_many_dft(fft3d_info%backward_plan_omp(2),2, &
               fft3d_info%iwork(2:3,1),1, &
               fft3d_info%zwork(1,1),fft3d_info%iwork(5:6,1),1,0, &
               fft3d_info%zwork(1,2),fft3d_info%iwork(5:6,1),1,0, &
               FFTW_BACKWARD,FFTW_MEASURE)
          call dfftw_plan_many_dft_c2r(fft3d_info%backward_plan_omp(1),1, &
               fft3d_info%iwork(1,1),n2,fft3d_info%dwork,fft3d_info%iwork(7,1),1, &
               ld1/2,fft3d_info%dwork,fft3d_info%iwork(4,1),1,ld1,FFTW_MEASURE)
      end if
#endif
#endif

#ifdef ACML
    call dzfft(ACML_MODE_PLAN,n1,fft3d_info%dwork(1,1),fft3d_info%tbl(1,1),ierr)
    call utils_assert(ierr == 0, 'Error in internal_parallel_init &
            &(fourier_mod.F90): dzfft failed with code ',ierr)
    call zfft2dx(ACML_MODE_PLAN,1.0_DP,.true.,.true.,n3,n2,fft3d_info%zwork(1,1), &
         1,ld3,fft3d_info%zwork(1,1),1,ld3,fft3d_info%tbl(1,2),ierr)
    call utils_assert(ierr == 0, 'Error in internal_parallel_init &
            &(fourier_mod.F90): dzfft failed with code ',ierr)
    call zdfft(ACML_MODE_PLAN,n1,fft3d_info%dwork(1,1),fft3d_info%tbl(1,3),ierr)
    call utils_assert(ierr == 0, 'Error in internal_parallel_init &
            &(fourier_mod.F90): zdfft failed with code ',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
    ierr = dfft_init(n1,fft3d_info%iwork(1,1),.true.)
    call utils_assert(ierr==0, 'Error in internal_parallel_init &
            &(fourier_mod.F90): dfft_init failed with code ',ierr)
    call utils_assert(ierr==0, 'Error in internal_parallel_init &
            &(fourier_mod.F90): zfft_init_2d failed with code ')
#endif
#ifdef SUN
    call dfftzm(0,n1,n2,1.0_DP,ddum,ld1,zdum,ld1/2,fft3d_info%dwork(:,1), &
         fft3d_info%iwork(:,1),fft3d_info%zwork(:,1),2*fft3d_info%zwork_len, &
         ierr)
    call utils_assert(ierr == 0,'Error in internal_parallel_init &
            &(fourier_mod.F90): dfftzm failed with code ',ierr)
    call zfftz2(0,n3,n2,1.0_DP,zdum,ld3,zdum,ld3,fft3d_info%dwork(:,2), &
         fft3d_info%iwork(:,2),fft3d_info%zwork(:,2),2*fft3d_info%zwork_len, &
         ierr)
    call utils_assert(ierr == 0, 'Error in internal_parallel_init &
            &(fourier_mod.F90): zfftz2 failed with code ',ierr)
    call zfftdm(0,n1,n2,1.0_DP,zdum,ld1/2,ddum,ld1,fft3d_info%dwork(:,3), &
         fft3d_info%iwork(:,3),fft3d_info%zwork(:,3),2*fft3d_info%zwork_len, &
         ierr)
    call utils_assert(ierr == 0, 'Error in internal_parallel_init &
            &(fourier_mod.F90): zfftdm failed with code ',ierr)
#endif
#endif

  end subroutine internal_parallel_init

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_fft3d_cell_forward(fft3d_info,rspc,gspc)

#ifdef VENDOR
    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_barrier
    use constants, only: DP
#else
    use comms, only:  pub_my_proc_id, pub_total_num_procs, comms_barrier, &
         comms_send, comms_recv
    use constants, only: DP
#endif
    ! jd: Some (older) MKL installation do not provide the .mod file, instead
    !     relying on a Fortran-style include. For these add
    !     '-DMKL_DONT_HAVE_SERVICE_MOD' to FFLAGS in the config file.
#if defined(MKL_FFTW3) && !defined(MKL_DONT_HAVE_SERVICE_MOD)
    ! jd: For environments which provide the .mod file.
    use mkl_service !! External dependency
#endif
    use rundat, only: pub_threads_cellfft, pub_threads_per_cellfft
    use utils, only: utils_assert, utils_use_var

    implicit none
#if defined(MKL_FFTW3) && defined(MKL_DONT_HAVE_SERVICE_MOD)
    ! jd: For environments which do not provide the .mod file.
    !     Ensure your -I path points to the directory with 'mkl_service.fi'
include "mkl_service.fi"
#endif

    type(parallel_fft3d_info), intent(inout) :: fft3d_info
    real(kind=DP), intent(in) :: &
         rspc(fft3d_info%ld1,fft3d_info%ld2,fft3d_info%max12slabs)
    complex(kind=DP), intent(out) :: &
         gspc(fft3d_info%ld3,fft3d_info%ld2,fft3d_info%max23slabs)

    ! External subroutines
#ifdef FFTW
    external :: rfftw_f77, fftwnd_f77_one
#endif
#ifdef FFTW3
    external :: dfftw_execute_dft_r2c, dfftw_execute_dft
#endif
#ifdef ACML
    external :: dzfft, zfft2dx
#endif
#ifdef VENDOR
#ifdef SUN
    external :: dfftzm, zfftz2
#endif
#endif

    ! Local variables
#ifdef ACML
    integer :: ierr                ! Error flag
    real(kind=DP) :: scale         ! Scale Factor
#endif
#ifdef VENDOR
    integer :: ierr                ! Error flag
#endif
#ifdef MKL_FFTW3
    integer :: ierr                ! Error flag
#endif
    integer :: i1,i2,i3            ! Grid loop counters
    integer :: islab12,islab23     ! Slab loop counters
    integer :: i1s,i1e             ! Start/end of 1-rods to send
    integer :: iproc               ! Proc loop counter
    integer :: send_id,recv_id     ! Send/receive proc

#ifdef MKL_FFTW3
    ierr = mkl_domain_set_num_threads(pub_threads_per_cellfft, MKL_DOMAIN_FFT)
#else
    call utils_use_var(pub_threads_per_cellfft)
#endif

    ! Zero output array
    gspc(:,:,:) = (0.0_DP,0.0_DP)

    !--------------------------------------------------------------------!
    !                                                                    !
    ! real   : rspc(ld1,ld2,num12slabs)            =  n(r1,r2,r3)        !
    !                       |                                            !
    !                       | 1D FFT along 1: r1 <-> g1                  !
    !                       V                                            !
    ! complex: buf12slab(ld1/2,n2)                 =  n(g1,r2)           !
    !                       |                         per 12-slab        !
    !                       | transpose                                  !
    !                       V                                            !
    ! complex: fsendbuf(max12slabs,n2,ld1/2)       =  n(r3,r2,g1)        !
    !                       |                                            !
    !                       | communicate                                !
    !                       V                                            !
    ! complex: frecvbuf(max12slabs,n2,num23slabs)  =  n(r3,r2,g1)        !
    !                       |                         per proc           !
    !                       | copy                                       !
    !                       V                                            !
    ! complex: gspc(ld3,ld2,num23slabs)            =  n(r3,r2,g1)        !
    !                       |                                            !
    !                       | 2D FFT in (2,3): r2,r3 <-> g2,g3           !
    !                       V                                            !
    ! complex: gspc(ld3,ld2,num23slabs)            =  n(g3,g2,g1)        !
    !                                                                    !
    !--------------------------------------------------------------------!

    ! Loop over 12-slabs on this proc

    do islab12=1,fft3d_info%num12slabs

       ! 1D FFT of 1-rods in this 12-slab on this proc

#ifdef FFTW
       call rfftw_f77(fft3d_info%forward_plan(1),fft3d_info%n2, &
            rspc(1,1,islab12),1,fft3d_info%ld1,fft3d_info%rfftwbuf, &
            1,fft3d_info%ld1)
       do i2=1,fft3d_info%n2
          fft3d_info%buf12slab(1,i2) = &
               cmplx(fft3d_info%rfftwbuf(1,i2),0.0_DP,kind=DP)
          do i1=2,fft3d_info%ld1/2-1
             fft3d_info%buf12slab(i1,i2) = &
                  cmplx(fft3d_info%rfftwbuf(i1,i2), &
                  fft3d_info%rfftwbuf(fft3d_info%ld1-i1,i2),kind=DP)
          end do
          fft3d_info%buf12slab(fft3d_info%ld1/2,i2) = &
               cmplx(fft3d_info%rfftwbuf(fft3d_info%ld1/2,i2),0.0_DP,kind=DP)
       end do
#endif

#ifdef FFTW3
#ifdef ParallelFFT
       if (pub_threads_cellfft) then
          call dfftw_execute_dft_r2c(fft3d_info%forward_plan_omp(1), &
               rspc(1,1,islab12),fft3d_info%buf12slab)
       else
#endif
          call dfftw_execute_dft_r2c(fft3d_info%forward_plan(1), &
               rspc(1,1,islab12),fft3d_info%buf12slab)
#ifdef ParallelFFT
       end if
#endif
#endif

#ifdef ACML
       do i2=1,fft3d_info%n2
          fft3d_info%rbuf(:) = rspc(:,i2,islab12)
          call dzfft(ACML_MODE_FORWARD,fft3d_info%n1,fft3d_info%rbuf(1), &
               fft3d_info%tbl(1,1),ierr)
          fft3d_info%buf12slab(1,i2) = &
               cmplx(fft3d_info%rbuf(1),0.0_DP,kind=DP)
          do i1=2,fft3d_info%ld1/2-1
             fft3d_info%buf12slab(i1,i2) = &
                  cmplx(fft3d_info%rbuf(i1), &
                  fft3d_info%rbuf(fft3d_info%ld1-i1),kind=DP)
          end do
          fft3d_info%buf12slab(fft3d_info%ld1/2,i2) = &
               cmplx(fft3d_info%rbuf(fft3d_info%ld1/2),0.0_DP,kind=DP)
       end do
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_forward &
            &(fourier_mod.F90): dzfft failed with code ',ierr)
       scale = sqrt(real(fft3d_info%n1,kind=DP))
       fft3d_info%buf12slab = fft3d_info%buf12slab * scale
#endif

#ifdef VENDOR
#ifdef ALPHA
       do i2=1,fft3d_info%n2
          ierr = dfft_apply('R','C','F',rspc(1,i2,islab12), &
               fft3d_info%buf12slab(1,i2),fft3d_info%iwork(1,1),1)
          call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_forward &
               &(fourier_mod.F90): dfft_apply failed with code ',ierr)
       end do
#endif
#ifdef SUN
       call dfftzm(-1,fft3d_info%n1,fft3d_info%n2,1.0_DP,rspc(:,:,islab12), &
            fft3d_info%ld1,fft3d_info%buf12slab,fft3d_info%ld1/2, &
            fft3d_info%dwork(:,1),fft3d_info%iwork(:,1), &
            fft3d_info%zwork(:,1),2*fft3d_info%zwork_len,ierr)
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_forward &
            &(fourier_mod.F90): dzfftm failed with code ',ierr)
#endif
#endif

       ! Transpose 1-rods in this 12-slab into 3-rods in 23-slab
       do i2=1,fft3d_info%n2
          do i1=1,fft3d_info%ld1/2
             fft3d_info%fsendbuf(islab12,i2,i1) = fft3d_info%buf12slab(i1,i2)
          end do
       end do

    end do

    ! Local communication phase
    i3 = fft3d_info%idx12slab(pub_my_proc_id)
    do islab12=1,fft3d_info%num12slabs
       do i2=1,fft3d_info%n2
          i1 = fft3d_info%idx23slab(pub_my_proc_id)
          do islab23=1,fft3d_info%num23slabs
             gspc(i3,i2,islab23) = fft3d_info%fsendbuf(islab12,i2,i1)
             i1 = i1 + 1
          end do
       end do
       i3 = i3 + 1
    end do

    ! Global communication phase
    call comms_barrier
    do iproc=1,pub_total_num_procs-1
       send_id = modulo(pub_my_proc_id + iproc,pub_total_num_procs)
       recv_id = modulo(pub_my_proc_id - iproc,pub_total_num_procs)

       ! Send packet to proc send_id
       i1s = fft3d_info%idx23slab(send_id)
       i1e = fft3d_info%idx23slab(send_id+1)-1
       call comms_send(send_id,fft3d_info%fsendbuf(:,:,i1s:i1e))

       ! Receive packet from proc recv_id
       call comms_recv(recv_id,fft3d_info%frecvbuf)

       ! Copy received data into gspc array
       do islab23=1,fft3d_info%num23slabs
          do i2=1,fft3d_info%n2
             i3 = fft3d_info%idx12slab(recv_id)
             do islab12=1,fft3d_info%idx12slab(recv_id+1)-i3
                gspc(i3,i2,islab23) = fft3d_info%frecvbuf(islab12,i2,islab23)
                i3 = i3 + 1
             end do
          end do
       end do

    end do

    ! Loop over 23-slabs on this proc

    do islab23=1,fft3d_info%num23slabs

       ! 2D FFT of this 23-slab on this proc

#ifdef FFTW
       call fftwnd_f77_one(fft3d_info%forward_plan(2),gspc(1,1,islab23), &
            gspc(1,1,islab23))
#endif

#ifdef FFTW3
#ifdef ParallelFFT
       if (pub_threads_cellfft) then
          call dfftw_execute_dft(fft3d_info%forward_plan_omp(2),gspc(1,1,islab23), &
               gspc(1,1,islab23))
       else
#endif
          call dfftw_execute_dft(fft3d_info%forward_plan(2),gspc(1,1,islab23), &
               gspc(1,1,islab23))
#ifdef ParallelFFT
       endif
#endif
#endif

#ifdef ACML
       call zfft2dx(ACML_MODE_BACKWARD,1.0_DP,.true.,.true.,fft3d_info%n3,fft3d_info%n2, &
            gspc(1,1,islab23),1,fft3d_info%ld3,gspc(1,1,islab23),1,fft3d_info%ld3, &
            fft3d_info%tbl(1,2),ierr)
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_forward &
            &(fourier_mod.F90): zfft2dx failed with code ',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
       ierr = zfft_apply_2d('C','C','F',gspc(1,1,islab23), &
            gspc(1,1,islab23),fft3d_info%ld3,fft3d_info%iwork(1,2),1,1)
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_forward &
            &(fourier_mod.F90): zfft_apply_2d failed with code ',ierr)
#endif
#ifdef SUN
       call zfftz2(-1,fft3d_info%n3,fft3d_info%n2,1.0_DP,gspc(:,:,islab23), &
            fft3d_info%ld3,gspc(:,:,islab23),fft3d_info%ld3, &
            fft3d_info%dwork(:,2),fft3d_info%iwork(:,2), &
            fft3d_info%zwork(:,2),2*fft3d_info%zwork_len,ierr)
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_forward &
            &(fourier_mod.F90): zfftz2 failed with code ',ierr)
#endif
#endif

    end do

#ifdef MKL_FFTW3
    ierr = mkl_domain_set_num_threads(1, MKL_DOMAIN_FFT)
#endif

    ! Synchronise all procs
    call comms_barrier

  end subroutine internal_fft3d_cell_forward

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_fft3d_cell_backward(fft3d_info,rspc,gspc)

#ifdef VENDOR
    use comms, only: pub_my_proc_id, pub_total_num_procs, comms_barrier
    use constants, only: DP
#else
    use comms, only:  pub_my_proc_id, pub_total_num_procs, comms_barrier, &
         comms_send, comms_recv
#endif
    ! jd: Some (older) MKL installation do not provide the .mod file, instead
    !     relying on a Fortran-style include. For these add
    !     '-DMKL_DONT_HAVE_SERVICE_MOD' to FFLAGS in the config file.
#if defined(MKL_FFTW3) && !defined(MKL_DONT_HAVE_SERVICE_MOD)
    ! jd: For environments which provide the .mod file.
    use mkl_service !! External dependency
#endif
    use rundat, only: pub_threads_cellfft, pub_threads_per_cellfft
    use utils, only: utils_assert, utils_use_var

    implicit none
#if defined(MKL_FFTW3) && defined(MKL_DONT_HAVE_SERVICE_MOD)
    ! jd: For environments which do not provide the .mod file.
    !     Ensure your -I path points to the directory with 'mkl_service.fi'
include "mkl_service.fi"
#endif

    type(parallel_fft3d_info), intent(inout) :: fft3d_info
    real(kind=DP), intent(out) :: &
         rspc(fft3d_info%ld1,fft3d_info%ld2,fft3d_info%max12slabs)
    complex(kind=DP), intent(in) :: &
         gspc(fft3d_info%ld3,fft3d_info%ld2,fft3d_info%max23slabs)

    ! External subroutines
#ifdef FFTW
    external :: fftwnd_f77_one, rfftw_f77
#endif
#ifdef FFTW3
     external :: dfftw_execute_dft, dfftw_execute_dft_c2r
#endif
#ifdef ACML
     external :: zfft2dx, zdfft
#endif
#ifdef VENDOR
#ifdef SUN
     external :: zfftz2, zfftdm
#endif
#endif

    ! Local variables
    integer :: i1,i2,i3            ! Grid loop counters
    integer :: islab12,islab23     ! Slab loop counters
    integer :: i3s,i3e             ! Start/end of 3-rods to send
    integer :: iproc               ! Proc loop counter
    integer :: send_id,recv_id     ! Send/receive proc
    integer :: ierr
    real(kind=DP) :: scale

#ifdef MKL_FFTW3
    ierr = mkl_domain_set_num_threads(pub_threads_per_cellfft, MKL_DOMAIN_FFT)
#else
    call utils_use_var(pub_threads_per_cellfft)
    ierr = 0  ! jme: this prevents 'unused variable' compilation warnings
#endif

    ! Zero output array
    rspc(:,:,:) = 0.0_DP

    !--------------------------------------------------------------------!
    !                                                                    !
    ! complex: gspc(ld3,ld2,num23slabs)            =  n(g3,g2,g1)        !
    !                       |                                            !
    !                       | 2D FFT in (2,3): g2,g3 <-> r2,r3           !
    !                       V                                            !
    ! complex: buf23slab(ld3,n2)                   =  n(r3,r2)           !
    !                       |                         per 23-slab        !
    !                       | transpose                                  !
    !                       V                                            !
    ! complex: bsendbuf(max23slabs,n2,n3)          =  n(g1,r2,r3)        !
    !                       |                                            !
    !                       | communicate                                !
    !                       V                                            !
    ! complex: brecvbuf(max23slabs,n2,num12slabs)  =  n(g1,r2,r3)        !
    !                       |                         per proc           !
    !                       | copy                                       !
    !                       V                                            !
    ! real   : rspc(ld1,ld2,num12slabs)            =  n(g1,r2,r3)        !
    !                       |                                            !
    !                       | 1D FFT along 1: g1 <-> r1                  !
    !                       V                                            !
    ! real   : rspc(ld1,ld2,num12slabs)            =  n(r1,r2,r3)        !
    !                                                                    !
    !--------------------------------------------------------------------!

    ! Loop over 23-slabs on this proc

    do islab23=1,fft3d_info%num23slabs

       ! 2D FFT of this 23-slab on this proc

#ifdef FFTW
       call fftwnd_f77_one(fft3d_info%backward_plan(2),gspc(1,1,islab23), &
            fft3d_info%buf23slab)
       scale = 1.0_DP / (fft3d_info%n2*fft3d_info%n3)
       fft3d_info%buf23slab = fft3d_info%buf23slab * scale
#endif

#ifdef FFTW3
#ifdef ParallelFFT
       if (pub_threads_cellfft) then
          call dfftw_execute_dft(fft3d_info%backward_plan_omp(2), &
               gspc(1,1,islab23),fft3d_info%buf23slab)
       else
#endif
          call dfftw_execute_dft(fft3d_info%backward_plan(2), &
               gspc(1,1,islab23),fft3d_info%buf23slab)
#ifdef ParallelFFT
       endif
#endif
       scale = 1.0_DP / (fft3d_info%n2*fft3d_info%n3)
       fft3d_info%buf23slab = fft3d_info%buf23slab * scale
#endif

#ifdef ACML
       scale = 1.0_DP / (fft3d_info%n2*fft3d_info%n3)
       call zfft2dx(ACML_MODE_FORWARD,scale,.true.,.false.,fft3d_info%n3,fft3d_info%n2, &
            gspc(1,1,islab23),1,fft3d_info%ld3,fft3d_info%buf23slab(1,1),1,fft3d_info%ld3, &
            fft3d_info%tbl(1,2),ierr)
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_backward &
            &(fourier_mod.F90): dzfft2dx failed with code ',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
       ierr = zfft_apply_2d('C','C','B',gspc(:,:,islab23), &
            fft3d_info%buf23slab,fft3d_info%ld3,fft3d_info%iwork(:,2),1,1)
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_backward &
            &(fourier_mod.F90): zfft_apply_2d failed with code ',ierr)
#endif
#ifdef SUN
       scale = 1.0_DP / (fft3d_info%n2*fft3d_info%n3)
       call zfftz2(1,fft3d_info%n3,fft3d_info%n2,scale,gspc(:,:,islab23), &
            fft3d_info%ld3,fft3d_info%buf23slab,fft3d_info%ld3, &
            fft3d_info%dwork(:,2),fft3d_info%iwork(:,2), &
            fft3d_info%zwork(:,2),2*fft3d_info%zwork_len,ierr)
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_backward &
            &(fourier_mod.F90): zfftz2 failed with code ',ierr)
#endif
#endif

       ! Transpose 3-rods in this 23-slab into 1-rods in 12-slab
       do i2=1,fft3d_info%n2
          do i3=1,fft3d_info%n3
             fft3d_info%bsendbuf(islab23,i2,i3) = fft3d_info%buf23slab(i3,i2)
          end do
       end do

    end do

    ! Local communication phase
    i3 = fft3d_info%idx12slab(pub_my_proc_id)
    do islab12=1,fft3d_info%num12slabs
       do i2=1,fft3d_info%n2
          i1 = fft3d_info%idx23slab(pub_my_proc_id)
          do islab23=1,fft3d_info%num23slabs
             rspc(2*i1-1,i2,islab12) = &
                  real(fft3d_info%bsendbuf(islab23,i2,i3),kind=DP)
             rspc(2*i1,i2,islab12) = &
                  aimag(fft3d_info%bsendbuf(islab23,i2,i3))
             i1 = i1 + 1
          end do
       end do
       i3 = i3 + 1
    end do

    ! Global communication phase
    call comms_barrier
    do iproc=1,pub_total_num_procs-1
       send_id = modulo(pub_my_proc_id + iproc,pub_total_num_procs)
       recv_id = modulo(pub_my_proc_id - iproc,pub_total_num_procs)

       ! Send packet to proc send_id
       i3s = fft3d_info%idx12slab(send_id)
       i3e = fft3d_info%idx12slab(send_id+1)-1
       call comms_send(send_id,fft3d_info%bsendbuf(:,:,i3s:i3e))

       ! Receive packet from proc recv_id
       call comms_recv(recv_id,fft3d_info%brecvbuf)

       ! Copy received data into rspc array
       do islab12=1,fft3d_info%num12slabs
          do i2=1,fft3d_info%n2
             i1 = fft3d_info%idx23slab(recv_id)
             do islab23=1,fft3d_info%idx23slab(recv_id+1)-i1
                rspc(2*i1-1,i2,islab12) = &
                     real(fft3d_info%brecvbuf(islab23,i2,islab12),kind=DP)
                rspc(2*i1,i2,islab12) = &
                     aimag(fft3d_info%brecvbuf(islab23,i2,islab12))
                i1 = i1 + 1
             end do
          end do
       end do

    end do

    ! Loop over 12-slabs on this proc

    do islab12=1,fft3d_info%num12slabs

       ! 1D FFT of 1-rods in this 12-slab on this proc

#ifdef FFTW
       do i2=1,fft3d_info%n2
          fft3d_info%rfftwbuf(1,i2) = rspc(1,i2,islab12)
          do i1=2,fft3d_info%ld1/2-1
             fft3d_info%rfftwbuf(i1,i2) = rspc(2*i1-1,i2,islab12)
             fft3d_info%rfftwbuf(fft3d_info%ld1-i1,i2) = &
                  rspc(2*i1,i2,islab12)
          end do
          fft3d_info%rfftwbuf(fft3d_info%ld1/2,i2) = &
               rspc(fft3d_info%ld1-1,i2,islab12)
       end do
       call rfftw_f77(fft3d_info%backward_plan(1),fft3d_info%n2, &
            fft3d_info%rfftwbuf,1,fft3d_info%ld1,rspc(1,1,islab12), &
            1,fft3d_info%ld1)
       scale = 1.0_DP / fft3d_info%n1
       rspc(:,:,islab12) = rspc(:,:,islab12) * scale
#endif

#ifdef FFTW3
#ifdef ParallelFFT
       if (pub_threads_cellfft) then
          call dfftw_execute_dft_c2r(fft3d_info%backward_plan_omp(1), &
               rspc(1,1,islab12),rspc(1,1,islab12))
       else
#endif
          call dfftw_execute_dft_c2r(fft3d_info%backward_plan(1), &
               rspc(1,1,islab12),rspc(1,1,islab12))
#ifdef ParallelFFT
       endif
#endif
       scale = 1.0_DP / fft3d_info%n1
       rspc(:,:,islab12) = rspc(:,:,islab12) * scale
#endif

#ifdef ACML
       do i2=1,fft3d_info%n2
          fft3d_info%rbuf(1) = rspc(1,i2,islab12)
          do i1=2,fft3d_info%ld1/2-1
             fft3d_info%rbuf(i1) = rspc(2*i1-1,i2,islab12)
             fft3d_info%rbuf(fft3d_info%ld1-i1) = &
                  -rspc(2*i1,i2,islab12)
          end do
          fft3d_info%rbuf(fft3d_info%ld1/2) = &
               rspc(fft3d_info%ld1-1,i2,islab12)
          call zdfft(ACML_MODE_FORWARD,fft3d_info%n1,fft3d_info%rbuf(1), &
               fft3d_info%tbl(1,3),ierr)
          scale = sqrt(1.0_DP / fft3d_info%n1)
          rspc(:,i2,islab12) = fft3d_info%rbuf(:) * scale
       end do
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_backward &
            &(fourier_mod.F90): zdfft failed with code ',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
       do i2=1,fft3d_info%n2
          ierr = dfft_apply('C','R','B',rspc(1,i2,islab12), &
               rspc(1,i2,islab12),fft3d_info%iwork(1,1),1)
          call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_backward &
               &(fourier_mod.F90): dfft_apply failed with code ',ierr)
       end do
#endif
#ifdef SUN
       scale = 1.0_DP / fft3d_info%n1
       call zfftdm(1,fft3d_info%n1,fft3d_info%n2,scale,rspc(:,:,islab12), &
            fft3d_info%ld1/2,rspc(:,:,islab12),fft3d_info%ld1, &
            fft3d_info%dwork(:,3),fft3d_info%iwork(:,3), &
            fft3d_info%zwork(:,3),2*fft3d_info%zwork_len,ierr)
       call utils_assert(ierr == 0, 'Error in internal_fft3d_cell_backward &
            &(fourier_mod.F90): zfftdm failed with code ',ierr)
#endif
#endif

    end do

#ifdef MKL_FFTW3
    ierr = mkl_domain_set_num_threads(1, MKL_DOMAIN_FFT)
#endif

    ! Synchronise all procs
    call comms_barrier

  end subroutine internal_fft3d_cell_backward

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_fft3d_box(fft3d_info,dir,data,apply_scale,padbox,omp)

    ! jd: Some (older) MKL installation do not provide the .mod file, instead
    !     relying on a Fortran-style include. For these add
    !     '-DMKL_DONT_HAVE_SERVICE_MOD' to FFLAGS in the config file.
#if defined(MKL_FFTW3) && !defined(MKL_DONT_HAVE_SERVICE_MOD)
    ! jd: For environments which provide the .mod file.
    use mkl_service !! External dependency
#endif
#ifdef MKL_FFTW3
    use rundat, only: pub_threads_per_fftbox
#endif
    use utils, only: utils_assert

    implicit none
#if defined(MKL_FFTW3) && defined(MKL_DONT_HAVE_SERVICE_MOD)
    ! jd: For environments which do not provide the .mod file.
    !     Ensure your -I path points to the directory with 'mkl_service.fi'
include "mkl_service.fi"
#endif

    type(serial_fft3d_info), intent(inout) :: fft3d_info
    character, intent(in) :: dir
    complex(kind=DP), intent(inout) :: &
         data(fft3d_info%ld1,fft3d_info%ld2,fft3d_info%n3)
    logical, intent(in) :: apply_scale
    logical, intent(in), optional :: omp
    logical, intent(in), optional :: padbox

    ! External subroutines and functions
#ifdef FFTW
    external :: fftwnd_f77_one
#endif
#ifdef FFTW3
    external :: dfftw_execute_dft
#endif
#ifdef ACML
    external :: zfft3dy
#endif
#ifdef VENDOR
#ifdef SUN
    external :: zfftz3
#endif
#endif
!$  integer, external :: omp_get_thread_num

    ! Local variables
#ifdef ACML
    integer :: ierr
#endif
#ifdef VENDOR
    integer :: ierr             ! Error flag
#endif
#ifdef MKL_FFTW3
    integer :: ierr             ! Error flag
#endif
    real(kind=DP) :: scale
#ifdef FFTW
    !   complex(kind=DP) :: zdum
    complex(kind=DP) :: zdum(1) !cks, 13/10/2004: Hack to compile with NAG f95
#endif
    integer :: tid
    integer :: i1,i2,j1
    logical :: do_padbox, open_mp

#ifdef MKL_FFTW3
    ierr = mkl_domain_set_num_threads(pub_threads_per_fftbox, MKL_DOMAIN_FFT)
#endif
    tid = 0
!$  tid = omp_get_thread_num()

    ! Optional execution of data parallel transforms
    if (present(omp)) then
       open_mp = omp
    else
       open_mp = .false.
    end if

    do_padbox = .false.
    if (present(padbox)) do_padbox = padbox

    ! Platform-dependent call

#ifdef FFTW
    if (dir == 'F' .or. dir == 'f') then
       call fftwnd_f77_one(fft3d_info%forward_plan,data,zdum)
    else
       call fftwnd_f77_one(fft3d_info%backward_plan,data,zdum)
       scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
       if (apply_scale) then
          data = data * scale
       end if
    end if
#endif

#ifdef FFTW3
    if (open_mp) then
       if (dir == 'F' .or. dir == 'f') then
          if (.not.do_padbox) then
             call dfftw_execute_dft(fft3d_info%forward_omp_plan,data,data)
          else

             ! ndmh: padded version: assumes G-vectors where:
             ! ndmh:      n1/4+2 <= i1 <= 3*n1/4+1
             ! ndmh:  or  n2/4+2 <= i2 <= 3*n2/4+1
             ! ndmh:  or  n3/4+2 <= i3 <= 3*n3/4+1
             ! ndmh: are not required, so avoids as many of the 1D transforms
             ! ndmh: as possible

             ! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
             call dfftw_execute_dft(fft3d_info%n1_f_omp_plan,data(1,1,1),data(1,1,1))
             ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
             do i1=1,fft3d_info%n1/4+1
                call dfftw_execute_dft(fft3d_info%n2_f_omp_plan,data(i1,1,1),data(i1,1,1))
             end do
             ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
             do i1=3*fft3d_info%n1/4+2,fft3d_info%n1
                call dfftw_execute_dft(fft3d_info%n2_f_omp_plan,data(i1,1,1),data(i1,1,1))
             end do
             ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
             ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
             j1 = 3*fft3d_info%n1/4+2
             do i2=1,fft3d_info%n2/4+1
                 call dfftw_execute_dft(fft3d_info%n3_f_omp_plan_2,data(1,i2,1),data(1,i2,1))
                 call dfftw_execute_dft(fft3d_info%n3_f_omp_plan_3,data(j1,i2,1),data(j1,i2,1))
             end do
             ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
             ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
             do i2=3*fft3d_info%n2/4+2,fft3d_info%n2
                 call dfftw_execute_dft(fft3d_info%n3_f_omp_plan_2,data(1,i2,1),data(1,i2,1))
                 call dfftw_execute_dft(fft3d_info%n3_f_omp_plan_3,data(j1,i2,1),data(j1,i2,1))
             end do

#if 0
             ! ndmh: full multi-1D-FFT version, equivalent to full 3D transform
             call dfftw_execute_dft(fft3d_info%n3_f_omp_plan,data(1,1,1),data(1,1,1))
             do i1=1,fft3d_info%n1
                call dfftw_execute_dft(fft3d_info%n2_f_omp_plan,data(i1,1,1),data(i1,1,1))
             end do
             call dfftw_execute_dft(fft3d_info%n1_f_omp_plan,data(1,1,1),data(1,1,1))
#endif

          end if
       else
          if (.not.do_padbox) then
             ! ndmh: normal full 3D FFT version
             call dfftw_execute_dft(fft3d_info%backward_omp_plan,data,data)
          else

             ! ndmh: padded version: assumes G-vectors where:
             ! ndmh:      n1/4+2 <= i1 <= 3*n1/4+1
             ! ndmh:  or  n2/4+2 <= i2 <= 3*n2/4+1
             ! ndmh:  or  n3/4+2 <= i3 <= 3*n3/4+1
             ! ndmh: are all zero, so avoids as many of the 1D transforms
             ! ndmh: as possible

             ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
             ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
             j1 = 3*fft3d_info%n1/4+2
             do i2=1,fft3d_info%n2/4+1
                 call dfftw_execute_dft(fft3d_info%n3_b_omp_plan_2,data(1,i2,1),data(1,i2,1))
                 call dfftw_execute_dft(fft3d_info%n3_b_omp_plan_3,data(j1,i2,1),data(j1,i2,1))
             end do
             ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
             ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
             do i2=3*fft3d_info%n2/4+2,fft3d_info%n2
                 call dfftw_execute_dft(fft3d_info%n3_b_omp_plan_2,data(1,i2,1),data(1,i2,1))
                 call dfftw_execute_dft(fft3d_info%n3_b_omp_plan_3,data(j1,i2,1),data(j1,i2,1))
             end do
             ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
             do i1=1,fft3d_info%n1/4+1
                call dfftw_execute_dft(fft3d_info%n2_b_omp_plan,data(i1,1,1),data(i1,1,1))
             end do
             ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
             do i1=3*fft3d_info%n1/4+2,fft3d_info%n1
                call dfftw_execute_dft(fft3d_info%n2_b_omp_plan,data(i1,1,1),data(i1,1,1))
             end do
             ! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
             call dfftw_execute_dft(fft3d_info%n1_b_omp_plan,data(1,1,1),data(1,1,1))

#if 0
             ! ndmh: full multi-1D-FFT version, equivalent to full 3D transform
             call dfftw_execute_dft(fft3d_info%n3_b_omp_plan,data(1,1,1),data(1,1,1))
             do i1=1,fft3d_info%n1
                call dfftw_execute_dft(fft3d_info%n2_b_omp_plan,data(i1,1,1),data(i1,1,1))
             end do
             call dfftw_execute_dft(fft3d_info%n1_b_omp_plan,data(1,1,1),data(1,1,1))
#endif
          end if
          scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
          if (apply_scale) then
             data = data * scale
          end if
       end if

    else

       if (dir == 'F' .or. dir == 'f') then
          if (.not.do_padbox) then
             call dfftw_execute_dft(fft3d_info%forward_plan,data,data)
          else

             ! ndmh: padded version: assumes G-vectors where:
             ! ndmh:      n1/4+2 <= i1 <= 3*n1/4+1
             ! ndmh:  or  n2/4+2 <= i2 <= 3*n2/4+1
             ! ndmh:  or  n3/4+2 <= i3 <= 3*n3/4+1
             ! ndmh: are not required, so avoids as many of the 1D transforms
             ! ndmh: as possible

             ! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
             call dfftw_execute_dft(fft3d_info%n1_f_plan,data(1,1,1),data(1,1,1))
             ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
             do i1=1,fft3d_info%n1/4+1
                call dfftw_execute_dft(fft3d_info%n2_f_plan,data(i1,1,1),data(i1,1,1))
             end do
             ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
             do i1=3*fft3d_info%n1/4+2,fft3d_info%n1
                call dfftw_execute_dft(fft3d_info%n2_f_plan,data(i1,1,1),data(i1,1,1))
             end do
             ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
             ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
             j1 = 3*fft3d_info%n1/4+2
             do i2=1,fft3d_info%n2/4+1
                 call dfftw_execute_dft(fft3d_info%n3_f_plan_2,data(1,i2,1),data(1,i2,1))
                 call dfftw_execute_dft(fft3d_info%n3_f_plan_3,data(j1,i2,1),data(j1,i2,1))
             end do
             ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
             ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
             do i2=3*fft3d_info%n2/4+2,fft3d_info%n2
                 call dfftw_execute_dft(fft3d_info%n3_f_plan_2,data(1,i2,1),data(1,i2,1))
                 call dfftw_execute_dft(fft3d_info%n3_f_plan_3,data(j1,i2,1),data(j1,i2,1))
             end do

#if 0
             ! ndmh: full multi-1D-FFT version, equivalent to full 3D transform
             call dfftw_execute_dft(fft3d_info%n3_f_plan,data(1,1,1),data(1,1,1))
             do i1=1,fft3d_info%n1
                call dfftw_execute_dft(fft3d_info%n2_f_plan,data(i1,1,1),data(i1,1,1))
             end do
             call dfftw_execute_dft(fft3d_info%n1_f_plan,data(1,1,1),data(1,1,1))
#endif

          end if
       else
          if (.not.do_padbox) then
             ! ndmh: normal full 3D FFT version
             call dfftw_execute_dft(fft3d_info%backward_plan,data,data)
          else

             ! ndmh: padded version: assumes G-vectors where:
             ! ndmh:      n1/4+2 <= i1 <= 3*n1/4+1
             ! ndmh:  or  n2/4+2 <= i2 <= 3*n2/4+1
             ! ndmh:  or  n3/4+2 <= i3 <= 3*n3/4+1
             ! ndmh: are all zero, so avoids as many of the 1D transforms
             ! ndmh: as possible

             ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
             ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
             j1 = 3*fft3d_info%n1/4+2
             do i2=1,fft3d_info%n2/4+1
                 call dfftw_execute_dft(fft3d_info%n3_b_plan_2,data(1,i2,1),data(1,i2,1))
                 call dfftw_execute_dft(fft3d_info%n3_b_plan_3,data(j1,i2,1),data(j1,i2,1))
             end do
             ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
             ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
             do i2=3*fft3d_info%n2/4+2,fft3d_info%n2
                 call dfftw_execute_dft(fft3d_info%n3_b_plan_2,data(1,i2,1),data(1,i2,1))
                 call dfftw_execute_dft(fft3d_info%n3_b_plan_3,data(j1,i2,1),data(j1,i2,1))
             end do
             ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
             do i1=1,fft3d_info%n1/4+1
                call dfftw_execute_dft(fft3d_info%n2_b_plan,data(i1,1,1),data(i1,1,1))
             end do
             ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
             do i1=3*fft3d_info%n1/4+2,fft3d_info%n1
                call dfftw_execute_dft(fft3d_info%n2_b_plan,data(i1,1,1),data(i1,1,1))
             end do
             ! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
             call dfftw_execute_dft(fft3d_info%n1_b_plan,data(1,1,1),data(1,1,1))

#if 0
             ! ndmh: full multi-1D-FFT version, equivalent to full 3D transform
             call dfftw_execute_dft(fft3d_info%n3_b_plan,data(1,1,1),data(1,1,1))
             do i1=1,fft3d_info%n1
                call dfftw_execute_dft(fft3d_info%n2_b_plan,data(i1,1,1),data(i1,1,1))
             end do
             call dfftw_execute_dft(fft3d_info%n1_b_plan,data(1,1,1),data(1,1,1))
#endif
          end if
          scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
          if (apply_scale) then
             data = data * scale
          end if
       end if
    end if
#endif

#ifdef ACML
    ! ndmh: Note ACML's definition of forward and backward seem to be reversed
    ! ndmh: with respect to everyone elses' definitions, hence the flags are
    ! ndmh: passed in seemingly the wrong way round.
    if (dir == 'F' .or. dir == 'f') then
       call zfft3dy(ACML_MODE_BACKWARD,1.0_DP,.true.,fft3d_info%n1, &
            fft3d_info%n2,fft3d_info%n3,data(1,1,1),1,fft3d_info%ld1, &
            fft3d_info%ld1*fft3d_info%ld2,data(1,1,1),1,fft3d_info%ld1, &
            fft3d_info%ld1*fft3d_info%ld2,fft3d_info%tbl(1,tid), &
            fft3d_info%tbl_size,ierr)
    else
       if (apply_scale) then
          scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
       else
          scale = 1.0_DP
       end if
       call zfft3dy(ACML_MODE_FORWARD,scale,.true.,fft3d_info%n1, &
            fft3d_info%n2,fft3d_info%n3,data(1,1,1),1,fft3d_info%ld1, &
            fft3d_info%ld1*fft3d_info%ld2,data(1,1,1),1,fft3d_info%ld1, &
            fft3d_info%ld1*fft3d_info%ld2,fft3d_info%tbl(1,tid), &
            fft3d_info%tbl_size,ierr)
    end if

    call utils_assert(ierr == 0, 'Error in internal_fft3d_box &
         &(fourier_mod.F90): zfft3dy failed with code ',ierr)

#endif

#ifdef VENDOR
#ifdef ALPHA
    if (dir == 'F' .or. dir == 'f') then
       ierr = zfft_apply_3d('C','C','F',data,data,fft3d_info%ld1, &
            fft3d_info%ld2,fft3d_info%iwork,1,1,1)
    else
       ierr = zfft_apply_3d('C','C','B',data,data,fft3d_info%ld1, &
            fft3d_info%ld2,fft3d_info%iwork,1,1,1)
    end if

    call utils_assert(ierr == 0, 'Error in internal_fft3d_box &
         &(fourier_mod.F90): zfft_apply_3d failed with code ',ierr)

#endif
#ifdef SUN
    if (dir == 'F' .or. dir == 'f') then
       call zfftz3(-1,fft3d_info%n1,fft3d_info%n2,fft3d_info%n3,1.0_DP, &
            data,fft3d_info%ld1,fft3d_info%ld2,data,fft3d_info%ld1, &
            fft3d_info%ld2,fft3d_info%dwork,fft3d_info%iwork, &
            fft3d_info%zwork,2*fft3d_info%zwork_len,ierr)
    else
       if (apply_scale) then
          scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
       else
          scale = 1.0_DP
       end if
       call zfftz3(1,fft3d_info%n1,fft3d_info%n2,fft3d_info%n3,scale, &
            data,fft3d_info%ld1,fft3d_info%ld2,data,fft3d_info%ld1, &
            fft3d_info%ld2,fft3d_info%dwork,fft3d_info%iwork, &
            fft3d_info%zwork,2*fft3d_info%zwork_len,ierr)
    end if
    call utils_assert(ierr == 0, 'Error in internal_fft3d_box &
         &(fourier_mod.F90): zfftz3 failed with code ',ierr)
#endif
#endif

#ifdef MKL_FFTW3
    ierr = mkl_domain_set_num_threads(1, MKL_DOMAIN_FFT)
#endif

  end subroutine internal_fft3d_box

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_serial_exit(fft3d_info)

    use rundat, only: pub_threads_fftbox
    use utils, only : utils_dealloc_check, utils_assert

    implicit none

    type(serial_fft3d_info), intent(inout) :: fft3d_info

    ! External subroutines
#ifdef FFTW
    external :: fftwnd_f77_destroy_plan
#endif
#ifdef FFTW3
    external :: dfftw_destroy_plan
#endif

    ! Local variables

    integer :: ierr             ! Error flag

    ! Platform-dependent finalisation

#ifdef FFTW
    call fftwnd_f77_destroy_plan(fft3d_info%forward_plan)
    call fftwnd_f77_destroy_plan(fft3d_info%backward_plan)
#endif

#ifdef FFTW3

    call dfftw_destroy_plan(fft3d_info%n1_f_plan)
    call dfftw_destroy_plan(fft3d_info%n1_b_plan)
    call dfftw_destroy_plan(fft3d_info%n2_f_plan)
    call dfftw_destroy_plan(fft3d_info%n2_b_plan)
    call dfftw_destroy_plan(fft3d_info%n3_f_plan)
    call dfftw_destroy_plan(fft3d_info%n3_b_plan)
    call dfftw_destroy_plan(fft3d_info%n3_f_plan_1)
    call dfftw_destroy_plan(fft3d_info%n3_b_plan_1)
    call dfftw_destroy_plan(fft3d_info%n3_f_plan_2)
    call dfftw_destroy_plan(fft3d_info%n3_b_plan_2)
    call dfftw_destroy_plan(fft3d_info%n3_f_plan_3)
    call dfftw_destroy_plan(fft3d_info%n3_b_plan_3)

    call dfftw_destroy_plan(fft3d_info%forward_plan)
    call dfftw_destroy_plan(fft3d_info%backward_plan)

#ifdef ParallelFFT
    if (pub_threads_fftbox) then

       call dfftw_destroy_plan(fft3d_info%n1_f_omp_plan)
       call dfftw_destroy_plan(fft3d_info%n1_b_omp_plan)
       call dfftw_destroy_plan(fft3d_info%n2_f_omp_plan)
       call dfftw_destroy_plan(fft3d_info%n2_b_omp_plan)
       call dfftw_destroy_plan(fft3d_info%n3_f_omp_plan)
       call dfftw_destroy_plan(fft3d_info%n3_b_omp_plan)
       call dfftw_destroy_plan(fft3d_info%n3_f_omp_plan_1)
       call dfftw_destroy_plan(fft3d_info%n3_b_omp_plan_1)
       call dfftw_destroy_plan(fft3d_info%n3_f_omp_plan_2)
       call dfftw_destroy_plan(fft3d_info%n3_b_omp_plan_2)
       call dfftw_destroy_plan(fft3d_info%n3_f_omp_plan_3)
       call dfftw_destroy_plan(fft3d_info%n3_b_omp_plan_3)

       call dfftw_destroy_plan(fft3d_info%forward_omp_plan)
       call dfftw_destroy_plan(fft3d_info%backward_omp_plan)
    end if
#endif

#endif

#ifdef ACML
    deallocate(fft3d_info%tbl,stat=ierr)
    call utils_dealloc_check('internal_serial_exit (fourier_mod.F90)', &
         'fft3d_info%tbl',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
    ierr = zfft_exit_3d(fft3d_info%iwork)
    call utils_assert(ierr == 0, 'Error in internal_serial_exit &
         &(fourier_mod.F90): zfft_exit_3d failed with code ',ierr)
#endif
#endif

    ! Workspace deallocation

    if (allocated(fft3d_info%zwork)) then
       deallocate(fft3d_info%zwork,stat=ierr)
       call utils_dealloc_check('internal_serial_exit (fourier_mod.F90)', &
            'fft3d_info%zwork',ierr)
    end if
    if (allocated(fft3d_info%dwork)) then
       deallocate(fft3d_info%dwork,stat=ierr)
       call utils_dealloc_check('internal_serial_exit (fourier_mod.F90)', &
            'fft3d_info%dwork',ierr)
    end if
    if (allocated(fft3d_info%iwork)) then
       deallocate(fft3d_info%iwork,stat=ierr)
       call utils_dealloc_check('internal_serial_exit (fourier_mod.F90)', &
            'fft3d_info%iwork',ierr)
    end if

  end subroutine internal_serial_exit

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_parallel_exit(fft3d_info)

    use utils, only : utils_dealloc_check, utils_assert
    use rundat, only: pub_threads_cellfft

    implicit none

    type(parallel_fft3d_info), intent(inout) :: fft3d_info

    ! External subroutines
#ifdef FFTW
    external :: rfftw_f77_destroy_plan, fftwnd_f77_destroy_plan
#endif
#ifdef FFTW3
    external :: dfftw_destroy_plan
#endif

    ! Local variables

    integer :: ierr             ! Error flag

    ! Platform-dependent finalisation

#ifdef FFTW
    call rfftw_f77_destroy_plan(fft3d_info%backward_plan(1))
    call fftwnd_f77_destroy_plan(fft3d_info%backward_plan(2))
    call fftwnd_f77_destroy_plan(fft3d_info%forward_plan(2))
    call rfftw_f77_destroy_plan(fft3d_info%forward_plan(1))

    deallocate(fft3d_info%rfftwbuf,stat=ierr)
    call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
         'fft3d_info%rfftwbuf',ierr)
#endif

#ifdef FFTW3
#ifdef ParallelFFT
    if (pub_threads_cellfft) then
       call dfftw_destroy_plan(fft3d_info%backward_plan_omp(1))
       call dfftw_destroy_plan(fft3d_info%backward_plan_omp(2))
       call dfftw_destroy_plan(fft3d_info%forward_plan_omp(2))
       call dfftw_destroy_plan(fft3d_info%forward_plan_omp(1))
    end if
#endif
    call dfftw_destroy_plan(fft3d_info%backward_plan(1))
    call dfftw_destroy_plan(fft3d_info%backward_plan(2))
    call dfftw_destroy_plan(fft3d_info%forward_plan(2))
    call dfftw_destroy_plan(fft3d_info%forward_plan(1))
#endif

#ifdef ACML
    deallocate(fft3d_info%rbuf,stat=ierr)
    call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
         'fft3d_info%rbuf',ierr)
    deallocate(fft3d_info%tbl,stat=ierr)
    call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
         'fft3d_info%tbl',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
    ierr = zfft_exit_2d(fft3d_info%iwork(1,2))
    call utils_assert(ierr == 0, 'Error in internal_parallel_exit &
         &(fourier_mod.F90): zfft_exit_2d failed with code ',ierr)

    ierr = dfft_exit(fft3d_info%iwork(1,1))
    call utils_assert(ierr == 0, 'Error in internal_parallel_exit &
         &(fourier_mod.F90): dfft_exit_2d failed with code ',ierr)
#endif
#endif

    ! Workspace deallocation
    if (allocated(fft3d_info%zwork)) then
       deallocate(fft3d_info%zwork,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%zwork',ierr)
    end if

    if (allocated(fft3d_info%dwork)) then
       deallocate(fft3d_info%dwork,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%dwork',ierr)
    end if

    if (allocated(fft3d_info%iwork)) then
       deallocate(fft3d_info%iwork,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%iwork',ierr)
    end if

    ! General workspace deallocation

    if (allocated(fft3d_info%brecvbuf)) then
       deallocate(fft3d_info%brecvbuf,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%brecvbuf',ierr)
    end if

    if (allocated(fft3d_info%frecvbuf)) then
       deallocate(fft3d_info%frecvbuf,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%frecvbuf',ierr)
    end if

    if (allocated(fft3d_info%bsendbuf)) then
       deallocate(fft3d_info%bsendbuf,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%bsendbuf',ierr)
    end if

    if (allocated(fft3d_info%fsendbuf)) then
       deallocate(fft3d_info%fsendbuf,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%fsendbuf',ierr)
    end if

    if (allocated(fft3d_info%buf23slab)) then
       deallocate(fft3d_info%buf23slab,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%buf23slab',ierr)
    end if

    if (allocated(fft3d_info%buf12slab)) then
       deallocate(fft3d_info%buf12slab,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%buf12slab',ierr)
    end if

    if (allocated(fft3d_info%idx23slab)) then
       deallocate(fft3d_info%idx23slab,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%idx23slab',ierr)
    end if

    if (allocated(fft3d_info%idx12slab)) then
       deallocate(fft3d_info%idx12slab,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%idx12slab',ierr)
    end if

  end subroutine internal_parallel_exit

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================!
  ! The following routines are provided for debugging purposes !
  !============================================================!
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_write_real_2(baseunit,darray)

    use comms, only: pub_my_proc_id

    implicit none

    integer, intent(in) :: baseunit
    real(kind=DP), intent(in) :: darray(:,:)

    integer :: i1,i2

    write(baseunit+pub_my_proc_id,'(i6)') 2
    write(baseunit+pub_my_proc_id,'(2i6)') size(darray,1),size(darray,2)
    do i2=1,size(darray,2)
       do i1=1,size(darray,1),2
          write(baseunit+pub_my_proc_id,'(2e24.16)') darray(i1,i2), &
               darray(i1+1,i2)
       end do
    end do

  end subroutine fourier_write_real_2

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_write_real_3(baseunit,darray)

    use comms, only: pub_my_proc_id

    implicit none

    integer, intent(in) :: baseunit
    real(kind=DP), intent(in) :: darray(:,:,:)

    integer :: i1,i2,i3

    write(baseunit+pub_my_proc_id,'(i6)') 3
    write(baseunit+pub_my_proc_id,'(2i6)') size(darray,1),size(darray,2), &
         size(darray,3)
    do i3=1,size(darray,3)
       do i2=1,size(darray,2)
          do i1=1,size(darray,1),2
             write(baseunit+pub_my_proc_id,'(2e24.16)') darray(i1,i2,i3), &
                  darray(i1+1,i2,i3)
          end do
       end do
    end do

  end subroutine fourier_write_real_3

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_write_complex_2(baseunit,zarray)

    use comms, only: pub_my_proc_id

    implicit none

    integer, intent(in) :: baseunit
    complex(kind=DP), intent(in) :: zarray(:,:)

    integer :: i1,i2

    write(baseunit+pub_my_proc_id,'(i6)') 2
    write(baseunit+pub_my_proc_id,'(2i6)') 2*size(zarray,1),size(zarray,2)
    do i2=1,size(zarray,2)
       do i1=1,size(zarray,1)
          write(baseunit+pub_my_proc_id,'(2e24.16)') real(zarray(i1,i2), &
               kind=DP),aimag(zarray(i1,i2))
       end do
    end do

  end subroutine fourier_write_complex_2

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_write_complex_3(baseunit,zarray)

    use comms, only: pub_my_proc_id

    implicit none

    integer, intent(in) :: baseunit
    complex(kind=DP), intent(in) :: zarray(:,:,:)

    integer :: i1,i2,i3

    write(baseunit+pub_my_proc_id,'(i6)') 3
    write(baseunit+pub_my_proc_id,'(2i6)') 2*size(zarray,1),size(zarray,2), &
         size(zarray,3)
    do i3=1,size(zarray,3)
       do i2=1,size(zarray,2)
          do i1=1,size(zarray,1)
             write(baseunit+pub_my_proc_id,'(2e24.16)') real(zarray(i1,i2,i3), &
                  kind=DP),aimag(zarray(i1,i2,i3))
          end do
       end do
    end do

  end subroutine fourier_write_complex_3

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  real(kind=DP) function calculate_flops(info)
    ! jd: Arguments
    type(parallel_fft3d_info), intent(in) :: info

    ! ------------------------------------------------------------------------
    calculate_flops = (2.5_DP * info%num12slabs * info%n2 * info%n1 * &
         log(real(info%n1,kind=DP)) + 5.0_DP * info%num23slabs * info%n2 * &
         info%n3 * log(real(info%n2*info%n3,kind=DP))) / log(2.0_DP)

  end function calculate_flops

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine make_license_heartbeat()

#ifdef ACCELRYS
    use license, only: license_sendHeartbeat,LIC_SUCCESS
    use constants, only: DP
    use utils, only: utils_assert
    integer :: ierr

    ! Make a license heartbeat call
    call license_sendHeartbeat(ierr)
    call utils_assert(ierr == LIC_SUCCESS, &
           'Error in fourier_apply_cell_{forward|backward}: &
            &license_sendHeartbeat failed with code ',ierr)
#endif

  end subroutine make_license_heartbeat

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef GPU_PGI
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef GPU_SP_TEST
  subroutine fourier_gpu_interp_prod_SP(in1, in2, out1, out2, &
                  gpu_sum, gpu_copy_off_current, gpu_copy_off_previous, is)
    !=========================================================================!
    ! Uses Nvidia GPUs to perform a Fast Fourier Transform from real to       !
    ! reciprocal space, upscale to a fine grid, return to real space and take !
    ! the product of the real and imaginary components of the resulting data. !
    ! Uses a combination of the PGI accelerator pragmas and CUDA Fortran.     !
    ! Input is double precision and then converted to single precision in     !
    ! to enhance performance throught the use of more cores on the GPU and the!
    ! reduction in data transfer resulting from the lower number of bits      !
    ! required to describe SP variables.
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   in1 (input): Sum of overlapping row functions.                        !
    !   in2 (input): Current col function.                                    !
    !   gpu_sum (input): Current col function is on the same atom as the      !
    !                    previous function so accumulate the density on the   !
    !                    GPU in order to minimize data transfer.              !
    !   gpu_copy_off_current (input): Current col function is the last of the !
    !                                 col functions for this batch so it needs!
    !                                 to be copied off after calculation.     !
    !   gpu_copy_off_previous (input): The previous col function was the last !
    !                                  one on its atom in this batch and needs!
    !                                  to be copied off the GPU.              !
    !   out1 (output): Array for the col function density for the previous    !
    !                  call to this routine.                                  !
    !   out2 (output): Array for the current col functions contribution to the!
    !                  density. Only used if this is the last col function of !
    !                  a batch.                                               !
    !-------------------------------------------------------------------------!
    ! Written by Karl Wilkinson 2012                                          !
    !                                                                         !
    ! Single precision routine with increased performance relative to the     !
    ! double precision equivalent. The single precision code is curently      !
    ! limited to ths routine and should be extended.                          !
    !                                                                         !
    ! In addition to obvious changes relating to the pragmas that control     !
    ! execution on the GPU the most significant differences relative to the   !
    ! CPU code are:                                                           !
    !      - Addition of CPU<>GPU data transfer steps.                        !
    !      - "Reordering" of output in order to maximise data transfer        !
    !          efficiency.                                                    !
    !      - Direct calls to the FFT routines within this routine rather      !
    !          than using a wrapper.                                          !
    ! The flags used to control the timing of data movement are controlled    !
    ! within an adapted version of the calling routine:                       !
    ! density_gpu_batch_interp_dep.                                           !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: DP, SP
    use fourier_gpu_wrapper_mod
    use utils, only: utils_assert, utils_dealloc_check
    use accel_lib
    use cudafor

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out1(:,:,:)
    real(kind=DP), intent(out) :: out2(:,:,:)

    real(kind=SP), allocatable :: in1SP(:,:,:)
    real(kind=SP), allocatable :: in2SP(:,:,:)
    real(kind=SP), allocatable :: out1SP(:,:,:)
    real(kind=SP), allocatable :: out2SP(:,:,:)


    logical, intent(in) :: gpu_sum, gpu_copy_off_current, gpu_copy_off_previous
    integer, intent(in) :: is

    ! Local variables
    integer :: n1,n2,n3
    real(kind=SP) :: scalefac
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3


    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3
    allocate(in1SP(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_interp_prod_SP','in1SP',ierr)
    allocate(in2SP(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_interp_prod_SP','in2SP',ierr)
    allocate(out1SP(2*n1,2*n2,2*n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_interp_prod_SP','out1SP',ierr)
    allocate(out2SP(2*n1,2*n2,2*n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_interp_prod_SP','out2SP',ierr)

    ! Check array sizes
    call utils_assert(&
         size(in1,1) >= n1 .and. size(in1,2) >= n2 .and. size(in1,3) >= n3, &
         'Error in fourier_gpu_interp_prod_SP: invalid dimensions for array in1')
    call utils_assert(&
         size(in2,1) >= n1 .and. size(in2,2) >= n2 .and. size(in2,3) >= n3, &
         'Error in fourier_gpu_interp_prod_SP: invalid dimensions for array in2')

    call utils_assert(&
         size(out1,1) >= 2*n1 .and. size(out1,2) >= 2*n2 .and. size(out1,3) >= 2*n3, &
         'Error in fourier_gpu_interp_prod_SP: invalid dimensions for array out1')
    call utils_assert(&
         size(out2,1) >= 2*n2 .and. size(out2,2) >= 2*n2 .and. size(out2,3) >= 2*n3, &
         'Error in fourier_gpu_interp_prod_SP: invalid dimensions for array out2')

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scalefac = 8.0 / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)

    ! kaw: Convert to SP
    in1SP = in1
    in2SP = in2

    istat = cudaMemcpyAsync(fftbox_coarse_1_gpuSP, in1SP, n1*n2*n3, stream1)
    istat = cudaMemcpyAsync(fftbox_coarse_2_gpuSP, in2SP, n1*n2*n3, stream1)

    ! kaw: Zero the fine work array.
    fine_work_gpuSP = 0.0

    ! kaw: Only proceed once the input arrays have been copied to the GPU
    istat = cudaStreamSynchronize(stream1)

    ! kaw: Copy off density from previous call to the routine if necessary
    if (gpu_copy_off_previous) then
       istat = cudaMemcpyAsync(out1SP, out_gpuSP(:,:,:,is), 8*n1*n2*n3, stream2)
       out1 = out1SP                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PROBLEM HERE
    end if

!$acc data region
!$acc kernels
    coarse_work_gpuSP(1:n1,1:n2,1:n3) = scalefac * cmplx(fftbox_coarse_1_gpuSP(1:n1,1:n2,1:n3), &
         fftbox_coarse_2_gpuSP(1:n1,1:n2,1:n3),kind=SP)
!$acc end kernels

    ! Fourier transform to reciprocal space on coarse grid
    call cufftExec(cufftplan_coarse_SP,planType_SP,coarse_work_gpuSP,coarse_work_gpuSP,CUFFT_FORWARD)

!$acc kernels
    fine_work_gpuSP(1:m1,1:m2,1:m3)          = coarse_work_gpuSP(1:m1,1:m2,1:m3)
    fine_work_gpuSP(k1:2*n1,1:m2,1:m3)       = coarse_work_gpuSP(l1:n1,1:m2,1:m3)
    fine_work_gpuSP(1:m1,k2:2*n2,1:m3)       = coarse_work_gpuSP(1:m1,l2:n2,1:m3)
    fine_work_gpuSP(k1:2*n1,k2:2*n2,1:m3)    = coarse_work_gpuSP(l1:n1,l2:n2,1:m3)
    fine_work_gpuSP(1:m1,1:m2,k3:2*n3)       = coarse_work_gpuSP(1:m1,1:m2,l3:n3)
    fine_work_gpuSP(k1:2*n1,1:m2,k3:2*n3)    = coarse_work_gpuSP(l1:n1,1:m2,l3:n3)
    fine_work_gpuSP(1:m1,k2:2*n2,k3:2*n3)    = coarse_work_gpuSP(1:m1,l2:n2,l3:n3)
    fine_work_gpuSP(k1:2*n1,k2:2*n2,k3:2*n3) = coarse_work_gpuSP(l1:n1,l2:n2,l3:n3)
!$acc end kernels

    ! Execute FFT.
    call cufftExec(cufftplan_fine_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_INVERSE)

    ! kaw: Should only proceed if the preceeding copy off is complete...
    istat = cudaStreamSynchronize(stream2)

    ! ndmh: unpack data from real and imaginary parts of complex array
    ! ndmh: and multiply them together to find product
    if (gpu_sum) then
!$acc kernels
       out_gpuSP(1:2*n1,1:2*n2,1:2*n3,is) = out_gpuSP(1:2*n1,1:2*n2,1:2*n3,is) + &
                                    (real(fine_work_gpuSP(1:2*n1,1:2*n2,1:2*n3),SP)*&
                                    aimag(fine_work_gpuSP(1:2*n1,1:2*n2,1:2*n3)))
!$acc end kernels
    else
!$acc kernels
       out_gpuSP(1:2*n1,1:2*n2,1:2*n3,is) = real(fine_work_gpuSP(1:2*n1,1:2*n2,1:2*n3),SP)*&
                                       aimag(fine_work_gpuSP(1:2*n1,1:2*n2,1:2*n3))
!$acc end kernels
    end if
!$acc end data region


    ! kaw: Copy off this density if necessary
    if (gpu_copy_off_current) then
       out2SP = out_gpuSP(:,:,:,is)
       ! kaw: convert back to DP
       out2 = out2SP
    end if


  deallocate(in1SP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_interp_prod_SP','in1SP',ierr)
  deallocate(in2SP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_interp_prod_SP','in2SP',ierr)

  end subroutine fourier_gpu_interp_prod_SP

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#ifdef GPU_SP_TEST_POT
  subroutine fourier_gpu_interpolate_SP(in1, in2, potential_box_dbl_pair, apply_pot_2)
    !=========================================================================!
    ! Uses Nvidia GPUs to perform a Fast Fourier Transform from real to       !
    ! reciprocal space, upscale to a fine grid, return to real space and take !
    ! the product of the resulting data and the corrosponding region of the   !
    ! potential as extracted from the simulation cell.                        !
    !                                                                         !
    ! Uses a combination of the PGI accelerator pragmas and CUDA Fortran.     !
    !                                                                         !
    ! Input is double precision and then converted to single precision in     !
    ! to enhance performance throught the use of more cores on the GPU and the!
    ! reduction in data transfer resulting from the lower number of bits      !
    ! required to describe SP variables.                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   in1 (input): Current col function.                                    !
    !   in2 (input): An optional second col function (usually present).       !
    !   potential_box_dbl_pair (input): Region of the potential extracted from!
    !           the simulation cell that is to be applied to the col          !
    !           function(s).                                                  !
    !   apply_pot_2 (input): Indicates the presence of the second col function!
    !-------------------------------------------------------------------------!
    ! Written by Karl Wilkinson 2012                                          !
    !                                                                         !
    ! In addition to obvious changes relating to the pragmas that control     !
    ! execution on the GPU the most significant differences relative to the   !
    ! CPU code are:                                                           !
    !      - Addition of CPU<>GPU data transfer steps.                        !
    !      - Direct calls to the FFT routines within this routine rather      !
    !          than using a wrapper.                                          !
    !      - Calculation of col function * potential within this routine      !
    !          rather than within the calling routine.                        !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: DP, SP
    use fourier_gpu_wrapper_mod
    use utils, only: utils_assert
    use accel_lib !! External dependency
    use cudafor !! External dependency

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(in) :: potential_box_dbl_pair(:,:,:,:)
    real(kind=SP), allocatable :: in1SP(:,:,:)
    real(kind=SP), allocatable :: in2SP(:,:,:)
    real(kind=SP), allocatable :: potential_box_dbl_pairSP(:,:,:,:)

    ! Local variables
    integer :: n1,n2,n3,i2,i1
    integer :: n1fine, n2fine, n3fine
    real(kind=SP) :: scalefac
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3,xj1,xj2,xj3
    logical :: apply_pot_2
    integer :: ierr

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3
    allocate(in1SP(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_interpolate_SP','in1SP',ierr)
    allocate(in2SP(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_interpolate_SP','in2SP',ierr)

    allocate(potential_box_dbl_pairSP(2*n1,2*n2,2*n3,2),stat=ierr)
    call utils_alloc_check('fourier_gpu_interpolate_SP','potential_box_dbl_pairSP',ierr)

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scalefac = 8.0 / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)

    ! kaw: Convert inputs to SP
    in1SP = in1
    in2SP = in2
    potential_box_dbl_pairSP = potential_box_dbl_pair

    istat = cudaMemcpyAsync(fftbox_coarse_1_gpuSP, in1SP, n1*n2*n3, stream1)
    istat = cudaMemcpyAsync(fftbox_coarse_2_gpuSP, in2SP, n1*n2*n3, stream2)

    ! kaw: Zero the fine work array. This is automatically done on the GPU as it has been declared as
    !      a CUDA Fortran array (device clause in declaration)
    fine_work_gpuSP = 0.0

    ! kaw: Only proceed once the input arrays have been copied to the GPU
    istat = cudaStreamSynchronize(stream1)
    istat = cudaStreamSynchronize(stream2)

    istat = cudaMemcpyAsync(potential_box_dbl_gpuSP, potential_box_dbl_pairSP, 16*n1*n2*n3, stream1)
    ! kaw: Transfer potential data to GPU

!$acc data region
!$acc kernels
    coarse_work_gpuSP(1:n1,1:n2,1:n3) = scalefac * cmplx(fftbox_coarse_1_gpuSP(1:n1,1:n2,1:n3), &
         fftbox_coarse_2_gpuSP(1:n1,1:n2,1:n3),kind=SP)
!$acc end kernels

    ! Fourier transform to reciprocal space on coarse grid
    call cufftExec(cufftplan_coarse_SP,planType_SP,coarse_work_gpuSP,coarse_work_gpuSP,CUFFT_FORWARD)

!$acc kernels
             fine_work_gpuSP(1:m1,1:m2,1:m3)    = coarse_work_gpuSP(1:m1,1:m2,1:m3)
          fine_work_gpuSP(k1:2*n1,1:m2,1:m3)   = coarse_work_gpuSP(l1:n1,1:m2,1:m3)
          fine_work_gpuSP(1:m1,k2:2*n2,1:m3)   = coarse_work_gpuSP(1:m1,l2:n2,1:m3)
       fine_work_gpuSP(k1:2*n1,k2:2*n2,1:m3)  = coarse_work_gpuSP(l1:n1,l2:n2,1:m3)
          fine_work_gpuSP(1:m1,1:m2,k3:2*n3)   = coarse_work_gpuSP(1:m1,1:m2,l3:n3)
       fine_work_gpuSP(k1:2*n1,1:m2,k3:2*n3)  = coarse_work_gpuSP(l1:n1,1:m2,l3:n3)
       fine_work_gpuSP(1:m1,k2:2*n2,k3:2*n3)  = coarse_work_gpuSP(1:m1,l2:n2,l3:n3)
    fine_work_gpuSP(k1:2*n1,k2:2*n2,k3:2*n3) = coarse_work_gpuSP(l1:n1,l2:n2,l3:n3)
!$acc end kernels

    ! Execute FFT.
!    call cufftExec(cufftplan_fine_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_INVERSE)

    n1fine = n1*2
    n2fine = n2*2
    n3fine = n3*2
    j1 = 3*n1fine/4+2
    ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
    do i2=1,n2fine/4+1
        call cufftExec(cufft_n3_2_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_INVERSE,(n1fine*(i2-1)))
        call cufftExec(cufft_n3_3_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_INVERSE,(n1fine*(i2-1))+j1-1)
    end do
    ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
    do i2=3*n2fine/4+2,n2fine
        call cufftExec(cufft_n3_2_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_INVERSE,(n1fine*(i2-1)))
        call cufftExec(cufft_n3_3_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_INVERSE,(n1fine*(i2-1))+j1-1)
    end do
    ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    do i1=1,n1fine/4+1
       call cufftExec(cufft_n2_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_INVERSE,i1-1)
    end do
    ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    do i1=3*n1fine/4+2,n1fine
       call cufftExec(cufft_n2_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_INVERSE,i1-1)
    end do
    ! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
    call cufftExec(cufft_n1_SP, planType_SP, fine_work_gpuSP, fine_work_gpuSP, CUFFT_INVERSE)

    ! ndmh: unpack data from real and imaginary parts of complex array
    ! ndmh: and multiply them together to find product
    ! -------- Apply potential for NGWFs as needed ------------
    if (apply_pot_2) then
       !$acc kernels
       do k3 = 1, 2*n3
          do k2 = 1, 2*n2
             fftbox_fine_1_gpuSP(:,k2,k3) = potential_box_dbl_gpuSP(:,k2,k3,1)*real(fine_work_gpuSP(:,k2,k3),SP)
             fftbox_fine_2_gpuSP(:,k2,k3) = potential_box_dbl_gpuSP(:,k2,k3,2)*aimag(fine_work_gpuSP(:,k2,k3))
             !fine_work_gpu(:,k2,k3) = cmplx(potential_box_dbl_gpu(:,k2,k3,1)*real(fine_work_gpu(:,k2,k3),
             !                               potential_box_dbl_gpu(:,k2,k3,2)*aimag(fine_work_gpu(:,k2,k3))))
          end do
       end do
       !$acc end kernels
    else
       !$acc kernels
       do k3 = 1, 2*n3
          do k2 = 1, 2*n2
             fftbox_fine_1_gpuSP(:,k2,k3) = potential_box_dbl_gpuSP(:,k2,k3,1)*real(fine_work_gpuSP(:,k2,k3),SP)
             !fine_work_gpu(:,k2,k3) = cmplx(potential_box_dbl_gpu(:,k2,k3,1)*real(fine_work_gpu(:,k2,k3),
             !                               aimag(fine_work_gpu(:,k2,k3))))
          end do
       end do
       !$acc end kernels
    end if

    ! kaw: Only proceed once the input arrays have been copied to the GPU
    istat = cudaStreamSynchronize(stream1)

    deallocate(in1SP,stat=ierr)
    call utils_dealloc_check('fourier_gpu_interpolate_SP','in1SP',ierr)
    deallocate(in2SP,stat=ierr)
    call utils_dealloc_check('fourier_gpu_interpolate_SP','in2SP',ierr)
    deallocate(potential_box_dbl_pairSP,stat=ierr)
    call utils_dealloc_check('fourier_gpu_interpolate_SP','potential_box_dbl_pairSP',ierr)

  end subroutine fourier_gpu_interpolate_SP

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_gpu_filter_SP(batch_count,multibatch)
    !=========================================================================!
    ! Uses Nvidia GPUs to perform a Fast Fourier Transform from real to       !
    ! reciprocal space, downscales to a coarse grid then returns to real space!
    !                                                                         !
    ! Uses a combination of the PGI accelerator pragmas and CUDA Fortran.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   batch_count (input): Which col function of the batch is being         !
    !           processed.                                                    !
    !   multibatch (input): Indicates the presence of the second col function !
    !-------------------------------------------------------------------------!
    ! Written by Karl Wilkinson 2012                                          !
    !                                                                         !
    ! In addition to obvious changes relating to the pragmas that control     !
    ! execution on the GPU the most significant differences relative to the   !
    ! CPU code are:                                                           !
    !      - Direct calls to the FFT routines within this routine rather      !
    !          than using a wrapper.                                          !
    !      - Removal of the out1 and out2 output arrays. These are dealt with !
    !          at a later stage in the calculation to allow the local         !
    !          potential integrals to be calculated on the GPU itself         !
    !          (depending upon) the relative data transfer costs.             !
    ! This routine is called from potential_gpu_app2_ngwf_batch and the   !
    !
    !=========================================================================!

    use constants, only: DP, SP
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor !! External dependency

    implicit none

    ! Arguments
    integer,intent(in) :: batch_count
    logical,intent(in) :: multibatch

    ! Local variables
    integer :: n1,n2,n3
    integer :: n1fine,n2fine,n3fine,i1,i2,j1
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3
    real(kind=SP) :: scale

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3


    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1 + l1 ; k2 = n2 + l2 ; k3 = n3 + l3

    ! Pack data into real and imaginary parts of complex array
!$acc data region
!$acc kernels
  fine_work_gpuSP(1:2*n1,1:2*n2,1:2*n3) = cmplx(fftbox_fine_1_gpuSP(1:2*n1,1:2*n2,1:2*n3), &
                                          fftbox_fine_2_gpuSP(1:2*n1,1:2*n2,1:2*n3),SP)
!$acc end kernels

    ! Fourier transform to reciprocal space on fine grid
    call cufftExec(cufftplan_fine_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_FORWARD)
    !n1fine = n1*2
    !n2fine = n2*2
    !n3fine = n3*2
    !! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
    !call cufftExec(cufft_n1_SP, planType_SP, fine_work_gpuSP, fine_work_gpuSP, CUFFT_FORWARD)
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=1,n1fine/4+1
    !   call cufftExec(cufft_n2_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_FORWARD,i1-1)
    !end do
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=3*n1fine/4+2,n1fine
    !   call cufftExec(cufft_n2_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_FORWARD,i1-1)
    !end do
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
    !j1 = 3*n1fine/4+2
    !do i2=1,n2fine/4+1
    !    call cufftExec(cufft_n3_2_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_FORWARD,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_FORWARD,(n1fine*(i2-1))+j1-1)
    !end do
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
    !do i2=3*n2fine/4+2,n2fine
    !    call cufftExec(cufft_n3_2_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_FORWARD,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3_SP,planType_SP,fine_work_gpuSP,fine_work_gpuSP,CUFFT_FORWARD,(n1fine*(i2-1))+j1-1)
    !end do


    ! Copy Fourier components to coarse reciprocal space grid
!$acc kernels
    coarse_work_gpuSP(1:m1,1:m2,1:m3) =    fine_work_gpuSP(1:m1,1:m2,1:m3)
    coarse_work_gpuSP(l1:n1,1:m2,1:m3) =   fine_work_gpuSP(k1:2*n1,1:m2,1:m3)
    coarse_work_gpuSP(1:m1,l2:n2,1:m3) =   fine_work_gpuSP(1:m1,k2:2*n2,1:m3)
    coarse_work_gpuSP(l1:n1,l2:n2,1:m3) =  fine_work_gpuSP(k1:2*n1,k2:2*n2,1:m3)
    coarse_work_gpuSP(1:m1,1:m2,l3:n3) =   fine_work_gpuSP(1:m1,1:m2,k3:2*n3)
    coarse_work_gpuSP(l1:n1,1:m2,l3:n3) =  fine_work_gpuSP(k1:2*n1,1:m2,k3:2*n3)
    coarse_work_gpuSP(1:m1,l2:n2,l3:n3) =  fine_work_gpuSP(1:m1,k2:2*n2,k3:2*n3)
    coarse_work_gpuSP(l1:n1,l2:n2,l3:n3) = fine_work_gpuSP(k1:2*n1,k2:2*n2,k3:2*n3)
!$acc end kernels

    ! Fourier transform to real space on coarse grid
    call cufftExec(cufftplan_coarse_SP,planType_SP,coarse_work_gpuSP,coarse_work_gpuSP,CUFFT_INVERSE)

    ! Normalise
    scale = 0.125_SP / (internal_box_coarse_info%ld1*&
         internal_box_coarse_info%ld2*internal_box_coarse_info%n3)

    ! Unpack data from real and imaginary parts of complex array
    ! ndmh: a loop saves on cache misses when repeating for aimag
    if (multibatch) then
       !$acc kernels
       do k3 = 1,n3
          do k2 = 1,n2
             do k1 = 1,n1
                potential_fftbox_batch_gpuSP(k1,k2,k3,batch_count-1) = real(scale*coarse_work_gpuSP(k1,k2,k3),SP)
                potential_fftbox_batch_gpuSP(k1,k2,k3,batch_count) = aimag(scale*coarse_work_gpuSP(k1,k2,k3))
             end do
          end do
       end do
       !$acc end kernels
    else
       !$acc kernels
       do k3 = 1,n3
          do k2 = 1,n2
             do k1 = 1,n1
                potential_fftbox_batch_gpuSP(k1,k2,k3,batch_count) = real(scale*coarse_work_gpuSP(k1,k2,k3),SP)
             end do
          end do
       end do
       !$acc end kernels
    end if
!$acc end data region

  end subroutine fourier_gpu_filter_SP
#endif
  !kaw: End of GPU_SP_TEST_POT ifdef

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#endif

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  subroutine fourier_gpu_interp_prod(in1, in2, out1, out2, &
                  gpu_sum, gpu_copy_off_current, gpu_copy_off_previous, is)
    !=========================================================================!
    ! Uses Nvidia GPUs to perform a Fast Fourier Transform from real to       !
    ! reciprocal space, upscale to a fine grid, return to real space and take !
    ! the product of the real and imaginary components of the resulting data. !
    ! Uses a combination of the PGI accelerator pragmas and CUDA Fortran.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   in1 (input): Sum of overlapping row functions.                        !
    !   in2 (input): Current col function.                                    !
    !   gpu_sum (input): Current col function is on the same atom as the      !
    !           previous function so accumulate the density on the GPU in     !
    !           order to minimize data transfer.                              !
    !   gpu_copy_off_current (input): Current col function is the last of the !
    !           col functions for this batch so it needs to be copied off     !
    !           after calculation.                                            !
    !   gpu_copy_off_previous (input): The previous col function was the last !
    !           one on its atom in this batch and needs to be copied off the  !
    !           GPU.                                                          !
    !   out1 (output): Array for the col function density for the previous    !
    !           call to this routine.                                         !
    !   out2 (output): Array for the current col functions contribution to the!
    !           density. Only used if this is the last col function of a batch!
    !-------------------------------------------------------------------------!
    ! Written by Karl Wilkinson 2012                                          !
    ! heavily edited by Karl Wilkinson and Levi Barnes in 2014/2015 to        !
    ! substantially improve performance through the use of CUDA streams for   !
    ! the concurrent use of multiple execution and data transfer engines.     !
    !                                                                         !
    ! In addition to obvious changes relating to the pragmas that control     !
    ! execution on the GPU the most significant differences relative to the   !
    ! CPU code are:                                                           !
    !      - Addition of CPU<>GPU data transfer steps.                        !
    !      - "Reordering" of output in order to maximise data transfer        !
    !          efficiency.                                                    !
    !      - Direct calls to the FFT routines within this routine rather      !
    !          than using a wrapper.                                          !
    ! The flags used to control the timing of data movement are controlled    !
    ! within an adapted version of the calling routine:                       !
    ! density_gpu_batch_interp_dep.                                           !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: DP, I4B
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor !! External dependency
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out1(:,:,:)
    real(kind=DP), intent(out) :: out2(:,:,:)
    logical, intent(in) :: gpu_sum, gpu_copy_off_current, gpu_copy_off_previous
    integer, intent(in) :: is
    type(DIM3) :: grid, block
    integer :: dim1

! Local variables
    integer :: n1,n2,n3
    real(kind=DP) :: scalefac
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3,xj1,xj2,xj3
    integer :: i,j,k,z

    integer :: istat
    ! Stream variables
    integer(kind=I4B) :: streamA, streamB

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scalefac = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)

!Start of OpenACC data region
!$acc data region

    ! kaw: Zero the fine work array. This is automatically done on the GPU as it has been declared as
    !      a device resident array. This should be streamed...
    fine_work_gpu = 0.0_DP

    ! kaw: Set up the cudastreams
    streamA = 4
    streamB = 5
    istat = cudaMemcpyAsync(fftbox_coarse_1_gpu, in1, n1*n2*n3, stream1)
    istat = cudaMemcpyAsync(fftbox_coarse_2_gpu, in2, n1*n2*n3, stream1)

    ! kaw: Copy off density from previous call to the routine if necessary
    if (gpu_copy_off_previous) then
       istat = cudaMemcpyAsync(out1, out_gpu(:,:,:,is), 8*n1*n2*n3, stream2)
    end if

    istat = cudaStreamSynchronize(stream1)
    !$acc kernels async(streamA)
    coarse_work_gpu(1:n1,1:n2,1:n3) = scalefac * cmplx(fftbox_coarse_1_gpu(1:n1,1:n2,1:n3), &
         fftbox_coarse_2_gpu(1:n1,1:n2,1:n3),kind=DP)
    !$acc end kernels

    ! Fourier transform to reciprocal space on coarse grid
    !$acc wait(streamA)
    call cufftSetStream(cufftplan_coarse, stream1)
    call cufftExec(cufftplan_coarse,planType,coarse_work_gpu,coarse_work_gpu,CUFFT_FORWARD)
    call cufftSetStream(cufftplan_coarse, 0)

    ! kaw: Should only proceed if the fine work array has been zeroed
    istat = cudaStreamSynchronize(stream1)

    !$acc parallel do async(streamA)
    do z=0,m1*m2*m3-1
       i = mod(z,m1)
       j = mod(z/m1,m2)
       k = z/m1/m2
       fine_work_gpu(1+i,1+j,1+k) = coarse_work_gpu(1+i,1+j,1+k)
       if(l1+i.le.n1) fine_work_gpu(k1+i,1+j,1+k) = coarse_work_gpu(l1+i,1+j,1+k)
       if(l2+j.le.n2) fine_work_gpu(1+i,k2+j,1+k) = coarse_work_gpu(1+i,l2+j,1+k)
       if(l1+i.le.n1.and.l2+j.le.n2) fine_work_gpu(k1+i,k2+j,1+k) = coarse_work_gpu(l1+i,l2+j,1+k)
       if(l3+k.le.n3) fine_work_gpu(1+i,1+j,k3+k) = coarse_work_gpu(1+i,1+j,l3+k)
       if(l1+i.le.n1.and.l3+k.le.n3) fine_work_gpu(k1+i,1+j,k3+k) = coarse_work_gpu(l1+i,1+j,l3+k)
       if(l2+j.le.n2.and.l3+k.le.n3) fine_work_gpu(1+i,k2+j,k3+k) = coarse_work_gpu(1+i,l2+j,l3+k)
       if(l1+i.le.n1.and.l2+j.le.n2.and.l3+k.le.n3) fine_work_gpu(k1+i,k2+j,k3+k) = coarse_work_gpu(l1+i,l2+j,l3+k)
    end do

    ! Execute FFT.
    !$acc wait(streamA)
    call cufftSetStream(cufftplan_fine, stream1)
    call cufftExec(cufftplan_fine,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE)
    call cufftSetStream(cufftplan_fine, 0)

    ! kaw: Should only proceed if the preceeding copy off is complete...
    istat = cudaStreamSynchronize(stream1)
    istat = cudaStreamSynchronize(stream2)

    ! ndmh: unpack data from real and imaginary parts of complex array
    ! ndmh: and multiply them together to find product
    if (gpu_sum) then
       !$acc kernels async(streamA)
       out_gpu(1:2*n1,1:2*n2,1:2*n3,is) = out_gpu(1:2*n1,1:2*n2,1:2*n3,is) + &
                                   (real(fine_work_gpu(1:2*n1,1:2*n2,1:2*n3),DP)*&
                                    real(cmplx(0.0d0, -1.0d0, DP)*fine_work_gpu(1:2*n1,1:2*n2,1:2*n3)))
!$acc end kernels
    else
       !$acc kernels async(streamA)
       out_gpu(1:2*n1,1:2*n2,1:2*n3,is) = real(fine_work_gpu(1:2*n1,1:2*n2,1:2*n3),DP)*&
                                         real(cmplx(0.0d0, -1.0d0, DP)*fine_work_gpu(1:2*n1,1:2*n2,1:2*n3))
!$acc end kernels
    end if


!$acc end data region
!End of OpenACC data region

    ! kaw: Synchronise streams.
    !$acc wait(streamA)

    ! kaw: Copy off this density if necessary, this cannot be done asynchronously.
    if (gpu_copy_off_current) then
       out2 = out_gpu(:,:,:,is)
    end if


  end subroutine fourier_gpu_interp_prod

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_gpu_interpolate(in1, in2, potential_box_dbl_pair, apply_pot_2)
    !=========================================================================!
    ! Uses Nvidia GPUs to perform a Fast Fourier Transform from real to       !
    ! reciprocal space, upscale to a fine grid, return to real space and take !
    ! the product of the resulting data and the corrosponding region of the   !
    ! potential as extracted from the simulation cell.                        !
    !                                                                         !
    ! Uses a combination of the PGI accelerator pragmas and CUDA Fortran.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   in1 (input): Current col function.                                    !
    !   in2 (input): An optional second col function (usually present).       !
    !   potential_box_dbl_pair (input): Region of the potential extracted from!
    !           the simulation cell that is to be applied to the col          !
    !           function(s).                                                  !
    !   apply_pot_2 (input): Indicates the presence of the second col function!
    !-------------------------------------------------------------------------!
    ! Written by Karl Wilkinson 2012                                          !
    !                                                                         !
    ! In addition to obvious changes relating to the pragmas that control     !
    ! execution on the GPU the most significant differences relative to the   !
    ! CPU code are:                                                           !
    !      - Addition of CPU<>GPU data transfer steps.                        !
    !      - Direct calls to the FFT routines within this routine rather      !
    !          than using a wrapper.                                          !
    !      - Calculation of col function * potential within this routine      !
    !          rather than within the calling routine.                        !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: DP, I4B
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor !! External dependency
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(in) :: potential_box_dbl_pair(:,:,:,:)

    ! Local variables
    integer :: n1,n2,n3,i2,i1
    integer :: n1fine, n2fine, n3fine
    real(kind=DP) :: scalefac
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3
    logical :: apply_pot_2
    integer :: istat
    integer :: z,i,j,k,batch_count
    real(kind=DP) :: arr_real,arr_imag

    ! Stream variables for PGI streams (CUDA streams are initialised in gpu_init)
    integer(kind=I4B) :: streamA, streamB

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3
    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)



    ! ndmh: and normalise in advance as it uses fewer operations
    scalefac = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)

!$acc data region


    ! kaw: Set up the cudastreams
    streamA = 4
    streamB = 5

    ! kaw: Copy data from host to device
 !Potential optimisation would skip this intermediate step and move directly
 !from a CPU array to the real and imaginary components of a complex array.
 !Alternatively, the complexation could be done on the CPU and the data transfer
 !could simply be that of a complex array.
    istat = cudaMemcpyAsync(fftbox_coarse_1_gpu, in1, n1*n2*n3, stream1)
    istat = cudaMemcpyAsync(fftbox_coarse_2_gpu, in2, n1*n2*n3, stream2)


    ! kaw: zero the fine grid array
    !$acc kernels async(streamB)
    do z=1,8*n1*n2*n3
       fine_work_gpu(z,1,1) = 0.0d0
    enddo
   !$acc end kernels


    istat = cudaStreamSynchronize(stream1)
    istat = cudaStreamSynchronize(stream2)


    ! Pack data into real and imaginary parts of coarse, complex array
    !$acc kernels async(streamA)
    coarse_work_gpu(1:n1,1:n2,1:n3) = scalefac * cmplx(fftbox_coarse_1_gpu(1:n1,1:n2,1:n3), &
         fftbox_coarse_2_gpu(1:n1,1:n2,1:n3),kind=DP)
    !$acc end kernels


    ! kaw: Transfer potential data to GPU
    istat = cudaMemcpyAsync(potential_box_dbl_gpu, potential_box_dbl_pair, 16*n1*n2*n3, stream2)
    !$acc wait(streamA)


    ! Fourier transform to reciprocal space on coarse grid
    call cufftSetStream(cufftplan_coarse, stream1)
    call cufftExec(cufftplan_coarse,planType,coarse_work_gpu,coarse_work_gpu,CUFFT_FORWARD)
    call cufftSetStream(cufftplan_coarse, 0)


    istat = cudaStreamSynchronize(stream1)
    !$acc wait(streamB)


    ! Populate fine grid array
    !$acc parallel do async(streamA)
    do z=0,m1*m2*m3-1
       i = mod(z,m1)
       j = mod(z/m1,m2)
       k = (z/m1)/m2
       fine_work_gpu(1+i,1+j,1+k) = coarse_work_gpu(1+i,1+j,1+k)
       if(l1+i.le.n1) fine_work_gpu(k1+i,1+j,1+k) = coarse_work_gpu(l1+i,1+j,1+k)
       if(l2+j.le.n2) fine_work_gpu(1+i,k2+j,1+k) = coarse_work_gpu(1+i,l2+j,1+k)
       if(l1+i.le.n1.and.l2+j.le.n2) fine_work_gpu(k1+i,k2+j,1+k) = coarse_work_gpu(l1+i,l2+j,1+k)
       if(l3+k.le.n3) fine_work_gpu(1+i,1+j,k3+k) = coarse_work_gpu(1+i,1+j,l3+k)
       if(l1+i.le.n1.and.l3+k.le.n3) fine_work_gpu(k1+i,1+j,k3+k) = coarse_work_gpu(l1+i,1+j,l3+k)
       if(l2+j.le.n2.and.l3+k.le.n3) fine_work_gpu(1+i,k2+j,k3+k) = coarse_work_gpu(1+i,l2+j,l3+k)
       if(l1+i.le.n1.and.l2+j.le.n2.and.l3+k.le.n3) fine_work_gpu(k1+i,k2+j,k3+k) = coarse_work_gpu(l1+i,l2+j,l3+k)
    end do
    !$acc end parallel

    !$acc wait(streamA)


    ! Execute inverse FFT.
    call cufftSetStream(cufftplan_fine, stream1)
    call cufftExec(cufftplan_fine,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE)
    call cufftSetStream(cufftplan_fine, 0)


!   ! kaw: In theory the following should be faster than the 3D FFT above.
!   !      This is not the case in practice but will be revisited with input
!   !      from NVIDIA.
    !n1fine = n1*2
    !n2fine = n2*2
    !n3fine = n3*2
    !j1 = 3*n1fine/4+2
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
    !do i2=1,n2fine/4+1
    !    call cufftExec(cufft_n3_2,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,(n1fine*(i2-1))+j1-1)
    !end do
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
    !do i2=3*n2fine/4+2,n2fine
    !    call cufftExec(cufft_n3_2,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,(n1fine*(i2-1))+j1-1)
    !end do
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=1,n1fine/4+1
    !   call cufftExec(cufft_n2,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,i1-1)
    !end do
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=3*n1fine/4+2,n1fine
    !   call cufftExec(cufft_n2,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,i1-1)
    !end do
    !! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
    !call cufftExec(cufft_n1, planType, fine_work_gpu, fine_work_gpu, CUFFT_INVERSE)


    ! ndmh: unpack data from real and imaginary parts of complex array
    ! ndmh: and multiply them together to find product


    ! -------- Apply potential to NGWFs as needed ------------
    !Synchronise all streams
    istat = cudaStreamSynchronize(stream2)
    istat = cudaStreamSynchronize(stream1)


    if (apply_pot_2) then
       !$acc parallel loop collapse(3) async(streamA)
       do k3 = 1,2*n3
          do k2 = 1,2*n2
             do k1 = 1,2*n1
                   arr_real = potential_box_dbl_gpu(k1,k2,k3,1)*real(fine_work_gpu(k1,k2,k3),DP)
                   arr_imag = potential_box_dbl_gpu(k1,k2,k3,2)*aimag(fine_work_gpu(k1,k2,k3))
                   fine_work_gpu(k1,k2,k3) = cmplx(arr_real,arr_imag,DP)
             end do
          end do
       end do
       !$acc end parallel loop
    else
       !$acc parallel loop collapse(3) async(streamA)
       do k3 = 1, 2*n3
          do k2 = 1, 2*n2
             do k1 = 1, 2*n1
                arr_real = potential_box_dbl_gpu(k1,k2,k3,1)*real(fine_work_gpu(k1,k2,k3),DP)
                arr_imag = 0.0d0
                fine_work_gpu(k1,k2,k3) = cmplx(arr_real,arr_imag,DP)
             end do
          end do
       end do
       !$acc end parallel loop
    end if


!$acc end data region


    ! kaw: Synchronise streams before exiting.
    !$acc wait(streamA)

  end subroutine fourier_gpu_interpolate

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_gpu_filter(batch_count,multibatch)
    !=========================================================================!
    ! Uses Nvidia GPUs to perform a Fast Fourier Transform from real to       !
    ! reciprocal space, downscales to a coarse grid then returns to real space!
    !                                                                         !
    ! Uses a combination of the PGI accelerator pragmas and CUDA Fortran.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   batch_count (input): Which col function of the batch is being         !
    !           processed.                                                    !
    !   multibatch (input): Indicates the presence of the second col function !
    !-------------------------------------------------------------------------!
    ! Written by Karl Wilkinson 2012                                          !
    !                                                                         !
    ! In addition to obvious changes relating to the pragmas that control     !
    ! execution on the GPU the most significant differences relative to the   !
    ! CPU code are:                                                           !
    !      - Direct calls to the FFT routines within this routine rather      !
    !          than using a wrapper.                                          !
    !      - Removal of the out1 and out2 output arrays. These are dealt with !
    !          at a later stage in the calculation to allow the local         !
    !          potential integrals to be calculated on the GPU itself         !
    !          (depending upon) the relative data transfer costs.             !
    ! This routine is called from potential_gpu_app2_ngwf_batch and the   !
    !
    !=========================================================================!

    use constants, only: DP
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor !! External dependency
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer,intent(in) :: batch_count
    logical,intent(in) :: multibatch

    ! Local variables
    integer :: n1,n2,n3
    integer :: n1fine,n2fine,n3fine,i1,i2,j1
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3
    real(kind=DP) :: scale
    integer :: istat
    integer :: i,j,k,z

    ! FFT box dimensions and bounds for converting from fine to coarse grids
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3
    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1 + l1 ; k2 = n2 + l2 ; k3 = n3 + l3


!$acc data region


    ! Fourier transform to reciprocal space on fine grid
    call cufftSetStream(cufftplan_coarse, stream1)
    call cufftExec(cufftplan_fine,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD)
    call cufftSetStream(cufftplan_coarse, 0)

    istat = cudaDeviceSynchronize()
    !This is critical, although I am not clear why. It should be possible to remove any
    !references to streams in this routine but even if this is done for the above call
    !to cufftExec a synchronisation issue emerges.
    !Even if streams are used, it should be possible to synchronise just stream1, but
    !this doesn't work so something more complex is going on.
    !Streams are definitely the issue as forcing sequential execution with environment
    !variables fixes the problem.


!   ! kaw: In theory the following should be faster than the 3D FFT above.
!   !      This is not the case in practice but will be revisited with input
!   !      from NVIDIA.
    !n1fine = n1*2
    !n2fine = n2*2
    !n3fine = n3*2
    !! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
    !call cufftExec(cufft_n1, planType, fine_work_gpu, fine_work_gpu, CUFFT_FORWARD)
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=1,n1fine/4+1
    !   call cufftExec(cufft_n2,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,i1-1)
    !end do
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=3*n1fine/4+2,n1fine
    !   call cufftExec(cufft_n2,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,i1-1)
    !end do
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
    !j1 = 3*n1fine/4+2
    !do i2=1,n2fine/4+1
    !    call cufftExec(cufft_n3_2,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,(n1fine*(i2-1))+j1-1)
    !end do
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
    !do i2=3*n2fine/4+2,n2fine
    !    call cufftExec(cufft_n3_2,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,(n1fine*(i2-1))+j1-1)
    !end do


    ! Copy Fourier components to coarse reciprocal space grid
    !$acc parallel do
    do z=0,m1*m2*m3-1
       i = mod(z,m1)
       j = mod(z/m1,m2)
       k = z/m1/m2
       coarse_work_gpu(1+i,1+j,1+k) = fine_work_gpu(1+i,1+j,1+k)
       if (l1+i.le.n1) coarse_work_gpu(l1+i,1+j,1+k) = fine_work_gpu(k1+i,1+j,1+k)
       if (l2+j.le.n2) coarse_work_gpu(1+i,l2+j,1+k) = fine_work_gpu(1+i,k2+j,1+k)
       if (l1+i.le.n1.and.l2+j.le.n2) coarse_work_gpu(l1+i,l2+j,1+k) = fine_work_gpu(k1+i,k2+j,1+k)
       if (l3+k.le.n3) coarse_work_gpu(1+i,1+j,l3+k) = fine_work_gpu(1+i,1+j,k3+k)
       if (l1+i.le.n1.and.l3+k.le.n3) coarse_work_gpu(l1+i,1+j,l3+k) = fine_work_gpu(k1+i,1+j,k3+k)
       if (l2+j.le.n2.and.l3+k.le.n3) coarse_work_gpu(1+i,l2+j,l3+k) = fine_work_gpu(1+i,k2+j,k3+k)
       if (l1+i.le.n1.and.l2+j.le.n2.and.l3+k.le.n3)coarse_work_gpu(l1+i,l2+j,l3+k) = fine_work_gpu(k1+i,k2+j,k3+k)
    end do
    !$acc end parallel


    ! Fourier transform to real space on coarse grid
    call cufftSetStream(cufftplan_coarse, stream1)
    call cufftExec(cufftplan_coarse,planType,coarse_work_gpu,coarse_work_gpu,CUFFT_INVERSE)
    call cufftSetStream(cufftplan_coarse, 0)

    istat = cudaDeviceSynchronize()
    !Again, this is critical, although I am not clear why. see above for further notes on this issue.


    ! Normalise
    scale = 0.125_DP / (internal_box_coarse_info%ld1*&
         internal_box_coarse_info%ld2*internal_box_coarse_info%n3)


    ! Unpack data from real and imaginary parts of complex array
    ! ndmh: a loop saves on cache misses when repeating for aimag
    if (multibatch) then
       !$acc parallel loop collapse(3)
       do k3 = 1,n3
          do k2 = 1,n2
             do k1 = 1,n1
                potential_fftbox_batch_gpu(k1,k2,k3,batch_count-1) = real(scale*coarse_work_gpu(k1,k2,k3),DP)
                potential_fftbox_batch_gpu(k1,k2,k3,batch_count) = aimag(scale*coarse_work_gpu(k1,k2,k3))
             end do
          end do
       end do
       !$acc end parallel loop
    else
       !$acc parallel loop collapse(3)
       do k3 = 1,n3
          do k2 = 1,n2
             do k1 = 1,n1
                potential_fftbox_batch_gpu(k1,k2,k3,batch_count) = real(scale*coarse_work_gpu(k1,k2,k3),DP)
             end do
          end do
       end do
       !$acc end parallel loop
    end if
!$acc end data region

  end subroutine fourier_gpu_filter


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_gpu_interpolate_grad(potential_box_dbl_pair, apply_pot_2, &
                                          is, batch_count)
    !=========================================================================!
    ! Uses Nvidia GPUs to perform a Fast Fourier Transform from real to       !
    ! reciprocal space, upscale to a fine grid, return to real space and take !
    ! the product of the resulting data and the corrosponding region of the   !
    ! potential as extracted from the simulation cell.                        !
    !                                                                         !
    ! Uses a combination of the PGI accelerator pragmas and CUDA Fortran.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   in1 (input): Current col function.                                    !
    !   in2 (input): An optional second col function (usually present).       !
    !   potential_box_dbl_pair (input): Region of the potential extracted from!
    !           the simulation cell that is to be applied to the col          !
    !           function(s).                                                  !
    !   apply_pot_2 (input): Indicates the presence of the second col function!
    !-------------------------------------------------------------------------!
    ! Written by Karl Wilkinson 2012                                          !
    !                                                                         !
    ! In addition to obvious changes relating to the pragmas that control     !
    ! execution on the GPU the most significant differences relative to the   !
    ! CPU code are:                                                           !
    !      - Addition of CPU<>GPU data transfer steps.                        !
    !      - Direct calls to the FFT routines within this routine rather      !
    !          than using a wrapper.                                          !
    !      - Calculation of col function * potential within this routine      !
    !          rather than within the calling routine.                        !
    !=========================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: DP, I4B
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor !! External dependency
    use timer, only: timer_clock

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: potential_box_dbl_pair(:,:,:,:)
    integer :: is, batch_count

    ! Local variables
    integer :: n1,n2,n3,i2,i1
    integer :: n1fine, n2fine, n3fine
    real(kind=DP) :: scalefac
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3
    logical :: apply_pot_2

    real(kind=DP) :: test_r(3)
    integer :: npt
    integer :: istat

    ! Stream variables
    integer(kind=I4B) :: streamA, streamB

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scalefac = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)


!$acc data region
    ! kaw: Zero the fine work array. This is automatically done on the GPU as it has been declared as
    !      a CUDA Fortran array
    fine_work_gpu = 0.0_DP

    ! kaw: Transfer potential data to GPU
    streamA = 4
    streamB = 5
    istat = cudaMemcpyAsync(potential_box_dbl_gpu, potential_box_dbl_pair, 16*n1*n2*n3, stream2)

!$acc kernels async(streamA)
    coarse_work_gpu(1:n1,1:n2,1:n3) = scalefac * cmplx(fftbox_batch_gpu(:,:,:,is,3,batch_count), &
         fftbox_batch_gpu(:,:,:,is,4,batch_count),kind=DP)
!$acc end kernels

!$acc wait(streamA)
    call cufftSetStream(cufftplan_coarse, stream1)
    ! Fourier transform to reciprocal space on coarse grid
    call cufftExec(cufftplan_coarse,planType,coarse_work_gpu,coarse_work_gpu,CUFFT_FORWARD)
    call cufftSetStream(cufftplan_coarse, 0)

    ! Populate fine grid array
    istat = cudaStreamSynchronize(stream1)
!$acc kernels async(streamA)
             fine_work_gpu(1:m1,1:m2,1:m3)    = coarse_work_gpu(1:m1,1:m2,1:m3)
          fine_work_gpu(k1:2*n1,1:m2,1:m3)   = coarse_work_gpu(l1:n1,1:m2,1:m3)
          fine_work_gpu(1:m1,k2:2*n2,1:m3)   = coarse_work_gpu(1:m1,l2:n2,1:m3)
       fine_work_gpu(k1:2*n1,k2:2*n2,1:m3)  = coarse_work_gpu(l1:n1,l2:n2,1:m3)
          fine_work_gpu(1:m1,1:m2,k3:2*n3)   = coarse_work_gpu(1:m1,1:m2,l3:n3)
       fine_work_gpu(k1:2*n1,1:m2,k3:2*n3)  = coarse_work_gpu(l1:n1,1:m2,l3:n3)
       fine_work_gpu(1:m1,k2:2*n2,k3:2*n3)  = coarse_work_gpu(1:m1,l2:n2,l3:n3)
    fine_work_gpu(k1:2*n1,k2:2*n2,k3:2*n3) = coarse_work_gpu(l1:n1,l2:n2,l3:n3)
!$acc end kernels

    ! Execute inverse FFT.
    !$acc wait(streamA)
    call cufftSetStream(cufftplan_fine, stream1)
    call cufftExec(cufftplan_fine,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE)
    call cufftSetStream(cufftplan_fine, 0)

!   ! kaw: In theory the following should be faster than the 3D FFT above.
!   !      This is not the case in practice but will be revisited with input
!   !      from the library engineers at NVIDIA.
    !n1fine = n1*2
    !n2fine = n2*2
    !n3fine = n3*2
    !j1 = 3*n1fine/4+2
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
    !do i2=1,n2fine/4+1
    !    call cufftExec(cufft_n3_2,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,(n1fine*(i2-1))+j1-1)
    !end do
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
    !do i2=3*n2fine/4+2,n2fine
    !    call cufftExec(cufft_n3_2,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,(n1fine*(i2-1))+j1-1)
    !end do
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=1,n1fine/4+1
    !   call cufftExec(cufft_n2,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,i1-1)
    !end do
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=3*n1fine/4+2,n1fine
    !   call cufftExec(cufft_n2,planType,fine_work_gpu,fine_work_gpu,CUFFT_INVERSE,i1-1)
    !end do
    !! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
    !call cufftExec(cufft_n1, planType, fine_work_gpu, fine_work_gpu, CUFFT_INVERSE)


    ! ndmh: unpack data from real and imaginary parts of complex array
    ! ndmh: and multiply them together to find product
    ! -------- Apply potential for NGWFs as needed ------------
    istat = cudaStreamSynchronize(stream2)
    istat = cudaStreamSynchronize(stream1)
    if (apply_pot_2) then
       !$acc kernels async(streamA)
       do k3 = 1, 2*n3
          do k2 = 1, 2*n2
             fftbox_fine_1_gpu(:,k2,k3) = real(fine_work_gpu(:,k2,k3),DP)
             fftbox_fine_2_gpu(:,k2,k3) = aimag(fine_work_gpu(:,k2,k3))
             !fine_work_gpu(:,k2,k3) = cmplx(potential_box_dbl_gpu(:,k2,k3,1)*real(fine_work_gpu(:,k2,k3),
             !                               potential_box_dbl_gpu(:,k2,k3,2)*aimag(fine_work_gpu(:,k2,k3))))
          end do
       end do
       !$acc end kernels

       !$acc kernels async(streamA)
       do k3 = 1, 2*n3
          do k2 = 1, 2*n2
             fftbox_fine_1_gpu(:,k2,k3) = potential_box_dbl_gpu(:,k2,k3,1)*real(fine_work_gpu(:,k2,k3),DP)
             fftbox_fine_2_gpu(:,k2,k3) = potential_box_dbl_gpu(:,k2,k3,2)*aimag(fine_work_gpu(:,k2,k3))
             !fine_work_gpu(:,k2,k3) = cmplx(potential_box_dbl_gpu(:,k2,k3,1)*real(fine_work_gpu(:,k2,k3),
             !                               potential_box_dbl_gpu(:,k2,k3,2)*aimag(fine_work_gpu(:,k2,k3))))
          end do
       end do
       !$acc end kernels

    else
       !$acc kernels async(streamA)
       do k3 = 1, 2*n3
          do k2 = 1, 2*n2
             fftbox_fine_1_gpu(:,k2,k3) = potential_box_dbl_gpu(:,k2,k3,1)*real(fine_work_gpu(:,k2,k3),DP)
             !fine_work_gpu(:,k2,k3) = cmplx(potential_box_dbl_gpu(:,k2,k3,1)*real(fine_work_gpu(:,k2,k3),
             !                               aimag(fine_work_gpu(:,k2,k3))))
          end do
       end do
       !$acc end kernels
    end if
!$acc end data region

    !$acc wait(streamA)

  end subroutine fourier_gpu_interpolate_grad

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_gpu_filter_grad(is, batch_count, multibatch)
    !=========================================================================!
    ! Uses Nvidia GPUs to perform a Fast Fourier Transform from real to       !
    ! reciprocal space, downscales to a coarse grid then returns to real space!
    !                                                                         !
    ! Uses a combination of the PGI accelerator pragmas and CUDA Fortran.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   batch_count (input): Which col function of the batch is being         !
    !           processed.                                                    !
    !   multibatch (input): Indicates the presence of the second col function !
    !-------------------------------------------------------------------------!
    ! Written by Karl Wilkinson 2012                                          !
    !                                                                         !
    ! In addition to obvious changes relating to the pragmas that control     !
    ! execution on the GPU the most significant differences relative to the   !
    ! CPU code are:                                                           !
    !      - Direct calls to the FFT routines within this routine rather      !
    !          than using a wrapper.                                          !
    !      - Removal of the out1 and out2 output arrays. These are dealt with !
    !          at a later stage in the calculation to allow the local         !
    !          potential integrals to be calculated on the GPU itself         !
    !          (depending upon) the relative data transfer costs.             !
    ! This routine is called from potential_gpu_app2_ngwf_batch and the   !
    !
    !=========================================================================!

    use constants, only: DP
    use fourier_gpu_wrapper_mod
    use accel_lib !! External dependency
    use cudafor !! External dependency
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer,intent(in) :: batch_count
    logical,intent(in) :: multibatch
    integer :: is

    ! Local variables
    integer :: n1,n2,n3
    integer :: n1fine,n2fine,n3fine,i1,i2,j1
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3
    real(kind=DP) :: scale
    integer :: istat

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3


    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1 + l1 ; k2 = n2 + l2 ; k3 = n3 + l3


!$acc data region
!$acc kernels
    ! Pack data into real and imaginary parts of complex array
    fine_work_gpu(1:2*n1,1:2*n2,1:2*n3) = &
                                 cmplx(fftbox_fine_1_gpu(1:2*n1,1:2*n2,1:2*n3), &
                                 fftbox_fine_2_gpu(1:2*n1,1:2*n2,1:2*n3),DP)
!$acc end kernels

    ! Fourier transform to reciprocal space on fine grid
! fourier_filter: fine forward - dfftw_forwards batches
    call cufftExec(cufftplan_fine,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD)

!   ! kaw: In theory the following should be faster than the 3D FFT above.
!   !      This is not the case in practice but will be revisited with input
!   !      from the library engineers at NVIDIA.
    !n1fine = n1*2
    !n2fine = n2*2
    !n3fine = n3*2
    !! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
    !call cufftExec(cufft_n1, planType, fine_work_gpu, fine_work_gpu, CUFFT_FORWARD)
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=1,n1fine/4+1
    !   call cufftExec(cufft_n2,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,i1-1)
    !end do
    !! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
    !do i1=3*n1fine/4+2,n1fine
    !   call cufftExec(cufft_n2,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,i1-1)
    !end do
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
    !j1 = 3*n1fine/4+2
    !do i2=1,n2fine/4+1
    !    call cufftExec(cufft_n3_2,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,(n1fine*(i2-1))+j1-1)
    !end do
    !! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
    !! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
    !do i2=3*n2fine/4+2,n2fine
    !    call cufftExec(cufft_n3_2,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,(n1fine*(i2-1)))
    !    call cufftExec(cufft_n3_3,planType,fine_work_gpu,fine_work_gpu,CUFFT_FORWARD,(n1fine*(i2-1))+j1-1)
    !end do


    ! Copy Fourier components to coarse reciprocal space grid
!$acc kernels
    coarse_work_gpu(1:m1,1:m2,1:m3) =    fine_work_gpu(1:m1,1:m2,1:m3)
    coarse_work_gpu(l1:n1,1:m2,1:m3) =   fine_work_gpu(k1:2*n1,1:m2,1:m3)
    coarse_work_gpu(1:m1,l2:n2,1:m3) =   fine_work_gpu(1:m1,k2:2*n2,1:m3)
    coarse_work_gpu(l1:n1,l2:n2,1:m3) =  fine_work_gpu(k1:2*n1,k2:2*n2,1:m3)
    coarse_work_gpu(1:m1,1:m2,l3:n3) =   fine_work_gpu(1:m1,1:m2,k3:2*n3)
    coarse_work_gpu(l1:n1,1:m2,l3:n3) =  fine_work_gpu(k1:2*n1,1:m2,k3:2*n3)
    coarse_work_gpu(1:m1,l2:n2,l3:n3) =  fine_work_gpu(1:m1,k2:2*n2,k3:2*n3)
    coarse_work_gpu(l1:n1,l2:n2,l3:n3) = fine_work_gpu(k1:2*n1,k2:2*n2,k3:2*n3)
!$acc end kernels

    ! Fourier transform to real space on coarse grid
    call cufftExec(cufftplan_coarse,planType,coarse_work_gpu,coarse_work_gpu,CUFFT_INVERSE)

    ! Normalisation factor
    scale = 0.125_DP / (internal_box_coarse_info%ld1*&
         internal_box_coarse_info%ld2*internal_box_coarse_info%n3)


    ! Unpack data from real and imaginary parts of complex array
    ! ndmh: a loop saves on cache misses when repeating for aimag
    if (multibatch) then
       !$acc kernels
       do k3 = 1,n3
          do k2 = 1,n2
             do k1 = 1,n1
                fftbox_batch_gpu(k1,k2,k3,is,3,batch_count) = real(scale*coarse_work_gpu(k1,k2,k3),DP)
                fftbox_batch_gpu(k1,k2,k3,is,4,batch_count) = aimag(scale*coarse_work_gpu(k1,k2,k3))
             end do
          end do
       end do
       !$acc end kernels
    else
       !$acc kernels
       do k3 = 1,n3
          do k2 = 1,n2
             do k1 = 1,n1
                fftbox_batch_gpu(k1,k2,k3,is,3,batch_count) = real(scale*coarse_work_gpu(k1,k2,k3),DP)
             end do
          end do
       end do
       !$acc end kernels
    end if
!$acc end data region
  end subroutine fourier_gpu_filter_grad

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_gpu_init()
  ! kaw: Subroutine to initialise the GPU related resources

    use constants, only: DP, stdout
    use comms, only: pub_my_proc_id, pub_on_root, pub_total_num_procs, &
         comms_barrier
    use fourier_gpu_wrapper_mod
    use rundat, only: pub_fftbox_batch_size, pub_num_spins
    use utils, only: utils_assert
    use accel_lib !! External dependency
    use cudafor !! External dependency

    implicit none

    integer :: n1, n2, n3, n1dbl, n2dbl, n3dbl
    integer :: istat, ierr
    integer :: total_gpu_mem, gpu_total_memory, gpu_after_alloc_free_memory

    if (pub_on_root) write(stdout,*) '                 '
    if (pub_on_root) write(stdout,*) &
    '============================== Initialising GPUs ==============================='

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3
    n1dbl = n1*2
    n2dbl = n2*2
    n3dbl = n3*2

    ! Set up the devices
    num_devices = acc_get_num_devices(acc_device_nvidia)
    call utils_assert(.not.num_devices.eq.0, &
         'Attempting to run GPU code but no GPUs detected.')

    my_device = mod(pub_my_proc_id,num_devices)

! kaw: Changed this to the equivalent cuda fortran function as the pgi function
!      does not seem to work on the iridis GPU procs.
    call acc_set_device_num(my_device,acc_device_nvidia)
!    istat = cudaSetDevice(my_device)

    call comms_barrier
    write(stdout,*) '     MPI thread',pub_my_proc_id,'allocated to GPU',my_device, &
               'on proc:'
    call comms_barrier


!Allocate GPU arrays
    if (pub_on_root) write(stdout,*) '     Allocating DP Arrays on GPU'

    total_gpu_mem = 0
    allocate(coarse_work_gpu(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','coarse_work_gpu',ierr)
    total_gpu_mem = total_gpu_mem + ((2*n1*n2*n3*8)/(1024*1024))
    allocate(fine_work_gpu(n1dbl,n2dbl,n3dbl),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','fine_work_gpu',ierr)
    total_gpu_mem = total_gpu_mem + ((2*n1dbl*n2dbl*n3dbl*8)/(1024*1024))
    allocate(fftbox_coarse_1_gpu(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','fftbox_coarse_1_gpu',ierr)
    total_gpu_mem = total_gpu_mem + ((n1*n2*n3*8)/(1024*1024))
    allocate(fftbox_coarse_2_gpu(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','fftbox_coarse_2_gpu',ierr)
    total_gpu_mem = total_gpu_mem + ((n1*n2*n3*8)/(1024*1024))
!KAW, removed to reduce memory usage    allocate(fftbox_fine_1_gpu(n1dbl,n2dbl,n3dbl))
!KAW, removed to reduce memory usage    allocate(fftbox_fine_2_gpu(n1dbl,n2dbl,n3dbl))
    allocate(potential_box_dbl_gpu(n1dbl,n2dbl,n3dbl,2),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','potential_box_dbl_gpu',ierr)
    total_gpu_mem = total_gpu_mem + ((n1dbl*n2dbl*n3dbl*8)/(1024*1024))
    allocate(potential_fftbox_batch_gpu(n1,n2,n3,pub_fftbox_batch_size),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','potential_fftbox_batch_gpu',ierr)
    total_gpu_mem = total_gpu_mem + ((n1*n2*n3*8*pub_fftbox_batch_size)/(1024*1024))
    allocate(out_gpu(n1dbl,n2dbl,n3dbl,pub_num_spins),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','out_gpu',ierr)
    total_gpu_mem = total_gpu_mem + ((n1dbl*n2dbl*n3dbl*8)/(1024*1024))
!KAW, removed to reduce memory usage    allocate(pub_fftbox_recip_grid_gpu(n1,n2,n3))


!create streams
    if (pub_on_root)write(stdout,*) '     Creating GPU Streams'
    istat=cudaStreamCreate(stream1)
    istat=cudaStreamCreate(stream2)


!create GPU plans
    if (pub_on_root)write(stdout,*) '     Creating Plans for DP GPU FFTs'
    call cufftPlan3D(cufftplan_coarse,n3,n2,n1,planType)
    call cufftPlan3D(cufftplan_fine,n3dbl,n2dbl,n1dbl,planType)

!Get size of FFT plans
    call cufftGetSize(cufftplan_coarse, coarsesize)
    call cufftGetSize(cufftplan_fine, finesize)
    if (pub_on_root)write(stdout,*) '        Coarse Plan size (Mb):', &
                                     coarsesize/(1024*1024)
    total_gpu_mem = total_gpu_mem + coarsesize/(1024*1024)
    if (pub_on_root)write(stdout,*) '        Fine Plan size (Mb):', &
                                     finesize/(1024*1024)
    total_gpu_mem = total_gpu_mem + finesize/(1024*1024)

!Feedback about GPU memopry usage
    if (pub_on_root)write(stdout,*) ' Estimated GPU RAM usage per MPI rank:', total_gpu_mem ,'(Mb):'


    !call cufftPlanMany(cufft_n1,1,n1dbl,n1dbl,1,n1dbl,n1dbl,1,n1dbl,planType,n2dbl*n3dbl)
    !call cufftPlanMany(cufft_n2,1,n2dbl,n2dbl,n1dbl,n1dbl*n2dbl,n2dbl,n1dbl,n1dbl*n2dbl,planType,n3dbl)
    !call cufftPlanMany(cufft_n3_2,1,n3dbl,n3dbl,n1dbl*n2dbl,1,n3dbl,n1dbl*n2dbl,1,planType,n1dbl/4+1)
    !call cufftPlanMany(cufft_n3_3,1,n3dbl,n3dbl,n1dbl*n2dbl,1,n3dbl,n1dbl*n2dbl,1,planType,n1dbl/4)


#ifdef GPU_SP_TEST
    GPU_SP_switched = .false.
    if (pub_on_root) write(stdout,*) '     Allocating SP Arrays on GPU'
    allocate(coarse_work_gpuSP(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','coarse_work_gpuSP',ierr)
    allocate(fine_work_gpuSP(n1dbl,n2dbl,n3dbl),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','fine_work_gpuSP',ierr)
    allocate(fftbox_coarse_1_gpuSP(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','fftbox_coarse_1_gpuSP',ierr)
    allocate(fftbox_coarse_2_gpuSP(n1,n2,n3),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','fftbox_coarse_2_gpuSP',ierr)
    allocate(fftbox_fine_1_gpuSP(n1dbl,n2dbl,n3dbl),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','fftbox_fine_1_gpuSP',ierr)
    allocate(fftbox_fine_2_gpuSP(n1dbl,n2dbl,n3dbl),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','fftbox_fine_1_gpuSP',ierr)
    allocate(potential_box_dbl_gpuSP(n1dbl,n2dbl,n3dbl,2),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','potential_box_dbl_gpuSP',ierr)
    allocate(potential_fftbox_batch_gpuSP(n1,n2,n3,pub_fftbox_batch_size),stat=ierr)
    call utils_alloc_check('fourier_gpu_init','potential_fftbox_batch_gpuSP',ierr)

    if (pub_on_root)write(stdout,*) '     Creating Plans for SP GPU FFTs'
    call cufftPlan3D(cufftplan_coarse_SP,n3,n2,n1,planType_SP)
    call cufftPlan3D(cufftplan_fine_SP,n3*2,n2*2,n1*2,planType_SP)
    !call cufftPlanMany(cufft_n1_SP,1,n1dbl,n1dbl,1,n1dbl,n1dbl,1,n1dbl,planType_SP,n2dbl*n3dbl)
    !call cufftPlanMany(cufft_n2_SP,1,n2dbl,n2dbl,n1dbl,n1dbl*n2dbl,n2dbl,n1dbl,n1dbl*n2dbl,planType_SP,n3dbl)
    !call cufftPlanMany(cufft_n3_2_SP,1,n3dbl,n3dbl,n1dbl*n2dbl,1,n3dbl,n1dbl*n2dbl,1,planType_SP,n1dbl/4+1)
    !call cufftPlanMany(cufft_n3_3_SP,1,n3dbl,n3dbl,n1dbl*n2dbl,1,n3dbl,n1dbl*n2dbl,1,planType_SP,n1dbl/4)
#endif

    if (pub_on_root) write(stdout,*) &
    '=============================== GPUs Initialised ==============================='
    if (pub_on_root) write(stdout,*) ' '
  end subroutine fourier_gpu_init
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine fourier_gpu_exit
  ! kaw: Subroutine to clean up the GPU related resources
    use constants, only: stdout
    use comms, only : pub_on_root
    use accel_lib !! External dependency
    use fourier_gpu_wrapper_mod
    integer :: istat, ierr
    if (pub_on_root) write(stdout,*) ' '
    if (pub_on_root) write(stdout,*) &
    '================================= Clearing GPUs ================================'

#ifdef GPU_SP_TEST
  if (pub_on_root) write(stdout,*) 'SP arrays'
  deallocate(coarse_work_gpuSP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','coarse_work_gpuSP',ierr)
  deallocate(fine_work_gpuSP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','fine_work_gpuSP',ierr)
  deallocate(fftbox_coarse_1_gpuSP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','fftbox_coarse_1_gpuSP',ierr)
  deallocate(fftbox_coarse_2_gpuSP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','fftbox_coarse_2_gpuSP',ierr)
  deallocate(fftbox_fine_1_gpuSP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','fftbox_fine_1_gpuSP',ierr)
  deallocate(fftbox_fine_2_gpuSP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','fftbox_fine_2_gpuSP',ierr)
  deallocate(potential_box_dbl_gpuSP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','potential_box_dbl_gpuSP',ierr)
  deallocate(potential_fftbox_batch_gpuSP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','potential_fftbox_batch_gpuSP',ierr)
  deallocate(out_gpuSP,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','out_gpuSP',ierr)

  if (pub_on_root) write(stdout,*) 'Destroying SP FFT plans'
  call cufftDestroy(cufftplan_coarse_SP)
  call cufftDestroy(cufftplan_fine_SP)
  !call cufftDestroy(cufft_n3_3_SP)
  !call cufftDestroy(cufft_n3_2_SP)
  !call cufftDestroy(cufft_n2_SP)
  !call cufftDestroy(cufft_n1_SP)
#endif

  if (pub_on_root) write(stdout,*) 'Deallocating DP arrays'
!KAW, removed to reduce mempry usage!  deallocate(pub_fftbox_recip_grid_gpu)
  deallocate(potential_fftbox_batch_gpu,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','potential_fftbox_batch_gpu',ierr)
  deallocate(potential_box_dbl_gpu,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','potential_box_dbl_gpu',ierr)
!KAW, removed to reduce mempry usage  deallocate(fftbox_fine_2_gpu)
!KAW, removed to reduce mempry usage  deallocate(fftbox_fine_1_gpu)
  deallocate(out_gpu,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','out_gpu',ierr)
  deallocate(fftbox_coarse_2_gpu,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','fftbox_coarse_2_gpu',ierr)
  deallocate(fftbox_coarse_1_gpu,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','fftbox_coarse_1_gpu',ierr)
  deallocate(fine_work_gpu,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','fine_work_gpu',ierr)
  deallocate(coarse_work_gpu,stat=ierr)
  call utils_dealloc_check('fourier_gpu_exit','coarse_work_gpu',ierr)

  if (pub_on_root) write(stdout,*) 'Destroying DP FFT plans'
  call cufftDestroy(cufftplan_coarse)
  call cufftDestroy(cufftplan_fine)
  !call cufftDestroy(cufft_n3_3)
  !call cufftDestroy(cufft_n3_2)
  !call cufftDestroy(cufft_n2)
  !call cufftDestroy(cufft_n1)

  if (pub_on_root) write(stdout,*) 'Destroying CUDA Streams'
  istat=cudaStreamDestroy(stream1)
  istat=cudaStreamDestroy(stream2)


    if (pub_on_root) write(stdout,*) &
    '================================= GPUs Cleared ================================='
    if (pub_on_root) write(stdout,*) ' '
  end subroutine fourier_gpu_exit
#else
  subroutine dummy_gpu
  use fourier_gpu_wrapper_mod
  !kaw: This is a dummy use statement to make the check_dependencies script happy.
  end subroutine dummy_gpu
#endif
end module fourier

