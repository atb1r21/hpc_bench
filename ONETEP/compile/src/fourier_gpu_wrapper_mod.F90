module fourier_gpu_wrapper_mod
! kaw: This module is currently slightly haphazard. The contents should mostly
!      be distributed into the fourier and wrapper modules. There are also
!      data structures that are essentially reused and a clean up is needed to
!      reduce this issue.

  use constants, only: DP, SP
  use utils
#ifdef GPU_PGI
  use accel_lib !! External dependency
  use constants, only: I4B
  use cudafor !! External dependency
#endif

#ifdef GPU_PGI
  ! kaw: CUFFT variables, both DP and SP.
  integer, parameter, public :: CUFFT_FORWARD = -1
  integer, parameter, public :: CUFFT_INVERSE = 1

  ! kaw: Used for full 3D FFT box FFTs
  integer, parameter, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex
  integer :: planType = CUFFT_Z2Z
  integer :: cufftplan_coarse, cufftplan_fine
  ! kaw: Used for batched 1D FFT box FFTs
  integer :: cufft_n1, cufft_n2, cufft_n3_2, cufft_n3_3 ! Functional but slower than 3D FFTs
  integer :: coarsesize, finesize !Memory used by plans on GPU

#ifdef GPU_SP_TEST
  integer, parameter, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
  ! kaw: Used for full 3D FFT box FFTs
  integer :: planType_SP = CUFFT_C2C
  integer :: cufftplan_coarse_SP, cufftplan_fine_SP
  ! kaw: Used for batched 1D FFT box FFTs
  integer :: cufft_n1_SP, cufft_n2_SP, cufft_n3_2_SP, cufft_n3_3_SP
#endif

  ! kaw: Dynamic precision based variables.
  logical :: GPU_SP     !kaw: Why is this outside the ifdef?
#ifdef GPU_SP_TEST
  logical :: GPU_SP_toggle
  logical :: GPU_SP_switched
  real(kind=DP) :: GPU_SP_switch = 5.0E-2_DP
  real(kind=DP) :: deltaE_gpu = 0.0_DP
  real(kind=DP) :: E_gpu_previous = 0.0_DP
  !kaw: General arrays, the number of these should be reduced to make more
  !     efficient use of memory.
  real(kind=SP), allocatable, device :: out_gpuSP(:,:,:,:)
  real(kind=SP), allocatable, device :: fftbox_coarse_1_gpuSP(:,:,:)
  real(kind=SP), allocatable, device :: fftbox_coarse_2_gpuSP(:,:,:)
  real(kind=SP), allocatable, device :: potential_box_dbl_gpuSP(:,:,:,:)
  real(kind=SP), allocatable, device :: potential_fftbox_batch_gpuSP(:,:,:,:)
  real(kind=SP), allocatable, device :: fftbox_fine_1_gpuSP(:,:,:)
  real(kind=SP), allocatable, device :: fftbox_fine_2_gpuSP(:,:,:)
  !kaw: Work arrays, consider adding another dimension to these to allow
  !     the unrolling of loops.
  complex(kind=SP), allocatable, device :: coarse_work_gpuSP(:,:,:)
  complex(kind=SP), allocatable, device :: fine_work_gpuSP(:,:,:)
#endif

  ! kaw: For the control of procs containing multiple GPUs
  integer :: num_devices, my_device
  !kaw: DP work arrays, consider adding another dimension to these to allow
  !     the unrolling of loops.
  complex(kind=DP), allocatable, device :: coarse_work_gpu(:,:,:)
  complex(kind=DP), allocatable, device :: fine_work_gpu(:,:,:)
  !kaw: General arrays, the number of these should be reduced to make more
  !     efficient use of memory.
  real(kind=DP), allocatable, device :: fftbox_coarse_1_gpu(:,:,:)
  real(kind=DP), allocatable, device :: fftbox_coarse_2_gpu(:,:,:)
  real(kind=DP), allocatable, device :: potential_fftbox_batch_gpu(:,:,:,:)
  real(kind=DP), allocatable, device :: potential_box_dbl_gpu(:,:,:,:)
  real(kind=DP), allocatable, device :: fftbox_fine_1_gpu(:,:,:)
  real(kind=DP), allocatable, device :: fftbox_fine_2_gpu(:,:,:)
  real(kind=DP), allocatable, device :: out_gpu(:,:,:,:)
  real(kind=DP), allocatable, dimension(:,:,:,:,:,:), device :: fftbox_batch_gpu
  real(kind=DP), allocatable, dimension(:,:,:,:,:,:), device :: nonloc_batch_gpu
  real(kind=DP), allocatable, dimension(:,:,:), device :: pub_fftbox_recip_grid_gpu
  integer(kind=I4B) :: stream1, stream2, stream3

  ! kaw: Check if this is used.
  interface cufftPlan1d
     subroutine cufftPlan1d(plan, nx, type, batch) bind(C,name='cufftPlan1d')
       use iso_c_binding !! External dependency
       integer(c_int):: plan
       integer(c_int),value:: nx, batch,type
     end subroutine cufftPlan1d
  end interface

  ! kaw: Check if this is used.
  interface cufftPlan2d
     subroutine cufftPlan2d(plan, nx, ny, type) bind(C,name='cufftPlan2d')
       use iso_c_binding !! External dependency
       integer (c_int) :: plan
       integer (c_int), value :: nx, ny, type
     end subroutine cufftPlan2d
  end interface

  interface cufftPlan3d
     subroutine cufftPlan3d(plan, nx, ny, nz, type) bind(C,name='cufftPlan3d')
       use iso_c_binding !! External dependency
       integer (c_int) :: plan
       integer (c_int), value :: nx, ny, nz, type
     end subroutine cufftPlan3d
  end interface

  interface cufftPlanMany
     subroutine cufftPlanMany(plan, rank, n, &
                              iembed, istride, idist, &
                              oembed, ostride, odist, &
                              type, nbatch) bind(C,name='cufftPlanMany')
       use iso_c_binding !! External dependency
       !pgi$ ignore_tkr n, iembed, oembed
       integer (c_int) :: plan            !Plan
       integer (c_int),value :: rank      !Dimensions of data to transform
       integer (c_int) :: n(*)            !array of size rank, contains size of each dimension
       integer (c_int) :: iembed(*)       !Pointer of size rank indicates input storage dimensions
       integer (c_int),value :: istride   !Distance between two successive input elements in
       integer (c_int),value :: idist     !Indicates the distance between the first element of two
       integer (c_int) :: oembed(*)       !See above, but for output.
       integer (c_int),value :: ostride
       integer (c_int),value :: odist
       integer (c_int),value :: type      !Transform data types
       integer (c_int),value :: nbatch    !Number of batches
     end subroutine cufftPlanMany
  end interface

  ! ------------
  ! cufftPlanSize
  ! ------------

  interface cufftGetSize
     subroutine cufftGetSize(plan,plansize) bind(C,name='cufftGetSize')
       use iso_c_binding !! External dependency
       integer(c_int),value:: plan
       integer(c_int):: plansize
     end subroutine cufftGetSize
  end interface

  ! ------------
  ! cufftDestroy
  ! ------------

  interface cufftDestroy
     subroutine cufftDestroy(plan) bind(C,name='cufftDestroy')
       use iso_c_binding !! External dependency
       integer(c_int),value:: plan
     end subroutine cufftDestroy
  end interface

  ! ------------
  ! cufftSetStream
  ! ------------
  interface cufftSetStream
     subroutine cufftSetStream(plan,stream) bind(C,name='cufftSetStream')
       use iso_c_binding !! External dependency
       integer(c_int),value:: plan,stream
     end subroutine cufftSetStream
  end interface

  ! A generic interface such as "cufftExec" which contained all the R2C, C2C, ...
  ! subroutine variants would work fine as long as transforms were not done in place.
  ! The types and kinds associated with idata and odata would determine the
  ! correct routine.
  ! However, this appoach would break down when performing in-place transforms.
  ! For example, suppose one wanted to perform a real to complex transform in-place.
  ! The soubroutine call:
  !    call cufftExec(plan, data, data)
  ! would result in a compile time error if "data" were of type real or
  ! the cufftExecC2C would be executed if "data" were of type complex
  ! So, we define a separate interface for each R2C, C2C, ... variant and
  ! use the "!dir$ ignore_tkr" directive to ignore compile-time checks of
  ! data type, kind, and ranks.

  ! ----------------------------------
  ! cufftExec[C2C|R2C|C2R|Z2Z|Z2D|D2Z]
  ! ----------------------------------
  interface cufftExecC2C
     subroutine cufftExecC2C(plan, idata, odata, direction) &
          & bind(C,name='cufftExecC2C')
       use iso_c_binding !! External dependency
       integer(c_int),value:: plan
       !pgi$ ignore_tkr idata, odata
       complex(kind=kind(1.0d0)), device :: idata(*), odata(*)
       integer(c_int),value:: direction
     end subroutine cufftExecC2C
  end interface

  interface cufftExecZ2Z
     subroutine cufftExecZ2Z(plan, idata, odata, direction) &
          & bind(C,name='cufftExecZ2Z')
       use iso_c_binding !! External dependency
       integer(c_int),value:: plan
       !pgi$ ignore_tkr idata, odata
       complex(kind=kind(1.0d0)), device :: idata(*), odata(*)
       integer(c_int),value:: direction
     end subroutine cufftExecZ2Z
  end interface

  interface cufftExecD2Z
     subroutine cufftExecD2Z(plan, idata, odata) &
          & bind(C,name='cufftExecD2Z')
       use iso_c_binding !! External dependency
       integer(c_int),value:: plan
       !pgi$ ignore_tkr idata, odata
       complex(kind=kind(1.0d0)), device :: idata(*), odata(*)
     end subroutine cufftExecD2Z
  end interface

  interface cufftExec
     subroutine cufftExec(plan, transform, idata, odata, direction, offset)
       integer :: plan, transform
       !pgi$ ignore_tkr idata, odata
       real, device :: idata(*), odata(*)
       integer, optional :: direction, offset
     end subroutine cufftExec
  end interface
 !kaw: old version
     !subroutine cufftExec(plan, transform, idata, odata, direction)
     !  integer :: plan, transform
     !  !pgi$ ignore_tkr idata, odata
     !  real, device :: idata(*), odata(*)
     !  integer, optional :: direction
     !end subroutine cufftExec

end module fourier_gpu_wrapper_mod

! A general subroutine that uses the transform type passed in "transform" to
! determine which routine to call.  The interfaces for the routines specified
! above use the "!dir$ ignore_tkr" directive to remove checking of
! type, kind, and rank to facilitate certain in-place transforms (eg. cufftExecR2C)
! and multidimensional arrays

 subroutine cufftExec(plan, transform, idata, odata, direction, offset)
  use fourier_gpu_wrapper_mod !! External dependency
  use utils, only: utils_abort
  implicit none
  integer :: plan, transform
  real, device :: idata(0:*), odata(0:*)
  integer :: direction
  integer, optional :: offset
  integer :: off

  off = 0
  if( present(offset) ) off = offset
  ! The values of off are scaled by 2 or 4 below since the double and
  ! double complex data types are 2 or 4 times bigger than the real
  ! type used to declare idata/odata.

  select case (transform)
  case (CUFFT_Z2Z)
     call cufftExecZ2Z(plan, idata(4*off), odata(4*off), direction)
  case default
     call utils_abort('Invalid transform type passed to cufftExec.')
  end select
 end subroutine cufftExec

#else

#ifdef BASIC_GPU_FFT
  !kaw: This is used to compile a version that simply offloads the FFT operations
  !     but none of the surrounding code. Should be discarded as it had no
  !     function other than to provide a comparison for the paper.
  integer, parameter, public :: CUFFT_FORWARD = -1
  integer, parameter, public :: CUFFT_INVERSE = 1
  integer, parameter, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
  integer, parameter, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex
  integer, parameter, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex
  integer, parameter, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double
  integer :: planType = CUFFT_Z2Z
  integer :: cufftplan_coarse, cufftplan_fine
  complex(kind=DP), allocatable, device :: coarse_work_gpu(:,:,:)
  complex(kind=DP), allocatable, device :: fine_work_gpu(:,:,:)
  interface cufftPlan3d
     subroutine cufftPlan3d(plan, nx, ny, nz, type) bind(C,name='cufftPlan3d')
       use iso_c_binding !! External dependency
       integer (c_int) :: plan
       integer (c_int), value :: nx, ny, nz, type
     end subroutine cufftPlan3d
  end interface
  interface cufftDestroy
     subroutine cufftDestroy(plan) bind(C,name='cufftDestroy')
       use iso_c_binding !! External dependency
       integer(c_int),value:: plan
     end subroutine cufftDestroy
  end interface
  interface cufftExecZ2Z
     subroutine cufftExecZ2Z(plan, idata, odata, direction) &
          & bind(C,name='cufftExecZ2Z')
       use iso_c_binding !! External dependency
       integer(c_int),value:: plan
       !pgi$ ignore_tkr idata, odata
       complex(kind=kind(1.0d0)), device :: idata(*), odata(*)
       integer(c_int),value:: direction
     end subroutine cufftExecZ2Z
  end interface
  interface cufftExec
     subroutine cufftExec(plan, transform, idata, odata, direction, offset)
       integer :: plan, transform
       !pgi$ ignore_tkr idata, odata
       real, device :: idata(*), odata(*)
       integer, optional :: direction, offset
     end subroutine cufftExec
  end interface
#endif


end module fourier_gpu_wrapper_mod


#ifdef BASIC_GPU_FFT
 subroutine cufftExec(plan, transform, idata, odata, direction, offset)
  use fourier_gpu_wrapper_mod !! External dependency
  use utils, only: utils_abort
  implicit none
  integer :: plan, transform
  real, device :: idata(0:*), odata(0:*)
  integer :: direction,
  integer, optional :: offset
  integer :: off

  off = 0
  if( present(offset) ) off = offset
  ! The values of off are scaled by 2 or 4 below since the double and
  ! double complex data types are 2 or 4 times bigger than the real
  ! type used to declare idata/odata.

  select case (transform)
  case (CUFFT_Z2Z)
     call cufftExecZ2Z(plan, idata(4*off), odata(4*off), direction)
  case default
     call utils_abort('Invalid transform type passed to cufftExec.')
  end select
 end subroutine cufftExec
#endif

#endif






