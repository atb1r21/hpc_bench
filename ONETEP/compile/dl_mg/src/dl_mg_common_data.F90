!> \brief Store of common data used by other modules
!!
!! \todo In order to protect data the variable should be declare private
!! and use get/set functions to access and set values
!
!! Lucian Anton
!! February 2015
!!

module dl_mg_common_data
  use dl_mg_params, only : wp
  use dl_mg_types
  implicit none


  integer, save  :: nx, ny, nz      !< fine grid size (inner points)
  integer, save  :: nxb, nyb, nzb   !< fine grid size (inner + boundary points)
  integer, save  :: nxc, nyc, nzc   !< coarse grid size
  integer, save  :: npx, npy, npz   !< number of MPI ranks
  integer, save  :: isx, isy, isz, iex, iey, iez !< global indices range for the received data
  integer, save  :: mx, my, mz                   !< per-MPI rank fine grid size

  integer, save  :: bc(3)              !< boundary condition type
  logical        :: fullpbc = .false.  !< short hand for full PBC

  real(wp), save :: dx, dy, dz   !< fine grid cell units
  real(wp), save :: grid_weight  !< = dx * dy *dz

  integer, save  :: mg_levels     !< number of multigrid levels

! weight for restriction
  integer, save  :: restriction_weight !< set in dl_mg_init
  integer, save  :: mg_comm !< communicator to be used by MG for general messages.
                     !! For data transfer during computation one must use
                     !! the communicator provided in mg_t derived type

  integer, save  :: report_unit !< unit for log messages
  character(128), save :: report_file !< output for solver info

!> aggregation parameters
!> full_aggregation_size takes precedence over full_aggregation_level if both presen in dl_mg_init
  integer :: full_agg_level = 0 !< level at which full aggregation is done
  integer :: full_agg_size  = 0 !< minum size for distributed MG

!> number of OpenMP threads
  integer :: nthreads = 1

!> array of block starts and block sizes
  type(block_list_t), allocatable, target, save :: blk(:)


end module dl_mg_common_data
