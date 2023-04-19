!> switch between 'use mpi' and 'include mpif.h'
!!
!! Lucian Anton
!! February 2012

module dl_mg_mpi_header
#ifdef MPI
#ifdef USE_INCLUDE_MPIF
  !implicit none
  include 'mpif.h'
#else
  use mpi
#endif
#endif

end module dl_mg_mpi_header
