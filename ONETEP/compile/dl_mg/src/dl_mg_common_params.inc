! public parameters which are needed across solver
! putting them in a  module won't work because I want the user application
! to access dl_mg only with one 'use' statement, i.e.
!
! use  dl_mg
!
! If dl_mg would use a separate module for these parameters,
! some compilers will require that module in the include directory

! Lucian Anton August 2013
!
! update December 2014

! boundary conditions
integer, parameter :: DL_MG_BC_PERIODIC = 0, DL_MG_BC_NEWMANN = 1, DL_MG_BC_DIRICHLET = 2

! error codes
integer, parameter :: DL_MG_SUCCESS = 0



  integer, parameter :: DL_MG_ERR_UNSPECIFIED    = 1, & !< use for cases which
       !! have no errors codes defined
       !! or when aborting on error.
       !! Extended error message can be added via msg
       !! argument of handle_error, which is printed
       !! only if execution is aborted
       DL_MG_ERR_INIT            = 2, &
       DL_MG_ERR_NEUMANN_BC      = 3, &
       DL_MG_ERR_UNKNOWN_BC      = 4, &
       DL_MG_ERR_MPI_TOPO        = 5, &
       DL_MG_ERR_MPI_COMM_DUP    = 6, &
       DL_MG_ERR_MPI_TOPO_BC     = 7, &
       DL_MG_ERR_GRID_EMPTY      = 8, &
       DL_MG_ERR_PROLONG_MAP     = 9, &
       DL_MG_ERR_DEBYE           = 10, &
       DL_MG_ERR_NOINIT_POISSON  = 11, &
       DL_MG_ERR_NOINIT_PBE      = 12, &
       DL_MG_ERR_DERPOT          = 13, &
       DL_MG_ERR_MOD_BC          = 14, &
       DL_MG_ERR_RES_RATIO       = 15, &
       DL_MG_ERR_NITER           = 16, &
       DL_MG_ERR_NEWTON_TYPE     = 17, &
       DL_MG_ERR_NEWTON_DAMP     = 18, &
       DL_MG_ERR_NEWTON_DAMP_DERFUN = 19, &
       DL_MG_ERR_ASSERTION       = 20, & ! JCW
       DL_MG_ERR_ALLOCATION      = 21, & ! JCW
       DL_MG_ERR_DEALLOCATION    = 22, & ! JCW
       DL_MG_ERR_IO              = 23, & ! JCW
       DL_MG_ERR_GRID_GEOM       = 24, & ! JCW
       DL_MG_ERR_MPI_FAIL        = 25, & ! JCW
       DL_MG_ERR_EQTYPE          = 26, &
       DL_MG_ERR_FD_ORDER        = 27, & ! JCW
       DL_MG_ERR_DEFCO_ITER      = 28, & ! JCW
       DL_MG_ERR_DEFCO_UNPHYSICAL= 29, & ! JCW
       DL_MG_ERR_DEFCO_DAMPING   = 30, &
       DL_MG_ERR_DEFCO_DERPOT    = 31, &
       DL_MG_ERR_NOMPI           = 32, &
       DL_MG_ERR_CION_ZERO       = 33, &
       DL_MG_ERR_CION_NEG        = 34

! error codes that could allow the code to use the (partial) solution
! if dl_mg returns before convergence is reached
! made public in dl_mg
!integer, parameter :: DL_MG_ERR_NITER              = 16, &
!                      DL_MG_ERR_NEWTON_DAMP        = 18, &
!                      DL_MG_ERR_NEWTON_DAMP_DERFUN = 19, &
!                      DL_MG_ERR_DEFCO_ITER         = 28, &
!                      DL_MG_ERR_DEFCO_UNPHYSICAL   = 29, &
!                      DL_MG_ERR_DEFCO_DAMPING      = 30, &
!                      DL_MG_ERR_DEFCO_DERPOT       = 31

! lengths of error/info strings
integer, parameter :: DL_MG_MAX_ERROR_STRING = 511, DL_MG_VERSION_STRING_LENGTH=30
