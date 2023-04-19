! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!
!========================================================================!
!                                                                        !
!              The ONETEP code is written and maintained by              !
!        Chris-Kriton Skylaris, Arash A. Mostofi, Nicholas D.M. Hine,    !
!                  Jacek Dziedzic and Peter D. Haynes.                   !
!                                                                        !
!========================================================================!
!
! Module to hold keyword list. This must be updated as
! new keywords are brought into existence.
!
! The 'LABEL' is the LABEL as used in calling the esdf routines
! 'TYP' defines the TYPe, with the following syntax. It is 3 characters
! long.
! The first indicates:
!  I - integer
!  S - single
!  D - double
!  P - physical
!  T - string (text)
!  E - defined (exists)
!  L - boolean (logical)
!  B - block
! The second is always a colon (:)
! The third indicates the "level" of the keyword
!  B - Basic
!  I - Intermediate
!  E - Expert
!  D - Dummy
!
! 'DSCRPT' is a description of the variable. It should contain a (short) title
! enclosed between *! ... !*, and then a more detailed description of the
! variable.
!
module esdf_key

  implicit none

  type KW_TYPE
     Character(32)   :: LABEL
     Character(3)    :: TYP
     Character(3000) :: DSCRPT
     Character(30)   :: KWGRP=""
  end type KW_TYPE


  integer, parameter :: maxnumkw=2000
  integer, save      :: numkw
  type(KW_TYPE), save :: KW(maxnumkw)

contains

 subroutine esdf_key_init

  integer :: i
  i=1

  ! Now define the keywords
  KW(i)%LABEL = "LATTICE_CART"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! The simulation cell lattice vectors !*"
  KW(i)%KWGRP = "CELLDATA"
  i=i+1

  KW(i)%LABEL = "PSINC_SPACING"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! PSINC grid spacing in atomic units !*"
  KW(i)%KWGRP = "CELLDATA"
  i=i+1

  KW(i)%LABEL = "GEOMETRY"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! input geometry !*"
  KW(i)%KWGRP = "CELLDATA"
  i=i+1

  KW(i)%LABEL = "PPD_NPOINTS"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! PPD edge length in grid points for &
       &each lattice direction !*"
  i=i+1

  KW(i)%LABEL = "XC_FUNCTIONAL"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Exchange-correlation functional !*"
  KW(i)%KWGRP = "XC"
  i=i+1

  KW(i)%LABEL = "KERNEL_CUTOFF"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Density kernel radius!*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "MAXIT_NGWF_CG"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of NGWF conjugate gradients (CG) &
       &iterations !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MAXIT_LNV"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of LNV iterations !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "LNV_THRESHOLD_ORIG"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! LNV convergence threshold !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "NGWF_THRESHOLD_ORIG"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! NGWF convergence threshold !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "OVLP_FOR_NONLOCAL"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Overlap sparsity for nonlocal !*"
  i=i+1

  KW(i)%LABEL = "EXACT_LNV"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Use original LNV algorithm !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "PRECOND_RECIP"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Recip-space KE preconditioning !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "PRECOND_SCHEME"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Recip-space preconditioning scheme &
       & BG = Bowler-Gillan method;&
       & MAURI = Mauri method;&
       & TETER = Teter-Allen-Payne method !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "K_ZERO"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! KE preconditioning parameter !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "PRECOND_REAL"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Real-space KE preconditioning !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "OUTPUT_DETAIL"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Level of output detail &
       &BRIEF, NORMAL or VERBOSE !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "KERNEL_UPDATE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Update density kernel during NGWF line search !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MAXIT_PEN"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of penalty functional iterations !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "PEN_PARAM"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Penalty functional parameter !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "R_PRECOND"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Radial cut-off for real-space preconditioner !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "WRITE_NGWF_PLOT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write NGWFs in plotting format !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_NGWF_GRAD_PLOT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write NGWF Gradients in plotting format !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "MAXIT_HOTELLING"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of Hotelling iteration per NGWF change !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "USE_SPACE_FILLING_CURVE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Re-arrange atoms according to space-filling curve !*"
  i=i+1

  KW(i)%LABEL = "COREHAM_DENSKERN_GUESS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Initial guess for density kernel from core Hamiltonian !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "WRITE_DENSKERN"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write density kernel restart information !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "READ_DENSKERN"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Read density kernel restart information !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "NGWF_CG_TYPE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Type of CG coefficient for NGWF optimisation &
      & NGWF_POLAK = Polak-Ribbiere formula; &
      & NGWF_FLETCHER = Fletcher-Reeves formula. !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "LNV_CG_TYPE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Type of CG coefficient for LNV denskern optimisation &
       & LNV_POLAK = Polak-Ribbiere formula; &
       & LNV_FLETCHER = Fletcher-Reeves formula. !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "WRITE_TIGHTBOX_NGWFS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write in universal tightbox NGWFs restart information !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "READ_TIGHTBOX_NGWFS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Read in universal tightbox NGWFs restart information !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_DENSITY_PLOT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write the charge density in plotting format !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "FFTBOX_PREF"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Preferred FFT box dimensions !*"
  i=i+1

  KW(i)%LABEL = "TIMINGS_LEVEL"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Level of timings output: 0(none) 1(procs summary), 2(proc details) !*"
  i=i+1

  KW(i)%LABEL = "CONSTANT_EFIELD"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Cartesian coordinates of constant electric field vector !*"
  i=i+1

  KW(i)%LABEL = "CHARGE"
  KW(i)%TYP = "D:B" ! pa: changed from I to D
  KW(i)%DSCRPT = "*! The total charge of the system !*"
  i=i+1

  KW(i)%LABEL = "WRITE_FORCES"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write ionic forces !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_POSITIONS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write ionic positions each geometry or MD step !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "TASK"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Type of calculation !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "CUTOFF_ENERGY"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Plane wave kinetic energy cutoff !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "NUM_EIGENVALUES"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of energy and occupancy eigenvalues to print &
       &below and above the Fermi level !*"
  i=i+1

  KW(i)%LABEL = "OLD_LNV"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Use LNV algorithm backwards compatible pre Dec 2004 !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "SMOOTH_PROJECTORS"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Halfwidth of Gaussian filter for nonlocal projectors !*"
  i=i+1

  KW(i)%LABEL = "SMOOTH_LOC_PSPOT"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Halfwidth of Gaussian filter for local pseudopotential !*"
  KW(i)%KWGRP = "PSEUDO"
  i=i+1

  KW(i)%LABEL = "OCC_MIX"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Mix fraction of occupancy preconditioned NGWF cov grad !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "NGWF_CG_MAX_STEP"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Maximum length of trial step for NGWF optimisation line search !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "LNV_CG_MAX_STEP"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Maximum length of trial step for kernel optimisation line search !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "GEOM_MODULUS_EST"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! The estimated bulk modulus !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_FREQUENCY_EST"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! The estimated average phonon frequency at the gamma point !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_ENERGY_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Geometry optimization energy convergence tolerance.&
       & The difference between max and min energies over&
       & geom_convergence_win iterations must be less than this !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_FORCE_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Geometry optimization force convergence tolerance !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_DISP_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Geometry optimization displacement convergence tolerance !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_MAX_ITER"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Maximum number of geometry optimization iterations !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_CONVERGENCE_WIN"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Geometry optimization convergence tolerance window.&
       & The geometry optimization convergence criteria must all be met for&
       & geom_convergence_win iterations before acceptance !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_CONTINUATION"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Read information for continuation of a previous geometry optimisation !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_OUTPUT_DETAIL"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Level of output detail for GEOM: BRIEF, NORMAL, VERBOSE, PROLIX !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_BACKUP_ITER"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of geometry optimisation iterations between backups of all data for continuation !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "DELTA_E_CONV"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Use consecutive energy gains as a criterion for NGWF convergence !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "PRINT_QC"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Print Quality Control information !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "MINIT_LNV"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Min number of LNV iterations !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MAX_RESID_HOTELLING"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Max allowed value in Hotelling residual !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "CUBE_FORMAT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Allow .cube format for plot outputs !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "GRD_FORMAT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Allow .grd format for plot outputs !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "POSITIONS_ABS"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Cartesian positions for each atom !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "SPECIES"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Species information (symbol, atomic number, number of NGWFs, NGWF radius) !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "SPECIES_POT"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Pseudopotential name for each species !*"
  KW(i)%KWGRP = "PSEUDO"
  i=i+1

  KW(i)%LABEL = "SPECIES_ATOMIC_SET"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Atomic set name for each species !*"
  i=i+1

  KW(i)%LABEL = "SPECIES_CONSTRAINTS"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Ionic constraints for each species !*"
  i=i+1

  KW(i)%LABEL = "ODD_PSINC_GRID"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Force odd number of points in simcell psinc grid !*"
  i=i+1

  KW(i)%LABEL = "GEOM_PRINT_INV_HESSIAN"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Write inverse Hessian to standard output !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_METHOD"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Method for geometry optimization &
       & CARTESIAN = BFGS from CASTEP; or &
       & BFGS = BFGS from CASTEP; &
       & LBFGS = limited memory BFGS; &
       & TPSD = two-point steepest descent; &
       & DELOCALIZED = DELOCALIZED INTERNALS from CASTEP. !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_LBFGS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Whether to perform LBFGS rather than BFGS in a Geometry Optimization !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "NNHO"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Initialise NGWFs to nonorthogonal natural hybrid orbitals !*"
  KW(i)%KWGRP = "BASIS"
  i=i+1

  KW(i)%LABEL = "ELEC_CG_MAX"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of NGWF iterations to reset CG !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "CHECK_ATOMS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Check atoms on top of each other !*"
  i=i+1

  KW(i)%LABEL = "SPECIES_NGWF_PLOT"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Species whose NGWFs to plot !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "TSSEARCH_METHOD"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Specifies method to be used for TS search (e.g., LSTQST !*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "TSSEARCH_LSTQST_PROTOCOL"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Specifies LSTQST protocol !*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "TSSEARCH_QST_MAX_ITER"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Specifies maximum number of QST steps !*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "TSSEARCH_CG_MAX_ITER"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Specifies maximum number of CG steps !*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "TSSEARCH_FORCE_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Force tolerance for TS search !*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "TSSEARCH_DISP_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Displacement tolerance for TS search !*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "POSITIONS_ABS_PRODUCT"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Cartesian positions for each atom in the product (TS search)!*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "POSITIONS_ABS_INTERMEDIATE"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Cartesian positions for each atom in the intermediate structure (TS search)!*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "NGWF_HALO"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Halo extension to NGWF radii !*"
  KW(i)%KWGRP = "BASIS"
  i=i+1

  KW(i)%LABEL = "HOMO_DENS_PLOT"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Number of squared MOs to plot from HOMO and lower !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "LUMO_DENS_PLOT"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Number of squared MOs to plot from LUMO and higher !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "HOMO_PLOT"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Number of MOs to plot from HOMO and lower !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "LUMO_PLOT"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Number of MOs to plot from LUMO and higher !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "DOS_SMEAR"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Half width of smearing Gaussians for DOS !*"
  i=i+1

  KW(i)%LABEL = "PDOS_MAX_L"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! The maximum azimuthal angular momentum channel to project on to in a pDOS calculation !*"
  i=i+1

  KW(i)%LABEL = "SPIN_POLARIZED"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Switch for spin polarisation !*"
  KW(i)%KWGRP = "SPIN"
  i=i+1

  KW(i)%LABEL = "SPIN_POLARISED"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Switch for spin polarisation !*"
  KW(i)%KWGRP = "SPIN"
  i=i+1

  KW(i)%LABEL = "SPIN"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! Total spin of system !*"
  KW(i)%KWGRP = "SPIN"
  i=i+1

  KW(i)%LABEL = "MAXIT_PALSER_MANO"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of iterations for Palser-Manolopoulos scheme !*"
  i=i+1

  KW(i)%LABEL = "MAXIT_KERNEL_FIX"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum # iterations of Penalty Functional idempotency correction per LNV step !*"
  i=i+1

  KW(i)%LABEL = "DO_PROPERTIES"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Allow calculation of properties !*"
  i=i+1

  KW(i)%LABEL = "POPN_CALCULATE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Allow population analysis !*"
  i=i+1

  KW(i)%LABEL = "POPN_BOND_CUTOFF"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Bond length cutoff for population analysis !*"
  i=i+1

  KW(i)%LABEL = "DEVEL_CODE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! For development code only !*"
  i=i+1

  KW(i)%LABEL = "WRITE_XYZ"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Output coordinates in .xyz file !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_PARAMS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Output runtime parameters at startup !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "DENSE_THRESHOLD"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Threshold for matrix segment filling for segment  to be dense!*"
  i=i+1

  KW(i)%LABEL = "NGWF_ANALYSIS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Perform NGWF analysis !*"
  i=i+1

  KW(i)%LABEL = "SPREAD_CALCULATE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Calculate spread of NGWFs !*"
  i=i+1

  KW(i)%LABEL = "DISPERSION"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Select dispersion correction !*"
  i=i+1

  KW(i)%LABEL = "PADDED_LATTICE_CART"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! The simulation cell lattice vectors for the padded cell !*"
  i=i+1

  KW(i)%LABEL = "COULOMB_CUTOFF_TYPE"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Type of cutoff coulomb interaction: NONE, SPHERE, CYLINDER, SLAB, WIRE !*"
  KW(i)%KWGRP = "CHARGE"
  i=i+1

  KW(i)%LABEL = "COULOMB_CUTOFF_RADIUS"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Radius of sphere or cylinder for cutoff coulomb interaction !*"
  KW(i)%KWGRP = "CHARGE"
  i=i+1

  KW(i)%LABEL = "COULOMB_CUTOFF_LENGTH"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Length of cylinder or width of slab for cutoff coulomb interaction !*"
  KW(i)%KWGRP = "CHARGE"
  i=i+1

  KW(i)%LABEL = "COULOMB_CUTOFF_WRITE_INT"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Write real-space cutoff Coulomb interaction scalarfield !*"
  KW(i)%KWGRP = "CHARGE"
  i=i+1

  KW(i)%LABEL = "LOCPOT_SCHEME"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Scheme for evaluating local potential matrix elements &
      & FULL = Calculate matrix and symmetrize;&
      & LOWER = Calculate lower triangle only and expand;&
      & ALTERNATE = Calculate alternating elements from both triangles and expand !*"
  i=i+1

  KW(i)%LABEL = "LNV_CHECK_TRIAL_STEPS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Check stability of kernel at each trial step during LNV !*"
  i=i+1

  KW(i)%LABEL = "BS_KPOINT_PATH"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! K-point path for bandstructure calculation !*"
  i=i+1

  KW(i)%LABEL = "BS_KPOINT_PATH_SPACING"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! K-point spacing along bandstructure path !*"
  i=i+1

  KW(i)%LABEL = "BS_NUM_EIGENVALUES"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Num of energy and occ. eigenvalues to print below and above the Fermi level from a bs calc !*"
  i=i+1

  KW(i)%LABEL = "BS_UNFOLD"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Number of times to unfold Brillouin zone in each lattice direction !*"
  i=i+1

  KW(i)%LABEL = "BS_METHOD"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Method to use: 'KP' or 'TB' !*"
  i=i+1

  KW(i)%LABEL = "GEOM_REUSE_DK_NGWFS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Re-use density kernel and NGWFs during geometry optimisation steps !*"
  i=i+1

  KW(i)%LABEL = "WRITE_CONVERGED_DK_NGWFS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Only write Density Kernel and NGWFs upon convergence of NGWF optimisation !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_SW_NGWFS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write NGWFs restart information in spherical waves representation !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "READ_SW_NGWFS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Read NGWFs restart information in spherical waves representation !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_MAX_L"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum angular momentum number when writing in SW representation !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "READ_MAX_L"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum angular momentum number when reading in SW representation !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "EXTRA_N_SW"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of zeros of the sph Bessel function for the SW !*"
  i=i+1

  KW(i)%LABEL = "CLASSICAL_INFO"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Include classical atoms Coulomb interaction in the calculation !*"
  i=i+1

  KW(i)%LABEL = "HUBBARD_MAX_ITER"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of DFT+U projector optimisation steps, 0 for none !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "HUBBARD_ENERGY_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Energy tolerance when using DFT+U projector optimisation !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "HUBBARD_CONV_WIN"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Energy convergence window when using DFT+U projector optimisation !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "HUBBARD_PROJ_MIXING"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Proportion of old Hubbard projector to mix with new !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "HUBBARD"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Hubbard species info (symb., ang. mom., U parameter (eV), effective charge, alpha parameter (eV)) !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "HUBBARD_FUNCTIONAL"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! DFT+U energy functional to use !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "HUBBARD_TENSOR_CORR"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! DFT+U projector tensorial correction to use !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "HUBBARD_NGWF_SPIN_THRESHOLD"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! NGWF RMS  gradient threshold at which to switch off DFT+U spin-splitting !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "POLARISATION_CALCULATE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Allow calculation of polarisation !*"
  i=i+1

  KW(i)%LABEL = "KERFIX"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Density kernel fixing approach !*"
  i=i+1

  KW(i)%LABEL = "PAW"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Uses a PAW construction to find correct core densities/wavefunctions !*"
  KW(i)%KWGRP = "PAW"
  i=i+1

  KW(i)%LABEL = "DO_TDDFT"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Allow Time-Dependent DFT calculation !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_MAXIMUM_ENERGY"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Desired maximum of spectrum from TDDFT !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_RESOLUTION"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Desired resolution of spectrum from TDDFT (in Hartree) !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_PROPAGATION_METHOD"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Method used to integrate von Neumann equation eg. RUNGEKUTTA or CRANKNICHOLSON!*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_SPARSITY_LEVEL"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Matrix sparsity when computing propagators e.g. 0(recommended),1,2,3 !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_TAMMDANCOFF"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Invoke Tamm-Dancoff decoupling approximation !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_DIPOLE_KICK_STRENGTH"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Maximum allowed phase shift in TDDFT delta-kick, units of PI !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_XC_FUNCTIONAL"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Exchange-correlation functional for TDDFT &
       & LDA = Adiabatic Perdew-Zunger LDA;&
       & NONE = Random Phase Approximation.!*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_HAMILTONIAN_MIXING"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Order of polynomial extrapolation to H(t + half Delta t) 0,1,2 !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_DAMPING"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Energy smearing when Fourier transforming for frequency-dependent dipole moment !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_ENFORCED_IDEMPOTENCY"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Project out at each timestep that part of change to denskern not respecting idempotency to 1st order !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_MAXIT_HOTELLING"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of Hotelling iteration per propagation step in Crank-Nicholson propagator !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_MAX_RESID_HOTELLING"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Max allowed value in Hotelling residual for Crank-Nicholson propagator !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "TDDFT_INV_OVERLAP_EXACT"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Renew inverse overlap with O N^3 algorithm before beginning TDDFT !*"
  KW(i)%KWGRP = "TDDFT"
  i=i+1

  KW(i)%LABEL = "PBC_CORRECTION_CUTOFF"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! alpha*L cutoff parameter for Martyna-Tuckerman PBC correction !*"
  i=i+1

  KW(i)%LABEL = "SPECIES_LDOS_GROUPS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Species groups for Local density of states calculation !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "DX_FORMAT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Allow .dx format for plot outputs !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "DX_FORMAT_DIGITS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of significant digits in .dx output !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "DX_FORMAT_COARSE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Output only points on the coarse grid !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "PAW_OUTPUT_DETAIL"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Level of output detail for PAW: BRIEF, NORMAL or VERBOSE !*"
  KW(i)%KWGRP = "PAW"
  i=i+1

  KW(i)%LABEL = "FINE_GRID_SCALE"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! Ratio of size of fine grid to standard grid !*"
  i=i+1

  KW(i)%LABEL = "SPECIES_COND"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Species information for Conduction NGWFs (symbol, atomic number, number of NGWFs, NGWF radius) !*"
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "LIBXC_X_FUNC_ID"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Functional identifier for which exchange functional to use with LIBXC - see LIBXC documentation !*"
  KW(i)%KWGRP = "XC"
  i=i+1

  KW(i)%LABEL = "LIBXC_C_FUNC_ID"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Functional identifier for which correlation functional to use with LIBXC - see LIBXC documentation !*"
  KW(i)%KWGRP = "XC"
  i=i+1

  KW(i)%LABEL = "IS_DENSITY_THRESHOLD"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Implicit solvent: rho_0 parameter in Fattebert-Gygi functional !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SOLVATION_BETA"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Implicit solvent: beta parameter in Fattebert-Gygi functional !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_BULK_PERMITTIVITY"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! Implicit solvent: eps_inf parameter(relative permittivity) in Fattebert-Gygi functional !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_MULTIGRID_ERROR_TOL"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! [OBSOLETE] Implicit solvent: stop criterion for the multigrid solver !*"
  KW(i)%KWGRP = "SOLVATION" ! JCW: Keyword kept so that we can detect old-syntax inputs
  i=i+1

  KW(i)%LABEL = "IS_SMEARED_ION_WIDTH"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Implicit solvent: Width of Gaussian smearing !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_IMPLICIT_SOLVENT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Use implicit solvent? !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_CHECK_SOLV_ENERGY_GRAD"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Implicit solvent: sanity check the energy gradient !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SMEARED_ION_REP"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Implicit solvent: use smeared ions for electrostatics !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SOLVATION_METHOD"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Implicit solvent: direct or corrective method !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_DIELECTRIC_MODEL"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Implicit solvent: how the cavity is determined !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_DIELECTRIC_FUNCTION"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Implicit solvent: how the cavity is determined !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SOLVATION_OUTPUT_DETAIL"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Implicit solvent: controls extra output !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SURFACE_THICKNESS"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Implicit solvent: thickness used for SA calculation !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SOLVENT_SURFACE_TENSION"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! [OBSOLETE] Use is_solvent_surf_tension_instead, specify unscaled value !*"
  KW(i)%KWGRP = "SOLVATION" ! jd: Keyword kept so that we can detect old-syntax inputs
  i=i+1

  KW(i)%LABEL = "IS_MULTIGRID_MAX_ITERS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! [OBSOLETE] Implicit solvent: max number of iterations in multigrid solver !*"
  KW(i)%KWGRP = "SOLVATION" ! JCW: Keyword kept so that we can detect old-syntax inputs
  i=i+1

  KW(i)%LABEL = "IS_MULTIGRID_NLEVELS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! [OBSOLETE] Implicit solvent: number of levels in multigrid solver !*"
  KW(i)%KWGRP = "SOLVATION" ! JCW: Keyword kept so that we can detect old-syntax inputs
  i=i+1

  KW(i)%LABEL = "IS_MULTIGRID_DEFECT_ERROR_TOL"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! [OBSOLETE] Implicit solvent: stop criterion for defect correction !*"
  KW(i)%KWGRP = "SOLVATION" ! JCW: Keyword kept so that we can detect old-syntax inputs
  i=i+1

  KW(i)%LABEL = "IS_DISCRETIZATION_ORDER"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! [OBSOLETE] Implicit solvent: discretization order (2nd, 4th, ...) for the PB solver !*"
  KW(i)%KWGRP = "SOLVATION" ! JCW: Keyword kept so that we can detect old-syntax inputs
  i=i+1

  KW(i)%LABEL = "IS_INCLUDE_APOLAR"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Implicit solvent: include apolar terms in solvation energy !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_APOLAR_SASA_DEFINITION"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Implicit solvent: sets the method used to define the surface area for the nonpolar term !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_APOLAR_METHOD"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Implicit solvent: the method by which the apolar contribution is calculated !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_DENSITY_MIN_THRESHOLD"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Implicit solvent: parameter in Andreussi functional !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_DENSITY_MAX_THRESHOLD"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Implicit solvent: parameter in Andreussi functional !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_BC_COARSENESS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Open BCs: controls boundary condition coarse-graining !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_BC_SURFACE_COARSENESS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Open BCs: controls boundary condition coarse-graining !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_BC_THRESHOLD"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Open BCs: controls boundary condition coarse-graining !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "OPENBC_PSPOT"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Force open BCs in local pseudopotential !*"
  KW(i)%KWGRP = "PSEUDO"
  i=i+1

  KW(i)%LABEL = "OPENBC_PSPOT_FINETUNE_F"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Open BCs in local pseudo, fineness parameter!*"
  KW(i)%KWGRP = "PSEUDO"
  i=i+1

  KW(i)%LABEL = "OPENBC_PSPOT_FINETUNE_NPTSX"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Open BCs in local pseudo, npts_x parameter!*"
  KW(i)%KWGRP = "PSEUDO"
  i=i+1

  KW(i)%LABEL = "OPENBC_PSPOT_FINETUNE_ALPHA"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Open BCs in local pseudo, alpha parameter!*"
  KW(i)%KWGRP = "PSEUDO"
  i=i+1

  KW(i)%LABEL = "OPENBC_ION_ION"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Force open BCs in ion-ion energy !*"
  i=i+1

  KW(i)%LABEL = "COND_NUM_STATES"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Number of conduction states to be optimised for !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_READ_DENSKERN"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Read in conduction density kernel !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_READ_TIGHTBOX_NGWFS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Read in universal tightbox conduction NGWFs !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_MAXIT_LNV"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of LNV iterations during conduction NGWF optimisation !*" ! ndmh
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_MINIT_LNV"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Min number of LNV iterations during conduction NGWF optimisation !*" ! ndmh
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_KERNEL_CUTOFF"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Conduction density kernel radius !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_INIT_SHIFT"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Initial shifting factor for projected conduction Hamiltonian !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_SHIFT_BUFFER"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Additional buffer for updating projected Hamiltonian shift !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_FIXED_SHIFT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Fixed projected conduction Hamiltonian shift !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_CALC_MAX_EIGEN"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Calculate maximum conduction Hamiltonian eigenvalue !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_CALC_OPTICAL_SPECTRA"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Calculate matrix elements for optical absorption spectra !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_SPEC_CALC_MOM_MAT_ELS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Calculate momentum matrix elements (default true otherwise use position) !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_SPEC_CALC_NONLOC_COMM"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Calculate nonlocal commutator for momentum matrix elements (default true) !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_SPEC_CONT_DERIV"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Calculate non-local commuator using continuous deriv in k-space (default true) !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_SPEC_NONLOC_COMM_SHIFT"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! Finite difference shift for non-local commutator if calculating using finite difference !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "KERNEL_DIIS_LSHIFT"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! The initial value of Beta in the level shifting matrix !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "KERNEL_DIIS_LS_ITER"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum Level shifting iteration !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "KERNEL_DIIS_SIZE"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of kernels saved during kernel DIIS !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "KERNEL_DIIS_MAXIT"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of kernel DIIS iterations !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "KERNEL_DIIS_SCHEME"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Select scheme for kernel-diis !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "KERNEL_DIIS_THRESHOLD"
  KW(i)%TYP = "D:A"
  KW(i)%DSCRPT = "*! Density mixing convergence threshold !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "KERNEL_DIIS_LINEAR_ITER"
  KW(i)%TYP = "I:A"
  KW(i)%DSCRPT = "*! Number of linear mixing iterations !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "KERNEL_DIIS_COEFF"
  KW(i)%TYP = "D:A"
  KW(i)%DSCRPT = "*! Coefficient for the input kernel in linear mixing DIIS !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "KERNEL_DIIS_CONV_CRITERIA"
  KW(i)%TYP = "T:A"
  KW(i)%DSCRPT = "*! Density mixing convergence criteria !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "KERNEL_CHRISTOFFEL_UPDATE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Update density kernel during NGWF line search !*"
  i=i+1

  KW(i)%LABEL = "ELEC_ENERGY_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Tolerance on total energy change during NGWF optimisation !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "ELEC_FORCE_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Tolerance on max force change during NGWF optimisation !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "NGWF_MAX_GRAD"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Maximum permissible value of NGWF Gradient for convergence !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "NONSC_FORCES"
  KW(i)%TYP = "L:A"
  KW(i)%DSCRPT = "*! Include non self-consistent forces due to NGWF optimisation !*"
  i=i+1

  KW(i)%LABEL = "GEOM_RESET_DK_NGWFS_ITER"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of geom iterations between resets of kernel and NGWFs !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "OPENBC_HARTREE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Force open BCs in Hartree potential !*"
  i=i+1

  KW(i)%LABEL = "VDW_PARAMS"
  KW(i)%TYP = "B:E"
  KW(i)%DSCRPT = "*! Replacement VDW parameters (atomic number, c6coeff, radzero, neff) !*"
  KW(i)%KWGRP = "VDW"
  i=i+1

  KW(i)%LABEL = "VDW_DCOEFF"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Replacement VDW damping coefficient !*"
  KW(i)%KWGRP = "VDW"
  i=i+1

  KW(i)%LABEL = "VDW_BC"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! 3 character string defining BCs for van der Waals interactions along &
                   &each lattice vector. 'O' for open, 'P' for periodic. !*"
  KW(i)%KWGRP = "VDW"
  i=i+1

  KW(i)%LABEL = "VDW_RADIAL_CUTOFF"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Radial cutoff for van der Waals interactions !*"
  KW(i)%KWGRP = "VDW"
  i=i+1

  KW(i)%LABEL = "COMMS_GROUP_SIZE"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of procs in a group (determines comms efficiency) !*"
  i=i+1

  KW(i)%LABEL = "DBL_GRID_SCALE"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! Ratio of charge density / potential working grid to standard grid (1 or 2 only) !*"
  i=i+1

  KW(i)%LABEL = "REALSPACE_PROJECTORS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Whether to evaluate and store projectors in real space !*"
  i=i+1

  KW(i)%LABEL = "EVEN_PSINC_GRID"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Force even number of points in simcell psinc grid !*"
  i=i+1

  KW(i)%LABEL = "COND_PLOT_JOINT_ORBITALS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Plot orbitals in the joint basis following a conduction calculation !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_PLOT_VC_ORBITALS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Plot orbitals in the val and cond NGWF basis sets after a cond calc !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_NUM_EXTRA_STATES"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of extra conduction states for initial 'preconditioning' !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_NUM_EXTRA_ITS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of NGWF iterations with extra conduction states for 'preconditioning' !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "WRITE_NGWF_RADIAL"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Write NGWFs radial distributions !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_NGWF_GRAD_RADIAL"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Write NGWFs gradients radial distributions !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_RADIAL_STEP"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Define the grid step used for writing radial distributions !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_RADIAL_SMEAR"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Define the gaussian smearing used for writing radial distributions !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "SMOOTH_SCHEME"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Smoothing scheme for the NGWF gradients at the edges !*"
  i=i+1

  KW(i)%LABEL = "R_SMOOTH"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Radius of the unshaved NGWF gradients !*"
  i=i+1

  KW(i)%LABEL = "K_SMOOTH"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Characteristic wavevector of the smoothing function(NGWF gradient in reciprocal space!*"
  i=i+1

  KW(i)%LABEL = "AUG_FUNCS_RECIP"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Construct Augmentation functions in recip space (T) or real (F) !*"
  i=i+1

  KW(i)%LABEL = "INITIAL_DENS_REALSPACE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Construct initial density in real space from atomsolver density !*"
  i=i+1

  ! lpl: NB62O output
  KW(i)%LABEL = "WRITE_NBO"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Performs Natural Population Analysis and writes &
       &a FILE.47 input for GENNBO !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NBO_WRITE_SPECIES"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! List of atoms to be included in the output to &
       &GENNBO FILE.47 !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NBO_SPECIES_NGWFLABEL"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! User-specified NGWF (false) lm-label and &
       &NMB/NRBs for GENNBO !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NBO_INIT_LCLOWDIN"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Performs atom-centered Lowdin symmetric &
       &orthogonalization in generating the NAOs. !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NBO_WRITE_LCLOWDIN"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Write a GENNBO FILE.47 containing all the atoms &
       &in the atom-centered Lowdin basis to satisfy the strict &
       &lm-orthogonality requirement in GENNBO !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NBO_WRITE_NPACOMP"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Writes individual NAO population into the &
       &standard output. !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NBO_SCALE_DM"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Scales density matrix in the FILE.47 output to &
       &achieve charge integrality (Required for proper GENNBO &
       &functionality). !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NLPP_FOR_EXCHANGE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Give exchange matrix same sparsity as non-local pseudopotential matrix !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "SWRI"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Defines spherical-wave resolutions of identity !*"
  KW(i)%KWGRP = "SWRI"
  i=i+1

  KW(i)%LABEL = "HFX_USE_RI"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! ID of the SWRI to use for HFx !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "HFX_MAX_L"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Maximum order of HFx expansion into SWs !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "NBO_AOPNAO_SCHEME"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! AO to PNAO scheme to use in generating NAOs (for testing purposes). !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NBO_LIST_PLOTNBO"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! List of NBOs to be plotted according to GENNBO &
          &output indices. !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "PLOT_NBO"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Plot NBO's orbitals from FILE.xx as defined by nbo_plot_orbtype. !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NBO_PLOT_ORBTYPE"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Type of GENNBO-generated orbital to plot. !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "NBO_PNAO_ANALYSIS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! S/P/D/F CHARACTER ANALYSIS ON PNAO. !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "PHONON_FINITE_DISP"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Amplitude of the ionic perturbation to be used in a finite displacement phonon calculation !*"
  i=i+1

  KW(i)%LABEL = "PHONON_FMAX"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Maximum force allowed on the unperturbed system for a phonon calculation !*"
  i=i+1

  KW(i)%LABEL = "PHONON_FARMING_TASK"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Operation to perform for phonon calc (for task farming or post-proc. of dynamical matrix) !*"
  i=i+1

  KW(i)%LABEL = "PHONON_DISP_LIST"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! List of displacements to perform for phonon calculation !*"
  i=i+1

  KW(i)%LABEL = "PHONON_TMIN"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Lower bound of temperature range for computation of vibrational thermodynamic quantities !*"
  i=i+1

  KW(i)%LABEL = "PHONON_TMAX"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Upper bound of temperature range for computation of vibrational thermodynamic quantities !*"
  i=i+1

  KW(i)%LABEL = "PHONON_DELTAT"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Temperature step for computation of vibrational thermodynamic quantities !*"
  i=i+1

  KW(i)%LABEL = "PHONON_MIN_FREQ"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Discard phonon frequencies smaller than this value for computation of vibrational &
          &thermodynamic quantities !*"
  i=i+1

  KW(i)%LABEL = "GEOM_LBFGS_MAX_UPDATES"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of LBFGS update vectors to store !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_LBFGS_BLOCK_LENGTH"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! How many updates to store before reallocation in an unbounded LBFGS calculation !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "SPECIES_CORE_WF"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Core Wavefunction filename for each species !*"
  i=i+1

  KW(i)%LABEL = "COND_SPEC_PRINT_MAT_ELS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write optical matrix elements to file !*"
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_SPEC_SCISSOR_OP"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Scissor operator for JDOS and imag. diel. fn. !*"
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_SPEC_OPT_SMEAR"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Half width of smearing Gaussians for JDOS and imag. diel. fn. !*"
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_CALC_EELS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Calculate matrix elements for electron energy loss spectra (EELS) !*" ! lr408
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "EXTERNAL_PRESSURE"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! External applied pressure !*"
  i=i+1

  KW(i)%LABEL = "SMOOTHING_FACTOR"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! Smoothing factor in volume term !*"
  i=i+1

  KW(i)%LABEL = "ISOSURFACE_CUTOFF"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! Isosurface cutoff in volume term !*"
  i=i+1

  ! gibo: cDFT keyword:START
  KW(i)%LABEL = "CONSTRAINED_DFT"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Constrained_DFT species info [1:symb., 2:l=angul. mom., 3. Ionic_charge, 4: U_occ. (eV), 5:U_q_up (eV),&
   & 6:U_q_down (eV), 7:U_spin (eV), 8: Targeted N_up, 9: Targeted N_down, 10:targeted N_up - N_down !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL = "CI_CDFT"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Perform a CONFIGURATION-INTERACTION CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CI_CDFT_NUM_CONF"
  KW(i)%TYP  = "I:I"
  KW(i)%DSCRPT  = "*! Number of cDFT-configurations for CI_CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_ATOM_CHARGE"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Perform an ATOM-CHARGE-constrained CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_ATOM_SPIN"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Perform an ATOM-SPIN-constrained CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_CHARGE_ACCEPTOR"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Perform an ACCEPTOR GROUP-CHARGE-constrained CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_CHARGE_DONOR"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Perform a DONOR GROUP-CHARGE-constrained CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_SPIN_ACCEPTOR"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Perform an ACCEPTOR GROUP-SPIN-constrained CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_SPIN_DONOR"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Perform a DONOR GROUP-SPIN-constrained CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_CHARGE_DIFF"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Perform a GROUP-CHARGE-DIFFERENCE-constrained CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_SPIN_DIFF"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Perform a GROUP-SPIN-DIFFERENCE-constrained CDFT simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_HUBBARD"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Perform a constrained-DFT+U simulation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_CHARGE_ACCEPTOR_TARGET"
  KW(i)%TYP  = "D:I"
  KW(i)%DSCRPT  = "*! Targeted group-CHARGE for GROUP-CHARGE-ACCEPTOR-constrained cDFT !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_CHARGE_DONOR_TARGET"
  KW(i)%TYP  = "D:I"
  KW(i)%DSCRPT  = "*! Targeted group-CHARGE for GROUP-CHARGE-DONOR-constrained cDFT !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_SPIN_ACCEPTOR_TARGET"
  KW(i)%TYP  = "D:I"
  KW(i)%DSCRPT  = "*! Targeted group-SPIN for GROUP-SPIN-ACCEPTOR-constrained cDFT !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_SPIN_DONOR_TARGET"
  KW(i)%TYP  = "D:I"
  KW(i)%DSCRPT  = "*! Targeted group-SPIN for GROUP-SPIN-DONOR-constrained cDFT !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_CHARGE_DIFF_TARGET"
  KW(i)%TYP  = "D:I"
  KW(i)%DSCRPT  = "*! Targeted CHARGE difference (acceptor-donor) for  GROUP-CHARGE-DIFFERENCE-constrained cDFT !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_SPIN_DIFF_TARGET"
  KW(i)%TYP  = "D:I"
  KW(i)%DSCRPT  = "*! Targeted SPIN difference (acceptor-donor) for  GROUP-SPIN-DIFFERENCE-constrained cDFT !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "MAXIT_CDFT_U_CG"
  KW(i)%TYP  = "I:I"
  KW(i)%DSCRPT  = "*! Max number of cdFT-U conjugate gradients (CG)&
   & iterations !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_CG_TYPE"
  KW(i)%TYP  = "T:E"
  KW(i)%DSCRPT  = "*! Type of CG coefficient for CDFT U-optimisation &
   & NGWF_POLAK = Polak-Ribbiere formula; &
   & NGWF_FLETCHER = Fletcher-Reeves formula.!*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_CG_THRESHOLD"
  KW(i)%TYP  = "D:I"
  KW(i)%DSCRPT  = "*! RMS gradient convergence threshold for U-pot. in CDFT !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_CG_MAX"
  KW(i)%TYP  = "I:E"
  KW(i)%DSCRPT  = "*! Number of U-opt iterations to reset CG !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_MAX_GRAD"
  KW(i)%TYP  = "D:I"
  KW(i)%DSCRPT  = "*! Maximum permissible value of CDFT U-Gradient for convergence !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_ELEC_ENERGY_TOL"
  KW(i)%TYP  = "P:I"
  KW(i)%DSCRPT  = "*! Tolerance on total energy change during CDFT optimisation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_CG_MAX_STEP"
  KW(i)%TYP  = "D:E"
  KW(i)%DSCRPT  = "*! Maximum length of trial step for cDFT optimisation line search !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GURU"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Let the user signal she/he does not need helpt with the cDFT U-initialisation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_CONTINUATION"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Logical to restart (from the *.cdft file) a cDFT U-optimisation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_TRIAL_LENGTH"
  KW(i)%TYP  = "D:I"
  KW(i)%DSCRPT  = "*! Trial length for cDFT line-search !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_PRINT_ALL_OCC"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Logical to have the occupancy-matrix of all the CDFT-atoms printed in stdout !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_TIGHT"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Logical to activate tight NGWFs-cDFT optimisation !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_READ_PROJECTORS"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Logical to read cDFT-projectors from file !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_CHARGE_UP_ONLY"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Logical to constrain only UP electrons !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_GROUP_CHARGE_DOWN_ONLY"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Logical to constrain only UP electrons !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL  = "CDFT_MULTI_PROJ"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Logical to use mutiple angular-momentum projectors on one cDFT-site !*"
  KW(i)%KWGRP = "CDFT"
  ! gibo: cDFT keyword:END
  i=i+1

  !ep: Mermin Keywords section
  KW(i)%LABEL = "MERMIN"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Use mermin method to optimise the kernel !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "CHECK_MERMIN"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Check kernel search direction during mermin optimisation!*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  !ep: still to be tested confirmed
  KW(i)%LABEL = "MERMIN_TEMP"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Use easy annealing during mermin !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MERMIN_SMEARING_WIDTH"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Occupancy smearing width in MERMIN calculations !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MERMIN_CHEB"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Chebyshev expansion to be used in mermin mod !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MAXIT_MERMIN"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maxit Mermin cycle iteraions !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MERMIN_THRESHOLD_ORIG"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Threshold for the lnv loop when doing extrapolation !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MERMIN_CG_TYPE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Type of CG coefficient for MERMIN denskern optimisation &
       & MERMIN_POLAK = Polak-Ribbiere formula; &
       & MERMIN_FLETCHER = Fletcher-Reeves formula. !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MERMIN_CG_MAX_STEP"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Maximum length of trial step for kernel optimisation &
       & line search !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MERMIN_CHECK_TRIAL_STEPS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Check stability of kernel at each trial step during LNV !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MERMIN_ROUND_EVALS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Round MERMIN eigenvalues to N decimal figures !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MERMIN_FREE_ENERGY_THRES"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Convergence threshold for the free energy &
       & in MERMIN calculations !*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1

  KW(i)%LABEL = "MERMIN_MU_SQ"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Check kernel search direction during mermin &
       & optimisation!*"
  KW(i)%KWGRP = "MERMIN"
  i=i+1
  !ep: Mermin Keywords section

  KW(i)%LABEL = "SPECIES_ATOMIC_SET_COND"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Atomic set description for each species, for initialising Conduction NGWFs !*"
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "SPECIES_ATOMIC_SET_AUX"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Atomic set description for each species, for initialising Auxiliary NGWFs !*"
  i=i+1

  KW(i)%LABEL = "SPECIES_AUX"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Species information for Auxiliary NGWFs (symbol, atomic number, number of NGWFs, NGWF radius) !*"
  i=i+1

  KW(i)%LABEL = "HFX_NLPP_FOR_EXCHANGE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Give exchange matrix same sparsity as non-local pseudopotential matrix !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "HFX_MAX_Q"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Maximum number of Bessel zeros in HFx''s SW expansion !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "DMA_BESSEL_AVERAGING"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Whether DMA expansion should average over even-odd Bessels !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "DMA_USE_RI"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! ID of the SWRI to use for DMA !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "THOLE_POLARISABILITIES"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Thole polarisabilities of all QM atoms !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_SMEARING_A"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Thole big A for short-range smearing of undamped !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "SWRI_CHEB_BATCHSIZE"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of SW pot expansions buffered !*"
  KW(i)%KWGRP = "SWRI"
  i=i+1

  KW(i)%LABEL = "SWRI_SWOP_SMOOTHING"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Apply SW/SWpot smoothing for more accurate overlaps? !*"
  KW(i)%KWGRP = "SWRI"
  i=i+1

  KW(i)%LABEL = "SWRI_VERBOSE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Verbose output for spherical wave resolution of identity? !*"
  KW(i)%KWGRP = "SWRI"
  i=i+1

  KW(i)%LABEL = "SWRI_PROXIMITY_SORT_POINT"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! SWRI atomblocks will be processed closest-first to this point !*"
  KW(i)%KWGRP = "SWRI"
  i=i+1

  KW(i)%LABEL = "POL_EMB_FIXED_CHARGE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Is the embedding a fixed-charge (non-polarisable) FF? !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_PERM_SCALING"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Scaling factor applied to interactions with MM perm. mpoles !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "HFX_READ_XMATRIX"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Read restart information for the X matrix !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "HFX_WRITE_XMATRIX"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Write restart information for the X matrix !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "SWRI_OVERLAP_INDIRECT"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Inversions done for overlap metric !*"
  KW(i)%KWGRP = "SWRI"
  i=i+1

  KW(i)%LABEL = "SWRI_IMPROVE_INVERSE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Use Hotelling improvement when calculating inverses !*"
  KW(i)%KWGRP = "SWRI"
  i=i+1

  KW(i)%LABEL = "HFX_METRIC"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Electrostatic or overlap metric !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "SWX_C_THRESHOLD"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Absolute magnitude below which expansion coefficients will be zeroed !*"
  KW(i)%KWGRP = "SWX"
  i=i+1

  KW(i)%LABEL = "CACHE_LIMIT_FOR_SWOPS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Max cache size for SWs or SWpots in PPDs (in MiB) !*"
  KW(i)%KWGRP = "SWX"
  i=i+1

  KW(i)%LABEL = "DMA_MULTIPOLE_SCALING"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Scaling factor applied to all DMA multipoles !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "CACHE_LIMIT_FOR_EXPANSIONS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Max cache size for expanded potentials (in MiB) !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "MG_CONTINUE_ON_ERROR"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! If true, solutions DL_MG will not abort on errors !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "IS_CORE_WIDTH"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Implicit solvent: Radius where eps is set to unity !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "HFX_CUTOFF"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! HFx cutoff radius!*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  ! lpl: More NBO
  KW(i)%LABEL = "NBO_WRITE_DIPOLE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Writes dipole matrix to FILE.47 !*"
  i=i+1

  KW(i)%LABEL = "NBO_SCALE_SPIN"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Whether or not partial density matrices for different spins are scaled independently !*"
  KW(i)%KWGRP = "NBO"
  i=i+1

  KW(i)%LABEL = "PROJECTORS_PRECALCULATE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Whether to pre-calculate the nonlocal projectors in FFTboxes rather than on-the-fly !*"
  i=i+1

  KW(i)%LABEL = "EIGENSOLVER_ORFAC"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Precision to which the parallel eigensolver will orthogonalise evecs !*"
  i=i+1

  KW(i)%LABEL = "EIGENSOLVER_ABSTOL"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Precision to which the parallel eigensolver will calculate the eigenvalues !*"
  i=i+1

  KW(i)%LABEL = "READ_HAMILTONIAN"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Read current Hamiltonian matrix from a file !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "WRITE_HAMILTONIAN"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Save current Hamiltonian matrix in a file !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "THERMOSTAT"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Thermostat for MD in NVT ensemble !*"
  i=i+1

  KW(i)%LABEL = "VELOCITIES"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! Initial velocities for each atom !*"
  i=i+1

  KW(i)%LABEL = "MD_DELTA_T"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Molecular dynamics time step !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_NUM_ITER"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Maximum number of molecular dynamics iterations !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_RESTART"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Restart MD from backup files !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_PROPERTIES"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Compute vibrational and IR spectra from MD !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_RESET_HISTORY"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Reset mixing scheme for initial guess of NGWFs and density kernel !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "KERNEL_FORCE_CONV"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Force density kernel convergence on last NGWFs optimization step !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_DKN_NUM"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of coefficients used to build new guess for dkn !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_NGWFS_NUM"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of coefficients used to build new guess for NGWFS !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_DKN_TYPE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Type of mixing used to build new guess for dkn !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_NGWFS_TYPE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Type of mixing used to build new guess for NGWFS !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_LOCAL_LENGTH"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Max radius for local mixing of NGWFs !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_LOCAL_SMEAR"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Radial smearing for local mixing of NGWFs !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_DKN_INIT_TYPE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Type of initialization before MIX_DKN_TYPE  !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_NGWFS_INIT_TYPE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Type of initialization before MIX_DKN_TYPE  !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_DKN_INIT_NUM"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of steps before extrapolation of the desity kernel !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_NGWFS_INIT_NUM"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of steps before extrapolation of the NGWFs !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_DKN_RESET"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of extrapolation steps between two resets of the density kernel !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_NGWFS_RESET"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of extrapolation steps between two resets of the NGWFs !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MTS_XI"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Internal thermostat in the multiple time-step scheme !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MTS_NSTEP"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of time steps in the multiple time-step scheme !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MTS_NGWF_THRESHOLD"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! NGWF convergence threshold for the mts correction !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MTS_LNV_THRESHOLD"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! LNV convergence threshold for the mts correction !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "ETRANS_SETUP"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Transport setup description !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_LCR"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Compute electronic transport coefficients (LCR config) !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_BULK"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Compute electronic transport coefficients (Bulk config) !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_LEADS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Define the transport leads !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_SAME_LEADS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Use same description for all leads !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_ECMPLX"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Imaginary part of the energy !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_EMAX"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Highest energy for the calculation of transmission coefficients !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_EMIN"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Lowest energy for the calculation of transmission coefficients !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_ENUM"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of energy points for the calculation of transmission coefficients !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_WRITE_HS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Write hamiltonian corresponding to transport setup !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "NGWF_CG_ROTATE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Rotate density kernel during NGWF optimization !*"
  i=i+1

  ! ars: ensemble DFT
  KW(i)%LABEL = "EDFT"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Use Ensemble-DFT method to optimise the kernel !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_MAXIT"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of EDFT iterations to optimise the kernel !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_MAX_STEP"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Maximum step during EDFT line search !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_SMEARING_WIDTH"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Occupancy smearing width in EDFT calculations !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_FREE_ENERGY_THRES"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Convergence threshold for the free energy in EDFT calculations !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_ENERGY_THRES"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Convergence threshold for the energy in EDFT calculations !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_ENTROPY_THRES"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Convergence threshold for the entropy in EDFT calculations !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_RMS_GRADIENT_THRES"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Convergence threshold for the RMS gradient in EDFT calculations !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_COMMUTATOR_THRES"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Convergence threshold for the RMS commutator in EDFT calculations !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_FERMI_THRES"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Convergence threshold for the Fermi level in EDFT calculations !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_NELEC_THRES"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Convergence threshold for number of electrons in grand canonical edft !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1


  KW(i)%LABEL = "EDFT_WRITE_OCC"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Generate file with smeared occupancies !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_EXTRA_BANDS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of additional MOs with non-zero occupancy in EDFT !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_ROUND_EVALS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Round EDFT eigenvalues to N decimal figures !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "MAXIT_KERNEL_OCC_CHECK"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of kernel resets after occupancy checks !*"
  i=i+1

 ! tjz07: TDDFT params
  KW(i)%LABEL = "LR_TDDFT_CALCULATE"    !tjz07
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT =  "*! enables LR-TDDFT calculation !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_NUM_STATES" !tjz07
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Sets the number of excitation energies we want to solve for!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_CG_THRESHOLD" !tjz07
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! sets convergence tolerance for CG routine !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_MAXIT_CG" !tjz07
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! sets the maximum number of iterations for CG routine !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_MAXIT_PEN" !tjz07
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! sets the maximum number of iterations for penalty functional routine !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_PENALTY_TOL" !tjz07
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! sets the convergence tolerance for the Penalty functional routine !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_RESET_CG" !tjz07
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! sets the number of iterations after which the search direction gets reset !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_KERNEL_CUTOFF" ! tjz07
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! sets cutoff on the effective response density kernel!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_WRITE_DENSITIES" ! tjz07
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Determines whether to write out TDDFT response densities!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_WRITE_KERNELS"
  KW(i)%TYP ="L:E"
  KW(i)%DSCRPT = "*! writes out response kernels after each iteration!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_RESTART"
  KW(i)%TYP ="L:E"
  KW(i)%DSCRPT = "*! Restart flag for the LR_TDDFT option !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_TRIPLET"
  KW(i)%TYP ="L:E"
  KW(i)%DSCRPT ="*! DEtermines if triplet states are calculated or singlets!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_PROJECTOR"
  KW(i)%TYP ="L:E"
  KW(i)%DSCRPT = "*! Use projector onto unoccupied subspace !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_JOINT_SET"
  KW(i)%TYP ="L:E"
  KW(i)%DSCRPT ="*! Use joint set to represent cond states!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_PREOPT"
  KW(i)%TYP ="L:E"
  KW(i)%DSCRPT ="*! Refine starting guess in preoptimisation routine!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_PREOPT_ITER"
  KW(i)%TYP ="I:E"
  KW(i)%DSCRPT ="*! Number of preoptimisation iterations!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_ANALYSIS"
  KW(i)%TYP ="L:E"
  KW(i)%DSCRPT ="*! Do a full O(N^3) analysis of TDDFT transitions!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="THREADS_MAX"
  KW(i)%TYP ="I:E"
  KW(i)%DSCRPT = "*! Number of threads in outer loops !*"
  KW(i)%KWGRP = "THREADS"
  i=i+1

  KW(i)%LABEL ="THREADS_PER_FFTBOX"
  KW(i)%TYP ="I:E"
  KW(i)%DSCRPT = "*! Number of nested threads used for FFT box operations. !*"
  KW(i)%KWGRP = "THREADS"
  i=i+1

  KW(i)%LABEL = "FFTBOX_BATCH_SIZE"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of NGWFs in fftbox batches !*"
  KW(i)%KWGRP = "THREADS"
  i=i+1

  KW(i)%LABEL ="THREADS_NUM_FFTBOXES"
  KW(i)%TYP ="I:E"
  KW(i)%DSCRPT = "*! Number of threads to use in OpenMP-parallel FFTs  !*"
  KW(i)%KWGRP = "THREADS"
  i=i+1

  KW(i)%LABEL ="THREADS_PER_CELLFFT"
  KW(i)%TYP ="I:E"
  KW(i)%DSCRPT = "*! Number of threads to use in OpenMP-parallel FFTs on simulation cell !*"
  KW(i)%KWGRP = "THREADS"
  i=i+1

  KW(i)%LABEL ="POSITIONS_FRAC"
  KW(i)%TYP ="B:B"
  KW(i)%DSCRPT = "*! Fractional positions of atomic species !*"
  i=i+1

  KW(i)%LABEL ="POSITIONS_XYZ_FILE"
  KW(i)%TYP ="T:B"
  KW(i)%DSCRPT = "*! .xyz file to read positional data from !*"
  KW(i)%KWGRP = "CELLDATA"
  i=i+1

  KW(i)%LABEL ="LATTICE_ABC"
  KW(i)%TYP ="B:B"
  KW(i)%DSCRPT = "*! The simulation cell vectors !*"
  KW(i)%KWGRP = "CELLDATA"
  i=i+1

  KW(i)%LABEL = "PADDED_LATTICE_ABC"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! The simulation cell lattice vectors for the padded cell !*"
  KW(i)%KWGRP = "CELLDATA"
  i=i+1

  KW(i)%LABEL = "ESDF_DUMP"     !ncor
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Dump all runtime parameters at startup !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL ="RMS_KERNEL_MEASURE"
  KW(i)%TYP ="L:E"
  KW(i)%DSCRPT = "*! Use root mean squared measure of [K,H] commutator and delta K !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "WRITE_POLARISATION_PLOT"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Write the polarisation density in plotting format !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "HUBBARDSCF_ON_THE_FLY"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Carry out on-the-fly HUBBARDSCF with new projectors for new NGWFs !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "THREADS_NUM_MKL"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of threads to use in OpenMP-parallel MKL operations !*"
  KW(i)%KWGRP = "THREADS"
  i=i+1

  KW(i)%LABEL = "MTS_MIX_INC"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Include the mts correction step in the NGWFs and dkn mixing scheme !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MTS_ELEC_ENERGY_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Tolerance on total energy change during NGWF optimisation !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MTS_ELEC_FORCE_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Tolerance on max force change during NGWF optimisation !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MTS_NGWF_MAX_GRAD"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Maximum permissible value of NGWF Gradient for convergence !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MIX_NGWFS_COEFF"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Mix the propagated NGWFs with the new NGWFs !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MTS_MAXIT_LNV"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of LNV iterations !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "MTS_MINIT_LNV"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Min number of LNV iterations !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MTS_MAXIT_PEN"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of penalty iterations !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MTS_MAXIT_NGWF_CG"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of conjugate gradients iterations !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "PHONON_SAMPLING"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Default number of sampling points for finite difference calculation (1 or 2) !*"
  i=i+1

  KW(i)%LABEL = "PHONON_VIB_FREE"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Default allowed vibrational degrees of freedom for all ions !*"
  i=i+1

  KW(i)%LABEL = "PHONON_EXCEPTION_LIST"
  KW(i)%TYP = "B:E"
  KW(i)%DSCRPT = "*! List of ionic degrees of freedom with modified properties !*"
  i=i+1

  KW(i)%LABEL = "PHONON_ENERGY_CHECK"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Check total energy of system doesn't decrease upon ionic displacement !*"
  i=i+1

  KW(i)%LABEL = "PHONON_WRITE_EIGENVECS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Write phonon mode eigenvectors to file !*"
  i=i+1

  KW(i)%LABEL = "PHONON_ANIMATE_LIST"
  KW(i)%TYP = "B:E"
  KW(i)%DSCRPT = "*! List of phonon modes to write as xyz animations !*"
  i=i+1

  KW(i)%LABEL = "PHONON_ANIMATE_SCALE"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Scaling factor for phonon mode animations !*"
  i=i+1

  KW(i)%LABEL = "SUPERCELL"
  KW(i)%TYP = "B:E"
  KW(i)%DSCRPT = "*! Definition of supercell !*"
  i=i+1

  KW(i)%LABEL = "PHONON_GRID"
  KW(i)%TYP = "B:E"
  KW(i)%DSCRPT = "*! Regular grid of q points used for calculating vibrational thermodynamic quantities !*"
  i=i+1

  KW(i)%LABEL = "PHONON_QPOINTS"
  KW(i)%TYP = "B:E"
  KW(i)%DSCRPT = "*! List of additional q points to calculate !*"
  i=i+1

  KW(i)%LABEL = "PHONON_DOS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Calculate phonon DOS from MP grid and write to file !*"
  i=i+1

  KW(i)%LABEL = "PHONON_DOS_MIN"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Lower bound of frequency range (in cm^-1) for phonon DOS calculation !*"
  i=i+1

  KW(i)%LABEL = "PHONON_DOS_MAX"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Upper bound of frequency range (in cm^-1) for phonon DOS calculation !*"
  i=i+1

  KW(i)%LABEL = "PHONON_DOS_DELTA"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Frequency step (in cm^-1) for phonon DOS calculation !*"
  i=i+1

  KW(i)%LABEL = "PHONON_SK"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Use Slater-Koster style interpolation for q points instead of real-space cutoff !*"
  i=i+1

  KW(i)%LABEL  = "HUBBARD_READ_PROJECTORS"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Logical to read Hubbard-projectors from file !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL  = "CDFT_WRITE_POTENTIALS"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Logical to write cDFT-potentials into file !*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL = "RUN_TIME"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! The maximum allocated run time for this job (in seconds) !*"
  i=i+1

  KW(i)%LABEL = "PSEUDO_PATH"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Path to pseudopotentials !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "ETRANS_LEAD_NKPOINTS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of kpoints to use to determine lead chemical potential !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_CALCULATE_LEAD_MU"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Calculate the lead chemical potentials !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_EREF"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Reference energy for electronic transport calculation !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ZERO_TOTAL_FORCE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Subtract avg force to ensure Newton's 3rd law holds !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "WRITE_VELOCITIES"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write ionic velocities each MD step !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "MD_OUTPUT_DETAIL"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Level of output detail for MD: BRIEF, NORMAL or VERBOSE !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "POSITIONS_FRAC_PRODUCT"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT= "*! Fractional positions for each atom in the product (TS search)!*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "POSITIONS_FRAC_INTERMEDIATE"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Fractional positions for each atom in the intermediate structure (TS search)!*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "POSITIONS_INTERMEDIATE_XYZ_FILE"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! .xyz file to read positional data for intermediate structure (TS search)!*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL = "POSITIONS_PRODUCT_XYZ_FILE"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! .xyz file to read positional data for product structure (TS search)!*"
  KW(i)%KWGRP = "TS"
  i=i+1

  ! lpl: DDEC keywords
  KW(i)%LABEL = "DDEC_CALCULATE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Performs DDEC charge analysis !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_WRITE_RAD"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Write DDEC partial radial density for each atom !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_CLASSICAL_HIRSHFELD"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! DDEC range of IH (+/-) ionic reference states to be generated !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_RAD_NPTS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of spherical shells per atom !*"
  i=i+1

  KW(i)%LABEL = "DDEC_RAD_RCUT"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Radius of largest spherical shell !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_MAXIT"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of DDEC iterations !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_CORE_MAXIT"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of DDEC core density iterations !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_CONV_THRESHOLD"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! DDEC charge convergence threshold !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_USE_COREDENS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Whether to include core densities in calculation !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_IH_FRACTION"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! DDEC IH fraction !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_ZERO_THRESHOLD"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! DDEC threshold to neglect 1/density !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_CORE_CORRECTION"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Whether to correct the integrated number of core electrons on the Cartesian grid !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_RESHAPE_DENS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Whether to reshape the partitioned AIM density after each iteration !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_CORE_CORR_MAXIT"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Number of core correction iterations !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_IH_IONIC_RANGE"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! DDEC range of IH (+/-) ionic reference states to be generated !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_REFDENS_PATH"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Path to DDEC reference densities !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_REFDENS_INIT"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Initialize ISA guess densities as neutral reference densities !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_RENORMALIZE_REFDENS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Renormalize reference densities !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_C3_REFDENS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Reshape reference densities to produce c3 reference densities !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_RAD_SHELL_MODE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! The effective radius of each shell !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_MULTIPOLE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Compute DDEC dipoles and quadrupoles !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_MOMENT"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Compute DDEC moment !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_INTERP_RAD_DENS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Interpolate converged radial densities for a smoother profile !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_AVG_RAD"
  KW(i)%TYP = "L:D"
  KW(i)%DSCRPT = "*! Compute expected radius of each atom based on the partitoned density !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_MIN_SHELL_DENS"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Minimum number of points per shell for grid point binning !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_REF_SHELL_MODE"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Mode for initializing reference densities from fine radial grid !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_RCOMP"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! DDEC rcomp block paramters !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_FORMAT_DENS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Whether to format the input densities to shave off density spikes !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_EFF_DECAY_EXP"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Calculate AIM effective decay exponents for r > ddec_eff_decay_rmin !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_EFF_DECAY_RMIN"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Minumum radius of AIM density to which effective decay exponents are fitted !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1
  ! lpl: END DDEC keywords

  ! vv: anharmonic keywords
  KW(i)%LABEL = "ANHARMONIC_CALCULATE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Active the calculation of the IR spectrum !*"
  i=i+1

  KW(i)%LABEL = "ANH_QC_FACTOR"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Quantum correction factor in IR spectrum !*"
  i=i+1

  KW(i)%LABEL = "ANH_ACF_FACTOR"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Prefactor for the autocorrelation function !*"
  i=i+1

  KW(i)%LABEL = "ANH_FIRST_ITER"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! First md iteration to include in the autocorrelation !*"
  i=i+1

  KW(i)%LABEL = "ANH_LAST_ITER"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Last md iteration to include in the autocorrelation !*"
  i=i+1

  KW(i)%LABEL = "ANH_MD_TEMP"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Temperature in the md simulation !*"
  i=i+1

  KW(i)%LABEL = "ANH_PLOT_FIRSTFREQ"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! First freq to be shown in the IR table !*"
  i=i+1

  KW(i)%LABEL = "ANH_PLOT_LASTFREQ"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Last freq to be shown in the IR table !*"
  i=i+1

  KW(i)%LABEL = "ANH_PLOT_ALL"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Plot the whole IR spectrum !*"
  i=i+1

  KW(i)%LABEL = "ANH_APPLY_FILTER"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Apply the gaussian filter !*"
  i=i+1

  KW(i)%LABEL = "ANH_TYPE"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Describe the type of calculation to perform !*"
 ! vv: end anharmonic keywords
  i=i+1

  KW(i)%LABEL = "ETRANS_LEAD_DISP_TOL"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! The maximum acceptable difference in atomic positions between a lead and its first principle layer !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "CACHE_LIMIT_FOR_DKNBLKS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Max cache size for remote DKN blocks (in MiB) !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "TIMINGS_ORDER"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Sorting order of timings !*"
  i=i+1

  KW(i)%LABEL = "DMA_CALCULATE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Compute distributed multipole analysis !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "DMA_MAX_L"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Maximum order of DMA multipoles for properties !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "DMA_OUTPUT_POTENTIAL"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Output distributed multipole potential on cell faces !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "DMA_OUTPUT_POTENTIAL_REFERENCE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Output reference (pointwise) potential on cell faces !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "DMA_METRIC"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Electrostatic or overlap metric !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "DMA_MAX_Q"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Maximum number of Bessel zeros in DMA''s SW expansion !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_PRECOND"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! If set to T, TDDFT search direction gets preconditioned !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_PRECOND_TOL"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Convergence tolerance in applying the preconditioner !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_PRECOND_ITER"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Max number of iterations in applying the preconditioner !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_NUM_CONV_STATES"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT ="*! Sets the number of already converged states.!*"
  KW(i)%KWGRP ="LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_CHECK_CONV_ITER"
  KW(i)%TYP ="I:E"
  KW(i)%DSCRPT ="*! Num of iterations at which conv of states is checked !*"
  KW(i)%KWGRP ="LRTDDFT"
  i=i+1

  KW(i)%LABEL = "SPECIES_TDDFT_KERNEL"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Species groups defining the region in which the TDDFT kernel is defined !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "ETRANS_NUM_EIGCHAN"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! The number of transmission eigenchannels to calculate !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_PLOT_EIGCHAN"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Plot the transmission eigenchannels defined in etrans_eigenchannel_energies block !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_PLOT_EIGCHAN_ENERGIES"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! The energies at which the transmission eigenchannels are to be plotted !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_WRITE_XYZ"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write the lead and device co-ordinates to .xyz files !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_LEAD_SIZE_CHECK"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Check if leads define a full principle layer !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL = "ETRANS_SEED_LEAD"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! The seed lead for determining the tri-diagonal partitioning !*"
  KW(i)%KWGRP = "ETRANS"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_HOMO_NUM"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT ="*! Defines number of occ KS transitions considered in TDDFT analysis.!*"
  KW(i)%KWGRP ="LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_LUMO_NUM"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT ="*! Defines number of unocc KS transitions considered in TDDFT analysis.!*"
  KW(i)%KWGRP ="LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_XC_FINITE_DIFF"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT ="*! Evaluate fxc using finite difference technique.!*"
  KW(i)%KWGRP ="LRTDDFT"
  i=i+1

  KW(i)%LABEL = "QNTO_ANALYSIS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT =  "*! Runs QNTO analysis !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "QNTO_NUM_TRANSITION"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT =  "*! Sets the number of NTO pairs to output or compare !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "QNTO_WRITE_ORBITALS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT =  "*! Enables plotting NTOs !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "QNTO_SVD_METHOD"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT =  "*! Sets the SVD method !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "IS_PBE"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Implicit solvent: include ionic strengths? !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_PBE_TEMPERATURE"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Implicit solvent: temperature for Boltzmann term !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SC_STERIC_MAGNITUDE"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Implicit solvent: softcore steric potential prefactor !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SC_STERIC_SMOOTHING_ALPHA"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Implicit solvent: softcore steric pot erf parameter !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SC_STERIC_CUTOFF"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Implicit solvent: cutoff rad for softcore steric pot !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_STERIC_WRITE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Implicit solvent: write steric pot to file? !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "SOL_IONS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Ions in solvent: name, charge, concentration !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "LOWDIN_POPN_CALCULATE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Allow Lowdin population analysis !*"
  i=i+1

  KW(i)%LABEL = "TURN_OFF_HARTREE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Omit Hartree terms !*"
  i=i+1

  KW(i)%LABEL = "HUBBARD_UNIFY_SITES"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Combine all projectors into one Hubbard site !*"
  i=i+1

  KW(i)%LABEL = "DFT_NU"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Constrained DFT+nu Species info!*"
  KW(i)%KWGRP = "CDFT"
  i=i+1

  KW(i)%LABEL = "HUBBARD_COMPUTE_U_OR_J"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Compute U or J correction !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "KERNEL_TRACK_MID_OCC"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Print middle occupancy after LNV convergence !*"
  i=i+1

  KW(i)%LABEL = "KERNEL_CHECK_ALL"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Turn on all kernel checking parameters !*"
  i=i+1

  KW(i)%LABEL = "HUBBARD_J_MINORITY_TERM"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Include minority-only energy term in DFT+U+J !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "POLARISATION_SIMCELL_CALCULATE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Calculate simcell polarisation (also, quadrupoles) !*"
  i=i+1

  KW(i)%LABEL = "POLARISATION_SIMCELL_REFPT"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Reference point for simcell dipoles, quadrupoles !*"
  i=i+1

  KW(i)%LABEL = "POLARISATION_LOCAL"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Allow the calculation of local polarisation !*"
  i=i+1

  KW(i)%LABEL = "IS_PBE_BC_DEBYE_SCREENING"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Should solvated BCs use Debye lambda exp() factor !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_PBE_USE_FAS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! [OBSOLETE] Should FAS be used for non-linear PBE !*"
  KW(i)%KWGRP = "SOLVATION" ! JCW: Keyword kept so that we can detect old-syntax inputs
  i=i+1

  KW(i)%LABEL = "IS_PBE_EXP_CAP"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Cap of exponential argument in Boltzmann term !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "PRINT_POTENTIAL_NOXC"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Print Local potential withouth XC (only Hartree+Ion) !*"
  i=i+1

  KW(i)%LABEL = "IS_STERIC_POT_TYPE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Implicit solvent: steric potential type !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_MULTIGRID_VERBOSE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Implicit solvent: verbose multigrid output? !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_MULTIGRID_VERBOSE_Y"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Implicit solvent: x-section Y for verbose output !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_MULTIGRID_VERBOSE_Z"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Implicit solvent: x-section Z for verbose output !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_MULTIGRID_ERROR_DAMPING"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! [OBSOLETE] Implicit solvent: error damping in defect correction? !*"
  KW(i)%KWGRP = "SOLVATION" ! JCW: Keyword kept so that we can detect old-syntax inputs
  i=i+1

  KW(i)%LABEL = "PDOS_OPTADOS_OUTPUT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Output angular momentum projected density of states weights for input into OptaDOS !*"
  i=i+1

  KW(i)%LABEL = "COUPLINGS_STATES"
  KW(i)%TYP = "B:E"
  KW(i)%DSCRPT = "*! Allow the calculation of electronic couplings !*"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_SUBSYSTEM_COUPLING"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Compute coupling between subsystems in LRTDDFT !*"
  KW(i)%KWGRP ="LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_OPTICAL_PERMITTIVITY"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Optical permittivity of solvent used in SCF response in TDDFT !*"
  KW(i)%KWGRP ="LRTDDFT"
  i=i+1

  ! agreco: add keyword to use extended NGWFs
  KW(i)%LABEL = "EXTEND_NGWF"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! directions along which NGWFs are extended !*"
  i=i+1

  ! agreco: add keyword to request NGWFs with initial random values
  KW(i)%LABEL = "FULL_RAND_NGWF"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! request NGWFs initialised with random values !*"
  i=i+1

  ! agreco: whether or not write NGWFs for plotting at every iteration
  KW(i)%LABEL = "WRITE_NGWF_PLOT_EVERY_IT"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Write NGWFs in plotting format at every iteration !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "DMA_SCALE_CHARGE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Scale DMA monopoles? !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "WRITE_INITIAL_RADIAL_NGWFS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Controls output of radial NGWF plots from atomsolver !*"
  i=i+1

  KW(i)%LABEL = "ETRANS_EREF_METHOD"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! The method used to determine the reference energy for electronic transport !*"
  i=i+1

  KW(i)%LABEL = "IS_AUTO_SOLVATION"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! If true, vacuum calculation will automatically precede solvated calculation !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SEPARATE_RESTART_FILES"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! If true, dielectric cavity will be constructed from a second set of restart files !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_RPA"
  KW(i)%TYP ="L:E"
  KW(i)%DSCRPT = "*! If true, perform a full LR_TDDFT calculation rather than the Tamm-Dancoff approx.!*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "SPECIES_LOCDIPOLE_GROUPS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Species groups for calculation of dipole moments of subsystems !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "DMA_TARGET_NUM_VAL_ELEC"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Target number of valence electrons in DMA region (for scaling) !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_RESTART_FROM_TDA"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Option to restart RPA calculation from Tamm-Dancoff response kernels !*"
  KW(i)%KWGRP ="LRTDDFT"
  i=i+1

  KW(i)%LABEL = "COND_EELS_FINE_PROJECTORS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Directly generate core wavefunctions on the fine grid !*" !ewt23
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "COND_EELS_REALSPACE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Compute matrix elements for EELS spectra using realspace method !*" !ewt23
  KW(i)%KWGRP = "COND"
  i=i+1

  KW(i)%LABEL = "POL_EMB_POT_FILENAME"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! File with multipoles and energy terms for polarisable embedding potential !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "DMA_PRECISE_GDMA_OUTPUT"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Extra precision in GDMA output at the cost of breaking format compat !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_MPOLE_EXCLUSION_RADIUS"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Exclusion radius for point multipole singularities !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  ! ja531 --> Fermi Operator Expansion
  KW(i)%LABEL = "FOE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Use the Fermi Operator expansion to evaluate density kernels !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  !gom
  KW(i)%LABEL  = "DFT_NU_OPT_U1_ONLY"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Selective optimisazion of U1 potential !*"
  KW(i)%KWGRP = "DFTnu"
  i=i+1

  KW(i)%LABEL  = "DFT_NU_OPT_U2_ONLY"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Selective optimisazion of U2 potential !*"
  KW(i)%KWGRP = "DFTnu"
  i=i+1

  KW(i)%LABEL  = "DFT_NU_CONTINUATION"
  KW(i)%TYP  = "L:I"
  KW(i)%DSCRPT  = "*! Restart U1/2 optimisation from .dft_nu file !*"
  KW(i)%KWGRP = "DFTnu"
  i=i+1

  KW(i)%LABEL = "SPECIES_PDOS_GROUPS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Species groups for Local, angular momentum projected density of states calculation !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  ! gcc32 : confined ngwfs
  KW(i)%LABEL = "CONFINED_NGWFS_BARRIER"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! The barrier potential (in Ha) for NGWF confinement !*"
  i=i+1

  KW(i)%LABEL = "CONFINED_NGWFS"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Whether to use the hybrid confinement method !*"
  i=i+1

  KW(i)%LABEL = "MAXIT_NGWF_CG_CONFINED"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Number of iterations for which the NGWFs are explicitly &
                       &confined !*"
  i=i+1

  KW(i)%LABEL = "POPN_MULLIKEN_PARTIAL"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Enable Mulliken partial charge analysis !*"
  i=i+1

  ! tjz21: New keywords controlling automatic cond initialisation
  KW(i)%LABEL = "COND_ENERGY_RANGE"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Energy range of optimised cond states measured from HOMO !*"
  i=i+1

  KW(i)%LABEL ="COND_ENERGY_GAP"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT ="*! Energy gap between highest optimised and lowest unoptimised cond state !*"
  i=i+1

  KW(i)%LABEL = "MD_WRITE_HISTORY"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Write mixing scheme for initial guess of NGWFs and density kernel !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_GLOBAL_RESTART"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Option to restart the md calculation with electronic history !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_WRITE_OUT"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Makes MD restart output cleaner !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "POL_EMB_PAIRWISE_POLARISABILITY"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Pairwise polarisability for emulating Thole damping !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_THOLE_A"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Thole constant (a) for emulating Thole damping !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "SWRI_ASSEMBLY_PREFIX"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Directory+rootname for assembling [VO]matrix blocks !*"
  KW(i)%KWGRP = "SWRI"
  i=i+1

  ! ja531 -> New pDOS approximation to speed things up.
  KW(i)%LABEL = "PDOS_REDUCE_SWS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Reduce bessel functions in pDOS before projection !*"
  i=i+1

  ! mjsp: EDA keywords
  KW(i)%LABEL = "EDA_IATM"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! List of the number atoms in each monomer. Eg. '3 4' " &
        // "would indicate the first 3 atoms in " &
        // "monomer one, and the final 4 atoms in monomer two. !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_DELTADENS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Write the delta densities (electron density differences) " &
        // "of the EDA energy components !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_FRAG_ISOL_POL"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Enables fragment polarisation calculations to be " &
        // "obtained for each fragment independently. !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_FRAG_ISOL_CT"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Calculate delocalisation energies " &
        // "for each fragment pair (state initialised from the frozen stage). " &
        // "NOTE: This calculation also allows relaxation of the surrounding fragments.!*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_RESET_NGWFS_POL"
  KW(i)%TYP = "L:A"
  KW(i)%DSCRPT = "*! Resets the NGWFs for the (full) polarisation stage of " &
        // "the EDA to the initial guess NGWFs. !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_RESET_NGWFS_CT"
  KW(i)%TYP = "L:A"
  KW(i)%DSCRPT = "*! Resets the NGWFs for the charge transfer stage of the " &
        // "EDA to the initial guess NGWFs. !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_READ_FRAGS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Enables the reading of fragments specified using " &
        // "the EDA_FRAGS block. !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_FRAGS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! List of the fragment density kernels and NGWFs prefix " &
        // "(e.g. 'frag1' for frag1.dkn and frag1.tightbox_ngwfs) !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_READ_SUPER"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Enables the reading of supermolecule dkn and NGWFs using " &
        // "the EDA_SUPER block. !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_SUPER"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! The supermolecule data filename prefix (e.g. 'super1' for " &
        // "super1.dkn_supermolecule and super1.tightbox_ngwfs_supermolecule) !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_DELOC"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! List of fragment pairs to calculate delocalisations for !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_SPLIT_ATOMS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Enables partitioning of fragments by atom-splitting !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_FRAG_ATOMS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT =  "*! (Developmental) Fragment to be split for defining " &
        // "bond-splitted fragments in EDA (frag#, atom#(not in use), NGWF#(not in use)) !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

 ! gcc32: lr_phonons keywords
  KW(i)%LABEL = "LR_PHONONS_CALCULATE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT =  "*! enables LR-PHONONS calculation !*"
  KW(i)%KWGRP = "LRPHONONS"
  i=i+1

  ! tjz21: TDDFT keywords specifying charge transfer in tddft kernel
  KW(i)%LABEL = "LR_TDDFT_CT_LENGTH"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT =  "*! Charge-tranfer length for definition of the TDDFT" &
       // "response matrix !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "SPECIES_TDDFT_CT"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Species groups defining the region in which the TDDFT" &
       // "ct is defined !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "MD_AUX_DKN_T"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Target temperature of the auxiliary density kernel !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_AUX_BEREN_TC"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Set the value for the relaxation time in the Berendsen coupling !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_RESTART_THERMO"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Restart MD from .thermo.restart file !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_LNV_THRESHOLD"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Threshold for the lnv loop when doing extrapolation !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_NGWF_THRESHOLD"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Threshold for the outer loop when doing extrapolation !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_INIT_RANDOM"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT =  "*! Set whether initial TDDFT vectors are set to random matrices " &
       // "or pure KS transitions with minimum energies !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_INIT_MAX_OVERLAP"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! If set to T, initialise to KS transitions that maximise the " &
       // "elec-hole overlap, ie. not charge transfer states !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "CHECK_STACK_SIZE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Check if stack size is sufficient? Set to F for valgrind. !*"
  KW(i)%KWGRP = ""
  i=i+1

  KW(i)%LABEL = "EDFT_INIT_MAXIT"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Number of iterations of Ensemble DFT (cubic scaling) to run before LNV. !*"
  KW(i)%KWGRP = ""
  i=i+1

  KW(i)%LABEL = "MW_TOTAL_FORCE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Subtract mass-weighted average force to ensure Newton's 3rd law holds !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "HUBBARD_TENSOR_FORCES"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Calculate force contributions due to Hubbard metric tensor !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "EDA_WRITE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Write EDA continuation data !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EDA_CONTINUATION"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Read information for continuation of a previous EDA calculation !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  ! ebl: DMFT keywords
  KW(i)%LABEL = "DMFT_FULLY_SC"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! If true uses the self energy in the ONETEP kernel NGWF optimization, in the energy functional !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_KERNEL"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Writes dmft density kernel, 0:does not calculate it, -1:writes the purified dmft density kernel !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_FULLY_SC_H"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! If true the H used in the DFT is the KS H built from 1 shot DMFT !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_NBO"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! ---Missing description--- !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_NKPOINTS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of K points for averaging the lattice Green Function !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_OPTICS"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Calculate the optical conductivity from the DMFT Green function !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_ORDER_PROJ"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Small number to enforce the orbital order of the NGWFS when the overlap with projectors are identical !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_PLOT_REAL_SPACE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Write out real space DMFT quantities !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_POINTS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of DMFT energy points on real or matsubara axis !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_PURIFY_SC"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Use the purified kernel for DFT_DMFT in the property module !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SC"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! if true will run the self consistent ONETEP+DMFT, if false, one shot ONETEP+DMFT is used !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SPOIL_KERNEL"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! If false will NOT update LHXC potential !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SWITCH_OFF_PROJ_ORDER"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! For compatibility with older version of ONETEP-DMFT, if true switches" &
     // "off the natural order of the projections !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "HUBBARD_PROJ_READ_ONLY"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Read, but do not write hubbard projectors !*"
  KW(i)%KWGRP = "HUBBARD"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_MLWF_ANALYSIS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! If set to T, a maximally localised wannier function " &
       // "analysis of the converged TDDFT evecs is performed !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_MOM_MAT_ELS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Compute oscillator strengths in momentum rather than " &
       // "position space  !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "LR_TDDFT_PENALTY_FUNC"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! If set to F, the idempotency violation through " &
       // "a penalty functional is not computed, and no iterative improvements are performed !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "IS_INCLUDE_CAVITATION"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! [OBSOLETE] Use is_include_apolar instead !*"
  KW(i)%KWGRP = "SOLVATION" ! jd: Keyword kept so that we can detect old-syntax inputs
  i=i+1

  KW(i)%LABEL = "DMA_DIPOLE_SCALING"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Scaling factor applied to all DMA dipoles !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "DMA_QUADRUPOLE_SCALING"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Scaling factor applied to all DMA quadrupoles !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "RAND_SEED"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Seed for generating velocities in MD from Maxwell-Boltzmann distribution !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "MD_AUTOCORR"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Keyword to output real and auxiliary kernels to external files during&
       & Niklasson/Berendsen propagation !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "SWX_OUTPUT_DETAIL"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Level of output detail for SWX: BRIEF, NORMAL or VERBOSE !*"
  KW(i)%KWGRP = "SWX"
  i=i+1

  KW(i)%LABEL = "TURN_OFF_EWALD"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Omit Ewald term from energies and forces !*"
  i=i+1

  KW(i)%LABEL = "EDA_NODIAG"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Set to avoid the EDA diagonalisation bottlenecks. " &
       // "0=Off 1=Frozen, 2=Polarisation, 3=Frozen and Polarisation. !*"
  KW(i)%KWGRP = "EDA"
  i=i+1

  KW(i)%LABEL = "EFIELD_CALCULATE"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Calculate and output electric field during properties !*"
  i=i+1

  KW(i)%LABEL = "IS_SOLVENT_SURF_TENSION"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! [OBSOLETE] Use is_solvent_surf_tension_instead, specify unscaled value !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SOLVENT_PRESSURE"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! Adjust the pressure used for the SAV apolar model !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_APOLAR_SCALING_FACTOR"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Implicit solvent: Scaling factor for apolar term !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  ! agreco: specify use of complex NGWFs from input file
  KW(i)%LABEL = "USE_CMPLX_NGWFS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! specify to use complex valued NGWFs !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  ! agreco: position of the discontinuity of the sawtooth potential
  ! associated to the constant electric field
  KW(i)%LABEL = "EFIELD_ORIGIN"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Cartesian coordinates of the origin of the sawtooth potential !*"
  i=i+1

  ! agreco: list of k-points to be used for BZ sampling in scf calculation
  KW(i)%LABEL = "KPOINT_LIST"
  KW(i)%TYP = "B:B"
  KW(i)%DSCRPT = "*! K-point list to be used for BZ sampling in scf calculation !*"
  i=i+1

   ! agreco: seed to be used for random number generation
  KW(i)%LABEL = "RAND_SEED_NGWF"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! User specified seed to start pseudo-random generator !*"
  i=i+1

  ! agreco: sigma of the normal distribution used to generate random numbers
  ! this will be used to apply a small perturbation to the initial NGWFs,
  ! breaking the initial symmetry but potentially allowing to check we do not
  ! get stuck in local minima
  KW(i)%LABEL = "RAND_NORMAL_SIGMA"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! User specified sigma for normal distributed random numbers !*"
  i=i+1

  ! agreco: if true, convolute fully random initialised NGWFs with erfc(r),
  ! where r is the distance from the NGWF centre; potentially useful to smooth
  ! random oscillations at the edges of the localisation region
  KW(i)%LABEL = "CONVOLUTE_RAND"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! If true, convolute randomly initialised NGWFs to smooth edges !*"
  i=i+1

  ! agreco: convolution function to use when convolute_rand is true
  KW(i)%LABEL = "CONVOLUTE_FUNC"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Convolute function to use with randomly initialised NGWFs !*"
  i=i+1

  ! agreco: save the overlap matrix in a format suitable for plotting
  ! this keyword is added mostly for some testing of overlap patterns when
  ! using partially extended NGWFs, may be removed in the future.....
  KW(i)%LABEL = "SHOW_OVERLAP"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Save to file overlap matrix at the end of the calculation !*"
  i=i+1

  ! agreco: possibility to apply convolution to smooth oscillations from randomly
  ! initialised NGWFs only in a small region starting from the sphere edge
  KW(i)%LABEL = "CONV_REGION_WIDTH"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Specify width of region (from the sphere edge) in which apply convolution !*"
  i=i+1

  ! agreco: single k-point for testing, in fractional coordinates
  ! to be removed when full k-point sampling is working
  KW(i)%LABEL = "SINGLE_KPOINT"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Single k-point for testing of k-point sampling !*"
  i=i+1

  ! agreco: keyword to specify k-point sample method to be used
  ! essentially a duplicate of BS_METHOD, kept to avoid confusion;
  ! consider to merge them at the end
  KW(i)%LABEL = "KPOINT_METHOD"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Method for sampling of BZ !*"
  i=i+1

  ! agreco: keyword to check a quantity is real when using complex NGWFs;
  ! this is the threshold against which the imaginary part of the quantity
  ! must be checked
  KW(i)%LABEL = "IMAG_THR"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Threshold to accept value of imaginary part of a quantity !*"
  i=i+1

  ! agreco: keyword to multiply all the initial NGWFs by the same complex phase,
  ! used to test implementation of complex NGWFs. The parameter specifies the
  ! phase theta such that the NGWFs are multiplied by exp(i * theta)
  KW(i)%LABEL = "MULT_NGWF_BY_PHASE"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Phase theta added to initial complex NGWFs !*"
  i=i+1

  ! agreco: keyword to multiply all the initial NGWFs by a random complex phase,
  ! between 0 and 2PI, used to test implementation of complex NGWFs.
  KW(i)%LABEL = "MULT_NGWF_BY_RANDOM_PHASE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! If true, random phase in [0,2PI] is added to initial complex NGWFs !*"
  i=i+1

  ! kkbd: Keyword to control the number of ONETEP images.
  KW(i)%LABEL = "NUM_IMAGES"
  KW(i)%TYP   = "I:I"
  KW(i)%DSCRPT = "*! Control the number of ONETEP images !*"
  i=i+1

  ! kkbd: Keyword to control the size of ONETEP images.
  KW(i)%LABEL  = "IMAGE_SIZES"
  KW(i)%TYP    = "T:E"
  KW(i)%DSCRPT = "*! MPI Process size of each ONETEP image separated by pipes | !*"
  i=i+1

  KW(i)%LABEL = "XC_MINTAU" ! JCW
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! The minimum threshold for tau (the kinetic energy " &
       // "density): this is typically used to determine the cutoff where " &
       // "expressions containing 1/tau in the XC energy and potential " &
       // "functions are set to zero, to avoid numerical issues !*"
  KW(i)%KWGRP = "XC"
  i=i+1

  KW(i)%LABEL = "XC_INITIAL_FUNCTIONAL" ! JCW
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Use an alternative XC functional when calculating " &
       // "XC energy for initial guess. Only available when XC_FUNCTIONAL " &
       // "is a meta-GGA. Default: 'PBE'. To omit XC when computing initial " &
       // "guess set to 'NONE'. !*"
  KW(i)%KWGRP = "XC"
  i=i+1

  KW(i)%LABEL = "MD_AUX_REP"
  KW(i)%TYP   = "T:E"
  KW(i)%DSCRPT = "*! Specify the representation of the auxiliary density matrix "&
       // " for the XLBOMD !*"
  KW(i)%KWGRP = "MD"
  i=i+1

  KW(i)%LABEL = "BS_PERTURBATIVE_SOC"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Add perturbative spin-orbit couplings " &
       // "to the bandstructure calculation. !*"
  KW(i)%KWGRP = "BS"
  i=i+1

  KW(i)%LABEL = "HHF_NSTATES"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Number of extra occupied states for HHF calculation !*"
  KW(i)%KWGRP = "HHF"
  i=i+1

  KW(i)%LABEL = "POL_EMB_QMSTAR"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Should the QM* rep be used for any QM/MM interaction !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POLARISATION_BERRY"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Allow calculation of polarisation using Berry phase!*"
  i=i+1

  KW(i)%LABEL = "LR_PHONONS_ZERO_DIM"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Ensure correction of dynamical matrix for molecules!*"
  i=i+1

  KW(i)%LABEL = "LR_PHONONS_KERNEL_CUTOFF" ! gcc32
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! sets cutoff on the effective response density kernel!*"
  i=i+1

  KW(i)%LABEL = "LR_PHONONS_RESTART"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Restart from previously written force constants!*"
  i=i+1

  KW(i)%LABEL = "SWX_DBL_GRID"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Should spherical-wave expansion use double grid? !*"
  KW(i)%KWGRP = "SWX"
  i=i+1

  KW(i)%LABEL = "CACHE_LIMIT_FOR_SWOPS2"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Max cache size for SWs or SWpots at points (in MiB) !*"
  KW(i)%KWGRP = "SWX"
  i=i+1

  KW(i)%LABEL = "POL_EMB_DBL_GRID"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Should polarisable embedding do gradients on double grid? !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "READ_REAL_NGWFS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Read real NGWFs from file into complex NGWFs !*"
  i=i+1

  KW(i)%LABEL = "CHECK_HERMITIAN_MATS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Check hermitian character of complex H/S/K matrices !*"
  i=i+1

  KW(i)%LABEL = "KPOINT_GRID_SIZE"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! size of the k-point grid !*"
  i=i+1

  KW(i)%LABEL = "KPOINT_GRID_SHIFT"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! shift of the k-point grid !*"
  i=i+1

  KW(i)%LABEL = "POL_EMB_REPULSIVE_MM_POT_ALPHA"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! OBSOLETE Pol-emb repulsive MM pot: alpha parameter !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_REPULSIVE_MM_POT_BETA"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! OBSOLETE Pol-emb repulsive MM pot: beta parameter !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_REPULSIVE_MM_POT_A"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! OBSOLETE Pol-emb repulsive MM pot: a parameter !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_REPULSIVE_MM_POT_B"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! OBSOLETE Pol-emb repulsive MM pot: b parameter !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_REPULSIVE_MM_POT_CUTOFF"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Pol-emb repulsive MM pot: cutoff rad around MM atom !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_REPULSIVE_MM_POT_WRITE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Pol-emb repulsive MM pot: write to file? !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_REPULSIVE_MM_POT_C"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! OBSOLETE Pol-emb repulsive MM pot: c parameter !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "ENERGY_COMPONENTS_INTERVAL"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! How often to print out energy components !*"
  i=i+1

  KW(i)%LABEL = "CHECK_DENSITY"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Check density is real when using complex NGWFs !*"
  i=i+1

  ! ja531 --> Fermi Operator Expansion (dense mode)
  KW(i)%LABEL = "DENSE_FOE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Use the dense matrix version of the FOE !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "POL_EMB_DMA_MIN_L"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Minimum order of DMA multipoles for polarisable embedding !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "POL_EMB_DMA_MAX_L"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Maximum order of DMA multipoles for polarisable embedding !*"
  KW(i)%KWGRP = "DMA"
  i=i+1

  KW(i)%LABEL = "POL_EMB_VACUUM_DMA_MIN_L"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Minimum order of DMA multipoles for polarisable embedding (vacuum calc) !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_VACUUM_DMA_MAX_L"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Maximum order of DMA multipoles for polarisable embedding (vacuum calc) !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_WRITE_VACUUM_RESTART"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Should QM/MM restart files be written out !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_VACUUM_QMSTAR"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Should the QM* rep be used for any QM/MM interaction !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  ! gcc32: keywords for spectral function (band-structure) unfolding
  KW(i)%LABEL = "SPECIES_BSUNFLD_GROUPS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Species groups for spectral function unfolding &
       &calculation !*"
  KW(i)%KWGRP = "BS_UNFOLDING"
  i=i+1

  KW(i)%LABEL = "BSUNFLD_KPOINT_PATH"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Primitive-cell k-point path for bandstructure &
       &unfolding calculation !*"
  KW(i)%KWGRP = "BS_UNFOLDING"
  i=i+1

  KW(i)%LABEL = "BSUNFLD_TRANSFORMATION"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Transformation matrix (flattened) between &
       &primitive-cell and supercell lattice vectors !*"
  KW(i)%KWGRP = "BS_UNFOLDING"
  i=i+1

  KW(i)%LABEL = "BSUNFLD_CALCULATE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Logical flag for bandstructure unfolding calculation !*"
  KW(i)%KWGRP = "BS_UNFOLDING"
  i=i+1

  KW(i)%LABEL = "BSUNFLD_NUM_KPTS_PATH"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of primitive-cell kpts sampled along each path !*"
  KW(i)%KWGRP = "BS_UNFOLDING"
  i=i+1

  KW(i)%LABEL = "BSUNFLD_NUM_ATOMS_PRIM"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of atoms in implicit primitive-cell !*"
  KW(i)%KWGRP = "BS_UNFOLDING"
  i=i+1

  KW(i)%LABEL = "BSUNFLD_RESTART"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Restart a bs-unfolding calculation !*"
  KW(i)%KWGRP = "BS_UNFOLDING"
  i=i+1

  KW(i)%LABEL = "BSUNFLD_NUM_EIGENVALUES"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Enforce provided number of kpts per path !*"
  KW(i)%KWGRP = "BS_UNFOLDING"
  i=i+1

  KW(i)%LABEL = "HT_STASH_SIZE"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! OMP stash size for hash tables (in MiB) !*"
  KW(i)%KWGRP = "HT"
  i=i+1

  ! ebl: keywords for DMFT functionality
  KW(i)%LABEL = "DMFT_COMPLEX_FREQ"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Perform DMFT using complex frequencies. &
                    &Real frequencies are required for DOS and optics &
                    &calculations !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_TEMP"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Temperature (in Hartree) for DFT+DMFT!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_EMIN"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Miniumum energy on real axis (Ha) for DFT+DMFT!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_EMAX"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Maximum energy on real axis (Ha) for DFT+DMFT!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SMEAR_SHIFT"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Frequency dependent smearing (Ha) for DFT+DMFT!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SMEAR"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Smearing (in Ha) for DFT+DMFT!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_PARAMAGNETIC"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Imposes the paramagnetic state in DFT+DMFT!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_CUTOFF_SMALL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Atomic energy threshold for double-precision calculations of the DMFT &
                   & Green-function at low frequencies. (Useful for GPU calculations.) !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_ROTATE_GREEN"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Write out the Green function in the rotated local atomic basis !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_OPTICS_WINDOW"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*!Window of energy around Fermi energy for optical conductivity !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_OPTICS_I1"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! First direction for optical conductivity current-current correlator !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_OPTICS_I2"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Second direction for optical conductivity current-current correlator !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_DOS_MIN"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Minimum of energy window on real axis (Ha) for DMFT DOS calculations !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_DOS_MAX"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Maximum of energy window on real axis (Ha) for DMFT DOS calculations !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SMEAR_T"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Frequency dependent smearing parameter T (Ha) for DFT+DMFT. &
                   &See documentation for details of the smearing!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SMEAR_W"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Frequency dependent smearing parameter w (Ha) for DFT+DMFT. &
                   &See documentation for details of the smearing!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SMEAR_ETA"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Frequency dependent smearing parameter eta (Ha) for DFT+DMFT. &
                   &See documentation for details of the smearing!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_NMU_LOOP"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of iterations of the Newtons method for finding the chemical potential!*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_KERNEL_MIX"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Mixing of DMFT and DFT kernels for DFT+DMFT !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_MU_DIFF_MAX"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! The threshold for the difference between the current and target occupancies &
                   & when updating the chemical potential in DMFT !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_MU_ORDER"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Order of the Newton Method used to fix the chemical potential in DMFT &
                   &(HouseHolder general form) !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SKIP_ENERGY"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Do not compute the energy along the Newton steepest descent used to find &
                   &the chemical potential !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_KPOINTS_SYM"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! If true and using k points, additional cubic symmetry " &
        // "is used to reduce the number of k points, and not only the k " &
        // "inversion symmetry. Note: this is NOT valid when computing the " &
        // "DMFT density kernel !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_WIN"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Window of energies around the chemical potential considered for the Green's &
                   & function inversion in DMFT !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_NVAL"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of eigenstates in this energy window used for the Green's &
                   &function inversion !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SCALING_CUTOFF"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SCALING_NMPI"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SCALING_METH"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_SCALING_TAIL"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_KS_SHIFT"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! if true adds the correction to the energy in the SC dmft minimization " &
        // "during kernel optimization where the re-occupation of the energy level by the DMFT "&
        // "density kernel is taken into account !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "DMFT_CUTOFF_TAIL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Atomic energy threshold for double-precision calculations of the DMFT &
                   & Green-function at high frequencies. (Useful for GPU calculations.) !*"
  KW(i)%KWGRP = "DMFT"
  i=i+1

  KW(i)%LABEL = "IS_DIELECTRIC_EXCLUSIONS"
  KW(i)%TYP = "B:E"
  KW(i)%DSCRPT = "*! Dielectric exclusion regions !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "MULTIGRID_BC"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! 3 character string defining BCs for multigrid solver along &
                   &each lattice vector. 'O' for open, 'P' for periodic, 'Z' for zero &
                   &(i.e. open, but with the potential assumed to be zero on cell &
                   &boundaries). !*"
  KW(i)%KWGRP = "BC"
  i=i+1

  KW(i)%LABEL = "ION_ION_BC"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! 3 character string defining BCs for ion-ion interaction along &
                   &each lattice vector. 'O' for open, 'P' for periodic. !*"
  KW(i)%KWGRP = "BC"
  i=i+1

  KW(i)%LABEL = "PSPOT_BC"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! 3 character string defining BCs for local pseudopotential along &
                   &each lattice vector. 'O' for open, 'P' for periodic. !*"
  KW(i)%KWGRP = "BC"
  i=i+1

  KW(i)%LABEL = "SMEARED_ION_BC"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! 3 character string defining BCs for smeared ion representation along &
                   &each lattice vector. 'O' for open, 'P' for periodic. !*"
  KW(i)%KWGRP = "BC"
  i=i+1

  KW(i)%LABEL = "MG_GRANULARITY_POWER"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Power of 2 which gives multigrid granularity, i.e. granularity = 2**N &
                   &where N is MG_GRANULARITY_POWER. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_DEFCO_FD_ORDER"
  KW(i)%TYP = "I:B"
  KW(i)%DSCRPT = "*! Order of finite differences to use in the high-order defect correction &
                   &component of the multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_RES_REL"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Relative tolerance in norm of residual for defect correction procedure in &
                   &multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_RES_ABS"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Absolute tolerance in norm of residual for defect correction procedure in &
                   &multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_POT_REL"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Relative tolerance in norm of potential for defect correction procedure &
                   &in multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_POT_ABS"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Absolute tolerance in norm of potential for defect correction procedure &
                   &in multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_VCYC_REL"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Relative tolerance for norm of residual in multigrid V-cycle iterations. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_VCYC_ABS"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Absolute tolerance for norm of residual in multigrid V-cycle iterations. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_NEWTON_REL"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Relative tolerance for norm of residual in Newton method iterations in &
                   &multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_NEWTON_ABS"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Absolute tolerance for norm of residual in Newton method iterations in &
                   &multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_MAX_ITERS_DEFCO"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of iterations for the high-order defect correction &
                   &procedure in the multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_MAX_ITERS_VCYCLE"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of multigrid V-cycle iterations. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_MAX_ITERS_NEWTON"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of iterations for the Newton method in the &
                   &multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_USE_ERROR_DAMPING"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Should we use error damping in the high-order defect correction &
                   &procedure of the multigrid solver? !*"
  KW(i)%KWGRP = " MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_USE_FAS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Should FAS be used for non-linear PB equation &
                   &in the multigrid solver (instead of Newton method)? !*"
  KW(i)%KWGRP = " MULTIGRID"
  i=i+1

  KW(i)%LABEL = "FINITE_DIFFERENCE_ORDER"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Order of finite differences to use outside of high-order &
                   &defect correction, e.g. for computing the electric field in &
                   &open boundary conditions. !*"
  KW(i)%KWGRP = " FD"
  i=i+1

  KW(i)%LABEL = "EDFT_SPIN_FIX"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of iterations to fix the spin, negative is forever !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "POL_EMB_REPULSIVE_MM_POT_R0"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Pol-emb repulsive MM pot: R0 parameter !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL = "POL_EMB_POLSCAL"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Pol-emb: QM polarisability scaling factor !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  ! ja531 --> Fermi Operator Expansion

  KW(i)%LABEL = "CONTRACOHAM_RADMULT"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Sparsity pattern for Contra-Covariant &
       &Ham radius multiplier !*"
  i=i+1

  KW(i)%LABEL = "H2DENSKERN_SPARSITY"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Use the H^2 sparsity pattern for K in FOE !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "FOE_MU_TOL"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Tolerance for stopping in FOE &
       &chemical potential search. !*"
  i=i+1

  KW(i)%LABEL = "FOE_TEST_SPARSITY"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Test the quality of the H^2 sparsity pattern for K &
       &in FOE !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "SPECIES_BSUNFLD_PROJATOMS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Species projected atoms for spectral function unfolding &
       &calculation !*"
  KW(i)%KWGRP = "BS_UNFOLDING"
  i=i+1

  KW(i)%LABEL = "AUGBOX_PREF"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Preferred Augmentation box dimensions !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL  = "NEB_CI_DELAY"
  KW(i)%TYP    = "I:I"
  KW(i)%DSCRPT = "*! Number of NEB chain LBFGS steps before climging image &
       &kicks in.  Negative number for no climbing image (default). !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  ! lam81
  KW(i)%LABEL = "GEOM_PRECOND_TYPE"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! type of preconditioner &
       & NONE = LBFGS from CASTEP; &
       & ID = ID/LBFGS should be identical to NONE; &
       & EXP = EXP/LBFGS. !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_PRECOND_SCALE_CELL"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Scaling cell in variable cell optimisation. !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_PRECOND_EXP_C_STAB"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! stabilization constant of EXP preconditioner !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_PRECOND_EXP_A"
  KW(i)%TYP = "D:B"
  KW(i)%DSCRPT = "*! A value of EXP preconditioner !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_PRECOND_EXP_R_NN"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! nearest neighbor distance for EXP preconditioner !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_PRECOND_EXP_R_CUT"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! cutoff distance for EXP preconditioner !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_PRECOND_EXP_MU"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! mu value for EXP preconditioner !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_PRECOND_FF_C_STAB"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! stabilization constant of FF preconditioner !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "GEOM_PRECOND_FF_R_CUT"
  KW(i)%TYP = "P:B"
  KW(i)%DSCRPT = "*! cutoff distance for FF preconditioner !*"
  KW(i)%KWGRP = "GEOM"
  i=i+1

  KW(i)%LABEL = "QNTO_NBO_PROJ"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT =  "*! Enables projection to NBOs !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "QNTO_REF_DIR"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT =  "*! Sets the reference directory for NTO projection !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "QNTO_NUM_REF_STATES"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT =  "*! Sets the number of states to reference !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL = "QNTO_NUM_CORE_ATOMS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT =  "*! Sets the number of atoms for projection !*"
  KW(i)%KWGRP = "LRTDDFT"
  i=i+1

  KW(i)%LABEL ="LR_TDDFT_SPECTRUM_SMEAR"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Gaussian smearing half-width for LR-TDDFT spectrum !*"
  KW(i)%KWGRP ="LRTDDFT"
  i=i+1

  KW(i)%LABEL = "IS_DIELECTRIC_EXCLUSIONS_SMEAR"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Implicit solvent: smoothing on boundaries of &
       &dielectric exclusion regions defined in the &
       &IS_DIELECTRIC_EXCLUSIONS block.  This is a smearing distance. !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "MM_REP_PARAMS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Repulsive MM potential params of all MM species !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL  = "FOE_CHEBY_THRES"
  KW(i)%TYP    = "D:I"
  KW(i)%DSCRPT = "*! The maximum error threshold on the Chebyshev expansions &
       &in the FOE !*"
  i=i+1

  KW(i)%LABEL  = "FOE_AVOID_INVERSIONS"
  KW(i)%TYP    = "L:I"
  KW(i)%DSCRPT = "*! Avoid performing any inversions or using any inverses &
       &in the FOE !*"
  i=i+1

  KW(i)%LABEL  = "FOE_CHECK_ENTROPY"
  KW(i)%TYP    = "L:I"
  KW(i)%DSCRPT = "*! Validate the FOE entropy approximation &
       &against a simple  quadratic form !*"
  i=i+1

  KW(i)%LABEL = "TSSEARCH_ENERGY_TOL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Energy tolerance for TS search !*"
  KW(i)%KWGRP = "TS"
  i=i+1

  KW(i)%LABEL  = "NEB_PRINT_SUMMARY"
  KW(i)%TYP    = "L:I"
  KW(i)%DSCRPT = "*! Flag to print NEB pathway and convergence information &
       &to stdout !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "REACTANT_ROOTNAME"
  KW(i)%TYP    = "T:I"
  KW(i)%DSCRPT = "*! Specification of reactant rootname for energy &
       &calculation.  User must also include .tightbox_ngwf and .dkn files &
       &in this directory. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "PRODUCT_ROOTNAME"
  KW(i)%TYP    = "T:I"
  KW(i)%DSCRPT = "*! Specification of product rootname for energy &
       &calculation.  User must also include .tightbox_ngwf and .dkn files &
       &in this directory. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "REACTANT_ENERGY"
  KW(i)%TYP    = "P:I"
  KW(i)%DSCRPT = "*! Direct specification of reactant energy. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "PRODUCT_ENERGY"
  KW(i)%TYP    = "P:I"
  KW(i)%DSCRPT = "*! Direct specification of product energy. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "DMFT_WRITE"
  KW(i)%TYP    = "L:I"
  KW(i)%DSCRPT = "*! Write DMFT data (Green's functions, Hamiltonians, etc.) !*"
  KW(i)%KWGRP  = "DMFT"
  i=i+1

  KW(i)%LABEL  = "DMFT_READ"
  KW(i)%TYP    = "L:I"
  KW(i)%DSCRPT = "*! Read DMFT data (Green's functions, Hamiltonians, etc.) &
       &if they are present!*"
  KW(i)%KWGRP  = "DMFT"
  i=i+1

  KW(i)%LABEL  = "NEB_SPRING_CONSTANT"
  KW(i)%TYP    = "P:I"
  KW(i)%DSCRPT = "*! Spring constant for the NEB chain. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "NEB_READ_XYZ"
  KW(i)%TYP    = "L:E"
  KW(i)%DSCRPT = "*! Read XYZ file for each image. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "NEB_CONVERGE_ALL"
  KW(i)%TYP    = "L:E"
  KW(i)%DSCRPT = "*! Use energy and displacement convergence criteria &
                  &for NEB as well as forces. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "NEB_UPDATE_METHOD"
  KW(i)%TYP    = "T:I"
  KW(i)%DSCRPT = "*! Update method for NEB. Currently supported: &
                    &FIRE (default), GLBFGS. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "NEB_MAX_ITER"
  KW(i)%TYP    = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of NEB iterations. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL = "HFX_OUTPUT_DETAIL"
  KW(i)%TYP = "T:B"
  KW(i)%DSCRPT = "*! Level of output detail for HFx: BRIEF, NORMAL, VERBOSE, PROLIX !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "HFX_DEBUG"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Perform extra sanity checks for HFx? !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "SWRI_PRINT_EIGENVALUES"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Print debugging SW metric matrix eigenvalue info? !*"
  KW(i)%KWGRP = "SWRI"
  i=i+1

  KW(i)%LABEL = "WRITE_OVERlAP"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Save current Overlap matrix in a file !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "POL_EMB_REPULSIVE_MM_POT_VERBOSE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Pol-emb repulsive MM pot: verbose output? !*"
  KW(i)%KWGRP = "QMMM"
  i=i+1

  KW(i)%LABEL  = "USE_SPH_HARM_ROT"
  KW(i)%TYP    = "L:E"
  KW(i)%DSCRPT = "*! Initialize spherical harmonic rotation module. Not &
       &needed in normal use, as this should be done automatically. !*"
  KW(i)%KWGRP  = "SPH_HARM_ROT"
  i=i+1

  KW(i)%LABEL  = "NEB_GLBFGS_HISTORY_SIZE"
  KW(i)%TYP    = "I:E"
  KW(i)%DSCRPT = "*! History size of GLBFGS for NEB. !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL  = "NEB_CONTINUATION"
  KW(i)%TYP    = "L:I"
  KW(i)%DSCRPT = "*! Continue NEB run from .neb_cont files !*"
  KW(i)%KWGRP  = "TS"
  i=i+1

  KW(i)%LABEL = "SPECIES_NGWF_REGIONS"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Regions that each atom is allocated to !*"
  KW(i)%KWGRP = "GENERAL"
  i=i+1

  KW(i)%LABEL = "DO_FANDT"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Perform a freeze-and-thaw optimisation of the NGWFs !*"
  i=i+1

  KW(i)%LABEL = "USE_EMFT"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Do an embedding mean field theory calculation !*"
  i=i+1

  KW(i)%LABEL = "ACTIVE_XC_FUNCTIONAL"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Defines the xc functional used for the higher level "&
                     // "calculation within an embedding mean field theory calculation !*"
  i=i+1

  KW(i)%LABEL = "FREEZE_SWITCH_STEPS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! No. of CG steps to perform before switching F+T !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "FREEZE_ENVIR_NGWFS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Never optimise the environment NGWFs !*"
  i=i+1

  KW(i)%LABEL = "ACTIVE_REGION"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Defines which region is the active region used for the higher level "&
                     // "calculation within an embedding mean field theory calculation !*"
  i=i+1

  KW(i)%LABEL = "PARALLEL_SCHEME"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Types of parallel strategies that can be used in subsystem calculations !*"
  i=i+1

  KW(i)%LABEL = "USE_EMFT_FOLLOW"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Do an EMFT calculation after a regular NGWF optimisation!*"
  i=i+1

  KW(i)%LABEL = "USE_EMFT_LNV_ONLY"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Do an LNV only EMFT calculation !*"
  i=i+1

  KW(i)%LABEL = "BLOCK_ORTHOGONALISE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Orthogonalise environment NGWFs wrt active subsystem !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "READ_SUB_DENSKERN"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Read density kernel restart information from subsystem kernels. !*"
  KW(i)%KWGRP = "IO"
  i=i+1

  KW(i)%LABEL = "EMBED_DEBUG"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Verbose printing for embedding if debugging !*"
  KW(i)%KWGRP = "I"
  i=i+1

  KW(i)%LABEL = "PROJECT_EMBED"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Do an embedding calculation using the projector method !*"
  i=i+1

  KW(i)%LABEL = "EMFT_LNV_STEPS"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Number of LNV iterations during EMFT kernel optimisation !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "EXTERNAL_BC_FROM_CUBE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Read external potential for boundary conditions from cube file !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  ! aeaa - DDEC anisotropy terms
  KW(i)%LABEL = "DDEC_ANISO"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Calculates off center point charges !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_ANISO_MAX_DIS"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Sets the maximum distance from the atom center !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_ANISO_MAX_DIS_HALOGEN"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Sets the maximum distance from the atom center &
                       &for the halogens!*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_ANISO_ERROR_THRES"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Sets the threshold above which off center charges &
       &are added !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL = "DDEC_ANISO_ERROR_REDUCE"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Sets the improve in ESP needed before off center charges &
       &are added !*"
  KW(i)%KWGRP = "DDEC"
  i=i+1

  KW(i)%LABEL  = "SPECIES_SCISSOR"
  KW(i)%TYP    = "B:I"
  KW(i)%DSCRPT = "*! Apply energy shift to species hamiltonian eigenvalues !*"
  i=i+1

  KW(i)%LABEL = "PDOS_SUM_MAG"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Sum over magnetic quantum number in PDOS (true by default) !*"
  i=i+1

  KW(i)%LABEL = "PDOS_CONSTRUCT_BASIS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Compute PDOS using a SW fit to NGWFs !*"
  i=i+1

  KW(i)%LABEL = "PDOS_LOWDIN"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Compute PDOS by taking the Lowdin factorization of the SW overlap !*"
  i=i+1

  KW(i)%LABEL = "PDOS_LCAO_OPTIMIZE"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Compute PDOS by solving Kohn-Sham equations with an LCAO basis !*"
  i=i+1

  KW(i)%LABEL = "PDOS_ORTH_ATOM_BLOCKS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Orthogonalize the LCAOs on same centres in PDOS !*"
  i=i+1

  KW(i)%LABEL = "PDOS_OUTPUT_BASIS"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Write the pDOS SW basis to disk in tightbox_ngwf format !*"
  i=i+1

  KW(i)%LABEL = "PDOS_OUTPUT_SWOPT_KERNHAM"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Write the kernel / Hamiltonian in a SW optimized pDOS !*"
  i=i+1

  KW(i)%LABEL = "PDOS_PSEUDOATOMIC"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Use pseudoatomic functions as the PDOS AM resolved basis !*"
  i=i+1

  KW(i)%LABEL = "PDOS_MAX_N"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! The maximum number of Bessels to use in SW expansion in a pDOS calc !*"
  i=i+1

  KW(i)%LABEL = "IS_HC_STERIC_SMEARING"
  KW(i)%TYP = "P:E"
  KW(i)%DSCRPT = "*! Implicit solvent: smearing distance for smoothed hard-core potential !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_HC_STERIC_DENS_ISOVALUE"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Implicit solvent: n_0 parameter in electrolyte accessibility !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL  = "ELD_CALCULATE"
  KW(i)%TYP    = "L:B"
  KW(i)%DSCRPT = "*! Calculate electron localisation descriptors !*"
  KW(i)%KWGRP  = "ELD"
  i=i+1

  KW(i)%LABEL  = "ELD_FUNCTION"
  KW(i)%TYP    = "T:B"
  KW(i)%DSCRPT = "*! Choose which electron localisation descriptor to use &
       &during the properties calculation, either ELF or LOL !*"
  KW(i)%KWGRP  = "ELD"
  i=i+1

  KW(i)%LABEL  = "KE_DENSITY_CALCULATE"
  KW(i)%TYP    = "L:B"
  KW(i)%DSCRPT = "*! Calculate kinetic energy density !*"
  KW(i)%KWGRP  = "ELD"
  i=i+1

  KW(i)%LABEL = "CACHE_LIMIT_FOR_PRODS"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Max cache size for Aa-Dd NGWF products (in MiB) in HFx !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "MG_VCYC_SMOOTHER_ITER_PRE"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! V cycle smoother iterations pre-smoothing: passed to DL_MG !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_VCYC_SMOOTHER_ITER_POST"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! V cycle smoother iterations post-smoothing: passed to DL_MG !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_MAX_RES_RATIO"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! Residual ratio threshold for giving up on MG convergence: passed to DL_MG !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_MU_REL"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Relative tolerance for chemical potential in multigrid &
       &calculations in PBC with Boltzmann ions. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_MU_ABS"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Absolute tolerance for chemical potential in multigrid &
       &calculations in PBC with Boltzmann ions. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "IS_SOLVATION_PROPERTIES"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Produce scalarfields of solvation inputs and outputs in properties? !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "SPECIES_SOLVENT_RADIUS"
  KW(i)%TYP = "B:E"
  KW(i)%DSCRPT = "*! Implicit solvent: solvent radius around ions !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_PBE_ENERGY_TOLERANCE"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Absolute tolerance for difference between energy &
       &expressions in Boltzmann solvation !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "HFX_BESSEL_RAD_NPTSX"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! HFx: Number of points in Bessel radial interpolation !*"
  KW(i)%KWGRP = "HFx"
  i=i+1

  KW(i)%LABEL = "IS_PBE_NEUTRALISATION_SCHEME"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Neutralisation scheme for PBE solvation in PBCs !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "HFX_MEMORY_LIMIT"
  KW(i)%TYP = "I:E"
  KW(i)%DSCRPT = "*! Max cache size (per MPI rank) for all of HFx (in MiB) !*"
  KW(i)%KWGRP = "HFX"
  i=i+1

  KW(i)%LABEL = "HFX_MEMORY_WEIGHTS"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Three weights for HFx memory use (SWOP, EXPA, PROD). !*"
  i=i+1

  KW(i)%LABEL = "IS_RESTART_VAC_FROM_VAC"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Decides whether the vacuum calculation in an IS autosolvation &
       &calculation should be restarted from vacuum_* files or not !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_EMFT_CAVITY"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Decides whether the IS cavity is determined using the &
       &EMFT-optimised kernel or the normal kernel !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "EDFT_TRIAL_STEP"
  KW(i)%TYP = "D:E"
  KW(i)%DSCRPT = "*! User input mixing parameter - replaces default line search. !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_UPDATE_SCHEME"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! Define the update scheme used in the inner loop of EDFT. !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_HAM_DIIS_SIZE"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Max number of Hamiltonians saved during EDFT Pulay DIIS. !*"
  KW(i)%KWGRP = "CONV"
  i=i+1

  KW(i)%LABEL = "EDFT_GRAND_CANONICAL"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! use grand canonical ensemble in ensemble DFT !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_REFERENCE_POTENTIAL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! electrochemical potential of reference electrode  !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "EDFT_ELECTRODE_POTENTIAL"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! electrode potential used to set the Fermi level !*"
  KW(i)%KWGRP = "EDFT"
  i=i+1

  KW(i)%LABEL = "IS_SOFT_SPHERE_RADII"
  KW(i)%TYP = "B:I"
  KW(i)%DSCRPT = "*! Implicit solvent: Block of Van der Waals radii used to &
       &define the cavity size of elements in the soft sphere continuum model. !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SOFT_SPHERE_DELTA"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Implicit solvent: Value delta used to define the smoothing function of the soft spheres. !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_SOFT_SPHERE_SCALE"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Implicit solvent: Value used to parametrise the library of solvation radii. !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "IS_BC_ALLOW_FRAC_CHARGE"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Implicit solvent: Don't check for total charge being an integer in MG BCs. !*"
  KW(i)%KWGRP = "SOLVATION"
  i=i+1

  KW(i)%LABEL = "DFTB"
  KW(i)%TYP = "L:B"
  KW(i)%DSCRPT = "*! Perform a DFTB calculation? !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_METHOD"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Choice of a method for DFTB !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_METHOD_PARAM_FILE"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Name of the DFTB method parameter definition file !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_COMMON_PARAM_FILE"
  KW(i)%TYP = "T:I"
  KW(i)%DSCRPT = "*! Name of the DFTB common parameter definition file !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_BC"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! 3 character string defining BCs IN DFTB calculations along &
                   &each lattice vector. 'O' for open, 'P' for periodic. !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_COORD_CUTOFF"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! distance cutoff for calculation of coordination number. !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_REP_CUTOFF"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! distance cutoff for calculation of DFTB repulsion term. !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_SRB_CUTOFF"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! distance cutoff for calculation of DFTB SRB term. !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_OVERLAP_CUTOFF"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! distance cutoff for calculation of overlap matrix. !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_OVERLAP_ANALYTICAL"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! In DFTB, is the overlap matrix to be calculated analytically? !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_EWALD_REPLICATE_XTB"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! In GFN0 DFTB, replicate Ewald summation in xTB? !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "DFTB_EWALD_PARAMETER"
  KW(i)%TYP = "P:I"
  KW(i)%DSCRPT = "*! Convergence parameter for the Ewald summation in DFTB. !*"
  KW(i)%KWGRP = "DFTB"
  i=i+1

  KW(i)%LABEL = "MG_USE_CG"
  KW(i)%TYP = "L:I"
  KW(i)%DSCRPT = "*! Implicit solvent: Use conjugate gradients in DL_MG. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_CG_RES_REL"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Relative tolerance in norm of residual for defect correction procedure in &
                   &multigrid solver for CG. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_TOL_CG_RES_ABS"
  KW(i)%TYP = "D:I"
  KW(i)%DSCRPT = "*! Absolute tolerance in norm of residual for defect correction procedure in &
                   &multigrid solver for CG. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "MG_MAX_ITERS_CG"
  KW(i)%TYP = "I:I"
  KW(i)%DSCRPT = "*! Maximum number of iterations for conjugate gradients &
                   &in the multigrid solver. !*"
  KW(i)%KWGRP = "MULTIGRID"
  i=i+1

  KW(i)%LABEL = "PERMIT_UNUSUAL_NGWF_COUNT"
  KW(i)%TYP = "L:E"
  KW(i)%DSCRPT = "*! Allows continuing the calc with suspect number of NGWFs. !*"
  i=i+1

  KW(i)%LABEL = "HFX_BC"
  KW(i)%TYP = "T:E"
  KW(i)%DSCRPT = "*! 3 character string defining BCs for HFx along &
                   &each lattice vector. 'O' for open, 'P' for periodic. !*"
  KW(i)%KWGRP = "VDW"
  i=i+1


  ! **********************************************************************
  ! jd: Add new keywords here, using this template:

!  KW(i)%LABEL = "..."
!  KW(i)%TYP = "..."
!  KW(i)%DSCRPT = "*! ... !*"
!  KW(i)%KWGRP = "..."
!  i=i+1
  ! **********************************************************************

  ! jd: This should follow the last keyword
  numkw = i-1

 end subroutine esdf_key_init

end module esdf_key
