! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                    G E O M E T R Y                                          !
!=============================================================================!
!                                                                             !
! $Id: geometry_optimiser_mod.F90,v 1.37 2009/09/23 16:56:54 cks22 Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! This module performs a geometry optimisation on the specified model.        !
! It may be used to modify ionic positions and/or cell parameters.            !
!                                                                             !
!-----------------------------------------------------------------------------!
! Written from  "Module Specification for New Plane Wave Code",  M. Segall,   !
!    P. Lindan, M. Probert, C. Pickard, P. Hasnip, S. Clark and M. Payne      !
!                           Copyright 1999/2000                               !
!-----------------------------------------------------------------------------!
! Written by Matt Probert, v1.0, 11/05/2002                                   !
!-----------------------------------------------------------------------------!
! Modified for ONETEP by Arash A Mostofi                                      !
! Cleaned up by Alvaro Ruiz Serrano in September 2009.                        !
! lBFGS implemented by Jolyon Aarons in 2011.                                 !
! Subsequent modifications and contributions by Nicholas Hine, Fabiano        !
! Corsetti, Jose M Escartin, Andrea Greco, and Loukas Kollias.                !
!-----------------------------------------------------------------------------!
!
! Revision 1.38  2022/08/11 12:59:20  lk
! lk_2022_08_11: Introduced pub_geom_output_detail to control the level of
! of detail of the output of the geometry optimization section.
! Subsequently, replaced pub_output_detail with pub_geom_output_detail
! throughout. Comments on the affected parts of the code start with "lk:".
!
! kkbd,2018_01_30: Adjusted to allow the minimization of a NEB chain-of-states
! that contains several beads with springs connecting neighbors.
!
! ars_2012_11_08: added the ability to read/write the Hamiltonian matrix during
! geometry optimization. It applies to the .continuation and .bfgs_backup files.
! So far, only active if we are doing ensemble-DFT calculations.
!
! $Log: geometry_optimiser_mod.F90,v $
! Revision 1.37  2009/09/23 16:56:54  cks22
!
! cks_2009_09_23:
! -> Alvaro's clean up of the geometry_optimiser
! -> Fix geometry optimisation continuations
! -> Add XLYP dispersion parameters
! -> Faster calculation of spherical harmonics
! -> Various minor fixes and cleanups
!
! Revision 1.36  2009/08/14 10:39:03  ndmh3
!
! ndmh_2009_08_14: Addition of SPAM3, function_basis and new function sum/integrals
! routines. Cleanup of whitespace, capitalisation and indentation. Bugfixes to
! TDDFT and DFT+U. Minor format changes, rearrangement of atom counting information
! to show projector count. New approach to nonlocal potential matrices in
! norm_conserv_pseudo, to go with SPAM3. Improvements to Peano Space Filling Curve
! distribution in parallel_strategy_mod. Rearrangement of ngwf_gradient_mod to
! group routines relating to construction of coefficients for function sums.
!
! Revision 1.35  2009/06/23 09:28:13  cks22
!
! cks_2009_06_23: Bug fixes for spin polarisation, an additional kernel fix approach,
! BLYP parameters and periodic boundary conditions for the VdW module, clean-up of
! various subroutines.
!
! Revision 1.34  2009/05/11 14:32:40  ndmh3
!
! ndmh_2009_05_11: Removed old-style usage of output_detail, replaced with
! pub_output_detail throughout. Fixed memory allocation problem when using
! NNHO's with dense_matrices: T. Tidied up some formatting in properties_mod.
! Fixed a badly formatted error message in energy_and_force_mod.
! qoh: fixed bug in spherical_wave_mod introduced in 2.3.13
!
! Revision 1.33  2009/04/22 15:11:54  cks22
!
! cks_2009_04_22: several minor bug fixes and the inclusion of "classical" point
! charges to the external potential.
!
! Revision 1.32  2009/02/05 03:20:59  cks22
!
! cks_2009_02_05: Hartree-Fock exchange (slow algorithm, using the FFT box equal
! to the simulation cell). Also a bug fix by Nick Hine in
! restart_ngwfs_tightbox_input to follow the F95 standard and fix wrong reading
! of NGWFs with pg compilers.
!
! Revision 1.19  2009/02/02 17:53:02  cks
!
! cks_2009_02_02: Version 2.2.16
!
! Revision 1.31  2009/01/26 16:49:33  ndmh3
! ndmh_2009_01_26: fixed step sizes in geometry_optimiser_mod by fixing bug
! in inverse Hessian initialisation involving units of mass (1 m_e /= 1 amu!)
!
! Revision 1.5  2008/12/01 14:08:43  vmilman
! Merge with the updated ODG version 2.2.14
!
! Revision 1.29  2008/10/28 11:09:56  cks22
!
! cks_2008_10_28: Trivial changes to
! 1) Eliminate "cell" as an argument and replace it with "cell".
! 2) Fix a recent bug in the initialisation of tight boxes which made the code
! crash when the number of ppds was not the same along a1 and a2.
! 3) Make the code run without crashing when the error bounds chacking is used
! with pgf90.
!
! Revision 1.28  2008/09/22 18:43:45  ndmh3
!
! ndmh_22_09_2008: addendums to previous checkin, cache-optimisation of
! norm_conserv_pseudo_mod structure factor code, various cleanups and minor
! fixes, adaptable output during penalty and lnv for very large systems (no
! change for small systems).
!
! Revision 1.27  2008/08/27 15:08:48  cks22
!
! cks_2008_08_27: Restored working of empirical dispersion in the code.
!Also added the input flag geom_reuse_dk_ngwfs (true by default) which when
!set to false it prevents reading the density kernel and ngwfs from the
!previous geometry step as this can sometimes cause (mysterious) ngwf
!convergence problems.
!
! Revision 1.26  2008/08/13 15:34:14  ndmh3
!
! ndmh_2008_08_13: file omitted from previous check-in
!
! Revision 1.25  2008/07/31 17:48:24  cks22
!
! cks_2008_07_31: (1) Trivial changes in spherical_mod and (2) incorporation of
! empirical dispersion forces and energy in order to be able to use dispersion
! contributions in geometry optimisation.
!
! Revision 1.12  2008/07/31 17:41:16  cks
!
! cks_2008_07_31: VdW forces and geopt.
!
! Revision 1.24  2008/07/02 16:39:22  pdh1001
! Minor changes to satisfy gfortran.
!
! Revision 1.23  2008/03/30 14:56:57  cks22
!
! cks_2008_03_30: Memory deallocation bug fixes and write_xyz now writes xyz
! history during optimisation.
!
! Revision 1.22  2008/02/15 10:38:17  pdh1001
! Several minor bug fixes.
! Added NGWF analysis functionality.
! Added change from Victor to population analysis output.
!
! Revision 1.21  2007/08/19 08:54:54  cks22
!
! cks_2007_08_19: Extension to use simulation cells smaller than the 6 NGWF
! radii rule, along
! one or more lattice vector directions. The minimum simulation cell is now
!defined by the twice
! the maximum NGWF radius.
!
! Revision 1.3  2007/08/17 17:17:14  cks
!
! cks_2007_08_17: Merged flex and current versions.
!
! Revision 1.20  2007/02/17 00:08:40  pdh1001
! Added capability to write sparse matrices in SPAM3 format to a binary file
! in a versioned, extensible format.
!
! Revision 1.19  2007/01/31 14:10:45  pdh1001
! Victor's changes to TS search and delocalised internals, plus Windows fix.
!
! Revision 1.18  2006/12/11 11:08:38  pdh1001
! Minor changes to avoid continuation bug: number of spins now included in
! .continuation file and SPAM%matrix initialised to zero after allocation.
!
! Revision 1.17  2006/11/16 13:23:55  pdh1001
! Bug fix for spin polarisation in backup, restore and continuation file routine
!
! Revision 1.16  2006/09/08 16:32:46  pdh1001
! FFTw version 3 interface added (compile with -DFFTW3 instead of -DFFTW)
! Energy components printed in BRIEF mode
! Other superficial changes
!
! Revision 1.15  2006/05/15 08:20:35  cks22
!
! cks_2006_05_15: added properties module, brief output level, Victor's
! latest TS changes, ngwf initialisation up to element 53.
!
! Revision 1.14  2006/04/21 15:42:09  cks22
!
! cks_2004_04_21: Added Victor's transition state search files.
! Also added ngwf_halo, a new experimentlal feature which can take into account
! small but non-zero matrix elements.
!
! Revision 1.13  2006/02/10 16:45:07  pdh1001
! Merged with source from Victor Milman (mainly to incorporate delocalized
! internals) and some minor bug fixes.
!
! Revision 1.12  2005/09/30 13:22:06  aam24
! Small change to supress inverse Hessian printout.
!
! Revision 1.11  2005/09/27 10:23:45  aam24
! New logical keywords "geom_print_inv_hessian" and "verbose_ewald_forces" for
! controlling whether to write to standard output the inverse hessian
! during a geometry optimisation run and debugging information on the
! Ewald forces, respectively.
!
! Revision 1.10  2005/09/13 15:17:54  aam24
! Parallel bug fix in subroutine geom_opt_continuation_read for geometry
! optimiser continuation facility.
!
! Revision 1.9  2005/08/31 14:39:29  aam24
! Small change to make pub_task=GeometryOptimization work with old_input=TRUE.
!
! Revision 1.8  2005/08/04 14:10:32  aam24
! Output absolute cartesian (rather than fractional) coordinates in geometry
! optimisation.
!
! Revision 1.7  2005/08/04 11:02:58  aam24
! First implementation of ionic constraints.
!
! Revision 1.6  2005/07/26 10:24:18  aam24
! Geometry optimiser now able to cope with fractional positions outside
! the interval [0,1]. Minor changes to norm_conserve_pseudo, electronic,
! energy_and_force and onetep.F90.
!
! Revision 1.5  2005/06/08 10:16:10  aam24
! Minor structural changes to geometry_optimiser and onetep
! modules. Average residual force on atoms is now explicitly
! removed. Addition of new "pub_task" and subroutines for testing
! forces. Addition of a function in energy_and_force module to count the
! number of psinc functions enclosed in NGWF spheres (you need to
! uncomment the relevant call).
!
! Revision 1.4  2005/02/28 17:10:26  aam24
! Fixing changes in onetep.F90 that were accidentally overwritten in the
! last revision. Very minor change in geometry_optimiser module.
!
! Revision 1.3  2005/02/23 18:32:52  aam24
! Additional geometry optimisation parameters in rundat and esdf_key,
! including logical keyword "geom_continuation" which allows one to
! resume a previous geometry optimisation. Corresponding changes to
! geometry_optimiser module and onetep.F90. Addition of two "physical"
! types in esdf_mod.F90 to deal with units of "ha/bohr**3" for
! rundat parameter "geom_modulus_est" (estimate of the bulk modulus) and
! "ha/bohr" for "geom_force_tol" (convergence tolerance for the forces).
!
! Revision 1.2  2005/01/27 17:49:22  aam24
! Small changes to geometry_optimiser module to get around an apparent
! Intel compiler bug.
!
! Revision 1.1  2004/12/02 11:59:09  aam24
! Addition of geometry optimisation module with corresponding minor
! changes to some other modules. Geometry optimisation is invoked with
! the 'pub_task' rundat parameter. Addition of two periodic table arrays
! to constants_mod.F90.
!
! Revision 1.60  2003/11/20 17:20:22  mijp2
! catch case of non-default fixed_npw with non-varying cell [bug 03324vymq01]
!
! Revision 1.59  2003/11/03  17:29:24  mijp2
! made frac_symm_ops allocatable to fix bug 03304vymq05
!
! Revision 1.58  2003/10/31  14:58:10  mijp2
! added constraints bug-fix in deloc-internals by Victor
!
! Revision 1.57  2003/09/30  09:30:52  mijp2
! minor forchk update
!
! Revision 1.56  2003/08/06 13:53:27  mijp2
! fixed latest restart bug
!
! Revision 1.55  2003/08/05  14:55:41  mijp2
! fixed bug 03148vymq02
!
! Revision 1.54  2003/07/28  17:00:33  mijp2
! added fixed_npw
! forchk neatening up
!
! Revision 1.53  2003/07/16  16:55:37  mijp2
! fixed variable cell, variable #pw restart bug
! fixed nb=0 bug
!
! Revision 1.52  2003/05/24  18:13:03  mijp2
! fixed bug in backup/restore with const Ecut cell changes
!
! Revision 1.51  2003/05/22  09:57:14  mijp2
! better purification of noise in frac_symm_ops
! changed convergence path in variable cell with const_pw=.false.
!
! Revision 1.50  2003/05/12  21:01:29  mijp2
! corrected intent() errors
!
! Revision 1.49  2003/05/07  15:02:10  mijp2
! added missing initialisation for last_write_time
!
! Revision 1.48  2003/04/30  22:11:40  mijp2
! Added saves to module variables.
! Added delocalised internals minimiser (Victor)
! Changed to fixed Ecutoff not fixed #plane waves
!
! Revision 1.47  2003/04/11  16:52:41  mijp2
! fixed run_time restart bug
!
! Revision 1.46  2003/03/20  11:35:04  mijp1
! FORCHK caught a type error in unsymm_to_symm
!
! Revision 1.44  2003/02/07  09:37:02  mijp2
! fixed parallel bug in timed backups
!
! Revision 1.43  2003/01/29  21:29:00  mijp2
! added timed backups using backup_interval
!
! Revision 1.42  2002/11/15  11:06:21  mijp2
! renamed %nbands -> %nbands_max and %nbands_occ -> %nbands
!
! Revision 1.41  2002/08/23 13:22:52  mijp2
! added off-diagonal external pressure (ie shear) capability
!
! Revision 1.40  2002/07/29  15:12:42  mijp2
! tweaked previous RCS log message as it causes pgf90 problems!
!
! Revision 1.39  2002/07/29  11:04:48  mijp2
! added frequency_label,
! deleted off-diagonal external_pressures (bug 02128s1wq01),
! clarified cell_constraints warning,
! added some extra checks on lambda to ensure minimum cell vector lengths,
! fixed bug in hs_vec calculation in geom_update_inv_Hessian_params,
! minor neatening of output,
! explicit dE/ion in convergence table,
! only do bisection if trying to guarantee downhill only steps.
! Phew!
!
! Revision 1.38  2002/07/16  12:42:02  mijp2
! extended geom_symm_delta to operate on ionic displacements as well as cell
! strains
!
! Revision 1.37  2002/07/12  16:24:16  mijp2
! fixed bug 02159vymq01 - BFGS sometimes violates cell constraints
! by the addition of new private routine geom_symm_delta which
! ensures that the strain obeys symmetry before the cell is updated.
!
! Revision 1.36  2002/07/09  14:20:48  mijp2
! fixed bug 02168vymq01 (final output forces showed zeroed constraints)
! changed output to only print cell info if cell changed and v.v. ions
!
! Revision 1.35  2002/05/16  16:08:27  mijp2
! fixed subtle restart bug [02136vymq01]
!
! Revision 1.34  2002/05/14  16:08:09  mijp2
! corrected dE bug in geom_converged [bug 02134vymq01]
!
! Revision 1.33  2002/05/13  12:30:33  mijp2
! tweaked uphill / quad heuristics to fix bug
!                    02119g1eq03
!
! Revision 1.32  2002/05/12  00:09:42  mijp2
! numerous minor bug fixes and tweaks:
! improved output messages to reduce user confusion [bugs 02120s1wq01 and
! 02119g1eq03]
! revised geom_check_lambda [bug 02114s1wq01]
! updated force output with constraints [bug 02114vymq02]
! deleted fake_forces as no longer needed
! improved backtracking if need to go uphill
! tweaked quad minimisation heuristic
! no longer rationalise ionic coords in .castep output
! sanitised final messages
! removed windows on |F|max, |dR|max and Smax convergence
!
! Revision 1.31  2002/04/24  17:09:02  mijp2
! bug fix in random noise
!
! Revision 1.30  2002/04/19  12:26:40  mijp2
! updated FBSC routine to be in sync with castep.f90 v1.30
!
! Revision 1.29  2002/04/18  16:30:56  mijp2
! deleted cell_update_symmetry (now in cell module) and fixed up call-logic
! beefed-up constraints & symmetry checks for cell and ions
! minor FBSC bug fixed
! catch instability in inv_Hessian for long runs and correct it
!
! Revision 1.28  2002/04/16  16:51:01  mijp2
! random noise bug fix
! changed order of basis_initialise and model_cell_changed to sort out beta-phi
! bug
! correctly updating symmetry if cell vectors change
! tweaked termination conditions if try to bisect and in the force/stress noise
! floor
! catch system that is fully constrained by symmetry
! correctly update the |dR|max field if backtrack a step
!
! Revision 1.27  2002/04/03  13:05:02  mijp2
! undid fudge factor
! tweaked error messages
!
! Revision 1.26  2002/03/04  17:26:37  mijp2
! changed BFGS backup file to scratch type so not visible to user
! made FBSC more efficient
!
! Revision 1.25  2002/02/18  16:20:37  mijp2
! added geom_update_params to catch user changes to modulus/frequency
! estimates and update inverse Hessian appropriately
! tweaked o/p some more
!
! Revision 1.24  2002/02/15  01:25:11  mijp2
! added nullify methods
! added DM stuff
! changed output formats
! neatened output and verbosity
! neatened up convergence output
!
! Revision 1.23  2002/02/07  17:41:41  mijp2
! fixed format error
!
! Revision 1.22  2002/02/07  17:21:49  mijp2
! allowed for on-the-fly unit changes in output
!
! Revision 1.21  2002/02/07  12.38.53  mijp2
! added check for no relaxation required
!
! Revision 1.20  2002/02/04  15:06:51  mijp2
! changed LAPACK usage for portability between v2 and v3
!
! Revision 1.19  2002/01/30  19:00:54  mijp2
! changed FBSC usage as problems with doing it mid-minimization
!
! Revision 1.18  2002/01/30  17:20:03  mijp2
! changed force & stress output
! added fermi energy to fbsc
! changed wave functions to S versions
!
! Revision 1.17  2002/01/17  17:40:03  mijp2
! fixed transposed cell vectors, and added efermi
!
! Revision 1.16  2001/12/19 10:54:43  mijp2
! fixed <-- E formatting in trajectory file
!
! Revision 1.15  2001/12/11  15:50:32  mijp2
! rewrote perfectly good F90 in geom_get_forces to workaround pgf90 compiler
! bug with -fast optimisation
!
! Revision 1.14  2001/11/28  09:42:57  mijp2
! fixed weird SGI bug with anint function
!
! Revision 1.13  2001/11/24  14:33:39  mijp2
! added file_maxpath
! changed .geom file format - now includes a header, cell vectors even for
! fixed cell calculation, non-rationalized positions, and everything in AU
! tweaks to output
! restart bug fixed
!
! Revision 1.12  2001/11/08  23:56:56  mijp2
! re-introduced finite basis correction calculation as in main
! deleted geom_ucase
! only o/p modulus and frequency if different to initial values
!
! Revision 1.11  2001/10/31  22:55:14  mijp2
! commented out one offending line missed in last checkin
!
! Revision 1.10  2001/10/31  22:36:42  mijp2
! removed finite basis set correction calculation - now in main
!
! Revision 1.9  2001/10/12  22:27:19  mijp2
! updated format statement
!
! Revision 1.8  2001/09/29  00:00:56  mijp2
! parallel IO bugfix, trajectory output bugfix
! improved inv_Hess reinitialise
! improved check_lambda
! improved IO
!
! Revision 1.7  2001/07/23  17:02:33  mijp2
! improved output
! fixed bug in |dR| calculation
!
! Revision 1.6  2001/07/19  10:50:51  mijp1
! added bisection search for line minimization failure situation
!  replaced quadratic lambda calculation algorithm
! tidied up logic for backtracking/further minimization
! implemented opt_strategy for fast backtracking
!
! Revision 1.5  2001/07/15  00:41:56  mijp2
! minor tweaks to logic and backup/backtracking
! some improvements to output formatting
!
! Revision 1.4  2001/07/03  02:35:42  mijp2
! finite basis corrections
!
! Revision 1.2  2001/06/27  22:06:05  mijp2
! added trajectory file output
!
! Revision 1.1  2001/06/15  14:31:18  mijp2
! Initial revision
!
!                                                                             !
!=============================================================================!

module geometry_optimiser

  use constants, only: DP
  use rundat, only: pub_debug_on_root
  use simulation_cell, only: castep_cell
#ifdef HDF5
  use hdf5
#endif

  implicit none                                 !Impose strong typing

  private                                       !Everything is private ...

  public :: geometry_optimise

  ! Cell data in CASTEP format
  type(castep_cell)  :: current_cell

  ! Fake constants
  integer, parameter :: file_maxpath = 256

  ! Constants and parameters
  real(kind=DP), parameter :: geom_stress_tol      = 0.1_dp
#ifdef PGI_BUG
  real(kind=DP), dimension(6) :: external_pressure = 0.0_dp
#else
  !kaw: This causes a compiler error for pgi v12.6-13.1.
  real(kind=DP), dimension(6), parameter :: external_pressure = 0.0_dp
#endif
  logical, save :: constrain_ions
  logical, save :: fix_all_ions
  integer, save :: num_ionic_constraints
  integer, save :: ndim  !the number of dimensions in the search space
  !so for the augmented-Hessian BFGS method used here,
  !we have ndim=9+3*num_ions as we seek to optimize the
  !cell (9 vars) and the ions (3*num_ions d.of.f) simultaneously.

  integer, save :: ndof  !the number of degrees of freedom in the problem
  ! - may be less than ndim due to symmetry and/or constraints

  integer, parameter :: shuffle_forwards  =  1
  integer, parameter :: shuffle_none      =  0
  integer, parameter :: shuffle_backwards = -1
  integer, parameter :: toggle_on         =  1
  integer, parameter :: toggle_off        = -1

#ifdef PGI_BUG
  integer, dimension(3,3) :: unit_matrix
#else
  integer, parameter, dimension(1:3,1:3) :: &
       unit_matrix=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
#endif

  integer, save             :: reset_iter, previous_reset_iter
  integer, save             :: iwarn             !global count of warning messages s.t. can remind user at end

  real (kind=dp), save      :: old_frequency_est, old_modulus_est
  real (kind=dp), save      :: new_frequency_est, new_modulus_est

  !Set up the pretty unit labels
  ! lk: changed label of energy units from Ha to Eh.
  character(len=*), parameter :: length_label    = 'Bohr'
  character(len=*), parameter :: energy_label    = 'Eh'
  character(len=*), parameter :: force_label     = 'Eh/Bohr'
  character(len=*), parameter :: pressure_label  = 'Eh/Bohr**3'
  character(len=*), parameter :: frequency_label = '1/aut'
  character(len=*), parameter :: force_constant_label = 'Eh/Bohr**2' ! lam81

  character(len=file_maxpath) :: cont_file

  ! ja531 --> LBFGS

  ! A few intermediate values.
  real(kind=dp), save,dimension(:,:), allocatable     :: YY_mat
  real(kind=dp), save,dimension(:),   allocatable     :: D_mat   ! The diagonal of an m x m matrix in which the ith diagonal
  ! entry is the dot product of the ith correction pair.
  real(kind=dp), save,dimension(:,:), allocatable     :: R_mat   ! m x m matrix L containing the matrix product S'*Y.
  real(kind=dp), save,dimension(:,:), allocatable     :: lbfgs_init_inv_hess
  real(kind=dp), save,dimension(:,:), allocatable     :: lbfgs_init_hess
  real(kind=dp),      dimension(:),   allocatable     :: tmpn

  real(kind=dp), save,dimension(:,:), allocatable     :: hessian !For BFGS

  integer, save                                    :: lbfgs_block_length=0

  logical, save                                    :: do_lbfgs, do_bfgs, do_tpsd
  logical, save                                    :: restart=.false.
  !  character(LEN=5),save                            :: method_string
  real(kind=dp), dimension(:,:), allocatable    :: tmpmm
  real(kind=dp), dimension(:), allocatable      :: tmp2m
  ! ja531 --> \LBFGS

  integer,save                                     :: iteration

  integer                                          :: hess_init_method
  integer, parameter                               :: hess_init_pfrommer=1
  integer, parameter                               :: hess_init_householder=2
  integer, parameter                               :: hess_init_identity=3
  integer, parameter                               :: hess_init_scaled=4
  integer, parameter                               :: hess_init_cellion=5

  logical                                          :: debug_square_order
  logical                                          :: bulk_freq_progress
  integer                                          :: debug_reinit
  logical                                          :: do_reinit

  real (kind=dp)                                   :: init_scalar(2) !scalar factor for scaled identity inv_hessian initialisation
  real (kind=dp)                                   :: debug_init_scale !factor to scale the above factor by!

  ! LBFGS test :: update the initialisation when a pair is removed from the representation
  logical                                          :: update_init_on_remove
  real(kind=dp), dimension(:), allocatable         :: old_x_vec_remove
  real(kind=dp), dimension(:), allocatable         :: old_f_vec_remove

  character(len=5)                                 :: geom_string
! ja531 -> If this is a continuation, then this is used for holding the read geom_string, for comparison
  character(len=5)                                 :: geom_continuation_string

  ! kkbd: NEB stuff
  integer, save :: prev_ci = -1

  ! lam81
  ! Precond type
  type preconditioner
      logical                                    :: initialized = .false.
      ! general
      logical                                    :: scale_cell
      real(kind=dp), dimension(:,:), allocatable :: P
      real(kind=dp), dimension(:,:), allocatable :: rmatrix
      real(kind=dp), dimension(:,:), allocatable :: dmatrix
      character(len=5), dimension(:), allocatable :: ion_constraint_type
      ! ID
      ! EXP
      real(kind=dp)                              :: exp_c_stab
      real(kind=dp)                              :: exp_A
      real(kind=dp)                              :: exp_mu
      real(kind=dp)                              :: exp_r_NN
      real(kind=dp)                              :: exp_r_cut
      ! FF
      real(kind=dp)                              :: ff_c_stab
      real(kind=dp)                              :: ff_r_cut
      integer, dimension(:), allocatable         :: ff_row
      real(kind=dp), dimension(:,:), allocatable :: ff_alpha
      real(kind=dp), dimension(:,:), allocatable :: ff_r_ref
      real(kind=dp)                              :: ff_kbond
      real(kind=dp)                              :: ff_kangle
      real(kind=dp)                              :: ff_kdihedral
  end type preconditioner

  type(preconditioner), save                     :: precond

contains

  subroutine geometry_optimise(total_energy,total_forces,mdl, &
       output_file,path,is_cmplx)

       !=========================================================================!
       ! Find ground state ionic positions and/or cell vectors for given model   !
       ! Also sets the private flags constrain_ions, symm_ions, etc.             !
       !-------------------------------------------------------------------------!
       ! Arguments:                                                              !
       !   cmdl, intent=inout, the model to be optimised                         !
       !   output_file, intent=in, the filename to which results will be written.!
       !-------------------------------------------------------------------------!
       ! Parent module variables used:                                           !
       !   model_file is set to checkpoint parameter                             !
       !   constrain_ions, symm_ions, constrain_cell, and symm_cell are set here !
       !   for all minimizers.                                                   !
       !   ndof and ndim are set here for all minimizers                         !
       !   pretty printing labels are set for all minimizers                     !
       !   on_root flag is set as appropriate for output on parallel machines    !
       !-------------------------------------------------------------------------!
       ! Modules used:                                                           !
       !   md for damped-md optimiser                                            !
       !   model for type definitions                                            !
       !   parameters to define the job etc                                      !
       !   cell to setup constraints flags                                       !
       !   io for error messages                                                 !
       !   comms to set on_root flag                                             !
       !-------------------------------------------------------------------------!
       ! Key Internal Variables:                                                 !
       !-------------------------------------------------------------------------!
       ! Necessary conditions:                                                   !
       !   cmdl has already been properly read and initialised                   !
       !-------------------------------------------------------------------------!
       ! Written by Matt Probert, v1.0, 11/05/2002                               !
       !-------------------------------------------------------------------------!
       ! Modified for ONETEP by Arash A Mostofi                                  !
       ! Modified to deal with embedding structures by Joseph Prentice, May 2018 !
       !=========================================================================!

       use comms, only: pub_on_root, comms_barrier
       use constants, only: dp, stdout
       use geometry, only: operator(.DOT.)
       use energy_and_force, only: energy_and_force_calculate
       use forces, only: forces_apply_constraints
       use model_type, only: MODEL
       use parallel_strategy, only: PARAL_INFO ! jcap
       use rundat, only: pub_write_denskern, pub_write_tightbox_ngwfs, pub_read_denskern, &
            pub_read_tightbox_ngwfs, pub_geom_reuse_dk_ngwfs, pub_maxit_pen, &
            pub_geom_method, pub_rootname, pub_geom_continuation, &
            pub_write_positions, pub_maxit_ngwf_cg, pub_edft, pub_write_hamiltonian, &
            pub_read_hamiltonian, pub_geom_lbfgs, &
            pub_geom_lbfgs_block_length, &
            pub_ngwf_regions_ngroups ! jcap
       use simulation_cell, only: castep_model, &
            copy_to_castep_cell, castep_cell_copy, &
            castep_model_alloc, castep_cell_2_elements, castep_model_dealloc, &
            castep_cell_dealloc
       use utils, only: utils_abort, utils_assert, utils_alloc_check
       use neb, only: neb_path, neb_energy_and_force

       implicit none

       character(len=*), intent(in) :: output_file
       type(MODEL), intent(inout)   :: mdl
       real(kind=DP), intent(inout) :: total_energy
       real(kind=DP), intent(inout) :: total_forces(1:3,mdl%nat)
       ! kkbd: Optional, if we're minimizing a chain of states
       type(NEB_PATH),optional,intent(inout) :: path
       ! agrecocmplx: optional argument for complex NGWFs
       logical, intent(in), optional :: is_cmplx

       ! <<< local variables >>>
       type(castep_model) :: cmdl
       integer :: ion_i, ierr, ion_in_species
       ! agrecocmplx
       logical :: loc_cmplx
       ! jcap: Further local variables
       integer :: iregion,jat(pub_ngwf_regions_ngroups)
       ! lk: string that holds the form of
       !     energy that is being minimized.
       character(len=80) :: energy_str

       !Set up flags for turning on/off experimental bits
       hess_init_method=hess_init_pfrommer
       debug_square_order=.false. ! for testing order of phonon frequency sqrt / sum
       bulk_freq_progress=.false.
       debug_reinit=0
       debug_init_scale=1.0_dp
       do_reinit=.true.
       update_init_on_remove=.false.

       ! agrecocmplx
       if (present(is_cmplx)) then
          loc_cmplx = is_cmplx
       else
          loc_cmplx = .false.
       end if

       ! aam: total number of ionic constraints
       num_ionic_constraints = 0
       do ion_i=1,mdl%nat
          select case (mdl%elements(ion_i)%ion_constraint_type)
          case ('NONE')  ; continue
          case ('LINE')  ; num_ionic_constraints = num_ionic_constraints + 1
          case ('PLANE') ; num_ionic_constraints = num_ionic_constraints + 2
          case ('FIXED') ; num_ionic_constraints = num_ionic_constraints + 3
          case default
             call utils_abort('Error in geometry_optimise: illegal &
                  &value for mdl%elements(i)%ion_constraint_type, i=', ion_i)
          end select
       enddo

       ! ars 0 <= num_ionic_constraints <= 3*mdl%nat
       ! aam: Analyse constraints
       num_ionic_constraints=min(3*mdl%nat,num_ionic_constraints)  !catch any sillies
       num_ionic_constraints=max(0,num_ionic_constraints)          !catch any sillies
       if (num_ionic_constraints>0) constrain_ions=.true.
       ! ars : correct initialisation of fix_all_ions
       if (num_ionic_constraints==3*mdl%nat) then
          fix_all_ions=.true.
       else
          fix_all_ions=.false.
       end if

       ndof = 3*mdl%nat - num_ionic_constraints ! aam: cell is fixed ; no symmetry

          ndim = 9 + 3*mdl%nat  !we always allow for both cell & ions here regardless
                                !of constraints to allow for changes upon restart, etc.

       !Make sure there is something here to do
       if (fix_all_ions) then
          if (pub_on_root) write (stdout,'(a)') &
          'WARNING - there is nothing to optimise - skipping relaxation'
          return
       end if

       ! aam: Write a little banner
       if (pub_on_root) then
          write(stdout,'(a)') ''
          if (pub_geom_continuation) then
             write(stdout,1)
        write(stdout,2) 'Resuming previous ONETEP Geometry Optimisation'
             write(stdout,1)
          else
             write(stdout,1)
             write(stdout,3) 'Starting ONETEP Geometry Optimisation'
             write(stdout,1)
          endif
          write(stdout,'(/a,a,a)') &
               ' Geometry Optimisation: output file is "',trim(output_file),'"'
       endif

       if (pub_on_root) write(stdout,*) ''

        ! fill CASTEP-style cell data structures from mdl%cell and "elements"
        ! kaw: Directly analogous to cks changes for energy calculation
       if (mdl%nat_classical > 0) then
        call copy_to_castep_cell(current_cell,mdl%cell,mdl%nat,mdl%elements, &
             mdl%nat_classical,mdl%classical_elements)
       else
        call copy_to_castep_cell(current_cell,mdl%cell,mdl%nat,mdl%elements, &
             mdl%nat_classical)
       end if

       ! Consistency check
       call utils_assert(current_cell%num_ions>0, 'Error in geometry_optimise: &
               &current_cell uninitialised in geometry_optimise')

       ! decide whether to do bfgs or lbfgs...
       select case(trim(adjustl(pub_geom_method)))
       case('CARTESIAN')
          if(pub_geom_lbfgs) then
             do_lbfgs=.true.
             do_bfgs=.false.
             do_tpsd=.false.
          else
             do_lbfgs=.false.
             do_bfgs=.true.
             do_tpsd=.false.
          end if
       case('DELOCALIZED','BFGS')
          do_lbfgs=.false.
          do_tpsd=.false.
          do_bfgs=.true.
       case('LBFGS')
          do_lbfgs=.true.
          do_bfgs=.false.
          do_tpsd=.false.
       case('TPSD')
          do_lbfgs=.false.
          do_bfgs=.false.
          do_tpsd=.true.
       end select

       ! set flags determining whether there is l/bfgs data in the cmdl...
       ! there is none for the moment...
       cmdl%bfgs_optimisation=.false.
       cmdl%lbfgs_optimisation=.false.
       cmdl%lbfgs_num_updates=0
       ! Allocate model data
       call castep_model_alloc(cmdl,mdl%nat+mdl%nat_classical,ndim,constrain_ions,do_lbfgs)

       if(do_lbfgs) then ! allocate the LBFGS intermediate matrices and vectors...
          if (.not.allocated(R_mat)) then
             allocate(R_mat(1:pub_geom_lbfgs_block_length,1:pub_geom_lbfgs_block_length),stat=ierr)
             call utils_alloc_check('geometry_optimise','R_mat',ierr)
          end if
          if (.not.allocated(YY_mat)) then
             allocate(YY_mat(1:pub_geom_lbfgs_block_length,1:pub_geom_lbfgs_block_length),stat=ierr)
             call utils_alloc_check('geometry_optimise','YY_mat',ierr)
          end if

          if (.not.allocated(tmpmm)) then
             allocate(tmpmm(1:pub_geom_lbfgs_block_length*2,1:pub_geom_lbfgs_block_length*2),stat=ierr)
             call utils_alloc_check('geometry_optimise','tmpmm',ierr)
          end if
          if (.not.allocated(tmp2m)) then
             allocate(tmp2m(1:pub_geom_lbfgs_block_length*2),stat=ierr)
             call utils_alloc_check('geometry_optimise','tmp2m',ierr)
          end if

          if (.not.allocated(lbfgs_init_inv_hess)) then
             allocate(lbfgs_init_inv_hess(1:3,1:ndim),stat=ierr)
             call utils_alloc_check('geometry_optimise','lbfgs_init_inv_hess',ierr)
          end if
          !       if (.not.allocated(lbfgs_init_hess)) then
          !          if (constrain_cell.and..not.fix_all_cell) allocate(lbfgs_init_hess(1:3,1:ndim))
          !       end if

          if (.not.allocated(tmpn)) then
             allocate(tmpn(1:ndim),stat=ierr)
             call utils_alloc_check('geometry_optimise','tmpn',ierr)
          end if
          if (.not.allocated(D_mat)) then
             allocate(D_mat(1:pub_geom_lbfgs_block_length),stat=ierr)
             call utils_alloc_check('geometry_optimise','D_mat',ierr)
          end if
       end if


       ! Over-ride the density kernel and tightbox_ngwf write flags if necessary
       if (.not.pub_write_denskern .and. pub_geom_reuse_dk_ngwfs.and..not.pub_edft) then
          pub_write_denskern = .true.
          if (pub_on_root) write(stdout,*)  'Geometry Optimisation: &
               &write_denskern parameter over-ridden to T'
       endif
       if (.not.pub_write_tightbox_ngwfs .and. pub_geom_reuse_dk_ngwfs .and. &
            (pub_maxit_ngwf_cg > 0)) then
          pub_write_tightbox_ngwfs = .true.
          if (pub_on_root) write(stdout,*) 'Geometry Optimisation: &
               &write_tightbox_ngwfs parameter over-ridden to TRUE'
       endif
       ! ars: if this is a EDFT calculation, we need to save the Hamiltonian
       if (pub_edft.and.(.not.pub_write_hamiltonian).and.pub_geom_reuse_dk_ngwfs) then
          pub_write_hamiltonian = .true.
          if (pub_on_root) write(stdout,*)  'Geometry Optimisation: &
               &write_hamiltonian parameter over-ridden to T'
       endif
       ! End over-rides


       !**************************************************************!
       !                     Copy model data                          !
       !**************************************************************!

       ! Name of continuation file
       write(cont_file,'(a,a)') trim(pub_rootname),'.continuation'
#ifdef HDF5
       cont_file=trim(cont_file)//'.h5'
#endif

       geom_continuation_string='NONE '

       if (pub_geom_continuation) then

          ! Read continuation data into model and current_cell.
          ! NB This subroutine overwrites current_cell%ionic_positions
          call geom_opt_continuation_read(cmdl,mdl,cont_file,loc_cmplx)

          call castep_cell_2_elements(current_cell,mdl%elements,mdl%nat)

          ! jcap: copy the continuation data into the appropriate subregions
          jat = 0
          do ion_i=1,mdl%nat
             ! Get the region counter
             iregion = mdl%elements(ion_i)%region
             ! Count the no. of atoms allocated to each region
             jat(iregion) = jat(iregion) + 1
             ! Copy the co-ordinates over
             mdl%regions(iregion)%elements(jat(iregion))%centre = &
                  &mdl%elements(ion_i)%centre
          end do

       else

          call castep_cell_copy(current_cell,cmdl%cell)

          if (pub_write_positions) call geom_output(cmdl,mdl)

          ! If this is not a continuation of a previous geometry optimisation:
          ! Find ground state energy and forces for initial configuration
          ! agrecocmplx: use optional argument to specify usage of complex NGWFs
          if (present(path)) then
             ! Part of a chain-of-states minimization
             call neb_energy_and_force(total_energy,total_forces,mdl,path, &
                  is_cmplx=loc_cmplx)
          else
             call energy_and_force_calculate(total_energy,total_forces,mdl, &
                  is_cmplx=loc_cmplx)
          end if

          ! Total energy
          cmdl%total_energy = total_energy

          ! Forces
          do ion_i=1,mdl%nat                                   ! Each ion has its own species
             do ion_in_species=1,current_cell%max_ions_in_species   ! ie max_ions_species = 1
                cmdl%forces(1,ion_in_species,ion_i) = total_forces(1,ion_i)
                cmdl%forces(2,ion_in_species,ion_i) = total_forces(2,ion_i)
                cmdl%forces(3,ion_in_species,ion_i) = total_forces(3,ion_i)
                if (constrain_ions) then ! aam: keep copy of unconstrained forces
                   cmdl%orig_forces(1,ion_in_species,ion_i) = total_forces(1,ion_i)
                   cmdl%orig_forces(2,ion_in_species,ion_i) = total_forces(2,ion_i)
                   cmdl%orig_forces(3,ion_in_species,ion_i) = total_forces(3,ion_i)
                endif
             enddo
          enddo

          ! BFGS data
          cmdl%bfgs_iteration  = 0
          cmdl%lbfgs_iteration = 0
          cmdl%tpsd_iteration = 0
          if (pub_on_root) then
             if(do_bfgs) then
                cmdl%bfgs_inv_Hessian = 0.0_dp
             else if(do_lbfgs) then
                cmdl%lbfgs_position_updates=0.0_dp
                cmdl%lbfgs_gradient_updates=0.0_dp
             end if
          end if
       endif

       cmdl%found_ground_state = .true.
       cmdl%stress             = 0.0_dp
       cmdl%strain             = 0.0_dp

       ! aam: apply ionic constraints
       if (constrain_ions) call forces_apply_constraints(cmdl%forces,mdl)

       call castep_cell_copy(current_cell,cmdl%cell)
       call castep_cell_copy(current_cell,cmdl%orig_cell)
       call castep_cell_copy(current_cell,cmdl%ref_cell)

       !**************************************************************!
       !                    Model data copy completed                 !
       !**************************************************************!


       ! Over-ride the density kernel and tightbox_ngwf read flags if necessary
       ! This is required for back-tracking in the BFGS algorithm
       ! cks, 2008/08/26: Added the ability to prevent the usage of the
       ! cks: density kernel and NGWFs from the previous geometry step
       ! cks: as this can occasionally cause convergence problems in the energy.
       if (.not.pub_read_denskern.and.pub_geom_reuse_dk_ngwfs.and..not.pub_edft) then
          pub_read_denskern = .true.
          if (pub_on_root) write(stdout,*)  'Geometry Optimisation: &
               &read_denskern parameter over-ridden to T'
       endif
       if (.not.pub_read_tightbox_ngwfs .and.  pub_geom_reuse_dk_ngwfs) then
          pub_read_tightbox_ngwfs = .true.
          if (pub_on_root) write(stdout,*) 'Geometry Optimisation: &
               &read_tightbox_ngwfs parameter over-ridden to TRUE'
       endif
       ! ars: if this is an EDFT calculation, we need to read the Hamiltonian
       if (pub_edft.and.(.not.pub_read_hamiltonian).and. pub_geom_reuse_dk_ngwfs) then
          pub_read_hamiltonian = .true.
          if (pub_on_root) write(stdout,*)  'Geometry Optimisation: &
               &read_hamiltonian parameter over-ridden to T'
       endif
       ! cks: do not do penalty if NGWFs and kernel are re-used
       if (pub_geom_reuse_dk_ngwfs) pub_maxit_pen = 0
       ! End over-rides

       ! aam: output initial atomic positions in fractional coordinates
       if (pub_write_positions) call geom_output(cmdl,mdl)
       call comms_barrier

       !Now select geometry minimizer (might want to add others in the future)
       select case (pub_geom_method)
       case ('CARTESIAN','BFGS','LBFGS')

          ! agrecocmplx
          call geom_BFGS(cmdl,mdl,output_file,path,loc_cmplx)

       case ('TPSD')
          call geom_tpsd(cmdl,mdl,output_file,path,loc_cmplx)

       case ('DELOCALIZED')

          ! agrecocmplx
          if (present(path)) then
             call geom_DELOC(cmdl,mdl,output_file,path,loc_cmplx,ierr)
          else
             call geom_DELOC(cmdl,mdl,output_file,is_cmplx=loc_cmplx,status=ierr)
          end if

          if (ierr/=0) then
             pub_geom_method = 'CARTESIAN'
             ! agrecocmplx
             if (present(path)) then
                call geom_BFGS(cmdl,mdl,output_file,path,loc_cmplx)
             else
                call geom_BFGS(cmdl,mdl,output_file,is_cmplx=loc_cmplx)
             end if
          endif

       case default
          call utils_abort('Error in geom_optimise - unrecognised geometry&
               &optimisation request='//trim(pub_geom_method)//'.')
       end select

       ! Deallocate memory

       call castep_model_dealloc(cmdl,do_lbfgs)

       call castep_cell_dealloc(current_cell)

       ! lk: print total energy after geometry optimization.
       !     Note that this is not printed in case of Nudged
       !     Elastic Band and Ensemble DFT calculations.
       if (pub_on_root.and..not.present(path)) then
          if(pub_edft) then
             energy_str='Free Energy'
          else
             energy_str='Total Energy'
          end if
          write(stdout,'(a)') '================================&
               &================================================'
          write(stdout,'(a,i10,a)') '           &
               &Finished geometry optimization after', &
               iteration,' iteration(s).'
          write(stdout,'(a,a12,es21.12e3,a,a)') '           ',&
               &energy_str,cmdl%total_energy,' ',trim(energy_label)
          write(stdout,'(a)') '================================&
               &================================================'
       end if

1      format(80('='))
2      format(16('<'),1x,a,1x,16('>'))
3      format(20('<'),1x,a,1x,21('>'))

       return

     end subroutine geometry_optimise

     !=========================================================================!
     ! Find ground state ionic positions for given model.                      !
     ! Based on 'Two-Point Step Size Gradient Methods' by Barzilai and Borwein !
     ! 1988.                                                                   !
     !-------------------------------------------------------------------------!
     ! Arguments:                                                              !
     !   cmdl, intent=inout, the model to be optimised                         !
     !   output_file, intent=in, the filename to which results will be written.!
     !-------------------------------------------------------------------------!
     ! Parent module variables used:                                           !
     !   on_root used for diagnostics                                          !
     !-------------------------------------------------------------------------!
     ! Modules used:                                                           !
     !   model for type definitions                                            !
     !   parameters                                                            !
     !   cell for cell vectors and constraints                                 !
     !   io for output and unit conversion routines                            !
     !-------------------------------------------------------------------------!
     ! Key Internal Variables:                                                 !
     !   x_vec       =  ndim   vector of strains & fractional coords           !
     !   f_vec       =  ndim   vector of stresses & forces (without the metric)!
     !   delta_vec   =  ndim   vector of x_vec update: new_x=old_x+lambda*delta!
     !   lambda      =  scalar = line minimization parameter along the proposed!
     !                  BFGS update direction delta                            !
     !   trial_OK    =  logical flag to say whether or not trial step has gone !
     !                  up or down in enthalpy. Used if need to backtrack.     !
     !                  ditto for line_OK and quad_OK.                         !
     !                                                                         !
     !-------------------------------------------------------------------------!
     ! Necessary conditions:                                                   !
     !   cmdl has already been properly read and initialised                   !
     !-------------------------------------------------------------------------!
     ! Written by Jolyon Aarons, v1.0, 23/05/2022.                             !
     !=========================================================================!
     subroutine geom_tpsd(cmdl,mdl,output_file,path,is_cmplx)
       use constants, only: dp, stdout, NORMAL, VERBOSE
       use comms, only: pub_on_root, pub_root_proc_id, comms_barrier, &
            comms_bcast
       !       use image_comms, only: pub_my_image
       use linalg, only: linalg_invert_sym_matrix
       use model_type, only: MODEL
       use rundat, only: pub_geom_energy_tol, pub_geom_output_detail, &
            pub_geom_max_iter,pub_geom_backup_iter,pub_rootname, &
            pub_write_positions,pub_geom_reset_dk_ngwfs_iter,pub_geom_reuse_dk_ngwfs, &
            pub_read_denskern,pub_read_tightbox_ngwfs,pub_maxit_ngwf_cg, &
            pub_edft, pub_read_hamiltonian, pub_hubbard,&
            pub_edft_spin_fix, pub_edft_spin_fix_orig, pub_geom_frequency_est
       use simulation_cell, only: castep_cell, castep_model, &
            castep_cell_alloc, castep_cell_copy, castep_cell_dealloc
       use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
            utils_open_unit_check, utils_close_unit_check, utils_assert, utils_abort
       use neb, only: neb_path, neb_sync_path, neb_converged
       use services, only: services_flush

       implicit none
       type(castep_model), intent(inout) :: cmdl
       type(MODEL),        intent(inout) :: mdl
       character(len=*),   intent(in)    :: output_file
       type(neb_path),optional,intent(inout) :: path
       logical,            intent(in)    :: is_cmplx

       real (kind=dp), allocatable, dimension(:)     :: x_vec, f_vec
       real (kind=dp), allocatable, dimension(:)     :: x_vec_old, f_vec_old
       real (kind=dp), allocatable, dimension(:)     :: s_dir
       real (kind=dp), allocatable, dimension(:)     :: x_diff, g_diff

       logical :: done_trial

       integer :: loc_iteration, iter ! jd: renamed 'iteration' to 'loc_iteration' to expose module-wide var of the same
       integer :: ierr                !     name which we need to set to get 'geom_write_trajectory' to work correctly
       real(kind=DP) :: enthalpy
       real(kind=DP) :: old_enthalpy, trial_enthalpy

       real(kind=DP) :: dE
       real(kind=DP) :: Fmax   !max Force
       real(kind=DP) :: dRmax  !max displacement
       real(kind=DP) :: Smax   !max Stress
       logical       :: converged, converged_dE, converged_Fmax, converged_dRmax, converged_Smax
       integer :: geom_tpsd_step_type
       real(kind=dp) :: old_f_dot_delta, trial_f_dot_delta

       real(kind=dp) :: metric(3,3)
       integer :: ispec
       real(kind=dp) :: avg_mass
       integer :: SD_init
       logical :: doing_restart

       ! JA: Allocate and initialize vectors
       ! JA: positions
       allocate (x_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_tpsd','x_vec',ierr)
       x_vec=0.0_dp
       allocate (x_vec_old(1:ndim),stat=ierr)
       call utils_alloc_check('geom_tpsd','x_vec_old',ierr)
       x_vec_old=0.0_dp
       ! JA: forces
       allocate (f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_tpsd','f_vec',ierr)
       f_vec=0.0_dp
       allocate (f_vec_old(1:ndim),stat=ierr)
       call utils_alloc_check('geom_tpsd','f_vec_old',ierr)
       f_vec_old=0.0_dp
       ! JA: search_direction
       allocate (s_dir(1:ndim),stat=ierr)
       call utils_alloc_check('geom_tpsd','s_dir',ierr)
       s_dir=0.0_dp
       ! JA: difference vectors
       allocate (x_diff(1:ndim),stat=ierr)
       call utils_alloc_check('geom_tpsd','x_diff',ierr)
       x_diff=0.0_dp
       allocate (g_diff(1:ndim),stat=ierr)
       call utils_alloc_check('geom_tpsd','g_diff',ierr)
       g_diff=0.0_dp

       ! JA: type of tpsd curvature estimate (1 or 2)
       geom_tpsd_step_type=1

       ! JA: set the geom string
       geom_string = 'TPSD'
       call comms_bcast(pub_root_proc_id,geom_string)

       ! JA: initialized all convergence checks to false
       converged=.false.; converged_dE=.false.; converged_Fmax=.false.; converged_dRmax=.false.; converged_Smax=.false.

       ! JA: get the vectors from the model
       call geom_mdl_to_xvec(cmdl,x_vec_old,mdl)
       call geom_mdl_to_fvec(cmdl,f_vec_old,enthalpy,mdl)

       ! JA: are we doing a restart?
       if(cmdl%tpsd_iteration==0) then
          doing_restart = .false.
       else
          doing_restart = .true.
       end if

       iteration = cmdl%tpsd_iteration ! jd: Needed for get geom_write_trajectory() to work

       if(doing_restart) then
          s_dir = cmdl%tpsd_alpha*f_vec_old
       else

          ! jd: Write .geom file
          call geom_write_trajectory(cmdl,enthalpy,output_file,mdl)

          ! JA: do a step of steepest descent first
          if (pub_on_root) write (stdout,4) 'Starting init '//trim(geom_string)//' iteration ',0,' (SD) ...'

          ! JA: Get an approximate inverse Hessian (from Pfrommer)
          ! JA: this is so we can appropriately scale the first step of SD.
          metric=matmul(cmdl%cell%real_lattice,transpose(cmdl%cell%real_lattice))
          call linalg_invert_sym_matrix(metric,3)

          avg_mass=0.0_dp
          do ispec=1,cmdl%cell%num_species - mdl%nat_classical
             avg_mass=avg_mass+cmdl%cell%species_mass(ispec)*cmdl%cell%num_ions_in_species(ispec)
          end do
          avg_mass=1822.888485_dp*avg_mass/cmdl%cell%num_ions

          metric=metric/(avg_mass*pub_geom_frequency_est**2) !(metric^-1)/(mass*omega^2)

          ! JA: choose between Pfrommer initialized SD and simple scalar init SD

          SD_init = 2 ! 1=simple, 2=Pfrommer

          if(SD_init==1) then
             ! JA: with step size of 0.001*fvec
             cmdl%tpsd_alpha = 0.001_dp

             ! JA: set initial SD search direction
             s_dir = cmdl%tpsd_alpha*f_vec_old

          elseif(SD_init==2) then
             do ispec=1,cmdl%cell%num_species - mdl%nat_classical
                s_dir(10+(ispec-1)*3:9+ispec*3) = matmul(f_vec_old(10+(ispec-1)*3:9+ispec*3),metric)
             end do
          end if
       end if

       ! JA: set new SD positions
       x_vec = x_vec_old + s_dir

       ! JA: Put new positions in model, run a DFT calc & retrieve forces.
       call geom_xvec_to_mdl(x_vec,cmdl,mdl)
       call geom_get_forces(cmdl,mdl,path,is_cmplx,en_sync=.true.)
       call geom_mdl_to_fvec(cmdl,f_vec,enthalpy,mdl)

       old_f_dot_delta = dot_product(f_vec,f_vec)

       ! JA: Check convergence
       call geom_converged(cmdl,x_vec_old,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
            & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

       loc_iteration = 0

       ! if(do_lbfgs) iteration=cmdl%lbfgs_iteration
       call geom_output_converged(loc_iteration,enthalpy,dE,Fmax,dRmax, &
            & converged_dE,converged_Fmax,converged_dRmax)

       old_enthalpy = enthalpy

       if (present(path)) converged = neb_converged(converged, loc_iteration, dE, Fmax, dRmax, path)

       ! JA: All we will ever do in tpsd is a trial... there is no option.
       done_trial = .true.

       trial_f_dot_delta = old_f_dot_delta

       ! JA: Do the TPSD loop
       tpsd_loop: do loc_iteration = 1,pub_geom_max_iter

          cmdl%tpsd_iteration=cmdl%tpsd_iteration+1

          if (pub_on_root) write (stdout,4) 'Starting '//trim(geom_string)//' iteration ',cmdl%tpsd_iteration,' ...'

          ! JA: Now we start tpsd by setting up the difference vectors
          x_diff = x_vec - x_vec_old
          g_diff = f_vec_old - f_vec ! Difference in gradient, not force

          ! JA: Do one of the 2 options for tpsd curvature estimate.
          if(geom_tpsd_step_type == 1) then
             cmdl%tpsd_alpha = dot_product(x_diff,g_diff)/dot_product(g_diff,g_diff)
          elseif(geom_tpsd_step_type == 2) then
             cmdl%tpsd_alpha = dot_product(x_diff,x_diff)/dot_product(x_diff,g_diff)
          end if

          ! JA: Backup for the next iteration
          x_vec_old = x_vec
          f_vec_old = f_vec

          old_f_dot_delta = trial_f_dot_delta
          trial_f_dot_delta = dot_product(f_vec,f_vec)

          ! JA: Make a new search direction
          s_dir = cmdl%tpsd_alpha*f_vec_old

          ! JA: Update positions based on TPSD search direction
          x_vec = x_vec_old + s_dir

          ! JA: Put new positions in model, run a DFT calc & retrieve forces.
          old_enthalpy = enthalpy
          call geom_xvec_to_mdl(x_vec,cmdl,mdl)
          call geom_get_forces(cmdl,mdl,path,is_cmplx,en_sync=.true.)
          call geom_mdl_to_fvec(cmdl,f_vec,enthalpy,mdl)
          trial_enthalpy = enthalpy

          ! JA: Write an update to screen
          call geom_tpsd_write_summary

          if (mod(cmdl%tpsd_iteration,pub_geom_backup_iter)==0) then
             call geom_opt_continuation_write(cmdl, mdl, cont_file, is_cmplx)
          end if

          iteration = cmdl%tpsd_iteration ! jd: Needed to get geom_write_trajectory() to work
          call geom_write_trajectory(cmdl,enthalpy,output_file,mdl)

          ! JA: test convergence
          call geom_converged(cmdl,x_vec_old,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
               & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

          !JA: check if we are done yet.
          if (converged) then
             !JA: tell user the calculation just converged
             call geom_output_converged(cmdl%tpsd_iteration,enthalpy,dE,Fmax,dRmax,&
                  & converged_dE,converged_Fmax,converged_dRmax)

             ! jd: Following (l)bfgs and deloc logic here
             call geom_write_trajectory(cmdl,enthalpy,output_file,mdl)

             if (present(path)) then
                converged = neb_converged(converged, cmdl%tpsd_iteration, dE, Fmax, dRmax, path)
                if (converged) then
                   exit tpsd_loop
                else
                   cycle tpsd_loop
                end if
             else
                exit tpsd_loop
             end if
          end if


       end do tpsd_loop


       ! JA: Clean up memory for vectors.
       deallocate (x_vec,stat=ierr)
       call utils_dealloc_check('geom_tpsd','x_vec',ierr)
       deallocate (x_vec_old,stat=ierr)
       call utils_dealloc_check('geom_tpsd','x_vec_old',ierr)
       deallocate (f_vec,stat=ierr)
       call utils_dealloc_check('geom_tpsd','f_vec',ierr)
       deallocate (f_vec_old,stat=ierr)
       call utils_dealloc_check('geom_tpsd','f_vec_old',ierr)
       deallocate (s_dir,stat=ierr)
       call utils_dealloc_check('geom_tpsd','s_dir',ierr)
       deallocate (x_diff,stat=ierr)
       call utils_dealloc_check('geom_tpsd','x_diff',ierr)
       deallocate (g_diff,stat=ierr)
       call utils_dealloc_check('geom_tpsd','g_diff',ierr)

4   format(/,80('='),/,1x,a,i4,a,/,80('='))

     contains

       subroutine geom_tpsd_write_summary()
         !=========================================================================!
         ! contain'ed subroutine to write out the TPSD minimisation summary.       !
         !-------------------------------------------------------------------------!
         ! Arguments:                                                              !
         !    none                                                                 !
         !-------------------------------------------------------------------------!
         ! Parent module variables used:                                           !
         !    none                                                                 !
         !-------------------------------------------------------------------------!
         ! Modules used:                                                           !
         !    none                                                                 !
         !-------------------------------------------------------------------------!
         ! Key Internal Variables:                                                 !
         !    none                                                                 !
         !-------------------------------------------------------------------------!
         ! Necessary conditions:                                                   !
         !-------------------------------------------------------------------------!


         use constants, only: dp, stdout
         use services, only: services_flush

         !local vars
         character(len=80) :: data_string
         character(len=80) :: divider_string
         character(len=80) :: label_string
         integer           :: string_index
         integer           :: len_label, len_lambda, len_fdelta, len_enthalpy
         real (kind=dp)    :: enthalpy, f_dot_delta

         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Entering geom_tpsd_write_summary'

         !initialise strings
         data_string     = repeat(' ',len(data_string))
         divider_string  = repeat(' ',len(divider_string))
         label_string    = repeat(' ',len(label_string))
         !                  123456789 2345678901 23456789012345678901 23456789012345678
         divider_string  = '+--------+----------+--------------------+-----------------+'
         label_string    = '|  Step  |  lambda  | F.delta (Ha/a0)**2 |  enthalpy (Ha)  |'
         !e.g.              | trial step   -100.12345    -100.12345    -123456.123456

         if(pub_on_root) then

            !header
            write(stdout,*) ' '
            select case (trim(geom_string))
            case ('TPSD')
               write(stdout,77) divider_string
               write(stdout,77) label_string
               write(stdout,77) divider_string
            end select

!             !previous                                                    '123456789012'
!             string_index=1
             len_label   =10
             len_lambda  =11
             len_fdelta  =21
             len_enthalpy=18
!             write(data_string(string_index:string_index+len_label),   1) '  previous  '
!             string_index=string_index+len_label
!             write(data_string(string_index:string_index+len_lambda),  2) 0.0_dp
!             string_index=string_index+len_lambda
!             f_dot_delta=old_f_dot_delta
!             if (abs(f_dot_delta)>1.0E4.or.abs(f_dot_delta)<1.0E-5) then
!                write(data_string(string_index:string_index+len_fdelta), 33) f_dot_delta
!             else
!                write(data_string(string_index:string_index+len_fdelta),  3) f_dot_delta
!             end if
!             string_index=string_index+len_fdelta
!             enthalpy=old_enthalpy
!             if (abs(enthalpy)>1.0E7.or.abs(enthalpy)<1.0E-5) then
!                write(data_string(string_index:string_index+len_enthalpy),34) enthalpy
!             else
!                write(data_string(string_index:string_index+len_enthalpy),4) enthalpy
!             end if
!             string_index=string_index+len_enthalpy
!             !cross check all is well (NB strings start from 1 not 0 ...)
!             if (string_index>71) write (stdout,*) 'geom_tpsd_write_summary: format problem?'
!             if(trim(geom_string)=='TPSD')  write(stdout,77) data_string

            !trial                                                          '123456789012'
            if (done_trial) then
               string_index=1
               write(data_string(string_index:string_index+len_label),   1) '  TPSD'
               string_index=string_index+len_label
               write(data_string(string_index:string_index+len_lambda),  2) cmdl%tpsd_alpha
               string_index=string_index+len_lambda
                f_dot_delta=trial_f_dot_delta
                if (abs(f_dot_delta)>1.0E4.or.abs(f_dot_delta)<1.0E-5) then
                   write(data_string(string_index:string_index+len_fdelta), 33) f_dot_delta
                else
                   write(data_string(string_index:string_index+len_fdelta),  3) f_dot_delta
                end if
               string_index=string_index+len_fdelta
               enthalpy=trial_enthalpy
               if (abs(enthalpy)>1.0E7.or.abs(enthalpy)<1.0E-5) then
                  write(data_string(string_index:string_index+len_enthalpy),34) enthalpy
               else
                  write(data_string(string_index:string_index+len_enthalpy),4) enthalpy
               end if
               string_index=string_index+len_enthalpy
               if(trim(geom_string)=='TPSD') write(stdout,77) data_string
            end if

            !footer
            if(trim(geom_string)=='TPSD') write(stdout,77) divider_string
            write(stdout,*) ' '

         end if          !pub_on_root
         call services_flush

         if (pub_debug_on_root) write(stdout,'(a)') &
              'DEBUG: Leaving geom_tpsd_write_summary'

1        format('|',a6,2x,'|')         !1+6+1+1=9
2        format(1x,f8.6,1x,'|')        !1+8+1+1=11
3        format(1x,f18.6,1x,'|')       !1+18+1+1=21
33       format(1x,es18.3e3,1x,'|')    !1+18+1+1=21
4        format(1x,f15.6,1x,'|')       !1+15+1+1=18
34       format(1x,es15.6e3,1x,'|')    !1+15+1+1=18
         !14+14+14+18=60
77       format(1x,a60,' <-- min TPSD')   !1+60+14 => ends col 75

       end subroutine geom_tpsd_write_summary

     end subroutine geom_tpsd




     !---------------------------------------------------------------------------!
     !                P R I V A T E   R O U T I N E S                            !
     !---------------------------------------------------------------------------!

     subroutine geom_BFGS(cmdl,mdl,output_file,path,is_cmplx)
       !=========================================================================!
       ! Find ground state ionic positions and/or cell vectors for given model.  !
       ! Based upon "Relaxation of Crystals with the Quasi-Newton Method" by     !
       !   B.G. Pfrommer et al, J.Comp.Phys. v131, p233-240 (1997).              !
       !                                                                         !
       ! We are in the constant pressure, adiabatic ensemble, and so work with   !
       !   the enthalpy: H=Etot+p^ext.V and seek to minimize this.               !
       ! Minimizing wrt ionic coords -> equilibrium condition: force=0           !
       !   and min. wrt cell strains -> equilibrium condition: stress+pressure=0 !
       ! That is, we use the following sign convention for stress & pressure:    !
       !   A positive external pressure tends to compress the system, and        !
       !   a positive internal stress   tends to expand   the system.            !
       !                                                                         !
       ! NB We use the ref_cell as the strain reference if doing variable cell.  !
       !                                                                         !
       ! NB As of version 1.48 the minimiser can be either constant Ecut or      !
       !    constant #PW, and as of version 1.54 this can be set at run-time via !
       !    the fixed_npw parameter.                                             !
       !    Constant Ecut is most physical and works best if long way from eqm   !
       !   - but can get trapped into finding wrong minimum (eg Fe).             !
       !    Constant PW can give smoother convergence near the minimum but is    !
       !    physical. Might be tempted to try to hybridise these schemes, but I  !
       !    found you can get nasty discontinuities if switch mid-minimisation   !
       !     => better to set it once before minimisation starts and not change. !
       !                                                                         !
       ! NB If variable cell with fixed_npw=F then we must abandon our strictly  !
       !    downhill in energy strategy. As we have a changing number of plane   !
       !    waves then there will be discontinuities in the energy as we step    !
       !    around the cell vectors and so an E-based search will go wild. Hence !
       !    the search direction must be based upon reducing F & S rather than E.!
       !-------------------------------------------------------------------------!
       ! Arguments:                                                              !
       !   cmdl, intent=inout, the model to be optimised                          !
       !   output_file, intent=in, the filename to which results will be written.!
       !-------------------------------------------------------------------------!
       ! Parent module variables used:                                           !
       !   on_root used for diagnostics                                          !
       !-------------------------------------------------------------------------!
       ! Modules used:                                                           !
       !   model for type definitions                                            !
       !   parameters                                                            !
       !   cell for cell vectors and constraints                                 !
       !   io for output and unit conversion routines                            !
       !-------------------------------------------------------------------------!
       ! Key Internal Variables:                                                 !
       !   x_vec       =  ndim   vector of strains & fractional coords           !
       !   f_vec       =  ndim   vector of stresses & forces (without the metric)!
       !   delta_vec   =  ndim   vector of x_vec update: new_x=old_x+lambda*delta!
       !   lambda      =  scalar = line minimization parameter along the proposed!
       !                  BFGS update direction delta                            !
       !   trial_OK    =  logical flag to say whether or not trial step has gone !
       !                  up or down in enthalpy. Used if need to backtrack.     !
       !                  ditto for line_OK and quad_OK.                         !
       !                                                                         !
       !-------------------------------------------------------------------------!
       ! Necessary conditions:                                                   !
       !   cmdl has already been properly read and initialised                   !
       !-------------------------------------------------------------------------!
       ! Written by Matt Probert, v1.0, 11/05/2002                               !
       !-------------------------------------------------------------------------!
       ! Modified for ONETEP by Arash A Mostofi, 2004                            !
       ! Modified for LBFGS by Jolyon Aarons, v2.0, 26/06/2011                   !
       !=========================================================================!

       use constants, only: dp, stdout, NORMAL, VERBOSE
       use comms, only: pub_on_root, pub_root_proc_id, comms_barrier, &
            comms_bcast
       use image_comms, only: pub_my_image
       use hubbard_init, only: h_species
       use linalg, only: linalg_invert_sym_matrix
       use model_type, only: MODEL
       use rundat, only: pub_geom_energy_tol,pub_geom_frequency_est, &
            pub_geom_output_detail, pub_geom_max_iter,pub_geom_modulus_est, &
            pub_geom_backup_iter,pub_rootname, pub_write_positions, &
            pub_geom_reset_dk_ngwfs_iter,pub_geom_reuse_dk_ngwfs, &
            pub_read_denskern,pub_read_tightbox_ngwfs,pub_maxit_ngwf_cg, &
            pub_geom_lbfgs_max_updates, pub_edft, pub_read_hamiltonian, &
            pub_hubbard, pub_edft_spin_fix, pub_edft_spin_fix_orig
       use simulation_cell, only: castep_cell, castep_model, &
            castep_cell_alloc, castep_cell_copy, castep_cell_dealloc
       use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
            utils_open_unit_check, utils_close_unit_check, utils_assert, utils_abort
       use neb, only: neb_path, neb_sync_path, neb_converged
       use services, only: services_flush

       implicit none

       type(castep_model), intent(inout) :: cmdl
       type(MODEL),        intent(inout) :: mdl
       character(len=*),   intent(in)    :: output_file
       type(neb_path),optional,intent(inout) :: path
       ! agrecocmplx: whether to use complex NGWFs
       logical,            intent(in)    :: is_cmplx
       type(castep_cell) :: init_cell


       type(castep_cell)                            :: backup_cell
       real(kind=DP)                                :: backup_total_energy,backup_enthalpy
       real(kind=DP), allocatable, dimension(:,:,:) :: backup_forces
       real(kind=DP), dimension(1:6)                :: backup_stress
       real(kind=DP), dimension(1:3,1:3)            :: backup_strain
       real(kind=DP), allocatable, dimension(:)     :: backup_x_vec, backup_f_vec
       integer                                      :: backup_dk_unit,&
            backup_ngwf_unit(mdl%nsub), backup_ham_unit !backup_unit
       character (len=file_maxpath)                 :: backup_dk_file,&
            backup_ngwf_file(mdl%nsub), backup_ham_file !backup_file

       real (kind=dp), allocatable, dimension(:)     :: x_vec, old_x_vec !strains and fractional coords
       real (kind=dp), allocatable, dimension(:)     :: f_vec, old_f_vec !derivatives of enthalpy wrt x_vec
       real (kind=dp), allocatable, dimension(:)     :: delta_vec, old_delta_vec !update x_vec
       integer                                      :: io_status

       real(kind=DP), allocatable, dimension(:)     :: trial_f_vec, line_f_vec, quad_f_vec

       real(kind=DP) :: enthalpy, old_enthalpy, trial_enthalpy, line_enthalpy, quad_enthalpy

       logical        :: monotonic_E !can be switched off to allow uphill steps if desperate!
       logical       :: abort_BFGS  !have we got stuck?
       !    logical        :: found_shear
       real (kind=dp) :: fudge       !amount of uphill ignored if monotonic_E=false

       logical        :: trial_OK, line_OK, quad_OK
       logical        :: do_line, do_quad
       logical        :: rescale_lambda, lambda_bad
       logical        :: done_trial, done_line, done_quad
       logical       :: revert_old
       logical        :: do_update
       logical        :: re_initialise

       real(kind=DP) :: dE
       real(kind=DP) :: Fmax   !max Force
       real(kind=DP) :: dRmax  !max displacement
       real(kind=DP) :: Smax   !max Stress
       logical       :: converged, converged_dE, converged_Fmax, converged_dRmax, converged_Smax

       real (kind=dp) :: lambda_trial,lambda_line,lambda_quad,lambda_backtrack
       logical       :: lambda_diff

       real (kind=dp) :: old_f_dot_delta, trial_f_dot_delta, line_f_dot_delta, quad_f_dot_delta
       real (kind=dp) :: df_dot_delta, df_dot_delta_dlambda

       integer       :: i,ierr, isub
!       logical       :: re_initialise

       real (kind=dp) :: switch_Fmax

       real(kind=dp) :: geom_linmin_tol
    logical       :: use_enthalpy
    integer       :: monotonic_E_steps, use_enthalpy_steps
    integer       :: sp

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_BFGS'

    !-------------------------------------------------------------------------!
    !Calculate and print initial memory estimate                              !
    !-------------------------------------------------------------------------!
    if (do_lbfgs) then
       if(pub_geom_lbfgs_max_updates.eq.0) then
          lbfgs_block_length = 30  !need to set this before call to geom_calc_storage
       else
          lbfgs_block_length = pub_geom_lbfgs_max_updates
       end if
    end if

    !-------------------------------------------------------------------------!
    !Initialize data structures required                                      !
    !-------------------------------------------------------------------------!
    monotonic_E=.false.
    geom_linmin_tol=0.6_dp
    use_enthalpy=.false.
    monotonic_E_steps=10
    use_enthalpy_steps=10

       abort_BFGS=.false.
    do_update=.false.
       fudge=pub_geom_energy_tol*real(cmdl%cell%num_ions,kind=dp) !absolute fudge in enthalpy for uphill steps

       old_frequency_est=pub_geom_frequency_est             !we store the starting values as 'old' to test o/p at end
    old_modulus_est  =pub_geom_modulus_est
       new_frequency_est=pub_geom_frequency_est             !we store the working value as 'new' s.t. can cope with
       new_modulus_est  =pub_geom_modulus_est               !automatic and user updates to the estimates


    converged=.false.

    if(do_bfgs) geom_string="BFGS "
    if(do_lbfgs) geom_string="LBFGS"

    call comms_bcast(pub_root_proc_id,geom_string)

    ! ja531 -> March 2015 : Changed mind about this : DI inv Hessian is Cartesian!
!    call utils_assert(trim(geom_continuation_string)/="DI","Error: Cannot continue DI geometry optimisation as non-DI.")

    !We need a local cell as that at which the inv_Hessian was last initialised ...
    call castep_cell_alloc(init_cell,mdl%nat+mdl%nat_classical)
    !We need a local cell as that at which the inv_Hessian was last initialised ...
    call castep_cell_copy(cmdl%cell,init_cell)
    call castep_cell_copy(cmdl%cell,init_cell)


    if(pub_on_root) then
       if(do_bfgs) then
          if (.not.allocated(cmdl%bfgs_inv_Hessian)) then
             allocate (cmdl%bfgs_inv_Hessian(1:ndim,1:ndim),stat=ierr)
             call utils_alloc_check('geom_BFGS','cmdl%bfgs_inv_Hessian',ierr)
             cmdl%bfgs_inv_Hessian=0.0_dp
          end if
       else if(do_lbfgs) then
!          if(allocated(cmdl%bfgs_inv_Hessian)) then
!             deallocate(cmdl%bfgs_inv_Hessian,stat=ierr)
!             call utils_alloc_check('geom_BFGS','cmdl%bfgs_inv_Hessian',ierr)
!          end if

          if (.not.allocated(R_mat)) then
             allocate(R_mat(1:lbfgs_block_length,1:lbfgs_block_length),stat=ierr)
             call utils_alloc_check('geom_BFGS','R_mat',ierr)
          end if
          if (.not.allocated(YY_mat)) then
             allocate(YY_mat(1:lbfgs_block_length,1:lbfgs_block_length),stat=ierr)
             call utils_alloc_check('geom_BFGS','YY_mat',ierr)
          end if

          if (.not.allocated(tmpmm)) then
             allocate(tmpmm(1:lbfgs_block_length*2,1:lbfgs_block_length*2),stat=ierr)
             call utils_alloc_check('geom_BFGS','tmpmm',ierr)
          end if
          if (.not.allocated(tmp2m)) then
             allocate(tmp2m(1:lbfgs_block_length*2),stat=ierr)
             call utils_alloc_check('geom_BFGS','tmp2m',ierr)
          end if

          if (.not.allocated(cmdl%lbfgs_position_updates)) then
             allocate(cmdl%lbfgs_position_updates(1:ndim,1:lbfgs_block_length),stat=ierr)
             call utils_alloc_check('geom_BFGS','=>cmdl%lbfgs_position_updates',ierr)
          end if
          if (.not.allocated(cmdl%lbfgs_gradient_updates)) then
             allocate(cmdl%lbfgs_gradient_updates(1:ndim,1:lbfgs_block_length),stat=ierr)
             call utils_alloc_check('geom_BFGS','=>cmdl%lbfgs_gradient_updates',ierr)
          end if
          if (.not.allocated(lbfgs_init_inv_hess)) then
             allocate(lbfgs_init_inv_hess(1:3,1:ndim),stat=ierr)
             call utils_alloc_check('geom_BFGS','lbfgs_init_inv_hess',ierr)
          end if
          if (.not.allocated(lbfgs_init_hess)) then
             !             if (constrain_cell.and..not.fix_all_cell) then
             allocate(lbfgs_init_hess(1:3,1:ndim),stat=ierr)
             call utils_alloc_check('geom_BFGS','lbfgs_init_hess',ierr)
             !             end if
          end if

          if (.not.allocated(tmpn)) then
             allocate(tmpn(1:ndim),stat=ierr)
             call utils_alloc_check('geom_BFGS','tmpn',ierr)
          end if
          if (.not.allocated(D_mat)) then
             allocate(D_mat(1:lbfgs_block_length),stat=ierr)
             call utils_alloc_check('geom_BFGS','D_mat',ierr)
          end if

       !We need a local cell as that at which the inv_Hessian was last initialised ...

          if(update_init_on_remove) then
             allocate(old_x_vec_remove(1:ndim),stat=ierr)
             call utils_alloc_check('geom_BFGS','old_x_vec_remove',ierr)
             allocate(old_f_vec_remove(1:ndim),stat=ierr)
             call utils_alloc_check('geom_BFGS','old_f_vec_remove',ierr)
          end if

       end if
    end if
    !    call comms_barrier()
    call comms_bcast(pub_root_proc_id,lbfgs_block_length)

    !can sometimes loose strain if SCF convergence fails but it is recoverable ...
    !  if (maxval(abs(cmdl%strain))<=1.0E-6_dp) then
    !     !calculate strain from cells
    !     call algor_invert(3,cmdl%ref_cell%real_lattice,inverse=cmdl%strain) !using cmdl%strain as temp ref_cell^-1
    !     !h'=(1+e)h => (h^T)'=(h^T)(1+e)^T => (1+e)^T=(h^T)^-1.(h^T)' => e=[(h^T)^-1.(h^T)']^T-1
    !     cmdl%strain=transpose(matmul(cmdl%strain,cmdl%cell%real_lattice))-unit_matrix
    !  end if

       !Then we set up the ndim arrays
       allocate (x_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','x_vec',ierr)
       x_vec=0.0_dp

       allocate (old_x_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','old_x_vec',ierr)
       old_x_vec=0.0_dp

       allocate (f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','f_vec',ierr)
       f_vec=0.0_dp

       allocate (old_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','old_f_vec',ierr)
       old_f_vec=0.0_dp

       allocate (delta_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','delta_vec',ierr)
       delta_vec=0.0_dp

    allocate (old_delta_vec(1:ndim),stat=ierr)
    call utils_alloc_check('geom_BFGS','old_delta_vec',ierr)
    old_delta_vec=0.0_dp

    !then we setup the backup arrays as needed


    call castep_cell_alloc(backup_cell,mdl%nat+mdl%nat_classical)


    ! AAM: Set file names for NGWF and density kernel backup files
    ! jcap: loop over regions
    do isub=1,mdl%nsub
       if (mdl%nsub==1) then
          write(backup_ngwf_file(isub),'(2a)')&
               trim(pub_rootname),'.tightbox_ngwfs.bfgs_back'
       else
          write(backup_ngwf_file(isub),'(a,a,i0,a)')&
               trim(pub_rootname),'.tightbox_ngwfs_',isub,'.bfgs_back'
       end if
#ifdef HDF5
       backup_ngwf_file(isub)=trim(backup_ngwf_file(isub))//'.h5'
#endif

    end do
    write(backup_dk_file,'(2a)') trim(pub_rootname),'.dkn.bfgs_back'
#ifdef HDF5
    backup_dk_file=trim(backup_dk_file)//'.h5'
#endif

    write(backup_ham_file,'(2a)') trim(pub_rootname),'.ham.bfgs_back'
#ifdef HDF5
    backup_ham_file=trim(backup_ham_file)//'.h5'
#endif


       if (pub_on_root) write(stdout,'(/3a)') &
         &trim(geom_string)//': continuation file name is "',trim(cont_file),'"'


       !kaw: We only need forces for the QM ions so the final dimension of this was too large
       !allocate(backup_forces(1:3,1:cmdl%cell%max_ions_in_species,1:cmdl%cell%num_species),stat=ierr)
       allocate(backup_forces(1:3,1:cmdl%cell%max_ions_in_species,1:mdl%nat),stat=ierr)
       call utils_alloc_check('geom_BFGS','backup_forces',ierr)
       backup_forces=0.0_dp

       allocate (backup_x_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','backup_x_vec',ierr)
       backup_x_vec=0.0_dp

       allocate (backup_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','backup_f_vec',ierr)
       backup_f_vec=0.0_dp

       allocate (trial_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','trial_f_vec',ierr)
       trial_f_vec=0.0_dp

       allocate (line_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','line_f_vec',ierr)
       line_f_vec=0.0_dp

       allocate (quad_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','quad_f_vec',ierr)
       quad_f_vec=0.0_dp

       !-------------------------------------------------------------------------!
       !Sort out start from restart                                              !
       !-------------------------------------------------------------------------!
       call comms_bcast(pub_root_proc_id,cmdl%bfgs_iteration)
       call comms_bcast(pub_root_proc_id,cmdl%lbfgs_iteration)

       if (((do_bfgs).and.(cmdl%bfgs_iteration==0)).or.&
            &((do_lbfgs).and.(cmdl%lbfgs_iteration==0))) then
          restart=.false.
          if(do_bfgs) then
             if(pub_on_root) then
                if((hess_init_method.eq.hess_init_identity).or.(hess_init_method.eq.hess_init_scaled)) then
                   call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
                        init_method=hess_init_identity)
                else
                   call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl)
                end if
              end if
             cmdl%strain=0.0_dp
          else
             cmdl%lbfgs_num_updates=0
             call comms_bcast(pub_root_proc_id,cmdl%lbfgs_num_updates)
             if(pub_on_root) then
                if (.not.allocated(R_mat)) then
                   allocate(R_mat(1:lbfgs_block_length,1:lbfgs_block_length),stat=ierr)
                   call utils_alloc_check('geom_BFGS','R_mat',ierr)
                end if
                if (.not.allocated(YY_mat)) then
                   allocate(YY_mat(1:lbfgs_block_length,1:lbfgs_block_length),stat=ierr)
                   call utils_alloc_check('geom_BFGS','YY_mat',ierr)
                end if
                if (.not.allocated(tmpmm)) then
                   allocate(tmpmm(1:lbfgs_block_length*2,1:lbfgs_block_length*2),stat=ierr)
                   call utils_alloc_check('geom_BFGS','tmpmm',ierr)
                end if
                if (.not.allocated(tmp2m)) then
                   allocate(tmp2m(1:lbfgs_block_length*2),stat=ierr)
                   call utils_alloc_check('geom_BFGS','tmp2m',ierr)
                end if
                if((hess_init_method.eq.hess_init_identity)) then !.or.(hess_init_method.eq.hess_init_scaled)) then
                   init_scalar=1.0_dp
                   call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
                        init_method=hess_init_identity)
                else
                   call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl)
                end if

             end if
          end if

          cmdl%strain=0.0_dp
          !Put the strains wrt original cell and fractional atomic coords wrt current_cell into x_vec
          call geom_mdl_to_xvec(cmdl,x_vec,mdl)

          !Relax electrons in starting structure to get forces & stresses ...
          !     call geom_get_forces(cmdl,elements)

          !Get stresses and forces from cmdl and put into f_vec
          call geom_mdl_to_fvec(cmdl,f_vec,enthalpy,mdl)

          !start up the output file
          call geom_write_trajectory(cmdl,enthalpy,output_file,mdl)

          if(pub_on_root.and.(hess_init_method.eq.hess_init_scaled)) then
             init_scalar = (dot_product(x_vec-old_x_vec,f_vec-old_f_vec)&
                  /dot_product(f_vec-old_f_vec,f_vec-old_f_vec))
             call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
                  init_method=hess_init_scaled)
          end if

          if(pub_on_root.and.update_init_on_remove) then
             old_x_vec_remove=x_vec
             old_f_vec_remove=f_vec
          end if

          if(pub_on_root.and.(hess_init_method.eq.hess_init_cellion)) then
             init_scalar(1) = (dot_product(x_vec(1:9)-old_x_vec(1:9),f_vec(1:9)-old_f_vec(1:9))/&
                  &dot_product(f_vec(1:9)-old_f_vec(1:9),f_vec(1:9)-old_f_vec(1:9)))
             init_scalar(2) = (dot_product(x_vec(10:ndim)-old_x_vec(10:ndim),f_vec(10:ndim)-old_f_vec(10:ndim))/&
                  &dot_product(f_vec(10:ndim)-old_f_vec(10:ndim),f_vec(10:ndim)-old_f_vec(10:ndim)))
             call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
                  init_method=hess_init_cellion)
          end if
       else
          if (pub_on_root) write (stdout,*)  geom_string//': Restarting previous &
               & geometry optimization.'

          !Catch the possibility that model_continuation/reuse has decided to nuke the inv_Hessian
          re_initialise = .false.
          if (pub_on_root) then
             if(cmdl%bfgs_optimisation) then
                re_initialise=all(abs(cmdl%bfgs_inv_Hessian)<tiny(0.0_dp))
             end if
             if(cmdl%lbfgs_optimisation) then
                re_initialise=(all(abs(cmdl%lbfgs_position_updates)<tiny(0.0_dp)).or.&
                     all(abs(cmdl%lbfgs_gradient_updates)<tiny(0.0_dp)))
             end if
          end if
          call comms_bcast(pub_root_proc_id,re_initialise)

       if(pub_on_root) then
          if(do_bfgs) then
             if (re_initialise) then
                call castep_cell_copy(cmdl%cell,init_cell)
                if((hess_init_method.eq.hess_init_identity).or.(hess_init_method.eq.hess_init_scaled)) then
                   call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
                        reset_updates=.true.,init_method=hess_init_identity)
                else
                   call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
                        reset_updates=.true.)
                end if
             else
                if(cmdl%lbfgs_optimisation) then
                   if((hess_init_method.eq.hess_init_identity).or.(hess_init_method.eq.hess_init_scaled)) then
                      call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
                           init_method=hess_init_identity)
                   else
                      call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl)
                   end if

                   call lbfgs_restart(toggle_on,init_cell,cmdl,mdl)
                   if (.not.allocated(tmpmm)) then
                      allocate(tmpmm(1:lbfgs_block_length*2,1:lbfgs_block_length*2),stat=ierr)
                      call utils_alloc_check('geom_BFGS','tmpmm',ierr)
                   end if
                   call linalg_invert_sym_matrix(R_mat,cmdl%lbfgs_num_updates)
                   tmpmm(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)=0.0_dp
                   tmpmm(1:cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)=-R_mat
                   tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,&
                        &1:cmdl%lbfgs_num_updates)= transpose(tmpmm(1:cmdl%lbfgs_num_updates,&
                        &cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2))
                   tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)=&
                        &matmul(matmul(tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,&
                        &1:cmdl%lbfgs_num_updates),YY_mat(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)),&
                        &tmpmm(1:cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2))
                   deallocate(R_mat,stat=ierr)
                   call utils_dealloc_check('geom_BFGS','R_mat',ierr)
                   deallocate(YY_mat,stat=ierr)
                   call utils_dealloc_check('geom_BFGS','YY_mat',ierr)
                   cmdl%bfgs_inv_Hessian=cmdl%bfgs_inv_Hessian + &
                        &matmul(matmul(cmdl%lbfgs_position_updates(:,1:cmdl%lbfgs_num_updates),&
                        &tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,1:cmdl%lbfgs_num_updates)),&
                        &transpose(cmdl%lbfgs_gradient_updates(:,1:cmdl%lbfgs_num_updates)))
                   cmdl%bfgs_inv_Hessian=cmdl%bfgs_inv_Hessian + &
                        &matmul(matmul(cmdl%lbfgs_gradient_updates(:,1:cmdl%lbfgs_num_updates),&
                        &tmpmm(1:cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)),&
                        &transpose(cmdl%lbfgs_position_updates(:,1:cmdl%lbfgs_num_updates)))
                   cmdl%bfgs_inv_Hessian=cmdl%bfgs_inv_Hessian + &
                        &matmul(matmul(cmdl%lbfgs_gradient_updates(:,1:cmdl%lbfgs_num_updates),&
                        &tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2, &
                        &cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)),&
                        &transpose(cmdl%lbfgs_gradient_updates(:,1:cmdl%lbfgs_num_updates)))
                   deallocate(tmpmm,stat=ierr)
                   call utils_dealloc_check('geom_BFGS','tmpmm',ierr)
                   deallocate(cmdl%lbfgs_position_updates,stat=ierr)
                   call utils_dealloc_check('geom_BFGS','=>cmdl%lbfgs_position_updates',ierr)
                   deallocate(cmdl%lbfgs_gradient_updates,stat=ierr)
                   call utils_dealloc_check('geom_BFGS','=>cmdl%lbfgs_gradient_updates',ierr)

                   ! ja531 -> March 2015 : We probably allocated this in lbfgs_restart
                   if(allocated(lbfgs_init_inv_hess)) then
                      deallocate(lbfgs_init_inv_hess,stat=ierr)
                      call utils_dealloc_check('geom_BFGS','lbfgs_init_inv_hess',ierr)
                   end if
                   if(allocated(lbfgs_init_hess)) then
                      deallocate(lbfgs_init_hess,stat=ierr)
                      call utils_dealloc_check('geom_BFGS','lbfgs_init_hess',ierr)
                   end if
                   if(allocated(tmpn)) then
                      deallocate(lbfgs_init_hess,stat=ierr)
                      call utils_dealloc_check('geom_BFGS','tmpn',ierr)
                   end if
                   if(allocated(D_mat)) then
                      deallocate(D_mat,stat=ierr)
                      call utils_dealloc_check('geom_BFGS','D_mat',ierr)
                   end if

                   cmdl%lbfgs_optimisation=.false.
                   cmdl%bfgs_optimisation=.true.
                end if
             end if
          else if(do_lbfgs) then
          if (re_initialise) then
             call castep_cell_copy(cmdl%cell,init_cell)
                if (.not.allocated(R_mat)) then
                   allocate(R_mat(1:lbfgs_block_length,1:lbfgs_block_length),stat=ierr)
                   call utils_alloc_check('geom_BFGS','R_mat',ierr)
                end if
                if (.not.allocated(YY_mat)) then
                   allocate(YY_mat(1:lbfgs_block_length,1:lbfgs_block_length),stat=ierr)
                   call utils_alloc_check('geom_BFGS','YY_mat',ierr)
                end if
                if (.not.allocated(tmpmm)) then
                   allocate(tmpmm(1:lbfgs_block_length*2,1:lbfgs_block_length*2),stat=ierr)
                   call utils_alloc_check('geom_BFGS','tmpmm',ierr)
                end if
                if (.not.allocated(tmp2m)) then
                   allocate(tmp2m(1:lbfgs_block_length*2),stat=ierr)
                   call utils_alloc_check('geom_BFGS','tmp2m',ierr)
                end if
                !call lbfgs_restart(toggle_on,init_cell,cmdl)
                if((hess_init_method.eq.hess_init_identity).or.(hess_init_method.eq.hess_init_scaled)) then
                   call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
                        reset_updates=.true.,init_method=hess_init_identity)
                else
                   call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
                        reset_updates=.true.)
                end if
                cmdl%lbfgs_num_updates=0
             else
                call lbfgs_restart(toggle_on,init_cell,cmdl,mdl)
             end if
          end if
          end if

       call comms_bcast(pub_root_proc_id,cmdl%bfgs_optimisation)
       call comms_bcast(pub_root_proc_id,cmdl%lbfgs_optimisation)
       call comms_bcast(pub_root_proc_id,cmdl%lbfgs_num_updates)

          !Put the strains wrt original cell and fractional atomic coords wrt current_cell into x_vec
          call geom_mdl_to_xvec(cmdl,x_vec,mdl)

          !Get stresses and forces from cmdl and put into f_vec
          call geom_mdl_to_fvec(cmdl,f_vec,enthalpy,mdl)

       end if

       call geom_precond_initialize(cmdl,mdl,is_cmplx) ! lam81
       !now prepare to relax the structure ...
       !calculate delta_vec (the update step)
       !NB Here and elsewhere we use matmul rather than BLAS for clarity,
       !   as efficiency is not an issue here (matrix sizes are either 3*3 or ndim*ndim)
       !   and electronic_minimisation always dominates.
       if (pub_on_root) then
       !       delta_vec=matmul(cmdl%bfgs_inv_Hessian,f_vec)
       delta_vec=geom_apply_vec(cmdl,f_vec,'y')
       end if
       call comms_bcast(pub_root_proc_id, delta_vec)

       !store old x_vec, f_vec and enthalpy
       old_x_vec=x_vec
       old_f_vec=f_vec
    old_delta_vec=delta_vec
       old_enthalpy=enthalpy

       !NB We reset bfgs_iteration here as the geom_max_iter refers to number of iterations
       !   THIS pass and not in total.
    if(do_bfgs) then
       cmdl%bfgs_iteration=0
       reset_iter=cmdl%bfgs_iteration
    else if(do_lbfgs) then
       cmdl%lbfgs_iteration=0
       reset_iter=cmdl%lbfgs_iteration
    end if
    previous_reset_iter=-10

    call comms_bcast(pub_root_proc_id, reset_iter)
    call comms_bcast(pub_root_proc_id, previous_reset_iter)

    if(do_bfgs)  iteration=cmdl%bfgs_iteration
    if(do_lbfgs) iteration=cmdl%lbfgs_iteration

       !initialise convergence tests (impossible to converge here with convergence windows)
       call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
            & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)
       !29/09/05: Actually it is possible to converge here now - we no longer use windows on
       !          Fmax, dRmax and Smax, and if the initial configuration passes all those but
       !          not dE(as that still has a window) then I think it is OK to call this
       !          structure converged and save some CPU time!

       !tell user where we are starting from
    if(do_bfgs) iteration=cmdl%bfgs_iteration
    if(do_lbfgs) iteration=cmdl%lbfgs_iteration
    call geom_output_converged(iteration,enthalpy,dE,Fmax,dRmax, &
         & converged_dE,converged_Fmax,converged_dRmax)

    if (present(path)) converged = neb_converged(converged, iteration, dE, Fmax, dRmax, path)

    if (pub_on_root) then
       switch_Fmax=0.10_dp*51.421_dp !'eV/Ang'
    end if
    call comms_bcast(pub_root_proc_id, switch_Fmax)


    !-------------------------------------------------------------------------!
    !BEGIN main loop over successive configurations as we minimize cell & ions!
    !-------------------------------------------------------------------------!
    main_loop: do

       !skip main_loop? i.e. already converged?
       if (converged) then
!          if (pub_on_root) write(stdout,*) "Main Loop Exit 1"
          exit main_loop
       end if

       !check to see if we have done enough iterations yet?
       !or perhaps already converged?
       !NB Test is done this way around s.t. we don't need to undo anything at end
       if (do_bfgs.and.(cmdl%bfgs_iteration>=pub_geom_max_iter)) then
!          if (pub_on_root) write(stdout,*) "Main Loop Exit 2"
          exit main_loop
       end if
       if (do_lbfgs.and.(cmdl%lbfgs_iteration>=pub_geom_max_iter)) then
!          if (pub_on_root) write(stdout,*) "Main Loop Exit 3"
          exit main_loop
       end if

       if(do_bfgs) then
          cmdl%bfgs_iteration=cmdl%bfgs_iteration+1
          iteration=cmdl%bfgs_iteration
       else if(do_lbfgs) then
          cmdl%lbfgs_iteration=cmdl%lbfgs_iteration+1
          iteration=cmdl%lbfgs_iteration
       end if

       !banner to output
       if (pub_on_root) write (stdout,4) 'Starting '//trim(geom_string)//' iteration ',iteration,' ...'
       if (present(path)) call neb_sync_path(mdl, path)

       ! kkbd: Climbing image counter for NEB runs
       if (present(path)) then
        if (path%ci_delay > 0) then
          path%ci_delay = path%ci_delay - 1
          if (pub_on_root) then
            if (path%ci_delay == 0 ) then
              write(stdout,*)"NEB: enabling climbing image"
            else if (pub_debug_on_root) then
              write(stdout,*)"NEB: neb_ci_delay reduced to ",path%ci_delay
            end if
          end if
          call services_flush
          else if (path%ci_delay == 0) then
          ! kkbd: If we're just now starting to be the climbing image
          !       OR we were last time and aren't now then reset the Hessian
          ! kkbd: TODO: this is in testing
          if (((path%climbing_image == pub_my_image).and.(prev_ci /= pub_my_image)) &
              .or. ((path%climbing_image /= prev_ci) .and. (prev_ci == pub_my_image))) then
            write(stdout,*) "NEB: Resetting inv Hessian because of a change in climbing image."
            call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl, &
                 reset_updates=.true.)
          end if

          prev_ci = path%climbing_image
        end if
      end if

       ! aam: write geometry optimiser continuation data?
       if (mod(iteration,pub_geom_backup_iter)==0) &

            call geom_opt_continuation_write(cmdl, mdl, cont_file, is_cmplx)

       ! ndmh: reset dk and NGWFs this iteration?
       ! kkbd: Reset if we've just calculated a reactant or product in NEB
       if ((pub_geom_reset_dk_ngwfs_iter>0).and.pub_geom_reuse_dk_ngwfs) then
          if (mod(iteration,pub_geom_reset_dk_ngwfs_iter)==0) then
             if (pub_read_denskern) then
                pub_read_denskern=.FALSE.
                if (pub_on_root) write(stdout,*) &
                     'Geometry Optimisation: Resetting: read_denskern &
                     &parameter set to FALSE'
             endif
             if (pub_read_tightbox_ngwfs .and. (pub_maxit_ngwf_cg>0)) then
                pub_read_tightbox_ngwfs=.FALSE.
                if (pub_on_root) write(stdout,*) &
                     'Geometry Optimisation: Resetting: read_tightbox_ngwfs &
                     &parameter set to FALSE'
             endif
             if (pub_edft) then
                if (pub_read_hamiltonian) then
                   pub_read_hamiltonian=.FALSE.
                   if (pub_on_root) write(stdout,*) &
                        'Geometry Optimisation: Resetting: read_hamiltonian &
                        &parameter set to FALSE'
                end if

                if (pub_edft_spin_fix_orig > 0) then
                   ! kkbd: Reset the fixed spin duration if we're doing a free-
                   !       spin run with a fixed-spin start
                   pub_edft_spin_fix = pub_edft_spin_fix_orig
                   if (pub_on_root) write(stdout,*) &
                           "Geometry Optimisation: Resetting: edft_spin_fix &
                           &parameter set to",pub_edft_spin_fix
                end if
             endif

             ! ebl: Restore spin-splitting if it has been turned off
             if (pub_hubbard) then
                do sp=1,mdl%regions(1)%par%num_hub_species
                   if (abs(h_species(sp)%hub_spin_splitting &
                        - h_species(sp)%orig_hub_spin_splitting) > 1E-12_DP) then
                      h_species(sp)%hub_spin_splitting = h_species(sp)%orig_hub_spin_splitting
                      if (pub_on_root) write(stdout,*) &
                           'Geometry Optimisation: Resetting: Hubbard spin splitting for ', &
                            trim(adjustl(h_species(sp)%hub_species)), ' turned back on'
                   end if
                end do
             end if

             !tell user what's about to happen
             if (pub_on_root)&
                  write (stdout,'(1x,a,i5)')  trim(geom_string)//': resetting NGWFs and &
                  &kernel at iteration', iteration

             ! re-evaluate step=0 structure with reset NGWFs/kernels
             call geom_xvec_to_mdl(x_vec,cmdl,mdl)
             ! agrecocmplx
             if(present(path)) call neb_sync_path(mdl,path)
             call geom_get_forces(cmdl,mdl,path,is_cmplx,en_sync=.true.)
             call geom_mdl_to_fvec(cmdl,f_vec,enthalpy,mdl)

             !apply inv_Hessian to get new delta after reset
             if (pub_on_root) then
                !              delta_vec=matmul(cmdl%bfgs_inv_Hessian,f_vec)
                delta_vec=geom_apply_vec(cmdl,f_vec,'y')
             end if
             call comms_bcast(pub_root_proc_id,delta_vec)

             !test convergence
             call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_none,dE,Fmax,dRmax,Smax, &
                  & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

             !are we done yet?
             if (converged) then
                ! kkbd: If we're minimizing a neb chain, only converged if all beads are converged.
                !       We should exit only if all images are converged, but we shouldn't keep going
                !       with this image's geomopt if it's converged - so cycle to the next iteration.
                !       This looks a little janky but it's necessary so we don't comms lock in neb_converged

                ! tell user what just happened
                call geom_output_converged(cmdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax,&
                    & converged_dE,converged_Fmax,converged_dRmax)

                if (present(path)) then
                  converged = neb_converged(converged, iteration, dE, Fmax, dRmax, path)
                  if (converged) then
                    exit main_loop
                  else
                    cycle main_loop
                  end if
                else
                  !tell user what just happened
                  exit main_loop
                end if
             end if

             !store old x_vec, f_vec and enthalpy
             old_x_vec=x_vec
             old_f_vec=f_vec
             old_enthalpy=enthalpy

             if (.not.pub_read_denskern.and.(.not.pub_edft)) then
                pub_read_denskern=.TRUE.
                if (pub_on_root) write(stdout,*) &
                     'Geometry Optimisation: read_denskern parameter &
                     &set to TRUE'
             endif
             if (.not.pub_read_tightbox_ngwfs .and. (pub_maxit_ngwf_cg>0)) then
                pub_read_tightbox_ngwfs=.TRUE.
                if (pub_on_root) write(stdout,*) &
                     'Geometry Optimisation: read_tightbox_ngwfs parameter &
                     &set to TRUE'
             endif
             if (pub_edft.and.(.not.pub_read_hamiltonian)) then
                pub_read_hamiltonian=.TRUE.
                if (pub_on_root) write(stdout,*) &
                     'Geometry Optimisation: pub_read_hamiltonian parameter &
                     &set to TRUE'
             endif
           end if
       end if


       !store key data in case of backtrack (including old x_vec, f_vec etc)
       call geom_BFGS_backup

       !set accept/reject flags
       done_trial=.false.
       trial_OK  =.false.
       done_line =.false.
       line_OK   =.false.
       done_quad =.false.
       quad_OK   =.false.

!       !set monotonic_E flag if it looks like things are going OK now
       if ((.not.monotonic_E).and.(iteration>reset_iter+monotonic_E_steps)) then
          if (pub_on_root) write (stdout,*) trim(geom_string)//': forbidding uphill steps from now on'
          !          if (on_root) write (stdout,*) 'BFGS: forbidding uphill steps from now on'
          monotonic_E=.true.
       end if
       if ((.not.use_enthalpy).and.(iteration>reset_iter+use_enthalpy_steps)) then
          if (pub_on_root) write (stdout,*) trim(geom_string)//': checking enthalpy rather than F.dx from now on'
          use_enthalpy=.true.
       end if
!
!       !catch case of external applied shear - we are bound to go uphill to get to equilibrium ...
!       found_shear=(any(abs(external_pressure(4:6))>epsilon(1.0_dp)))
!       call comms_bcast(pub_root_proc_id,found_shear)
!       if (found_shear.and..not.fix_all_cell) then
!          monotonic_E=.false.
!       end if

       !----------------------------------------------------------------------!
       !BEGIN trial move                                                      !
       !----------------------------------------------------------------------!

       !store F^old.delta for lambda line minimisation
       old_f_dot_delta=dot_product(f_vec,delta_vec)
       call geom_BFGS_write_summary

       !propose a trial move with lambda=1
       lambda_trial=1.0_dp

       !check trial lambda does not lead to too large a displacement wrt old structure ...
       call geom_check_lambda(lambda_trial,delta_vec,1.0_DP)


       !generate trial structure
       x_vec=old_x_vec+lambda_trial*delta_vec             !x=x+dx - includes strains & ionic coords
       call comms_bcast(pub_root_proc_id,x_vec)

       !       if (constrain_cell.and..not.fix_all_cell) then
       !          call geom_constrain_strain(cmdl,lambda_trial,old_x_vec,x_vec,old_delta_vec,delta_vec,old_f_vec,f_vec)
       !          old_f_vec=f_vec !ensure constraints are built into old as well as current f_vec for backtracking
       !       end if

       !diagnostic
       ! lk: changed output detail from > NORMAL to > VERBOSE.
       if ((pub_geom_output_detail > VERBOSE) .and. pub_on_root) then
          write (stdout,7)    trim(geom_string)//': trial lambda =',lambda_trial
          write (stdout,10)   trim(geom_string)//': trial:','i','xvec','old_xvec','delta_vec'
          do i=1,ndim
             write (stdout,1) trim(geom_string)//': trial:',i,x_vec(i),old_x_vec(i),delta_vec(i)
          end do
       end if

       !update model and current_cell for new cell & positions
       call geom_xvec_to_mdl(x_vec,cmdl,mdl)

       !tell user what's about to happen
       if (pub_on_root)&
            write (stdout,5) trim(geom_string)//': starting iteration',iteration, &
            & ' with trial guess (lambda=',lambda_trial
       if (pub_write_positions) call geom_output(cmdl,mdl)  ! info on current_cell to stdout


       !evaluate trial structure
       ! agrecocmplx
       call geom_get_forces(cmdl,mdl,path,is_cmplx)
       !get stresses and forces from cmdl and put into f_vec
       call geom_mdl_to_fvec(cmdl,f_vec,enthalpy,mdl)
       trial_f_vec=f_vec

       !store F^trial.delta for lambda line minimisation
       trial_f_dot_delta=dot_product(f_vec,delta_vec)

       !setup backtrack logic
       trial_enthalpy=enthalpy
       done_trial=.true.
       trial_OK=.false.

! Ja531/20131028: enthalpy should be a better test.
!      if (abs(trial_f_dot_delta)<=abs(old_f_dot_delta)) trial_OK=.true.
       if (use_enthalpy) then
          if (monotonic_E) then
             if(trial_enthalpy<old_enthalpy) trial_OK=.true.
          else
             trial_OK=.true.
          end if
       else
          if (abs(trial_f_dot_delta)<=abs(old_f_dot_delta)) then
             trial_OK=.true.
          else
             if(trial_enthalpy<old_enthalpy) trial_OK=.true.
          end if
       end if
       !inform user about what is going on ...
       call geom_BFGS_write_summary

       !test convergence
       call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
            & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

        !are we done yet?
        if (converged) then
          !tell user what just happened
          call geom_output_converged(iteration,enthalpy,dE,Fmax,dRmax,&
          & converged_dE,converged_Fmax,converged_dRmax)

          if (present(path)) then
            converged = neb_converged(converged, iteration, dE, Fmax, dRmax, path)
            if (converged) then
              exit main_loop
            else
              cycle main_loop
            end if
          else
            exit main_loop
          end if
        end if

       !line minimisation for best lambda is safe -> lambda_line
       !NB If the energy surface is quadratic, then lambda=1 is optimal, but
       !   if we start some way from the minimum this is unlikely, so we do
       !   a line minimization instead.
       !BUT for efficiency we keep the trial structure if lambda_line~1
       !   (unless it has gone uphill) to prevent needless calls to geom_get_forces!!!
       if (abs(trial_f_dot_delta-old_f_dot_delta)<epsilon(old_f_dot_delta)) then
          if (pub_on_root) then
             write (stdout,*) trim(geom_string)//': Looks like this system is as converged as possible.'
             write (stdout,*) '      Maybe your geometry convergence tolerances are too tight?'
          end if

          !test convergence
          call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_none,dE,Fmax,dRmax,Smax, &
               & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)
          !spoof the main flag
          converged=.true.

          if (present(path)) then
            !tell user what just happened
            call geom_output_converged(iteration,enthalpy,dE,Fmax,dRmax,&
            & converged_dE,converged_Fmax,converged_dRmax)

            converged = neb_converged(converged, iteration, dE, Fmax, dRmax, path)
            if (converged) then
              exit main_loop
            else
              cycle main_loop
            end if
          else
            exit main_loop
          end if

       end if

       !so we know that abs(df_dot_delta)>epsilon ...
       df_dot_delta=trial_f_dot_delta-old_f_dot_delta
       lambda_line=-lambda_trial*old_f_dot_delta/df_dot_delta

       !parallel synch
       !        call comms_gcopy(lambda_line,1)

       !check line minimisation does not lead to too large a displacement wrt old structure ...
       if ((1.0_dp-lambda_trial)>epsilon(lambda_trial).and.lambda_line>lambda_trial) then
          !lambda_trial has been truncated and line>trial => line will be truncated too => nothing will change
          !so this is a time NOT to check the proposed new lambda and let a big displacment go through
          if (pub_geom_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
               &trim(geom_string)//': trial lambda was truncated for safety so NOT truncating new&
             & lambda'
          call geom_check_lambda(lambda_line,delta_vec,1.5_DP)
       else
          call geom_check_lambda(lambda_line,delta_vec,1.0_DP)
       end if

          !do some simple error checking on lambda
          !       lambda_bad=(monotonic_E.and.fixed_npw.and.(lambda_line<=0.0_dp))
          lambda_bad=(monotonic_E.and.use_enthalpy.and.(lambda_line<=0.0_dp))
          if (lambda_bad) then
             !best line lambda<0 => stepping backwards => not good!
             if (trial_OK) then
                !backtrack to trial structure as that was OK
                ! lk: changed output detail from > NORMAL to > VERBOSE.
                if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
                   write (stdout,*) trim(geom_string)//': Warning - lambda < 0'
                   write (stdout,*) trim(geom_string)//': Warning - reset to trial'
                end if
                lambda_line=lambda_trial
             else
                !backtrack to old structure as that was OK
                ! lk: changed output detail from > NORMAL to > VERBOSE.
                if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
                   write (stdout,*) trim(geom_string)//': Warning - lambda < 0'
                   write (stdout,*) trim(geom_string)//':         - reverting to earlier configuration and reseting inverse Hessian'
                end if

                abort_BFGS=.false.
                call geom_BFGS_revert_old !calculate new value for abort_BFGS
                if (abort_BFGS) then
!                   if (pub_on_root) write(stdout,*) "Main Loop Exit 7"
                    ! kkbd: TODO: Work out if we can just cycle safely here
                    if (present(path)) call utils_abort("NEB: BFGS Aborted.")
                   exit main_loop     !no point continuing
                else
!                   if (pub_on_root) write(stdout,*) "Main Loop Cycle 1"
                   ! kkbd: If we're minimizing a neb chain, consider the other images
                   if (present(path)) converged = neb_converged(.false., iteration, 0.0_dp, 0.0_dp, 0.0_dp, path)
                   cycle main_loop    !successful backtrack, do fresh iteration
                end if

             end if                   !trial_OK=.false.
          end if                      !lambda_line<0

       !best line lambda not close enough to trial lambda (Pfrommer tolerance)
       !or trial structure has gone uphill from end of previous iteration and lambda_line<>lambda_trial
       ! => evaluate new trial structure
          do_line=.false.
       lambda_diff=abs(lambda_line-lambda_trial)>epsilon(lambda_line)
             do_line=   (abs(lambda_line-lambda_trial)>geom_linmin_tol)  &
                  & .or.((monotonic_E.and..not.trial_OK).and.lambda_diff) &
                  & .or.((.not.monotonic_E.and.(trial_enthalpy-fudge)>old_enthalpy).and.lambda_diff)
             if (do_line) then

                !-------------------------------------------------------------------!
                !BEGIN line-minimisation move                                       !
                !-------------------------------------------------------------------!

                !store key data in case of backtrack (including trial x_vec, f_vec etc)
                if (trial_OK) call geom_BFGS_backup

                !setup new state of system
                x_vec=old_x_vec+lambda_line*delta_vec

                !diagnostic
                ! lk: changed output detail from > NORMAL to > VERBOSE.
                if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
                   write (stdout,7)    trim(geom_string)//': line minimisation lambda=',lambda_line
                   write (stdout,10)   trim(geom_string)//': line :','i','xvec','old_xvec','delta_vec'
                   do i=1,ndim
                      write (stdout,1) trim(geom_string)//': line :',i,x_vec(i),old_x_vec(i),delta_vec(i)
                   end do
                end if

                !update model and current_cell for new cell & positions
                call geom_xvec_to_mdl(x_vec,cmdl,mdl)

                !tell user what's about to happen
                if (pub_on_root) write (stdout,5) trim(geom_string)//': improving iteration',iteration, &
                     & ' with line minimization (lambda=',lambda_line
                if (pub_write_positions) call geom_output(cmdl,mdl) !info on current_cell to stdout

                !evaluate line minimisation structure
                ! agrecocmplx
                call geom_get_forces(cmdl,mdl,path,is_cmplx)

                !get stresses and forces from cmdl and put into f_vec
                call geom_mdl_to_fvec(cmdl,f_vec,enthalpy,mdl)
                line_f_vec=f_vec

                !store F^line.delta for quad minimisation
                line_f_dot_delta=dot_product(f_vec,delta_vec)

                !setup backtrack logic
                line_enthalpy=enthalpy
                done_line=.true.
                line_OK=.false.

                ! Ja531/20131028: enthalpy should be a better test.
                !             if (abs(line_f_dot_delta)<=abs(trial_f_dot_delta)) then
                if(use_enthalpy) then
                   if (monotonic_E) then
                      if (trial_OK) then
                         if (line_enthalpy<trial_enthalpy) line_OK=.true.
                      else
                         if (line_enthalpy<old_enthalpy) line_OK=.true.
                      end if
                   else
                      if (line_enthalpy<trial_enthalpy) then
                         line_OK=.true.
                         trial_OK =.false.
                      else
                         trial_OK=.true.
                         line_OK=.false.
                      end if
                   end if
                else
                   if (abs(line_f_dot_delta)<=abs(trial_f_dot_delta)) then
                      trial_OK=.false.
                      line_OK=.true.
                   else
                      if (line_enthalpy<trial_enthalpy) then
                         line_OK=.true.
                         trial_OK =.false.
                      end if
                   end if
                end if


          !inform user about what is going on ...
          call geom_BFGS_write_summary

          !test convergence
          call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_none,dE,Fmax,dRmax,Smax, &
               & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

          !are we done yet?
          if (converged) then
            !tell user what just happened
            call geom_output_converged(iteration,enthalpy,dE,Fmax,dRmax,&
            & converged_dE,converged_Fmax,converged_dRmax)

            if (present(path)) then
              converged = neb_converged(converged, iteration, dE, Fmax, dRmax, path)
              if (converged) then
                exit main_loop
              else
                cycle main_loop
              end if
            else
              exit main_loop
            end if
          end if

          !Is quadratic interpolation required?
          !We actually do this using simple Newton-Raphson and assuming quadratic form for (F.delta) vs lambda
          !This rarely makes a big difference, except when the line minimisation step has failed to go downhill
          ! whereupon it is worth doing everything possible to avoid backtracking to old structure
          ! and resetting the inv_Hessian.
          !Thats also why I've added the monotonic_E flag - sometimes going uphill gets around a barrier
          !if all else has failed!

          !we know that lambda_trial-lambda_line>epsilon and lambda_line<>0 from if blocks above
          df_dot_delta_dlambda=(trial_f_dot_delta-line_f_dot_delta)/(lambda_trial-lambda_line)
          if (lambda_line<lambda_trial) then
             !can do two-sided estimate (averaging gradients) to get a best value
             df_dot_delta_dlambda=0.5_dp*(df_dot_delta_dlambda+(line_f_dot_delta-old_f_dot_delta)/lambda_line)
          else
             !cannot do any more as one sided estimate is all that is possible I'm afraid!
          end if
          if (abs(df_dot_delta_dlambda)>tiny(df_dot_delta_dlambda)) then
             lambda_quad=lambda_line-line_f_dot_delta/df_dot_delta_dlambda
          else
             !this implies that the function is linear and so lambda_line must be the best answer!
             lambda_quad=lambda_line
          end if

          !check quadratic lambda does not lead to too large a displacement wrt old structure ...
             rescale_lambda=((1.0_dp-lambda_trial)>epsilon(lambda_trial)).and.(lambda_quad>lambda_trial)
             ! lk: changed output detail from > NORMAL to > VERBOSE.
             if (rescale_lambda) then
                !lambda_trial has been truncated and quad>trial => quad will be truncated too
                !=> nothing will change unless we rescale the allowed displacement limit
                if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                     & trim(geom_string)//': trial lambda was truncated so increasing allowed displacement this time'
                call geom_check_lambda(lambda_quad,delta_vec,1.0_DP)
             else
             !compare actual to expected value of lamba_line (not saved)
                rescale_lambda=((lambda_line+lambda_trial*old_f_dot_delta/(trial_f_dot_delta-old_f_dot_delta))&
                     &<epsilon(lambda_line) .and.lambda_quad>lambda_line)
                if (rescale_lambda) then
                   !lambda_line has been truncated and quad>line => quad will be truncated too
                   !=> nothing will change unless we rescale the allowed displacement limit
                   if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                        & trim(geom_string)//': line lambda was truncated so increasing allowed displacement this time'
                   call geom_check_lambda(lambda_quad,delta_vec,1.0_DP) !scale=2.0 for quad
                else
                   !all is OK so do normal check
                   call geom_check_lambda(lambda_quad,delta_vec,1.0_DP) !scale=1.0 normally
                end if
             end if

             !do some simple error checking on lambda_quad
             if (lambda_bad) then
                if (line_OK) then
                   if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
                      write (stdout,*) trim(geom_string)//': Warning - lambda < 0'
                      write (stdout,*) trim(geom_string)//': Warning - reset to line minimisation'
                   end if
                   lambda_quad=lambda_line
                else
                   if (trial_OK) then
                      if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
                         write (stdout,*) trim(geom_string)//': Warning - lambda < 0'
                         write (stdout,*) trim(geom_string)//': Warning - reset to trial'
                      end if
                      lambda_quad=lambda_trial
             else
                      if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
                         write (stdout,*) trim(geom_string)//': Warning - lambda < 0'
                         write (stdout,*) trim(geom_string)//&
                              &':         - reverting to earlier configuration and reseting inverse Hessian'
             end if

                      abort_BFGS=.false.
                      call geom_BFGS_revert_old !calculate new value for abort_BFGS
                      if (abort_BFGS) then
!                         if (pub_on_root) write(stdout,*) "Main Loop Exit 9"
                         ! kkbd: TODO: Work out if we can just cycle safely here
                         if (present(path)) call utils_abort("NEB: BFGS Aborted.")
                         exit main_loop     !no point continuing
                      else
!                         if (pub_on_root) write(stdout,*) "Main Loop Cycle 2"
                         ! kkbd: If minimizing a neb chain, consider the other images
                         if (present(path)) converged = neb_converged(.false., iteration, 0.0_dp, 0.0_dp, 0.0_dp, path)
                         cycle main_loop    !successful backtrack, do fresh iteration
          end if

                   end if                   !trial_OK=.false.
                end if                      !line_OK=.false.
             end if                         !lambda_quad<0

             !best quadratic lambda not close enough to line minimisation lambda (use geom_linmin_tol as fractional)
          !or line minimisation structure has gone uphill from trial structure and lamba_quad<>lambda_line
             do_quad=.false.
          lambda_diff=abs(lambda_quad-lambda_line)>epsilon(lambda_quad)
                do_quad = ((monotonic_E.and..not.line_OK) .and. &
                     &     (abs((lambda_quad-lambda_line)/lambda_line)>geom_linmin_tol)) &
                     &.or.((monotonic_E.and..not.trial_OK.and..not.line_OK).and.lambda_diff) &
                     &.or.((.not.monotonic_E.and.(trial_enthalpy-fudge)>old_enthalpy.and. &
                     &     (line_enthalpy-fudge)>old_enthalpy).and.lambda_diff)

             if (do_quad) then

             !----------------------------------------------------------------!
             !BEGIN quadratic-minimisation move                               !
             !----------------------------------------------------------------!

             !store key data in case of backtrack (including line minimisation x_vec, f_vec etc)
             if (line_OK) call geom_BFGS_backup

             !set up new state of system
             x_vec=old_x_vec+lambda_quad*delta_vec

             !diagnostic
             ! lk: changed output detail from > NORMAL to > VERBOSE.
             if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
                write (stdout,7)    trim(geom_string)//': quad lambda =',lambda_quad
                write (stdout,10)   trim(geom_string)//': quad :','i','xvec','old_xvec','delta_vec'
                do i=1,ndim
                   write (stdout,1) trim(geom_string)//': quad :',i,x_vec(i),old_x_vec(i),delta_vec(i)
                end do
             end if

             !update model and current_cell for new cell & positions
             call geom_xvec_to_mdl(x_vec,cmdl,mdl)

             !tell user what's about to happen
             if (pub_on_root)write(stdout,5) trim(geom_string)//': improving iteration',iteration, &
                  & ' with quad minimization (lambda=',lambda_quad
             if (pub_write_positions) call geom_output(cmdl,mdl)              !info on current_cell to stdout

             !evaluate quad structure
             ! agrecocmplx
             call geom_get_forces(cmdl,mdl,path,is_cmplx)

             !get stresses and forces from cmdl and put into f_vec
             call geom_mdl_to_fvec(cmdl,f_vec,enthalpy,mdl)
             quad_f_vec=f_vec

             !store F^quad.delta for neat output
             quad_f_dot_delta=dot_product(f_vec,delta_vec)

             !setup backtrack logic
             quad_enthalpy=enthalpy
             done_quad=.true.
             quad_OK=.false.

! Ja531/20131028: enthalpy should be a better test.
!                if (abs(quad_f_dot_delta)<=abs(line_f_dot_delta)) then
             if(use_enthalpy) then
                if (monotonic_E) then
                   if (line_OK) then
                      if (quad_enthalpy<line_enthalpy) quad_OK=.true.
                   else
                      if (trial_OK) then
                         if (quad_enthalpy<trial_enthalpy) quad_OK=.true.
                      else
                         if (quad_enthalpy<old_enthalpy) quad_OK=.true.
                      end if
                   end if
                else
                   !still only want one flag true, rest false! choose the lowest of the attempts made
                   if ((quad_enthalpy<line_enthalpy).and.(quad_enthalpy<trial_enthalpy)) then
                      quad_OK=.true.
                      line_OK=.false.
                      trial_OK=.false.
                   end if
                end if
             else
                if (abs(quad_f_dot_delta)<=abs(line_f_dot_delta)) then
                   trial_OK=.false.
                   line_OK=.false.
                   quad_OK=.true.
                else
                   if ((quad_enthalpy<line_enthalpy).and.(quad_enthalpy<trial_enthalpy)) then
                      quad_OK=.true.
                      line_OK=.false.
                      trial_OK=.false.
                   end if
                end if
             end if


             !inform user about what is going on ...
             call geom_BFGS_write_summary

             !test convergence
             call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_none,dE,Fmax,dRmax,Smax, &
                  & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

             !are we done yet?
             if (converged) then
                !tell user what just happened
                call geom_output_converged(iteration,enthalpy,dE,Fmax,dRmax,&
                & converged_dE,converged_Fmax,converged_dRmax)

                if (present(path)) then
                  converged = neb_converged(converged, iteration, dE, Fmax, dRmax, path)
                  if (converged) then
                    exit main_loop
                  else
                    cycle main_loop
                  end if
                else
                  exit main_loop
                end if
              end if

             !----------------------------------------------------------------!
             !END quadratic-minimisation move                                 !
             !----------------------------------------------------------------!

          else !lambda_quad close enough to lambda_line so no need to do any more
             if (pub_geom_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                  &trim(geom_string)//': line minimization close enough - no need to evaluate &
                &quadratic step'

             !----------------------------------------------------------------!
             !END line-minimisation move                                      !
             !----------------------------------------------------------------!

          end if

       else !lambda_line close enough to trial (lambda=1) so no need to change anything
          if (pub_geom_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
               &trim(geom_string)//': trial structure close enough - no need to evaluate line step'

          !-------------------------------------------------------------------!
          !END trial move                                                     !
          !-------------------------------------------------------------------!

       end if


       !----------------------------------------------------------------------!
       !BEGIN backtrack code                                                  !
       !----------------------------------------------------------------------!
       ! lk: changed output detail from > NORMAL to > VERBOSE.
       !check if we need to backtrack the proposed structure
       lambda_backtrack=-999.9_dp    !set to impossible value
       revert_old =.false.
       if (done_quad) then           !we may or may not have done a quad step, so check first
          if (quad_OK) then
             if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                  &trim(geom_string)//': accepting quad step - no need to backtrack'
          else                       !quad   step has gone uphill so reject and revert to line step
             if (line_OK) then       !line step went downhill so accept line step
                if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                     &trim(geom_string)//': backtracking from quad to line'
                lambda_backtrack=lambda_line
             else                    !line step has gone uphill so reject and revert to trial step
                if (trial_OK) then   !trial  step went downhill so accept trial step
                   if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                        &trim(geom_string)//': backtracking from quad to trial'
                   lambda_backtrack=lambda_trial
                else                 !trial  step has gone uphill so reject and revert to old structure
                   if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                        &trim(geom_string)//': backtracking from quad to old'
                   lambda_backtrack=0.0_dp
                   revert_old=.true.
                end if
             end if
          end if
       else if (done_line) then    !we may or may not have done a line step, so check next
          if (line_OK) then
             if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                  &trim(geom_string)//': accepting line step - no need to backtrack'
          else                       !line step has gone uphill so reject and revert to trial step
             if (trial_OK) then      !trial  step went downhill so accept trial step
                if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                     &trim(geom_string)//': backtracking from line to trial'
                lambda_backtrack=lambda_trial
             else                    !trial  step has gone uphill so reject and revert to old structure
                if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                     &trim(geom_string)//': backtracking from line to old'
                lambda_backtrack=0.0_dp
                revert_old=.true.
             end if
          end if
       else                          !we always do a trial step, so check last.
          if (trial_OK) then
             if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                  &trim(geom_string)//': accepting trial step - no need to backtrack'
          else                       !trial  step has gone uphill so reject and revert to old structure
             if (pub_geom_output_detail>VERBOSE.and.pub_on_root) write (stdout,*) &
                  &trim(geom_string)//': backtracking from trial to old'
             lambda_backtrack=0.0_dp
             revert_old=.true.
          end if
       end if

       !warn user we are about to backtrack
       if (revert_old) then

          if (pub_on_root) write (stdout,*) trim(geom_string)//': reverting to earlier &
               &configuration and reseting inverse Hessian'
          abort_BFGS=.false.
          call geom_BFGS_revert_old !calculate new value for abort_BFGS
          if (abort_BFGS) then
!             if (pub_on_root) write(stdout,*) "Main Loop Exit 11"
             ! kkbd: TODO: Work out if we can just cycle safely here
             if (present(path)) call utils_abort("NEB: BFGS Aborted.")
             exit main_loop     !no point continuing
          else
             ! kkbd: if we're doing a NEB run we need to consider the other images.
             if (present(path)) then
                converged = neb_converged(.false., iteration, 0.0_dp, 0.0_dp, &
                     0.0_dp, path)
             end if
!             if (pub_on_root) write(stdout,*) "Main Loop Cycle 3"
             cycle main_loop    !successful backtrack, do fresh iteration
          end if

       else                     !not revert_old

          if (lambda_backtrack>-999.0_dp) then !ie we have assigned a definite value in check-logic above

             x_vec=old_x_vec+lambda_backtrack*delta_vec
             call geom_xvec_to_mdl(x_vec,cmdl,mdl)

             !tell user what's about to happen
             if (pub_on_root) write (stdout,*) trim(geom_string)//': reverting to earlier configuration'
             if (pub_write_positions) call geom_output(cmdl,mdl)              !info on current_cell to stdout

             !NB Presume above logic is consistent with geom_BFGS_backup logic s.t. we get the right restore ...
             call geom_BFGS_restore

          end if                !lambda_backtrack>0

       end if                   !revert_old

       !----------------------------------------------------------------------!
       !END backtrack code                                                    !
       !----------------------------------------------------------------------!

       !----------------------------------------------------------------------!
       !BEGIN finish off this iteration                                       !
       !----------------------------------------------------------------------!

       !update inv_Hessian using BFGS
       call geom_update_inv_Hessian(old_x_vec,x_vec,old_f_vec,f_vec,init_cell,cmdl,mdl)

       ! ja531 --> bugfix 2011/09/18 10:43 : The following variables need to be distributed, as they are
       !           set within geom_update_inv_Hessian on the root proc.
       call comms_bcast(pub_root_proc_id,reset_iter)
       !apply inv_Hessian to get new delta
       if (pub_on_root) then
          !          delta_vec=matmul(cmdl%bfgs_inv_Hessian,f_vec)
          delta_vec=geom_apply_vec(cmdl,f_vec,'y')
       end if
       call comms_bcast(pub_root_proc_id,delta_vec)

       !test convergence
       call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_none,dE,Fmax,dRmax,Smax, &
            & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

       !tell user what just happened
       call geom_output_converged(iteration,enthalpy,dE,Fmax,dRmax,&
            & converged_dE,converged_Fmax,converged_dRmax)

       !are we done yet?
       if (present(path)) converged = neb_converged(converged, iteration, dE, &
            Fmax, dRmax, path)

       if (converged) then
         exit main_loop
       end if

       !update output file
       call geom_write_trajectory(cmdl,enthalpy,output_file,mdl)

       !store old x_vec, f_vec and enthalpy
       old_x_vec=x_vec
       old_f_vec=f_vec
       old_delta_vec=delta_vec
       old_enthalpy=enthalpy


       !----------------------------------------------------------------------!
       !END finish off this iteration                                         !
       !----------------------------------------------------------------------!

      end do main_loop
      !-------------------------------------------------------------------------!
      !END main loop over successive configurations as we minimize cell & ions  !
      !-------------------------------------------------------------------------!

      !-------------------------------------------------------------------------!
      !BEGIN final output and tidy up                                           !
      !-------------------------------------------------------------------------!

      !write final entry to the output file
      if (converged) call geom_write_trajectory(cmdl,enthalpy,output_file,mdl)


      !Actually, we need to UNSET the found_forces flag if we have any ionic
      !constraints s.t. check_forces_stresses in main will recalculate the
      !unconstrained forces and so produce non-zero outputs with firstd_output_forces
      !and then reset the found_forces flag before the call to model_write in main.

      if (pub_on_root) then
       if (converged) then
          write(stdout,*) geom_string//': Geometry optimization completed successfully.'
       else
          write(stdout,*) geom_string//': Geometry optimization failed to converge after ', &
               & iteration,' steps'
       end if
       write (stdout,11) geom_string//': Final Configuration:'
       call geom_output(cmdl,mdl)
       write (stdout,2) geom_string//': Final Enthalpy     =',enthalpy,trim(energy_label)
      end if

    if (iteration-reset_iter>1) then
       !       if(house_init.and.do_lbfgs) then
       !          call lbfgs_update_bandmatrix(cmdl,update_vecs,lbfgs_recurs_mat,lbfgs_init_inv_hess)
       !       else
       if(hess_init_method.lt.3) then
          if(pub_on_root) call geom_update_inv_Hessian_params(init_cell,cmdl,mdl,final=.true.)
       end if
    end if
      if (pub_on_root) then
         if (.not.fix_all_ions) then
            if (abs((new_frequency_est-old_frequency_est)/old_frequency_est)>0.01_dp) then !warn if 1% change
             write (stdout,2) geom_string//': Final <frequency>  =',new_frequency_est,trim(frequency_label)
            end if
         end if
      end if

    ! aam: store as cont_file whether converged or not
    call geom_opt_continuation_write(cmdl,mdl,cont_file,is_cmplx)

    !clean up at the end
    call castep_cell_dealloc(backup_cell)

    ! ===============================================================================
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<< Remove backup files >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ===============================================================================
    if (pub_on_root) then
       ! NGWF backup file first....
       ! jcap: loop over regions
       do isub=1,mdl%nsub
          ! Find next available unit specifier
          backup_ngwf_unit(isub) = utils_unit()
          ! tell user what's going on
          write(stdout,'(/3a)',advance='no') &
               ' Removing NGWF backup file "',trim(backup_ngwf_file(isub)),'" ... '
          ! open file
          open(unit=backup_ngwf_unit(isub),file=backup_ngwf_file(isub),iostat=io_status,&
               form='UNFORMATTED',action='READWRITE')
          call utils_open_unit_check('geom_BFGS','backup_ngwf_unit',io_status)

          ! delete file
          close(unit=backup_ngwf_unit(isub),iostat=io_status,status='DELETE')
          call utils_close_unit_check('geom_BFGS','backup_ngwf_unit',io_status)

          ! success
          write(stdout,'(a)') 'done'
       end do

       ! ... then density kernel backup file...
       ! Find next available unit specifier
       backup_dk_unit = utils_unit()
       ! tell user what's going on
       write(stdout,'(3a)',advance='no') &
            ' Removing density kernel backup file "',trim(backup_dk_file),'" ... '
       ! open file
       open(unit=backup_dk_unit,file=backup_dk_file,iostat=io_status,&
            form='UNFORMATTED',action='READWRITE')
       call utils_open_unit_check('geom_BFGS','backup_dk_unit',io_status)

       ! delete file
       close(unit=backup_dk_unit,iostat=io_status,status='DELETE')
       call utils_close_unit_check('geom_BFGS','backup_dk_unit', io_status)

       ! success
       write(stdout,'(a)') 'done'
    endif
    ! =================================================================================
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<< End remove backup files >>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! =================================================================================


    deallocate (trial_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','trial_f_vec',ierr)

    deallocate (line_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','line_f_vec',ierr)

    deallocate (quad_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','quad_f_vec',ierr)

    deallocate(backup_forces,stat=ierr)
    call utils_dealloc_check('geom_BFGS','backup_forces',ierr)

    deallocate (backup_x_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','backup_x_vec',ierr)

    deallocate (backup_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','backup_f_vec',ierr)

    call castep_cell_dealloc(init_cell)

    deallocate (x_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','x_vec',ierr)

    deallocate (old_x_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','old_x_vec',ierr)

    deallocate (f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','f_vec',ierr)

    deallocate (old_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','old_f_vec',ierr)

    deallocate (delta_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','delta_vec',ierr)


    ! ja531 -> March 2015: l-BFGS Cleanup
    if(allocated(lbfgs_init_inv_hess)) then
       deallocate(lbfgs_init_inv_hess,stat=ierr)
       call utils_dealloc_check('geom_BFGS','lbfgs_init_inv_hess',ierr)
    end if
    if(allocated(lbfgs_init_hess)) then
       deallocate(lbfgs_init_hess,stat=ierr)
       call utils_dealloc_check('geom_BFGS','lbfgs_init_hess',ierr)
    end if
    if(allocated(tmpn)) then
       deallocate(tmpn,stat=ierr)
       call utils_dealloc_check('geom_BFGS','tmpn',ierr)
    end if
    if(allocated(D_mat)) then
       deallocate(D_mat,stat=ierr)
       call utils_dealloc_check('geom_BFGS','D_mat',ierr)
    end if
    if(allocated(R_mat)) then
       deallocate(R_mat,stat=ierr)
       call utils_dealloc_check('geom_BFGS','R_mat',ierr)
    end if
    if(allocated(YY_mat)) then
       deallocate(YY_mat,stat=ierr)
       call utils_dealloc_check('geom_BFGS','YY_mat',ierr)
    end if
    if(allocated(tmpmm)) then
       deallocate(tmpmm,stat=ierr)
       call utils_dealloc_check('geom_BFGS','tmpmm',ierr)
    end if
    if(allocated(tmp2m)) then
       deallocate(tmp2m,stat=ierr)
       call utils_dealloc_check('geom_BFGS','tmp2m',ierr)
    end if
    if(allocated(tmpmm)) then
       deallocate(tmpmm,stat=ierr)
       call utils_dealloc_check('geom_BFGS','tmpmm',ierr)
    end if
    if(allocated(tmp2m)) then
       deallocate(tmp2m,stat=ierr)
       call utils_dealloc_check('geom_BFGS','tmp2m',ierr)
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_BFGS'

    !-------------------------------------------------------------------------!
    !Its the end of the BFGS code as we know it and I feel fine ...           !
    !-------------------------------------------------------------------------!

1   format(1x,a,i4,4f15.5)
2   format(1x,a,es17.8e3,1x,a)
4   format(/,80('='),/,1x,a,i4,a,/,80('='))
5   format(/,80('-'),/,1x,a,i4,a,f12.6,')',/,80('-'))
7   format(1x,a,f11.6)
10  format(1x,a,a4,4a15)
11  format(/,80('='),/,1x,a,/,80('='))

    return

  contains

    subroutine geom_BFGS_revert_old
      !=========================================================================!
      ! 'contain'ed subroutine to implement backtracking from current to        !
      ! previous configuration as done in several places and in same way        !
      !-------------------------------------------------------------------------!
      ! Arguments:                                                              !
      !   none - access all geom_BFGS variables directly as in same scope       !
      !-------------------------------------------------------------------------!
      ! Parent module variables used:                                           !
      !   same as geom_BFGS as in same scope                                    !
      !-------------------------------------------------------------------------!
      ! Modules used:                                                           !
      !   same as geom_BFGS as in same scope                                    !
      !-------------------------------------------------------------------------!
      ! Key Internal Variables:                                                 !
      !   reset_iter = value of cmdl%bfgs_iteration last time we did this        !
      !                so can trap consecutive backtracks and abort if required !
      !   monotonic_E= can turn off if we get desperate                         !
      !   abort_BFGS = give up altogether                                       !
      !-------------------------------------------------------------------------!
      ! Necessary conditions:                                                   !
      !-------------------------------------------------------------------------!
      ! Written by Matt Probert, v1.0, 11/05/2002                               !
      !-------------------------------------------------------------------------!
      ! Modified for ONETEP by Arash A Mostofi, 2004                            !
      !=========================================================================!


      use constants, only: stdout
      use rundat, only: pub_task
      implicit none

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering geom_BFGS_revert_old'

      !backtrack model
      x_vec=old_x_vec
      f_vec=old_f_vec
      enthalpy=old_enthalpy
      call geom_xvec_to_mdl(old_x_vec,cmdl,mdl)
      call geom_fvec_to_mdl(old_f_vec,cmdl,mdl)

      !NB geom_xvec_to_mdl calls model_cell_changed
      !=> invalidates model flags so we need to reset them here
      !=> but we have not restored the wavefunction etc
!      if (.not.fix_all_ions) cmdl%found_forces = .true.
!      if (.not.fix_all_cell) cmdl%found_stress = .true.
      if(do_bfgs) then
         cmdl%bfgs_optimisation=.true.
      else if(do_lbfgs) then
         cmdl%lbfgs_optimisation=.true.
      end if
      !before we throw away the old inv_Hessian, lets try to extract some useful stuff out of it
      if(pub_on_root) then
         if(hess_init_method.lt.3) then
            if (iteration-reset_iter>1) call geom_update_inv_Hessian_params(init_cell,cmdl,mdl)
         end if
      end if


      call castep_cell_copy(cmdl%cell,init_cell)

      if(pub_on_root) then
         call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
              reset_updates=.true.,init_method=hess_init_method)
      end if



      !apply inv_Hessian to get new delta
      if (pub_on_root) then
         !         delta_vec=matmul(cmdl%bfgs_inv_Hessian,f_vec)
         delta_vec=geom_apply_vec(cmdl,f_vec,'y')
      end if
      call comms_bcast(pub_root_proc_id,delta_vec)

      !test convergence
      call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_backwards,dE,Fmax,dRmax,Smax, &
           & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

      !tell user what just happened
      call geom_output_converged(iteration,enthalpy,dE,Fmax,dRmax,&
           converged_dE,converged_Fmax,converged_dRmax)

      !abort or continue?
      if (iteration-reset_iter==1) then !whoops - we've reseting after a previous reset!
         !Test if we are able to toggle off-diag element in inv_Hessian on/off ..

            !We're in trouble ...
            if (pub_on_root) then
               write (stdout,*) trim(geom_string)//': Warning - Repeated consecutive reset of inverse Hessian'
               write (stdout,*) trim(geom_string)//':           without satisfying convergence criteria which'
               write (stdout,*) trim(geom_string)//':           looks like '//trim(geom_string)// &
                    & ' has run out of search directions.'
            end if
            if (monotonic_E) then
               if (pub_on_root) then
                  write (stdout,*) trim(geom_string)//': Warning - Lets try allowing some uphill steps and see if'
                  write (stdout,*) trim(geom_string)//':           we can get around this barrier.'
                  write (stdout,*) trim(geom_string)//': Warning - It is possible that the system may now converge to'
                  write (stdout,*) trim(geom_string)//':           a stationary point OTHER than the desired minimum.'
                  write (stdout,*) trim(geom_string)//': Hint    - this may be an indication that either:'
                  write (stdout,*) trim(geom_string)//':           a) you are using a poor guess at geom_frequency_est'
                  write (stdout,*) trim(geom_string)//':           and/or geom_modulus_est, or'
                  write (stdout,*) trim(geom_string)//':           b) you are using unrealistic convergence criteria.'
                  write (stdout,*) trim(geom_string)//':           Suggest therefore that you consider changing them!'
               end if
               monotonic_E=.false.
               use_enthalpy=.false.
               !rather than simply repeat the last step, lets take the best of the (trial/line/quad) steps
               !and use that as our starting point

               if (done_quad) then
                  quad_OK=((quad_enthalpy<line_enthalpy).and.(quad_enthalpy<trial_enthalpy))
               else
                  quad_OK=.false.
               end if
               call comms_bcast(pub_root_proc_id,quad_OK)
               if (quad_OK) then
                  !revert to quad step
                  lambda_backtrack=lambda_quad
                  enthalpy=quad_enthalpy
                  f_vec=quad_f_vec
               else
                  if (done_line) then
                     line_OK=(line_enthalpy<trial_enthalpy)
                  else
                     line_OK=.false.
                  end if
                  call comms_bcast(pub_root_proc_id,line_OK)
                  if (line_OK) then
                     !revert to line step
                     lambda_backtrack=lambda_line
                     enthalpy=line_enthalpy
                     f_vec=line_f_vec
                  else
                     !revert to trial step
                     lambda_backtrack=lambda_trial
                     enthalpy=trial_enthalpy
                     f_vec=trial_f_vec
                  end if
            end if

               !restore the structure
               x_vec=old_x_vec+lambda_backtrack*delta_vec
               call comms_bcast(pub_root_proc_id,x_vec)
               call geom_xvec_to_mdl(x_vec,cmdl,mdl)
               call geom_fvec_to_mdl(f_vec,cmdl,mdl)
               old_x_vec=x_vec
               old_f_vec=f_vec
               old_enthalpy=enthalpy

               !we don't bother with the wavefunction etc as it has not been backed up
               !and it is probably pretty close in all directions anyway
               if(do_bfgs) then
                  cmdl%bfgs_optimisation=.true.
               else if(do_lbfgs) then
                  cmdl%lbfgs_optimisation=.true.
               end if
            else                      !not monotonic_E
               iwarn=iwarn+1
               if (pub_task == "TRANSITIONSTATESEARCH") then
                  ! kkbd: don't abort if one image happens to be well optimized
                  if (pub_on_root) write (stdout,*) "NEB: Bypassing "//trim(geom_string)//" termination."
               else
                  if (pub_on_root) write (stdout,*) trim(geom_string)//': Warning - Terminating '//trim(geom_string)//' loop.'
                  abort_BFGS=.true.
               end if
            end if                    !monotonic_E
      else                            !not consecutive reset
         previous_reset_iter=reset_iter
         reset_iter=iteration
         call comms_bcast(pub_root_proc_id,reset_iter)

      end if

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving geom_BFGS_revert_old'

      return
    end subroutine geom_BFGS_revert_old

    subroutine geom_BFGS_backup
      !=======================================================================!
      ! 'contain'ed subroutine to store key data to speed up backtracking     !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      !   none - access all geom_BFGS variables directly as in same scope     !
      !-----------------------------------------------------------------------!
      ! Parent module variables used:                                         !
      !   same as geom_BFGS as in same scope                                  !
      !-----------------------------------------------------------------------!
      ! Modules used:                                                         !
      !   same as geom_BFGS as in same scope                                  !
      !-----------------------------------------------------------------------!
      ! Key Internal Variables:                                               !
      !-----------------------------------------------------------------------!
      ! Necessary conditions:                                                 !
      !-----------------------------------------------------------------------!
      ! Written by Matt Probert, v1.0, 11/05/2002                             !
      !-----------------------------------------------------------------------!
      ! Modified for ONETEP by Arash A Mostofi, 2004                          !
      !=======================================================================!

      use comms, only: comms_barrier, pub_on_root
      use rundat, only: pub_geom_reuse_dk_ngwfs, pub_maxit_ngwf_cg

      implicit none

      ! Local Variables
#ifdef ACCELRYS
      integer :: row           ! vm: required for Intel bug workaround
#endif
      character(len=file_maxpath) :: read_file
#ifdef HDF5
      integer :: ierr, isub
#endif

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering geom_BFGS_backup'

      !always backup cmdl%cell in memory - small arrays so not an issue here
      call castep_cell_copy(cmdl%cell,backup_cell)

      !backup energies
      backup_total_energy = cmdl%total_energy
      backup_enthalpy = enthalpy

      if (pub_on_root) then

#ifdef HDF5
         call h5open_f(ierr)
#endif

         !**********************************************************!
         ! AAM: read current tightbox_ngwfs and copy to backup file !
         !**********************************************************!

         if (pub_geom_reuse_dk_ngwfs.and.(pub_maxit_ngwf_cg>0)) then

            do isub=1,mdl%nsub

               write(stdout,'(3a)',advance='no') ' Backing up NGWFs to file "', &
                    trim(backup_ngwf_file(isub)),'" ...'

               ! Name of current NGWF file to be read
               if (mdl%nsub==1) then
                  write(read_file,'(2a)') trim(pub_rootname),'.tightbox_ngwfs'
               else
                  write(read_file,'(a,a,i0)') trim(pub_rootname),'.tightbox_ngwfs_',isub
               end if
#ifdef HDF5
               read_file=trim(read_file)//'.h5'

               ! Copy file.
               call geom_opt_copy_file_h5(read_file, backup_ngwf_file(isub), &
                    'tightbox_ngwfs')
#else
               ! Copy file.
               call geom_opt_copy_file('geom_BFGS_backup', read_file, &
                    backup_ngwf_file(isub), mdl%regions(isub)%elements, &
                    mdl%regions(isub)%par, is_cmplx)
#endif
            end do

            write(stdout,'(a)') ' done'

         end if

         !***************************************************************!
         ! AAM: read current density kernel from file and backup to file !
         !***************************************************************!

         if (.not.pub_edft.and.pub_geom_reuse_dk_ngwfs) then

            write(stdout,'(3a)',advance='no') &
                 ' Backing up density kernel to file "', &
                 trim(backup_dk_file), '" ...'

            ! Name of current density kernel file to be read
            write(read_file,'(2a)') trim(pub_rootname),'.dkn'
#ifdef HDF5
            read_file=trim(read_file)//'.h5'

            ! Copy file.
            call geom_opt_copy_file_h5(read_file, backup_dk_file, 'sparse')
#else
            ! Copy file.
            call geom_opt_copy_file('geom_BFGS_backup', read_file, &
                 backup_dk_file)
#endif

            write(stdout,'(a)') ' done'

         end if

         !***************************************************************!
         ! ARS: read current Hamiltonian from file and backup to file    !
         !***************************************************************!

         if (pub_edft.and.pub_geom_reuse_dk_ngwfs) then

            write(stdout,'(3a)',advance='no') &
                 ' Backing up Hamiltonian to file "', &
                 trim(backup_ham_file), '" ...'

            ! Name of current density kernel file to be read
            write(read_file,'(2a)') trim(pub_rootname),'.ham'
#ifdef HDF5
            read_file=trim(read_file)//'.h5'

            ! Copy file.
            call geom_opt_copy_file_h5(read_file, backup_ham_file, 'sparse')
#else
            ! Copy file.
            call geom_opt_copy_file('geom_BFGS_backup', read_file, &
                 backup_ham_file)
#endif

            write(stdout,'(a)') ' done'

         end if

#ifdef HDF5
         call h5close_f(ierr)
#endif

      end if

      !backup ionic stuff
      backup_forces = cmdl%forces
      backup_stress = cmdl%stress
      backup_strain = cmdl%strain
      backup_x_vec = x_vec
      backup_f_vec = f_vec

      call comms_barrier

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving geom_BFGS_backup'

      return
    end subroutine geom_BFGS_backup

    subroutine geom_BFGS_restore
      !=======================================================================!
      ! 'contain'ed subroutine to restore key data to speed up backtracking   !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      !   none - access all geom_BFGS variables directly as in same scope     !
      !-----------------------------------------------------------------------!
      ! Parent module variables used:                                         !
      !   same as geom_BFGS as in same scope                                  !
      !-----------------------------------------------------------------------!
      ! Modules used:                                                         !
      !   same as geom_BFGS as in same scope                                  !
      !-----------------------------------------------------------------------!
      ! Key Internal Variables:                                               !
      !-----------------------------------------------------------------------!
      ! Necessary conditions:                                                 !
      !-----------------------------------------------------------------------!
      ! Written by Matt Probert, v1.0, 11/05/2002                             !
      !-----------------------------------------------------------------------!
      ! Modified for ONETEP by Arash A Mostofi, 2004                          !
      !=======================================================================!

      use constants, only: stdout
      use rundat, only: pub_geom_reuse_dk_ngwfs, pub_maxit_ngwf_cg

      implicit none

      ! Local Variables
#ifdef ACCELRYS
      integer :: row           ! vm: required for Intel bug workaround
#endif
      character(len=file_maxpath) :: write_file
#ifdef HDF5
      integer :: ierr
#endif

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering geom_BFGS_restore'

      !restore cmdl%cell from memory
      call castep_cell_copy(backup_cell,cmdl%cell)
      !and make sure that current_cell is in sync
      call castep_cell_copy(cmdl%cell,current_cell)

      ! AAM: bring orig_cell in sync too
      call castep_cell_copy(cmdl%cell,cmdl%orig_cell)

      !restore energies
      cmdl%total_energy = backup_total_energy
      enthalpy = backup_enthalpy

      if (pub_on_root) then

#ifdef HDF5
         call h5open_f(ierr)
#endif

         !*******************************************************************!
         !                AAM: start with tightbox_ngwfs                     !
         !*******************************************************************!

         if ((pub_geom_reuse_dk_ngwfs).and.(pub_maxit_ngwf_cg>0)) then

            do isub=1,mdl%nsub
               ! Name of current NGWF file to be written
               if (mdl%nsub==1) then
                  write(write_file,'(2a)') trim(pub_rootname),'.tightbox_ngwfs'
               else
                  write(write_file,'(a,a,i0)') trim(pub_rootname),'.tightbox_ngwfs_',isub
               end if
#ifdef HDF5
               write_file=trim(write_file)//'.h5'
#endif
               write(stdout,'(3a)',advance='no') &
                    'Restoring NGWFs to file "', trim(write_file),'" ...'

               ! Copy file.
#ifdef HDF5
               call geom_opt_copy_file_h5(backup_ngwf_file(isub), write_file, &
                    'tightbox_ngwfs')
#else
               call geom_opt_copy_file('geom_BFGS_restore', backup_ngwf_file(isub), &
                    write_file, mdl%regions(isub)%elements, &
                    mdl%regions(isub)%par, is_cmplx)
#endif
            end do

            write(stdout,'(a)') ' done'

         end if

         !********************************************************************!
         !        AAM: read density kernel from backup file and restore       !
         !********************************************************************!

         if (.not.pub_edft.and.pub_geom_reuse_dk_ngwfs) then

            ! Name of current density kernel file to be written
            write(write_file,'(2a)') trim(pub_rootname),'.dkn'
#ifdef HDF5
            write_file=trim(write_file)//'.h5'
#endif
            write(stdout,'(3a)',advance='no') &
                 'Restoring density kernel to file "', trim(write_file),'" ...'

            ! Copy file.
#ifdef HDF5
            call geom_opt_copy_file_h5(backup_dk_file, write_file, 'sparse')
#else
            call geom_opt_copy_file('geom_BFGS_restore', backup_dk_file, &
                 write_file)
#endif

            write(stdout,'(a)') ' done'

         end if

         !********************************************************************!
         !        ARS: read Hamiltonian from backup file and restore          !
         !********************************************************************!

         if (pub_edft.and.pub_geom_reuse_dk_ngwfs) then

            ! Name of current Hamiltonian file to be written
            write(write_file,'(2a)') trim(pub_rootname),'.ham'
#ifdef HDF5
            write_file=trim(write_file)//'.h5'
#endif
            write(stdout,'(3a)',advance='no') &
                 'Restoring density kernel to file "', trim(write_file),'" ...'

            ! Copy file.
#ifdef HDF5
            call geom_opt_copy_file_h5(backup_ham_file, write_file, 'sparse')
#else
            call geom_opt_copy_file('geom_BFGS_restore', backup_ham_file, &
                 write_file)
#endif

            write(stdout,'(a)') ' done'

         end if

#ifdef HDF5
         call h5close_f(ierr)
#endif

      end if

      !restore ionic stuff
      cmdl%forces = backup_forces
      cmdl%stress = backup_stress
      cmdl%strain = backup_strain
      x_vec = backup_x_vec
      f_vec = backup_f_vec

      !restore model flags
      cmdl%found_ground_state = .true.

      !NB No need to call wave_Sorthonormalise as wvfn was Sorthonormal when stored

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving geom_BFGS_restore'

      return
    end subroutine geom_BFGS_restore

    subroutine geom_BFGS_write_summary()
      !=========================================================================!
      ! contain'ed subroutine to write out the BFGS minimisation summary.       !
      !-------------------------------------------------------------------------!
      ! Arguments:                                                              !
      !    none                                                                 !
      !-------------------------------------------------------------------------!
      ! Parent module variables used:                                           !
      !    none                                                                 !
      !-------------------------------------------------------------------------!
      ! Modules used:                                                           !
      !    none                                                                 !
      !-------------------------------------------------------------------------!
      ! Key Internal Variables:                                                 !
      !    none                                                                 !
      !-------------------------------------------------------------------------!
      ! Necessary conditions:                                                   !
      !-------------------------------------------------------------------------!
      ! Written by Matt Probert, v1.0, 11/05/2002                               !
      !-------------------------------------------------------------------------!
      ! Modified for ONETEP by Arash A Mostofi, 2004                            !
      !=========================================================================!


      use constants, only: dp, stdout
      use services, only: services_flush

      !local vars
      character(len=80) :: data_string
      character(len=80) :: divider_string
      character(len=80) :: label_string
      integer           :: string_index
      integer           :: len_label, len_lambda, len_fdelta, len_enthalpy
      real (kind=dp)    :: enthalpy, f_dot_delta

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering geom_BFGS_write_summary'

      !initialise strings
      data_string     = repeat(' ',len(data_string))
      divider_string  = repeat(' ',len(divider_string))
      label_string    = repeat(' ',len(label_string))
      !                  12345678901234 2345678901234 2345678901234 23456789012345678
      divider_string  = '+------------+-------------+-------------+-----------------+'
      label_string    = '|    Step    |   lambda    |   F.delta   |    enthalpy     |'
      !e.g.              | trial step    -100.12345    -100.12345    -123456.123456

      if(pub_on_root) then

         !header
         write(stdout,*) ' '
         select case (trim(geom_string))
         case ('BFGS')
            write(stdout,5) divider_string
            write(stdout,5) label_string
            write(stdout,5) divider_string
         case ('LBFGS')
            write(stdout,6) divider_string
            write(stdout,6) label_string
            write(stdout,6) divider_string
         end select

         !previous                                                    '123456789012'
         string_index=1
         len_label   =14
         len_lambda  =14
         len_fdelta  =14
         len_enthalpy=18
         write(data_string(string_index:string_index+len_label),   1) '  previous  '
         string_index=string_index+len_label
         write(data_string(string_index:string_index+len_lambda),  2) 0.0_dp
         string_index=string_index+len_lambda
         f_dot_delta=old_f_dot_delta
         if (abs(f_dot_delta)>1.0E4.or.abs(f_dot_delta)<1.0E-5) then
            write(data_string(string_index:string_index+len_fdelta), 33) f_dot_delta
         else
            write(data_string(string_index:string_index+len_fdelta),  3) f_dot_delta
         end if
         string_index=string_index+len_fdelta
         enthalpy=old_enthalpy
         if (abs(enthalpy)>1.0E7.or.abs(enthalpy)<1.0E-5) then
            write(data_string(string_index:string_index+len_enthalpy),34) enthalpy
         else
            write(data_string(string_index:string_index+len_enthalpy),4) enthalpy
         end if
         string_index=string_index+len_enthalpy
         !cross check all is well (NB strings start from 1 not 0 ...)
         if (string_index>71) write (stdout,*) 'geom_BFGS_write_summary: format problem?'
         if(trim(geom_string)=='BFGS')  write(stdout,5) data_string
         if(trim(geom_string)=='LBFGS') write(stdout,6) data_string

         !trial                                                          '123456789012'
         if (done_trial) then
            string_index=1
            write(data_string(string_index:string_index+len_label),   1) ' trial step '
            string_index=string_index+len_label
            write(data_string(string_index:string_index+len_lambda),  2) lambda_trial
            string_index=string_index+len_lambda
            f_dot_delta=trial_f_dot_delta
            if (abs(f_dot_delta)>1.0E4.or.abs(f_dot_delta)<1.0E-5) then
               write(data_string(string_index:string_index+len_fdelta), 33) f_dot_delta
            else
               write(data_string(string_index:string_index+len_fdelta),  3) f_dot_delta
            end if
            string_index=string_index+len_fdelta
            enthalpy=trial_enthalpy
            if (abs(enthalpy)>1.0E7.or.abs(enthalpy)<1.0E-5) then
               write(data_string(string_index:string_index+len_enthalpy),34) enthalpy
            else
               write(data_string(string_index:string_index+len_enthalpy),4) enthalpy
            end if
            string_index=string_index+len_enthalpy
            if(trim(geom_string)=='BFGS') write(stdout,5) data_string
            if(trim(geom_string)=='LBFGS') write(stdout,6) data_string
         end if

         !line                                                           '123456789012'
         if (done_line) then
            string_index=1
            write(data_string(string_index:string_index+len_label),   1) '  line step '
            string_index=string_index+len_label
            write(data_string(string_index:string_index+len_lambda),  2) lambda_line
            string_index=string_index+len_lambda
            f_dot_delta=line_f_dot_delta
            if (abs(f_dot_delta)>1.0E4.or.abs(f_dot_delta)<1.0E-5) then
               write(data_string(string_index:string_index+len_fdelta), 33) f_dot_delta
            else
               write(data_string(string_index:string_index+len_fdelta),  3) f_dot_delta
            end if
            string_index=string_index+len_fdelta
            enthalpy=line_enthalpy
            if (abs(enthalpy)>1.0E7.or.abs(enthalpy)<1.0E-5) then
               write(data_string(string_index:string_index+len_enthalpy),34) enthalpy
            else
               write(data_string(string_index:string_index+len_enthalpy),4) enthalpy
            end if
            string_index=string_index+len_enthalpy
            if(trim(geom_string)=='BFGS') write(stdout,5) data_string
            if(trim(geom_string)=='LBFGS') write(stdout,6) data_string
         end if

         !quad                                                           '123456789012'
         if (done_quad) then
            string_index=1
            write(data_string(string_index:string_index+len_label),   1) '  quad step '
            string_index=string_index+len_label
            write(data_string(string_index:string_index+len_lambda),  2) lambda_quad
            string_index=string_index+len_lambda
            f_dot_delta=quad_f_dot_delta
            if (abs(f_dot_delta)>1.0E4.or.abs(f_dot_delta)<1.0E-5) then
               write(data_string(string_index:string_index+len_fdelta), 33) f_dot_delta
            else
               write(data_string(string_index:string_index+len_fdelta),  3) f_dot_delta
            end if
            string_index=string_index+len_fdelta
            enthalpy=quad_enthalpy
            if (abs(enthalpy)>1.0E7.or.abs(enthalpy)<1.0E-5) then
               write(data_string(string_index:string_index+len_enthalpy),34) enthalpy
            else
               write(data_string(string_index:string_index+len_enthalpy),4) enthalpy
            end if
            string_index=string_index+len_enthalpy
            if(trim(geom_string)=='BFGS') write(stdout,5) data_string
            if(trim(geom_string)=='LBFGS') write(stdout,6) data_string

         end if

         !footer
         if(trim(geom_string)=='BFGS') write(stdout,5) divider_string
         if(trim(geom_string)=='LBFGS') write(stdout,6) divider_string
         write(stdout,*) ' '

      end if          !pub_on_root
      call services_flush

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving geom_BFGS_write_summary'

1     format('|',a11,1x,'|')        !1+11+1+1=14
2     format(1x,f11.6,1x,'|')       !1+11+1+1=14
3     format(1x,f11.6,1x,'|')       !1+11+1+1=14
33    format(1x,es11.3e3,1x,'|')    !1+11+1+1=14
4     format(1x,f15.6,1x,'|')       !1+15+1+1=18
34    format(1x,es15.6e3,1x,'|')    !1+15+1+1=18
      !14+14+14+18=60
5     format(1x,a60,' <-- min BFGS')   !1+60+14 => ends col 75
6     format(1x,a60,' <-- min LBFGS')  !1+60+15 => ends col 76

    end subroutine geom_BFGS_write_summary
  end subroutine geom_BFGS


  subroutine geom_DELOC(cmdl,mdl,output_file,path,is_cmplx,status)
    !=========================================================================!
    ! Find ground state ionic positions for given model using delocalized     !
    ! internal coordinates. This technology is based on the DMol3             !
    ! implementation which is described in:                                   !
    ! J. Andzelm, R. D. King-Smith, G. Fitzgerald,                            !
    ! Chem. Phys. Lett. 335 (2001) 321-326                                    !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=inout, the model to be optimised                          !
    !   output_file, intent=in, the filename to which results will be written.!
    !   status, intent=out                                                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model                                                                 !
    !   parameters                                                            !
    !   comms                                                                 !
    !   geometry_deloc_utils                                                  !
    !   geometry_deloc_algor                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   cmdl has already been properly read and initialised                    !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !=========================================================================!

    use simulation_cell, only: castep_model, castep_cell_alloc, &
         castep_cell_copy, castep_cell_dealloc
    use comms, only: comms_bcast, pub_root_proc_id, pub_on_root
    use constants, only: dp, stdout
    use geometry_deloc_utils, only: deloc_utils_initialize,&
         deloc_utils_mdl_to_internals, deloc_utils_update_params, &
         deloc_utils_energy_gradient, deloc_utils_output_converged, &
         deloc_utils_internals_to_mdl, deloc_utils_deallocate,deloc_utils_nullify
    use geometry_deloc_algor, only: geom_algor_deloc
    use linalg, only: linalg_invert_sym_matrix
    use model_type, only: MODEL
    use rundat, only: pub_geom_backup_iter, pub_geom_output_detail, &
         pub_geom_frequency_est, pub_geom_modulus_est, pub_geom_max_iter, &
         pub_write_positions
    use services, only: services_flush
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    use neb, only: neb_path, neb_sync_path, neb_converged

    implicit none

    type(castep_model), intent(inout)               :: cmdl
    type(MODEL),        intent(inout)               :: mdl
    character(len=*), intent(in)                    :: output_file ! trajectory
    type(neb_path),optional,intent(inout)           :: path
    ! agrecocmplx
    logical, intent(in)                             :: is_cmplx
    integer,intent(out)                             :: status
    !-------------------------------------------------------------------------!

    type (castep_cell) :: init_cell !the current_cell when last initialised inv_Hessian
    integer :: ierr,fake_iter
    logical :: converged

    real(kind=DP) :: enthalpy
    real(kind=DP) :: dE
    real(kind=DP) :: Fmax   !max Force
    real(kind=DP) :: dRmax  !max displacement
    real(kind=DP) :: Smax   !max Stress
    logical       :: first_step=.true.
    logical       :: converged_dE, converged_Fmax, converged_dRmax, converged_Smax
    real(kind=DP), allocatable, dimension(:)     :: old_x_vec !old strains and fractional coords
    logical       :: hessian_nuked ! ja531 -> is the Hessian bad?
    !only to be able to use geometry_converged
    !-------------------------------------------------------------------------!
    status = 0

    old_frequency_est=pub_geom_frequency_est             !we store the starting values as 'old' to test o/p at end
    new_frequency_est=pub_geom_frequency_est             !we store the working value as 'new' s.t. can cope with
    new_modulus_est  =pub_geom_modulus_est               !automatic and user updates to the estimates

    geom_string = "DI"
    call comms_bcast(pub_root_proc_id,geom_string)

    ! ja531 -> March 2015 : Changed mind about this : DI inv Hessian is Cartesian!
!    call utils_assert(trim(geom_string)==trim(geom_continuation_string), &
!   "Error: Cannot continue non-DI geometry optimisation as DI.")

    allocate(old_x_vec(1:ndim),stat=ierr)
    call utils_alloc_check('geom_DELOC','old_x_vec',ierr)

    if (pub_on_root) write(stdout,'(/a,a,a)') &
         ' DI:   continuation file name is "',trim(cont_file),'"'

    call castep_cell_alloc(init_cell,mdl%nat+mdl%nat_classical)

    call castep_cell_copy(cmdl%cell,init_cell)

    !-------------------------------------------------------------------------!
    !Sort out start from restart                                              !
    !-------------------------------------------------------------------------!

    !If this is not a restart, then need to initialise the inv_Hessian
    if (cmdl%bfgs_iteration==0) then

       if(pub_on_root) call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,&
            reset_updates=.true.)
       cmdl%strain=0.0_dp

       ! agrecocmplx
       if (present(path)) call neb_sync_path(mdl, path)
       call geom_get_forces(cmdl,mdl,path,is_cmplx,en_sync=.true.)

       enthalpy = cmdl%total_energy
       !start up the output file
       call geom_write_trajectory(cmdl,enthalpy,output_file,mdl)

    else

       if (pub_on_root)write (stdout,*) 'Delocalized Internals: Restarting previous geometry optimization.'

       !Catch the possibility that model_continuation/reuse has decided to nuke the inv_Hessian
       ! ja531 -> mdl%bfgs_inv_Hessian is only allocated on root
       if(pub_on_root) then
          if(cmdl%bfgs_optimisation) then
             hessian_nuked=all(abs(cmdl%bfgs_inv_Hessian)<tiny(0.0_dp))
          end if
          if(cmdl%lbfgs_optimisation) then
             hessian_nuked=(all(abs(cmdl%lbfgs_position_updates)<tiny(0.0_dp)).or.&
                  &all(abs(cmdl%lbfgs_gradient_updates)<tiny(0.0_dp)))
          end if
       end if
       call comms_bcast(pub_root_proc_id,hessian_nuked)
       if (hessian_nuked) then
          call castep_cell_copy(cmdl%cell,init_cell)
          if(pub_on_root) call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,&
               mdl,reset_updates=.true.)
       else
          ! ja531 -> March 2015: In the case we are reading from an l-BFGS continuation file, convert to an inverse Hessian.
          if(cmdl%lbfgs_optimisation) then
             if (.not.allocated(R_mat)) then
                allocate(R_mat(1:lbfgs_block_length,1:lbfgs_block_length),stat=ierr)
                call utils_alloc_check('geom_DELOC','R_mat',ierr)
             end if
             if (.not.allocated(YY_mat)) then
                allocate(YY_mat(1:lbfgs_block_length,1:lbfgs_block_length),stat=ierr)
                call utils_alloc_check('geom_DELOC','YY_mat',ierr)
             end if

             if (.not.allocated(tmpmm)) then
                allocate(tmpmm(1:lbfgs_block_length*2,1:lbfgs_block_length*2),stat=ierr)
                call utils_alloc_check('geom_DELOC','tmpmm',ierr)
             end if
             if (.not.allocated(tmp2m)) then
                allocate(tmp2m(1:lbfgs_block_length*2),stat=ierr)
                call utils_alloc_check('geom_DELOC','tmp2m',ierr)
             end if

             if (.not.allocated(cmdl%lbfgs_position_updates)) then
                allocate(cmdl%lbfgs_position_updates(1:ndim,1:lbfgs_block_length),stat=ierr)
                call utils_alloc_check('geom_DELOC','=>cmdl%lbfgs_position_updates',ierr)
             end if
             if (.not.allocated(cmdl%lbfgs_gradient_updates)) then
                allocate(cmdl%lbfgs_gradient_updates(1:ndim,1:lbfgs_block_length),stat=ierr)
                call utils_alloc_check('geom_DELOC','=>cmdl%lbfgs_gradient_updates',ierr)
             end if
             if (.not.allocated(lbfgs_init_inv_hess)) then
                allocate(lbfgs_init_inv_hess(1:3,1:ndim),stat=ierr)
                call utils_alloc_check('geom_DELOC','lbfgs_init_inv_hess',ierr)
             end if
             if (.not.allocated(lbfgs_init_hess)) then
                !             if (constrain_cell.and..not.fix_all_cell) then
                allocate(lbfgs_init_hess(1:3,1:ndim),stat=ierr)
                call utils_alloc_check('geom_DELOC','lbfgs_init_hess',ierr)
                !             end if
             end if

             if (.not.allocated(tmpn)) then
                allocate(tmpn(1:ndim),stat=ierr)
                call utils_alloc_check('geom_DELOC','tmpn',ierr)
             end if
             if (.not.allocated(D_mat)) then
                allocate(D_mat(1:lbfgs_block_length),stat=ierr)
                call utils_alloc_check('geom_DELOC','D_mat',ierr)
             end if

             if(pub_on_root) then
                if(.not.allocated(cmdl%bfgs_inv_Hessian)) then
                   allocate(cmdl%bfgs_inv_Hessian(1:ndim,1:ndim),stat=ierr)
                   call utils_alloc_check('geom_DELOC','=>cmdl%bfgs_inv_Hessian', ierr)
                end if
             end if

             if((hess_init_method.eq.hess_init_identity).or.(hess_init_method.eq.hess_init_scaled)) then
                if(pub_on_root) call geom_inv_Hessian_initialize(toggle_on,init_cell,&
                     cmdl,mdl,init_method=hess_init_identity)
             else
                if(pub_on_root) call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl)
             end if

             call lbfgs_restart(toggle_on,init_cell,cmdl,mdl)
             if (.not.allocated(tmpmm)) then
                allocate(tmpmm(1:lbfgs_block_length*2,1:lbfgs_block_length*2),stat=ierr)
                call utils_alloc_check('geom_DELOC','tmpmm',ierr)
             end if
             call linalg_invert_sym_matrix(R_mat,cmdl%lbfgs_num_updates)
             tmpmm(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)=0.0_dp
             tmpmm(1:cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)=-R_mat
             tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,&
                  1:cmdl%lbfgs_num_updates)= transpose(tmpmm(1:cmdl%lbfgs_num_updates,&
                  cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2))
             tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)=&
                  matmul(matmul(tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,&
                  1:cmdl%lbfgs_num_updates),YY_mat(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)),&
                  tmpmm(1:cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2))
             deallocate(R_mat,stat=ierr)
             call utils_dealloc_check('geom_DELOC','R_mat',ierr)
             deallocate(YY_mat,stat=ierr)
             call utils_dealloc_check('geom_DELOC','YY_mat',ierr)
             cmdl%bfgs_inv_Hessian=cmdl%bfgs_inv_Hessian+matmul(matmul(cmdl%lbfgs_position_updates(:,1:cmdl%lbfgs_num_updates),&
                  tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,1:cmdl%lbfgs_num_updates)),&
                  transpose(cmdl%lbfgs_gradient_updates(:,1:cmdl%lbfgs_num_updates)))
             cmdl%bfgs_inv_Hessian=cmdl%bfgs_inv_Hessian+matmul(matmul(cmdl%lbfgs_gradient_updates(:,1:cmdl%lbfgs_num_updates),&
                  tmpmm(1:cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)),&
                  transpose(cmdl%lbfgs_position_updates(:,1:cmdl%lbfgs_num_updates)))
             cmdl%bfgs_inv_Hessian=cmdl%bfgs_inv_Hessian+matmul(matmul(cmdl%lbfgs_gradient_updates(:,1:cmdl%lbfgs_num_updates),&
                  tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)),&
                  transpose(cmdl%lbfgs_gradient_updates(:,1:cmdl%lbfgs_num_updates)))
             deallocate(tmpmm,stat=ierr)
             call utils_dealloc_check('geom_DELOC','tmpmm',ierr)
             deallocate(cmdl%lbfgs_position_updates,stat=ierr)
             call utils_dealloc_check('geom_DELOC','=>cmdl%lbfgs_position_updates',ierr)
             deallocate(cmdl%lbfgs_gradient_updates,stat=ierr)
             call utils_dealloc_check('geom_DELOC','=>cmdl%lbfgs_gradient_updates',ierr)

             ! ja531 -> March 2015 : We probably allocated this in lbfgs_restart
             if(allocated(lbfgs_init_inv_hess)) then
                deallocate(lbfgs_init_inv_hess,stat=ierr)
                call utils_dealloc_check('geom_DELOC','lbfgs_init_inv_hess',ierr)
             end if
             if(allocated(lbfgs_init_hess)) then
                deallocate(lbfgs_init_hess,stat=ierr)
                call utils_dealloc_check('geom_DELOC','lbfgs_init_hess',ierr)
             end if
             if(allocated(tmpn)) then
                deallocate(tmpn,stat=ierr)
                call utils_dealloc_check('geom_DELOC','tmpn',ierr)
             end if
             if(allocated(D_mat)) then
                deallocate(D_mat,stat=ierr)
                call utils_dealloc_check('geom_DELOC','D_mat',ierr)
             end if

             cmdl%lbfgs_optimisation=.false.
             cmdl%bfgs_optimisation=.true.
          end if
       end if

       ! ja531 -> March 2015: Need to update ethalpy!
       enthalpy = cmdl%total_energy
    end if

    ! initialize geometry_deloc_utils module
    call deloc_utils_initialize(cmdl,pub_geom_output_detail)

    !NB We reset bfgs_iteration here as the geom_max_iter refers to number of iterations
    !   THIS pass and not in total.
    cmdl%bfgs_iteration=0
    reset_iter=cmdl%bfgs_iteration

    iteration=cmdl%bfgs_iteration
    call comms_bcast(pub_root_proc_id,iteration)

    ! fill symmetry information
    call deloc_utils_mdl_to_internals(cmdl,mdl%elements)

    ! generate initial delocalized coordinates
    fake_iter = 0
    if (pub_on_root) then
       call geom_algor_DELOC(.true.,fake_iter,status)
    end if
    call comms_bcast(pub_root_proc_id,status)
    if (status/=0) go to 999

    !Put the strains wrt original cell and fractional atomic coords wrt current_cell into old_x_vec
    call geom_mdl_to_xvec(cmdl,old_x_vec,mdl)

    main_loop : do

       ! check to see if we have done enough iterations yet?
       !NB Test is done this way around s.t. we don't need to undo anything at end
       if (iteration>=pub_geom_max_iter) exit main_loop

       !checkpoint and update as appropriate
       call deloc_utils_update_params

       ! get energy and gradients
       call deloc_utils_energy_gradient                ! SCF energy/gradient

       if (first_step) then
          !initialise convergence tests (impossible to converge here with convergence windows)
          call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
               & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

          !tell user where we are starting from
          call deloc_utils_output_converged(iteration,enthalpy,dE,Fmax,dRmax, &
                converged_dE,converged_Fmax,converged_dRmax)

          first_step = .false.

       endif

       cmdl%bfgs_iteration = cmdl%bfgs_iteration + 1             ! next geometry cycle
       iteration=cmdl%bfgs_iteration

       !banner to output
       if (pub_on_root) write (stdout,4) 'Starting DI iteration ',iteration,' ...'
       if (present(path)) call neb_sync_path(mdl, path)

       !store the current model?

       !Put the strains wrt original cell and fractional atomic coords wrt current_cell into old_x_vec
       call geom_mdl_to_xvec(cmdl,old_x_vec,mdl)

       call deloc_utils_mdl_to_internals(cmdl,mdl%elements)

       !info on current_cell to stdout
       ! make a step of geometry optimization

       if (pub_on_root) then
          ! Do all the things on the root proc (i.e. look for  new coordinates):
          call geom_algor_DELOC(.false.,iteration,status)
       end if
       call comms_bcast(pub_root_proc_id,status)
       if (status/=0) go to 999

       ! update model data
       call deloc_utils_internals_to_mdl(cmdl)
       call castep_cell_copy(cmdl%cell,current_cell)

       !tell user what's happening
       if (pub_on_root) write (stdout,5) 'DI: structure has been updated on iteration',iteration

       !tell user what's the current structure
       if (pub_write_positions) call geom_output(cmdl,mdl)

       ! find the forces and energy
       ! agrecocmplx
       call geom_get_forces(cmdl,mdl,path,is_cmplx,en_sync=.true.)

       call services_flush

       ! save trajectory data
       enthalpy = cmdl%total_energy
       call geom_write_trajectory(cmdl,enthalpy,output_file,mdl)

       !test convergence
       call geom_converged(cmdl,old_x_vec,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
            & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)

       !tell user what just happened (stop iterations if converged)
       call deloc_utils_output_converged(iteration,enthalpy,dE,Fmax,dRmax, &
             converged_dE,converged_Fmax,converged_dRmax)

       if (present(path)) converged = neb_converged(converged, iteration, dE, Fmax, dRmax, path)

       if (converged) exit main_loop


       !checkpoint and update as appropriate
       call deloc_utils_update_params
       if (mod(iteration,pub_geom_backup_iter)==0) &
            call geom_opt_continuation_write(cmdl, mdl, cont_file, is_cmplx)

    enddo main_loop

    if (pub_on_root) then
       if (converged) then
          write(stdout,*) 'DI: Geometry optimization completed successfully.'
       else
          write(stdout,*) 'DI: Geometry optimization failed to converge after',iteration,'steps.'
       end if
       write (stdout,11) 'DI: Final Configuration:'
       call geom_output(cmdl,mdl)
       write (stdout,2) 'DI: Final Enthalpy     =',enthalpy,trim(energy_label)
    end if


999 continue

    call deloc_utils_deallocate
    call deloc_utils_nullify

    deallocate(old_x_vec,stat=ierr)
    call utils_dealloc_check('geom_DELOC','old_x_vec',ierr)

    call castep_cell_dealloc(init_cell)

2   format(1x,a,es17.8e3,1x,a)
4   format(/,80('='),/,1x,a,i4,a,/,80('='))
5   format(/,80('-'),/,1x,a,i4,/,80('-'))
11  format(/,80('='),/,1x,a,/,80('='))

    return

  end subroutine geom_DELOC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine geom_xvec_to_mdl(x_vec,cmdl,mdl)
    !=========================================================================!
    ! Updates cmdl and current_cell for current strains and fractional ionic   !
    ! positions in x_vec                                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   x_vec, intent=in, the ndim vector of search parameters                !
    !   cmdl,intent=inout, the model to be updated for new cell & coords       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !   fixed_npw to see which kind of basis re-initialisation is required    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model for type definitions                                            !
    !   cell to update components of cell                                     !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use constants, only: DP, stdout
    use model_type, only: MODEL
    use simulation_cell, only: castep_model, castep_cell_copy, &
         castep_model_cell_changed

    implicit none

    real(kind=DP), dimension(1:ndim), intent(in)  :: x_vec
    type(MODEL), intent(in)                       :: mdl  ! jcap
    type(castep_model), intent(inout)             :: cmdl

    !local variables
    integer                                       :: i,j,iatom,ispec

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_xvec_to_mdl'

    !Only need to update the fractional atomic coords if allowed to move
    if (.not.fix_all_ions) then
       !NB We no longer apply PBCs here and so cmdl%cell is allowed to become
       !   unrationalized, i.e. fractional coords>1 etc.
       i=10
       do ispec=1,cmdl%cell%num_species - mdl%nat_classical
          do iatom=1,cmdl%cell%num_ions_in_species(ispec)
             do j=1,3
                cmdl%cell%ionic_positions(j,iatom,ispec)=x_vec(i+j-1)
             end do
             i=i+3
          end do
       end do
    end if

    !Not much to do for the basis here:
    !Just synchronize cmdl%cell with current_cell ...
    call castep_cell_copy(cmdl%cell,current_cell)

    !... and update model flags and S-orthogonalise wvfn etc
    call castep_model_cell_changed(cmdl)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_xvec_to_mdl'

    return
  end subroutine geom_xvec_to_mdl

  subroutine geom_fvec_to_mdl(f_vec,cmdl,mdl)
    !=========================================================================!
    ! Updates cmdl and current_cell for current stresses and fractional forces !
    ! positions in f_vec                                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   f_vec, intent=in, the ndim vector of fractional forces                !
    !   cmdl,intent=inout, the model to be updated                             !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model for type definitions                                            !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use constants, only: dp, stdout
    use linalg, only: linalg_invert_sym_matrix
    use model_type, only: MODEL
    use simulation_cell, only: castep_model, castep_cell_frac_to_cart

    implicit none

    type(MODEL), intent(in)                       :: mdl ! jcap
    real(kind=DP), dimension(1:ndim), intent(in)  :: f_vec
    type(castep_model), intent(inout)             :: cmdl

    !local variables
    real(kind=DP), dimension(1:3,1:3)             :: metric, inv_metric
    real(kind=DP), dimension(1:3)                 :: frac_forces
    integer                                       :: i,iatom,ispec

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_fvec_to_mdl'

    !NB WE REQUIRE that geom_xvec_to_mdl has been called immediately before this routine
    !   s.t. cmdl%cell and cmdl%strain are in sync with x_vec

    !calculate metric in same form as inv_Hessian, metric=h^T.h
    metric=matmul(cmdl%cell%real_lattice,transpose(cmdl%cell%real_lattice))
    inv_metric=metric ; call linalg_invert_sym_matrix(inv_metric,3)

    if (.not.fix_all_ions) then
       i=10
       do ispec=1,cmdl%cell%num_species - mdl%nat_classical
          do iatom=1,cmdl%cell%num_ions_in_species(ispec)

             !convert forces from fractional form
             frac_forces=matmul(inv_metric,f_vec(i:i+2))

             call castep_cell_frac_to_cart(current_cell,frac_forces,cmdl%forces(:,iatom,ispec))

             !next ion
             i=i+3
          end do
       end do
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_fvec_to_mdl'

    return
  end subroutine geom_fvec_to_mdl

  subroutine geom_mdl_to_xvec(cmdl,x_vec,mdl)
    !=========================================================================!
    ! Extract strains wrt original cell and fractional ionic coords from cmdl  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl,   intent=in, the model to extract strain & coords from           !
    !   x_vec, intent=out,the ndim vector of search parameters                !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   io for error messages                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    ! lk: added VERBOSE to constants
    use constants, only: dp, stdout, NORMAL, VERBOSE
    use comms, only: pub_on_root
    use model_type, only: MODEL
    use rundat, only: pub_geom_output_detail
    use simulation_cell, only: castep_model

    implicit none

    type(MODEL), intent(in)                       :: mdl
    type(castep_model), intent(in)                :: cmdl
    real(kind=DP), dimension(1:ndim), intent(out) :: x_vec

    !local variables
    integer                                       :: i,j,iatom,ispec

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_mdl_to_xvec'

    !Pack strain into x_vec using same order as stresses into f_vec
    x_vec(1:9)=reshape(cmdl%strain,(/9/))

    !... and pack the fractional atomic coords into x_vec
    i=10



    do ispec=1,cmdl%cell%num_species - mdl%nat_classical
       do iatom=1,cmdl%cell%num_ions_in_species(ispec)
          do j=1,3
             x_vec(i+j-1)=cmdl%cell%ionic_positions(j,iatom,ispec)
          end do
          i=i+3
       end do
    end do
    ! lk: changed output detail from > NORMAL to > VERBOSE.
    if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
       do i=1,ndim
          write (stdout,99) i,x_vec(i)
       end do
    end if
99  format('x_vec(',i4,')=',f10.5)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_mdl_to_xvec'

    return
  end subroutine geom_mdl_to_xvec

  subroutine geom_mdl_to_fvec(cmdl,f_vec,enthalpy,mdl)
    !=========================================================================!
    ! Extract stresses and forces from cmdl                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl,   intent=in, the model to extract strain & coords from           !
    !   f_vec,         intent=out,the ndim vector of gradients                !
    !   enthalpy,      intent=out,the enthalpy=Etot+pV of the system          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   cell for external pressure and cart_to_frac conversion                !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   stress = total unsymmetrised stress, including external pressure      !
    !   strain = difference between current cell & original cell              !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout, NORMAL, VERBOSE, two_pi
    use model_type, only: MODEL
    use rundat, only: pub_geom_output_detail
    use simulation_cell, only: castep_model

    implicit none

    type(MODEL), intent(in)                       :: mdl ! jcap
    type(castep_model), intent(in)                :: cmdl
    real(kind=DP), dimension(1:ndim), intent(out) :: f_vec
    real(kind=DP),    intent(out)                 :: enthalpy

    !local variables
    real(kind=DP), dimension(1:3,1:3)             :: metric
    real(kind=DP), dimension(1:3)                 :: frac_forces
    real(kind=DP)                                 :: pressure
    real(kind=DP)                                 :: shear_energy

    integer                                       :: i,iatom,ispec

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_mdl_to_fvec'

    pressure=(external_pressure(1)+external_pressure(2)+external_pressure(3))/3.0_dp

    !and now include the ionic forces (Pfrommer eq 5) ...
    !calculate metric in same form as inv_Hessian, metric=h^T.h
    metric=matmul(cmdl%cell%real_lattice,transpose(cmdl%cell%real_lattice))

    ! jd: Takes care of the part of f_vec that would otherwise remain
    !     uninitialized. Since inv_Hessian is not necessarily diagonal,
    !     leaving garbage in f_vec(1:9) is risky.
    f_vec(1:9) = 0.0_DP

    if (.not.fix_all_ions) then
       i=10
       do ispec=1,cmdl%cell%num_species - mdl%nat_classical
          do iatom=1,cmdl%cell%num_ions_in_species(ispec)

             !convert forces from model to fractional form
             !NB we DO NOT use cell_cart_to_frac as this imposes PBCs on frac
             !   which is disastrous for forces!!!
             frac_forces=matmul(cmdl%cell%recip_lattice,cmdl%forces(:,iatom,ispec))/two_pi

             !store metric*fractional_forces in f_vec
             f_vec(i:i+2)=matmul(metric,frac_forces)

             !next ion
             i=i+3
          end do
       end do
    end if

    enthalpy=cmdl%total_energy+pressure*cmdl%cell%volume

    !Calculate elastic shear energy contribution (off diagonal pressure) according to Kittel (3rd ed) eq 4.14
    !NB this requires strain wrt EQUILIBRIUM => only correct if ref_cell is equilibrated ...
    !NB if external pressure is shear, then we ought to have no symmetry ops and so cannot presume
    !   that strain is symmetric => cannot do factor of 0.5 & 2 cancellation ...
    shear_energy=0.0_dp
    !pdh: remove any() structure which upsets gfortran
    if (abs(external_pressure(4))>epsilon(1.0_dp).or.abs(external_pressure(5))>epsilon(1.0_dp).or. &
         &abs(external_pressure(6))>epsilon(1.0_dp)) then
       shear_energy=shear_energy+(cmdl%stress(4)+external_pressure(4))*(cmdl%strain(2,3)+cmdl%strain(3,2))
       shear_energy=shear_energy+(cmdl%stress(5)+external_pressure(5))*(cmdl%strain(3,1)+cmdl%strain(1,3))
       shear_energy=shear_energy+(cmdl%stress(6)+external_pressure(6))*(cmdl%strain(1,2)+cmdl%strain(2,1))
       shear_energy=0.5_dp*cmdl%cell%volume*shear_energy
       enthalpy=enthalpy+shear_energy
    end if

    if (pub_geom_output_detail>NORMAL.and.pub_on_root) then
       write (stdout,*)  'Enthalpy contributions:'
       write (stdout,98) 'total energy',cmdl%total_energy,trim(energy_label)
       write (stdout,98) 'diagonal pressure term',pressure*cmdl%cell%volume,trim(energy_label)
       write (stdout,98) 'off-diagonal pressure term',shear_energy,trim(energy_label)
       if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
          do i=1,ndim
             write (stdout,99) i,f_vec(i)
          end do
       end if
    end if

98  format(a30,'=',1x,es17.8e3,1x,a)
99  format('f_vec(',i4,')=',f10.5)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_mdl_to_fvec'

    return
  end subroutine geom_mdl_to_fvec

  subroutine geom_check_lambda(lambda,delta_vec,scale)
    !=========================================================================!
    ! Check proposed value of lambda does not lead to too large a displacement!
    ! and if it does, reduce it until tolerance satisfied.                    !
    ! Large values of lambda in themselves are not a problem - only if cause  !
    ! an unreasonably large change in cell or ionic coords.                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   lambda, intent=inout, the lambda to be tested/updated                 !
    !   delta_vec, intent=in, the displacement vector to check against        !
    !   old_x_vec, intent=in, the reference structure to check delta against  !
    !   cmdl, intent=in, the model for the proposed cell vector change         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   on_root for diagnostics                                               !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout, NORMAL, VERBOSE
    use rundat, only: pub_geom_output_detail
    use simulation_cell, only: castep_cell_frac_to_cart

    implicit none
    real(kind=DP), intent(inout)                 :: lambda
    real(kind=DP), intent(in), dimension(1:ndim) :: delta_vec
    real(kind=DP), intent(in)                    :: scale

    !local vars
    real(kind=DP), parameter                     :: delta_ions_tol=0.5_dp !max allowed mag. change in Bohr (au)
    real(kind=DP)                                :: lambda_new
    real(kind=DP)                                :: delta_mag
    real(kind=DP), dimension(1:3)                :: delta_frac, delta_cart
    integer                                      :: i,j
    logical                                      :: warn_ions

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_check_lambda'

    !convert lambda*delta into 3D ionic change
    if (.not.fix_all_ions) then

       delta_mag=-1.0_dp
       do i=10,ndim,3     !check appropriate displacements
          do j=1,3
             delta_frac(j)=lambda*delta_vec(i+j-1)  !this is the proposed displacement wrt that at end of last iter
          end do                                    !NOT wrt last minimization !?

          !convert from fractional to cartesian
          call castep_cell_frac_to_cart(current_cell,delta_frac,delta_cart)

          !compare cartesian length to max allowed displacement (delta_ions_tol)
          delta_mag=max(delta_mag,sqrt(dot_product(delta_cart,delta_cart)))
       end do

       ! lk: changed output detail from > NORMAL to > VERBOSE.
       if (delta_mag>scale*delta_ions_tol) then
          !NB delta_mag >=0 but lambda may (if .not.monotonic_E) be < 0 so must be careful here
          lambda_new=lambda*scale*delta_ions_tol/delta_mag
          if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
             write (stdout,1) 'geom_check_lambda: Warning - BFGS step would &
                  &result in ionic displacement >',scale*delta_ions_tol,' au'
             write (stdout,2) 'geom_check_lambda: Warning - reducing lambda &
                  &from ',lambda,' to ',lambda_new,' for safety'
          end if
          lambda=lambda_new
       end if

       !warn user if overall ionic coords have changed by a lot
       warn_ions=.false.
       do i=10,ndim,3     !check appropriate displacements
          do j=1,3
             delta_frac(j)=lambda*delta_vec(i+j-1)  !this is the checked & scaled displacement
          end do
          if (maxval(delta_frac)>0.1_dp) warn_ions=.true.
       end do
       if (warn_ions .and. pub_on_root) then
          write (stdout,*) 'geom_check_lambda: Warning - fractional ionic displacement > 10%'
          write (stdout,*) 'geom_check_lambda:         - suggest you reconsider the initial ionic positions'
       end if

    end if

1   format(1x,a,f10.3,a)
2   format(1x,a,f10.6,a,f10.6,a)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_check_lambda'

    return
  end subroutine geom_check_lambda

  subroutine geom_inv_Hessian_initialize(toggle,init_cell, cmdl, mdl, reset_updates, init_method, force_opt_method)
    !=========================================================================!
    ! Initialize the inverse Hessian matrix for good convergence to gs ions   !
    ! as in B.G. Pfrommer et al, J.Comp.Phys. v131, p233-240 (1997).          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   toggle,    intent=in,  switch to determine symm break or not          !
    !   ^^^^^^^ -- jd (retreat2014): This is apparently ignored. I'm adding   !
    !                                a utils_use_var() to silence it, but     !
    !                                I'm not removing it in case it's more    !
    !                                than a vestige.                          !
    !   init_cell, intent=in,  the cell to derive inv_Hessian from            !
    !   cmdl,     intent=inout, model containing suitably initialized inv_Hess !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   parameters to define the bulk modulus and average phonon frequency    !
    !   io for error messages and unit conversion routines                    !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   init_cell must be valid                                               !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified for LBFGS by Jolyon Aarons, v2.0, 26/06/2011                   !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout, NORMAL, VERBOSE
    use linalg, only: linalg_invert_sym_matrix
    use model_type, only: MODEL
    use rundat, only: pub_geom_print_inv_hessian, pub_geom_output_detail
    use simulation_cell, only: castep_cell, castep_model
    use utils, only: utils_alloc_check, utils_assert, utils_use_var, utils_abort

    implicit none

    type(MODEL), intent(in)              :: mdl  ! jcap
    integer,            intent(in)    :: toggle
    type(castep_cell), intent(in)        :: init_cell
    type(castep_model), intent(inout)    :: cmdl
    logical, optional,  intent(in)    :: reset_updates
    integer, optional,  intent(in)    :: init_method
    character, optional,intent(in)    :: force_opt_method

    !local variables
    real(kind=DP), dimension(1:3,1:3) :: metric
    real(kind=DP) :: avg_mass
    real(kind=DP) :: inv_3VB
    integer :: i,j,k,l,ispec
    integer :: ierr
    logical :: loc_do_bfgs, loc_do_lbfgs

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_inv_Hessian_initialize'

    call utils_use_var(toggle)

    ! ja531 -> March 2015: Added the option to force initialising the l-BFGS stuff.
    if(present(force_opt_method)) then
       select case(force_opt_method)
       case("B")
          loc_do_bfgs=.true.
          loc_do_lbfgs=.false.
       case("L")
          loc_do_lbfgs=.true.
          loc_do_bfgs=.false.
       case default
          call utils_abort("Unrecognised force_opt_method in geom_inv_Hessian_initialize : "//force_opt_method)
       end select
    else
       loc_do_bfgs=do_bfgs
       loc_do_lbfgs=do_lbfgs
    end if

    if (.not. pub_on_root) return

    !Initialize most of the matrix with a lot of zeros
    if(loc_do_bfgs) then
    cmdl%bfgs_inv_Hessian=0.0_dp
    else if(loc_do_lbfgs) then
       lbfgs_init_inv_hess=0.0_dp
       if(present(reset_updates)) then
          if(reset_updates) then
             cmdl%lbfgs_num_updates = 0

             !   Initialise correction pairs
             cmdl%lbfgs_position_updates = 0.0_dp
             cmdl%lbfgs_gradient_updates = 0.0_dp

             !   Initialise intermediate matrices.
             YY_mat = 0.0_dp
             D_mat  = 0.0_dp
             R_mat  = 0.0_dp
          end if
       end if
       if(present(init_method)) then
          call utils_assert(init_method.ne.hess_init_householder,&
               &'Can not initialise with Householder transform in geom_inv_Hessian_initialize')

          if((init_method.ge.hess_init_identity).and.(init_method.le.hess_init_cellion)) then
             lbfgs_init_inv_hess=0.0_dp
             lbfgs_init_inv_hess(3,1:9)=abs(init_scalar(1))*debug_init_scale
             lbfgs_init_inv_hess(3,10:ndim)=abs(init_scalar(2))*debug_init_scale
             ! lk: show these only if output_detail > VERBOSE.
             if ((pub_geom_output_detail > VERBOSE) .and. pub_on_root) then
                write(stdout,*) "BFGS : (re)initialising inverse Hessian &
                                 &with Identity matrix"
                write(stdout,*) "BFGS : scaling identity initial &
                                 &inverse Hessian by : ", &
                                 abs(init_scalar(1))*debug_init_scale
                write(stdout,*) "BFGS : ", &
                                 abs(init_scalar(1))*debug_init_scale! &
                                 !,pressure_unit), &trim(pressure_label)
                write(stdout,*) "BFGS : ", &
                                 abs(init_scalar(1))*debug_init_scale! &
                                 !,frequency_unit), &trim(frequency_label)
             end if

             if (pub_debug_on_root) write(stdout,'(a)') &
                  'DEBUG: Leaving geom_inv_Hessian_initialize'
             return
          end if
       end if
    end if

    !Initialize stress/strain part of inv_Hessian using bulk modulus
    !first 9 diagonal elements=1/(3*Vol*B)
    inv_3VB=1.0_dp/(3.0_dp*init_cell%volume*new_modulus_est)
    if(loc_do_bfgs) then
    do i=1,9
       cmdl%bfgs_inv_Hessian(i,i)=inv_3VB
    end do
    else if(loc_do_lbfgs) then
       lbfgs_init_inv_hess(3,1:9)=inv_3VB
    endif

    !Initialize ionic part of inv_Hessian using average phonon frequency
    !3x3 blocks=(metric^-1)/(mass*omega^2), where metric=h^T.h
    metric=matmul(init_cell%real_lattice,transpose(init_cell%real_lattice))
!    call algor_invert(3,metric)       !metric -> metric^-1
    call linalg_invert_sym_matrix(metric,3)

    avg_mass=0.0_dp
    do ispec=1,init_cell%num_species - mdl%nat_classical
       avg_mass=avg_mass+init_cell%species_mass(ispec)*init_cell%num_ions_in_species(ispec)
    end do

    ! ndmh: 26/01/09 very important fix: atomic units of mass are m_e so
    ! ndmh: multiply by (1 a.m.u.) / (1 m_e) = 1822.888485
    avg_mass=1822.888485_dp*avg_mass/init_cell%num_ions

    if (pub_geom_output_detail>NORMAL .and. pub_on_root) then
       write(stdout,*) 'avg_mass= ',avg_mass
    end if

    !(metric^-1)/(mass*omega^2)
    metric=metric/(avg_mass*new_frequency_est**2)

    if(loc_do_bfgs) then
       !do i=1,init_cell%num_ions
    do i=1,init_cell%num_ions - mdl%nat_classical
       k=10+3*(i-1)
       cmdl%bfgs_inv_Hessian(k:k+2,k:k+2)=metric(1:3,1:3)
    end do
    else if(loc_do_lbfgs) then
       lbfgs_init_inv_hess(2:3,11) = metric(1:2,2)
       lbfgs_init_inv_hess(1,11) = 0.0_dp
       lbfgs_init_inv_hess(3,10) = metric(1,1)
       lbfgs_init_inv_hess(1:2,10) = 0.0_dp
       lbfgs_init_inv_hess(1:3,12) = metric(1:3,3)
       do i=2,init_cell%num_ions
          k=10+3*(i-1)
          lbfgs_init_inv_hess(1:3,k:k+2) = lbfgs_init_inv_hess(1:3,10:12)
       end do

       !       if(want_lbfgs_inverse) then
       if(.not.allocated(lbfgs_init_hess)) then
          allocate(lbfgs_init_hess(1:3,1:ndim),stat=ierr)
          call utils_alloc_check('geom_inv_Hessian_initialise', &
                                 'lbfgs_init_hess',ierr)
       end if
       lbfgs_init_hess=0.0_dp
       inv_3VB=(3.0_dp*init_cell%volume*new_modulus_est)
       ! call comms_bcast(pub_root_proc_id,inv_3VB) !parallel paranoia
       lbfgs_init_hess(3,1:9)=inv_3VB
       call linalg_invert_sym_matrix(metric,3)
       lbfgs_init_hess(2:3,11) = metric(1:2,2)
       lbfgs_init_hess(1,11) = 0.0_dp
       lbfgs_init_hess(3,10) = metric(1,1)
       lbfgs_init_hess(1:2,10) = 0.0_dp
       lbfgs_init_hess(1:3,12) = metric(1:3,3)
       do i=2,init_cell%num_ions
          k=10+3*(i-1)
          lbfgs_init_hess(1:3,k:k+2)=lbfgs_init_hess(1:3,10:12)
       end do
       ! end if
    end if

    !Finally, if symmetry is not to be strictly enforced, then we may add terms
    !to mix forces & stresses so ensure symmetry breaking is possible if required
    !NB We add these terms asymmetrically - we want to break symm for ions but not cell!


    if (pub_geom_print_inv_hessian.and.loc_do_bfgs.and.pub_on_root) then
       l=ndim/5
       do i=1,ndim
          do j=1,5*l,5
             write (stdout,99) (i,j+k,cmdl%bfgs_inv_Hessian(i,j+k:j+k),k=0,4)
          end do
          if (5*l<ndim) write (stdout,99) (i,j,cmdl%bfgs_inv_Hessian(i,j:j),j=5*l,ndim)
       end do
    else if(loc_do_lbfgs) then
       !           if (on_root.and.(verbose.or.iprint>1)) then
       !              write (stdout,*) 'LBFGS: WARNING - adding off-diagonal terms into inv_Hessian to break symmetry'
       !              write (stdout,*) '                 was not possible in this version of LBFGS.'
                      ! need general sparse form for initial H_0.
       !           end if
    end if

99  format('inv_Hessian:',5('(',i2,',',i2,')=',f15.10,:))

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_inv_Hessian_initialize'

    return
  end subroutine geom_inv_Hessian_initialize

  subroutine geom_update_inv_Hessian(old_x_vec,x_vec,old_f_vec,f_vec,init_cell,cmdl,mdl)
    !=========================================================================!
    ! Update the inverse Hessian matrix using the standard BFGS algorithm     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   inv_Hessian, intent=inout, BFGS updated (Hessian)^-1                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !   reset_iter and previous_reset_iter may be updated                     !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified for LBFGS by Jolyon Aarons, v2.0, 26/06/2011                   !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout, NORMAL, VERBOSE
    use model_type, only: MODEL ! jcap
    use rundat, only: pub_geom_output_detail, pub_geom_lbfgs_max_updates
    use simulation_cell, only: castep_cell, castep_model, castep_cell_copy
    use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    type(MODEL),                             intent(in)    :: mdl
    real(kind=DP), dimension(1:ndim),        intent(in)    :: old_x_vec
    real(kind=DP), dimension(1:ndim),        intent(in)    :: x_vec
    real(kind=DP), dimension(1:ndim),        intent(in)    :: old_f_vec
    real(kind=DP), dimension(1:ndim),        intent(in)    :: f_vec
    type(castep_cell),                       intent(inout) :: init_cell
    type(castep_model),                      intent(inout) :: cmdl

    !local variables
    real (kind=dp), dimension(:), allocatable :: dx_vec, df_vec, u_vec, h_df_vec, dx_vec_remove, df_vec_remove
    real (kind=dp)                    :: df_dot_h_df, dx_dot_df, dx_dot_df_remove
    real (kind=dp)                    :: symm_tmp, symm_err
    integer                           :: i,j
    integer                           :: ierr
    logical                           :: dx_df_bad, debug_doit

    ! BLAS subroutines
    external :: dgemv, dsbmv

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_update_inv_Hessian'

    allocate(dx_vec(1:ndim),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian','dx_vec',ierr)
    allocate(df_vec(1:ndim),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian','df_vec',ierr)

    dx_vec=x_vec-old_x_vec
    df_vec=f_vec-old_f_vec
    dx_dot_df=dot_product(dx_vec,df_vec)

    if (.not.pub_on_root) return
    !check norms
    if (pub_geom_output_detail>NORMAL) then
       write (stdout,*) trim(geom_string)//': norms of important vectors in update_inv_Hessian'
       write (stdout,*) ' ||dx||1  =',sum(abs(dx_vec))
       write (stdout,*) ' ||dx||2  =',sqrt(dot_product(dx_vec,dx_vec))
       write (stdout,*) ' ||dx||inf=',maxval(abs(dx_vec))
       write (stdout,*) ' ||df||1  =',sum(abs(df_vec))
       write (stdout,*) ' ||df||2  =',sqrt(dot_product(df_vec,df_vec))
       write (stdout,*) ' ||df||inf=',maxval(abs(df_vec))
       write (stdout,*) '  dx.df   =',dx_dot_df
    end if

    !parallel synch
    dx_df_bad=(dx_dot_df>(-1.0_dp*tiny(dx_dot_df)))
    !    call comms_gcopy(dx_df_bad,1)

    !diagnostic
    if (pub_geom_output_detail>NORMAL) then
       if (dx_df_bad) then
          write (stdout,*) trim(geom_string)//': dx.df=',dx_dot_df,' NOT OK!'
       else
          write (stdout,*) trim(geom_string)//': dx.df=',dx_dot_df,' OK'
       end if
    end if

    if(bulk_freq_progress) call geom_update_inv_Hessian_params(init_cell,cmdl,mdl,printonly=.true.)

    !To ensure the inv_Hessian stays positive-definite we may need to reset it rather than update it
    !if we get too much error accumulation due to non-exact line minimisation,
    !which will happen after too many steps or a bad update, so:

    if(debug_reinit.gt.0) then
       debug_doit=(mod(iteration,debug_reinit).eq.0)
    else
       debug_doit=.false.
    end if

    if ((((iteration-reset_iter)>ndof).and.do_reinit).or.dx_df_bad.or.debug_doit) then   !dx.df ought to be <0
       if (pub_geom_output_detail>NORMAL) then
          write (stdout,*) trim(geom_string)//': resetting inv_Hessian to prevent error accumulation'
          if ((iteration-reset_iter)>ndof) write (stdout,*) '      due to number of steps taken'
          if (dx_df_bad)                            write (stdout,*) '      due to bad update'
       end if
       !       if(cmdl%lbfgs_num_updates.gt.0) then !moved into routine
       if(hess_init_method.lt.3) then
       call geom_update_inv_Hessian_params(init_cell,cmdl,mdl)
       end if
       !       end if
       call castep_cell_copy(cmdl%cell,init_cell)
       !       call lbfgs_init(toggle_on,init_cell,cmdl)
       !     if(on_root) then
       if(do_lbfgs) then
          if (hess_init_method.eq.hess_init_scaled) then
             init_scalar=(-dx_dot_df/dot_product(df_vec,df_vec))!*0.1_dp
          else if (hess_init_method.eq.hess_init_cellion) then
             init_scalar(1)=(-dot_product(dx_vec(1:9),df_vec(1:9))/dot_product(df_vec(1:9),df_vec(1:9)))
             init_scalar(2)=(-dot_product(dx_vec(10:ndim),df_vec(10:ndim))/dot_product(df_vec(10:ndim),df_vec(10:ndim)))
          end if
       end if
       !     end if
       call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,reset_updates=.true.,init_method=hess_init_method)
       previous_reset_iter=reset_iter
       reset_iter=iteration
       !     call trace_exit('geom_update_inv_Hessian',status)
       return
    end if


    allocate(h_df_vec(1:ndim),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian','h_df_vec',ierr)

    !All is OK so we can continue to do the update:
    h_df_vec=geom_apply_vec(cmdl,df_vec,'y')
    !    call comms_gcopy(h_df_vec(1),ndim) !parallel paranoia
    df_dot_h_df=dot_product(df_vec,h_df_vec)

    !this ought to be impossible with the above positive-definite test, but for safety:
    call utils_assert(abs(dx_dot_df) >= tiny(dx_dot_df), &
          'Error in geometry_optimise: in geom_update_inv_Hessian - dx_dot_df->0')
    call utils_assert(abs(df_dot_h_df) >= tiny(df_dot_h_df), &
          'Error in geometry_optimise: in geom_update_inv_Hessian - df_dot_h_df->0')

    allocate(u_vec(1:ndim),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian','u_vec',ierr)
    if(do_bfgs) then
    u_vec=(dx_vec/dx_dot_df) - (h_df_vec/df_dot_h_df)

       !check norms
       ! lk: changed output detail from > NORMAL to > VERBOSE.
       if (pub_geom_output_detail>VERBOSE) then
          write (stdout,*) trim(geom_string)//': norms of important vectors in update_inv_Hessian'
          write (stdout,*) '||Hdf||1  =',sum(abs(h_df_vec))
          write (stdout,*) '||Hdf||2  =',sqrt(dot_product(h_df_vec,h_df_vec))
          write (stdout,*) '||Hdf||inf=',maxval(abs(h_df_vec))
          write (stdout,*) ' ||u||1   =',sum(abs(u_vec))
          write (stdout,*) ' ||u||2   =',sqrt(dot_product(u_vec,u_vec))
          write (stdout,*) ' ||u||inf =',maxval(abs(u_vec))
          write (stdout,*) ' df.Hdf   =',df_dot_h_df
       end if

       do i=1,ndim
          do j=1,ndim
             cmdl%bfgs_inv_Hessian(j,i)=cmdl%bfgs_inv_Hessian(j,i) - dx_vec(i)*dx_vec(j)/dx_dot_df &
                  &                   - h_df_vec(i)*h_df_vec(j)/df_dot_h_df &
                  &                   + df_dot_h_df*u_vec(i)*u_vec(j)
          end do
       end do
       !not a good idea to gcopy a 3N*3N array ...

       !force update to be symmetric
       symm_err=0.0_dp
       do i=1,ndim
          do j=i+1,ndim
             symm_tmp=(cmdl%bfgs_inv_Hessian(j,i)+cmdl%bfgs_inv_Hessian(i,j))/2.0_dp
             symm_err=max(symm_err,abs((cmdl%bfgs_inv_Hessian(j,i)-cmdl%bfgs_inv_Hessian(i,j))/2.0_dp))
             cmdl%bfgs_inv_Hessian(j,i)=symm_tmp
             cmdl%bfgs_inv_Hessian(i,j)=symm_tmp
          end do
       end do
       if (pub_geom_output_detail>NORMAL) then
          if(do_bfgs) write (stdout,10) symm_err
          if(do_lbfgs) write (stdout,11) symm_err
       end if
10     format(1x,'BFGS: max inv_Hessian symmetry error',es17.8e3)
11     format(1x,'LBFGS: max inv_Hessian symmetry error',es17.8e3)

    else if(do_lbfgs) then
       df_vec=-df_vec
       dx_dot_df=-dx_dot_df

       !      if(hess_init_method.eq.hess_init_scaled) then
       !         if(iteration.eq.0) then
       !            lbfgs_init_inv_hess=lbfgs_init_inv_hess*(dx_dot_df/dot_product(df_vec,df_vec))
       !         end if
       !      end if

       ! Store the new correction pair (s,y). There are two possible
       ! situations: we haven't yet used up the maximum allowable storage, or
       ! we have used up all the storage, in which case we have to remove one
       ! of the correction pairs in order to store the new one. Note that a
       ! new correction pair is always placed at the end of the matrix.
       if(cmdl%lbfgs_num_updates.ne.0) then
          if (cmdl%lbfgs_num_updates == pub_geom_lbfgs_max_updates) then
             if(update_init_on_remove) then
                allocate(dx_vec_remove(1:ndim),stat=ierr)
                call utils_alloc_check('geom_update_inv_Hessian','dx_vec_remove',ierr)
                allocate(df_vec_remove(1:ndim),stat=ierr)
                call utils_alloc_check('geom_update_inv_Hessian','df_vec_remove',ierr)
                dx_vec_remove=cmdl%lbfgs_position_updates(:,1)-old_x_vec_remove !s
                df_vec_remove=cmdl%lbfgs_gradient_updates(:,1)-old_f_vec_remove !y
                dx_dot_df_remove=dot_product(dx_vec_remove,df_vec_remove) !sy
                init_scalar=(dx_dot_df_remove/dot_product(df_vec_remove,df_vec_remove))
                old_x_vec_remove=cmdl%lbfgs_position_updates(:,1)
                old_f_vec_remove=cmdl%lbfgs_gradient_updates(:,1)
                deallocate(dx_vec_remove,stat=ierr)
                call utils_dealloc_check('geom_update_inv_hessian','dx_vec_remove', ierr)
                deallocate(df_vec_remove,stat=ierr)
                call utils_dealloc_check('geom_update_inv_hessian','df_vec_remove', ierr)

                call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,reset_updates=.false.,init_method=hess_init_scaled)

                cmdl%lbfgs_num_updates = cmdl%lbfgs_num_updates - 1

                ! Remove the first column from S and Y.
                cmdl%lbfgs_position_updates = eoshift(cmdl%lbfgs_position_updates, shift=1, dim=2)
                cmdl%lbfgs_gradient_updates = eoshift(cmdl%lbfgs_gradient_updates, shift=1, dim=2)

                ! Update R = S'*Y.
                R_mat=0.0_dp
                do i=1,cmdl%lbfgs_num_updates
                   call dgemv('t',ndim,i,1.0_dp,cmdl%lbfgs_position_updates(1,1),&
                        &size(cmdl%lbfgs_position_updates,1),cmdl%lbfgs_gradient_updates(1,i),1,0.0_dp,R_mat(1,i),1)
          end do

                YY_mat(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)=0.0_dp
                call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_inv_hess,size(lbfgs_init_inv_hess,1),&
                     &cmdl%lbfgs_gradient_updates(1,1),1,0.0_dp,tmpn,1)
                YY_mat(1,1)=dot_product(tmpn,cmdl%lbfgs_gradient_updates(1:ndim,1))+R_mat(1,1)

                do i=2,cmdl%lbfgs_num_updates
                   call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_inv_hess,size(lbfgs_init_inv_hess,1),&
                        &cmdl%lbfgs_gradient_updates(1,i),1,0.0_dp,tmpn,1)
                   call dgemv('t',ndim,i-1,1.0_dp,cmdl%lbfgs_gradient_updates,&
                        &size(cmdl%lbfgs_gradient_updates,1),tmpn,1,0.0_dp,YY_mat(1,i),1)
                   YY_mat(i,1:i-1)=YY_mat(1:i-1,i)
                   YY_mat(i,i)=dot_product(tmpn,cmdl%lbfgs_gradient_updates(1:ndim,i)) + R_mat(i,i)
       end do
             else
                call lbfgs_removepair(cmdl)
             end if
          else if (mod(cmdl%lbfgs_num_updates,lbfgs_block_length).eq.0) then
             call lbfgs_resize(cmdl,cmdl%lbfgs_num_updates+lbfgs_block_length)
          end if
       end if
       cmdl%lbfgs_num_updates = cmdl%lbfgs_num_updates + 1
       ! Update D_mat.
       D_mat(cmdl%lbfgs_num_updates) = dx_dot_df
       ! Update R = S'*Y.
       R_mat(cmdl%lbfgs_num_updates,:) = 0.0_dp
       !      R_mat(1:cmdl%lbfgs_num_updates-1,cmdl%lbfgs_num_updates) = matmul(cmdl%S,y)
       call dgemv('t',ndim,cmdl%lbfgs_num_updates-1,1.0_dp,cmdl%lbfgs_position_updates,&
            &size(cmdl%lbfgs_position_updates,1),df_vec,1,0.0_dp,R_mat(1,cmdl%lbfgs_num_updates),1)
       R_mat(cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates) = dx_dot_df

       ! Update the collection of correction pairs (dx_vec,df_vec).
       cmdl%lbfgs_position_updates(1:ndim,cmdl%lbfgs_num_updates) = dx_vec(1:ndim)
       cmdl%lbfgs_gradient_updates(1:ndim,cmdl%lbfgs_num_updates) = df_vec(1:ndim)

       ! Update YY_mat_mat = Y'*Y
       u_vec=0.0_dp
       call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_inv_hess,size(lbfgs_init_inv_hess,1),&
            &df_vec,1,0.0_dp,u_vec,1)
       if(cmdl%lbfgs_num_updates.gt.1) then
          call dgemv('t',ndim,cmdl%lbfgs_num_updates-1,1.0_dp,cmdl%lbfgs_gradient_updates,&
               &size(cmdl%lbfgs_gradient_updates,1),u_vec,1,0.0_dp,YY_mat(1,cmdl%lbfgs_num_updates),1)
          YY_mat(cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates-1)=YY_mat(1:cmdl%lbfgs_num_updates-1,cmdl%lbfgs_num_updates)
       end if
       YY_mat(cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates)=dot_product(u_vec,df_vec) + dx_dot_df
    end if

    deallocate(dx_vec,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian','dx_vec',ierr)
    deallocate(df_vec,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian','df_vec',ierr)
    deallocate(u_vec,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian','u_vec',ierr)
    deallocate(h_df_vec,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian','h_df_vec',ierr)

! jd: Not referenced anywhere:
! 99  format('inv_Hessian:',5('(',i2,',',i2,')=',f15.10,:))

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_update_inv_Hessian'

    return
  end subroutine geom_update_inv_Hessian

  subroutine geom_update_inv_hessian_params(init_cell,cmdl,mdl,final,printonly)
    !=========================================================================!
    ! Update the parameters used to initialize the inverse Hessian matrix,    !
    ! i.e. geom_modulus_est and geom_frequency_est, by analysing the current  !
    ! inv_Hessian.                                                            !
    ! This follows Pfrommer section 3 closely and uses same notation.         !
    !                                                                         !
    ! NB Cannot call this routine WITHOUT then reseting inv_Hessian as it will!
    !    update geom_modulus_est and/or geom_frequency_est which on the next  !
    !    call will then be out of sync with what has actually been used!      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   init_cell,   intent=in, cell when inv_Hessian last initialized        !
    !   cmdl,         intent=in, current state of the system                   !
    !   final,       intent=in, final call of analysis => dont print results  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   geom_modulus_est & geom_frequency_est in parameters are recalculated  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Jolyon Aarons, v3.0, 30/06/2011, based on v2.0:              !
    ! Written by Matt Probert, v2.0, 21/02/2007                               !
    ! v2.0 now includes Householder transformation s.t. can handle mixed case !
    ! v3.0 now outsources the construction of a minimum spanning basis to     !
    ! geom_lbfgs_construct_basis, which calculates an augmented basis by      !
    ! combining that calculated for the cell and ion spaces independently in  !
    ! such a way as to span the whole space... instead of a Householder for   !
    ! LBFGS, in order to keep memory scaling linear.                          !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!
    use constants, only : dp, stdout, NORMAL, VERBOSE
    use comms, only: pub_on_root, comms_bcast !, pub_root_proc_id
    use linalg, only: linalg_invert_sym_matrix
    use model_type, only: MODEL ! jcap
    use rundat, only: pub_geom_output_detail
    use simulation_cell, only: castep_cell, castep_model
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    type(MODEL),             intent(in)    :: mdl
    type (castep_cell),  intent(inout) :: init_cell
    type(castep_model),           intent(inout) :: cmdl
    logical, optional,            intent(in)    :: final
    logical, optional, intent(in)    :: printonly

    !local variables
    real (kind=dp), dimension(:,:), allocatable :: tmp_Hess_matrix
    real (kind=dp), dimension(1:3,1:3)          :: tmp_strain, H0_block

    real (kind=dp), dimension(:,:), allocatable :: Abar        !Hessian in reduced coords, Pfrommer eqn 16

    real (kind=dp), dimension(:,:), allocatable :: tmpgradpos
    real (kind=dp), dimension(:,:), allocatable :: tmpmm
    real (kind=dp), dimension(:,:), allocatable :: R_inv
    real (kind=dp), dimension(:), allocatable   :: mass_vector
    real (kind=dp), dimension(:,:), allocatable :: basis_vecs
    integer                                     :: num_basis_vecs
    integer                                     :: num_cell_basis_vecs
    integer                                     :: num_ion_basis_vecs
    character(len=5)                            :: bfgs_method_local

    logical                                     :: update_params

    integer :: i,iatom,ispec
    integer :: ierr
    integer :: ndim3

    integer :: cell_size=9

    real (kind=dp) :: bulk, phonon
    real (kind=dp) :: dum_max,dum_min,max_B
    logical        :: update_OK, local_final

    !LAPACK requires
    integer :: info

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_update_inv_Hessian_params'

    update_params=.true.
    if(present(printonly)) then
       update_params=.not.printonly
    end if

    if(hess_init_method.eq.hess_init_identity) then
!       call trace_exit('geom_update_inv_hessian_params',status)
       return
    end if

    if (present(final)) then           !don't print analysis results here if final call
       local_final=final
    else
       local_final=.false.
    end if

    !initialize symmetric array size
    ndim3=ndim-3
    ! qoh: Initialize to avoid compiler warning
    update_ok = .false.

    !First get I+strain from init_cell to current_cell
    tmp_strain=init_cell%real_lattice
    call linalg_invert_sym_matrix(tmp_strain,3) !using tmp_strain as temp ref_cell^-1
    !h'=(1+e)h0 => (h^T)'=(h0^T)(1+e)^T => (1+e)^T=(h0^T)^-1.(h^T)'
    !tmp_strain=transpose(matmul(tmp_strain,cmdl%cell%real_lattice))
#ifdef PGI_BUG
    unit_matrix=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
#endif
    tmp_strain = unit_matrix

    info=0
    if(do_bfgs) then
       bfgs_method_local="BFGS"
       !now calculate original inv_Hessian in Hp_matrix
    !(backup current inv_Hessian in tmp_Hess_matrix)
    !NB must NOT have changed new_frequency_est or new_modulus_est between the 'real' call in main program
    !to get inv_Hessian and this call!
       allocate(tmp_Hess_matrix(ndim,ndim),stat=ierr)
       call utils_alloc_check('geom_update_inv_hessian_params','tmp_Hess_matrix',ierr)
    tmp_Hess_matrix=cmdl%bfgs_inv_Hessian
       tmp_Hess_matrix=cmdl%bfgs_inv_Hessian
       !ja531: with LBFGS we cannot do toggle_on/off stuff
       !       so for testing scaled identity we also do not do it
       !       call geom_inv_Hessian_initialize(toggle_off,init_cell,cmdl)
       if((hess_init_method.eq.hess_init_identity).or.(hess_init_method.eq.hess_init_scaled)) then
          !          call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,method=hess_init_identity)
          ! ja531 --> 16 August 2011 : toggle_off?
          call geom_inv_Hessian_initialize(toggle_off,init_cell,cmdl,mdl,init_method=hess_init_identity)
       else
          !          call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl)
          call geom_inv_Hessian_initialize(toggle_off,init_cell,cmdl,mdl)
       end if
       !       Hp_matrix=cmdl%bfgs_inv_Hessian
       !       cmdl%bfgs_inv_Hessian=tmp_Hess_matrix
       tmp_Hess_matrix(1:3,1:3) = matmul(transpose(tmp_strain),matmul(tmp_Hess_matrix(1:3,1:3),tmp_strain))
       tmp_Hess_matrix(4:6,4:6) = matmul(transpose(tmp_strain),matmul(tmp_Hess_matrix(4:6,4:6),tmp_strain))
       tmp_Hess_matrix(7:9,7:9) = matmul(transpose(tmp_strain),matmul(tmp_Hess_matrix(7:9,7:9),tmp_strain))
       cmdl%bfgs_inv_Hessian(1:3,1:3) = matmul(transpose(tmp_strain),matmul(cmdl%bfgs_inv_Hessian(1:3,1:3),tmp_strain))
       cmdl%bfgs_inv_Hessian(4:6,4:6) = matmul(transpose(tmp_strain),matmul(cmdl%bfgs_inv_Hessian(4:6,4:6),tmp_strain))
       cmdl%bfgs_inv_Hessian(7:9,7:9) = matmul(transpose(tmp_strain),matmul(cmdl%bfgs_inv_Hessian(7:9,7:9),tmp_strain))

       cmdl%bfgs_inv_Hessian=tmp_Hess_matrix-cmdl%bfgs_inv_Hessian
       !cmdl%bfgs_inv_Hessian is now the update

    else if(do_lbfgs) then
       if(cmdl%lbfgs_num_updates.lt.1) then
!          call trace_exit('geom_update_inv_hessian_params',status)
          return
       end if

       bfgs_method_local="LBFGS"
       allocate(tmpmm(2*cmdl%lbfgs_num_updates,2*cmdl%lbfgs_num_updates),stat=ierr)
       call utils_alloc_check('geom_update_inv_hessian_params','tmpmm',ierr)
       tmpmm=0.0_dp
       allocate(R_inv(cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates),stat=ierr)
       call utils_alloc_check('geom_update_inv_hessian_params','R_inv',ierr)
       R_inv=R_mat(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)
       call linalg_invert_sym_matrix(R_inv,cmdl%lbfgs_num_updates)
       tmpmm(1:cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)=&
            &-R_inv
       deallocate(R_inv,stat=ierr)
       call utils_dealloc_check('geom_update_inv_hessian_params','R_inv',ierr)
       tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,&
            &1:cmdl%lbfgs_num_updates)= transpose(tmpmm(1:cmdl%lbfgs_num_updates,&
            &cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2))
       tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)=&
            &matmul(matmul(tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,&
            &1:cmdl%lbfgs_num_updates),YY_mat(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)),&
            &tmpmm(1:cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2))


       allocate(tmpgradpos(ndim,2*cmdl%lbfgs_num_updates),stat=ierr)
       call utils_alloc_check('geom_update_inv_hessian_params','tmpgradpos',ierr)
       call lbfgs_initH_mul('n',cmdl%lbfgs_gradient_updates(:,1:cmdl%lbfgs_num_updates),lbfgs_init_inv_hess,&
            &tmpgradpos(:,1:cmdl%lbfgs_num_updates),info)
       tmpgradpos(:,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)=cmdl%lbfgs_position_updates(:,1:cmdl%lbfgs_num_updates)

       tmpgradpos(1:3,1:cmdl%lbfgs_num_updates*2)=matmul(transpose(tmp_strain),tmpgradpos(1:3,1:cmdl%lbfgs_num_updates*2))
       tmpgradpos(4:6,1:cmdl%lbfgs_num_updates*2)=matmul(transpose(tmp_strain),tmpgradpos(4:6,1:cmdl%lbfgs_num_updates*2))
       tmpgradpos(7:9,1:cmdl%lbfgs_num_updates*2)=matmul(transpose(tmp_strain),tmpgradpos(7:9,1:cmdl%lbfgs_num_updates*2))

       do i=1,9,3
          H0_block=0.0_dp
          H0_block(:,3)=lbfgs_init_inv_hess(1:3,i+2)
          H0_block(1:2,2)=lbfgs_init_inv_hess(2:3,i+1)
          H0_block(1,1)=lbfgs_init_inv_hess(3,i)
          H0_block(2:3,1)=H0_block(1,2:3)
          H0_block(3,2)=H0_block(2,3)

          H0_block = matmul(transpose(tmp_strain),matmul(H0_block,tmp_strain))
          !       H0_block = matmul(transpose(D_matrix(i:i+2,i:i+2)),matmul(H0_block,D_matrix(i:i+2,i:i+2)))

          lbfgs_init_inv_hess(1:3,i+2)=H0_block(:,3)
          lbfgs_init_inv_hess(2:3,i+1)=H0_block(1:2,2)
          lbfgs_init_inv_hess(3,i)=H0_block(1,1)
          lbfgs_init_inv_hess(1:2,i)=0.0_dp
          lbfgs_init_inv_hess(1,i+1)=0.0_dp
       end do

    end if

    if((hess_init_method.eq.hess_init_householder).and.do_lbfgs) then
       call lbfgs_update_bandmatrix(cmdl,tmpgradpos,tmpmm,lbfgs_init_inv_hess)
!       call trace_exit('geom_update_inv_hessian_params',status)
       return
    end if

    ! Create a vector of the masses of the atoms.
    allocate(mass_vector(1:ndim),stat=ierr)
    call utils_alloc_check('geom_update_inv_hessian_params','mass_vector',ierr)
    mass_vector=0.0_dp
    mass_vector(1:9)=1.0_dp
    ispec=1
    iatom=1
    do i=10,ndim,3
       ! VCA: begin alteration
       mass_vector(i)  =sqrt(cmdl%cell%species_mass(ispec))
       mass_vector(i+1)=sqrt(cmdl%cell%species_mass(ispec))
       mass_vector(i+2)=sqrt(cmdl%cell%species_mass(ispec))
       ! VCA: end alteration
       if (iatom<cmdl%cell%num_ions_in_species(ispec)) then
          iatom=iatom+1
       else
          iatom=1
          ispec=ispec+1
       end if
    end do

    call geom_calculate_update_basis(trim(bfgs_method_local),ndim,cell_size,num_cell_basis_vecs,num_ion_basis_vecs)

    ! ja531 --> 2011/09/14 01:17
    !           Fix suggested by Alexander Perlov -> num_basis_vecs needs to be set before printing the matrix.
    num_basis_vecs=num_cell_basis_vecs+num_ion_basis_vecs

    ! ja531 --> 2011/09/14 01:19
    !           Fix suggested by Alexander Perlov -> invert in geom_calculate_abar will crash if num_basis_vecs is zero.
    !           Shouldn't really happen anyway, since we expect at least one vector out from the SVD... further investigation?
    if(num_basis_vecs.gt.0) then

       allocate(Abar(1:num_basis_vecs,1:num_basis_vecs),stat=ierr)
       call utils_alloc_check('geom_update_inv_hessian_params','Abar',ierr)
       ! calculate the Hessian in the minimum-spanning space (A_bar).
       call geom_calculate_abar(trim(bfgs_method_local),ndim, num_basis_vecs)

       if(allocated(tmp_Hess_matrix)) then
          deallocate(tmp_Hess_matrix,stat=ierr)
          call utils_dealloc_check('geom_update_inv_hessian_params','tmp_Hess_matrix',ierr)
       end if

       ! Extract a bulk modulus (bulk) and average optical phonon frequency at the gamma point (phonon) from the
       ! reduced Hessian, A_bar in atomic units.
       call geom_calc_bulk_freq(ndim,cell_size,num_cell_basis_vecs,num_ion_basis_vecs,Abar,transpose(cmdl%cell%real_lattice),&
            &basis_vecs,bulk,phonon)

          end if

    !       call comms_bcast(pub_root_proc_id,phonon)
    !       call comms_bcast(pub_root_proc_id,bulk)

    if (num_ion_basis_vecs>0) then
       !now check phonon modes ...
       !finally, calculate <frequency> from generalized eigenvalues - is it really worth all this hassle???
       !NB Eigenvalues come from atomic units, where omega is an energy but when IO converts to/from energy it uses
       !   E=hbar.omega and so no two_pi factors here
       !NB We only update new_frequency_est if all is OK: freq(H2 vib)=132 THz so accept (0<freq<396 THz)
       if (pub_on_root) then
          dum_max=(396.0_dp*1.51982984550e-4_dp)   !no longer **2 as calc_bulk_freq does sqrt
          dum_min=(0.001_dp*1.51982984550e-4_dp)   !no longer **2 as calc_bulk_freq does sqrt
       end if
       !       call comms_gcopy(dummy,1)

       !       update_OK=((info==0) .and. (minval(omega2_vec)>0.0_dp) .and. (maxval(omega2_vec)<dummy))
       update_OK=((phonon>dum_min) .and. (phonon<dum_max))
       !       call comms_gcopy(update_OK,1)

       if (update_OK.and.update_params) then
          !          new_frequency_est=0.0_dp
          new_frequency_est=phonon
          !          do i=1,num_ion_basis_vecs
          !             new_frequency_est=new_frequency_est+sqrt(omega2_vec(i))
          !          end do
          !          new_frequency_est=new_frequency_est/num_ion_basis_vecs
          if (pub_on_root.and..not.local_final) write (stdout,12) trim(geom_string)//': updated estimated <frequency>  =', &
               & new_frequency_est*219474.6314431_dp,trim("cm-1")
       elseif(.not.update_params) then
          write (stdout,*) 'geom_update_inv_Hessian_params: phonon = ',phonon*219474.6314431_dp,trim("cm-1")
       else
          if ((pub_geom_output_detail>VERBOSE).and.pub_on_root) then
             write (stdout,*) 'geom_update_inv_Hessian_param: skipping <frequency> update because'
             if (info/=0) write (stdout,*)                   '                               dsygv failed to converge'
             if (phonon<=dum_min) write (stdout,*) '                               min omega**2<0'
             if (phonon>=dum_max) write (stdout,*) '                               max omega**2 too high'
          end if
       end if
    end if
    if (num_cell_basis_vecs>0) then
       !now check bulk modulus

       !NB we only update new_modulus_est if all is OK (0<B<3000 GPa)
       if (pub_on_root) then
          max_B=(3000.0_dp/29421.912_dp)!,'GPa')           !only works on root proc
          update_OK=(bulk>0.0_dp).and.(bulk<max_B)
       end if
       !       call comms_gcopy(update_OK,1)

       if (update_OK.and.update_params) then
          new_modulus_est=bulk
          if (pub_on_root.and..not.local_final) write (stdout,12) trim(geom_string)//': updated estimated bulk modulus =', &
               & new_modulus_est/29421.912_dp,"GPa"

       elseif(.not.update_params) then
          write (stdout,*) 'geom_update_inv_Hessian_params: bulk = ',bulk/29421.912_dp,trim("GPa")
       else
          ! lk: changed output detail from > NORMAL to > VERBOSE.
          if ((pub_geom_output_detail>VERBOSE).and.pub_on_root) then
             write (stdout,*) 'geom_update_inv_Hessian_param: skipping bulk modulus update because'
             if (bulk<0) then
                write (stdout,*) '                               update B <0'
             else
                write (stdout,*) '                               update B >3000 GPa'
             end if
          end if
       end if

          end if

    deallocate(mass_vector,stat=ierr)
    call utils_dealloc_check('geom_update_inv_hessian_params','mass_vector',ierr)

    if(allocated(Abar)) then
       deallocate (Abar,stat=ierr)
       call utils_dealloc_check('geom_update_inv_hessian_params','Abar',ierr)
       end if
    !parallel sync
    !    call comms_gcopy(new_frequency_est,1)
    !    call comms_gcopy(new_modulus_est,1)

12  format(1x,a,f10.5,1x,a)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_update_inv_hessian_params'

    return

  contains

    subroutine geom_calculate_update_basis(reduced_mem,ndim,cell_size,num_cell_basis_vecs,&
         &num_ion_basis_vecs)
      !=========================================================================!
      ! Calculate a basis from the inverse Hessian matrix, either using an SVD  !
      ! decomposition of the matrix itself in the case of BFGS, or through QR   !
      ! factorisation and SVD in the case of LBFGS.                             !
      !                                                                         !
      ! NB In the case of LBFGS, there is a contained routine which handles the !
      !    manipulations, however it cannot handle the case when ndim is        !
      !    smaller than cmdl%lbfgs_num_updates*2. For such a situation, a        !
      !    temporary matrix is created and a standard SVD is performed.         !
      !-------------------------------------------------------------------------!
      ! Arguments:                                                              !
      !   reduced_mem, intent=in, choose between "BFGS" and "LBFGS"             !
      !   ndim,        intent=in, number of dimensions of system                !
      !   cell_size,   intent=in, symmetrised (6) or unsymmetric (9) cell size  !
      !   num_cell_basis_vecs, intent=out, number of calculated cell DOFs       !
      !   num_ion_basis_vecs, intent=out, number of calculated ion DOFs         !
      !                                                                         !
      !-------------------------------------------------------------------------!
      ! Parent subroutine variables used:                                       !
      !   cmdl               the model, used only for lbfgs_num_updates          !
      !   basis_vecs        the basis vectors                                   !
      !   hessian_update    the hessian update matrix (BFGS)                    !
      !   tmpgradpos        the position and gradient update pairs (LBFGS)      !
      !   tmpmm             the middle, LBFGS recursion matrix.                 !
      !-------------------------------------------------------------------------!
      ! Modules used:                                                           !
      !                                                                         !
      !-------------------------------------------------------------------------!
      ! Key Internal Variables:                                                 !
      !-------------------------------------------------------------------------!
      ! Necessary conditions:                                                   !
      !-------------------------------------------------------------------------!
      ! Written by Jolyon Aarons, v1.0 26/05/2011                               !
      !=========================================================================!
      use comms, only: comms_bcast
      use constants, only : dp, stdout
      use simulation_cell, only: castep_model
      use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
      implicit none

      ! Arguments
      character(len=*),              intent(in)  :: reduced_mem
      integer,                       intent(in)  :: ndim
      integer,                       intent(in)  :: cell_size
      integer,                      intent(out)  :: num_cell_basis_vecs
      integer,                      intent(out)  :: num_ion_basis_vecs

      ! LAPACK subroutine
      external :: dgesvd

      ! Local variables
      real(kind=dp), dimension(:),   allocatable :: D_ion, D_cell
      real(kind=dp), dimension(:,:), allocatable :: Q_ion, Q_cell
      real(kind=dp), dimension(:,:), allocatable :: cell_update,ion_update
      real(kind=dp), dimension(:,:), allocatable :: hess_bkp
      real(kind=dp), dimension(:),   allocatable :: work
      integer                                    :: work_size
      integer                                    :: info
      integer                                    :: i
      real(kind=dp)                              :: eps=epsilon(1.0_dp)
      integer                                    :: num_basis_vecs
      logical                                    :: do_bfgs, do_lbfgs
      integer                                    :: ierr

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering geom_calculate_update_basis'

      !Find out whether we are doing BFGS or LBFGS
      call utils_assert((trim(reduced_mem).eq."BFGS").or.(trim(reduced_mem).eq."LBFGS"), &
           &'Error in geom_calculate_update_basis : reduced_mem argument must be either "BFGS" or "LBFGS"!')

      if(trim(reduced_mem).eq."BFGS") then
         do_bfgs=.true.
         do_lbfgs=.false.
      else if(trim(reduced_mem).eq."LBFGS") then
         do_lbfgs=.true.
         do_bfgs=.false.
      end if

      !Check we have everything we need
      num_cell_basis_vecs=0
      num_ion_basis_vecs=0


      !Check whether we're doing variable cell
      if(cell_size.gt.0) then
         allocate(Q_cell(1:cell_size,1:cell_size),stat=ierr)
         call utils_alloc_check('geom_calculate_update_basis','Q_cell',ierr)
         allocate(D_cell(1:cell_size),stat=ierr)
         call utils_alloc_check('geom_calculate_update_basis','D_cell',ierr)

         work_size=5*cell_size
         allocate(work(work_size),stat=ierr)
         call utils_alloc_check('geom_calculate_update_basis','work',ierr)
         !We will always do the cell part of the basis as a full SVD, because there are frequently more than 6 or 9
         !update vectors.
         if(do_bfgs) then
            allocate(hess_bkp(cell_size,cell_size),stat=ierr)
            call utils_alloc_check('geom_calculate_update_basis','hess_bkp',ierr)
            hess_bkp=cmdl%bfgs_inv_hessian(1:cell_size,1:cell_size)
            !Use LAPACK to perform the SVD
            call dgesvd('S','N',cell_size,cell_size,cmdl%bfgs_inv_hessian,size(cmdl%bfgs_inv_hessian,1),&
                 &D_cell, Q_cell, size(Q_cell,1), Q_cell, size(Q_cell,1), work, work_size, info)
            cmdl%bfgs_inv_hessian(1:cell_size,1:cell_size)=hess_bkp
            deallocate(hess_bkp,stat=ierr)
            call utils_dealloc_check('geom_calculate_update_basis','hess_bkp',ierr)
         else if(do_lbfgs) then
            !Create a temporary array, cell_update to hold the cell-hessian.
            allocate(cell_update(1:cell_size,1:cell_size),stat=ierr)
            call utils_alloc_check('geom_calculate_update_basis','cell_update',ierr)
            !form the cell-hessian from the lbfgs arrays.
            cell_update=matmul(tmpgradpos(1:cell_size,1:2*cmdl%lbfgs_num_updates),matmul(tmpmm,&
                 &transpose(tmpgradpos(1:cell_size,1:2*cmdl%lbfgs_num_updates))))
            ! SVD as before
            call dgesvd('S','N',cell_size,cell_size,cell_update,size(cell_update,1),&
                 &D_cell, Q_cell, size(Q_cell,1), Q_cell, size(Q_cell,1), work, work_size, info)
            deallocate(cell_update,stat=ierr)
            call utils_dealloc_check('geom_calculate_update_basis','cell_update',ierr)
         end if
         deallocate(work,stat=ierr)
         call utils_dealloc_check('geom_calculate_update_basis','work',ierr)
         !Count the singular values --> cell DOF.
         num_cell_basis_vecs=0
         do i=1,cell_size
            if(abs(D_cell(i))>eps) then
               num_cell_basis_vecs=num_cell_basis_vecs+1
            else
               exit
            end if
         end do
            deallocate(D_cell,stat=ierr)
            call utils_dealloc_check('geom_calculate_update_basis','D_cell',ierr)
      end if

      !Check we are doing a variable ion calculation
      if(ndim-cell_size.gt.0) then
         !for BFGS or LBFGS when max DOF in the hessian are fewer than the number of updates:
         !   normal SVD by the same process as cell part
         if((do_bfgs).or.(do_lbfgs.and.(2*cmdl%lbfgs_num_updates.ge.ndim-cell_size))) then
            allocate(Q_ion(1:ndim-cell_size,1:ndim-cell_size),stat=ierr)
            call utils_alloc_check('geom_calculate_update_basis','Q_ion',ierr)
            allocate(D_ion(1:ndim-cell_size),stat=ierr)
            call utils_alloc_check('geom_calculate_update_basis','D_ion',ierr)

            work_size=5*ndim-cell_size
            allocate(work(work_size),stat=ierr)
            call utils_alloc_check('geom_calculate_update_basis','work',ierr)

            if(do_bfgs) then
               allocate(hess_bkp(ndim-cell_size,ndim-cell_size),stat=ierr)
               call utils_alloc_check('geom_calculate_update_basis','hess_bkp',ierr)
               hess_bkp=cmdl%bfgs_inv_hessian(cell_size+1:,cell_size+1:)
               call dgesvd('S','N',ndim-cell_size,ndim-cell_size,cmdl%bfgs_inv_hessian(cell_size+1,cell_size+1),&
                    &size(cmdl%bfgs_inv_hessian,1),D_ion, Q_ion, size(Q_ion,1), Q_ion, size(Q_ion,1),&
                    &work, work_size, info)
               cmdl%bfgs_inv_hessian(cell_size+1:,cell_size+1:)=hess_bkp
               deallocate(hess_bkp,stat=ierr)
               call utils_dealloc_check('geom_calculate_update_basis','hess_bkp',ierr)

            else if(do_lbfgs) then
               allocate(ion_update(1:ndim-cell_size,1:ndim-cell_size),stat=ierr)
               call utils_alloc_check('geom_calculate_update_basis','ion_update',ierr)
               ion_update=matmul(tmpgradpos(1+cell_size:ndim,1:2*cmdl%lbfgs_num_updates),matmul(tmpmm,&
                    &transpose(tmpgradpos(1+cell_size:ndim,1:2*cmdl%lbfgs_num_updates))))
               call dgesvd('S','N',ndim-cell_size,ndim-cell_size,ion_update(1,1),&
                    &size(ion_update,1),D_ion, Q_ion, size(Q_ion,1), Q_ion, size(Q_ion,1),&
                    &work, work_size, info)
               deallocate(ion_update,stat=ierr)
               call utils_dealloc_check('geom_calculate_update_basis','ion_update',ierr)
            end if
            deallocate(work,stat=ierr)
            call utils_dealloc_check('geom_calculate_update_basis','work',ierr)

            num_ion_basis_vecs=0
            do i=1,ndim-cell_size
               if(abs(D_ion(i))>eps) then
                  num_ion_basis_vecs=num_ion_basis_vecs+1
               else
                  exit
               end if
            end do
            deallocate(D_ion,stat=ierr)
            call utils_dealloc_check('geom_calculate_update_basis','D_ion',ierr)

            ! For LBFGS in the general case:
         else if(do_lbfgs.and.(2*cmdl%lbfgs_num_updates.lt.ndim-cell_size)) then
            allocate(Q_ion(1:ndim-cell_size,1:cmdl%lbfgs_num_updates*2),stat=ierr)
            call utils_alloc_check('geom_calculate_update_basis','Q_ion',ierr)
            allocate(D_ion(1:cmdl%lbfgs_num_updates*2),stat=ierr)
            call utils_alloc_check('geom_calculate_update_basis','D_ion',ierr)
            ! Hand control over to the dedicated, contained routine.
            call lbfgs_calculate_basis(ndim-cell_size,cmdl,tmpgradpos(cell_size+1:ndim,:), tmpmm, &
                 & num_ion_basis_vecs, Q_ion, D_ion)
            do i=1,num_ion_basis_vecs
               if(abs(D_ion(i)).le.eps) then
                  exit
               end if
            end do
            deallocate(D_ion,stat=ierr)
            call utils_dealloc_check('geom_calculate_update_basis','D_ion',ierr)
            num_ion_basis_vecs=i-1;
         end if
      end if


      ! Total number of basis vectors, cell and ion...
      num_basis_vecs=num_cell_basis_vecs+num_ion_basis_vecs

      ! Deallocate Basis_vecs
      if(allocated(basis_vecs)) then
         deallocate(basis_vecs,stat=ierr)
         call utils_dealloc_check('geom_calculate_update_basis','basis_vecs',ierr)
      end if

      ! Allocate Basis with the correct sizes
      allocate(basis_vecs(1:ndim, 1:num_basis_vecs),stat=ierr)
      call utils_alloc_check('geom_calculate_update_basis','basis_vecs',ierr)
      basis_vecs=0.0_dp
      ! Fill it with our cell and ion basis vectors.
      if(cell_size.gt.0) then
         basis_vecs(1:cell_size,1:num_cell_basis_vecs)=Q_cell(1:cell_size,1:num_cell_basis_vecs)
         deallocate(Q_cell,stat=ierr)
         call utils_dealloc_check('geom_calculate_update_basis','Q_cell',ierr)
      end if
      if(ndim-cell_size.gt.0) then
         basis_vecs(cell_size+1:ndim,1+num_cell_basis_vecs:num_basis_vecs)=Q_ion(1:ndim-cell_size,1:num_ion_basis_vecs)
         deallocate(Q_ion,stat=ierr)
         call utils_dealloc_check('geom_calculate_update_basis','Q_ion',ierr)
      end if

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving geom_calculate_update_basis'

    end subroutine geom_calculate_update_basis

    subroutine lbfgs_calculate_basis(ndim,cmdl,update_matrix, lbfgs_W_matrix, data_size, basis_vecs, eigenvalues)
      !=========================================================================!
      ! This routine finds the eigenvalues and eigenvectors of an lbfgs matrix  !
      ! product.                                                                !
      !-------------------------------------------------------------------------!
      ! Arguments:                                                              !
      !   ndim, intent=in,      number of degrees of freedom in system          !
      !   cmdl, intent=in,       current state of the system                     !
      !   update_matrix, intent(in), matrix of position and gradient updates    !
      !      --> U=[P  G] where P=[p1 p2 p3...pndim] and G=[g1 g2 g3 ... gndim] !
      !          and p1,etc and g1,etc are position and gradient update vectors !
      !   lbfgs_W_matrix, intent(in), LBFGS middle-recursion matrix             !
      !                                (see Byrd-Nocedal (1984))                !
      !   data_size ,intent(out), calculated size of basis in number of vecs   !
      !   basis_vecs ,intent(out), the calculated basis vecs                    !
      !   eigenvalues,intent(out), and the eigenvalues returned in a vector     !
      !-------------------------------------------------------------------------!
      ! Parent module variables used:                                           !
      !   ndim and on_root                                                      !
      !-------------------------------------------------------------------------!
      ! Modules used:                                                           !
      !   algor : algor_sort                                                    !
      !-------------------------------------------------------------------------!
      ! Key Internal Variables:                                                 !
      !-------------------------------------------------------------------------!
      ! Necessary conditions:                                                   !
      ! Doing LBFGS                                                             !
      !-------------------------------------------------------------------------!
      ! Written by Jolyon Aarons, v1.0, 26/05/11                                !
      !=========================================================================!
      use comms, only: comms_bcast
      use simulation_cell, only: castep_model
      use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort, &
           utils_assert
      implicit none

      ! Arguments
      integer,                                    intent(in)    :: ndim
      type(castep_model),                          intent(in)    :: cmdl
      real(kind=dp), dimension(:,:),              intent(in)    :: update_matrix
      real(kind=dp), dimension(:,:),              intent(in)    :: lbfgs_W_matrix
      integer,                                    intent(out)   :: data_size
      real(kind=dp), dimension(:,:), allocatable, intent(out)   :: basis_vecs
      real(kind=dp), dimension(:),   allocatable, intent(out)   :: eigenvalues

      ! BLAS & LAPACK subroutines
      external :: dtrmm, dormqr, dsyevd, dgeqrf

      ! Local variables
      integer :: work_size, info, i, iwork_size, ierr
      real(kind=dp),dimension(:),allocatable :: work, qr_tau, inner_eigs
      real(kind=dp),dimension(:,:),allocatable :: qr_mat, rwr_mat, sorted_inner_vecs
      integer,dimension(:),allocatable :: iwork, pivot

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering lbfgs_calculate_basis'

      if(allocated(basis_vecs)) then
         deallocate(basis_vecs,stat=ierr)
         call utils_dealloc_check('lbfgs_calculate_basis','basis_vecs',ierr)
      end if
      if(allocated(eigenvalues)) then
         deallocate(eigenvalues,stat=ierr)
         call utils_dealloc_check('lbfgs_calculate_basis','eigenvalues',ierr)
      end if

      !QR factorise the update-matrix (does not explicitly compute Q)
      work_size=max(1,ndim)
      allocate(work(1:work_size),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','work',ierr)
      allocate(qr_mat(1:ndim,1:2*cmdl%lbfgs_num_updates),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','qr_mat',ierr)
      allocate(qr_tau(1:2*cmdl%lbfgs_num_updates),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','qr_tau',ierr)
      qr_mat=update_matrix
      call dgeqrf( ndim, cmdl%lbfgs_num_updates*2, qr_mat, size(qr_mat,1), qr_tau, work, work_size, info )
      call utils_assert(info==0,'Error in calculating QR-factorisation in lbfgs_calculate_basis')

      deallocate(work,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','work',ierr)

      !Fold R into the middle matrix --> RWR'
      allocate(rwr_mat(1:cmdl%lbfgs_num_updates*2,1:cmdl%lbfgs_num_updates*2),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','rwr_mat',ierr)
      rwr_mat=lbfgs_W_matrix
      call dtrmm('left','upper','notrans','nounitdiag',2*cmdl%lbfgs_num_updates,2*cmdl%lbfgs_num_updates,&
           &1.0_dp,qr_mat,size(qr_mat,1),rwr_mat,size(rwr_mat,1))
      call dtrmm('right','upper','trans','nounitdiag',2*cmdl%lbfgs_num_updates,2*cmdl%lbfgs_num_updates,&
           &1.0_dp,qr_mat,size(qr_mat,1),rwr_mat,size(rwr_mat,1))

      !Find the eigenvalues and eigenvectors of this new middle matrix
      work_size=-1
      iwork_size=-1
      allocate(work(1),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','work',ierr)
      allocate(iwork(1),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','iwork',ierr)
      allocate(inner_eigs(1:2*cmdl%lbfgs_num_updates),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','inner_eigs',ierr)
      allocate(pivot(1:2*cmdl%lbfgs_num_updates),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','pivot',ierr)
      call dsyevd('V','U',2*cmdl%lbfgs_num_updates,rwr_mat,size(rwr_mat,1),inner_eigs,work,work_size,iwork,iwork_size, info)
      call utils_assert(info==0, 'Error in performing specular decomposition in lbfgs_calculate_basis')

      work_size=INT(work(1))
      iwork_size=iwork(1)
      deallocate(work,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','work',ierr)
      deallocate(iwork,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','iwork',ierr)
      allocate(iwork(iwork_size),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','iwork',ierr)
      allocate(work(work_size),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','work',ierr)
      call dsyevd('V','U',2*cmdl%lbfgs_num_updates,rwr_mat,size(rwr_mat,1),inner_eigs,work,work_size,iwork,iwork_size, info)
      call utils_assert(info==0, 'Error in performing specular decomposition in lbfgs_calculate_basis')

      deallocate(work,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','work',ierr)
      deallocate(iwork,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','work',ierr)

      !Sort the eigenvalues by magnitude and produce a pivot-vector
      allocate(eigenvalues(1:2*cmdl%lbfgs_num_updates),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','eigenvalues',ierr)
      eigenvalues=abs(inner_eigs)
      !      call algor_sort(2*cmdl%lbfgs_num_updates,eigenvalues,pivot)
      call utils_heapsort(2*cmdl%lbfgs_num_updates,eigenvalues,pivot)
      eigenvalues=0.0_dp
      !work out the rank --> data_size
      do i=1,2*cmdl%lbfgs_num_updates
         !       if(abs(inner_eigs(pivot(i))).lt.epsilon(inner_eigs(1))) then
         !          data_size=i-1
         !          exit
         !       end if
         eigenvalues(i)=inner_eigs(pivot(i))
      end do
      data_size=2*cmdl%lbfgs_num_updates
      !
      !    !reduce size of eigenvalue vector to data_size
      !    inner_eigs=eigenvalues
      !    deallocate(eigenvalues,stat=ierr)
      !    if (ierr/=0) call io_abort('Error in deallocating eigenvalues in lbfgs_calculate_basis')
      !    allocate(eigenvalues(1:data_size),stat=ierr)
      !    if (ierr/=0) call io_allocate_abort('eigenvalues','lbfgs_calculate_basis')
      !    eigenvalues=inner_eigs(1:data_size)
      !    deallocate(inner_eigs,stat=ierr)
      !    if (ierr/=0) call io_abort('Error in deallocating inner_eigs in lbfgs_calculate_basis')
      !
      !reorder the eigenvectors by pivot vector
      allocate(sorted_inner_vecs(1:2*cmdl%lbfgs_num_updates,1:data_size),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','sorted_inner_vecs',ierr)
      sorted_inner_vecs=0.0_dp
      do i=1,data_size
         sorted_inner_vecs(:,i)=rwr_mat(:,pivot(i))
      end do
      deallocate(rwr_mat,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','rwr_mat',ierr)
      deallocate(pivot,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','pivot',ierr)

      !multiply the Q matrix by the eigenvectors to give the basis in the full space.
      allocate(basis_vecs(1:ndim,1:data_size),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','basis_vecs',ierr)
      basis_vecs=0.0_dp
      basis_vecs(1:2*cmdl%lbfgs_num_updates,1:data_size)=sorted_inner_vecs
      deallocate(sorted_inner_vecs,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','sorted_inner_vecs',ierr)

      work_size=-1
      allocate(work(1),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','work',ierr)
      call dormqr('l','n',ndim,data_size,2*cmdl%lbfgs_num_updates,qr_mat,size(qr_mat,1),qr_tau,basis_vecs,&
           &size(basis_vecs,1),work,work_size,info)
      call utils_assert(info==0, 'Error in multiplying by Q in lbfgs_calculate_basis')

      work_size=INT(work(1))
      deallocate(work,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','work',ierr)

      allocate(work(work_size),stat=ierr)
      call utils_alloc_check('lbfgs_calculate_basis','work',ierr)

      call dormqr('l','n',ndim,data_size,2*cmdl%lbfgs_num_updates,qr_mat,size(qr_mat,1),qr_tau,basis_vecs,&
           &size(basis_vecs,1),work,work_size,info)
      call utils_assert(info==0, 'Error in multiplying by Q in lbfgs_calculate_basis')

      deallocate(qr_mat,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','qr_mat',ierr)
      deallocate(qr_tau,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','qr_tau',ierr)
      deallocate(work,stat=ierr)
      call utils_dealloc_check('lbfgs_calculate_basis','work',ierr)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving lbfgs_calculate_basis'

    end subroutine lbfgs_calculate_basis

    subroutine geom_calc_bulk_freq(ndim,cell_size,num_cell_basis_vecs,num_ion_basis_vecs,Abar, &
         & real_lattice,basis_vecs,bulk,phonon)
      !=========================================================================!
      ! Update the parameters used to initialize the inverse Hessian matrix,    !
      ! i.e. geom_modulus_est and geom_frequency_est, by analysing the current  !
      ! inv_Hessian.                                                            !
      ! This follows Pfrommer section 3 closely and uses same notation.         !
      !                                                                         !
      ! NB Cannot call this routine WITHOUT then reseting inv_Hessian as it will!
      !    update geom_modulus_est and/or geom_frequency_est which on the next  !
      !    call will then be out of sync with what has actually been used!      !
      !-------------------------------------------------------------------------!
      ! Arguments:                                                              !
      !     ndim       , intent(in) :: number of degrees of freedom             !
      !     cell_size  , intent(in) :: whether to expect un-symmetrised cell    !
      !     --> which is 9x9 or Voigt notation which is 6x6.... given as 6 or 9 !
      !     num_cell_basis_vecs , intent(in) :: num of basis vecs in cell basis !
      !     num_ion_basis_vecs  , intent(in) :: num of basis vecs in ion basis  !
      !     Abar       , intent(in) :: reduced Hessian (as in Pfrommer)         !
      !     real_lattice , intent(in) ::   the real-space cell lattice vectors  !
      !     basis_vecs , intent(inout) ::  the Hessian basis... will be scaled  !
      !                                       on output                         !
      !     bulk  , intent(out) ::       calculated bulk modulus                !
      !     phonon  , intent(out) ::     calculated average phonon frequency at !
      !                                      gamma point                        !
      !-------------------------------------------------------------------------!
      ! Parent module variables used:                                           !
      !   ndim and on_root                                                      !
      !-------------------------------------------------------------------------!
      ! Modules used:                                                           !
      !   algor :: invert, constants :: dp                                      !
      !-------------------------------------------------------------------------!
      ! Key Internal Variables:                                                 !
      !-------------------------------------------------------------------------!
      ! Necessary conditions:                                                   !
      !-------------------------------------------------------------------------!
      ! Written by Jolyon Aarons, v1.0, 26/05/11                                !
      ! using some code from geom_update_inv_hessian_params, v2.0, written by   !
      ! Matt Probert, but                                                       !
      ! no longer performs Householder transformations, since the               !
      ! calculation of the basis is outsourced to geom_calculate_update_basis   !
      ! which performs this task differently.                                   !
      !=========================================================================!

      use comms, only: comms_bcast
      use constants, only : dp, stdout
      use linalg, only     : linalg_invert_sym_matrix
      use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

      implicit none

      ! Arguments
      integer,                                 intent(in)     :: ndim
      integer,                                 intent(in)     :: cell_size
      integer,                                 intent(in)     :: num_cell_basis_vecs
      integer,                                 intent(in)     :: num_ion_basis_vecs
      real(kind=dp), dimension(num_cell_basis_vecs+num_ion_basis_vecs,num_cell_basis_vecs+num_ion_basis_vecs), intent(in) :: Abar
      real(kind=dp), dimension(1:3,1:3),       intent(in)     :: real_lattice
      real(kind=dp), dimension(:,:),           intent(inout)  :: basis_vecs
      real(kind=dp),                           intent(out)    :: bulk
      real(kind=dp),                           intent(out)    :: phonon

      ! LAPACK subroutine
      external :: dpotrf

      ! Local variables
      integer                                              :: trace_max, trace_stride
      integer                                              :: num_basis_vecs
      integer                                              :: i,j,k
      integer                                              :: info

      real(kind=dp), dimension(:,:),           allocatable :: Bbar, Abar_cell_inv, Bbar_inv
      real(kind=dp), dimension(:,:),           allocatable :: S_matrix

      real(kind=dp)                                        :: sumtrace
      real(kind=dp)                                        :: vol

      real(kind=dp)                                        :: phonon_tol
      integer                                              :: phonon_count

      integer                                              :: ierr

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering geom_calc_bulk_freq'

      ! Perform det(real_lattice) to get the volume
      vol = real_lattice(1,1)*real_lattice(2,2)*real_lattice(3,3)  &
           - real_lattice(1,1)*real_lattice(2,3)*real_lattice(3,2)  &
           - real_lattice(1,2)*real_lattice(2,1)*real_lattice(3,3)  &
           + real_lattice(1,2)*real_lattice(2,3)*real_lattice(3,1)  &
           + real_lattice(1,3)*real_lattice(2,1)*real_lattice(3,2)  &
           - real_lattice(1,3)*real_lattice(2,2)*real_lattice(3,1)

      ! decide the stride lengths in the symmetrised, rank-2 strain tensor trace (cell_size==6)
      ! and the non-symmetrised, rank-4 strain tensor trace (cell_size==9).
      call utils_assert((cell_size==9).or.(cell_size==6),'Invalid cell_size in geom_calc_bulk_freq : must be either 6, or 9.')
      if(cell_size==9) then
         trace_max=9
         trace_stride=4
      else if(cell_size==6) then
         trace_max=3
         trace_stride=1
       end if

      !total basis size.
      num_basis_vecs=num_cell_basis_vecs+num_ion_basis_vecs

      if(num_cell_basis_vecs.gt.0) then
         if(num_ion_basis_vecs.gt.0) then
            ! We need an inverted ion-Abar for calculation of the bulk modulus.
            allocate(Abar_cell_inv(1:num_ion_basis_vecs,1:num_ion_basis_vecs),stat=ierr)
            call utils_alloc_check('geom_calc_bulk_freq','Abar_cell_inv',ierr)

            Abar_cell_inv=Abar(num_cell_basis_vecs+1:num_basis_vecs,num_cell_basis_vecs+1:num_basis_vecs)
            call linalg_invert_sym_matrix(Abar_cell_inv,num_ion_basis_vecs)
       end if

         allocate(Bbar(1:num_cell_basis_vecs,1:num_cell_basis_vecs),stat=ierr)
         call utils_alloc_check('geom_calc_bulk_freq','Bbar',ierr)

         !Construct the Elastic constants tensor from the symmetrised Hessian, Abar.
         Bbar=Abar(1:num_cell_basis_vecs,1:num_cell_basis_vecs)
         if(num_ion_basis_vecs.gt.0) then
            Bbar=Bbar-matmul(matmul(&
                 &Abar(1:num_cell_basis_vecs,num_cell_basis_vecs+1:num_basis_vecs),&
                 &Abar_cell_inv(1:num_ion_basis_vecs,1:num_ion_basis_vecs)),&
                 &Abar(num_cell_basis_vecs+1:num_basis_vecs,1:num_cell_basis_vecs))
            deallocate(Abar_cell_inv,stat=ierr)
            call utils_dealloc_check('geom_calc_bulk_freq','Abar_cell_inv',ierr)
       end if

         allocate(Bbar_inv(1:num_cell_basis_vecs,1:num_cell_basis_vecs),stat=ierr)
         call utils_alloc_check('geom_calc_bulk_freq','Bbar_inv',ierr)
         Bbar_inv=Bbar
         call linalg_invert_sym_matrix(Bbar_inv,num_cell_basis_vecs);
         deallocate(Bbar,stat=ierr)
         call utils_dealloc_check('geom_calc_bulk_freq','Bbar',ierr)

         ! Tensor trace to find the bulk modulus (note the special max and stride arrangements from before)!
         sumtrace=0.0_DP
         do i=1,num_cell_basis_vecs
            do j=1,num_cell_basis_vecs
               sumtrace=sumtrace + sum(basis_vecs(1:trace_max:trace_stride,i))*Bbar_inv(i,j)*&
                    &sum(basis_vecs(1:trace_max:trace_stride,j))
            end do
         end do
         deallocate(Bbar_inv,stat=ierr)
         call utils_dealloc_check('geom_calc_bulk_freq','Bbar_inv',ierr)

         bulk=1.0/(sumtrace*vol)

      end if

      if(num_ion_basis_vecs.gt.0) then
         ! cartesian transform
         do i=num_cell_basis_vecs+1,num_basis_vecs
            do k=cell_size+1,ndim,3
               basis_vecs(k:k+2,i)=matmul(real_lattice,basis_vecs(k:k+2,i))
            end do
         end do

         ! We have a generalised eigenvalue problem, with S_matrix the overlap!
         allocate(S_matrix(1:num_ion_basis_vecs,1:num_ion_basis_vecs),stat=ierr)
         call utils_alloc_check('geom_calc_bulk_freq','S_matrix',ierr)
         S_matrix=matmul(transpose(basis_vecs(cell_size+1:ndim,num_cell_basis_vecs+1:num_basis_vecs)),&
              &basis_vecs(cell_size+1:ndim,num_cell_basis_vecs+1:num_basis_vecs));
         call linalg_invert_sym_matrix(S_matrix,num_ion_basis_vecs);
         !But, we'll transform it into a standard eigenvalue problem.
         call dpotrf('u',num_ion_basis_vecs, S_matrix,size(S_matrix,1),info)
         S_matrix=matmul(transpose(S_matrix),matmul(Abar(num_cell_basis_vecs+1:num_basis_vecs,&
              &num_cell_basis_vecs+1:num_basis_vecs),(S_matrix)))
         !      S_matrix=matmul(Abar(num_cell_basis_vecs+1:num_basis_vecs,&
         !           &num_cell_basis_vecs+1:num_basis_vecs),S_matrix)
         ! Cholesky transform the matrix
         !         if(.not.debug_square_order) call dpotrf('u',num_ion_basis_vecs, S_matrix,size(S_matrix,1),info)
         ! and now we can trace it to get an average square-rooted eigenvalue
         phonon=0.0_dp
         phonon_tol=(396.0_dp*1.5198298455E-4_dp)**2
         phonon_count=0
         do i=1,num_ion_basis_vecs
            if((S_matrix(i,i).gt.0.0_dp).and.(S_matrix(i,i).lt.phonon_tol)) then
               phonon=phonon+sqrt(S_matrix(i,i))
               phonon_count=phonon_count+1
    end if
         end do
         !Which is the average phonon-mode at the Gamma-point.
         !         if(debug_square_order) phonon=sqrt(phonon)
         !         phonon=phonon/real(num_ion_basis_vecs,dp)
         if (phonon_count>0) phonon=phonon/real(phonon_count,dp)

         !  omega_squared=eig(Abar(num_cell_basis_vecs+1:num_basis_vecs,num_cell_basis_vecs+1:num_basis_vecs),S_matrix)

         !  deallocate(Abar)
         deallocate(S_matrix,stat=ierr)
         call utils_dealloc_check('geom_calc_bulk_freq','S_matrix',ierr)
      end if

    if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving geom_calc_bulk_freq'

    end subroutine geom_calc_bulk_freq

    !   subroutine geom_calculate_abar(reduced_mem,ndim, cmdl, num_basis_vecs, hessian, &
    !        &lbfgs_init_hess, tmpgradpos, tmpmm, mass, basis_vecs, Abar)
    subroutine geom_calculate_abar(reduced_mem,ndim, num_basis_vecs)
      !=========================================================================!
      ! Calculate a Hessian matrix in the reduced basis.                        !
      !-------------------------------------------------------------------------!
      ! Arguments:                                                              !
      !   reduced_mem, intent=in, choose between "BFGS" and "LBFGS"             !
      !   ndim,        intent=in, number of dimensions of system                !
      !   cell_size,   intent=in, symmetrised (6) or unsymmetric (9) cell size  !
      !   cmdl,         intent=in, the model, used only for lbfgs_num_updates    !
      ! hessian_update intent=in, the hessian update matrix (BFGS)              !
      !   tmpgradpos,  intent=in, the position and gradient update pairs (LBFGS)!
      !   tmpmm,       intent=in, the middle, LBFGS recursion matrix.           !
      !   num_cell_basis_vecs, intent=out, number of calculated cell DOFs       !
      !   num_ion_basis_vecs, intent=out, number of calculated ion DOFs         !
      !                                                                         !
      !-------------------------------------------------------------------------!
      ! Parent module variables used:                                           !
      !   basis                                                                 !
      !-------------------------------------------------------------------------!
      ! Modules used:                                                           !
      !                                                                         !
      !-------------------------------------------------------------------------!
      ! Key Internal Variables:                                                 !
      !-------------------------------------------------------------------------!
      ! Necessary conditions:                                                   !
      !-------------------------------------------------------------------------!
      ! Written by Jolyon Aarons, v1.0 26/05/2011                               !
      !=========================================================================!

      use comms, only: comms_bcast
      use constants, only    : dp, stdout
      use linalg, only     : linalg_invert_sym_matrix
      use simulation_cell, only: castep_model
      use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

      implicit none

      ! Arguments
      character(len=*),                intent(in)  :: reduced_mem
      integer,                       intent(in)    :: ndim
!      type(castep_model),             intent(in)    :: cmdl
      integer,                       intent(in)    :: num_basis_vecs
      !      real(kind=dp), dimension(ndim,ndim),                                   optional,intent(inout) :: hessian
      !      real(kind=dp), dimension(1:3, 1:ndim),                                 optional,intent(in)    :: lbfgs_init_hess
      !      real(kind=dp), dimension(1:ndim, 1:2*cmdl%lbfgs_num_updates),           optional,intent(inout) :: tmpgradpos
      !      real(kind=dp), dimension(1:2*cmdl%lbfgs_num_updates, 1:2*cmdl%lbfgs_num_updates), optional,intent(in) :: tmpmm
      !      real(kind=dp), dimension(1:ndim),                                               intent(in)    :: mass
      !      real(kind=dp), dimension(1:ndim, 1:num_basis_vecs),                             intent(in)    :: basis_vecs
      !      real(kind=dp), dimension(1:num_basis_vecs,1:num_basis_vecs),                    intent(out)   :: Abar

      ! BLAS subroutine
      external :: dgemm

      ! Local variables
      real(kind=dp), dimension(:,:), allocatable   :: mass_scaled_H0
      real(kind=dp), dimension(:,:), allocatable   :: update_space_H0_p1
      real(kind=dp), dimension(:,:), allocatable   :: update_space_H0
      real(kind=dp), dimension(:,:), allocatable   :: update_space_H_up
      integer                                      :: i
      integer                                      :: info
      integer                                      :: ierr
      logical                                      :: do_bfgs, do_lbfgs

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering geom_calculate_abar'

      call utils_assert((trim(reduced_mem).eq."BFGS").or.(trim(reduced_mem).eq."LBFGS"),&
           & 'geom_calculate_abar : reduced_mem argument must be either "BFGS" or "LBFGS"!')
      if(trim(reduced_mem).eq."BFGS") then
         do_bfgs=.true.
         do_lbfgs=.false.
      else if(trim(reduced_mem).eq."LBFGS") then
         do_lbfgs=.true.
         do_bfgs=.false.
      end if

      !      if(do_bfgs) then
      !         if(.not.present(hessian)) then
      !            call io_abort('Hessian is required to calculate Abar in BFGS : geom_calculate_abar')
      !         end if
      !      else if(do_lbfgs) then
      !         if((.not.present(tmpgradpos)).or.(.not.present(tmpmm)).or.(.not.present(lbfgs_init_hess))) then
      !            call io_abort('Update matrix, recursion matrix and initial hessian are required to &
      !                 &calculate Abar in LBFGS : geom_calculate_abar')
      !         end if
      !      end if

      Abar=0.0_dp

      if(do_bfgs) then
         do i=1,ndim
            cmdl%bfgs_inv_hessian(i,:)=cmdl%bfgs_inv_hessian(i,:)*mass_vector(i)
         end do
         do i=1,ndim
            cmdl%bfgs_inv_hessian(:,i)=cmdl%bfgs_inv_hessian(:,i)*mass_vector(i)
         end do
         Abar=matmul(matmul(transpose(basis_vecs),cmdl%bfgs_inv_hessian),basis_vecs)
      else if(do_lbfgs) then

         allocate(mass_scaled_H0(1:3,1:ndim),stat=ierr)
         call utils_alloc_check('geom_calculate_abar','mass_scaled_H0',ierr)
         mass_scaled_H0=0.0_dp
         !call outputdata(transpose(lbfgs_init_inv_hess))
         mass_scaled_H0(3,1)=lbfgs_init_inv_hess(3,1)*mass_vector(1)
         mass_scaled_H0(2,2)=lbfgs_init_inv_hess(2,2)*mass_vector(2)
         mass_scaled_H0(3,2)=lbfgs_init_inv_hess(3,2)*mass_vector(3)
         do i=3,ndim
            mass_scaled_H0(1,i)=lbfgs_init_inv_hess(1,i)*mass_vector(i-2)
            mass_scaled_H0(2,i)=lbfgs_init_inv_hess(2,i)*mass_vector(i-1)
            mass_scaled_H0(3,i)=lbfgs_init_inv_hess(3,i)*mass_vector(i)
         end do
         do i=1,ndim
            mass_scaled_H0(:,i)=mass_scaled_H0(:,i)*mass_vector(i)
         end do

         ! Y' H_0 Y
         allocate(update_space_H0_p1(1:ndim,num_basis_vecs),stat=ierr)
         call utils_alloc_check('geom_calculate_abar','update_space_H0_p1',ierr)
         call lbfgs_initH_mul('n',basis_vecs,mass_scaled_H0,update_space_H0_p1,info)
         deallocate(mass_scaled_H0,stat=ierr)
         call utils_dealloc_check('geom_calculate_abar','update_space_H0',ierr)
         allocate(update_space_H0(1:num_basis_vecs,1:num_basis_vecs),stat=ierr)
         call utils_alloc_check('geom_calculate_abar','update_space_H0',ierr)

         call dgemm('t','n',num_basis_vecs,num_basis_vecs,ndim,1.0_dp,basis_vecs,size(basis_vecs,1),&
              &update_space_H0_p1,size(update_space_H0_p1,1),0.0_dp,update_space_H0,&
              &size(update_space_H0,1))
         deallocate(update_space_H0_p1,stat=ierr)
         call utils_dealloc_check('geom_calculate_abar','update_space_H0_p1',ierr)

         do i=1,ndim
            tmpgradpos(i,:)=mass_vector(i)*tmpgradpos(i,:)
         end do
         allocate(update_space_H_up(1:num_basis_vecs,1:2*cmdl%lbfgs_num_updates),stat=ierr)
         call utils_alloc_check('geom_calculate_abar','update_space_H_up',ierr)
         update_space_H_up=matmul(transpose(basis_vecs),tmpgradpos)
         Abar=update_space_H0

         Abar=Abar+ matmul(update_space_H_up,matmul(tmpmm,transpose(update_space_H_up)))


      end if

      call linalg_invert_sym_matrix(Abar,num_basis_vecs)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving geom_calculate_abar'

    end subroutine geom_calculate_abar

  end subroutine geom_update_inv_hessian_params

  subroutine geom_get_forces(cmdl,mdl,path,is_cmplx,en_sync)
    !=========================================================================!
    ! Find ground state energy, forces and stresses for given model.          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=inout, the model to be evaluated                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   fix_all_ions, constrain_ions, symm_ions, fix_all_cell, constrain_cell,!
    !   and symm_cell are used to determine which firstd calls we need.       !
    !   pub_geom_output_detail, stdout and on_root for diagnostics                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model for type definitions and checkpointing                          !
    !   firstd for forces and stresses                                        !
    !   electronic for electron minimiser                                     !
    !   wave for adding random noise to highest occupied band                 !
    !   parameters to define the job etc                                      !
    !   cell for cell vectors and constraints                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   cmdl has already been properly read and initialised                   !
    !   current_cell must be valid                                            !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified to deal with embedding structures by Joseph Prentice, May 2018 !
    !=========================================================================!

    use constants, only: stdout
    use energy_and_force, only : energy_and_force_calculate
    use forces, only: forces_apply_constraints
    use model_type, only: MODEL
    use simulation_cell, only: castep_model, castep_cell_2_elements
    use neb, only: neb_path, neb_energy_and_force
    use rundat, only: pub_ngwf_regions_ngroups  ! jcap

    implicit none

    type(castep_model), intent(inout) :: cmdl
    type(MODEL),      intent(inout) :: mdl
    type(NEB_PATH),optional,intent(inout) :: path
    ! agrecocmplx
    logical,          intent(in) :: is_cmplx
    logical, optional, intent(in) :: en_sync

    !local variables
    integer       :: ii,is
    ! jcap: further local variables
    integer :: iregion,jat(pub_ngwf_regions_ngroups)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_get_forces'

    !-------------------------------------------------------------------------!
    !Now we must find the ground state wvfn for this model                    !
    !-------------------------------------------------------------------------!

    if (.not.cmdl%found_ground_state) then

       ! aam: Update ionic coordinates in elements with new ones from current_cell
       call castep_cell_2_elements(cmdl%cell,mdl%elements,mdl%nat)

       ! jcap: copy these ionic co-ordinates into the appropriate subregions
       jat = 0
       do ii=1,mdl%nat
          ! Get the region counter
          iregion = mdl%elements(ii)%region
          ! Count the no. of atoms allocated to each region
          jat(iregion) = jat(iregion) + 1
          ! Copy the co-ordinates over
          mdl%regions(iregion)%elements(jat(iregion))%centre = mdl%elements(ii)%centre
       end do

       ! aam: Calculate new ground state energy and forces
       ! agrecocmplx: use complex NGWFs if requested
       if (present(path)) then
         call neb_energy_and_force(cmdl%total_energy,cmdl%forces(:,1,:),mdl,path, &
              is_cmplx=is_cmplx,en_sync=en_sync)
       else
         call energy_and_force_calculate(cmdl%total_energy,cmdl%forces,mdl, &
              is_cmplx=is_cmplx)
       end if

       ! aam: In the presence of ionic constraints keep copy of unconstrained forces
       if (constrain_ions) then
          do ii=1,mdl%nat                          ! Each ion has its own species
             do is=1,current_cell%max_ions_in_species   ! ie max_ions_species = 1
                cmdl%orig_forces(1,is,ii) = cmdl%forces(1,is,ii)
                cmdl%orig_forces(2,is,ii) = cmdl%forces(2,is,ii)
                cmdl%orig_forces(3,is,ii) = cmdl%forces(3,is,ii)
             enddo
          enddo
       endif

       ! aam: Apply ionic constraints
       if (constrain_ions) call forces_apply_constraints(cmdl%forces,mdl)

       cmdl%found_ground_state = .true.
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_get_forces'

    return
  end subroutine geom_get_forces

  subroutine geom_converged(cmdl,old_x_vec,enthalpy,shuffle,dE,Fmax,dRmax,Smax, &
       & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax,mdl)
    !=========================================================================!
    ! Check convergence of current configuration of the system.               !
    ! Need to satisfy convergence in enthalpy for a given number of iterations!
    ! and SIMULTANEOUSLY get force, displacement and stress < tolerances.     !
    ! Complicated by possibility of user changing length of history during    !
    ! run and possibility of backtracking to old configurations.              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl,    intent=in, the current state of the system                    !
    !   old_x_vec,    =in, the previous fractional coords for calculating dR  !
    !   enthalpy,     =in, the current enthalpy                               !
    !                                                                         !
    !   shuffle, intent=in,integer flag to update, or not change histories    !
    !   dE,        intent=out, width of enthalpy window                       !
    !   Fmax,          "     , max force on any ion                           !
    !   dRmax,         "     , max displacement of any ion                    !
    !   Smax,          "     , max component of cell stress                   !
    !   converged,     "     , logical flag to signal model has converged     !
    !   converged_dE,  "     , logical flag to signal dE     "     "          !
    !   converged_Fmax,"     , logical flag to signal max force    "          !
    !   converged_dRmax,"    , logical flag to signal max disp     "          !
    !   converged_Smax, "    , logical flag to signal stress       "          !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_geom_output_detail, stdout and on_root for diagnostics                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !    parameters                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !  current_cell must be in sync with current model                        !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified for chain-of-states convergence by Kevin Duff, Jan 2018        !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    ! kaw: Gives access to the number of classical atoms.
    use comms, only: pub_on_root
    use constants, only: dp, stdout, NORMAL, VERBOSE
    use model_type, only: MODEL
    use rundat, only: pub_geom_convergence_win, pub_geom_disp_tol, &
         pub_geom_energy_tol, pub_geom_force_tol, pub_geom_output_detail, &
         pub_geom_continuation, pub_task, tssearch_method, tssearch_disp_tol, &
         tssearch_force_tol, tssearch_energy_tol, pub_task
    use simulation_cell, only: castep_model, castep_cell_frac_to_cart
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    type(MODEL),                      intent(in)  :: mdl
    type(castep_model),               intent(in)  :: cmdl
    real(kind=DP), dimension(1:ndim), intent(in)  :: old_x_vec
    real(kind=DP),                    intent(in)  :: enthalpy
    integer,                          intent(in)  :: shuffle

    real (kind=dp),                    intent(out) :: dE
    real (kind=dp),                    intent(inout) :: Fmax   !max Force
    real (kind=dp),                    intent(out) :: dRmax  !max displacement
    real (kind=dp),                    intent(out) :: Smax   !max Stress
    logical,                          intent(out) :: converged
    logical,                          intent(out) :: converged_dE
    logical,                          intent(out) :: converged_Fmax
    logical,                          intent(out) :: converged_dRmax
    logical,                          intent(out) :: converged_Smax

    !local variables
    integer, save                                        :: history_iter=0
    integer, parameter                                   :: max_history_length=10

    real(kind=DP), save, dimension(1:max_history_length) :: history_E    !NB E not dE
    real(kind=DP), save, dimension(1:max_history_length) :: history_Fmax
    real(kind=DP), save, dimension(1:max_history_length) :: history_dRmax
    real(kind=DP), save, dimension(1:max_history_length) :: history_Smax

    real(kind=DP), allocatable, dimension(:,:,:)         :: old_frac
    real(kind=DP), dimension(1:3)                        :: dR, dR_frac

    real(kind=DP) :: force_tol, energy_tol, disp_tol

    integer                                              :: i,j,iatom,ispec
    integer                                              :: num_shuffle,ierr

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_converged'

    ! ars: catch errors in tag 'shuffle'
    select case (shuffle)
    case (shuffle_forwards, shuffle_none, shuffle_backwards)
    case default
       call utils_abort('Error in geom_converged: &
            &unrecognised shuffle in geom_converged:',shuffle)
    end select

    ! kkbd: If we're relaxing a NEB chain-of-states, use tssearch tolerances
    if (pub_task == "TRANSITIONSTATESEARCH") then
       disp_tol   = tssearch_disp_tol
       force_tol  = tssearch_force_tol
       energy_tol = tssearch_energy_tol
    else
       disp_tol   = pub_geom_disp_tol
       force_tol  = pub_geom_force_tol
       energy_tol = pub_geom_energy_tol
    end if

    !cross-check history lengths
    !why ever would you want such a long convergence window? You must be joking ....
    if (pub_geom_convergence_win>max_history_length) then
       if (pub_on_root) write (stdout,*) 'geom_converged: truncating geom_convergence_win to ',max_history_length-1
       !NB we need the -1 to allow for backtracking
       pub_geom_convergence_win=max_history_length-1
    end if
    !... and a window of 1 will give all sorts of problems
    if (pub_geom_convergence_win<2) then
       if (pub_on_root) write (stdout,*) 'geom_converged: increasing geom_convergence_win to 2'
       pub_geom_convergence_win=2
    end if

    !initialize various bits
    allocate (old_frac(1:3,1:cmdl%cell%max_ions_in_species,1:cmdl%cell%num_species),stat=ierr)
    call utils_alloc_check('geom_converged','old_frac',ierr)

    i=10
    do ispec=1,cmdl%cell%num_species - mdl%nat_classical
       do iatom=1,cmdl%cell%num_ions_in_species(ispec)
          do j=1,3
             old_frac(j,iatom,ispec)=old_x_vec(i+j-1)
          end do
          i=i+3
       end do
    end do

    if (history_iter==0) then
       history_E    =0.0_dp
       history_Fmax =0.0_dp
       history_dRmax=0.0_dp
       history_Smax =0.0_dp
    end if

    !calculate convergence criteria (dE has to wait until done any shuffling)

    !calculate ionic max |F| and |dR| values
    Fmax =0.0_dp
    dRmax=0.0_dp
    if (.not.fix_all_ions) then
    ! kaw: Allows for classical atoms.
       do ispec=1,cmdl%cell%num_species - mdl%nat_classical
          do iatom=1,cmdl%cell%num_ions_in_species(ispec)
             Fmax=max(Fmax,dot_product(cmdl%forces(:,iatom,ispec),cmdl%forces(:,iatom,ispec)))
             dR_frac=geom_opt_min_image_frac(cmdl%cell%ionic_positions(:,iatom,ispec),old_frac(:,iatom,ispec))
             call castep_cell_frac_to_cart(current_cell,dR_frac,dR)
             dRmax=max(dRmax,dot_product(dR,dR))
          end do
       end do
       Fmax =sqrt(Fmax)
       dRmax=sqrt(dRmax)
    end if

    !calculate cell max |S| values
    Smax=0.0_dp

    !update history list? forwards/backwards/none?
    select case (shuffle)
    case (shuffle_forwards)       !first attempt at new configuation this iter

       !update history_iter
       history_iter=history_iter+1

       !now shuffle histories forwards, i.e. normal order in time
       num_shuffle=max(min(history_iter,max_history_length)-1,0)
       do i=num_shuffle,1,-1
          history_E(i+1)    =history_E(i)
          history_Fmax(i+1) =history_Fmax(i)
          history_dRmax(i+1)=history_dRmax(i)
          history_Smax(i+1) =history_Smax(i)
       end do

       !assign new values
       history_E(1)    =enthalpy/cmdl%cell%num_ions  !enthalpy    per atom
       history_Fmax(1) =Fmax                        !max |F|  on any atom
       history_dRmax(1)=dRmax                       !max |dR| of any atom
       history_Smax(1) =Smax                        !max S component

    case (shuffle_none)           !subsequent attempts at new configuration this iter
       !assign new values, overwriting previous attempt at this iter, without updating history_iter
       history_E(1)    =enthalpy/cmdl%cell%num_ions  !enthalpy    per atom
       history_Fmax(1) =Fmax                        !max |F|  on any atom
       history_dRmax(1)=dRmax                       !max |dR| of any atom
       history_Smax(1) =Smax                        !max S component

    case (shuffle_backwards)      !backtracking from current configuration to old one
       !now shuffle histories backwards as we are backtracking
       !and assign the old values to the working versions!
       !backtrack history_iter
       history_iter=history_iter-1

       num_shuffle=min(history_iter,max_history_length)-1
       do i=1,num_shuffle
          history_E(i)    =history_E(i+1)
          history_Fmax(i) =history_Fmax(i+1)
          history_dRmax(i)=history_dRmax(i+1)
          history_Smax(i) =history_Smax(i+1)
       end do

       !assign old values
       Fmax =history_Fmax(1)                        !max |F|  on any atom
       dRmax=history_dRmax(1)                       !max |dR| of any atom
       Smax =history_Smax(1)                        !max S component

    case default
       call utils_abort('Error in geom_converged: &
            &unrecognised shuffle in geom_converged: ',shuffle)
    end select

    !Calculate energy convergence
    !NB Energy is |max-min| enthalpy/atom over window, unlike others, as we don't know the min value!
    if (history_iter>=pub_geom_convergence_win) then
       dE=maxval(history_E(1:pub_geom_convergence_win))-minval(history_E(1:pub_geom_convergence_win))
    else
       dE=maxval(history_E(1:history_iter))-minval(history_E(1:history_iter))
    end if

    !test convergence
    if (history_iter>=pub_geom_convergence_win) then
       converged_dE=(dE<=energy_tol)
    else
       converged_dE=.false.
    end if
    converged_Fmax =( Fmax<=force_tol)
    if (history_iter>1) then
       converged_dRmax=(dRmax<=disp_tol)
    else
       converged_dRmax=.false.
    end if

    converged_Smax =( Smax<=geom_stress_tol)

    !spoof flags where irrevelant
    if (fix_all_ions) then
       converged_Fmax =.true.
       converged_dRmax=.true.
    end if
    converged_Smax =.true.

    !special case: initial configuration satisfies F and S so don't want to do
    !any more just to satisfy dE or dR
    ! qoh: However if we are doing a continuation don't stop just because the
    ! qoh: forces are converged.
    if (iteration==0 .and. converged_Fmax .and. converged_Smax &
         .and. .not. pub_geom_continuation) then
       converged_dE=.true.
       converged_dRmax=.true.
    end if

    !combine all the different parts into an overall convergence flag
    converged=(converged_dE.and.converged_Fmax.and.converged_dRmax.and.converged_Smax)

    !inform user
    ! lk: changed output detail from > NORMAL to > VERBOSE.
    if (pub_geom_output_detail>VERBOSE.and.pub_on_root) then
       write (stdout,99)    'geom_converged: history_iter=',history_iter
       do i=1,min(history_iter,max_history_length)
          write (stdout,98) 'geom_converged: history_E    (',i,')=',history_E(i),trim(energy_label)
          write (stdout,98) 'geom_converged: history_Fmax (',i,')=',history_Fmax(i),trim(force_label)
          write (stdout,98) 'geom_converged: history_dRmax(',i,')=',history_dRmax(i),trim(length_label)
          write (stdout,98) 'geom_converged: history_Smax (',i,')=',history_Smax(i),trim(pressure_label)
       end do
    end if

    !clean up at end
    deallocate (old_frac,stat=ierr)
    call utils_dealloc_check('geom_converged','old_frac',ierr)

98  format (1x,a,i4,a,1x,f20.6,1x,a)
99  format (1x,a,i4)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_converged'

    return
  end subroutine geom_converged

  subroutine geom_output_converged(bfgs_iteration,enthalpy,dE,Fmax,dRmax, &
        converged_dE,converged_Fmax,converged_dRmax)
    !=========================================================================!
    ! Output relevant convergence data to stdout                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   on_root                                                               !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified for chain-of-states convergence by Kevin Duff, Jan 2018        !
    ! Modified output detail and content by Loukas Kollias, Aug 2022          !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout
    use rundat, only: pub_geom_disp_tol, pub_geom_energy_tol, &
                      pub_geom_force_tol, tssearch_disp_tol, &
                      tssearch_force_tol, tssearch_energy_tol, pub_task
    use services, only: services_flush

    implicit none

    integer      , intent(in) :: bfgs_iteration
    real(kind=DP), intent(in) :: enthalpy

    real(kind=DP), intent(in) :: dE
    real(kind=DP), intent(in) :: Fmax   !max Force
    real(kind=DP), intent(in) :: dRmax  !max displacement

    logical,        intent(in) :: converged_dE
    logical,        intent(in) :: converged_Fmax
    logical,        intent(in) :: converged_dRmax

    !local vars
    character(len=80) :: data_string
    character(len=80) :: divider_string
    character(len=80) :: label_string
    integer           :: string_index
    integer           :: len_label, len_value, len_tol, len_unit, len_flag

    real(kind=DP) :: disp_tol, force_tol, energy_tol

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_output_converged'

    ! kkbd: If we're relaxing a NEB chain-of-states, use tssearch tolerances
    if (pub_task == "TRANSITIONSTATESEARCH") then
       disp_tol   = tssearch_disp_tol
       force_tol  = tssearch_force_tol
       energy_tol = tssearch_energy_tol
    else
       disp_tol   = pub_geom_disp_tol
       force_tol  = pub_geom_force_tol
       energy_tol = pub_geom_energy_tol
    end if

    !initialise strings
    data_string     = repeat(' ',len(data_string))
    divider_string  = repeat(' ',len(divider_string))
    label_string    = repeat(' ',len(label_string))
    !                  1234567890123 23456789012345678 23456789012345678 234567890123 23456
    divider_string  = '+-----------+-----------------+-----------------+------------+-----+'
    label_string    = '| Parameter |      value      |    tolerance    |    units   | OK? |'

    !write out all the relevant data (root proc only)
    if (pub_on_root) then

       !explain what we are testing ...
       write (stdout,1) trim(geom_string)//': finished iteration',bfgs_iteration,' with enthalpy=', &
            & enthalpy,trim(energy_label)

       !header
       write(stdout,*) ' '
       select case (geom_string)
       case ('BFGS')
          write(stdout,7) divider_string
          write(stdout,7) label_string
          write(stdout,7) divider_string
       case ('LBFGS')
          write(stdout,8) divider_string
          write(stdout,8) label_string
          write(stdout,8) divider_string
       case ('TPSD')
          write(stdout,9) divider_string
          write(stdout,9) label_string
          write(stdout,9) divider_string
       end select
       !dE                                                       '1234567890'
       string_index=1
       len_label   =13
       len_value   =18
       len_tol     =18
       len_unit    =13
       len_flag    = 6
       write(data_string(string_index:string_index+len_label),2) '  dE/ion  '
       string_index=string_index+len_label
       write(data_string(string_index:string_index+len_value),3) dE
       string_index=string_index+len_value
       write(data_string(string_index:string_index+len_tol),  4) energy_tol
       string_index=string_index+len_tol
       write(data_string(string_index:string_index+len_unit), 5) trim(energy_label)
       string_index=string_index+len_unit
       if (converged_dE) then
          write(data_string(string_index:string_index+len_flag), 6) 'Yes'
       else
          write(data_string(string_index:string_index+len_flag), 6) 'No '
       end if
       string_index=string_index+len_flag
       !cross check all is well (NB strings start from 1 not 0 ...)
       if (string_index>71) write (stdout,*) 'geom_output_converged: format problem?'
       if(trim(geom_string).eq.'BFGS') write(stdout,7) data_string
       if(trim(geom_string).eq.'LBFGS') write(stdout,8) data_string
       if(trim(geom_string).eq.'TPSD') write(stdout,9) data_string


       !|F|max                                                      '1234567890'
       if (.not.fix_all_ions) then
          string_index=1
          write(data_string(string_index:string_index+len_label),2) '  |F|max  '
          string_index=string_index+len_label
          write(data_string(string_index:string_index+len_value),3) Fmax
          string_index=string_index+len_value
          write(data_string(string_index:string_index+len_tol),  4) force_tol
          string_index=string_index+len_tol
          write(data_string(string_index:string_index+len_unit), 5) trim(force_label)
          string_index=string_index+len_unit
          if (converged_Fmax) then
             write(data_string(string_index:string_index+len_flag), 6) 'Yes'
          else
             write(data_string(string_index:string_index+len_flag), 6) 'No '
          end if
          if(trim(geom_string).eq.'BFGS') write(stdout,7) data_string
          if(trim(geom_string).eq.'LBFGS') write(stdout,8) data_string
          if(trim(geom_string).eq.'TPSD') write(stdout,9) data_string
       end if

       !|dR|max                                                     '1234567890'
       if (.not.fix_all_ions) then
          string_index=1
          write(data_string(string_index:string_index+len_label),2) '  |dR|max '
          string_index=string_index+len_label
          write(data_string(string_index:string_index+len_value),3) dRmax
          string_index=string_index+len_value
          write(data_string(string_index:string_index+len_tol),  4) disp_tol
          string_index=string_index+len_tol
          write(data_string(string_index:string_index+len_unit), 5) trim(length_label)
          string_index=string_index+len_unit
          if (converged_dRmax) then
             write(data_string(string_index:string_index+len_flag), 6) 'Yes'
          else
             write(data_string(string_index:string_index+len_flag), 6) 'No '
          end if
          if(trim(geom_string).eq.'BFGS') write(stdout,7) data_string
          if(trim(geom_string).eq.'LBFGS') write(stdout,8) data_string
          if(trim(geom_string).eq.'TPSD') write(stdout,9) data_string
       end if

       !footer
       if(trim(geom_string).eq.'BFGS') write(stdout,7) divider_string
       if(trim(geom_string).eq.'LBFGS') write(stdout,8) divider_string
       if(trim(geom_string).eq.'TPSD') write(stdout,9) divider_string
       write(stdout,*) ' '

       !back to all procs
    end if

    call services_flush

1   format(1x,a,i4,1x,a,es21.12e3,1x,a)

2   format('|',a10,    1x,'|')       !1+10+1+1=13
3   format(1x,es15.6e3,1x,'|')       !1+15+1+1=18
4   format(1x,es15.6e3,1x,'|')       !1+15+1+1=18
    !NB In IO the phys_unit is defined as 30 characters long but max. is 10 in practice
5   format(1x,a10,     1x,'|')       !1+10+1+1=13
6   format(1x,a3,      1x,'|')       !1+3+1+1=6
    !13+18+18+13+6=68
7   format(1x,a68,' <-- BFGS')
8   format(1x,a68,' <-- LBFGS')
9   format(1x,a68,' <-- TPSD')

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_output_converged'

    return
  end subroutine geom_output_converged

  subroutine geom_output(cmdl,mdl)
    !=========================================================================!
    ! Output relevant parts of model to stdout. Based on cell_output.         !
    ! NB unlike md_output this is does NOT output energies etc as called at   !
    !    intermediate stages of calculation as well as at end.                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=in, the model to be written                               !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   on_root                                                               !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

     use comms, only: pub_on_root
     use constants, only: dp, stdout, VERBOSE
     use model_type, only: MODEL
     use simulation_cell, only: castep_model, castep_cell_frac_to_cart
     use rundat, only: pub_geom_output_detail

     implicit none

     type(MODEL), intent(in)   :: mdl  ! jcap
     type(castep_model), intent(in) :: cmdl

     !local variables
     integer                        :: i,j,k
     real(kind=DP)                  :: cart(3)

     if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_output'

     ! Do it all on root proc
     ! lk: write out all coordinates in output file
     !     only if geom_output_detail >= VERBOSE
     if(pub_on_root .and. pub_geom_output_detail >= VERBOSE)then


       if (.not.fix_all_ions) then

          write(stdout,*)
          write(stdout,*)'                          -------------------&
     &------------'
          write(stdout,*)'                                    Cell Cont&
     &ents'
          write(stdout,*)'                          -------------------&
     &------------'
          write(stdout,*)

          ! write the ionic positions and species in cell
          write(stdout,7)
          write(stdout,*)'           x  Element    Atom      Absolute &
     &co-ordinates of atoms     x'
          write(stdout,*)'           x            Number        x      &
     &     y           z       x'
          write(stdout,8)


          do i=1,cmdl%cell%num_species - mdl%nat_classical
             do j=1,cmdl%cell%num_ions_in_species(i)
!             do j=1,1
                !NB Model is no longer rationalised and we don't bother to make it so here!
                !                write(stdout,3) cmdl%cell%species_symbol(i), j, (cmdl%cell%ionic_positions(k,j,i),k=1,3)
                !               write(stdout,3) cmdl%cell%species_symbol(i), i, (cmdl%cell%ionic_positions(k,j,i),k=1,3)
                ! aam: convert to absolute cartesian co-ordinates before outputting
                call castep_cell_frac_to_cart(current_cell, &
     &               cmdl%cell%ionic_positions(:,j,i),cart)
                write(stdout,3) cmdl%cell%species_symbol(i), i, &
     &(cart(k),k=1,3)
             end do
          end do
          write(stdout,7)
3         format(11x,' x',a6,1x,i8,4x,3f12.6,'   x')
7         format(1x,'           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx&
     &xxxxxxxxxxxxxxxxxxxxx')
8         format(1x,'           x--------------------------------------&
     &--------------------x')

       end if

       write(stdout,*)

       ! Back to all procs
     endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_output'

    return
  end subroutine geom_output

  subroutine geom_write_trajectory(cmdl,enthalpy,output_file,mdl)
    !=========================================================================!
    ! Output relevant geom_data in AU to specified trajectory file.           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=in, the model to be written                               !
    !   enthalpy, intent=in, the current enthalpy                             !
    !   output_file, intent=in, the filename to which results will be written.!
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   on_root                                                               !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    ! kaw: Gives access to the number of classical atoms.
    use comms, only: pub_on_root
    use constants, only: dp, stdout
    use model_type, only: MODEL
    use simulation_cell, only: castep_model, castep_cell_frac_to_cart
    use utils, only : utils_unit, utils_open_unit_check

    implicit none

    type(MODEL), intent(in)         :: mdl
    type(castep_model), intent(in)  :: cmdl
    real(kind=DP),    intent(in)    :: enthalpy
    character(len=*),  intent(in)   :: output_file

    !local variables
    integer                         :: i,iatom,ispec
    real(kind=DP), dimension(1:3)   :: frac_pos, abs_pos
    character(len=12)               :: action,form,stat,position,access

    integer                         :: kk

    !for filehandling we need
    integer                         :: out_unit, io_status

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_write_trajectory'

    action   = 'WRITE'
    form     = 'FORMATTED'
    stat     = 'UNKNOWN'
    position = 'REWIND'
    access   = 'SEQUENTIAL'

    if (iteration==0) then
       !create a new output file
       stat     ='REPLACE'
    else
       !append to existing output file
       position ='APPEND'
    end if

    if (pub_on_root) then

       ! Find next available unit specifier
       out_unit = utils_unit()

       open(unit=out_unit,iostat=io_status,file=output_file,status=stat,&
            access=access,form=form,position=position,action=action)
       call utils_open_unit_check('geom_write_trajectory','output_file',io_status)

       !now write out all the relevant data in atomic units (root proc only)
       write (out_unit,1) iteration
       write (out_unit,2) cmdl%total_energy,enthalpy

       !cell vectors
       do i=1,3
          write (out_unit,4) cmdl%cell%real_lattice(i,:)
       end do

       !ion stuff
       if (.not.fix_all_ions) then

          !ionic positions
          do ispec=1,cmdl%cell%num_species - mdl%nat_classical
             do iatom=1,cmdl%cell%num_ions_in_species(ispec)
                frac_pos(1:3)=cmdl%cell%ionic_positions(1:3,iatom,ispec)
                ! rationalise fractional co-ordinates before
                ! converting to cartesians and writing to file
                do kk=1,3
                   if (frac_pos(kk).lt.0.0_dp) then
                      frac_pos(kk)=frac_pos(kk)+1.0_dp
                   else if (frac_pos(kk).gt.1.0_dp) then
                      frac_pos(kk)=frac_pos(kk)-1.0_dp
                   else
                      continue
                   endif
                enddo
                call castep_cell_frac_to_cart(current_cell,frac_pos,abs_pos)
                write (out_unit,7) cmdl%cell%species_symbol(ispec),iatom, abs_pos(:)
             end do
          end do

          !ionic forces
    ! kaw: modified this to allow for classical atoms.
          do ispec=1,cmdl%cell%num_species - mdl%nat_classical
             do iatom=1,cmdl%cell%num_ions_in_species(ispec)
                write (out_unit,9)  cmdl%cell%species_symbol(ispec),iatom, cmdl%forces(:,iatom,ispec)
             end do
          end do

       end if

       !blank line to signal end of this iter (datablock size variable)
       write (out_unit,*) ' '

       close (unit=out_unit)

    end if

    !NB Try to keep format numbers/layout same here as in md and fit into 80-char page-width
1   format(12x,i18)
2   format(9x,2(3x,es18.8e3),21x,'  <-- E')
4   format(9x,3(3x,es18.8e3),    '  <-- h')
7   format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- R')
9   format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- F')

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_write_trajectory'

    return
  end subroutine geom_write_trajectory

  subroutine lbfgs_restart(toggle,init_cell,cmdl,mdl)
    !  jd (retreat2014): This 'toggle' thing is apparently ignored. I'm adding
    !                    a utils_use_var() to silence it, but I'm not removing
    !                    it in case it's more than a vestige from CASTEP.

    !=========================================================================!
    ! Prepare the LBFGS intermediate matrices after a restart                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=inout, the model to be optimised                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   algor for matrix inversion                                            !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !

    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !

    !-------------------------------------------------------------------------!
    ! Written by Jolyon Aarons, v1.0, 15/02/2010                              !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use comms, only: comms_bcast
    use constants, only: dp, stdout
    use linalg, only : linalg_invert_sym_matrix
    use model_type, only: MODEL ! jcap
    use rundat, only: pub_geom_lbfgs_max_updates
    use simulation_cell, only: castep_model, castep_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_use_var

    ! Arguments
    implicit none
    type(castep_model)                            :: cmdl ! @fixme intent?
    type(MODEL), intent(in)                       :: mdl  ! jcap
    type (castep_cell),   intent(in)                :: init_cell
    integer,            intent(in)                :: toggle

    ! BLAS subroutines
    external :: dgemv, dsbmv, dger

    ! Local variables
    integer                                       :: current_m
    integer                                       :: i
    integer,save                                  :: count
    !logical                                       :: inv=.false.
    integer                                       :: ierr

    !Householder variables
    real(kind=dp), dimension(:), allocatable      :: H_vector, HH_vector
    real(kind=dp)                                 :: alpha, r

    ! Initialize the scalars.
    !    cmdl%ndim     = ndim

    call utils_use_var(toggle)

    count=count+1

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering lbfgs_restart'

    if(cmdl%bfgs_optimisation) then
       ! Use the stored Inverse Hessian as an initialisation by tridiagonalising it with a Householder transformation.
       ! Could probably improve this by pentadiagonalising instead.
       allocate(H_vector(ndim),stat=ierr)
       call utils_alloc_check('lbfgs_restart','H_vector',ierr)
       allocate(HH_vector(ndim),stat=ierr)
       call utils_alloc_check('lbfgs_restart','HH_vector',ierr)

       do i=1,ndim-2
          H_vector=cmdl%bfgs_inv_Hessian(:,i)
          alpha=dot_product(H_vector(i+1:ndim),H_vector(i+1:ndim))
          alpha=sign(sqrt(alpha),-H_vector(i+1))
          r=sqrt(0.5_dp*(alpha**2 - H_vector(i+1)*alpha))

          ! ja531 -> March 2015 : avoid NaNs if we have a zero-length Householder vector
          if(abs(r)<epsilon(r)) cycle

          H_vector(1:i)=0.0_dp
          H_vector(i+1)=(H_vector(i+1)-alpha)/(2.0_dp*r)
          H_vector(i+2:ndim)=H_vector(i+2:ndim)/(2.0_dp*r)

          call dgemv('N', ndim, ndim, 1.0_dp, cmdl%bfgs_inv_Hessian, &
               size(cmdl%bfgs_inv_Hessian,1), H_vector, 1, 0.0_dp,HH_vector, 1)
          call dger(ndim, ndim, -2.0_dp, H_vector, 1, HH_vector, 1, &
               cmdl%bfgs_inv_Hessian, size(cmdl%bfgs_inv_Hessian,1))
          call dgemv('N', ndim, ndim, 1.0_dp, cmdl%bfgs_inv_Hessian, &
               size(cmdl%bfgs_inv_Hessian,1), H_vector, 1, 0.0_dp,HH_vector, 1)
          call dger(ndim, ndim, -2.0_dp, HH_vector, 1, H_vector, 1, &
               cmdl%bfgs_inv_Hessian, size(cmdl%bfgs_inv_Hessian,1))
       end do

       deallocate(H_vector,stat=ierr)
       call utils_dealloc_check('lbfgs_restart','H_vector',ierr)
       deallocate(HH_vector,stat=ierr)
       call utils_dealloc_check('lbfgs_restart','HH_vector',ierr)

       ! ja531 -> March 2015 : Get the order of indices the right way around.
       lbfgs_init_inv_hess(1,:)=0.0_dp
       lbfgs_init_inv_hess(2,1)=0.0_dp
       lbfgs_init_inv_hess(3,1)=cmdl%bfgs_inv_Hessian(1,1)
       do i=2,ndim
          lbfgs_init_inv_hess(2:3,i)=cmdl%bfgs_inv_Hessian(i,i-1:i)
       end do
       deallocate(cmdl%bfgs_inv_Hessian,stat=ierr)
       call utils_dealloc_check('lbfgs_restart','=>cmdl%bfgs_inv_Hessian',ierr)
       cmdl%bfgs_optimisation=.false.
       cmdl%lbfgs_optimisation=.true.
       cmdl%lbfgs_num_updates=0
    else
       ! Have to find out how much was stored in the model file and whether
       ! that fits with the history length parameters for this run.
       ! ja531 -> March 2015: I don't think we need these lbfgs_resizes.
       if(pub_geom_lbfgs_max_updates.le.0) then
          current_m=cmdl%lbfgs_num_updates+lbfgs_block_length
!          call lbfgs_resize(cmdl,current_m)
       else if(pub_geom_lbfgs_max_updates.ge.cmdl%lbfgs_num_updates) then
          current_m=pub_geom_lbfgs_max_updates
!          call lbfgs_resize(cmdl,current_m)
       else
          ! more updates in model file than allowed by new geom_lbfgs_max_updates --> delete some updates from the start.
          call lbfgs_removepair(cmdl,shift=(cmdl%lbfgs_num_updates-pub_geom_lbfgs_max_updates))
          current_m=pub_geom_lbfgs_max_updates
       end if

       if(allocated(R_mat)) then
          deallocate(R_mat,stat=ierr)
          call utils_dealloc_check('lbfgs_restart','R_mat',ierr)
       end if
       if(allocated(YY_mat)) then
          deallocate(YY_mat,stat=ierr)
          call utils_dealloc_check('lbfgs_restart','YY_mat',ierr)
       end if
       if(allocated(tmpmm)) then
          deallocate(tmpmm,stat=ierr)
          call utils_dealloc_check('lbfgs_restart','tmpmm',ierr)
       end if
       if(allocated(tmp2m)) then
          deallocate(tmp2m,stat=ierr)
          call utils_dealloc_check('lbfgs_restart','tmp2m',ierr)
       end if


       ! ja531 -> March 2015: We need this too...
       if(allocated(lbfgs_init_inv_hess)) then
          deallocate(lbfgs_init_inv_hess,stat=ierr)
          call utils_dealloc_check('lbfgs_restart','lbfgs_init_inv_hess',ierr)
       end if
       if(allocated(lbfgs_init_hess)) then
          deallocate(lbfgs_init_hess,stat=ierr)
          call utils_dealloc_check('lbfgs_restart','lbfgs_init_hess',ierr)
       end if
       if(allocated(tmpn)) then
          deallocate(tmpn,stat=ierr)
          call utils_dealloc_check('lbfgs_restart','tmpn',ierr)
       end if
       allocate(lbfgs_init_inv_hess(1:3,1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','lbfgs_init_inv_hess',ierr)
       allocate(lbfgs_init_hess(1:3,1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','lbfgs_init_hess',ierr)
       allocate(tmpn(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','tmpn',ierr)


       allocate(R_mat(1:current_m,1:current_m),stat=ierr)
       call utils_alloc_check('lbfgs_restart','R_mat',ierr)
       allocate(YY_mat(1:current_m,1:current_m),stat=ierr)
       call utils_alloc_check('lbfgs_restart','YY_mat',ierr)
       allocate(tmpmm(1:current_m*2,1:current_m*2),stat=ierr)
       call utils_alloc_check('lbfgs_restart','tmpmm',ierr)
       allocate(tmp2m(1:current_m*2),stat=ierr)
       call utils_alloc_check('lbfgs_restart','tmp2m',ierr)

       !       if (constrain_cell.and..not.fix_all_cell) inv=.true.
       !       call lbfgs_update_initial_hessian(toggle,init_cell,cmdl,inv=inv)
       !       call geom_inv_Hessian_initialize(toggle, init_cell, cmdl)

       ! ja531 -> March 2015: Force l-BFGS.
       if((hess_init_method.eq.hess_init_identity).or.(hess_init_method.eq.hess_init_scaled)) then
          call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,init_method=hess_init_identity,force_opt_method="L")
       else
          call geom_inv_Hessian_initialize(toggle_on,init_cell,cmdl,mdl,force_opt_method="L")
       end if

       ! ja531-> lbfgs restart bugfix 03/2018
       D_mat=0.0_dp
       do i=1,cmdl%lbfgs_num_updates
          D_mat(i)=-dot_product(cmdl%lbfgs_position_updates(1:ndim,i),cmdl%lbfgs_gradient_updates(1:ndim,i))
       end do

       ! Update R = S'*Y.
       R_mat=0.0_dp
       do i=1,cmdl%lbfgs_num_updates
          call dgemv('t',ndim,i,1.0_dp,cmdl%lbfgs_position_updates(1,1),&
               &size(cmdl%lbfgs_position_updates,1),cmdl%lbfgs_gradient_updates(1,i),1,0.0_dp,R_mat(1,i),1)
       end do

       YY_mat(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)=0.0_dp
       call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_inv_hess,size(lbfgs_init_inv_hess,1),&
            &cmdl%lbfgs_gradient_updates(1,1),1,0.0_dp,tmpn,1)
       YY_mat(1,1)=dot_product(tmpn,cmdl%lbfgs_gradient_updates(1:ndim,1))+R_mat(1,1)

       do i=2,cmdl%lbfgs_num_updates
          call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_inv_hess,size(lbfgs_init_inv_hess,1),&
               &cmdl%lbfgs_gradient_updates(1,i),1,0.0_dp,tmpn,1)
          call dgemv('t',ndim,i-1,1.0_dp,cmdl%lbfgs_gradient_updates,&
               &size(cmdl%lbfgs_gradient_updates,1),tmpn,1,0.0_dp,YY_mat(1,i),1)
          YY_mat(i,1:i-1)=YY_mat(1:i-1,i)
          YY_mat(i,i)=dot_product(tmpn,cmdl%lbfgs_gradient_updates(1:ndim,i)) + R_mat(i,i)
       end do

    end if
    !    deallocate(tmpm)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving lbfgs_restart'

    return

  end subroutine lbfgs_restart

  subroutine lbfgs_removepair (cmdl,shift)
    !=========================================================================!
    ! This function removes the first correction pair in storage, which       !
    ! effectively means that it removes the first column of S and Y.          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=inout, the model to be optimised                          !
    !   shift, optional, intent=in, the number of updates to remove           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    ! shift is not sign checked... but it probably ought to be positive!      !
    !-------------------------------------------------------------------------!
    ! Written by Jolyon Aarons, v1.0, 15/04/2011                              !
    !=========================================================================!
    use constants, only: stdout
    use utils, only: utils_alloc_check, utils_dealloc_check
    use comms, only: comms_bcast
    use simulation_cell, only: castep_model
    type(castep_model)                                 :: cmdl
    integer, optional, intent(in)                     :: shift
    integer                                           :: lshift

    integer,save :: count
    count=count+1

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering lbfgs_removepair'

    lshift=1
    if(present(shift)) lshift=shift

    cmdl%lbfgs_num_updates = cmdl%lbfgs_num_updates - 1

    ! Remove the first column from S and Y.
    cmdl%lbfgs_position_updates = eoshift(cmdl%lbfgs_position_updates, shift=lshift, dim=2)
    cmdl%lbfgs_gradient_updates = eoshift(cmdl%lbfgs_gradient_updates, shift=lshift, dim=2)

    ! Remove the first column and row from D_mat, SS and L.
    D_mat  = eoshift(D_mat, shift=lshift)
    YY_mat = eoshift(eoshift(YY_mat, shift=lshift, dim=2),shift=lshift)
    R_mat  = eoshift(eoshift(R_mat, shift=lshift, dim=2),shift=lshift)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving lbfgs_removepair'

    return

  end subroutine lbfgs_removepair

  subroutine lbfgs_resize(cmdl,new_size)
    !=========================================================================!
    ! Resizes all lbfgs matrices along cmdl%lbfgs_num_updates dimension to be  !
    ! new_size. Handles enlargement and reduction.                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=inout, the model to be optimised                          !
    !   new_size, intent=in, the new number of updates to store in lbfgs      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Jolyon Aarons, v1.0, 15/04/2011                              !
    !=========================================================================!
    use comms, only: comms_bcast
    use constants, only: dp, stdout
    use simulation_cell, only: castep_model
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    ! Arguments
    type(castep_model)                             :: cmdl ! @fixme intent?
    integer,intent(in)                            :: new_size

    ! BLAS subroutine
    external :: dcopy

    ! Local variables
    integer                                       :: old_size
    real(kind=dp), dimension(:,:), allocatable    :: old_s, old_y, old_YY_mat, old_R
    real(kind=dp), dimension(:), allocatable      :: old_D_mat
    integer                                       :: ierr

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering lbfgs_resize'
    old_size = size(cmdl%lbfgs_position_updates,2)

    if(new_size.lt.old_size) then
       if(cmdl%lbfgs_num_updates.gt.new_size) then
          call lbfgs_removepair(cmdl,shift=(cmdl%lbfgs_num_updates-new_size))
       end if
       old_size=cmdl%lbfgs_num_updates
    end if

    ! Start to re-allocate arrays, backing up along the way and keeping temporary arrays
    ! in memory for the least possible time to reduce high-water mark.

    ! re-allocate cmdl%lbfgs_position_updates and cmdl%lbfgs_gradient_updates
    allocate(old_s(ndim,old_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','old_s',ierr)
    allocate(old_y(ndim,old_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','old_y',ierr)
    !    call dcopy(ndim*old_size,cmdl%lbfgs_position_updates,1,old_s,1)
    !    call dcopy(ndim*old_size,cmdl%lbfgs_gradient_updates,1,old_y,1)
    old_s(1:ndim,1:old_size)=cmdl%lbfgs_position_updates(1:ndim,1:old_size)
    old_y(1:ndim,1:old_size)=cmdl%lbfgs_gradient_updates(1:ndim,1:old_size)
    deallocate(cmdl%lbfgs_position_updates,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','=>cmdl%lbfgs_position_updates',ierr)
    deallocate(cmdl%lbfgs_gradient_updates,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','=>cmdl%lbfgs_gradient_updates',ierr)
    allocate(cmdl%lbfgs_position_updates(ndim,new_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','=>cmdl%lbfgs_position_updates',ierr)
    allocate(cmdl%lbfgs_gradient_updates(ndim,new_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','=>cmdl%lbfgs_gradient_updates',ierr)
    cmdl%lbfgs_position_updates(1:ndim,1:old_size)=old_s(1:ndim,1:old_size)
    cmdl%lbfgs_gradient_updates(1:ndim,1:old_size)=old_y(1:ndim,1:old_size)
    cmdl%lbfgs_position_updates(1:ndim,old_size+1:new_size)=0.0_dp
    cmdl%lbfgs_gradient_updates(1:ndim,old_size+1:new_size)=0.0_dp
    deallocate(old_s,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','old_s',ierr)
    deallocate(old_y,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','old_y',ierr)

    ! re-allocate YY_mat
    allocate(old_YY_mat(old_size,old_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','old_YY_mat',ierr)
    !    call dcopy(old_size**2,YY_mat,1,old_YY_mat,1)
    old_YY_mat(1:old_size,1:old_size)=YY_mat(1:old_size,1:old_size)
    deallocate(YY_mat,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','YY_mat',ierr)
    allocate(YY_mat(new_size,new_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','YY_mat',ierr)
    !    call dcopy(old_size**2,old_YY_mat,1,YY_mat,1)
    YY_mat=0.0_dp
    YY_mat(1:old_size,1:old_size)=old_YY_mat(1:old_size,1:old_size)
    deallocate(old_YY_mat,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','old_YY_Mat',ierr)

    ! re-allocate tmpmm and tmp2m
    deallocate(tmpmm,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','tmpmm',ierr)
    allocate(tmpmm(new_size*2,new_size*2),stat=ierr)
    call utils_alloc_check('lbfgs_resize','tmpmm',ierr)
    deallocate(tmp2m,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','tmp2m',ierr)
    allocate(tmp2m(new_size*2),stat=ierr)
    call utils_alloc_check('lbfgs_resize','tmp2m',ierr)

    ! re-allocate D_mat
    allocate(old_D_mat(old_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','old_D_mat',ierr)
    call dcopy(old_size,D_mat,1,old_D_mat,1)
    old_D_mat(1:old_size)=D_mat(1:old_size)
    deallocate(D_mat,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','D_mat',ierr)
    allocate(D_mat(new_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','D_mat',ierr)
    call dcopy(old_size,old_D_mat,1,D_mat,1)
    D_mat=0.0_dp
    D_mat(1:old_size)=old_D_mat(1:old_size)
    deallocate(old_D_mat,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','old_D_mat',ierr)

    ! re-allocate R_mat
    allocate(old_R(old_size,old_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','old_R',ierr)
    !    call dcopy(old_size**2,R_mat,1,old_R,1)
    old_R(1:old_size,1:old_size)=R_mat(1:old_size,1:old_size)
    deallocate(R_mat,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','R_mat',ierr)
    allocate(R_mat(new_size,new_size),stat=ierr)
    call utils_alloc_check('lbfgs_resize','R_mat',ierr)
    !    call dcopy(old_size**2,old_R,1,R_mat,1)
    R_mat=0.0_dp
    R_mat(1:old_size,1:old_size)=old_R(1:old_size,1:old_size)
    deallocate(old_R,stat=ierr)
    call utils_dealloc_check('lbfgs_resize','old_R',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving lbfgs_resize'

    return
  end subroutine lbfgs_resize

  subroutine lbfgs_update_bandmatrix(cmdl,update_vecs,lbfgs_recurs_mat,H0,sd_bandwidth)
    !=========================================================================!
    ! Updates the initialisation to the inverse Hessian by using a            !
    ! Householder transformation applied to the lbfgs update space matrices   !
    ! to reflect into a band-diagonal form to sum with the previous pre-      !
    ! conditioner. (destroys update_vecs)                                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   update_vecs, the lbfgs matrix of updates to the positions and grads   !
    !   lbfgs_recurs_mat, the lbfgs recursion matrix (M)                      !
    !   H0, the band-diagonal symmetric matrix in LAPACK form                 !
    !       --> i.e. ndim*(1+sd_bandwidth)                                    !
    !  optional:
    !   sd_bandwidth, the desired superdiagonal bandwidth (i.e. 2=pentadiag)  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Jolyon Aarons, v1.0, 29/07/2011                              !
    !=========================================================================!
    use simulation_cell, only: castep_model
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    ! Arguments
    type(castep_model),             intent(in)     :: cmdl ! only for cmdl%lbfgs_num_updates
    real(kind=dp), dimension(:,:), intent(inout)  :: update_vecs
    real(kind=dp), dimension(:,:), intent(in)     :: lbfgs_recurs_mat
    real(kind=dp), dimension(:,:), intent(inout)  :: H0
    integer,             optional, intent(in)     :: sd_bandwidth ! super-diagonal bandwidth ... assumes pentadiagonal

    ! BLAS subroutines
    external :: dgemm, dsbmv, dgemv, dger

    ! Local variables
    integer                                              :: i, bw, ierr
    real(kind=dp)                                        :: alpha,r
    real(kind=dp), dimension(ndim)                       :: H_vector
    real(kind=dp), dimension(cmdl%lbfgs_num_updates*2)    :: m_vector
    real(kind=dp), dimension(:,:), allocatable           :: left_mat

!    call trace_entry('lbfgs_update_bandmatrix',status)

    bw=2
    if(present(sd_bandwidth)) bw=sd_bandwidth

    !Perform a partial diagonalisation through a Householder transformation.
    do i=1,ndim-bw
       call dgemv('N', cmdl%lbfgs_num_updates*2, cmdl%lbfgs_num_updates*2, 1.0_dp, lbfgs_recurs_mat, size(lbfgs_recurs_mat,1),&
            update_vecs(i,1), size(update_vecs,1), 0.0_dp,m_vector, 1)
       call dgemv('N', ndim, cmdl%lbfgs_num_updates*2, 1.0_dp, update_vecs, size(update_vecs,1), m_vector, 1, 0.0_dp,H_vector, 1)

       alpha=dot_product(H_vector(i+bw-1:ndim),H_vector(i+bw-1:ndim))
       alpha=sign(sqrt(alpha),-H_vector(i+bw-1))
       r=sqrt(0.5_dp*(alpha**2 - H_vector(i+bw-1)*alpha))

       H_vector(1:i+bw-2)=0.0_dp
       H_vector(i+bw-1)=(H_vector(i+bw-1)-alpha)/(2.0_dp*r)
       H_vector(i+bw:ndim)=H_vector(i+bw:ndim)/(2.0_dp*r)

       call dgemv('T', ndim, cmdl%lbfgs_num_updates*2, 1.0_dp, update_vecs, size(update_vecs,1), H_vector, 1, 0.0_dp,m_vector, 1)
       call dger(ndim,cmdl%lbfgs_num_updates*2, -2.0_dp, H_vector, 1, m_vector, 1, update_vecs, size(update_vecs,1))
    end do

    !Put the result into the band matrix.
    !       do i=1,ndim
    !          H0(:,i)=H0(:,i)+matmul(update_vecs(i-bw:i,:),matmul(lbfgs_recurs_mat,transpose(update_vecs(i,:))))
    !       end do
    allocate(left_mat(ndim,2*cmdl%lbfgs_num_updates),stat=ierr)
    call utils_alloc_check('lbfgs_update_bandmatrix','left_mat',ierr)
    left_mat=0.0_dp
    call dgemm('n','n',ndim,2*cmdl%lbfgs_num_updates,2*cmdl%lbfgs_num_updates,1.0_dp,update_vecs,size(update_vecs,1),&
         &lbfgs_recurs_mat,size(lbfgs_recurs_mat,1),0.0_dp,left_mat,size(left_mat,1))
    do i=1,ndim-bw
       call dgemv('n',3,2*cmdl%lbfgs_num_updates,1.0_dp,left_mat(i,1),size(left_mat,1),update_vecs(i+bw,1),&
            &2*cmdl%lbfgs_num_updates,1.0_dp,H0(1,i+bw),1)
    end do
    deallocate(left_mat,stat=ierr)
    call utils_dealloc_check('lbfgs_update_bandmatrix','left_mat',ierr)

    !Sort out the first 1:bw columns of H0
    do i=1,bw
       H0(bw-i+2:1+bw,i)=H0(bw-i+2:1+bw,i)+matmul(update_vecs(1:i,:),matmul(lbfgs_recurs_mat,update_vecs(i,:)))
    end do

!    call trace_exit('lbfgs_update_bandmatrix',status)

  end subroutine lbfgs_update_bandmatrix


  subroutine lbfgs_initH_mul(trans,A,B,C,info)
    !=========================================================================!
    ! Multiply a General matrix by a symmetric penta-diagonal band matrix and !
    ! put the result in C.                                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   trans, whether or not to transpose A, 't'/'T' or 'n'/'N'              !
    !   A, the general matrix --> not necessarily square.                     !
    !   B, the penta-diagonal symmetric band matrix in LAPACK form            !
    !       --> i.e. ndim*3                                                   !
    !   C, intent(out), the resulting general matrix                          !
    !   info, intent(out), successful completion of work --> info=1           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Jolyon Aarons, v1.0, 15/04/2011                              !
    !=========================================================================!
    use constants, only: dp, stdout
    use comms, only: comms_bcast
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none
    character,                   intent(in)       :: trans
    real(kind=dp),dimension(:,:),intent(in)       :: A,B
    integer,                     intent(out)      :: info
    real(kind=dp),dimension(:,:),intent(out)      :: C

    integer                                       :: i,j,n,m,which

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering lbfgs_initH_mul'

    n=size(A,1)
    m=size(A,2)
    if(size(A,1).ne.size(B,2)) then
       info=1

       if (pub_debug_on_root) write(stdout,'(a)') &
            'DEBUG: Leaving lbfgs_initH_mul'

       return
    end if

    C=0.0_dp
    !    if(present(trans)) then
    if((trans.eq.'t').or.(trans.eq.'T')) then
       which=2
    else if((trans.eq.'n').or.(trans.eq.'N')) then
       which=1
    else
       info=2

       if (pub_debug_on_root) write(stdout,'(a)') &
            'DEBUG: Leaving lbfgs_initH_mul'

       return
    end if
    !    else
    !       which=1
    !    end if

    if(which.eq.1) then ! No transpose --> C=AB
       do i=1,m
          c(1,i) = c(1,i) + a(1,i)*b(3,1) + a(2,i)*b(2,2) + a(3,i)*b(1,3)
          c(2,i) = c(2,i) + a(1,i)*b(2,2) + a(2,i)*b(3,2) + a(3,i)*b(2,3) + a(4,i)*b(1,4)
          do j=3,n-2
             c(j,i) = c(j,i) + a(j-2,i)*b(1,j) + a(j-1,i)*b(2,j) + a(j  ,i)*b(3,j) + a(j+1,i)*b(2,j+1) + a(j+2,i)*b(1,j+2)
          end do
          c(n-1,i) = c(n-1,i) + a(n-3,i)*b(1,n-1) + a(n-2,i)*b(2,n-1) + a(n-1,i)*b(3,n-1) + a(n  ,i)*b(2,n)
          c(n,i) = c(n,i) + a(n-2,i)*b(1,n) + a(n-1,i)*b(2,n) + a(n,i)*b(3,n)
       end do

       info=0
    else if(which.eq.2) then ! transpose --> C=A'B or C'=BA
       do i=1,m
          c(i,1) = c(i,1) + a(1,i)*b(3,1) + a(2,i)*b(2,2) + a(3,i)*b(1,3)
          c(i,2) = c(i,2) + a(1,i)*b(2,2) + a(2,i)*b(3,2) + a(3,i)*b(2,3) + a(4,i)*b(1,4)
       end do
       do j=3,n-2
          do i=1,m
             c(i,j) = c(i,j) + a(j-2,i)*b(1,j) + a(j-1,i)*b(2,j) + a(j  ,i)*b(3,j) + a(j+1,i)*b(2,j+1) + a(j+2,i)*b(1,j+2)
          end do
       end do
       do i=1,m
          c(i,n-1) = c(i,n-1) + a(n-3,i)*b(1,n-1) + a(n-2,i)*b(2,n-1) + a(n-1,i)*b(3,n-1) + a(n  ,i)*b(2,n)
          c(i,n) = c(i,n) + a(n-2,i)*b(1,n) + a(n-1,i)*b(2,n) + a(n,i)*b(3,n)
       end do
       info=0
    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving lbfgs_initH_mul'

  end subroutine lbfgs_initH_mul


  function geom_apply_vec(cmdl,f_vec,inv) result(p)
    !=========================================================================!
    ! This function multiplies a vector by the BFGS hessian using a BLAS for  !
    ! BFGS and a bespoke routine for LBFGS. The resulting vector is returned  !
    ! in p.                                                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent(in)                                                       !
    !   f_vec, intent(in) : the input vector to multiply be Hessian apprx.    !
    !   inv, intent(in)   : whether or not we want inverse Hessian 'y' or 'n' !
    !   result --> p=Hv                                                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! tmp2m, tmpmm                                                            !
    ! hessian if doing inverse and BFGS                                       !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Jolyon Aarons, v1.0, 15/04/2011                              !
    !=========================================================================!

    use comms, only: comms_bcast
    use constants, only: dp, stdout
    use linalg, only: linalg_invert_sym_matrix
    use simulation_cell, only: castep_model
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(castep_model), intent(inout)              :: cmdl
    real(kind=dp), dimension(:), intent(in)       :: f_vec
    character(len=1)                              :: inv
    real(kind=dp), dimension(ndim)                :: p, v
    !    type (unit_cell),intent(in),optional          :: init_cell

    ! BLAS & LAPACK subroutines
    external :: dgemv, dsbmv, dswap, dtrsm, dsymv, dpotrf, dpotrs

    ! Local variables
    integer                                       :: i !,lbfgs_num_pairs
    real(kind=dp), dimension(:,:), allocatable    :: HS_mat !, tmpmm
    real(kind=dp), dimension(:), allocatable      :: tmp2m_2 !,tmpn,tmpm,tmp2m
    integer                                       :: status
    integer                                       :: ierr
    !    logical                                       :: woodbury
    integer,save                                  :: count

    real(kind=dp), dimension(:), allocatable      :: alpha, beta      ! lam81
    real(kind=dp), dimension(:,:), allocatable    :: H0               ! lam81
    real(kind=dp)                                 :: c_cell           ! lam81
    integer                                       :: info             ! lam81

    integer :: loc_ndim ! kkbd TODO: Adapt lbfgs search direction finding

    count=count+1
    p=0.0_dp

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_apply_vec'

    v=f_vec

    ! Do the matrix vector product in BFGS
    if(do_bfgs) then
       if((inv.eq.'y').or.(inv.eq.'Y')) then
          call dsymv('u',ndim,1.0_dp,cmdl%bfgs_inv_Hessian,&
               &size(cmdl%bfgs_inv_Hessian,1),f_vec,1,0.0_dp,p,1)
       else if((inv.eq.'n').or.(inv.eq.'N')) then
          ! If we're doing a non-inverse we might need to set that up in hessian
          ! but then just do a dsymv as before.
          if(allocated(hessian)) then
             if(size(hessian).eq.ndim**2) then
                call dsymv('u',ndim,1.0_dp,hessian,&
                     &size(hessian,1),f_vec,1,0.0_dp,p,1)
             else
                deallocate(hessian,stat=ierr)
                call utils_dealloc_check('geom_apply_vec','hessian',ierr)
                allocate(hessian(ndim,ndim),stat=ierr)
                call utils_alloc_check('geom_apply_vec','hessian',ierr)
                hessian=cmdl%bfgs_inv_hessian
                call linalg_invert_sym_matrix(hessian,ndim)
                call dsymv('u',ndim,1.0_dp,hessian,&
                     &size(hessian,1),f_vec,1,0.0_dp,p,1)
             end if
          else
             allocate(hessian(ndim,ndim),stat=ierr)
             call utils_alloc_check('geom_apply_vec','hessian',ierr)
             hessian=cmdl%bfgs_inv_hessian
             call linalg_invert_sym_matrix(hessian,ndim)
             call dsymv('u',ndim,1.0_dp,hessian,&
                  &size(hessian,1),f_vec,1,0.0_dp,p,1)
          end if
       end if

       !LBFGS
    else if(do_lbfgs .and. .not. precond%initialized) then ! lam81
       ! if we're applying to the inverse Hessian:
       wb_inv : if((inv.eq.'y').or.(inv.eq.'Y')) then
          if(cmdl%lbfgs_num_updates.gt.0) then

             tmpmm=0.0_dp
             tmp2m=0.0_dp


             call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_inv_hess,size(lbfgs_init_inv_hess,1),v,1,0.0_dp,p,1)

             ! set up S'v and Y'H_0.v
             call dgemv('t',ndim,cmdl%lbfgs_num_updates,1.0_dp,cmdl%lbfgs_position_updates,&
                  &size(cmdl%lbfgs_position_updates,1),v(1),1,0.0_dp,tmp2m(1),1)
             call dgemv('t',ndim,cmdl%lbfgs_num_updates,1.0_dp,cmdl%lbfgs_gradient_updates,&
                  &size(cmdl%lbfgs_gradient_updates,1),p,1,0.0_dp,tmp2m(cmdl%lbfgs_num_updates+1),1)

             ! solve linear systems to find R'S'v and R'Y'H_0.v
             call dtrsm('l','u','n','n',cmdl%lbfgs_num_updates,1,1.0_dp,R_mat,&
                  &size(R_mat,1),tmp2m(1),size(tmp2m,1))
             call dtrsm('l','u','t','n',cmdl%lbfgs_num_updates,1,1.0_dp,R_mat,&
                  &size(R_mat,1),tmp2m(cmdl%lbfgs_num_updates+1),size(tmp2m,1))
             call dswap(cmdl%lbfgs_num_updates,tmp2m(1),1,tmp2m(cmdl%lbfgs_num_updates+1),1)

             tmpmm(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)=&
                  &YY_mat(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)

             ! solve linear system to find R'B and multiply out to find top half of vector
             call dtrsm('l','u','t','n',cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates,1.0_dp,R_mat,&
                  &size(R_mat,1),tmpmm,size(tmpmm,1))
             call dgemv('n',cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates,1.0_dp,tmpmm,&
                  &size(tmpmm,1),tmp2m(cmdl%lbfgs_num_updates+1),1,-1.0_dp,tmp2m(1),1)

             tmp2m(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)=-tmp2m(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)
             ! pre-multiply by U' and add sigma*H_0.v to finish product.
             !       call dcopy(ndim,v,1,p,1)

             call dgemv('n',ndim,cmdl%lbfgs_num_updates,1.0_dp,cmdl%lbfgs_position_updates,&
                  &size(cmdl%lbfgs_position_updates,1),tmp2m,1,1.0_dp,p,1)
             tmpn=0.0_dp
             call dgemv('n',ndim,cmdl%lbfgs_num_updates,1.0_dp,cmdl%lbfgs_gradient_updates,&
                  &size(cmdl%lbfgs_gradient_updates,1),tmp2m(cmdl%lbfgs_num_updates+1),1,0.0_dp,tmpn,1)
             call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_inv_hess,size(lbfgs_init_inv_hess,1),tmpn,1,1.0_dp,p,1)

          else
             p=0.0_dp
             call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_inv_hess,size(lbfgs_init_inv_hess,1),v,1,0.0_dp,p,1)
          end if
          !If we're applying to the Hessian:
       else if((inv.eq.'n').or.(inv.eq.'N')) then
          if(cmdl%lbfgs_num_updates.gt.0) then

             ! Same as above... but using the Woodbury inverse!

             tmpmm=0.0_dp
             tmp2m=0.0_dp

             allocate(HS_mat(1:ndim,1:cmdl%lbfgs_num_updates),stat=ierr)
             call utils_alloc_check('geom_apply_vec','HS_mat',ierr)
             HS_mat=0.0_dp
             call lbfgs_initH_mul('n',cmdl%lbfgs_position_updates(1:ndim,1:cmdl%lbfgs_num_updates),lbfgs_init_hess,&
                  &HS_mat(1:ndim,1:cmdl%lbfgs_num_updates),status)

             tmpmm(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)=&
                  &+matmul(transpose(cmdl%lbfgs_position_updates(1:ndim,1:cmdl%lbfgs_num_updates)),HS_mat)
             tmpmm(1:cmdl%lbfgs_num_updates,cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2)=&
                  &transpose(R_mat(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates))
             tmpmm(cmdl%lbfgs_num_updates+1:cmdl%lbfgs_num_updates*2,1:cmdl%lbfgs_num_updates)=&
                  &R_mat(1:cmdl%lbfgs_num_updates,1:cmdl%lbfgs_num_updates)

             do i=1,cmdl%lbfgs_num_updates
                tmpmm(cmdl%lbfgs_num_updates+i,i)=0.0_dp
                tmpmm(i,cmdl%lbfgs_num_updates+i)=0.0_dp
                tmpmm(cmdl%lbfgs_num_updates+i,cmdl%lbfgs_num_updates+i)=-D_mat(i)
             end do

             !call outputdata(tmpmm,unit=69)
             call linalg_invert_sym_matrix(tmpmm,2*cmdl%lbfgs_num_updates)

             !call outputdata(tmpmm,unit=69)

             call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_hess,size(lbfgs_init_hess,1),v,1,0.0_dp,p,1)

             ! set up S'v and Y'v
             call dgemv('t',ndim,cmdl%lbfgs_num_updates,1.0_dp,cmdl%lbfgs_position_updates,&
                  &size(cmdl%lbfgs_position_updates,1),p,1,0.0_dp,tmp2m,1)
             call dgemv('t',ndim,cmdl%lbfgs_num_updates,1.0_dp,cmdl%lbfgs_gradient_updates,&
                  &size(cmdl%lbfgs_gradient_updates,1),v,1,1.0_dp,tmp2m(cmdl%lbfgs_num_updates+1),1)

             allocate(tmp2m_2(1:cmdl%lbfgs_num_updates*2),stat=ierr)
             call utils_alloc_check('geom_apply_vec','tmp2m_2',ierr)
             tmp2m_2=0.0_dp
             call dgemv('n',cmdl%lbfgs_num_updates*2,cmdl%lbfgs_num_updates*2,1.0_dp,tmpmm,size(tmpmm,1),tmp2m,1,1.0_dp,tmp2m_2,1)
             call dgemv('n',ndim,cmdl%lbfgs_num_updates,-1.0_dp,HS_mat,size(HS_mat,1),tmp2m_2,1,1.0_dp,p,1)
             call dgemv('n',ndim,cmdl%lbfgs_num_updates,-1.0_dp,cmdl%lbfgs_gradient_updates,size(cmdl%lbfgs_gradient_updates,1),&
                  &tmp2m_2(cmdl%lbfgs_num_updates+1),1,1.0_dp,p,1)

             deallocate(HS_mat,stat=ierr)
             call utils_dealloc_check('geom_apply_vec','HS_mat',ierr)
             deallocate(tmp2m_2,stat=ierr)
             call utils_dealloc_check('geom_apply_vec','tmp2m',ierr)

             !             deallocate(tmpmm)
             !             deallocate(tmp2m)

          else
             p=0.0_dp
             call dsbmv('u',ndim,2,1.0_dp,lbfgs_init_hess,size(lbfgs_init_hess,1),v,1,0.0_dp,p,1)
          end if
       end if wb_inv
    else if (do_lbfgs .and. precond%initialized) then ! lam81

       allocate(alpha(cmdl%lbfgs_num_updates),stat=ierr)
       call utils_alloc_check('geom_apply_vec','alpha',ierr)
       allocate(beta(cmdl%lbfgs_num_updates),stat=ierr)
       call utils_alloc_check('geom_apply_vec','beta',ierr)
       allocate(H0(1:ndim, 1:ndim),stat=ierr)
       call utils_alloc_check('geom_apply_vec','H0',ierr)

       v = -v
       do i=cmdl%lbfgs_num_updates,1,-1
          alpha(i) = dot_product(cmdl%lbfgs_position_updates(:,i), v) / D_mat(i)
          v = v - alpha(i)*cmdl%lbfgs_gradient_updates(:,i)
       end do

       call geom_precond_calc(cmdl)

       H0 = 0.0_dp
       if ( precond%scale_cell .and. (cmdl%lbfgs_num_updates > 0) ) then
          c_cell = (dot_product(cmdl%lbfgs_gradient_updates(:,cmdl%lbfgs_num_updates), &
                              & cmdl%lbfgs_position_updates(:,cmdl%lbfgs_num_updates)) &
                    -dot_product(matmul(precond%P, cmdl%lbfgs_position_updates(10:,cmdl%lbfgs_num_updates)), &
                               & cmdl%lbfgs_position_updates(10:,cmdl%lbfgs_num_updates))) / &
                   dot_product(cmdl%lbfgs_position_updates(1:9,cmdl%lbfgs_num_updates), &
                             & cmdl%lbfgs_position_updates(1:9,cmdl%lbfgs_num_updates))
       else
          c_cell = 3.0_dp*cmdl%ref_cell%volume*new_modulus_est
       end if
       if (c_cell <= 0.0_dp) c_cell = 3.0_dp*cmdl%ref_cell%volume*new_modulus_est
       if (c_cell > 3.0_dp*cmdl%ref_cell%volume*new_modulus_est) c_cell = 3.0_dp*cmdl%ref_cell%volume*new_modulus_est
       do i=1,9
          H0(i,i) = c_cell
       end do
       H0(10:ndim,10:ndim) = precond%P

       p = v
       call dpotrf('u', ndim, H0, ndim, info)
       if (info /= 0) call utils_abort('Error dpotrf in geom_apply_vec')
       call dpotrs('u', ndim, 1, H0, ndim, p, ndim, info)
       if (info /= 0) call utils_abort('Error dpotrs in geom_apply_vec')

       do i=1,cmdl%lbfgs_num_updates
          beta(i) = dot_product(cmdl%lbfgs_gradient_updates(:,i),p) / D_mat(i)
          p = p + cmdl%lbfgs_position_updates(:,i)*(alpha(i)-beta(i))
       end do

       p = -p

       deallocate(alpha,stat=ierr)
       call utils_dealloc_check('geom_apply_vec','alpha',ierr)
       deallocate(beta,stat=ierr)
       call utils_dealloc_check('geom_apply_vec','beta',ierr)
       deallocate(H0,stat=ierr)
       call utils_dealloc_check('geom_apply_vec','H0',ierr)

    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_apply_vec'

    return

  end function geom_apply_vec


  function geom_opt_min_image_frac(f1,f2)
    !=========================================================================!
    ! Handles minimum image convention for periodic boundary conditions by    !
    ! applying translation vectors t of Bravais lattice to minimize r12, where!
    ! f12 = |f1 - f2 - t| for input f1 & f2 vectors, for any Bravais lattice. !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   f1,   intent=in, fractional input vector1                             !
    !   f2,   intent=in, fractional input vector2                             !
    ! Returns:                                                                !
    !   fractional vector f12=f1-f2                                           !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use constants, only: dp
    use simulation_cell, only:castep_cell_frac_to_cart, castep_cell_cart_to_frac

    implicit none

    real(kind=DP),dimension(3), intent(in) :: f1,f2    ! Fractional coords
    real(kind=DP), dimension(3) :: geom_opt_min_image_frac

    ! <<< local variables >>>
    real(kind=DP), dimension(3) :: r1,r2,r12,f12

    ! Get absolute Cartesian co-ordinates
    call castep_cell_frac_to_cart(current_cell,f1,r1)
    call castep_cell_frac_to_cart(current_cell,f2,r2)

    ! Get minimum image
    r12 = geom_opt_min_image_cart(r1,r2)

    ! Set to fractional
    call castep_cell_cart_to_frac(current_cell,r12,f12)

    ! Set function
    geom_opt_min_image_frac=f12

    return
  end function geom_opt_min_image_frac


  function geom_opt_min_image_cart(r1,r2)
    !=========================================================================!
    ! Handles minimum image convention for periodic boundary conditions by    !
    ! applying translation vectors t of Bravais lattice to minimize r12, where!
    ! r12 = |r1 - r2 - t| for input r1 & r2 vectors, for any Bravais lattice. !
    ! Routine taken from SCAMPI (mijp's path integral MD code)                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   r1,   intent=in, Cartesian input vector1                              !
    !   r2,   intent=in, Cartesian input vector2                              !
    ! Returns:                                                                !
    !   Cartesian vector r12=r1-r2                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! current_cell                                                            !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    ! parameters, io                                                          !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   transvec = the 8 different translation vectors = 't' in above         !
    !            = 0, a, b, c, a+b, a+c, b+c, a+b+c                           !
    !   r12      = the basic r1-r2 vector in the unit cell                    !
    !   mod_r12  = the 8 different values of |r12-t| in the unit cell         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.2, 05/01/2001                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use constants, only: dp
    use simulation_cell, only: castep_cell_cart_to_frac

    implicit none

    real(kind=DP), dimension(1:3), intent(in) :: r1, r2 !Cartesian
    real(kind=DP), dimension(1:3)             :: geom_opt_min_image_cart

    !local variables
    real(kind=DP), dimension(1:3)             :: r12
    real(kind=DP), dimension(1:3)             :: fract_r12
    real(kind=DP), dimension(1:8,1:3)         :: transvec
    real(kind=DP), dimension(1:8)             :: mod_r12
    real(kind=DP), dimension(1:3)             :: min_r12
    real(kind=DP)                             :: min_mod
    integer                                   :: i,j

    !Calculate basic r12=r1-r2 displacement vector ...
    r12=r1-r2

    call castep_cell_cart_to_frac(current_cell,r12,fract_r12)

    !Find origin of cell containing r12 ...
    do i=1,3
       if (fract_r12(i)>=0.0_dp) then
          fract_r12(i)=aint(fract_r12(i),kind=dp)
       else
          fract_r12(i)=aint(fract_r12(i)-1.0_dp,kind=dp)
       end if
    end do

    !... and map r12 s.t. new r12(i)= new fract_r12(i)*real_lattice(i),
    !where 0 <= new fract_r12 < 1 ...
    do i=1,3
       do j=1,3
          r12(i)=r12(i)-fract_r12(j)*current_cell%real_lattice(j,i)
       end do
    end do

    !... and now calculate translation vectors from each of eight corners
    !(of parallelepiped that contains r12) to r12 ...

    transvec(1,:)=r12(:)
    transvec(2,:)=r12(:)-current_cell%real_lattice(1,:)
    transvec(3,:)=r12(:)-current_cell%real_lattice(2,:)
    transvec(4,:)=r12(:)-current_cell%real_lattice(3,:)
    transvec(5,:)=r12(:)-current_cell%real_lattice(1,:)-current_cell%real_lattice(2,:)
    transvec(6,:)=r12(:)-current_cell%real_lattice(1,:)-current_cell%real_lattice(3,:)
    transvec(7,:)=r12(:)-current_cell%real_lattice(2,:)-current_cell%real_lattice(3,:)
    transvec(8,:)=r12(:)-current_cell%real_lattice(1,:)-current_cell%real_lattice(2,:)-current_cell%real_lattice(3,:)

    mod_r12=0.0_dp
    do i=1,3
       mod_r12(:)=mod_r12(:)+transvec(:,i)**2
    end do

    !... and finally chose the one with minimum length ...
    min_mod=mod_r12(1)
    min_r12(:)=transvec(1,:)
    do i=2,8
       if (mod_r12(i)<min_mod) then
          min_mod=mod_r12(i)
          min_r12(:)=transvec(i,:)
       end if
    end do

    ! Set function
    geom_opt_min_image_cart=min_r12

    return
  end function geom_opt_min_image_cart




  subroutine geom_opt_continuation_read(cmdl,mdl,filename,is_cmplx)
    !=========================================================================!
    ! Read all relevant data from file for geometry optimisation continuation !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=inout, the model                                         !
    !   mdl, intent=in, contains the element data array                       !
    !   filename, intent=in, the filename from which continuation             !
    !                           data will be read.                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !    model type                                                           !
    !    ndim                                                                 !
    !    current_cell                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !    rundat for geom_modulus_est, geom_frequency_est                      !
    !    comms for pub_on_root                                                !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   cmdl has already been properly read and initialised                   !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 08/02/2005                                  !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    ! Modified to remove elements argument by Robert Charlton, October 2018.  !
    !=========================================================================!

    use comms, only: pub_on_root, comms_barrier, comms_bcast, pub_root_proc_id
    use constants, only: stdout
    use model_type, only: MODEL
    use rundat, only: pub_geom_modulus_est,pub_geom_frequency_est, pub_rootname, &
         pub_maxit_ngwf_cg, pub_geom_reuse_dk_ngwfs,pub_edft, pub_num_spins, &
         pub_read_real_ngwfs
    use simulation_cell, only: castep_model, castep_model_alloc, &
         castep_model_dealloc
    use utils, only: utils_unit, utils_read_check, utils_open_unit_check, &
         utils_close_unit_check, utils_read_check, utils_binary_copy, &
         utils_abort

    implicit none

    ! Arguments
    type(castep_model), intent(inout)       :: cmdl
    type(MODEL), intent(in)                 :: mdl  ! jcap
    character(len=file_maxpath), intent(in) :: filename
    logical, intent(in)                     :: is_cmplx

    ! Local Variables
    integer :: read_unit,write_unit
    integer :: ios,ierr
    integer :: atom,ngwf,dim
    integer :: ionidx
    integer :: num,tn1,tn2,tn3,tn(3)
    integer :: tmp_num_spins
    integer :: isub
#ifdef ACCELRYS
    integer :: row           ! vm: required for Intel bug workaround
#endif
    real(kind=DP) :: orig1,orig2,orig3
    real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox
    character(len=file_maxpath) :: ngwf_file,dk_file,ham_file
    character(len=5) :: method_string
#ifdef HDF5
    integer(hid_t) :: filesystem_id
    integer(hid_t) :: file_id, file2_id
    integer(hid_t) :: group_id
    integer(hid_t) :: attr_id
    integer(hid_t) :: dset_id
    integer(hsize_t) :: adims(1), dims(2)
#else
    logical :: read_cmplx
#endif

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_opt_continuation_read'

    if (pub_on_root) then

#ifdef HDF5
       call h5open_f(ierr)
       call h5fopen_f(trim(filename),h5f_acc_rdonly_f,filesystem_id,ierr)
       call h5gopen_f(filesystem_id,'continuation',file_id,ierr)
#else
       ! Find available unit specifier for reading continuation file
       read_unit = utils_unit()

       ! Open continuation file
       open(unit=read_unit,file=filename,iostat=ios,&
            form='UNFORMATTED',action='READ',status='OLD')
       call utils_open_unit_check('geom_opt_continuation_read',filename,ios)
       rewind(read_unit)
#endif

       ! Decide whether continuation has lbfgs or bfgs data...
#ifdef HDF5
       adims=(/1/)
       call h5aopen_f(file_id,'geom_string',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_character,method_string(1:1), &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)
#else
       read(read_unit,iostat=ios) method_string
       call utils_read_check('geom_opt_continuation_read','method_string',ios)
#endif
       if((method_string(1:1).eq.'b').or.(method_string(1:1).eq.'B')) then
          geom_continuation_string='BFGS '
          cmdl%bfgs_optimisation=.true.
          cmdl%lbfgs_optimisation=.false.
          cmdl%tpsd_optimisation=.false.

       else if((method_string(1:1).eq.'d').or.(method_string(1:1).eq.'D')) then
          geom_continuation_string='DI   '
          cmdl%bfgs_optimisation=.true.
          cmdl%lbfgs_optimisation=.false.
          cmdl%tpsd_optimisation=.false.

       else if((method_string(1:1).eq.'t').or.(method_string(1:1).eq.'T')) then
          geom_continuation_string='TPSD '
          cmdl%bfgs_optimisation=.false.
          cmdl%lbfgs_optimisation=.false.
          cmdl%tpsd_optimisation=.true.

       else if((method_string(1:1).eq.'l').or.(method_string(1:1).eq.'L')) then
          geom_continuation_string='LBFGS'
          cmdl%bfgs_optimisation=.false.
          cmdl%lbfgs_optimisation=.true.
          cmdl%tpsd_optimisation=.false.
#ifdef HDF5
          call h5aopen_f(file_id,'lbfgs_num_updates',attr_id,ierr)
          call h5aread_f(attr_id,h5t_native_integer,cmdl%lbfgs_num_updates, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)
#else
          read(read_unit,iostat=ios) cmdl%lbfgs_num_updates
          call utils_read_check('geom_opt_continuation_read','mdl%lbfgs_num_updates',ios)
#endif
       else ! Must be old style model-file
          geom_continuation_string='NONE '
#ifndef HDF5
          backspace(read_unit)
#endif
          cmdl%bfgs_optimisation=.true.
          cmdl%lbfgs_optimisation=.false.
          cmdl%tpsd_optimisation=.false.
       end if

       ! Tell user what is going on
       write(stdout,'(3a)',advance ='no') &
            ' Reading continuation file "',trim(filename),'" ...'

       if(do_lbfgs.neqv.cmdl%lbfgs_optimisation) then
          call castep_model_dealloc(cmdl,do_lbfgs)
          call castep_model_alloc(cmdl,mdl%nat+mdl%nat_classical,ndim,constrain_ions, &
               cmdl%lbfgs_optimisation)
       end if

#ifdef HDF5
       if (cmdl%bfgs_optimisation) then
          call h5aopen_f(file_id,'bfgs_iteration',attr_id,ierr)
          call h5aread_f(attr_id,h5t_native_integer,cmdl%bfgs_iteration, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)
       else if(cmdl%lbfgs_optimisation) then
          call h5aopen_f(file_id,'lbfgs_iteration',attr_id,ierr)
          call h5aread_f(attr_id,h5t_native_integer,cmdl%lbfgs_iteration, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)
       else if(cmdl%tpsd_optimisation) then
          call h5aopen_f(file_id,'tpsd_iteration',attr_id,ierr)
          call h5aread_f(attr_id,h5t_native_integer,cmdl%tpsd_iteration, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)

       end if

       call h5aopen_f(file_id,'total_energy',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_double,cmdl%total_energy, &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)

       call h5aopen_f(file_id,'num_spins',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_integer,pub_num_spins, &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)

       dims=(/3,mdl%nat/)
       call h5dopen_f(file_id,'ionic_positions',dset_id,ierr)
       call h5dread_f(dset_id,h5t_native_double, &
            current_cell%ionic_positions(1:3,1,1:mdl%nat), &
            dims,ierr)
       call h5dclose_f(dset_id,ierr)

       dims=(/3,mdl%nat/)
       call h5dopen_f(file_id,'forces',dset_id,ierr)
       call h5dread_f(dset_id,h5t_native_double, &
            cmdl%forces(1:3,1,1:mdl%nat), &
            dims,ierr)
       call h5dclose_f(dset_id,ierr)
#else
       ! Read bfgs_iteration
       if(cmdl%bfgs_optimisation) then
          read(read_unit,iostat=ios) cmdl%bfgs_iteration
          call utils_read_check('geom_opt_continuation_read','cmdl%bfgs_iteration',ios)
       else if(cmdl%lbfgs_optimisation) then
          read(read_unit,iostat=ios) cmdl%lbfgs_iteration
          call utils_read_check('geom_opt_continuation_read','cmdl%lbfgs_iteration',ios)
       else if(cmdl%tpsd_optimisation) then
          read(read_unit,iostat=ios) cmdl%tpsd_iteration
          call utils_read_check('geom_opt_continuation_read','cmdl%tpsd_iteration',ios)
       end if

       ! Read total_energy
       read(read_unit,iostat=ios) cmdl%total_energy
       call utils_read_check('geom_opt_continuation_read','cmdl%total_energy',ios)

       ! Read number of spins
       read(read_unit,iostat=ios) tmp_num_spins
       call utils_read_check('geom_opt_continuation_read','pub_num_spins',ios)

       if (tmp_num_spins/=pub_num_spins) then
          call utils_abort('Error in geom_opt_continuations_read: continuation &
               &file value of num_spins does not match current value.')
       end if

       ! Read ionic positions in fractional co-ordinates into current_cell
       ! NB it is assumed that each ion has its own species
       do ionidx=1,mdl%nat
          read(read_unit,iostat=ios) current_cell%ionic_positions(1:3,1,ionidx)
          call utils_read_check('geom_opt_continuation_read','current_cell%ionic_positions',ios)
       end do

       ! Read ionic forces
       ! NB it is assumed that each ion has its own species
       do ionidx=1,mdl%nat
          read(read_unit,iostat=ios) cmdl%forces(1:3,1,ionidx)
          call utils_read_check('geom_opt_continuation_read','cmdl%forces',ios)
       end do
#endif

       ! If constrain_ions is TRUE, then cmdl%orig_forces will have been
       ! allocated (in castep_model_alloc) and we copy cmdl%forces to it
       ! NB it is assumed that each ion has its own species
       if (constrain_ions) then
          do ionidx=1,mdl%nat
             cmdl%orig_forces(1:3,1,ionidx) = cmdl%forces(1:3,1,ionidx)
          end do
       end if

#ifdef HDF5
       if(cmdl%bfgs_optimisation) then
          dims=(/ndim,ndim/)
          call h5dopen_f(file_id,'bfgs_inv_hessian',dset_id,ierr)
          call h5dread_f(dset_id,h5t_native_double,cmdl%bfgs_inv_hessian, &
               dims,ierr)
          call h5dclose_f(dset_id,ierr)
       else if(cmdl%lbfgs_optimisation) then
          dims=(/ndim,cmdl%lbfgs_num_updates/)
          call h5dopen_f(file_id,'lbfgs_position_updates',dset_id,ierr)
          call h5dread_f(dset_id,h5t_native_double,cmdl%lbfgs_position_updates, &
               dims,ierr)
          call h5dclose_f(dset_id,ierr)
          call h5dopen_f(file_id,'lbfgs_gradient_updates',dset_id,ierr)
          call h5dread_f(dset_id,h5t_native_double,cmdl%lbfgs_gradient_updates, &
               dims,ierr)
          call h5dclose_f(dset_id,ierr)
       else if(cmdl%tpsd_optimisation) then
          call h5aopen_f(file_id,'tpsd_alpha',attr_id,ierr)
          call h5aread_f(attr_id,h5t_native_double,cmdl%tpsd_alpha, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)

       end if

       call h5aopen_f(file_id,'modulus_est',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_double,pub_geom_modulus_est, &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)

       call h5aopen_f(file_id,'frequency_est',attr_id,ierr)
       call h5aread_f(attr_id,h5t_native_double,pub_geom_frequency_est, &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)
#else
       ! Read inverse Hessian
       !  read(read_unit,iostat=ios) cmdl%bfgs_inv_Hessian(1:ndim,1:ndim)


       if(cmdl%bfgs_optimisation) then
          do dim = 1,ndim
             read(read_unit,iostat=ios) cmdl%bfgs_inv_Hessian(dim,1:ndim)
          enddo
          call utils_read_check('geom_opt_continuation_read','cmdl%bfgs_inv_Hessian',ios)
       else if(cmdl%lbfgs_optimisation) then
          do dim=1,cmdl%lbfgs_num_updates
             read(read_unit,iostat=ios) cmdl%lbfgs_position_updates(1:ndim,dim)
          end do
          call utils_read_check('geom_opt_continuation_read','cmdl%lbfgs_position_updates',ios)
          do dim=1,cmdl%lbfgs_num_updates
             read(read_unit,iostat=ios) cmdl%lbfgs_gradient_updates(1:ndim,dim)
          end do
          call utils_read_check('geom_opt_continuation_read','cmdl%lbfgs_gradient_updates',ios)
       else if(cmdl%tpsd_optimisation) then
          read(read_unit,iostat=ios) cmdl%tpsd_alpha
          call utils_read_check('geom_opt_continuation_read','cmdl%tpsd_alpha',ios)
       end if

       ! Read geom_modulus_est
       read(read_unit,iostat=ios) pub_geom_modulus_est
       call utils_read_check('geom_BFGS_restore','geom_modulus_est',ios)

       ! Read geom_frequency_est
       read(read_unit,iostat=ios) pub_geom_frequency_est
       call utils_read_check('geom_BFGS_restore','geom_frequency_est',ios)
#endif

       !******************************************************************!
       ! AAM: read ngwfs from continuation file and write to              !
       !      .tightbox_ngwfs file                                        !
       !******************************************************************!

       if (pub_geom_reuse_dk_ngwfs.and.(pub_maxit_ngwf_cg>0)) then

          ! jcap: loop over regions
          do isub=1,mdl%nsub
             ! Set up name of current NGWF file to be written
             if(mdl%nsub==1)then
                write(ngwf_file,'(a,a)') trim(pub_rootname),'.tightbox_ngwfs'
             else
                write(ngwf_file,'(a,a,i0)') trim(pub_rootname),&
                     &'.tightbox_ngwfs_',isub
             end if

#ifdef HDF5
             ngwf_file=trim(ngwf_file)//'.h5'
             call h5fcreate_f(trim(ngwf_file),h5f_acc_trunc_f,file2_id,ierr)
             call h5ocopy_f(file_id,'tightbox_ngwfs',file2_id,'tightbox_ngwfs',ierr)
             call h5fclose_f(file2_id,ierr)
#else
             ! Find available unit for writing ngwfs
             write_unit = utils_unit()

             ! Open current NGWF file to be written
             open(unit=write_unit,file=ngwf_file,form='UNFORMATTED', &
                  action='WRITE',iostat=ios)
             call utils_open_unit_check('geom_opt_continuation_read',ngwf_file,ios)
             rewind(write_unit)

             ! jme: do we want to read a real or a complex file?
             read_cmplx = ( is_cmplx .and. (.not. pub_read_real_ngwfs) )

             ! jme: Copy NGWFs.
             call geom_opt_copy_ngwfs(read_unit, write_unit, &
                  mdl%regions(isub)%elements, mdl%regions(isub)%par, read_cmplx)

             ! Close current NGWF file
             close(write_unit,iostat=ios)
             call utils_close_unit_check('geom_opt_continuation_read','write_unit',ios)
          end do
#endif

       end if

       !******************************************************************!
       !================   End read/write tightbox_ngwfs  ================!
       !******************************************************************!

       !******************************************************************!
       ! AAM: read current density kernel from continuation file and      !
       !      write .dkn file                                             !
       !******************************************************************!

       if (.not.pub_edft.and.pub_geom_reuse_dk_ngwfs) then

          ! Name of current density kernel file to be written
          write(dk_file,'(2a)') trim(pub_rootname),'.dkn'

#ifdef HDF5
          dk_file=trim(dk_file)//'.h5'
          call h5fcreate_f(trim(dk_file),h5f_acc_trunc_f,file2_id,ierr)
          call h5ocopy_f(file_id,'sparse_dkn',file2_id,'sparse',ierr)
          call h5fclose_f(file2_id,ierr)
#else
          ! Find available unit for current density kernel file to be written
          write_unit = utils_unit()

          ! Open density kernel file
          open(unit=write_unit,file=dk_file,form='UNFORMATTED', &
               action='WRITE',iostat=ios)
          call utils_open_unit_check('geom_opt_continuation_read',dk_file,ios)
          rewind(write_unit)

          ! Read in density kernel file
          call utils_binary_copy(read_unit,write_unit)

          ! Close current density kernel file
          close(write_unit,iostat=ios)
          call utils_close_unit_check('geom_opt_continuation_read','write_unit',ios)
#endif

       end if

       !******************************************************************!
       !============= END read/write current density kernel ==============!
       !******************************************************************!

       !******************************************************************!
       ! ARS: read current Hamiltonian from continuation file and         !
       !      write .ham file                                             !
       !******************************************************************!

       if (pub_edft.and.pub_geom_reuse_dk_ngwfs) then

          ! Name of current hamiltonian file to be written
          write(ham_file,'(2a)') trim(pub_rootname),'.ham'

#ifdef HDF5
          ham_file=trim(ham_file)//'.h5'
          call h5fcreate_f(trim(ham_file),h5f_acc_trunc_f,file2_id,ierr)
          call h5ocopy_f(file_id,'sparse_ham',file2_id,'sparse',ierr)
          call h5fclose_f(file2_id,ierr)
#else
          ! Find available unit for current hamiltonian file to be written
          write_unit = utils_unit()

          ! Open density kernel file
          open(unit=write_unit,file=ham_file,form='UNFORMATTED', &
               action='WRITE',iostat=ios)
          call utils_open_unit_check('geom_opt_continuation_read',ham_file,ios)
          rewind(write_unit)

          ! Read in hamiltonian file
          call utils_binary_copy(read_unit, write_unit)

          ! Close current density kernel file
          close(write_unit,iostat=ios)
          call utils_close_unit_check('geom_opt_continuation_read','write_unit',ios)
#endif

       end if

       !******************************************************************!
       !============= END read/write current Hamiltonian =================!
       !******************************************************************!

#ifdef HDF5
       call h5gclose_f(file_id,ierr)
       call h5fclose_f(filesystem_id,ierr)
       call h5close_f(ierr)
#else
       ! Close read_unit
       close(read_unit,iostat=ios)
       call utils_close_unit_check('geom_opt_continuation_read','read_unit',ios)
#endif

       write(stdout,'(a)') 'done'

    end if


    call comms_barrier

    ! Broadcast data to all procs
    call comms_bcast(pub_root_proc_id,cmdl%bfgs_iteration)
    call comms_bcast(pub_root_proc_id,cmdl%lbfgs_iteration)
    call comms_bcast(pub_root_proc_id,cmdl%tpsd_iteration)
    call comms_bcast(pub_root_proc_id,cmdl%tpsd_alpha)
    call comms_bcast(pub_root_proc_id,cmdl%total_energy)
    call comms_bcast(pub_root_proc_id,cmdl%forces)
    call comms_bcast(pub_root_proc_id,geom_continuation_string)
    !call comms_bcast(pub_root_proc_id,cmdl%bfgs_inv_Hessian)
    call comms_bcast(pub_root_proc_id,current_cell%ionic_positions)
    call comms_bcast(pub_root_proc_id,pub_geom_modulus_est)
    call comms_bcast(pub_root_proc_id,pub_geom_frequency_est)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_opt_continuation_read'

    return
  end subroutine geom_opt_continuation_read

  subroutine geom_opt_continuation_write(cmdl,mdl,output_file,is_cmplx)
    !=========================================================================!
    ! Write all relevant data to file for geometry optimisation continuation  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   cmdl, intent=in, the model                                            !
    !   mdl, intent=in, contains the element data array                       !
    !   output_file, intent=in, the filename to which continuation            !
    !                           data will be written.                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !    model type                                                           !
    !    ndim                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !    comms for pub_on_root                                                !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   cmdl has already been properly read and initialised                   !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 08/02/2005                                  !
    ! Coordinate output in xyz format by Chris-Kriton Skylaris, 28/03/2008    !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    ! Modified to remove elements argument by Robert Charlton, October 2018.  !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use model_type, only: MODEL
    use rundat, only: pub_write_xyz, pub_rootname, pub_maxit_ngwf_cg, &
         pub_geom_reuse_dk_ngwfs, pub_edft, pub_num_spins
    use services, only: services_write_xyz
    use simulation_cell, only: castep_model
    use utils, only: utils_unit, utils_write_check, utils_open_unit_check, &
         utils_close_unit_check, utils_binary_copy

    implicit none

    ! Arguments
    type(castep_model), intent(in)          :: cmdl
    type(MODEL), intent(in)                 :: mdl
    character(len=file_maxpath), intent(in) :: output_file
    logical, intent(in)                     :: is_cmplx

    ! Local Variables
    integer :: output_unit,read_unit,ios,ierr,ionidx
    integer :: atom,ngwf,dim
    integer :: num,tn1,tn2,tn3,tn(3)
    integer :: isub
#ifdef ACCELRYS
    integer :: row           ! vm: required for Intel bug workaround
#endif
    real(kind=DP) :: orig1,orig2,orig3
    real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox
    character(len=file_maxpath) :: read_file
    character(len=file_maxpath) :: xyz_title_line ! title line of coords in xyz file
#ifdef HDF5
    integer(hid_t) :: filesystem_id
    integer(hid_t) :: file_id, file2_id
    integer(hid_t) :: group_id
    integer(hid_t) :: attr_id
    integer(hid_t) :: aspace_id1
    integer(hid_t) :: dspace_id
    integer(hid_t) :: dset_id
    integer(hsize_t) :: adims(1), dims(2)
#endif

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Entering geom_opt_continuation_write'

    if (pub_on_root) then

       ! cks: -------- XYZ OUTPUT ----------------------------------------
       if (pub_write_xyz) then
          if(do_lbfgs) then
             write(xyz_title_line,'(a, i5)')"Geometry optimisation iteration: ", cmdl%lbfgs_iteration
          else if(do_bfgs) then
             write(xyz_title_line,'(a, i5)')"Geometry optimisation iteration: ", cmdl%bfgs_iteration
          else if(do_tpsd) then
             write(xyz_title_line,'(a, i5)')"Geometry optimisation iteration: ", cmdl%tpsd_iteration
          end if
          call services_write_xyz(mdl%elements, pub_rootname, trim(xyz_title_line) )
       endif
       ! cks: ---- END XYZ OUTPUT ----------------------------------------



#ifdef HDF5
       call h5open_f(ierr)
       call h5fcreate_f(trim(output_file),h5f_acc_trunc_f,filesystem_id,ierr)
       call h5gcreate_f(filesystem_id,'continuation',file_id,ierr)

       adims=(/1/)
       call h5screate_simple_f(1,adims,aspace_id1,ierr)
       call h5acreate_f(file_id,'geom_string',h5t_native_character, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_character,geom_string, &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)
       if (do_lbfgs) then
          call h5acreate_f(file_id,'lbfgs_num_updates',h5t_native_integer, &
               aspace_id1,attr_id,ierr)
          call h5awrite_f(attr_id,h5t_native_integer,cmdl%lbfgs_num_updates, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)
       end if
#else
       ! Find available unit specifier for writing continuation file
       output_unit = utils_unit()

       ! Open output_file
       open(unit=output_unit,file=output_file,iostat=ios,&
            form='UNFORMATTED',action='WRITE')
       call utils_open_unit_check('geom_opt_continuation_write',output_file,ios)
       rewind(output_unit)

       if(do_bfgs.or.do_tpsd) then
          write(output_unit) geom_string!'BFGS' or 'TPSD'
       else if(do_lbfgs) then
          write(output_unit) geom_string!'LBFGS'
          write(output_unit) cmdl%lbfgs_num_updates
       end if
#endif


       ! Tell user what is going on
       write(stdout,'(/3a)',advance='no') ' Writing continuation file "', &
            trim(output_file),'" ...'

#ifdef HDF5
       if (do_bfgs) then
          call h5acreate_f(file_id,'bfgs_iteration',h5t_native_integer, &
               aspace_id1,attr_id,ierr)
          call h5awrite_f(attr_id,h5t_native_integer,cmdl%bfgs_iteration, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)
       else if(do_tpsd) then
          call h5acreate_f(file_id,'tpsd_iteration',h5t_native_integer, &
               aspace_id1,attr_id,ierr)
          call h5awrite_f(attr_id,h5t_native_integer,cmdl%tpsd_iteration, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)
       else
          call h5acreate_f(file_id,'lbfgs_iteration',h5t_native_integer, &
               aspace_id1,attr_id,ierr)
          call h5awrite_f(attr_id,h5t_native_integer,cmdl%lbfgs_iteration, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)
       end if

       call h5acreate_f(file_id,'total_energy',h5t_native_double, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_double,cmdl%total_energy, &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)

       call h5acreate_f(file_id,'num_spins',h5t_native_integer, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_integer,pub_num_spins, &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)

       dims=(/3,cmdl%cell%num_ions - mdl%nat_classical/)
       call h5screate_simple_f(2,dims,dspace_id,ierr)
       call h5dcreate_f(file_id,'ionic_positions',h5t_native_double, &
            dspace_id,dset_id,ierr)
       call h5dwrite_f(dset_id,h5t_native_double, &
            cmdl%cell%ionic_positions(1:3,1,1:cmdl%cell%num_ions-mdl%nat_classical), &
            dims,ierr)
       call h5dclose_f(dset_id,ierr)

       dims=(/3,mdl%nat/)
       call h5screate_simple_f(2,dims,dspace_id,ierr)
       call h5dcreate_f(file_id,'forces',h5t_native_double, &
            dspace_id,dset_id,ierr)
       if (constrain_ions) then
          call h5dwrite_f(dset_id,h5t_native_double, &
               cmdl%orig_forces(1:3,1,1:mdl%nat), &
               dims,ierr)
       else
          call h5dwrite_f(dset_id,h5t_native_double, &
               cmdl%forces(1:3,1,1:mdl%nat), &
               dims,ierr)
       end if
       call h5dclose_f(dset_id,ierr)

       if(do_bfgs) then
          dims=(/ndim,ndim/)
          call h5screate_simple_f(2,dims,dspace_id,ierr)
          call h5dcreate_f(file_id,'bfgs_inv_hessian',h5t_native_double, &
               dspace_id,dset_id,ierr)
          call h5dwrite_f(dset_id,h5t_native_double,cmdl%bfgs_inv_hessian, &
               dims,ierr)
          call h5dclose_f(dset_id,ierr)
       else if(do_lbfgs) then
          dims=(/ndim,cmdl%lbfgs_num_updates/)
          call h5screate_simple_f(2,dims,dspace_id,ierr)
          call h5dcreate_f(file_id,'lbfgs_position_updates',h5t_native_double, &
               dspace_id,dset_id,ierr)
          call h5dwrite_f(dset_id,h5t_native_double,cmdl%lbfgs_position_updates, &
               dims,ierr)
          call h5dclose_f(dset_id,ierr)
          call h5dcreate_f(file_id,'lbfgs_gradient_updates',h5t_native_double, &
               dspace_id,dset_id,ierr)
          call h5dwrite_f(dset_id,h5t_native_double,cmdl%lbfgs_gradient_updates, &
               dims,ierr)
          call h5dclose_f(dset_id,ierr)
       else if(do_tpsd) then
          call h5acreate_f(file_id,'tpsd_alpha',h5t_native_double, &
               aspace_id1,attr_id,ierr)
          call h5awrite_f(attr_id,h5t_native_double,cmdl%tpsd_alpha, &
               adims,ierr)
          call h5aclose_f(attr_id,ierr)

       end if

       call h5acreate_f(file_id,'modulus_est',h5t_native_double, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_double,new_modulus_est, &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)

       call h5acreate_f(file_id,'frequency_est',h5t_native_double, &
            aspace_id1,attr_id,ierr)
       call h5awrite_f(attr_id,h5t_native_double,new_frequency_est, &
            adims,ierr)
       call h5aclose_f(attr_id,ierr)
#else
       ! Write bfgs_iteration
       if(do_bfgs) then
          write(output_unit,iostat=ios) cmdl%bfgs_iteration
          call utils_write_check('geom_opt_continuation_write','cmdl%bfgs_iteration',ios)
       else if(do_lbfgs) then
          write(output_unit,iostat=ios) cmdl%lbfgs_iteration
          call utils_write_check('geom_opt_continuation_write','cmdl%lbfgs_iteration',ios)
       else if(do_tpsd) then
          write(output_unit,iostat=ios) cmdl%tpsd_iteration
          call utils_write_check('geom_opt_continuation_write','cmdl%tpsd_iteration',ios)
       end if

       ! Write total_energy
       write(output_unit,iostat=ios) cmdl%total_energy
       call utils_write_check('geom_opt_continuation_write','cmdl%total_energy',ios)

       ! Write number of spins
       write(output_unit,iostat=ios) pub_num_spins
       call utils_write_check('geom_opt_continuation_write','pub_num_spins',ios)

       ! Write ionic positions in fractional co-ordinates
       ! NB it is assumed that each ion has its own species
       do ionidx=1,cmdl%cell%num_ions - mdl%nat_classical
          write(output_unit,iostat=ios) cmdl%cell%ionic_positions(1:3,1,ionidx)
          call utils_write_check('geom_opt_continuation_write','cmdl%cell%ionic_positions',ios)
       end do

       ! Write the UNCONSTRAINED ionic forces so that on continuation different
       ! ionic constraints may be used
       ! NB it is assumed that each ion has its own species
       !do ionidx=1,cmdl%cell%num_ions
       ! kaw: We only need forces for the QM ions
       do ionidx=1,mdl%nat
          if (constrain_ions) then
             write(output_unit,iostat=ios) cmdl%orig_forces(1:3,1,ionidx)
             call utils_write_check('geom_opt_continuation_write','cmdl%orig_forces',ios)
          else
             write(output_unit,iostat=ios) cmdl%forces(1:3,1,ionidx)
             call utils_write_check('geom_opt_continuation_write','cmdl%forces',ios)
          end if
       end do

       ! Write inverse Hessian
       ! write(output_unit,iostat=ios) cmdl%bfgs_inv_Hessian(1:ndim,1:ndim)

       if(do_bfgs) then
          do dim=1,ndim
             write(output_unit,iostat=ios) cmdl%bfgs_inv_Hessian(dim,1:ndim)
          end do
          call utils_write_check('geom_opt_continuation_write','cmdl%bfgs_inv_Hessian',ios)
       else if(do_lbfgs) then
          do dim=1,cmdl%lbfgs_num_updates
             write(output_unit,iostat=ios) cmdl%lbfgs_position_updates(1:ndim,dim)
          end do
          call utils_write_check('geom_opt_continuation_write','cmdl%lbfgs_position_updates',ios)
          do dim=1,cmdl%lbfgs_num_updates
             write(output_unit,iostat=ios) cmdl%lbfgs_gradient_updates(1:ndim,dim)
          end do
          call utils_write_check('geom_opt_continuation_write','cmdl%lbfgs_gradient_updates',ios)
       else if(do_tpsd) then
          write(output_unit,iostat=ios) cmdl%tpsd_alpha
          call utils_write_check('geom_opt_continuation_write','cmdl%tpsd_alpha',ios)
       end if

       ! Write the latest estimate of the bulk modulus
       write(output_unit,iostat=ios) new_modulus_est
       call utils_write_check('geom_opt_continuation_write','new_modulus_est',ios)

       ! Write the latest estimate of the phonon frequency
       write(output_unit,iostat=ios) new_frequency_est
       call utils_write_check('geom_opt_continuation_write','new_frequency_est',ios)
#endif

       !******************************************************************!
       ! AAM: read current tightbox_ngwfs from file and write to          !
       !      continuation file                                           !
       !******************************************************************!

       if (pub_geom_reuse_dk_ngwfs.and.(pub_maxit_ngwf_cg>0)) then

          ! jcap: loop over regions
          do isub=1,mdl%nsub
             ! Set up name of current NGWF file to be read
             if(mdl%nsub==1)then
                write(read_file,'(a,a)') trim(pub_rootname),'.tightbox_ngwfs'
             else
                write(read_file,'(a,a,i0)') trim(pub_rootname),&
                     &'.tightbox_ngwfs_',isub
             end if

#ifdef HDF5
             read_file=trim(read_file)//'.h5'
             call h5fopen_f(trim(read_file),h5f_acc_rdonly_f,file2_id,ierr)
             call h5ocopy_f(file2_id,'tightbox_ngwfs',file_id,'tightbox_ngwfs',ierr)
             call h5fclose_f(file2_id,ierr)
#else
             ! Find available unit for reading current NGWFs
             read_unit = utils_unit()

             ! Open current NGWF file to be read
             open(unit=read_unit,file=read_file,form='UNFORMATTED',status='OLD', &
                  action='READ',iostat=ios)
             call utils_open_unit_check('geom_opt_continuation_write',read_file,ios)
             rewind(read_unit)

             ! jme: Copy NGWFs.
             call geom_opt_copy_ngwfs(read_unit, output_unit, &
                  mdl%regions(isub)%elements, mdl%regions(isub)%par, is_cmplx)

             ! Close current NGWF file
             close(read_unit,iostat=ios)
             call utils_close_unit_check('geom_opt_continuation_write','read_unit',ios)
#endif
          end do

       end if

       !******************************************************************!
       !================   End read/write tightbox_ngwfs  ================!
       !******************************************************************!

       !******************************************************************!
       ! AAM: read current density kernel from file and write to          !
       !      continuation file                                           !
       !******************************************************************!

       if (.not.pub_edft.and.pub_geom_reuse_dk_ngwfs) then

          ! Name of current density kernel file to be read
          write(read_file,'(2a)') trim(pub_rootname),'.dkn'

#ifdef HDF5
          read_file=trim(read_file)//'.h5'
          call h5fopen_f(trim(read_file),h5f_acc_rdonly_f,file2_id,ierr)
          call h5ocopy_f(file2_id,'sparse',file_id,'sparse_dkn',ierr)
          call h5fclose_f(file2_id,ierr)
#else
          ! Find available unit for current density kernel file to be read
          read_unit = utils_unit()

          ! Open current density kernel file
          open(unit=read_unit,file=read_file,form='UNFORMATTED', &
               status='OLD',action='READ',iostat=ios)
          call utils_open_unit_check('geom_opt_continuation_write',read_file,ios)
          rewind(read_unit)

          ! Read in density kernel file
          call utils_binary_copy(read_unit,output_unit)

          ! Close current density kernel file
          close(read_unit,iostat=ios)
          call utils_close_unit_check('geom_opt_continuation_write','read_unit',ios)
#endif

       end if

       !******************************************************************!
       !============= END read/write current density kernel ==============!
       !******************************************************************!

       !******************************************************************!
       ! ARS: read current Hamiltonian from file and write to             !
       !      continuation file                                           !
       !******************************************************************!

       if (pub_edft.and.pub_geom_reuse_dk_ngwfs) then

          ! Name of current hamiltonian file to be read
          write(read_file,'(2a)') trim(pub_rootname),'.ham'

#ifdef HDF5
          read_file=trim(read_file)//'.h5'
          call h5fopen_f(trim(read_file),h5f_acc_rdonly_f,file2_id,ierr)
          call h5ocopy_f(file2_id,'sparse',file_id,'sparse_ham',ierr)
          call h5fclose_f(file2_id,ierr)
#else
          ! Find available unit for current hamiltonian file to be read
          read_unit = utils_unit()

          ! Open current hamiltonian file
          open(unit=read_unit,file=read_file,form='UNFORMATTED', &
               status='OLD',action='READ',iostat=ios)
          call utils_open_unit_check('geom_opt_continuation_write',read_file,ios)
          rewind(read_unit)

          ! Read in hamiltonian file
          call utils_binary_copy(read_unit,output_unit)

          ! Close current hamiltonian file
          close(read_unit,iostat=ios)
          call utils_close_unit_check('geom_opt_continuation_write','read_unit',ios)
#endif

       end if

       !******************************************************************!
       !============= END read/write current Hamiltonian =================!
       !******************************************************************!

#ifdef HDF5
       call h5gclose_f(file_id,ierr)
       call h5fclose_f(filesystem_id,ierr)
       call h5close_f(ierr)
#else
       ! Close output_unit
       close(output_unit,iostat=ios)
       call utils_close_unit_check('geom_opt_continuation_write','output_unit',ios)
#endif

       write(stdout,'(a)') 'done'

    end if

    if (pub_debug_on_root) &
         write(stdout,'(a)') 'DEBUG: Leaving geom_opt_continuation_write'

    return
  end subroutine geom_opt_continuation_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef HDF5

  subroutine geom_opt_copy_ngwfs(read_unit, write_unit, elements, par, is_cmplx)
    !=========================================================================!
    ! Copy NGWFs from one file unit to another (non-HDF5 compilations).       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   read_unit:  source      unit id                                       !
    !   write_unit: destination unit id                                       !
    !   elements:   the element data array                                    !
    !   is_cmplx:   type of NGWF data to be copies                            !
    !-------------------------------------------------------------------------!
    ! Initial version extracted by Jose M Escartin from similar code in       !
    ! geom_opt_continuation_read & geom_opt_continuation_write, 30 Jan 2017.  !
    ! Modified by Jose M Escartin in February 2017 to incorporate new format  !
    ! for tightbox files.                                                     !
    ! Modified to remove pub_par by Joseph Prentice, September 2018           !
    !-------------------------------------------------------------------------!

    use constants, only: stdout
    use datatypes, only: FFTBOX_DATA, data_fftbox_alloc, data_fftbox_dealloc
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_devel_code
    use utils, only: utils_read_check, utils_write_check, utils_assert, &
         utils_alloc_check, utils_dealloc_check, utils_devel_code

    implicit none

    ! Arguments
    integer, intent(in) :: read_unit, write_unit
    type(PARAL_INFO), intent(in) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)
    logical, intent(in) :: is_cmplx

    ! Local variables
    integer :: num, tn1, tn2, tn3, tn(3)
    integer :: dim, atom, ngwf, ierr, ios
    real(kind=DP) :: orig1, orig2, orig3
    type(FFTBOX_DATA) :: uni_tbox
    logical :: file_iscmplx, use_old_tightbox_format
    integer :: tightbox_version
    character(len=len('VERSION_')) :: tag

    ! jme: The new tightbox format is enforced for complex tightboxes.
    ! For real tightboxes we also use the new format, unless the appropriate
    ! devel code is used.
    ios = 0
    use_old_tightbox_format = utils_devel_code(.false., 'RESTART', &
         'OLD_TIGHTBOX_FORMAT', pub_devel_code, no_bcast=.true.)
    if (is_cmplx.or.(.not.use_old_tightbox_format)) then
       if (pub_debug_on_root) then
          write(stdout,*) 'Attempting to read new tightbox format file.'
       end if
       read(read_unit, iostat=ios) tag

       if ((ios==0) .and. (tag == 'VERSION_')) then
          ! This file follows the new format.
          read(read_unit, iostat=ios) tightbox_version
          call utils_read_check('geom_opt_copy_ngwfs','tightbox_version',ios)

          if (pub_debug_on_root) then
             write(stdout,'(a,f4.2,a)') &
                  ' ONETEP tightbox format detected v', &
                  real(tightbox_version,DP)/100.0_DP, ' .'
          end if

          ! Write version headers to file.
          write(write_unit,iostat=ios) 'VERSION_'
          call utils_write_check('geom_opt_copy_ngwfs', 'VERSION_', ios)
          write(write_unit,iostat=ios) tightbox_version
          call utils_write_check('geom_opt_copy_ngwfs', 'tightbox_version', ios)

          ! Read and write complex flag.
          read(read_unit, iostat=ios) file_iscmplx
          call utils_read_check('geom_opt_copy_ngwfs','file_iscmplx',ios)
          write(write_unit,iostat=ios) file_iscmplx
          call utils_write_check('geom_opt_copy_ngwfs','file_iscmplx', ios)
       else
          write(*,*) 'Tightbox format version > v1.00 expected but not &
               &recognised.'
          ios = -1
          backspace read_unit
       end if
    end if
    if ( ((.not.is_cmplx).and.use_old_tightbox_format) .or. (ios.ne.0)) then
       ! File should follow the old format.
       write(stdout,*) 'Attempting to read old tightbox format file.'
       ! jme: Complex NGWFs were consolidated at the same time that the new
       ! tightbox file format was introduced (and enforced for them), so
       ! nobody should be using files with complex tightboxes
       ! in the old format.
       file_iscmplx = .false.
    end if

    ! jme: check that the file actually contains the type of data
    ! (real/complex) that we want to read
    call utils_assert(is_cmplx .eqv. file_iscmplx, 'Error in &
         &geom_opt_copy_ngwfs: complex character of data &
         &expected to be read differs from that detected in the &
         &source file. Their respective complex attributes are: ', &
         is_cmplx, file_iscmplx)

    ! Read number of NGWFs
    read(read_unit,iostat=ios) num
    call utils_read_check('geom_opt_copy_ngwfs','num',ios)

    ! Sanity check
    call utils_assert(num == par%num_ngwfs, &
         'Error in geom_opt_copy_ngwfs: ''num'' not equal to &
         &''par%num_ngwfs''. The two values were ',num, &
         par%num_ngwfs)

    ! Write number of NGWFs to file
    write(write_unit,iostat=ios) num
    call utils_write_check('geom_opt_copy_ngwfs','num',ios)

    ! Loop over dimensions
    do dim=1,3

       ! Read number of points in dimension dim of universal tightbox
       read(read_unit,iostat=ios) tn(dim)
       call utils_read_check('geom_opt_copy_ngwfs','tn',ios)

       ! Write number of points in dimension dim of universal tightbox
       write(write_unit,iostat=ios) tn(dim)
       call utils_write_check('geom_opt_copy_ngwfs','tn',ios)

    end do

    ! Allocate universal tightbox
    tn1 = tn(1) ; tn2 = tn(2) ; tn3 = tn(3)
    call data_fftbox_alloc(uni_tbox, tn1, tn2, tn3, iscmplx=file_iscmplx)

    atom_loop: do atom=1,par%nat

       ngwfs_on_atom_loop: do ngwf=1,elements(atom)%nfunctions

          ! Read origin
          read(read_unit,iostat=ios) orig1, orig2, orig3
          call utils_read_check('geom_opt_copy_ngwfs','orig1, orig2, orig3',ios)

          ! Write origin
          write(write_unit,iostat=ios) orig1, orig2, orig3
          call utils_write_check('geom_opt_copy_ngwfs','orig1, orig2, orig3',ios)

          if (file_iscmplx) then
             ! Read tightbox
             read(read_unit,iostat=ios) uni_tbox%z(1:tn1,1:tn2,1:tn3)
             call utils_read_check('geom_opt_copy_ngwfs','uni_tbox',ios)

             ! Write tightbox
             write(write_unit,iostat=ios) uni_tbox%z(1:tn1,1:tn2,1:tn3)
             call utils_write_check('geom_opt_copy_ngwfs','uni_tbox',ios)
          else
             ! Read tightbox
             read(read_unit,iostat=ios) uni_tbox%d(1:tn1,1:tn2,1:tn3)
             call utils_read_check('geom_opt_copy_ngwfs','uni_tbox',ios)

             ! Write tightbox
             write(write_unit,iostat=ios) uni_tbox%d(1:tn1,1:tn2,1:tn3)
             call utils_write_check('geom_opt_copy_ngwfs','uni_tbox',ios)
          end if

       end do ngwfs_on_atom_loop

    end do atom_loop

    ! Deallocate memory
    call data_fftbox_dealloc(uni_tbox)

  end subroutine geom_opt_copy_ngwfs

  subroutine geom_opt_copy_file(routine, src_file, dest_file, elements, par, &
       is_cmplx)
    !=========================================================================!
    ! Copy files (create backup or restore backup).                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   routine:   calling procedure name                                     !
    !   src_file:  source      file  name                                     !
    !   dest_file: destination file  name                                     !
    !   elements:  the element data array         (optional, reqd. for NGWFs) !
    !   is_cmplx:  type of NGWF data to be copies (optional, reqd. for NGWFs) !
    !-------------------------------------------------------------------------!
    ! Written by Jose M Escartin, 14/03/2017.                                 !
    ! Modified to remove pub_par by Joseph Prentice, September 2018           !
    !-------------------------------------------------------------------------!

    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use utils, only: utils_binary_copy, utils_unit, utils_open_unit_check, &
         utils_close_unit_check, utils_assert
    implicit none

    ! Arguments
    character(len=*), intent(in)           :: routine
    character(len=*), intent(in)           :: src_file
    character(len=*), intent(in)           :: dest_file
    type(PARAL_INFO), intent(in), optional :: par
    type(ELEMENT),    intent(in), optional :: elements(:)
    logical,          intent(in), optional :: is_cmplx

    ! Local variables
    integer :: read_unit, write_unit, ios
    logical :: is_ngwf_copy
    character(len=*), parameter :: myself = 'geom_opt_copy_file'
    character(len=:), allocatable :: reference

    ! Combination of name of subroutine and name of calling routine.
    reference = myself // '(' // trim(routine) // ')'

    ! Check that either two or zero optional arguments are present.
    call utils_assert(present(elements).eqv.present(is_cmplx), 'Error in ' &
         // myself // ' when called from ' // trim(routine) //': inconsistent &
         &presence of optional arguments "elements" and "is_cmplx":', &
         present(elements), present(is_cmplx))
    is_ngwf_copy =  present(is_cmplx)

    ! Open source file.
    read_unit = utils_unit()
    open(unit=read_unit, file=src_file, form='UNFORMATTED', status='OLD', &
         action='READ', iostat=ios)
    call utils_open_unit_check(reference, trim(src_file), ios)
    rewind(read_unit)

    ! Open destination file.
    write_unit = utils_unit()
    open(unit=write_unit, file=dest_file, form='UNFORMATTED', &
         status='REPLACE', action='WRITE', iostat=ios)
    call utils_open_unit_check(reference, trim(dest_file), ios)
    rewind(write_unit)

    ! Copy data from one file to the other.
    if (is_ngwf_copy) then
       call geom_opt_copy_ngwfs(read_unit, write_unit, elements, par, is_cmplx)
    else
       call utils_binary_copy(read_unit, write_unit)
    end if

    ! Close both files.
    close(read_unit, iostat=ios)
    call utils_close_unit_check(reference, trim(src_file), ios)
    close(write_unit, iostat=ios)
    call utils_close_unit_check(reference, trim(dest_file), ios)

  end subroutine geom_opt_copy_file

#else

  subroutine geom_opt_copy_file_h5(src_file, dest_file, obj)
    !=========================================================================!
    ! Copy HDF5 files with a single object.                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   src_file:  source      file       name                                !
    !   dest_file: destination file       name                                !
    !   obj:                   object     name                                !
    !-------------------------------------------------------------------------!
    ! Written by Jose M Escartin, 14/03/2017.                                 !
    !-------------------------------------------------------------------------!

    implicit none

    ! Arguments
    character(len=*), intent(in)           :: src_file
    character(len=*), intent(in)           :: dest_file
    character(len=*), intent(in)           :: obj

    ! Local variables
    integer :: ierr
    integer(hid_t) :: src_file_id, dest_file_id

    ! Open src file (read-only mode).
    call h5fopen_f(trim(src_file), h5f_acc_rdonly_f, src_file_id, ierr)

    ! Create dest file (rewrite if it already exists).
    call h5fcreate_f(trim(dest_file), h5f_acc_trunc_f, dest_file_id, ierr)

    ! Copy content of src file into dest file.
    call h5ocopy_f(src_file_id, obj, dest_file_id, obj, ierr)

    ! Close files
    call h5fclose_f(src_file_id, ierr)
    call h5fclose_f(dest_file_id, ierr)

  end subroutine geom_opt_copy_file_h5

#endif

  subroutine geom_precond_initialize(cmdl,mdl,is_cmplx)

    use simulation_cell, only: castep_model
    use model_type, only: MODEL
    use rundat, only: pub_geom_precond_type, pub_geom_precond_scale_cell
    use bibliography, only: bibliography_cite
    use utils, only: utils_abort

    type(castep_model), intent(inout) :: cmdl
    type(MODEL),        intent(inout) :: mdl
    ! ja531 -> added is_cmplx for geom_get_forces
    ! agrecocmplx
    logical, intent(in)               :: is_cmplx


    select case(pub_geom_precond_type)
    case('NONE')
       return
    case('ID')
       call geom_precond_id_initialize(cmdl,mdl)
    case('EXP')
       call geom_precond_exp_initialize(cmdl,mdl,is_cmplx)
       call bibliography_cite('Precond-EXP')
    case('FF')
       call geom_precond_ff_initialize(cmdl,mdl)
       call bibliography_cite('Precond-FF')
       call bibliography_cite('Lindh-FF')
    case default
       call utils_abort('Error in geom_precond_initialize - unrecognised preconditioner request='// &
               & trim(pub_geom_precond_type)//'.')
    end select

    precond%scale_cell = pub_geom_precond_scale_cell

  end subroutine geom_precond_initialize

  subroutine geom_precond_id_initialize(cmdl,mdl)

    use simulation_cell, only: castep_model
    use model_type, only: MODEL
    use comms, only: pub_on_root
    use utils, only: utils_abort, utils_alloc_check

    type(castep_model), intent(inout) :: cmdl
    type(MODEL),        intent(inout) :: mdl

    integer                          :: ierr, i

    if (.not. precond%initialized) then
       if (pub_on_root) then
          precond%initialized = .true.

          if (constrain_ions) then
             allocate(precond%ion_constraint_type(mdl%nat),stat=ierr)
             call utils_alloc_check('geom_precond_id_initialize','precond%ion_constraint_type',ierr)
             do i=1,mdl%nat
                precond%ion_constraint_type(i) = mdl%elements(i)%ion_constraint_type
                select case (precond%ion_constraint_type(i))
                case ('NONE')
                case ('FIXED')
                case default
                   call utils_abort('Error in geom_precond_id_initialize: &
                                   & illegal value for precond%ion_constraint_type(i)=', i)
                end select
             end do
          end if

          allocate(precond%P(3*cmdl%cell%num_ions, 3*cmdl%cell%num_ions),stat=ierr)
          call utils_alloc_check('geom_precond_id_initialize','precond%P',ierr)
          precond%P=0.0_dp

       end if
    end if

  end subroutine geom_precond_id_initialize

  subroutine geom_precond_exp_initialize(cmdl,mdl,is_cmplx)

    use simulation_cell, only: castep_model
    use model_type, only: MODEL
    use constants, only: stdout
    use comms, only: pub_on_root
    use utils, only: utils_alloc_check, utils_abort
    use rundat, only: pub_geom_precond_exp_c_stab, pub_geom_precond_exp_A, &
                    & pub_geom_precond_exp_r_NN, pub_geom_precond_exp_r_cut, pub_geom_precond_exp_mu

    type(castep_model), intent(inout) :: cmdl
    type(MODEL),        intent(inout) :: mdl
    ! ja531 -> added is_cmplx for geom_get_forces
    ! agrecocmplx
    logical, intent(in)               :: is_cmplx

    integer :: ierr, i

    if (.not. precond%initialized) then
       if (pub_on_root) then
          precond%initialized = .true.

          if (constrain_ions) then
             allocate(precond%ion_constraint_type(mdl%nat),stat=ierr)
             call utils_alloc_check('geom_precond_exp_initialize','precond%ion_constraint_type',ierr)
             do i=1,mdl%nat
                precond%ion_constraint_type(i) = mdl%elements(i)%ion_constraint_type
                select case (precond%ion_constraint_type(i))
                case ('NONE')
                case ('FIXED')
                case default
                   call utils_abort('Error in geom_precond_exp_initialize: &
                                   & illegal value for precond%ion_constraint_type(i)=', i)
                end select
             end do
          end if

          allocate(precond%P(3*cmdl%cell%num_ions, 3*cmdl%cell%num_ions),stat=ierr)
          call utils_alloc_check('geom_precond_exp_initialize','precond%P',ierr)
          precond%P=0.0_dp

          allocate(precond%rmatrix(3, cmdl%cell%num_ions),stat=ierr)
          call utils_alloc_check('geom_precond_exp_initialize','precond%rmatrix',ierr)
          precond%rmatrix=0.0_dp

          allocate(precond%dmatrix(cmdl%cell%num_ions, cmdl%cell%num_ions),stat=ierr)
          call utils_alloc_check('geom_precond_exp_initialize','precond%dmatrix',ierr)
          precond%dmatrix=0.0_dp

          precond%exp_c_stab = pub_geom_precond_exp_c_stab
          if (precond%exp_c_stab <= 0.0_dp) then
             call utils_abort('pub_geom_precond_exp_c_stab must be larger than 0!')
          end if

          precond%exp_A = pub_geom_precond_exp_A

          precond%exp_r_NN = pub_geom_precond_exp_r_NN
          if (precond%exp_r_NN < 0.0_dp) then
             call utils_abort('pub_geom_precond_exp_r_NN must be non-negatve!')
          end if
          if (precond%exp_r_NN == 0.0_dp ) then
             precond%exp_r_NN = geom_precond_exp_get_r_NN(cmdl)
             write(stdout, 1) 'EXP precond LBFGS => estimated nearest neighbor distance', &
                              & precond%exp_r_NN, trim(length_label)
          end if

          precond%exp_r_cut = pub_geom_precond_exp_r_cut
          if (precond%exp_r_cut < 0.0_dp) then
             call utils_abort('pub_geom_precond_exp_r_cut must be non-negatve!')
          end if
          if (precond%exp_r_cut == 0.0_dp ) then
             precond%exp_r_cut = 2.0 * precond%exp_r_NN
             write(stdout, 1) 'EXP precond LBFGS => estimated cutoff distance', &
                              & precond%exp_r_cut, trim(length_label)
          elseif (precond%exp_r_cut < precond%exp_r_NN ) then
             precond%exp_r_cut = precond%exp_r_NN
          end if
       end if

       precond%exp_mu = pub_geom_precond_exp_mu
       if (precond%exp_mu < 0.0_dp) then
          call utils_abort('pub_geom_precond_exp_mu must be non-negatve!')
       end if
       if (precond%exp_mu == 0.0_dp ) then
           precond%exp_mu = geom_precond_exp_estimate_mu(cmdl,mdl,is_cmplx)
           if (pub_on_root) write(stdout, 1) 'EXP precond LBFGS => estimated parameter mu', &
                                              & precond%exp_mu, trim(force_constant_label)
       end if
    end if

1   format(1x,A,' : ',G10.4,3x,A)

  end subroutine geom_precond_exp_initialize

  function geom_precond_exp_get_r_NN(cmdl) result(r_NN)

    use simulation_cell, only: castep_model, castep_cell_cart_lattice_to_abc
    use utils, only: utils_abort

    type(castep_model), intent(in) :: cmdl
    real (kind=dp)                :: r_NN

    real (kind=dp) :: r_cut, phi, extent
    real (kind=dp) :: a, b, c, alpha, beta, gamma
    real (kind=dp) :: d_NN
    integer :: i, j

    call  geom_precond_get_dmatrix(cmdl)

    r_cut = 1.0_dp                          ! starting with 1.0 Bohr
    phi = (1.0_dp + sqrt(5.0_dp)) / 2.0_dp  ! golden ratio

    call castep_cell_cart_lattice_to_abc(cmdl%cell,a,b,c,alpha,beta,gamma)
    extent = max(a, b, c)

    r_NN = 0.0_dp
    do i=1, cmdl%cell%num_ions
       do while (.true.)
          if (r_cut > 2.0 * extent) then
             call utils_abort('Error geom_precond_exp_get_r_NN: increased r_cut &
                             & to twice system extent without finding neighbours &
                             & for all atoms. This can happen if your system is too small; &
                             & try setting r_cut manually')
          end if
          d_NN = 0.0_dp
          ! find nearest neighbor of ith atom
          do j=1, cmdl%cell%num_ions
             if (i==j) cycle
             if (precond%dmatrix(i,j) <= r_cut) then
                if ( (d_NN == 0.0_dp) .or. (precond%dmatrix(i,j) < d_NN) ) then
                   d_NN = precond%dmatrix(i,j)
                end if
             end if
          end do
          ! maximum of nearest neigbour distances
          if (d_NN > 0.0_dp) then
             if (d_NN > r_NN) r_NN = d_NN
             exit
          else
             r_cut = r_cut*phi
          end if
       end do
    end do

  end function geom_precond_exp_get_r_NN

  function geom_precond_exp_estimate_mu(cmdl,mdl,is_cmplx) result(mu)

    use simulation_cell, only: castep_model, castep_cell_frac_to_cart, castep_cell_cart_to_frac
    use model_type, only: MODEL
    use comms, only: pub_on_root, pub_root_proc_id, comms_bcast
    use utils, only: utils_alloc_check, utils_dealloc_check

    type(castep_model), intent(inout) :: cmdl
    type(MODEL),        intent(inout) :: mdl
    ! ja531 -> added is_cmplx for geom_get_forces
    ! agrecocmplx
    logical, intent(in)               :: is_cmplx
    real (kind=dp) :: mu

    integer :: ierr
    integer :: i, k
    integer :: ispec, iatom
    real (kind=dp), allocatable :: x_vec(:)
    real (kind=dp), allocatable :: xcart(:,:), vcart(:,:), g0cart(:,:), g1cart(:,:)
    real (kind=dp) :: L(3), M(3,3)
    logical :: constrain_ions_orig

    constrain_ions_orig = constrain_ions
    constrain_ions = .false.
    if (constrain_ions_orig) call geom_get_forces(cmdl,mdl,is_cmplx=is_cmplx)

    allocate(x_vec(ndim),stat=ierr)
    call utils_alloc_check('geom_precond_exp_estimate_mu','x_vec',ierr)

    if (pub_on_root) then

       allocate(xcart(3, cmdl%cell%num_ions),stat=ierr)
       call utils_alloc_check('geom_precond_exp_estimate_mu','xcart',ierr)

       allocate(vcart(3, cmdl%cell%num_ions),stat=ierr)
       call utils_alloc_check('geom_precond_exp_estimate_mu','vcart',ierr)

       allocate(g0cart(3, cmdl%cell%num_ions),stat=ierr)
       call utils_alloc_check('geom_precond_exp_estimate_mu','g0cart',ierr)

       allocate(g1cart(3, cmdl%cell%num_ions),stat=ierr)
       call utils_alloc_check('geom_precond_exp_estimate_mu','g1cart',ierr)

       i = 1
       do ispec=1,cmdl%cell%num_species
          do iatom=1,cmdl%cell%num_ions_in_species(ispec)
             call castep_cell_frac_to_cart(cmdl%cell, cmdl%cell%ionic_positions(:,iatom,ispec), xcart(:,i))
             i=i+1
          end do
       end do

       do k=1,3
          L(k) = maxval(xcart(k,:)) - minval(xcart(k,:))
       end do

       M = 0.0_dp
       do k=1,3
          if (L(k) /= 0.0_dp) then
             M(k,k) = 0.01_dp * precond%exp_r_NN
             vcart(k,:) = sin(xcart(k,:) / L(k))
          else
             vcart(k,:) = xcart(k,:)
          end if
       end do

       vcart = matmul(M,vcart)

       precond%exp_mu = 1.0_dp

       call geom_precond_exp_calc_cart(cmdl)

       i = 1
       do ispec=1,cmdl%cell%num_species
          do iatom=1,cmdl%cell%num_ions_in_species(ispec)
             g0cart(:,i) = -cmdl%forces(:,iatom,ispec)
             i=i+1
          end do
       end do

       x_vec(1:9) = reshape(cmdl%strain,(/9/))
       i = 1
       do ispec=1,cmdl%cell%num_species
          do iatom=1,cmdl%cell%num_ions_in_species(ispec)
             call castep_cell_cart_to_frac(cmdl%cell, xcart(:,i)+vcart(:,i), x_vec(9+(i-1)*3+1:9+i*3))
             i=i+1
          end do
       end do

    end if

    call comms_bcast(pub_root_proc_id,x_vec,ndim)

    call geom_xvec_to_mdl(x_vec,cmdl,mdl)

    call geom_get_forces(cmdl,mdl,is_cmplx=is_cmplx)

    if (pub_on_root) then

       i = 1
       do ispec=1,cmdl%cell%num_species
          do iatom=1,cmdl%cell%num_ions_in_species(ispec)
             g1cart(:,i) = -cmdl%forces(:,iatom,ispec)
             i=i+1
          end do
       end do

       i = 1
       do ispec=1,cmdl%cell%num_species
          do iatom=1,cmdl%cell%num_ions_in_species(ispec)
             call castep_cell_cart_to_frac(cmdl%cell, xcart(:,i), cmdl%cell%ionic_positions(:,iatom,ispec))
             i=i+1
          end do
       end do

       mu = dot_product(reshape(vcart, (/3*cmdl%cell%num_ions/)), reshape(g1cart-g0cart, (/3*cmdl%cell%num_ions/))) / &
          & dot_product(reshape(vcart, (/3*cmdl%cell%num_ions/)), matmul(precond%P, reshape(vcart, (/3*cmdl%cell%num_ions/))))

       deallocate(xcart,stat=ierr)
       call utils_dealloc_check('geom_precond_exp_estimate_mu','xcart',ierr)

       deallocate(vcart,stat=ierr)
       call utils_dealloc_check('geom_precond_exp_estimate_mu','vcart',ierr)

       deallocate(g0cart,stat=ierr)
       call utils_dealloc_check('geom_precond_exp_estimate_mu','g0cart',ierr)

       deallocate(g1cart,stat=ierr)
       call utils_dealloc_check('geom_precond_exp_estimate_mu','g1cart',ierr)

    end if

    deallocate(x_vec,stat=ierr)
    call utils_dealloc_check('geom_precond_exp_estimate_mu','x_vec',ierr)

    constrain_ions = constrain_ions_orig

    call comms_bcast(pub_root_proc_id,mu)

  end function geom_precond_exp_estimate_mu

  subroutine geom_precond_ff_initialize(cmdl,mdl)

    use simulation_cell, only: castep_model
    use model_type, only: MODEL
    use comms, only: pub_on_root
    use utils, only: utils_abort, utils_alloc_check
    use rundat, only: pub_geom_precond_ff_c_stab, pub_geom_precond_ff_r_cut

    type(castep_model), intent(inout) :: cmdl
    type(MODEL),        intent(inout) :: mdl

    integer :: ierr
    integer :: i, ispec, iatom, row

    if (.not. precond%initialized) then
       if (pub_on_root) then
          precond%initialized = .true.

          if (constrain_ions) then
             allocate(precond%ion_constraint_type(mdl%nat),stat=ierr)
             call utils_alloc_check('geom_precond_ff_initialize','precond%ion_constraint_type',ierr)
             do i=1,mdl%nat
                precond%ion_constraint_type(i) = mdl%elements(i)%ion_constraint_type
                select case (precond%ion_constraint_type(i))
                case ('NONE')
                case ('FIXED')
                case default
                   call utils_abort('Error in geom_precond_ff_initialize: &
                                   & illegal value for precond%ion_constraint_type(i)=', i)
                end select
             end do
          end if

          allocate(precond%P(3*cmdl%cell%num_ions, 3*cmdl%cell%num_ions),stat=ierr)
          call utils_alloc_check('geom_precond_ff_initialize','precond%P',ierr)
          precond%P=0.0_dp

          allocate(precond%rmatrix(3, cmdl%cell%num_ions),stat=ierr)
          call utils_alloc_check('geom_precond_ff_initialize','precond%rmatrix',ierr)
          precond%rmatrix=0.0_dp

          allocate(precond%dmatrix(cmdl%cell%num_ions, cmdl%cell%num_ions),stat=ierr)
          call utils_alloc_check('geom_precond_ff_initialize','precond%dmatrix',ierr)
          precond%dmatrix=0.0_dp

          precond%ff_c_stab = pub_geom_precond_ff_c_stab
          if (precond%ff_c_stab <= 0.0_dp) then
             call utils_abort('pub_geom_precond_ff_c_stab must be larger than 0!')
          end if

          precond%ff_r_cut = pub_geom_precond_ff_r_cut
          if (precond%ff_r_cut <= 0.0_dp) then
             call utils_abort('pub_geom_precond_ff_r_cut must be larger than 0!')
          end if

          allocate(precond%ff_row(cmdl%cell%num_ions),stat=ierr)
          call utils_alloc_check('geom_precond_ff_initialize','precond%ff_row',ierr)
          i = 1
          do ispec=1,cmdl%cell%num_species
             select case(trim(cmdl%cell%species_symbol(ispec)))
             case('H','He')
                row = 1
             case('Li','Be','B ','C ','N ','O ','F ','Ne')
                row = 2
             case('Na','Mg','Al','Si','P ','S ','Cl','Ar')
                row = 3
             case default
                call utils_abort('row > 3 element is not supported in geom_precond_ff_initialize')
             end select
             do iatom=1,cmdl%cell%num_ions_in_species(ispec)
                precond%ff_row(i) = row
                i=i+1
             end do
          end do

          allocate(precond%ff_alpha(3,3),stat=ierr)
          call utils_alloc_check('geom_precond_ff_initialize','precond%ff_alpha',ierr)
          precond%ff_alpha = reshape( (/ 1.0000_dp, 0.3949_dp, 0.3949_dp, &
                                       & 0.3949_dp, 0.2800_dp, 0.2800_dp, &
                                       & 0.3949_dp, 0.2800_dp, 0.2800_dp  /), (/3,3/))

          allocate(precond%ff_r_ref(3,3),stat=ierr)
          call utils_alloc_check('geom_precond_ff_initialize','precond%ff_r_ref',ierr)
          precond%ff_r_ref = reshape( (/ 1.35_dp, 2.10_dp, 2.53_dp, &
                                       & 2.10_dp, 2.87_dp, 3.40_dp, &
                                       & 2.53_dp, 3.40_dp, 3.40_dp  /), (/3,3/))

          precond%ff_kbond = 0.45_dp
          precond%ff_kangle = 0.15_dp
          precond%ff_kdihedral = 0.005_dp
       end if
    end if

  end subroutine geom_precond_ff_initialize

  subroutine geom_precond_calc(cmdl)

    use simulation_cell, only: castep_model
    use utils, only: utils_abort
    use rundat, only: pub_geom_precond_type

    type(castep_model), intent(in) :: cmdl

    select case(pub_geom_precond_type)
    case('NONE')
       return
    case('ID')
       call geom_precond_id_calc_frac(cmdl)
       call geom_precond_cart_to_frac(cmdl)
    case('EXP')
       call geom_precond_exp_calc_cart(cmdl)
       call geom_precond_cart_to_frac(cmdl)
    case('FF')
       call geom_precond_ff_calc_cart(cmdl)
       call geom_precond_cart_to_frac(cmdl)
    case default
       call utils_abort('Error in geom_precond_apply_vec - unrecognised preconditioner request='// &
               & trim(pub_geom_precond_type)//'.')
    end select

  end subroutine geom_precond_calc

  subroutine geom_precond_id_calc_frac(cmdl)

    use simulation_cell, only: castep_model

    type(castep_model), intent(in)    :: cmdl

    integer                          :: ierr

    real (kind=dp)                   :: avg_mass
    integer                          :: i, ispec, icomp


    avg_mass=0.0_dp
    do ispec=1,cmdl%ref_cell%num_species
       avg_mass=avg_mass+cmdl%ref_cell%species_mass(ispec)*cmdl%ref_cell%num_ions_in_species(ispec)
    end do

    avg_mass=1822.888485_dp*avg_mass/cmdl%ref_cell%num_ions

    do i=1,cmdl%ref_cell%num_ions
       precond%P(3*(i-1)+1:3*(i-1)+3,3*(i-1)+1:3*(i-1)+3) = avg_mass*new_frequency_est**2
    end do

    if (constrain_ions) then
       do i=1,cmdl%ref_cell%num_ions
          select case (precond%ion_constraint_type(i))
          case ('NONE')
             continue
          case ('FIXED')
             do icomp=1,3
                precond%P(3*(i-1)+icomp,:3*(i-1)+icomp-1) = 0.0_dp
                precond%P(3*(i-1)+icomp,3*(i-1)+icomp+1:) = 0.0_dp
                precond%P(:3*(i-1)+icomp-1,3*(i-1)+icomp) = 0.0_dp
                precond%P(3*(i-1)+icomp+1:,3*(i-1)+icomp) = 0.0_dp
             end do
          end select
       enddo
    end if

  end subroutine geom_precond_id_calc_frac

  subroutine geom_precond_exp_calc_cart(cmdl)

    use simulation_cell, only: castep_model

    type(castep_model), intent(in) :: cmdl

    integer :: i, j, k, icomp
    real (kind=dp) :: coeff

    call geom_precond_get_dmatrix(cmdl)
    precond%P = 0.0_dp

    do i=1,cmdl%cell%num_ions
       do j=i+1,cmdl%cell%num_ions
          if (precond%dmatrix(i,j) <= precond%exp_r_cut) then
             coeff = precond%exp_mu * exp(-precond%exp_A * (precond%dmatrix(i,j)/precond%exp_r_NN - 1.0_dp))
             do k=1,3
                precond%P(3*(i-1)+k,3*(j-1)+k) = precond%P(3*(i-1)+k,3*(j-1)+k) - coeff
                precond%P(3*(j-1)+k,3*(i-1)+k) = precond%P(3*(j-1)+k,3*(i-1)+k) - coeff
                precond%P(3*(i-1)+k,3*(i-1)+k) = precond%P(3*(i-1)+k,3*(i-1)+k) + coeff
                precond%P(3*(j-1)+k,3*(j-1)+k) = precond%P(3*(j-1)+k,3*(j-1)+k) + coeff
             end do
          end if
       end do
    end do

    do i=1,cmdl%cell%num_ions
       do k=1,3
          precond%P(3*(i-1)+k,3*(i-1)+k) = precond%P(3*(i-1)+k,3*(i-1)+k) + precond%exp_mu * precond%exp_c_stab
       end do
    end do

    if (constrain_ions) then
       do i=1,cmdl%cell%num_ions
          select case (precond%ion_constraint_type(i))
          case ('NONE')
             continue
          case ('FIXED')
             do icomp=1,3
                precond%P(3*(i-1)+icomp,:3*(i-1)+icomp-1) = 0.0_dp
                precond%P(3*(i-1)+icomp,3*(i-1)+icomp+1:) = 0.0_dp
                precond%P(:3*(i-1)+icomp-1,3*(i-1)+icomp) = 0.0_dp
                precond%P(3*(i-1)+icomp+1:,3*(i-1)+icomp) = 0.0_dp
             end do
          end select
       enddo
    end if

  end subroutine geom_precond_exp_calc_cart

  subroutine geom_precond_ff_calc_cart(cmdl)

    use simulation_cell, only: castep_model

    type(castep_model), intent(in) :: cmdl

    integer :: i, j, k, l, icomp
    real (kind=dp) :: rhoij, rhoijk, rhoijkl, fk
    integer :: index6(6), index9(9), index12(12)
    real (kind=dp) :: grad6(6), grad9(9), grad12(12)

    call geom_precond_get_dmatrix(cmdl)
    precond%P = 0.0_dp

    do i=1,cmdl%cell%num_ions
       do j=i+1,cmdl%cell%num_ions
          if (precond%dmatrix(i,j) <= precond%ff_r_cut) then
             rhoij = geom_precond_ff_calc_rho(i,j)
             fk = precond%ff_kbond * rhoij
             index6 = (/ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                         3*(j-1)+1,3*(j-1)+2,3*(j-1)+3 /)
             grad6 = geom_precond_get_distance_grad(i,j,cmdl)
             precond%P(index6,index6) = precond%P(index6,index6) + fk*tensor_product(grad6, grad6, 6)
          end if
          do k=j+1,cmdl%cell%num_ions
             if (precond%dmatrix(i,j) <= precond%ff_r_cut) then
                if (precond%dmatrix(i,k) <= precond%ff_r_cut) then
                   rhoijk = geom_precond_ff_calc_rho(j,i) * geom_precond_ff_calc_rho(i,k)
                   fk = precond%ff_kangle * rhoijk
                   index9 = (/ 3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                               3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                               3*(k-1)+1,3*(k-1)+2,3*(k-1)+3 /)
                   grad9 = geom_precond_get_angle_grad(j,i,k,cmdl)
                   precond%P(index9,index9) = precond%P(index9,index9) + fk*tensor_product(grad9, grad9, 9)
                end if
                if (precond%dmatrix(j,k) <= precond%ff_r_cut) then
                   rhoijk = geom_precond_ff_calc_rho(i,j) * geom_precond_ff_calc_rho(j,k)
                   fk = precond%ff_kangle * rhoijk
                   index9 = (/ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                               3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                               3*(k-1)+1,3*(k-1)+2,3*(k-1)+3 /)
                   grad9 = geom_precond_get_angle_grad(i,j,k,cmdl)
                   precond%P(index9,index9) = precond%P(index9,index9) + fk*tensor_product(grad9, grad9, 9)
                end if
             end if
             if (precond%dmatrix(i,k) <= precond%ff_r_cut) then
                rhoijk = geom_precond_ff_calc_rho(i,k) * geom_precond_ff_calc_rho(k,j)
                index9 = (/ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                            3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                            3*(j-1)+1,3*(j-1)+2,3*(j-1)+3 /)
                grad9 = geom_precond_get_angle_grad(i,k,j,cmdl)
                precond%P(index9,index9) = precond%P(index9,index9) + fk*tensor_product(grad9, grad9, 9)
             end if
             do l=k+1,cmdl%cell%num_ions
                if (precond%dmatrix(i,j) <= precond%ff_r_cut) then
                   if ( (precond%dmatrix(i,k) <= precond%ff_r_cut) .and. (precond%dmatrix(j,l) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(k,i) * geom_precond_ff_calc_rho(i,j) * geom_precond_ff_calc_rho(j,l)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                                   3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                                   3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                                   3*(l-1)+1,3*(l-1)+2,3*(l-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(k,i,j,l,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   else if ( (precond%dmatrix(i,l) <= precond%ff_r_cut) .and. (precond%dmatrix(j,k) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(l,i) * geom_precond_ff_calc_rho(i,j) * geom_precond_ff_calc_rho(j,k)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(l-1)+1,3*(l-1)+2,3*(l-1)+3, &
                                   3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                                   3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                                   3*(k-1)+1,3*(k-1)+2,3*(k-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(l,i,j,k,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   end if
                end if
                if (precond%dmatrix(i,k) <= precond%ff_r_cut) then
                   if ( (precond%dmatrix(i,j) <= precond%ff_r_cut) .and. (precond%dmatrix(k,l) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(j,i) * geom_precond_ff_calc_rho(i,k) * geom_precond_ff_calc_rho(k,l)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                                   3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                                   3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                                   3*(l-1)+1,3*(l-1)+2,3*(l-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(j,i,k,l,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   else if ( (precond%dmatrix(i,l) <= precond%ff_r_cut) .and. (precond%dmatrix(j,k) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(l,i) * geom_precond_ff_calc_rho(i,k) * geom_precond_ff_calc_rho(k,j)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(l-1)+1,3*(l-1)+2,3*(l-1)+3, &
                                   3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                                   3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                                   3*(j-1)+1,3*(j-1)+2,3*(j-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(l,i,k,j,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   end if
                end if
                if (precond%dmatrix(i,l) <= precond%ff_r_cut) then
                   if ( (precond%dmatrix(i,j) <= precond%ff_r_cut) .and. (precond%dmatrix(k,l) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(j,i) * geom_precond_ff_calc_rho(i,l) * geom_precond_ff_calc_rho(l,k)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                                   3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                                   3*(l-1)+1,3*(l-1)+2,3*(l-1)+3, &
                                   3*(k-1)+1,3*(k-1)+2,3*(k-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(j,i,l,k,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   else if ( (precond%dmatrix(i,k) <= precond%ff_r_cut) .and. (precond%dmatrix(j,l) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(k,i) * geom_precond_ff_calc_rho(i,l) * geom_precond_ff_calc_rho(l,j)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                                   3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                                   3*(l-1)+1,3*(l-1)+2,3*(l-1)+3, &
                                   3*(j-1)+1,3*(j-1)+2,3*(j-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(k,i,l,j,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   end if
                end if
                if (precond%dmatrix(j,k) <= precond%ff_r_cut) then
                   if ( (precond%dmatrix(i,j) <= precond%ff_r_cut) .and. (precond%dmatrix(k,l) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(i,j) * geom_precond_ff_calc_rho(j,k) * geom_precond_ff_calc_rho(k,l)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                                   3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                                   3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                                   3*(l-1)+1,3*(l-1)+2,3*(l-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(i,j,k,l,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   else if ( (precond%dmatrix(j,l) <= precond%ff_r_cut) .and. (precond%dmatrix(i,k) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(l,j) * geom_precond_ff_calc_rho(j,k) * geom_precond_ff_calc_rho(k,i)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(l-1)+1,3*(l-1)+2,3*(l-1)+3, &
                                   3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                                   3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                                   3*(i-1)+1,3*(i-1)+2,3*(i-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(l,j,k,i,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   end if
                end if
                if (precond%dmatrix(j,l) <= precond%ff_r_cut) then
                   if ( (precond%dmatrix(i,j) <= precond%ff_r_cut) .and. (precond%dmatrix(k,l) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(i,j) * geom_precond_ff_calc_rho(j,l) * geom_precond_ff_calc_rho(l,k)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                                   3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                                   3*(l-1)+1,3*(l-1)+2,3*(l-1)+3, &
                                   3*(k-1)+1,3*(k-1)+2,3*(k-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(i,j,l,k,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   else if ( (precond%dmatrix(j,k) <= precond%ff_r_cut) .and. (precond%dmatrix(i,l) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(k,j) * geom_precond_ff_calc_rho(j,l) * geom_precond_ff_calc_rho(l,i)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                                   3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                                   3*(l-1)+1,3*(l-1)+2,3*(l-1)+3, &
                                   3*(i-1)+1,3*(i-1)+2,3*(i-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(k,j,l,i,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   end if
                end if
                if (precond%dmatrix(k,l) <= precond%ff_r_cut) then
                   if ( (precond%dmatrix(i,k) <= precond%ff_r_cut) .and. (precond%dmatrix(j,l) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(i,k) * geom_precond_ff_calc_rho(k,l) * geom_precond_ff_calc_rho(l,j)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3, &
                                   3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                                   3*(l-1)+1,3*(l-1)+2,3*(l-1)+3, &
                                   3*(j-1)+1,3*(j-1)+2,3*(j-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(i,k,l,j,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   else if ( (precond%dmatrix(j,k) <= precond%ff_r_cut) .and. (precond%dmatrix(i,l) <= precond%ff_r_cut) ) then
                      rhoijkl = geom_precond_ff_calc_rho(j,k) * geom_precond_ff_calc_rho(k,l) * geom_precond_ff_calc_rho(l,i)
                      fk = precond%ff_kdihedral * rhoijkl
                      index12 = (/ 3*(j-1)+1,3*(j-1)+2,3*(j-1)+3, &
                                   3*(k-1)+1,3*(k-1)+2,3*(k-1)+3, &
                                   3*(l-1)+1,3*(l-1)+2,3*(l-1)+3, &
                                   3*(i-1)+1,3*(i-1)+2,3*(i-1)+3 /)
                      grad12 = geom_precond_get_dihedral_grad(j,k,l,i,cmdl)
                      precond%P(index12,index12) = precond%P(index12,index12) + fk*tensor_product(grad12, grad12, 12)
                   end if
                end if
             end do
          end do
       end do
    end do

    do i=1,cmdl%cell%num_ions
       do k=1,3
          precond%P(3*(i-1)+k,3*(i-1)+k) = precond%P(3*(i-1)+k,3*(i-1)+k) + precond%ff_c_stab
       end do
    end do

    if (constrain_ions) then
       do i=1,cmdl%cell%num_ions
          select case (precond%ion_constraint_type(i))
          case ('NONE')
             continue
          case ('FIXED')
             do icomp=1,3
                precond%P(3*(i-1)+icomp,:3*(i-1)+icomp-1) = 0.0_dp
                precond%P(3*(i-1)+icomp,3*(i-1)+icomp+1:) = 0.0_dp
                precond%P(:3*(i-1)+icomp-1,3*(i-1)+icomp) = 0.0_dp
                precond%P(3*(i-1)+icomp+1:,3*(i-1)+icomp) = 0.0_dp
             end do
          end select
       enddo
    end if

  end subroutine geom_precond_ff_calc_cart

  function geom_precond_ff_calc_rho(i,j) result(rho)

    integer, intent(in) :: i, j
    real (kind=dp)      :: rho

    rho = exp( precond%ff_alpha(precond%ff_row(i), precond%ff_row(j)) * &
              (precond%ff_r_ref(precond%ff_row(i), precond%ff_row(j))**2-precond%dmatrix(i,j)**2) )

  end function geom_precond_ff_calc_rho

  subroutine geom_precond_cart_to_frac(cmdl)

    use simulation_cell, only: castep_model

    type(castep_model), intent(in) :: cmdl

    integer :: i, j
    real (kind=dp) :: block(3,3)

    do i=1,cmdl%cell%num_ions
       do j=1,cmdl%cell%num_ions
          block = precond%P(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3)
          if (all(block == 0.0_dp)) cycle
          block = matmul(cmdl%cell%real_lattice,matmul(block,transpose(cmdl%cell%real_lattice)))
          precond%P(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3) = block
       end do
    end do

  end subroutine geom_precond_cart_to_frac

  function geom_precond_get_distance_grad(i,j,cmdl) result(grad)

    use simulation_cell, only: castep_model

    integer, intent(in)           :: i, j
    type(castep_model), intent(in) :: cmdl
    real (kind=dp) :: grad(6)

    real (kind=dp) :: rij(3), dij, eij(3)

    rij = geom_opt_min_image_cart(precond%rmatrix(:,i), precond%rmatrix(:,j))
    dij = precond%dmatrix(i,j)
    eij = rij/dij

    grad(1:3) = eij
    grad(4:6) = -eij

  end function geom_precond_get_distance_grad

  function geom_precond_get_angle_grad(i,j,k,cmdl) result(grad)

    use simulation_cell, only: castep_model

    integer, intent(in)           :: i, j, k
    type(castep_model), intent(in) :: cmdl
    real (kind=dp) :: grad(9)

    integer :: n
    real (kind=dp) :: rij(3), dij, eij(3)
    real (kind=dp) :: rkj(3), dkj, ekj(3)
    real (kind=dp) :: eijekj
    real (kind=dp) :: a, sina
    real (kind=dp) :: eye(3,3), Pij(3,3), Qij(3,3), Pkj(3,3), Qkj(3,3)

    rij = geom_opt_min_image_cart(precond%rmatrix(:,i), precond%rmatrix(:,j))
    dij = precond%dmatrix(i,j)
    eij = rij/dij

    rkj = geom_opt_min_image_cart(precond%rmatrix(:,k), precond%rmatrix(:,j))
    dkj = precond%dmatrix(k,j)
    ekj = rkj/dkj

    eijekj = dot_product(eij, ekj)

    if (abs(eijekj) > 1.0_dp) then
       eijekj = sign(1.0_dp, eijekj)
    end if

    a = acos(eijekj)
    sina = sin(a)

    grad = 0.0_dp

    if (abs(sina) > 1.0e-3_dp) then
       eye = 0.0_dp
       do n=1,3
          eye(n,n) = 1.0_dp
       end do
       Pij = tensor_product(eij, eij, 3)
       Qij = eye - Pij
       Pkj = tensor_product(ekj, ekj, 3)
       Qkj = eye - Pkj

       grad(1:3) = -matmul(Qij,ekj)/sina/dij
       grad(7:9) = -matmul(Qkj,eij)/sina/dkj
       grad(4:6) = -grad(1:3)-grad(7:9)
    end if

  end function geom_precond_get_angle_grad

  function geom_precond_get_dihedral_grad(i,j,k,l,cmdl) result(grad)

    use simulation_cell, only: castep_model

    integer, intent(in)           :: i, j, k, l
    type(castep_model), intent(in) :: cmdl
    real (kind=dp) :: grad(12)

    integer :: n
    real (kind=dp) :: rij(3)
    real (kind=dp) :: rkj(3), dkj, dkj2
    real (kind=dp) :: rkl(3)
    real (kind=dp) :: rijrkj, rkjrkl
    real (kind=dp) :: rmj(3), dmj2
    real (kind=dp) :: rnk(3), dnk2
    real (kind=dp) :: dddri(3), dddrl(3)

    rij = geom_opt_min_image_cart(precond%rmatrix(:,i), precond%rmatrix(:,j))
    rkj = geom_opt_min_image_cart(precond%rmatrix(:,k), precond%rmatrix(:,j))
    dkj = precond%dmatrix(k,j)
    dkj2 = dkj*dkj
    rkl = geom_opt_min_image_cart(precond%rmatrix(:,k), precond%rmatrix(:,l))

    rijrkj = dot_product(rij, rkj)
    rkjrkl = dot_product(rkj, rkl)

    rmj = cross_product(rij, rkj)
    dmj2 = dot_product(rmj, rmj)
    rnk = cross_product(rkj, rkl)
    dnk2 = dot_product(rnk, rnk)

    dddri = dkj/dmj2*rmj
    dddrl = -dkj/dnk2*rnk

    grad(1:3) = dddri
    grad(4:6) = (rijrkj/dkj2-1.0_dp)*dddri-rkjrkl/dkj2*dddrl
    grad(7:9) = (rkjrkl/dkj2-1.0_dp)*dddrl-rijrkj/dkj2*dddri
    grad(10:12) = dddrl

  end function geom_precond_get_dihedral_grad

  subroutine geom_precond_get_dmatrix(cmdl)

    use simulation_cell, only: castep_model, castep_cell_frac_to_cart
    use utils, only: utils_alloc_check

    type(castep_model), intent(in) :: cmdl

    real (kind=dp) :: d(3)
    integer :: i, j
    integer :: ispec, iatom

    if (.not. allocated(precond%rmatrix)) then
       call utils_alloc_check('geom_precond_get_dmatrix','precond%rmatrix',1)
    end if

    if (.not. allocated(precond%dmatrix)) then
       call utils_alloc_check('geom_precond_get_dmatrix','precond%dmatrix',1)
    end if

    i = 1
    do ispec=1,cmdl%cell%num_species
       do iatom=1,cmdl%cell%num_ions_in_species(ispec)
          call castep_cell_frac_to_cart(cmdl%cell, cmdl%cell%ionic_positions(:,iatom,ispec), precond%rmatrix(:,i))
          i=i+1
       end do
    end do

    do i=1, cmdl%cell%num_ions
       do j=i+1, cmdl%cell%num_ions
          d = geom_opt_min_image_cart(precond%rmatrix(:,i), precond%rmatrix(:,j))
          precond%dmatrix(i,j) = sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
          precond%dmatrix(j,i) = precond%dmatrix(i,j)
       end do
    end do

  end subroutine geom_precond_get_dmatrix

  function tensor_product(a,b,n)

    integer, intent(in) :: n
    real(kind=dp), dimension(n), intent(in) :: a,b
    real(kind=dp), dimension(n,n) :: tensor_product

    ! Local variables
    integer :: i,j

    do j = 1, n
       do i = 1, n
          tensor_product(i,j) = a(i) * b(j)
       end do
    end do

  end function tensor_product

  function cross_product(a,b)

    real(kind=dp), dimension(3), intent(in) :: a,b
    real(kind=dp), dimension(3) :: cross_product

    cross_product(1) = a(2)*b(3) - a(3)*b(2)
    cross_product(2) = a(3)*b(1) - a(1)*b(3)
    cross_product(3) = a(1)*b(2) - a(2)*b(1)

  end function cross_product

end module geometry_optimiser
