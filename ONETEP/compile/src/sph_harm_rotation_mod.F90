! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!==============================================================================!
!                                                                              !
!                      Spherical harmonic rotation module                      !
!                                                                              !
!                         James C. Womack, 2018 - 20XX                         !
!                                                                              !
!------------------------------------------------------------------------------!
! # MODULE DESCRIPTION
! A module providing subroutines for the rotation of (real) spherical          !
! harmonics using Wigner D-matrices.                                           !
!------------------------------------------------------------------------------!
! # BACKGROUND                                                                 !
! A Wigner D-matrix is the representation of a rotation operator in the Euler  !
! angle parameterization in the basis of angular momentum eigenstates of a     !
! particular angular momentum quantum number, l.                               !
!                                                                              !
! The spherical harmonics, Y_{l,m} are position representations of angular     !
! momentum eigenstates with angular momentum l. For each l there are 2l+1      !
! different m values.                                                          !
!                                                                              !
! The real spherical harmonics Z_{l,m} are related to the usual (complex)      !
! spherical harmonics by a unitary transformation, i.e. they are constructed   !
! as linear combinations of the complex spherical harmonics. As their name     !
! suggests, Z_{l,m} are real functions.                                        !
!                                                                              !
! The Wigner D-matrix for angular momentum quantum number l, D_{mu,m}^{l},     !
! defines the transformation of spherical harmonics of angular momentum l      !
! between two coordinate systems related by an Euler rotation. This applies    !
! to the complex spherical harmonics.                                          !
!                                                                              !
! The real spherical harmonics (RSHs) can be rotated between coordinate        !
! systems (related by Euler rotations) by combining the unitary transformation !
! matrix, U_{mu,m}^{l}, which converts between real and complex spherical      !
! harmonics (CSHs), and the  Wigner D-matrix D_{mu,m}^{l} (for a given angular !
! momentum, l).                                                                !
!                                                                              !
! ## Definitions                                                               !
! The definitions used in this module are largely based on the following       !
! sources:                                                                     !
! * Arfken, "Mathematical Methods for Physicists", 7th ed., ch. 15             !
!   (spherical harmonics)                                                      !
! * M.A. Morrison and G.A. Parker, Aust. J. Phys. 40, 465 (1987)               !
!   (Wigner D-matrices and rotation conventions)                               !
! * J. Fernandez Rico, R. Lopez, and G. Ramirez, J. Chem. Phys. 91, 4204       !
!   (1989) (unitary transformation for RSH <--> CSH)                           !
!                                                                              !
! The complex spherical harmonics (CSHs) are defined as follows, including the !
! Condon-Shortley phase factor                                                 !
!                                                                              !
!   Y_{l,m}(theta,phi) = (-1)^m N_{l,m} P_{l,m}(cos(theta)) exp(i m phi)       !
!                                                                              !
! where                                                                        !
! * the coordinate system is spherical polar, with theta the polar angle       !
!   (0 <= theta <= pi) and phi the azimuthal angle (0 <= phi < 2*pi)           !
! * N_{l,m} is a normalization factor:                                         !
!     sqrt( (2*l+1)/(4*pi) * (l-m)!/(l+m)! )                                   !
! * P_{l,m}(x) is an associated Legendre function with argument x = cos(theta) !
! * (-1)^m is the Condon-Shortley phase factor                                 !
!                                                                              !
! Note that unlike in Arfken, 7th edn., we do not include the phase factor in  !
! the definition of P_{l,m} -- see p. 758 of Arfken 7th edn.                   !
!                                                                              !
! The RSHs are defined in terms of the CSHs as                                 !
!                                                                              !
!             | Y_{l,0}                                    m = 0               !
!   Z_{l,m} = | i/sqrt(2) (Y_{l,-|m|} - (-1)^m Y_{l,|m|})  m < 0               !
!             | 1/sqrt(2) (Y_{l,-m}   + (-1)^m Y_{l,m})    m > 0               !
!                                                                              !
! which corresponds to a 2l+1 x 2l+1 unitary transformation matrix for each l: !
!                                                                              !
!                     / ...                        ... \                       !
!                     |      i                  1      |                       !
!                     |         i            1         |                       !
!   U^{l} = 1/sqrt(2) |            sqrt(2)             |                       !
!                     |         i           -1         |                       !
!                     |     -i                  1      |                       !
!                     \ ...                        ... /                       !
!                                                                              !
! as defined in J. Fernandez Rico et al, JCP, 1989 (Eq. B8).                   !
!                                                                              !
! To rotate CSHs between coordinate systems related by an Euler rotation, we   !
! can apply the Wigner D-matrix for a given l:                                 !
!                                                                              !
!   Y^{l,m}(theta,phi)                                                         !
!     = \sum_{mu=-l}^{+l} Y_{l,mu}(theta',phi') D_{mu,m}^{l}(alpha,beta,gamma) !
!                                                                              !
! where theta, phi are angular coordinates in the new coordinate system and    !
! theta', phi' are in the old coordinate system, i.e. we are expressing a set  !
! of spherical harmonics in the unprimed coordinate system in terms of a set   !
! of spherical harmonics in the primed coordinate system.                      !
!                                                                              !
! The two coordinate systems are related by an Euler rotation defined by the   !
! three angles alpha, beta, gamma.                                             !
!                                                                              !
! To rotate the RSHs, we can use the rotation matrix                           !
!   R^{l} = (U^{l})* D^{l} U^{l}                                               !
! which is formed from the CSH <--> RSH transformation matrix and the          !
! Wigner D-matrix for a given l. (U^{l})* is the conjugate transpose of U^{l}. !
!                                                                              !
! The Wigner D-matrices in this module are constructed based on the definition !
! from Morrison and Parker, Aust. J. Phys. 1987 (appendix B, Eqs. B6 -- B10):  !
!                                                                              !
!   D_{mu,m}^{l}(alpha,beta,gamma)                                             !
!     = exp(-i mu alpha) d_{mu,m}^{l}(beta) exp(-i m gamma)                    !
!                                                                              !
! with the d_{mu,m}^{l} the so-called "small-d" matrix, which can be defined   !
! in terms of the hypergeometric function 2F1 (Eq. B7) or the Jacobi           !
! polynomials (Eq. B9).                                                        !
!                                                                              !
! ## Conventions                                                               !
! When dealing with Euler rotations of coordinate systems, it is important to  !
! clearly define the conventions used.                                         !
!                                                                              !
! The key conventions used in this module are:                                 !
! * active rotation                                                            !
! * about fixed coordinate axes                                                !
! * all Euler angles are positive, right-handed rotations wrt rotation axes    !
!                                                                              !
! These are the conventions used in Morrison and Parker, Aust. J. Phys. 1987,  !
! and are implicit in the definition of the Wigner D-matrix above.             !
!                                                                              !
! Using these conventions, the Euler rotation about the fixed Cartesian        !
! axes e_x, e_y, e_z proceeds as follows:                                      !
! 1. Rotate system by gamma about e_z (0 <= gamma <= 2*pi)                     !
! 2. Rotate system by beta about e_y  (0 <= beta <= pi)                        !
! 3. Rotate system by alpha about e_z (0 <= alpha <= 2*pi)                     !
!                                                                              !
! The rotation is active, so describes the rotation of a vector in a fixed     !
! coordinate frame.                                                            !
!                                                                              !
! ## Motivation                                                                !
! This module was originally developed to allow the rotation of RSHs           !
! associated with spherical wave (SW) basis functions (see the spherical_wave  !
! module), as used in the SW resolution-of-the-identity (SWRI) for evaluating  !
! metric matrix elements in the sw_resolution_of_identity module.              !
!                                                                              !
! If the metric matrix elements are expressed in a spherical polar coordinate  !
! system with the z-axis pointing between atomic centres (and the RSHs         !
! associated with SWs on each centre aligned in this coordinate system), it    !
! becomes possible to separate the 3-D integral over SWs into a 2-D            !
! integration over r,theta and a 1-D integration over phi, which can be done   !
! analytically.                                                                !
!                                                                              !
! The SW basis used in ONETEP is aligned along the Cartesian z-axis of the     !
! simulation cell, so in order to evaluate the metric matrix in this basis     !
! we must rotate (the RSH component of) the SWs between the "integral frame"   !
! (with z-axis pointing along the interatomic vector) and the "ONETEP frame".  !
!                                                                              !
! This can be achieved using the RSH rotation matrices defined above:          !
! 1. Determine active rotation which corresponds to rotating a vector pointing !
!    along z-axis of the integral coordinate frame into alignment with the     !
!    ONETEP coordinate frame z-axis.                                           !
! 2. Construct the RSH active rotation matrix from the Wigner D-matrix and     !
!    RSH <--> CSH transformation matrices, for a given l:                      !
!      R^{l} = (U^{l})* D^{l} U^{l}                                            !
! 3. Apply the rotation to the metric matrix (twice, one for each SW in the    !
!    integral:                                                                 !
!      (R^{lp})* M^{lp,lq} R^{lq}                                              !
!    where M^{lp,lq} is a block of a single atomblock of the metric matrix     !
!    containing matrix elements over SWs with angular momentum lp, lq.         !
!                                                                              !
! Each atomblock of the metric matrix corresponds to a pair of atoms, and thus !
! a single Euler rotation between integral frame and ONETEP frame. The         !
! atomblocks contain matrix elements for SWs with angular momentum up to       !
! lmax, so we need a R^{l} for each l value.                                   !
!                                                                              !
! The rotation of an entire atomblock of a metric matrix can be evaluated      !
! at once by constructing a single R matrix for all l. This is simply a block  !
! diagonal matrix with R^{l} submatrices along the diagonal.                   !
!                                                                              !
! ## Defining the Euler rotation                                               !
! As described in Morrison and Parker, Aust. J. Phys. 1987, the transformation !
! of CSHs between two coordinate systems can be understood in terms of the     !
! rotation of the axis of quantization for the CSHs. Since the same            !
! applies to RSHs, we can determine the Euler rotation needed to construct     !
! the RSH rotation matrix which transforms between the integral coordinate     !
! frame RSHs and ONETEP coordinate frame RSHs by finding the rotation which    !
! relates the axes of quantization (i.e. the z-axes of the two coordinate      !
! systems) of the two sets of RSHs.                                            !
!                                                                              !
! We want to define the Euler rotation for constructing the RSH rotation       !
! matrix to be used in the above process for evaluation of the ONETEP          !
! coordinate frame metric matrix in terms of the integral coordinate frame     !
! metric matrix. In this procedure, it is important to be clear about the      !
! definition of the fixed coordinate axes about which the active Euler         !
! rotation will occur.                                                         !
!                                                                              !
! In the implementation of spherical harmonic rotation in this module, we      !
! differ from the description in Morrison and Parker (III.6) by choosing the   !
! fixed axes of rotation to be those of the "final" coordinate system, i.e.    !
! the coordinate system of the spherical harmonics AFTER application of the    !
! rotation (Morrison and Parker use the "initial" coordinate system as the     !
! fixed axes of rotation, i.e. the coordinate system of the spherical          !
! harmonics BEFORE application of the rotation). In the specific case of       !
! rotating SWs in ONETEP, this is the ONETEP frame coordinate system, with     !
! atom-centered Cartesian axes parallel to the global Cartesian axes of the    !
! simulation cell.                                                             !
!                                                                              !
! When the ONETEP coordinate frame is selected as the fixed coordinate system  !
! for the rotation, the Euler rotation angles are simply the rotation angles   !
! for an active rotation which takes a vector aligned along the integral frame !
! coordinate system z-axis (defined in the fixed coordinate system) into       !
! alignment with the ONETEP frame coordinate system z-axis.                    !
!                                                                              !
! For the rotation between the ONETEP coordinate frame and integral coordinate !
! frame described above, only two Euler angles are needed: beta and gamma.     !
! This corresponds to rotating a vector aligned along the integral coordinate  !
! frame z-axis, e_z' to align with the ONETEP coordinate frame z-axis (with    !
! the ONETEP coordinate frame providing the fixed axes of rotation), i.e.      !
!   e_z = R^{a}(0,beta,gamma) e_z'                                             !
! with R^{a}(0,beta,gamma) an active Euler rotation operator.                  !
!                                                                              !
! This is a two step procedure:                                                !
! 1. Apply R^{a}_{z}(gamma) to rotate the vector anticlockwise about the       !
!    (fixed) e_z axis (0 <= gamma < 2*pi)                                      !
!   [ This results in a vector in the xz plane, with a negative projection     !
!   on the (fixed) e_x axis ]                                                  !
!                                                                              !
! 2. Apply R^{a}_{y}(beta) to rotate the vector (now in the xy-plane)          !
!    anticlockwise about the (fixed) e_y axis (0 <= beta <= pi)                !
!    [ This rotates the vector in the xz plane into alignment with the         !
!    (fixed) e_z axis ]                                                        !
!------------------------------------------------------------------------------!
! # NOTES                                                                      !
! ## TODO LIST                                                                 !
! * Consider implementing generic interfaces for routines to allow output of   !
!   matrices as Fortran arrays or dense matrix (DEM) objects (from dense       !
!   module)                                                                    !
!==============================================================================!
module sph_harm_rotation

  ! # MODULE-WIDE INCLUDES
  use comms, only : pub_on_root
  use constants, only: DP, stdout, stderr
  use rundat, only: pub_debug_on_root

  implicit none
  private

  ! # PUBLIC PROCEDURES
  ! Initialize module (check devel_code values, set module variables)
  public :: sph_harm_rot_initialise
  ! De-initialize module (reset module state to default)
  public :: sph_harm_rot_exit
  ! Return active Euler angles for a 'z to z' rotation between coord systems
  public :: sph_harm_rot_get_active_Euler_ang_z2z
  ! Return Wigner D-matrix for a given l and active Euler rotation (D)
  public :: sph_harm_rot_get_Wigner_Dmat
  ! Return unitary transformation matrix between RSHs and CSHs (U)
  public :: sph_harm_rot_get_RSH_CSH_Umat
  ! Return RSH rotation matrix for a given l and active Euler rotation
  ! (R = U* D U)
  public :: sph_harm_rot_get_RSH_Rmat
  ! Apply RSH rotation matrix for l between 0 and lmax to atomblock of a
  ! matrix with matrix elements containing RSHs (M' = R* M R)
  public :: sph_harm_rot_apply_RSH_Rmat_to_atomblock
  ! Unit tests for module procedures
  public :: sph_harm_rot_run_unit_tests

  ! # PUBLIC MODULE VARIABLES AND TYPES

  ! # PRIVATE MODULE VARIABLES AND TYPES
  ! Short code for reporting (useful for grepping output)
  ! This is also the block name for devel codes
  character(len=*), parameter :: module_short_code = "SHROT"
  ! Is the module initialized?
  logical, save :: sph_harm_rot_initialised = .false.
  ! ## Logical variables set by devel code values
  ! Should we run unit tests?
  !   devel_code: <module_short_code>:UNIT_TEST=[T/F]:<module_short_code>
  logical, save :: sph_harm_rot_unit_tests = .false.
  ! Should we output additional debugging information?
  !   devel_code: <module_short_code>:DEBUG=[T/F]:<module_short_code>
  logical, save :: sph_harm_rot_debug = .false.

  ! # INTERFACE BLOCKS

  ! # PROCEDURES
contains
  ! ## PUBLIC PROCEDURES

  ! Initialize module (check devel_code values, set module variables)
  subroutine sph_harm_rot_initialise()
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Initialize the sph_harm_rotation module:                                 !
    ! * If module has already been initialized, return immediately             !
    ! * If sph_harm_rotation module not needed, return immediately             !
    ! * Check relevant devel_codes and set module variables accordingly        !
    ! * If unit tests are requested, run these                                 !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! None                                                                     !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 07/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use constants, only: NORMAL
    use rundat, only: pub_devel_code, pub_output_detail, pub_use_sph_harm_rot
    use utils, only: utils_trace_in, utils_trace_out, &
         utils_abort, utils_devel_code

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "sph_harm_rot_initialise"

    ! Local variables
    integer :: ierror

    call utils_trace_in(myself)

    ! Check if sph_harm_rotation module is needed, if not return immediately
    if (.not.pub_use_sph_harm_rot .or. sph_harm_rot_initialised) then
      call utils_trace_out(myself)
      return
    end if

    if (pub_on_root .and. pub_output_detail > NORMAL) then
       write(stdout,"(a)") &
            "Initializing spherical harmonic rotation module."
    end if

    ! Check for relevant devel_code strings to set module variables
    sph_harm_rot_debug = &
         utils_devel_code(.false.,module_short_code,"DEBUG",pub_devel_code)
    if (pub_on_root .and. sph_harm_rot_debug) then
       write(stdout,"(a)") module_short_code//": debug mode active."
    end if
    sph_harm_rot_unit_tests = &
         utils_devel_code(.false.,module_short_code,"UNIT_TEST",pub_devel_code)
    if (pub_on_root .and. sph_harm_rot_unit_tests) then
       write(stdout,"(a)") module_short_code//": unit testing active."
    end if

    ! If requested, run unit tests for this module
    if (sph_harm_rot_unit_tests) then
       if (pub_on_root) then
          write(stdout,"(a)") module_short_code//": Running unit tests."
       end if
       call sph_harm_rot_run_unit_tests(ierror)
       if (ierror /= 0) then
          call utils_abort("Error in "//myself//": unit tests failed.")
       else
          if (pub_on_root) then
             write(stdout,"(a)") module_short_code//": Unit tests passed. "
          end if
       end if
    end if

    ! Module is now initialized
    sph_harm_rot_initialised = .true.

    call utils_trace_out(myself)

  end subroutine sph_harm_rot_initialise

  ! De-initialize module (reset module state to default)
  subroutine sph_harm_rot_exit()
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! De-initialize / clean up sph_harm_rotation module:                       !
    ! * If module not initialized, return immediately                          !
    ! * Restore module variables to default values                             !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! None                                                                     !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 07/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "sph_harm_rot_exit"

    ! Local variables

    call utils_trace_in(myself)

    ! Check if sph_harm_rotation module is initialized --> if not return now
    if (.not.sph_harm_rot_initialised) then
       call utils_trace_out(myself)
       return
    end if

    ! Restore module variables to default values
    sph_harm_rot_debug = .false.
    sph_harm_rot_unit_tests = .false.

    ! Module is not de-initialized
    sph_harm_rot_initialised = .false.

    call utils_trace_out(myself)

  end subroutine sph_harm_rot_exit

  subroutine sph_harm_rot_get_active_Euler_ang_z2z(beta_ang,gamma_ang,z_vec)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Compute the active Euler rotation angles for a rotation about fixed axes !
    ! which rotates a unit vector aligned along the z-axis of a rotated        !
    ! coordinate system into the z-axis of the coordinate system in which the  !
    ! unit vector is defined (i.e. the fixed coordinate axes of rotation).     !
    !                                                                          !
    ! This corresponds to the following rotation                               !
    !   e_z = R^{a}(0,beta,gamma) e_z'                                         !
    ! where e_z is the z-axis of the fixed coordinate system (in which z_vec   !
    ! is defined and which provides the fixed axes of rotation), e_z' is the   !
    ! z-axis of the rotated coordinate system (i.e. z_vec) and R^{a} is an     !
    ! active Euler rotation wrt fixed axes, as defined in Eq. 30 of            !
    !   M.A. Morrison and G.A. Parker, Aust. J. Phys. 40, 465 (1987)           !
    !                                                                          !
    ! Only two steps are needed (the alpha angle is arbitrary once the vector  !
    ! is aligned with e_z):                                                    !
    !                                                                          !
    ! 1. Apply R^{a}_{z}(gamma) to rotate the vector anticlockwise about the   !
    !    (fixed) e_z axis (0 <= gamma < 2*pi)                                  !
    !   [ This results in a vector in the xz plane, with a negative projection !
    !   on the (fixed) e_x axis ]                                              !
    !                                                                          !
    ! 2. Apply R^{a}_{y}(beta) to rotate the vector (now in the xz-plane)      !
    !    anticlockwise about the (fixed) e_y axis (0 <= beta <= pi)            !
    !    [ This rotates the vector in the xz plane into alignment with the     !
    !    (fixed) e_z axis ]                                                    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! beta_ang   out        Euler rotation angle about fixed y-axis (step 2)   !
    ! gamma_ang  out        Euler rotation angle about fixed z-axis (step 1)   !
    ! z_vec      in         Unit vector along z-direction of rotated coord.    !
    !                       system expressed in fixed coord. system of         !
    !                       rotation                                           !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 11/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use constants, only: PI
    use geometry, only: POINT
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = &
         "sph_harm_rot_get_active_Euler_ang_z2z"

    ! Arguments
    real(kind=DP), intent(out) :: beta_ang  ! Euler angle beta (about e_y)
    real(kind=DP), intent(out) :: gamma_ang ! Euler angle alpha (about e_z)
    type(POINT), intent(in)    :: z_vec     ! unit vector along z-direction of
                                            ! rotated coord system, expressed
                                            ! in coord system of rotation

    ! Local variables
    real(kind=DP) :: z_vec_xy_mag ! Magnitude of projection of z_vec on fixed
                                  ! xy-plane
    real(kind=DP) :: gamma_half   ! gamma angle in half range [0,pi]

    call utils_trace_in(myself)

    ! Compute gamma angle: +ve rotation about fixed z-axis, e_z, which takes
    ! the vector z_vec into the fixed xy-plane (with negative projection on
    ! the fixed e_x axis).
    !
    ! arccos has a range of [0, pi], but gamma should be in the range [0,2*pi),
    ! so we use the projection of z_vec on the fixed y-axis, e_y to determine
    ! whether gamma is in [0,pi] or (pi,2*pi), i.e.
    !
    !   gamma_half = arccos(R_xy.(-e_x)/|(R_xy|)
    !
    !                | gamma_half        ,  R_xy.e_y >= 0
    !        gamma = |
    !                | 2*pi - gamma_half ,  R_xy.e_y <  0
    ! where R = z_vec, R_xy is the projection of z_vec onto the fixed xy-plane.
    if (z_vec%x == 0.0_DP .and. z_vec%y == 0.0_DP) then
       ! Special case: if z_vec%x == z_vec%y == 0, then vector points along
       ! the fixed z-axis and no rotation is necessary. The gamma rotation is
       ! arbitrary, and we choose a value of 0 radians. If z_vec%x or z_vec%y
       ! are even slightly different to 0.0_DP, then a gamma rotation is needed.
       gamma_ang = 0.0_DP
    else
       ! Only evaluate gamma_half after checking for special case to avoid
       ! division by zero when z_vec_xy_mag == 0.0_DP
       z_vec_xy_mag = sqrt(z_vec%x**2+z_vec%y**2)
       gamma_half   = acos(-z_vec%x/z_vec_xy_mag)
       if (z_vec%y >= 0.0_DP) then
          ! gamma in [0,pi]
          gamma_ang = gamma_half
       else
          ! gamma in (pi,2*pi)
          gamma_ang = 2.0_DP*PI - gamma_half
       end if
    end if
    ! +ve/anticlockwise active rotation with 0 <= gamma <= 2*pi

    ! Compute beta angle: +ve rotation about fixed y-axis, e_y, which takes
    ! the vector aligned in the fixed xy-plane (with negative projection on the
    ! fixed e_x axis) from the gamma rotation into the fixed z-axis, e_z. Use
    !   beta = arccos( R.e_z / |R| )
    ! with R = z_vec, |R| = 1, i.e.
    !   beta = arccos( (z_vec)_z )
    ! where (z_vec)_z is the z-component of z_vec in the coordinate system of
    ! rotation.
    beta_ang = acos( z_vec%z )
    ! +ve/anticlockwise active rotation with 0 <= beta <= pi

    call utils_trace_out(myself)

  end subroutine sph_harm_rot_get_active_Euler_ang_z2z

  ! Return Wigner D-matrix for a given l and active Euler rotation (D)
  subroutine sph_harm_rot_get_Wigner_Dmat(Dmat,l,alpha_ang,beta_ang,gamma_ang)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Return Wigner-D matrix for rotation of CSHs (or equivalently, expressing !
    ! one set of CSHs in terms of another with a different axis of             !
    ! quantization)                                                            !
    !                                                                          !
    ! The Wigner D-matrix evaluated by this routine is as defined in           !
    ! appendix B of                                                            !
    !     M.A. Morrison and G.A. Parker, Aust. J. Phys. 40, 465 (1987)         !
    ! (MP1987), i.e. Eq. B6 from MP1987                                        !
    !   D_{mu,m}^{l}(alpha,beta,gamma)                                         !
    !     = exp(-i mu alpha) d_{mu,m}^{l}(beta) exp(-i m gamma)                !
    !                                                                          !
    ! with the d_{mu,m}^{l} the so-called "small-d" matrix.                    !
    !                                                                          !
    ! The small-d matrix is defined in terms of Jacobi polynomials (Eqs. B9    !
    ! in MP1987) --- see the internal_Wigner_smalld routine for details.       !
    !                                                                          !
    ! The full D-matrix is defined in terms of a rotation operator for         !
    ! active Euler rotations about fixed axes, i.e. Eq. 30 from MP1987:        !
    !   R^{a}(alpha,beta,gamma)                                                !
    !     = R^{a}_{z}(alpha)*R^{a}_{y}(beta)*R^{a}_{z}(gamma)                  !
    ! which corresponds to the following series of rotations about the fixed   !
    ! e_y and e_z axes:                                                        !
    !   1. Rotate system by gamma about e_z (0 <= gamma <= 2*pi)               !
    !   2. Rotate system by beta about e_y  (0 <= beta <= pi)                  !
    !   3. Rotate system by alpha about e_z (0 <= alpha <= 2*pi)               !
    !                                                                          !
    ! The returned D-matrix has elements D_{mu,m}^{l}(alpha,beta,gamma) with   !
    ! mu being the row index and m being the column index. Since Fortran uses  !
    ! a column-major ordering, this means that the fastest changing index of   !
    ! Dmat is mu.                                                              !
    !                                                                          !
    ! Dmat has intent inout to enable the D-matrix to be placed into a larger  !
    ! array without overwriting that array.                                    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! Dmat       inout      Wigner D-matrix (complex) with dimensions of       !
    !                       (2l+1) x (2l+1)                                    !
    ! l          in         Angular momentum of set of CSHs being rotated      !
    ! alpha_ang  in         Euler rotation angle about fixed z-axis (step 3)   !
    ! beta_ang   in         Euler rotation angle about fixed y-axis (step 2)   !
    ! gamma_ang  in         Euler rotation angle about fixed z-axis (step 1)   !
    !                       (angles are in radians)                            !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 10/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use constants, only: cmplx_0, cmplx_i
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "sph_harm_rot_get_Wigner_Dmat"

    ! Arguments
    complex(kind=DP), intent(inout) :: Dmat(:,:) ! Wigner D-matrix
    integer, intent(in)             :: l         ! angmom of SHs to rotate
    real(kind=DP), intent(in)       :: alpha_ang ! Euler angle alpha (about e_z)
    real(kind=DP), intent(in)       :: beta_ang  ! Euler angle beta (about e_y)
    real(kind=DP), intent(in)       :: gamma_ang ! Euler angle gamma (about e_z)

    ! Local variables
    integer :: irow, icol
    integer :: mu ! row index
    integer :: m  ! col index
    complex(kind=DP) :: exp_imu_alpha
    complex(kind=DP) :: exp_im_gamma
    real(kind=DP) :: smalld

    call utils_trace_in(myself)

    Dmat(:,:) = cmplx_0
    do m = -l, l
       do mu = -l, l
          irow = 1+l+mu
          icol = 1+l+m
          exp_imu_alpha = exp(-cmplx_i*mu*alpha_ang)
          exp_im_gamma = exp(-cmplx_i*m*gamma_ang)
          smalld = internal_Wigner_smalld(mu,m,l,beta_ang)
          Dmat(irow,icol) = exp_imu_alpha * smalld * exp_im_gamma
       end do
    end do

    call utils_trace_out(myself)

    contains

      elemental function internal_Wigner_smalld(mu,m,l,beta) result(smalld)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_get_Wigner_Dmat                              !
        ! Function:  Return an element of the Wigner small-d matrix, evaluated !
        !            in terms of Jacobi polynomials, as defined in Eq. B9 of   !
        !            Morrison and Parker Aust. J. Phys. 40, 465 (1987).        !
        !                                                                      !
        ! The Wigner small-d matrix elements have the form                     !
        !   d^{l}_{mu m}(beta) = f1(l,m,mu) * f2(beta,m,mu) *                  !
        !       P^(m-mu,m+mu)_{l-m}(cos(beta))                                 !
        ! with                                                                 !
        !   f1 = ( [ (l-m)!*(l+m)! ] / [ (l+mu)!*(l-mu)! ] )^{1/2}             !
        !   f2 = cos(beta/2)^{m+mu} * cos(beta/2)^{m-mu}                       !
        ! and where                                                            !
        !   P^(m-mu,m+mu)_{l-m}(cos(beta))                                     !
        ! is a Jacobi polynomial.                                              !
        !                                                                      !
        ! Note that in the context of the full Wigner D-matrix, defined using  !
        ! the convention of active rotation about fixed axes:                  !
        ! * l is the total angular momentum quantum number of the set of       !
        !   angular momentum eigenstates undergoing rotation;                  !
        ! * mu, m are the eigenvalues of the z-projection of the total         !
        !   angular momentum operator, J_z (with z the axis of quantization    !
        !   of the angular momentum eigenstates prior to application of the    !
        !   rotation operator);                                                !
        ! * beta is the active Euler rotation angle about the fixed y-axis     !
        !   of rotation.                                                       !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! mu          in             mu index of small-d matrix                !
        ! m           in             m index of small-d matrix                 !
        ! l           in             total angular momentum                    !
        ! beta        in             angle of rotation (radians)               !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 10/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!

        use constants, only: PI
        use services, only: services_factorial

        implicit none

        ! Parameters
        !character(len=*), parameter :: myself = &
        !     "sph_harm_rot_get_Wigner_Dmat::internal_Wigner_smalld"

        ! Arguments
        integer, intent(in) :: mu
        integer, intent(in) :: m
        integer, intent(in) :: l
        real(kind=DP), intent(in) :: beta

        ! Result
        real(kind=DP) :: smalld

        ! Local variables
        integer :: mu1, m1
        real(kind=DP) :: beta1
        real(kind=DP) :: f1, f2

        ! To avoid issues with division by zero, we need to deal with the following cases:
        ! 1. cos(beta/2) = 0 and m+mu < 0
        !    * cos(beta/2) = 0 when beta = pi
        !    * m+mu < 0 when m < -mu
        !    => Use Eq. B11 of Morrison and Parker, i.e.
        !         d^l_mu,m(beta) = d^l_-mu,-m(-beta),
        !       since -(m+mu) > 0 if m+mu < 0 and cos(-beta/2)^(-(m+mu)) is no longer a
        !       a division by zero.
        !
        ! 2. sin(beta/2) = 0 and m-mu < 0
        !    * sin(beta/2) = 0 when beta = 0
        !    * m-mu < 0 when m > mu
        !    => Use Eq. B11 of Morrison and Parker, i.e.
        !         d^l_mu,m(beta) = d^l_m,mu(-beta)
        !       since mu-m > 0 if m-mu < 0 and sin(-beta/2)^(mu-m) is no longer a division
        !       by zero
        !
        if (beta > PI/2.0_DP .and. m+mu < 0) then
           ! 1. cos(beta/2) = 0 and m+mu < 0
           beta1 = -beta
           mu1 = -mu
           m1 = -m
        else if (beta <= PI/2.0 .and. m-mu < 0) then
           ! 2. sin(beta/2) = 0 and m-mu < 0
           beta1 = -beta
           mu1 = m
           m1 = mu
        else
           ! No change from input required
           beta1 = beta
           mu1 = mu
           m1 = m
        end if

        f1 = sqrt(&
             real(services_factorial(l-m1)*services_factorial(l+m1),kind=DP)/&
             real(services_factorial(l+mu1)*services_factorial(l-mu1),kind=DP)&
             )
        f2 = cos(beta1/2.0_DP)**(m1+mu1)*sin(beta1/2.0_DP)**(m1-mu1)

        smalld = f1*f2*jacobi_polynomial_realx(l-m1,m1-mu1,m1+mu1,cos(beta1))

      end function internal_Wigner_smalld

  end subroutine sph_harm_rot_get_Wigner_Dmat

  subroutine sph_harm_rot_get_RSH_CSH_Umat(Umat,l)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Return unitary transformation matrix between RSHs and CSHs (U).          !
    !                                                                          !
    ! The RSHs are defined in terms of the CSHs as                             !
    !                                                                          !
    !             | Y_{l,0}                                    m = 0           !
    !   Z_{l,m} = | i/sqrt(2) (Y_{l,-|m|} - (-1)^m Y_{l,|m|})  m < 0           !
    !             | 1/sqrt(2) (Y_{l,-m}   + (-1)^m Y_{l,m})    m > 0           !
    !                                                                          !
    ! which corresponds to a 2l+1 x 2l+1 unitary transformation matrix for     !
    ! each l. See the module header for a definition of this matrix, which is  !
    ! based on Eq. B8 of                                                       !
    !   J. Fernandez Rico, R. Lopez, and G. Ramirez,                           !
    !                                         J. Chem. Phys. 91, 4204 (1989)   !
    !                                                                          !
    ! Umat has intent inout to enable the U-matrix to be placed into a larger  !
    ! array without overwriting that array.                                    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! Umat   inout          RSH-CSH unitary transformation matrix (complex)    !
    !                       with dimensions of (2l+1) x (2l+1)                 !
    ! l      in             Angular momentum of SHs being transformed          !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 07/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use constants, only: cmplx_0, cmplx_1, cmplx_i, SQRT_TWO
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "sph_harm_rot_get_RSH_CSH_Umat"
    real(kind=DP), parameter    :: OO_SQRT_TWO = 1.0_DP/SQRT_TWO

    ! Arguments
    complex(kind=DP), intent(inout) :: Umat(:,:) ! RSH-CSH transformation matrix
    integer, intent(in)             :: l         ! angmom of SHs to rotate

    ! Local variables
    integer :: nsh ! side length of Umat (2l+1)
    integer :: ish ! counter

    call utils_trace_in(myself)

    nsh = 2*l+1

    ! Matrix is zero everywhere except diagonal and antidiagonal
    Umat(1:nsh,1:nsh) = cmplx_0

    ! Central element
    Umat(1+l,1+l) = cmplx_1

    ! Fill diagonal and antidiagonal
    do ish = 1, l
       ! Top left diagonal
       !    U(ish,ish) = 2^(-1/2)*i, for ish = 1,l
       Umat(ish,ish)         = OO_SQRT_TWO * cmplx_i
       ! Bottom left antidiagonal
       !    U(nsh-ish+1,ish) = 2^(-1/2)*(-1)^(ish+l)*i, for ish = 1, l
       Umat(nsh-ish+1,ish)   = OO_SQRT_TWO * (-1.0_DP)**(ish+l) * cmplx_i
       ! Top right antidiagonal
       !    U(l+1-ish,l+1+ish) = 2^(-1/2)*1, for ish = 1, l
       Umat(l+1-ish,l+1+ish) = OO_SQRT_TWO * cmplx_1
       ! Bottom right diagonal
       !    U(l+1+ish,l+1+ish) = 2^(-1/2)*(-1)^(ish)
       Umat(l+1+ish,l+1+ish) = OO_SQRT_TWO * (-1.0_DP)**ish * cmplx_1
    end do

    call utils_trace_out(myself)

  end subroutine sph_harm_rot_get_RSH_CSH_Umat

  ! Return RSH rotation matrix for a given l and Euler rotation (R = U* D U)
  subroutine sph_harm_rot_get_RSH_Rmat(Rmat,l,alpha_ang,beta_ang,gamma_ang)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Return matrix for rotation of RSHs (or equivalently, expressing one set  !
    ! of RSHs in terms of another with a different axis of quantization).      !
    !                                                                          !
    ! The RSH rotation matrix evaluated by this routine is constructed using   !
    ! Wigner D-matrices and RSH <-> CSH transformation matrices using          !
    !   R^{l} = (U^{l})* D^{l} U^{l}                                           !
    ! i.e. Eq. B9 from                                                         !
    !   J. Fernandez Rico, R. Lopez, and G. Ramirez,                           !
    !                                         J. Chem. Phys. 91, 4204 (1989)   !
    ! where R^{l}, D^{l} and U^{l} are, respectively, the RSH rotation matrix, !
    ! the Wigner D-matrix and RSH <-> CSH transformation matrix for total      !
    ! angular momentum l ("*" denotes the complex conjugate operation).        !
    !                                                                          !
    ! The RSH rotation matrix is constructed from a Wigner D-matrix defined as !
    ! in sph_harm_rot_get_Wigner_Dmat, i.e. in terms of an active Euler        !
    ! rotation about fixed axes. The Euler rotations occur in the following    !
    ! order (about the fixed e_z and e_y axes):                                !
    !   1. Rotate system by gamma about e_z (0 <= gamma <= 2*pi)               !
    !   2. Rotate system by beta about e_y  (0 <= beta <= pi)                  !
    !   3. Rotate system by alpha about e_z (0 <= alpha <= 2*pi)               !
    !                                                                          !
    ! The returned matrix has elements R_{mu,m}^{l}(alpha,beta,gamma) with     !
    ! mu being the row index and m being the column index. Since Fortran uses  !
    ! a column-major ordering, this means that the fastest changing index of   !
    ! Rmat is mu.                                                              !
    !                                                                          !
    ! Rmat has intent inout to enable the R-matrix to be placed into a larger  !
    ! array without overwriting that array.                                    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! Rmat       inout      RSH rotation matrix (real) with dimensions of      !
    !                       (2l+1) x (2l+1)                                    !
    ! l          in         Angular momentum of set of RSHs being rotated      !
    ! alpha_ang  in         Euler rotation angle about fixed z-axis (step 3)   !
    ! beta_ang   in         Euler rotation angle about fixed y-axis (step 2)   !
    ! gamma_ang  in         Euler rotation angle about fixed z-axis (step 1)   !
    !                       (angles are in radians)                            !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 11/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "sph_harm_rot_get_RSH_Rmat"

    ! Arguments
    real(kind=DP), intent(inout) :: Rmat(:,:) ! RSH rotation matrix
    integer, intent(in)          :: l         ! angmom of SHs to rotate
    real(kind=DP), intent(in)    :: alpha_ang ! Euler angle alpha (about e_z)
    real(kind=DP), intent(in)    :: beta_ang  ! Euler angle beta (about e_y)
    real(kind=DP), intent(in)    :: gamma_ang ! Euler angle gamma (about e_z)

    ! Local variables
    complex(kind=DP), dimension(2*l+1,2*l+1) :: Dmat, Umat

    call utils_trace_in(myself)

    ! Get RSH <-> CSH transformation matrix for ang. mom. l
    call sph_harm_rot_get_RSH_CSH_Umat(Umat,l)
    ! Get Wigner D-matrix for ang. mom. l and Euler angles
    call sph_harm_rot_get_Wigner_Dmat(Dmat,l,alpha_ang,beta_ang,gamma_ang)

    !---------------------------------------------------------------------------
    ! @optimize Matrix multiplication
    ! The implementation of this routine currently uses the Fortran intrinsic
    ! function matmul to perform matrix multiplication. This must involve
    ! a copy, since the result of the function must be distinct from the
    ! arguments.
    !
    ! The matrices multiplied by this routine will have dimensions
    !   2*l+1 x 2*l+1
    ! probably not exceeding l values of ~6, i.e.
    !   13 x 13
    !
    ! It may be more efficient to use a hand-coded matrix multiplication or
    ! a routine from a linear-algebra library.
    !---------------------------------------------------------------------------

    ! Evaluate RSH rotation matrix as conj_trans(Umat).Dmat.Umat
    ! as described in Eq. B9 of J. Fernandez Rico, R. Lopez, and G. Ramirez,
    ! J. Chem. Phys. 91, 4204 (1989).
    Rmat = real(matmul(matmul(conjg(transpose(Umat)),Dmat),Umat),kind=DP)

    ! NOTE:
    ! For speed, we don't check whether the result of the matrix multiplication
    ! is real (it should be). We simply take the real part of the result and
    ! assume that the input produces a sane result.

    call utils_trace_out(myself)

  end subroutine sph_harm_rot_get_RSH_Rmat

  ! Apply RSH rotation matrix for l between 0 and lmax to atomblock of a
  ! matrix with matrix elements containing RSHs (M' = R* M R)
  subroutine sph_harm_rot_apply_RSH_Rmat_to_atomblock(Mmat, &
      lmax, num_q_per_l, alpha_ang, beta_ang, gamma_ang)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Apply a per-atomblock RSH rotation matrix to both sides of an input      !
    ! matrix Mmat and return the result in Mmat. The RSH rotation matrix is    !
    ! constructed in terms of an active Euler rotation about fixed axes (see   !
    ! the documentation for sph_harm_rot_get_RSH_Rmat).                        !
    !                                                                          !
    ! The RSH rotation matrix is applied to both sides of Mmat, i.e.           !
    !   M' = transpose(R) M R                                                  !
    ! where M is Mmat on entry (before rotation) and M' is Mmat on exit        !
    ! (after rotation).                                                        !
    !                                                                          !
    ! The per-atomblock RSH rotation matrix is constructed as a block-diagonal !
    ! matrix of RSH rotation matrices evaluated by sph_harm_rot_get_RSH_Rmat   !
    ! for l-values between 0 and lmax. As described in the documentation for   !
    ! sph_harm_rot_get_RSH_Rmat, the RSH rotation matrix for angular momentum  !
    ! l is defined in terms of Euler angles alpha, beta, and gamma, which      !
    ! describe an active rotation about fixed axes. The per-atomblock matrix   !
    ! is defined in terms of the same set of Euler angles, but for a range of  !
    ! angular momentum values from 0 to lmax.                                  !
    !                                                                          !
    ! For each l-value there are num_q_per_l(l) repetitions of the same        !
    ! rotation matrix along the diagonal of the per-atomblock RSH rotation     !
    ! matrix, e.g. for lmax = 1, num_q_per_l = [ 2, 2 ], the block diagonal    !
    ! matrix would have the following form:                                    !
    !                                                                          !
    !    l     0       1                                                       !
    !  l q   1   2   1   2                                                     !
    !      +---+---+---+---+                                                   !
    !    1 | R | 0 | 0 | 0 |   R with l = 0 is a 1 by 1 RSH rotation matrix    !
    !      |l=0|   |   |   |                                                   !
    !  0   +---+---+---+---+                                                   !
    !    2 | 0 | R | 0 | 0 |                                                   !
    !      |   |l=0|   |   |                                                   !
    !      +---+---+---+---+                                                   !
    !    1 | 0 | 0 | R | 0 |   R with l = 1 is a 3 by 3 RSH rotation matrix    !
    !      |   |   |l=1|   |                                                   !
    !  1   +---+---+---+---+                                                   !
    !    2 | 0 | 0 | 0 | R |                                                   !
    !      |   |   |   |l=1|                                                   !
    !      +---+---+---+---+                                                   !
    !                                                                          !
    ! The RSH rotation matrix therefore has dimensions nmat by nmat, with      !
    ! nmat = ( \sum_{l=0}^{lmax} num_q_per_l(l) * (2*l+1) )                    !
    !                                                                          !
    ! In a typical use case, M is an atomblock of matrix elements in a basis   !
    ! of spherical waves (products of RSHs and spherical Bessel functions).    !
    ! The row and column indices of matrix elements M_{ij} are compound indexes!
    ! of the total angular momentum, l, angular momentum component, m, and     !
    ! spherical Bessel index, q. These indices are incremented as with the     !
    ! per-atomblock RSH rotation matrix, i.e. (fastest changing index first)   !
    !                                                                          !
    ! Rows:    m_i, q_i, l_i                                                   !
    ! Columns: m_j, q_j, l_j                                                   !
    !                                                                          !
    ! Applying the per-atomblock RSH rotation matrix R to both sides of an M   !
    ! matrix constructed in this way is equivalent to applying RSH rotation    !
    ! matrices for each l-value to the appropriate subblocks of the M matrix,  !
    ! e.g. with the lmax = 1, num_q_per_l = [ 2, 2 ] example above, we have    !
    !                                                                          !
    !      / R0          \^T / M00 M00 M01 M01 \ / R0          \               !
    !      |    R0       |   | M00 M00 M01 M01 | |    R0       |               !
    !      |       R1    |   | M10 M10 M11 M11 | |       R1    |               !
    !      \          R1 /   \ M10 M10 M11 M11 / \          R1 /               !
    !                                                                          !
    !         / R0^T M00 R0   R0^T M00 R0   R0^T M01 R1   R0^T M01 R1 \        !
    !      =  | R0^T M00 R0   R0^T M00 R0   R0^T M01 R1   R0^T M01 R1 |        !
    !         | R1^T M10 R0   R1^T M10 R0   R1^T M11 R1   R1^T M11 R1 |        !
    !         \ R1^T M10 R0   R1^T M10 R0   R1^T M11 R1   R1^T M11 R1 /        !
    !                                                                          !
    ! R0 and R1 are the RSH rotation matrices for l = 0 and 1. M00, M01, M10,  !
    ! M11 are sub-blocks of the M matrix with different angular momentum       !
    ! spherical waves:                                                         !
    ! * M00 is a 1 by 1 matrix of elements M00_{ij} with l_i = l_j = 0         !
    ! * M10 is a 3 by 1 matrix of elements M10_{ij} with l_i = 1, l_j = 0      !
    !   with the 3 m_i components of l_i = 1 running with the compound row     !
    !   index i (see above).                                                   !
    ! * M11 is a 3 by 3 matrix of elements M11_{ij} with l_i = l_j = 1         !
    !   with the 3 m_i components of l_i = 1 and 3 m_j components of l_j       !
    !   running with the compound row index i and column index j,              !
    !   respectively.                                                          !
    !                                                                          !
    ! In summary, the effect is to apply the appropriate RSH rotation matrices !
    ! for each total angular momentum value to the appropriate regions of an   !
    ! atomblock matrix in the basis of RSHs / spherical waves.                 !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name>   <in/out/inout> <arg descrption>                                 !
    ! Mmat         inout    atomblock in RSH / spherical wave (SW) basis       !
    !                       on entry -- matrix to be rotated                   !
    !                       on exit  -- rotated matrix                         !
    ! lmax         in       Maximum angular momentum of RSHs in atomblock      !
    ! num_q_per_l  in       Number of identical submatrices per l (this is the !
    !                       number of spherical Bessels for a SW basis)        !
    ! alpha_ang    in       Euler rotation angle about fixed z-axis (step 3)   !
    ! beta_ang     in       Euler rotation angle about fixed y-axis (step 2)   !
    ! gamma_ang    in       Euler rotation angle about fixed z-axis (step 1)   !
    !                       (angles are in radians)                            !
    !                                                                          !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 11/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_trace_in, utils_trace_out

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = &
         "sph_harm_rot_apply_RSH_Rmat_to_atomblock"

    ! Arguments
    real(kind=DP), intent(inout) :: Mmat(:,:) ! on entry -- matrix to be rotated
                                              ! on exit  -- rotated matrix
    integer, intent(in)        :: lmax        ! max angmom of SHs in Mmat
    integer, intent(in)        :: num_q_per_l(0:lmax) ! num submatrices per l
    real(kind=DP), intent(in)  :: alpha_ang   ! Euler angle alpha (about e_z)
    real(kind=DP), intent(in)  :: beta_ang    ! Euler angle beta (about e_y)
    real(kind=DP), intent(in)  :: gamma_ang   ! Euler angle gamma (about e_z)

    ! Local variables
    real(kind=DP), allocatable :: Rmat(:,:)
    integer :: nmat
    integer :: ierr_alloc
    integer :: il

    call utils_trace_in(myself)

    ! Compute dimension of (square) Mmat and Rmat, i.e.
    ! nmat = ( \sum_{l=0}^{lmax} num_q_per_l(l) * (2*l+1) )                    !
    nmat = sum( [ ( num_q_per_l(il)*(2*il+1), il = 0, lmax ) ] )

    ! Allocate Rmat (per-atomblock RSH rotation matrix)
    allocate(Rmat(nmat,nmat),stat=ierr_alloc)
    call utils_alloc_check(myself,"Rmat",ierr_alloc)

    ! Build Rmat
    call internal_build_RSH_Rmat_atomblock(Rmat,lmax,num_q_per_l,&
         alpha_ang,beta_ang,gamma_ang)

    ! TEMPORARY CODE --> probably better to use LAPACK routines
    Mmat(1:nmat,1:nmat) = matmul(&
      matmul(transpose(Rmat(1:nmat,1:nmat)),Mmat(1:nmat,1:nmat)),&
      Rmat(1:nmat,1:nmat))
    !/TEMPORARY CODE

    ! Deallocate Rmat
    deallocate(Rmat,stat=ierr_alloc)
    call utils_dealloc_check(myself,"Rmat",ierr_alloc)

    call utils_trace_out(myself)

    contains

      subroutine internal_build_RSH_Rmat_atomblock(Rmat_atomblock,lmax,&
           num_q_per_l,alpha_ang,beta_ang,gamma_ang)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:     sph_harm_rot_apply_RSH_Rmat_to_atomblock                   !
        ! Function: Build the per-atomblock RSH rotation matrix for rotating   !
        !           the atomblock matrix Mmat.                                 !
        !                                                                      !
        ! The per-atomblock rotation matrix is constructed as a block-diagonal !
        ! matrix of per-l RSH rotation matrices for l-values between 0 and     !
        ! lmax. See the documentation for the host routine for further         !
        ! details.                                                             !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>    <in/out/inout> <arg descrption>                            !
        ! Rmat_atomblock out                                                   !
        ! lmax           in   Maximum angular momentum of RSHs in atomblock    !
        ! num_q_per_l    in   Number of identical submatrices per l            !
        ! alpha_ang      in   Euler rotation angle about fixed z-axis (step 3) !
        ! beta_ang       in   Euler rotation angle about fixed y-axis (step 2) !
        ! gamma_ang      in   Euler rotation angle about fixed z-axis (step 1) !
        !                     (angles are in radians)                          !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 11/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        implicit none

        ! Parameters
        !character(len=*), parameter :: myself = &
        !     "sph_harm_rot_apply_RSH_Rmat_to_atomblock::&
        !     &internal_build_RSH_Rmat_atomblock"

        ! Arguments
        real(kind=DP), intent(out) :: Rmat_atomblock(:,:)
        integer, intent(in)        :: lmax      ! max angmom of SHs in Mmat
        integer, intent(in)        :: num_q_per_l(0:lmax) ! submatrices per l
        real(kind=DP), intent(in)  :: alpha_ang ! Euler angle alpha (about e_z)
        real(kind=DP), intent(in)  :: beta_ang  ! Euler angle beta (about e_y)
        real(kind=DP), intent(in)  :: gamma_ang ! Euler angle gamma (about e_z)

        ! Local variables
        integer :: il
        integer :: iq
        integer :: istart_l, iend_l
        integer :: istart_q, iend_q
        integer :: nsh

        ! Zero Rmat_atomblock
        Rmat_atomblock(:,:) = 0.0_DP

        !----------------------------------------------------------------------
        ! Non-contiguous array arguments
        ! To avoid allocating extra memory for a work array to hold the per-l
        ! R-matrices, the result of sph_harm_rot_get_RSH_Rmat is directly
        ! placed in a section of Rmat_atomblock. Rmat_atomblock is thus passed
        ! to sph_harm_rot_get_RSH_Rmat as a non-contiguous array.
        !
        ! If we used a contiguous work-array to evaluate the per-l R-matrices,
        ! we would still need to copy these into the per-atomblock R-matrix in
        ! a non-contiguous way. This suggests that placing the result of
        ! sph_harm_rot_get_RSH_Rmat directly into a non-contiguous section of
        ! Rmat_atomblock should not be more expensive than placing the result
        ! into a contiguous work array and then copying this into the per-
        ! atomblock array.
        !----------------------------------------------------------------------
        istart_l = 1
        do il = 0, lmax
           nsh = 2*il + 1
           iend_l = istart_l + nsh - 1
           ! Get per-l RSH R-matrix
           call sph_harm_rot_get_RSH_Rmat(Rmat_atomblock(istart_l:iend_l,&
                istart_l:iend_l),il,alpha_ang,beta_ang,gamma_ang)
           ! The per-l RSH R-matrix is the same for all q-values, so simply
           ! copy this from the already-evaluated array
           istart_q = istart_l + nsh
           iend_q   = istart_q + nsh - 1
           do iq = 2, num_q_per_l(il)
              Rmat_atomblock(istart_q:iend_q,istart_q:iend_q) = &
                   Rmat_atomblock(istart_l:iend_l,istart_l:iend_l)
              istart_q = istart_q + nsh
              iend_q   = istart_q + nsh - 1
           end do
           ! Next l-value starts immediately after last q value for previous l
           istart_l = istart_q
        end do

      end subroutine internal_build_RSH_Rmat_atomblock

  end subroutine sph_harm_rot_apply_RSH_Rmat_to_atomblock

  subroutine sph_harm_rot_run_unit_tests(ierror)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Run unit tests for sph_harm_rotation module and output results to        !
    ! stdout.                                                                  !
    !                                                                          !
    ! This routine should be called by all MPI ranks. All MPI ranks will       !
    ! run the tests, but only the root MPI rank will output to stdout. The     !
    ! consistency of results across MPI ranks is checked and if discrepancies  !
    ! are found, this is reported.                                             !
    !                                                                          !
    ! An error code is returned. If this is zero, all tests passed and no      !
    ! discrepancies in outcomes were detected between MPI ranks. A non-zero    !
    ! error code indicates test failure and/or discrepancies between ranks.    !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS                                                              !
    ! <name>    <in/out/inout> <arg descrption>                                !
    ! ierror    out            error code (== 0 success, != 0 failure)         !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 07/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!
    use comms, only: pub_total_num_procs, comms_allgather
    use rundat, only: pub_print_qc
    use utils, only: utils_qc_print, utils_trace_in, utils_trace_out, &
         utils_alloc_check, utils_dealloc_check

    implicit none

    ! Parameters
    character(len=*), parameter :: myself = "sph_harm_rot_run_unit_tests"
    integer, parameter :: IFAIL = 1
    integer, parameter :: IPASS = 2
    integer, parameter :: ISKIP = 3
    integer, parameter :: TEST_UNDEFINED = 1
    integer, parameter :: TEST_PASSED = 0
    integer, parameter :: TEST_FAILED = -1
    integer, parameter :: TEST_SKIPPED = -2
    integer, parameter :: ERROR_MPI_DISAGREE = 100
    integer, parameter :: ERROR_FAILED_TEST  = 1000

    ! Arguments
    integer, intent(out) :: ierror ! 0 for success

    ! Local variables
    ! - General
    integer :: itest   ! count total number of tests
    integer :: outcome ! integer representation of test outcome

    integer :: my_results(3) ! test results
      ! results(IFAIL) counts fails
      ! results(IPASS) counts passes
      ! results(ISKIP) counts skips
    integer,allocatable :: all_results(:) ! results gathered for all MPI ranks
    integer :: irank           ! counter for MPI ranks
    integer :: iblock          ! index of block in all_results array
    logical :: all_ranks_agree ! do results for all MPI ranks agree?
    integer :: ierr_alloc      ! holds status returned by allocate/deallocate

    call utils_trace_in(myself)

    ! Assume the best, unless otherwise informed
    ierror = 0

    ! Initialize counters
    itest = 0
    my_results(1:3) = 0

    ! START OF UNIT TESTS
    call internal_test_sph_harm_rot_get_active_Euler_ang_z2z(outcome)
    call internal_update_results(itest,my_results,outcome)

    call internal_test_sph_harm_rot_get_RSH_CSH_Umat(outcome)
    call internal_update_results(itest,my_results,outcome)

    call internal_test_sph_harm_rot_get_Wigner_Dmat(outcome)
    call internal_update_results(itest,my_results,outcome)

    call internal_test_sph_harm_rot_get_RSH_Rmat(outcome)
    call internal_update_results(itest,my_results,outcome)

    call internal_test_sph_harm_rot_apply_RSH_Rmat_to_atomblock(outcome)
    call internal_update_results(itest,my_results,outcome)

    call internal_test_jacobi_polynomial_realx(outcome)
    call internal_update_results(itest,my_results,outcome)
    ! END OF UNIT TESTS

    ! All MPI processes run the unit tests
    ! At the end of the tests, we check the equality of the counters on each
    ! process by first gathering to the root MPI rank, then checking equality
    if (pub_on_root) then
       ! Allocate full-size array on root MPI rank (to be gathered to)
       allocate(all_results(pub_total_num_procs*3),stat=ierr_alloc)
       call utils_alloc_check(myself,"all_results",ierr_alloc)
    else
       ! Allocate single element array on non-root MPI ranks to satisfy
       ! comms_allgather_integer_1, which does not give the dummy argument
       ! i_array_dest the allocatable attribute (so the actual argument should
       ! be allocated).
       !
       ! Note: I have chosen not to use a zero-length array here since within
       ! comms_allgather_integer_1, the call to MPI_GATHER passes
       ! i_array_dest(1) as an actual argument. It is unclear to me whether
       ! it is legal Fortran to pass a zero-length array in this way.
       allocate(all_results(1),stat=ierr_alloc)
       call utils_alloc_check(myself,"all_results",ierr_alloc)
    end if
    call comms_allgather(all_results,my_results,3,gather_not_allgather=.true.)

    if (pub_on_root) then
       ! Check reports from all MPI ranks agree with root MPI rank
       all_ranks_agree = .true.
       do irank = 0, pub_total_num_procs-1
          iblock = irank * 3
          if (.not.all(all_results(iblock+1:iblock+3)==my_results(1:3))) then
             all_ranks_agree = .false.
          end if
       end do
       write(stdout,"(a)")      module_short_code//": [Unit test report]"

       if (all_ranks_agree) then
          ! All ranks agree, output a summary of results
          write(stdout,"(a)")      module_short_code//": All MPI ranks agree"
          write(stdout,"(a,i0,a)") module_short_code//": ", itest, " tests run"
          write(stdout,"(a,i0,a)") module_short_code//": * ", &
               my_results(IPASS), " passed"
          write(stdout,"(a,i0,a)") module_short_code//": * ", &
               my_results(IFAIL), " failed"
          write(stdout,"(a,i0,a)") module_short_code//": * ", &
               my_results(ISKIP), " skipped"
          ! If running QC tests, allow unit tests (pass, fail, skip) to be
          ! checked against benchmarks
          if (pub_print_qc) then
             call utils_qc_print(module_short_code//"_unit_test_pass",&
                  my_results(IPASS))
             call utils_qc_print(module_short_code//"_unit_test_fail",&
                  my_results(IFAIL))
             call utils_qc_print(module_short_code//"_unit_test_skip",&
                  my_results(ISKIP))
          end if
       else
          ! Disagreement between ranks, quit without summarizing results
          write(stdout,"(a)") module_short_code//": Disagreement &
               &between MPI ranks -- quitting tests with error."
          ierror = ierror + ERROR_MPI_DISAGREE
       end if
       if (sph_harm_rot_debug) then
          ! In debug mode, output results for all ranks
          write(stdout,"(a,i0,a)") module_short_code//" DEBUG: &
               &[per-MPI rank results]"
          do irank = 0, pub_total_num_procs-1
             iblock = irank * 3
             write(stdout,"(a,i0)")   module_short_code//" DEBUG: &
                  &Results for MPI rank ",irank
             write(stdout,"(a,i0,a)") module_short_code//" DEBUG: * ", &
                  all_results(iblock+IPASS), " passed"
             write(stdout,"(a,i0,a)") module_short_code//" DEBUG: * ", &
                  all_results(iblock+IFAIL), " failed"
             write(stdout,"(a,i0,a)") module_short_code//" DEBUG: * ", &
                  all_results(iblock+ISKIP), " skipped"
          end do
       end if
    end if

    ! Non-zero error code in event of failed tests (in root MPI rank)
    if (my_results(IFAIL) > 0) ierror = ierror + ERROR_FAILED_TEST

    ! Allocated on all MPI ranks
    deallocate(all_results,stat=ierr_alloc)
    call utils_dealloc_check(myself,"all_results",ierr_alloc)

    call utils_trace_out(myself)

    contains

      subroutine internal_update_results(itest,my_results,outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Update test and results counters based on outcome of      !
        !            a unit test                                               !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! itest       inout          test counter (increment by 1)             !
        ! my_results  inout          array counting of passes, fails and skips !
        ! outcome     inout          integer outcome code of preceding test    !
        !                            set to TEST_UNDEFINED on exit             !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 07/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use utils, only: utils_abort

        implicit none

        ! Parameters
        character(len=*), parameter :: myself = &
             "sph_harm_rot_run_unit_tests::internal_update_results"

        ! Arguments
        integer,intent(inout) :: itest
        integer,intent(inout) :: my_results(3)
        integer,intent(inout) :: outcome

        ! Increment test counter
        itest = itest + 1

        select case (outcome)

        case (TEST_PASSED)
           my_results(IPASS) = my_results(IPASS) + 1

        case (TEST_FAILED)
           my_results(IFAIL) = my_results(IFAIL) + 1

        case (TEST_SKIPPED)
           my_results(ISKIP) = my_results(ISKIP) + 1

        case default
           call utils_abort("Error in "//myself//": &
                &outcome code not recognized.")
        end select

        ! Reset outcome to undefined state
        outcome = TEST_UNDEFINED

      end subroutine internal_update_results

      subroutine internal_check_unitary(outcome,mat_dem,nmat,mat_label,err_thr)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Check whether a matrix is unitary and report via outcome  !
        ! (and to stdout, when in debug mode)                                  !
        !                                                                      !
        ! ## NOTE                                                              !
        ! This routine supports complex and real matrices mat_dem, and will    !
        ! automatically detect the type of the matrix based on the icmplx      !
        ! component of the DEM derived type.                                   !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        ! mat_dem     inout          matrix to check (as an instance of DEM)   !
        ! nmat        in             side length of (square) matrix            !
        ! mat_label   in             one-character label for matrix, e.g. "U"  !
        ! err_thr     in             maximum threshold for numerical error     !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 09/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use constants, only: cmplx_1
        use dense, only: DEM, dense_axpy, dense_create, dense_copy, &
             dense_destroy, dense_determinant, dense_invert, dense_norm, &
             dense_product, dense_put_element
        implicit none

        ! Parameters
        !character(len=*), parameter :: myself = &
        !     "sph_harm_rot_run_unit_tests::internal_check_unitary"

        ! Arguments
        integer,intent(inout) :: outcome
        type(DEM),intent(in) :: mat_dem
        integer,intent(in) :: nmat
        character(len=1),intent(in) :: mat_label
        real(kind=DP),intent(in) :: err_thr

        ! Local variables
        type(DEM) :: mat_dem_inv, mat_dem_conj, mat_dem_prod
        type(DEM) :: Imat_dem
        integer :: icol ! column counter for putting columns in DEM instance
        complex(kind=DP) :: mat_det_cmplx, mat_det_conj_cmplx, &
             mat_det_prod_cmplx ! complex determinants
        real(kind=DP)    :: mat_det_real, mat_det_conj_real, &
             mat_det_prod_real ! real determinants
        real(kind=DP)    :: norm
        logical          :: det_test_passed
        character(len=9) :: transpose_label

        ! Allocate workspace of inverse of mat_dem
        call dense_create(mat_dem_inv,nmat,nmat,iscmplx=mat_dem%iscmplx)
        ! Allocate workspace of conjugate transpose of mat_dem
        call dense_create(mat_dem_conj,nmat,nmat,iscmplx=mat_dem%iscmplx)
        ! Allocate workspace of product of mat with conj. trans. of mat_dem
        call dense_create(mat_dem_prod,nmat,nmat,iscmplx=mat_dem%iscmplx)
        ! Allocate workspace for identity matrix
        call dense_create(Imat_dem,nmat,nmat,iscmplx=mat_dem%iscmplx)

        if (mat_dem%iscmplx) then
           transpose_label = "conjtrans"
        else
           transpose_label = "transpose"
        end if

        if (mat_dem%iscmplx) then
           ! Build complex identity matrix Imat
           do icol = 1, nmat
              ! Contents are always zeroed by dense_create, so just place
              ! cmplx_1 along diagonal
              call dense_put_element(cmplx_1,Imat_dem,icol,icol)
           end do
        else
           ! Build real identity matrix Imat
           do icol = 1, nmat
              ! Contents are always zeroed by dense_create, so just place
              ! 1.0_DP along diagonal
              call dense_put_element(1.0_DP,Imat_dem,icol,icol)
           end do
        end if

        ! Get inverse of mat_dem
        ! * For unitary matrix M^-1 == conjugate_transpose(M)
        call dense_copy(mat_dem_inv,mat_dem)
        call dense_invert(mat_dem_inv) ! in place inversion

        ! Get (conjugate) transpose of mat_dem
        ! --> workaround for lack of conjugate transpose routine by
        !     using dense_product for form the conjugate transpose
        ! --> note that opB="C" results in conjugate transpose of
        !     argument B for complex matrix and transpose of argument B
        !     for real matrix (no need to use opB="T")
        call dense_product(mat_dem_conj,Imat_dem,mat_dem,opB="C")

        ! [ Check matrix is unitary ]
        ! (i)  M.conjugate_transpose(M) = I (complex M, I)
        !      M.transpose(M)           = I (real M, I)
        call dense_product(mat_dem_prod,mat_dem,mat_dem_conj)
        ! ... check for equality with identity
        if (mat_dem%iscmplx) then
           call dense_axpy(mat_dem_prod,Imat_dem,-cmplx_1)
        else
           call dense_axpy(mat_dem_prod,Imat_dem,-1.0_DP)
        end if
        ! ... using maximum column sum ("1") norm
        norm = dense_norm("1",mat_dem_prod)
        ! ... and comparing to err_thr
        if (sph_harm_rot_debug.and.pub_on_root) then
           write(stdout,'(a,es16.8)') module_short_code//" DEBUG: &
                &* 1-norm of "//mat_label//&
                " "//transpose_label//"("//mat_label//") - I = ",norm
        end if
        if (norm > err_thr) then
           outcome = TEST_FAILED
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   &Test failure!"
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   & --> 1-norm of "//mat_label//" "//transpose_label//"("&
                   //mat_label//") - I not within error threshold."
           end if
        end if

        ! (ii) det(M)det(conjugate(M)) = 1 (complex M)
        !      det(M)det(transpose(M)) = 1 (real M)
        det_test_passed = .true.
        if (mat_dem%iscmplx) then
           ! Check complex determinants
           call dense_determinant(amat=mat_dem,det=mat_det_cmplx)
           call dense_determinant(amat=mat_dem_conj,det=mat_det_conj_cmplx)
           ! ... check for equality with real 1
           mat_det_prod_cmplx = mat_det_cmplx*mat_det_conj_cmplx
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a,2es16.8)') module_short_code//" DEBUG: &
                   &* det("//mat_label//&
                   ") det("//transpose_label//"("//mat_label//&
                   "))     = ",mat_det_prod_cmplx
           end if
           if ( real(mat_det_prod_cmplx,kind=DP)-1.0_DP > err_thr.or.&
                aimag(mat_det_prod_cmplx) > err_thr ) det_test_passed = .false.
        else
           ! Check real determinants
           call dense_determinant(amat=mat_dem,det=mat_det_real)
           call dense_determinant(amat=mat_dem_conj,det=mat_det_conj_real)
           ! ... check for equality with real 1
           mat_det_prod_real = mat_det_real*mat_det_conj_real
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a,2es16.8)') module_short_code//" DEBUG: &
                   &* det("//mat_label//&
                   ") det("//transpose_label//"("//mat_label//&
                   "))     = ",mat_det_prod_real
           end if
           if ( mat_det_prod_real-1.0_DP > err_thr ) det_test_passed = .false.
        end if

        if (.not.det_test_passed) then
           outcome = TEST_FAILED
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   &Test failure!"
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   & --> det("//mat_label//") det("//transpose_label//"("&
                   //mat_label//")) - 1 not within &
                   &error threshold."
           end if
        end if

        ! Deallocate workspace of inverse of mat_dem
        call dense_destroy(mat_dem_inv)
        ! Deallocate workspace of conjugate transpose of mat_dem
        call dense_destroy(mat_dem_conj)
        ! Deallocate workspace of product of Umat with conj. trans. of mat_dem
        call dense_destroy(mat_dem_prod)
        ! Deallocate workspace for complex identity matrix
        call dense_destroy(Imat_dem)

      end subroutine internal_check_unitary

      subroutine internal_check_ref(outcome,mat_dem,mat_dem_ref,mat_label,&
           err_thr)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Check a matrix against reference data and report via      !
        ! outcome (and to stdout, when in debug mode)                          !
        !                                                                      !
        ! ## NOTE                                                              !
        ! This routine supports complex and real matrices mat_dem, and will    !
        ! automatically detect the type of the matrix based on the icmplx      !
        ! component of the DEM derived type.                                   !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        ! mat_dem     inout          matrix to check (as an instance of DEM)   !
        ! mat_dem_ref inout          reference data (as an instance of DEM)    !
        ! mat_label   in             one-character label for matrix, e.g. "U"  !
        ! err_thr     in             maximum threshold for numerical error     !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 09/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use constants, only: cmplx_1
        use dense, only: DEM, dense_axpy, dense_create, dense_copy, &
             dense_destroy, dense_norm
        use utils, only: utils_assert
        implicit none

        ! Parameters
        character(len=*), parameter :: myself = &
             "sph_harm_rot_run_unit_tests::internal_check_ref"

        ! Arguments
        integer,intent(inout) :: outcome
        type(DEM),intent(in) :: mat_dem
        type(DEM),intent(in) :: mat_dem_ref
        character(len=1),intent(in) :: mat_label
        real(kind=DP),intent(in) :: err_thr

        ! Local variables
        type(DEM) :: mat_dem_diff
        real(kind=DP) :: norm

        ! Check mat_dem and mat_dem_ref are the same type (real or complex)
        call utils_assert(mat_dem%iscmplx.eqv.mat_dem_ref%iscmplx, &
             "Error in "//myself//": mat_dem and mat_dem_ref are not of the &
             &same type (both must be real or both must be complex).")

        ! Create copy of mat_ref_dem (this will be overwritten by dense_axpy)
        ! as workspace for difference between computed and reference matrices
        call dense_create(mat_dem_diff,mat_dem_ref)
        call dense_copy(mat_dem_diff,mat_dem_ref)

        ! Check equality of computed mat_dem with reference mat_dem_ref by
        ! ... subtracting mat_dem from mat_dem_ref
        if (mat_dem%iscmplx) then
           call dense_axpy(mat_dem_diff,mat_dem,-cmplx_1)
        else
           call dense_axpy(mat_dem_diff,mat_dem,-1.0_DP)
        end if
        ! ... computing the maximum column sum ("1") norm
        norm = dense_norm("1",mat_dem_diff)
        ! ... and comparing to error threshold
        if (sph_harm_rot_debug.and.pub_on_root) then
           write(stdout,'(a,es16.8)') module_short_code//" DEBUG: &
                &* 1-norm of "//mat_label//&
                "ref - "//mat_label//"           = ",norm
        end if
        if (norm > err_thr) then
           outcome = TEST_FAILED
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   &Test failure!"
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   & --> 1-norm of "//mat_label//"ref - "//mat_label//&
                   " not within error threshold."
           end if
        end if

        ! Deallocate workspace for difference between computed and reference
        ! matrices
        call dense_destroy(mat_dem_diff)

      end subroutine internal_check_ref

      elemental function internal_isclose(A,Aref,relthr,absthr) result(isclose)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Check a for closeness in real values A and Aref, using    !
        !            an absolute and relative error threshold                  !
        !                                                                      !
        ! The function uses the test                                           !
        !   absolute(A - Aref) <= (absthr + relthr * absolute(Aref))           !
        ! to check closeness between real value A and a reference value Aref   !
        !                                                                      !
        ! This is the same expression used by numpy.isclose()                  !
        ! https://docs.scipy.org/doc/numpy-1.15.1/reference/                   !
        !         generated/numpy.isclose.html                                 !
        ! and has the effect that very small reference values Aref are         !
        ! checked against an error threshold close to absthr                   !
        ! while larger Aref values are checked against a threshold which       !
        ! scales with the size of Aref.                                        !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS AND RESULT                                               !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! A           in             real value to compare against reference   !
        ! Aref        in             reference value to compare A against      !
        ! relthr      in             relative threshold                        !
        ! absthr      in             absolute threshold                        !
        ! iclose      out (result)   true if test passes, false if test fails  !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 11/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        implicit none

        ! Arguments
        real(kind=DP), intent(in) :: A
        real(kind=DP), intent(in) :: Aref
        real(kind=DP), intent(in) :: relthr
        real(kind=DP), intent(in) :: absthr

        ! Result
        logical :: isclose

        isclose = ( abs(A-Aref) <= (absthr + relthr * abs(Aref)) )

      end function internal_isclose

      subroutine internal_test_sph_harm_rot_get_active_Euler_ang_z2z(outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Unit test for sph_harm_rot_get_active_Euler_ang_z2z       !
        !            routine                                                   !
        !----------------------------------------------------------------------!
        ! ## TESTS                                                             !
        ! * Test computed angles against some simple known input vectors       !
        ! * Test against some precomputed input/output combinations            !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 11/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use constants, only: PI
        use geometry, only: POINT, unit_vector
        implicit none

        ! Parameters
        character(len=*), parameter :: unit_test_routine = &
             "sph_harm_rot_get_active_Euler_ang_z2z"
        integer, parameter :: N_SIMPLE_REF = 9
        ! --> number of simple reference vectors
        real(kind=DP), parameter :: PI_O_2 = PI/2.0_DP
        real(kind=DP), parameter :: PI_O_4 = PI/4.0_DP
        real(kind=DP), parameter :: THREE_PI_O_2 = 3.0_DP*PI/2.0_DP
        real(kind=DP), parameter :: THREE_PI_O_4 = 3.0_DP*PI/4.0_DP
        real(kind=DP), parameter :: SEVEN_PI_O_4 = 7.0_DP*PI/4.0_DP
        real(kind=DP), parameter :: ONE_O_SQRT2  = 1.0_DP/sqrt(2.0_DP)
        real(kind=DP), parameter :: MINUS_ONE_O_SQRT2  = -ONE_O_SQRT2
        ! --> constants used for simple input vector tests
        integer, parameter :: N_ARB_REF = 4
        ! --> number of arbitrary input vectors
        real(kind=DP), parameter :: ERR_FACTOR = 1.0e2_DP
        real(kind=DP), parameter :: ERR_THRESHOLD = ERR_FACTOR * EPSILON(1.0_DP)
        ! --> rough error threshold for unit testing based on machine epsilon
        !     (to adjust, change ERR_FACTOR)

        ! Reference data (do not modify)
        ! Simple unnormalized input vectors for routine
        ! (these must be converted to unit vectors before being passed to
        ! sph_harm_rot_get_active_Euler_ang_z2z)
        real(kind=DP) :: Z_VEC_SIMPLE_UNNORM_REF(3,N_SIMPLE_REF)
        ! Euler angles corresponding to simple input vectors
        real(kind=DP) :: BETA_GAMMA_SIMPLE_REF(2,N_SIMPLE_REF)

        ! Vector points along fixed z-axis
        !   beta  = 0 deg = 0.0_DP
        !   gamma = 0 deg = 0.0_DP
        ! (special case, since gamma is arbitrary in this case)
        data Z_VEC_SIMPLE_UNNORM_REF(:,1) &
             / 0.0_DP, 0.0_DP, 1.0_DP /
        data BETA_GAMMA_SIMPLE_REF(:,1) &
             / 0.0_DP, 0.0_DP /

        ! Vector points along fixed negative of z-axis
        !   beta  = 180 deg = PI
        !   gamma = 0 deg   = 0.0_DP
        ! (special case, since gamma is arbitrary in this case)
        data Z_VEC_SIMPLE_UNNORM_REF(:,2) &
             / 0.0_DP, 0.0_DP, -1.0_DP /
        data BETA_GAMMA_SIMPLE_REF(:,2) &
             / PI, 0.0_DP /

        ! Vector points along fixed x-axis
        !   beta  = 90 deg  = PI/2
        !   gamma = 180 deg = PI
        data Z_VEC_SIMPLE_UNNORM_REF(:,3) &
             / 1.0_DP, 0.0_DP, 0.0_DP /
        data BETA_GAMMA_SIMPLE_REF(:,3) &
             / PI_O_2, PI /

        ! Vector points along negative of fixed x-axis
        !   beta  = 90 deg  = PI/2
        !   gamma = 0 deg   = 0.0_DP
        data Z_VEC_SIMPLE_UNNORM_REF(:,4) &
             / -1.0_DP, 0.0_DP, 0.0_DP /
        data BETA_GAMMA_SIMPLE_REF(:,4) &
             / PI_O_2, 0.0_DP /

        ! Vector points along fixed y-axis
        !   beta  = 90 deg  = PI/2
        !   gamma = 90 deg  = PI/2
        data Z_VEC_SIMPLE_UNNORM_REF(:,5) &
             / 0.0_DP, 1.0_DP, 0.0_DP /
        data BETA_GAMMA_SIMPLE_REF(:,5) &
             / PI_O_2, PI_O_2 /

        ! Vector points along negative of fixed y-axis
        !   beta  = 90 deg  = PI/2
        !   gamma = 270 deg = 3*PI/2
        data Z_VEC_SIMPLE_UNNORM_REF(:,6) &
             / 0.0_DP, -1.0_DP, 0.0_DP /
        data BETA_GAMMA_SIMPLE_REF(:,6) &
             / PI_O_2, THREE_PI_O_2 /

        ! Vector in xz plane, 45 deg from z axis and +ve projection on x-axis
        !   beta  = 45 deg  = PI/4
        !   gamma = 180 deg = PI
        data Z_VEC_SIMPLE_UNNORM_REF(:,7) &
             / 1.0_DP, 0.0_DP, 1.0_DP /
        data BETA_GAMMA_SIMPLE_REF(:,7) &
             / PI_O_4, PI /

        ! Vector (1/sqrt(2),1/sqrt(2),1)
        !   beta  = 45 deg  = PI/4
        !   gamma = 135 deg = 3*PI/4
        data Z_VEC_SIMPLE_UNNORM_REF(:,8) &
             / ONE_O_SQRT2, ONE_O_SQRT2, 1.0_DP /
        data BETA_GAMMA_SIMPLE_REF(:,8) &
             / PI_O_4, THREE_PI_O_4 /

        ! Vector (-1/sqrt(2),-1/sqrt(2),-1)
        !   beta  = 45 deg  = 3*PI/4
        !   gamma = 315 deg = 7*PI/4
        data Z_VEC_SIMPLE_UNNORM_REF(:,9) &
             / MINUS_ONE_O_SQRT2, MINUS_ONE_O_SQRT2, -1.0_DP /
        data BETA_GAMMA_SIMPLE_REF(:,9) &
             / THREE_PI_O_4, SEVEN_PI_O_4 /

        ! Arbitrary input vectors (with precomputed angles) for routine
        ! (these are already unit vectors, so no need to use unit_vector()
        ! before passing to sph_harm_rot_get_active_Euler_ang_z2z)
        real(kind=DP) :: Z_VEC_ARB_REF(3,N_ARB_REF)
        ! Euler angles corresponding to arbitrary input vectors
        real(kind=DP) :: BETA_GAMMA_ARB_REF(2,N_ARB_REF)

        ! Unit vector 0
        data Z_VEC_ARB_REF(:,1) &
             / 9.8574785881601901e-02_DP, &
             -3.8444715936400897e-01_DP, &
             -9.1786894121401696e-01_DP /
        data BETA_GAMMA_ARB_REF(:,1) &
             / 2.7334735456329469e+00_DP, &
             4.4613897322305691e+00_DP /

        ! Unit vector 1
        data Z_VEC_ARB_REF(:,2) &
             / -2.5402719632555198e-01_DP, &
             9.5277277727437804e-01_DP, &
             1.6641579976627000e-01_DP /
        data BETA_GAMMA_ARB_REF(:,2) &
             / 1.4036026675811963e+00_DP, &
             1.3102385638418190e+00_DP /

        ! Unit vector 2
        data Z_VEC_ARB_REF(:,3) &
             / -4.3347754445643601e-02_DP, &
             -7.4519080848054597e-01_DP, &
             -6.6544092986577696e-01_DP /
        data BETA_GAMMA_ARB_REF(:,3) &
             / 2.2988806914552127e+00_DP, &
             4.7704935092215219e+00_DP /

        ! Unit vector 3
        data Z_VEC_ARB_REF(:,4) &
             / 1.0337006774894200e-01_DP, &
             -7.6181814332273101e-01_DP, &
             -6.3949022322306504e-01_DP /
        data BETA_GAMMA_ARB_REF(:,4) &
             / 2.2646313276086558e+00_DP, &
             4.5775240001637494e+00_DP /

        ! Arguments
        integer,intent(out) :: outcome

        ! Local variables
        integer       :: iref
        type(POINT)   :: z_vec
        real(kind=DP) :: beta_ang_ref, gamma_ang_ref
        real(kind=DP) :: beta_ang, gamma_ang
        real(kind=DP) :: beta_abs_diff, gamma_abs_diff

        if (pub_on_root.and.sph_harm_rot_debug) then
           write(stdout,'(a)') module_short_code//&
                " DEBUG: Unit test for "//unit_test_routine
           write(stdout,'(a,es16.8)') module_short_code//&
                " DEBUG: Error threshold = ", ERR_THRESHOLD
        end if

        ! Assume everything is fine, unless it's not
        outcome = TEST_PASSED

        ! [ Check against some simple known input vectors ]
        if (sph_harm_rot_debug.and.pub_on_root) then
           write(stdout,'(a,i0)') module_short_code//" DEBUG: &
                &Testing computed beta and gamma angles against simple vectors"
        end if
        do iref = 1, N_SIMPLE_REF
           ! Get simple reference vector
           z_vec%x = Z_VEC_SIMPLE_UNNORM_REF(1,iref)
           z_vec%y = Z_VEC_SIMPLE_UNNORM_REF(2,iref)
           z_vec%z = Z_VEC_SIMPLE_UNNORM_REF(3,iref)
           ! Convert to a unit vector
           z_vec = unit_vector(z_vec)
           ! Get corresponding simple reference angles
           beta_ang_ref  = BETA_GAMMA_SIMPLE_REF(1,iref)
           gamma_ang_ref = BETA_GAMMA_SIMPLE_REF(2,iref)
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a,3f10.4)') module_short_code//" DEBUG: &
                   &Input z_vec            = ",z_vec%x,z_vec%y,z_vec%z
              write(stdout,'(a,2es16.8)') module_short_code//" DEBUG: &
                   &Expected beta, gamma   = ",beta_ang_ref, gamma_ang_ref
           end if
           call sph_harm_rot_get_active_Euler_ang_z2z(beta_ang,gamma_ang,z_vec)
           beta_abs_diff = abs(beta_ang-beta_ang_ref)
           gamma_abs_diff = abs(gamma_ang-gamma_ang_ref)
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a,2es16.8)') module_short_code//" DEBUG: &
                   &Computed beta, gamma   = ",beta_ang, gamma_ang
              write(stdout,'(a,2es16.8)') module_short_code//" DEBUG: &
                   &Abs. diff. beta, gamma = ",beta_abs_diff, gamma_abs_diff
           end if
           if (beta_abs_diff > ERR_THRESHOLD ) then
              outcome = TEST_FAILED
              if (sph_harm_rot_debug.and.pub_on_root) then
                 write(stdout,'(a)') module_short_code//" DEBUG: &
                      &Test failure!"
                 write(stdout,'(a)') module_short_code//" DEBUG: &
                      & --> Computed beta not within error &
                      &threshold of expected value."
              end if
           end if
           if (gamma_abs_diff > ERR_THRESHOLD ) then
              outcome = TEST_FAILED
              if (sph_harm_rot_debug.and.pub_on_root) then
                 write(stdout,'(a)') module_short_code//" DEBUG: &
                      &Test failure!"
                 write(stdout,'(a)') module_short_code//" DEBUG: &
                      & --> Computed gamma not within error &
                      &threshold of expected value."
              end if
           end if
        end do

        ! [ Check against a set of arbitrary input vectors ]
        if (sph_harm_rot_debug.and.pub_on_root) then
           write(stdout,'(a,i0)') module_short_code//" DEBUG: &
                &Testing computed beta and gamma angles against arbitrary &
                &vectors"
        end if
        do iref = 1, N_ARB_REF
           ! Get arbitrary unit vector
           z_vec%x = Z_VEC_ARB_REF(1,iref)
           z_vec%y = Z_VEC_ARB_REF(2,iref)
           z_vec%z = Z_VEC_ARB_REF(3,iref)
           ! Get corresponding reference angles
           beta_ang_ref  = BETA_GAMMA_ARB_REF(1,iref)
           gamma_ang_ref = BETA_GAMMA_ARB_REF(2,iref)
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a,3f10.4)') module_short_code//" DEBUG: &
                   &Input z_vec            = ",z_vec%x,z_vec%y,z_vec%z
              write(stdout,'(a,2es16.8)') module_short_code//" DEBUG: &
                   &Expected beta, gamma   = ",beta_ang_ref, gamma_ang_ref
           end if
           call sph_harm_rot_get_active_Euler_ang_z2z(beta_ang,gamma_ang,z_vec)
           beta_abs_diff = abs(beta_ang-beta_ang_ref)
           gamma_abs_diff = abs(gamma_ang-gamma_ang_ref)
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a,2es16.8)') module_short_code//" DEBUG: &
                   &Computed beta, gamma   = ",beta_ang, gamma_ang
              write(stdout,'(a,2es16.8)') module_short_code//" DEBUG: &
                   &Abs. diff. beta, gamma = ",beta_abs_diff, gamma_abs_diff
           end if
           if (beta_abs_diff > ERR_THRESHOLD ) then
              outcome = TEST_FAILED
              if (sph_harm_rot_debug.and.pub_on_root) then
                 write(stdout,'(a)') module_short_code//" DEBUG: &
                      &Test failure!"
                 write(stdout,'(a)') module_short_code//" DEBUG: &
                      & --> Computed beta not within error &
                      &threshold of expected value."
              end if
           end if
           if (gamma_abs_diff > ERR_THRESHOLD ) then
              outcome = TEST_FAILED
              if (sph_harm_rot_debug.and.pub_on_root) then
                 write(stdout,'(a)') module_short_code//" DEBUG: &
                      &Test failure!"
                 write(stdout,'(a)') module_short_code//" DEBUG: &
                      & --> Computed gamma not within error &
                      &threshold of expected value."
              end if
           end if
        end do

      end subroutine internal_test_sph_harm_rot_get_active_Euler_ang_z2z

      subroutine internal_test_sph_harm_rot_get_RSH_CSH_Umat(outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Unit test for sph_harm_rot_get_RSH_CSH_Umat routine       !
        !----------------------------------------------------------------------!
        ! ## TESTS                                                             !
        ! * Test unitarity of matrix up to LMAX                                !
        ! * Check matrix agrees with precomputed reference up to LMAX_REF      !
        !                                                                      !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 07/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use constants, only: garbage_complex
        use dense, only: DEM, dense_create, dense_copy, &
             dense_destroy, dense_put_col, dense_put_element
        use utils, only: utils_alloc_check, utils_dealloc_check
        implicit none

        ! Parameters
        character(len=*), parameter :: unit_test_routine = &
             "sph_harm_rot_get_RSH_CSH_Umat"
        integer, parameter :: LMAX = 3 ! max angular momentum to test
        integer, parameter :: LMAX_REF = 2
        ! --> max angular momentum in reference data
        integer, parameter :: NSH_TOT_MAX = LMAX*(LMAX+2)+2
        ! --> total number of angular momentum components from l = 0 .. lmax
        real(kind=DP), parameter :: ERR_FACTOR = 1.0e2_DP
        real(kind=DP), parameter :: ERR_THRESHOLD = ERR_FACTOR * EPSILON(1.0_DP)
        ! --> rough error threshold for unit testing based on machine epsilon
        !     (to adjust, change ERR_FACTOR)

        ! Reference data (do not modify)
        ! Store matrices sequentially in 1-D array up to LMAX_REF
        ! l = 0 (1x1 matrix, 1 element)
        ! l = 1 (3x3 matrix, 9 elements)
        ! l = 2 (5x5 matrix, 25 elements)
        complex(kind=DP) :: UMAT_REF(1+9+25)
        ! Pre-computed matrices from Sage Math implementation of U matrix
        ! l = 0  (1 element,   1:1+1-1)
        data UMAT_REF(1:1) &
             / (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP) /
        ! l = 1  (9 elements,  2:2+9-1)
        data UMAT_REF(2:2+9-1) &
             / (0.0000000000000000e+00_DP, 7.0710678118654746e-01_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 7.0710678118654746e-01_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (7.0710678118654746e-01_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-7.0710678118654746e-01_DP, 0.0000000000000000e+00_DP) /

        ! l = 2 (25 elements, 10:10+25-1)
        data UMAT_REF(11:11+25-1) &
             / (0.0000000000000000e+00_DP, 7.0710678118654746e-01_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, -7.0710678118654746e-01_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 7.0710678118654746e-01_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 7.0710678118654746e-01_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (7.0710678118654746e-01_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-7.0710678118654746e-01_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (7.0710678118654746e-01_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (7.0710678118654746e-01_DP, 0.0000000000000000e+00_DP) /

        ! Arguments
        integer,intent(out) :: outcome

        ! Local variables
        complex(kind=DP),allocatable :: Umat(:,:)
        type(DEM) :: Umat_dem, Umat_dem_ref
        integer :: il ! angular momentum counter
        integer :: nsh ! number of spherical harmonics (2*l + 1)
        integer :: icol ! column counter for putting columns in DEM instance
        integer :: istart_ref ! start index in UMAT_REF array
        integer :: ierr_alloc ! holds status returned by allocate/deallocate

        if (pub_on_root.and.sph_harm_rot_debug) then
           write(stdout,'(a)') module_short_code//&
                " DEBUG: Unit test for "//unit_test_routine
           write(stdout,'(a,es16.8)') module_short_code//&
                " DEBUG: Error threshold = ", ERR_THRESHOLD
        end if

        ! Assume everything is fine, unless it's not
        outcome = TEST_PASSED

        ! Counter for position in 1-D reference data array
        istart_ref = 1

        do il = 0, LMAX

           ! -------------------------------------------------------------------
           ! # Note on use of local allocatable array and DEM object
           ! We are storing the U matrix in a local complex array and an
           ! instance of DEM.
           !
           ! This is wasteful, since we have two copies of U, however, it
           ! serves a purpose, as it allows us
           !
           ! (i)  to ensure that Umat is always evaluated locally on an MPI
           !      process by sph_harm_rot_get_RSH_CSH_Umat; and
           ! (ii) to use the dense module, which provides an agnostic interface
           !      to LAPACK and ScaLAPACK.
           !
           ! In this case, we are not concerned about whether Umat is distributed
           ! as a DEM, because we expect all MPI ranks to hold the same Umat,
           ! so it is not important whether operations are performed locally or
           ! in a distributed fashion.
           !
           ! If we had a different Umat (e.g. different l values) on each MPI
           ! rank, then we would want to avoid using the dense module to do
           ! matrix operations, since when using ScaLAPACK, we would be
           ! attempting to distribute different data into a single distributed
           ! object.
           ! -------------------------------------------------------------------
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a,i0)') module_short_code//" DEBUG: &
                   &Testing U matrix with l = ",il
           end if

           nsh = 2*il + 1

           ! Allocate local array for U matrix
           allocate(Umat(nsh,nsh),stat=ierr_alloc)
           call utils_alloc_check(myself,"Umat",ierr_alloc)

           ! Initialize to garbage
           Umat(:,:) = garbage_complex

           ! Evaluate U matrix
           call sph_harm_rot_get_RSH_CSH_Umat(Umat,il)

           ! Allocate workspace for dense matrix operations (DEM instance)
           ! (this may be distributed, if using ScaLAPACK)
           call dense_create(Umat_dem,nsh,nsh,iscmplx=.true.)

           ! Put local Umat (identical on each MPI rank) into DEM
           do icol = 1, nsh
              call dense_put_col(Umat(:,icol),Umat_dem,icol)
           end do

           ! [ Check matrix is unitary ]
           call internal_check_unitary(outcome,Umat_dem,nsh,"U",ERR_THRESHOLD)

           ! [ Check matrix against known results ]
           ! Only perform check if reference data is available (il <= LMAX_REF)
           if (il <= LMAX_REF) then
              ! Allocate workspace for reference Umat with total ang mom il
              call dense_create(Umat_dem_ref,nsh,nsh,iscmplx=.true.)
              ! Check against pre-computed reference data

              ! Put local UMAT_REF (identical on each MPI rank) into DEM
              do icol = 1, nsh
                 ! This takes the 1-D array containing the refence data and
                 ! extracts columns to put in the Umat_dem_ref DEM instance
                 call dense_put_col(&
                      UMAT_REF(istart_ref+(icol-1)*nsh:istart_ref+icol*nsh-1),&
                      Umat_dem_ref,icol)
              end do

              call internal_check_ref(outcome,Umat_dem,Umat_dem_ref,"U",&
                   ERR_THRESHOLD)

              ! Deallocate workspace for reference Umat with total ang mom il
              call dense_destroy(Umat_dem_ref)

              ! Update start point in reference data array
              ! (Start point in reference data array is the sum over the squares
              ! of angular momentum components for for l = 0 .. il-1)
              istart_ref = istart_ref + nsh*nsh
           end if

           ! Deallocate workspace for dense matrix operations (DEM instance)
           call dense_destroy(Umat_dem)

           ! Deallocate local array for U matrix
           deallocate(Umat,stat=ierr_alloc)
           call utils_dealloc_check(myself,"Umat",ierr_alloc)
        end do

      end subroutine internal_test_sph_harm_rot_get_RSH_CSH_Umat

      subroutine internal_test_sph_harm_rot_get_Wigner_Dmat(outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Unit test for sph_harm_rot_get_Wigner_Dmat routine        !
        !----------------------------------------------------------------------!
        ! ## TESTS                                                             !
        ! * Test unitarity of matrix up to LMAX                                !
        ! * Check matrix agrees with precomputed reference up to LMAX_REF      !
        !                                                                      !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 09/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use constants, only: garbage_complex
        use dense, only: DEM, dense_create, dense_copy, &
             dense_destroy, dense_put_col, dense_put_element
        use utils, only: utils_alloc_check, utils_dealloc_check
        implicit none

        ! Parameters
        character(len=*), parameter :: unit_test_routine = &
             "sph_harm_rot_get_Wigner_Dmat"
        integer, parameter :: NROT_REF = 5
        integer, parameter :: LMAX_REF = 2
        ! --> max angular momentum in reference data
        real(kind=DP), parameter :: ERR_FACTOR = 1.0e2_DP
        real(kind=DP), parameter :: ERR_THRESHOLD = ERR_FACTOR * EPSILON(1.0_DP)
        ! --> rough error threshold for unit testing based on machine epsilon
        !     (to adjust, change ERR_FACTOR)

        ! Reference data (do not modify)
        ! Store matrices sequentially in 1-D array up to LMAX_REF
        ! We have 5 rotation sets, each containing Wigner D-matrices up
        ! to LMAX_REF (rotation sets are numbered 0 to 4)
        ! Each rotation set has:
        !   l = 0 (1x1 matrix, 1 element)
        !   l = 1 (3x3 matrix, 9 elements)
        !   l = 2 (5x5 matrix, 25 elements)
        ! Rotation angles for rotation sets (in radians)
        !   0: alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        !   1: alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        !   2: alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        !   3: alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        !   4: alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        ! Start indices for each reference matrix:
        !   Rotation set   l   Start index   End index
        ! +--------------+---+-------------+-----------+
        !   0              0   1             1
        !   0              1   2             10
        !   0              2   11            35
        !   1              0   36            36
        !   1              1   37            45
        !   1              2   46            70
        !   2              0   71            71
        !   2              1   72            80
        !   2              2   81            105
        !   3              0   106           106
        !   3              1   107           115
        !   3              2   116           140
        !   4              0   141           141
        !   4              1   142           150
        !   4              2   151           175

        ! Euler angles for each rotation set
        !   0: alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        !   1: alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        !   2: alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        !   3: alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        !   4: alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        real(kind=DP) :: DMAT_EULER_ANGLES_REF(5*3)

        ! Rotation set 0
        data DMAT_EULER_ANGLES_REF(1:3) &
             / 0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP /

        ! Rotation set 1
        data DMAT_EULER_ANGLES_REF(4:6) &
             / 5.4641481231202804e-01_DP, &
             8.1804498047957241e-01_DP, &
             6.1422738037708546e+00_DP /

        ! Rotation set 2
        data DMAT_EULER_ANGLES_REF(7:9) &
             / 3.2078413032109494e+00_DP, &
             1.9773138676292171e-01_DP, &
             2.2398638613670325e+00_DP /

        ! Rotation set 3
        data DMAT_EULER_ANGLES_REF(10:12) &
             / 1.0649532499842701e-01_DP, &
             1.2464827205876983e+00_DP, &
             2.8973157492962169e+00_DP /

        ! Rotation set 4
        data DMAT_EULER_ANGLES_REF(13:15) &
             / 1.9651774322985169e+00_DP, &
             5.1291187749676292e-01_DP, &
             5.0980521145061619e-01_DP /

        ! Pre-computed matrices from Sage Math implementation of D matrix
        ! Rotation set 0
        complex(kind=DP) :: DMAT_REF(5*(1+9+25))
        ! Dmat for il = 0, alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        data DMAT_REF(1:1) &
             / (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP) /

        ! Dmat for il = 1, alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        data DMAT_REF(2:10) &
             / (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP) /

        ! Dmat for il = 2, alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        data DMAT_REF(11:35) &
             / (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (-0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (0.0000000000000000e+00_DP, 0.0000000000000000e+00_DP), &
             (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP) /

        ! Rotation set 1
        ! Dmat for il = 0, alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        data DMAT_REF(36:36) &
             / (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP) /

        ! Dmat for il = 1, alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        data DMAT_REF(37:45) &
             / (7.7355601412781783e-01_DP, 3.3208408609360124e-01_DP), &
             (-5.1093915274122848e-01_DP, 7.2477546342137725e-02_DP), &
             (1.2226086784397174e-01_DP, -1.0035796149981616e-01_DP), &
             (4.4091297775202531e-01_DP, 2.6815584755425925e-01_DP), &
             (6.8364930698757953e-01_DP, 0.0000000000000000e+00_DP), &
             (-4.4091297775202537e-01_DP, 2.6815584755425931e-01_DP), &
             (1.2226086784397174e-01_DP, 1.0035796149981616e-01_DP), &
             (5.1093915274122848e-01_DP, 7.2477546342137725e-02_DP), &
             (7.7355601412781783e-01_DP, -3.3208408609360124e-01_DP) /

        ! Dmat for il = 2, alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        data DMAT_REF(46:70) &
             / (4.8810906675669430e-01_DP, 5.1377128398769045e-01_DP), &
             (-5.9299204215248991e-01_DP, -1.6066784076753038e-01_DP), &
             (3.1329686988951516e-01_DP, -9.0708563846620488e-02_DP), &
             (-7.8056337306768961e-02_DP, 8.5047946271468222e-02_DP), &
             (4.8759993695625375e-03_DP, -2.4539702936038851e-02_DP), &
             (3.5641094497562914e-01_DP, 5.0042477498087556e-01_DP), &
             (2.8412605182129613e-01_DP, 1.2197402454538717e-01_DP), &
             (-6.0501088563559025e-01_DP, 8.5821774013391208e-02_DP), &
             (2.8942798299043432e-01_DP, -2.3757726315988714e-01_DP), &
             (-3.8176453397227125e-02_DP, 1.0894266090357800e-01_DP), &
             (1.5002719942002787e-01_DP, 2.8961148401240477e-01_DP), &
             (5.2209181803117888e-01_DP, 3.1752745106095204e-01_DP), &
             (2.0106456241689677e-01_DP, 0.0000000000000000e+00_DP), &
             (-5.2209181803117877e-01_DP, 3.1752745106095204e-01_DP), &
             (1.5002719942002787e-01_DP, -2.8961148401240477e-01_DP), &
             (3.8176453397227125e-02_DP, 1.0894266090357800e-01_DP), &
             (2.8942798299043432e-01_DP, 2.3757726315988711e-01_DP), &
             (6.0501088563559025e-01_DP, 8.5821774013391208e-02_DP), &
             (2.8412605182129619e-01_DP, -1.2197402454538717e-01_DP), &
             (-3.5641094497562908e-01_DP, 5.0042477498087556e-01_DP), &
             (4.8759993695625375e-03_DP, 2.4539702936038851e-02_DP), &
             (7.8056337306768961e-02_DP, 8.5047946271468222e-02_DP), &
             (3.1329686988951511e-01_DP, 9.0708563846620488e-02_DP), &
             (5.9299204215248991e-01_DP, -1.6066784076753035e-01_DP), &
             (4.8810906675669430e-01_DP, -5.1377128398769045e-01_DP) /

        ! Rotation set 2
        ! Dmat for il = 0, alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        data DMAT_REF(71:71) &
             / (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP) /

        ! Dmat for il = 1, alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        data DMAT_REF(72:80) &
             / (6.6428611635506241e-01_DP, -7.3439337795903381e-01_DP), &
             (8.6158293063794467e-02_DP, -1.0895940488587019e-01_DP), &
             (5.5237418205800019e-03_DP, -8.0253926619393585e-03_DP), &
             (-1.3860318148067619e-01_DP, -9.1957305410949788e-03_DP), &
             (9.8051475932109300e-01_DP, 0.0000000000000000e+00_DP), &
             (1.3860318148067619e-01_DP, -9.1957305410949788e-03_DP), &
             (5.5237418205800019e-03_DP, 8.0253926619393585e-03_DP), &
             (-8.6158293063794467e-02_DP, -1.0895940488587019e-01_DP), &
             (6.6428611635506241e-01_DP, 7.3439337795903381e-01_DP) /

        ! Dmat for il = 2, alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        data DMAT_REF(81:105) &
             / (-9.8057589207988807e-02_DP, -9.7569464984256404e-01_DP), &
             (-3.2223290922449459e-02_DP, -1.9184415056332191e-01_DP), &
             (-5.4487680087543493e-03_DP, -2.2995212858201594e-02_DP), &
             (-5.6360065047790075e-04_DP, -1.8290275087770463e-03_DP), &
             (-3.3895203677885620e-05_DP, -8.8660394146660590e-05_DP), &
             (-1.3976029559396672e-01_DP, 1.3531289058593932e-01_DP), &
             (6.3839856664139261e-01_DP, -7.0577371451397897e-01_DP), &
             (1.4632274807031578e-01_DP, -1.8504590775959101e-01_DP), &
             (1.6355962584095714e-02_DP, -2.3763424570696828e-02_DP), &
             (9.7836544636681127e-04_DP, -1.6449280757643331e-03_DP), &
             (2.3424813897749606e-02_DP, 3.1220155449679815e-03_DP), &
             (-2.3538997448016599e-01_DP, -1.5617121874627478e-02_DP), &
             (9.4211378986975103e-01_DP, 0.0000000000000000e+00_DP), &
             (2.3538997448016596e-01_DP, -1.5617121874627474e-02_DP), &
             (2.3424813897749610e-02_DP, -3.1220155449679815e-03_DP), &
             (-9.7836544636681106e-04_DP, -1.6449280757643331e-03_DP), &
             (1.6355962584095714e-02_DP, 2.3763424570696828e-02_DP), &
             (-1.4632274807031576e-01_DP, -1.8504590775959101e-01_DP), &
             (6.3839856664139261e-01_DP, 7.0577371451397897e-01_DP), &
             (1.3976029559396672e-01_DP, 1.3531289058593932e-01_DP), &
             (-3.3895203677885620e-05_DP, 8.8660394146660590e-05_DP), &
             (5.6360065047790075e-04_DP, -1.8290275087770463e-03_DP), &
             (-5.4487680087543502e-03_DP, 2.2995212858201598e-02_DP), &
             (3.2223290922449459e-02_DP, -1.9184415056332191e-01_DP), &
             (-9.8057589207988807e-02_DP, 9.7569464984256404e-01_DP) /

        ! Rotation set 3
        ! Dmat for il = 0, alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        data DMAT_REF(106:106) &
             / (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP) /

        ! Dmat for il = 1, alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        data DMAT_REF(107:115) &
             / (-6.5308074464460597e-01_DP, 9.0556255790004966e-02_DP), &
             (6.5034714907320024e-01_DP, -1.6210196121135434e-01_DP), &
             (-3.1992662084989276e-01_DP, 1.1706238784128796e-01_DP), &
             (6.6644795549324654e-01_DP, 7.1243124265282953e-02_DP), &
             (3.1865824911272272e-01_DP, 0.0000000000000000e+00_DP), &
             (-6.6644795549324654e-01_DP, 7.1243124265282953e-02_DP), &
             (-3.1992662084989276e-01_DP, -1.1706238784128796e-01_DP), &
             (-6.5034714907320024e-01_DP, -1.6210196121135434e-01_DP), &
             (-6.5308074464460597e-01_DP, -9.0556255790004966e-02_DP) /

        ! Dmat for il = 2, alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        data DMAT_REF(116:140) &
             / (4.1831402356284814e-01_DP, -1.1828109392712773e-01_DP), &
             (-5.7989806439506730e-01_DP, 2.3300394891357895e-01_DP), &
             (4.8582489840724097e-01_DP, -2.5823145079969545e-01_DP), &
             (-2.6740985684212920e-01_DP, 1.8100791321078466e-01_DP), &
             (8.8649440081326927e-02_DP, -7.4902748341365663e-02_DP), &
             (-6.2465224401252539e-01_DP, 1.9549407883201121e-02_DP), &
             (2.3686161140923923e-01_DP, -3.2843259957511271e-02_DP), &
             (3.5894758329289211e-01_DP, -8.9469304673317804e-02_DP), &
             (-5.2382113453904611e-01_DP, 1.9166817893420654e-01_DP), &
             (2.8973637020666171e-01_DP, -1.4256483423593375e-01_DP), &
             (5.3775766473774822e-01_DP, 1.1630136762018339e-01_DP), &
             (3.6783413805334331e-01_DP, 3.9321379847211700e-02_DP), &
             (-3.4768538040862079e-01_DP, 0.0000000000000000e+00_DP), &
             (-3.6783413805334320e-01_DP, 3.9321379847211700e-02_DP), &
             (5.3775766473774822e-01_DP, -1.1630136762018339e-01_DP), &
             (-2.8973637020666171e-01_DP, -1.4256483423593377e-01_DP), &
             (-5.2382113453904622e-01_DP, -1.9166817893420654e-01_DP), &
             (-3.5894758329289211e-01_DP, -8.9469304673317804e-02_DP), &
             (2.3686161140923925e-01_DP, 3.2843259957511271e-02_DP), &
             (6.2465224401252550e-01_DP, 1.9549407883201121e-02_DP), &
             (8.8649440081326927e-02_DP, 7.4902748341365663e-02_DP), &
             (2.6740985684212920e-01_DP, 1.8100791321078466e-01_DP), &
             (4.8582489840724091e-01_DP, 2.5823145079969545e-01_DP), &
             (5.7989806439506741e-01_DP, 2.3300394891357895e-01_DP), &
             (4.1831402356284814e-01_DP, 1.1828109392712773e-01_DP) /

        ! Rotation set 4
        ! Dmat for il = 0, alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        data DMAT_REF(141:141) &
             / (1.0000000000000000e+00_DP, 0.0000000000000000e+00_DP) /

        ! Dmat for il = 1, alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        data DMAT_REF(142:150) &
             / (-7.3535577783286710e-01_DP, 5.7854201068142952e-01_DP), &
             (-3.0286570236253252e-01_DP, -1.6933312554065497e-01_DP), &
             (7.4099484707884564e-03_DP, -6.3912232594790980e-02_DP), &
             (-1.3332595607191797e-01_DP, 3.2035219760511452e-01_DP), &
             (8.7131929730390190e-01_DP, 0.0000000000000000e+00_DP), &
             (1.3332595607191797e-01_DP, 3.2035219760511452e-01_DP), &
             (7.4099484707884564e-03_DP, 6.3912232594790980e-02_DP), &
             (3.0286570236253257e-01_DP, -1.6933312554065497e-01_DP), &
             (-7.3535577783286710e-01_DP, -5.7854201068142952e-01_DP) /

        ! Dmat for il = 2, alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        data DMAT_REF(151:175) &
             / (2.0603726186886970e-01_DP, -8.5086842054726708e-01_DP), &
             (4.5351052994690194e-01_DP, -7.1701250101670846e-02_DP), &
             (7.7224972810712000e-02_DP, 1.2562256155890258e-01_DP), &
             (-1.8479072328660889e-02_DP, 2.5600197822101579e-02_DP), &
             (-4.0298661389109227e-03_DP, -9.4717270036089529e-04_DP), &
             (-1.2345404501653599e-01_DP, -4.4223519625637808e-01_DP), &
             (-5.4610358138652892e-01_DP, 4.2964762573402993e-01_DP), &
             (-4.5707561777078037e-01_DP, -2.5555235327672043e-01_DP), &
             (2.0322810660039495e-02_DP, -1.7528815578202459e-01_DP), &
             (3.0352364311776948e-02_DP, -8.6936887075500843e-03_DP), &
             (-1.0391925934524647e-01_DP, -1.0462080068213434e-01_DP), &
             (-2.0121143880301820e-01_DP, 4.8346569942512935e-01_DP), &
             (6.3879597678124822e-01_DP, 0.0000000000000000e+00_DP), &
             (2.0121143880301817e-01_DP, 4.8346569942512929e-01_DP), &
             (-1.0391925934524644e-01_DP, 1.0462080068213432e-01_DP), &
             (-3.0352364311776951e-02_DP, -8.6936887075500843e-03_DP), &
             (2.0322810660039499e-02_DP, 1.7528815578202461e-01_DP), &
             (4.5707561777078032e-01_DP, -2.5555235327672038e-01_DP), &
             (-5.4610358138652892e-01_DP, -4.2964762573402993e-01_DP), &
             (1.2345404501653599e-01_DP, -4.4223519625637808e-01_DP), &
             (-4.0298661389109227e-03_DP, 9.4717270036089529e-04_DP), &
             (1.8479072328660889e-02_DP, 2.5600197822101579e-02_DP), &
             (7.7224972810711986e-02_DP, -1.2562256155890258e-01_DP), &
             (-4.5351052994690189e-01_DP, -7.1701250101670846e-02_DP), &
             (2.0603726186886970e-01_DP, 8.5086842054726708e-01_DP) /

        ! Arguments
        integer,intent(out) :: outcome

        ! Local variables
        complex(kind=DP),allocatable :: Dmat(:,:)
        type(DEM) :: Dmat_dem, Dmat_dem_ref
        real(kind=DP) :: alpha_ang, beta_ang, gamma_ang
        integer :: irot ! rotation set counter
        integer :: il ! angular momentum counter
        integer :: nsh ! number of spherical harmonics (2*l + 1)
        integer :: icol ! column counter for putting columns in DEM instance
        integer :: istart_ref ! start index in DMAT_REF array
        integer :: ierr_alloc ! holds status returned by allocate/deallocate

        if (pub_on_root.and.sph_harm_rot_debug) then
           write(stdout,'(a)') module_short_code//&
                " DEBUG: Unit test for "//unit_test_routine
           write(stdout,'(a,es16.8)') module_short_code//&
                " DEBUG: Error threshold = ", ERR_THRESHOLD
        end if

        ! Assume everything is fine, unless it's not
        outcome = TEST_PASSED

        ! Counter for position in 1-D reference data array
        istart_ref = 1

        ! -------------------------------------------------------------------
        ! # Note on use of local allocatable array and DEM object
        ! We are storing the D matrix in a local complex array and an
        ! instance of DEM.
        !
        ! This is wasteful, since we have two copies of D, however, it
        ! serves a purpose, as it allows us
        !
        ! (i)  to ensure that Dmat is always evaluated locally on an MPI
        !      process by sph_harm_rot_get_Wigner_Dmat; and
        ! (ii) to use the dense module, which provides an agnostic interface
        !      to LAPACK and ScaLAPACK.
        !
        ! In this case, we are not concerned about whether Dmat is distributed
        ! as a DEM, because we expect all MPI ranks to hold the same Dmat,
        ! so it is not important whether operations are performed locally or
        ! in a distributed fashion.
        !
        ! If we had a different Dmat (e.g. different l values) on each MPI
        ! rank, then we would want to avoid using the dense module to do
        ! matrix operations, since when using ScaLAPACK, we would be
        ! attempting to distribute different data into a single distributed
        ! object.
        ! -------------------------------------------------------------------
        do irot = 0, NROT_REF-1
           do il = 0, LMAX_REF
              ! Get alpha, beta and gamma angles for rotation set
              alpha_ang = DMAT_EULER_ANGLES_REF(irot*3+1)
              beta_ang  = DMAT_EULER_ANGLES_REF(irot*3+2)
              gamma_ang = DMAT_EULER_ANGLES_REF(irot*3+3)

              if (sph_harm_rot_debug.and.pub_on_root) then
                 write(stdout,'(a,i0)') module_short_code//" DEBUG: &
                      &Testing D matrix with l = ",il
                 write(stdout,'(a,f10.4,a,f10.4,a,f10.4)') &
                      module_short_code//" DEBUG: &
                      &alpha = ",alpha_ang,", &
                      &beta  = ",beta_ang,", &
                      &gamma = ",gamma_ang
              end if

              nsh = 2*il + 1

              ! Allocate local array for D matrix
              allocate(Dmat(nsh,nsh),stat=ierr_alloc)
              call utils_alloc_check(myself,"Dmat",ierr_alloc)

              ! Initialize to garbage
              Dmat(:,:) = garbage_complex

              ! Evaluate D matrix for this angular momentum and set of angles
              call sph_harm_rot_get_Wigner_Dmat(Dmat,il,alpha_ang,beta_ang,&
                   gamma_ang)

              ! Allocate workspace for dense matrix operations (DEM instance)
              ! (this may be distributed, if using ScaLAPACK)
              call dense_create(Dmat_dem,nsh,nsh,iscmplx=.true.)

              ! Put local Dmat (identical on each MPI rank) into DEM
              do icol = 1, nsh
                 call dense_put_col(Dmat(:,icol),Dmat_dem,icol)
              end do

              ! [ Check matrix is unitary ]
              call internal_check_unitary(outcome,Dmat_dem,nsh,"D",ERR_THRESHOLD)

              ! [ Check matrix against known results ]
              ! Allocate workspace for reference Dmat with total ang mom il
              call dense_create(Dmat_dem_ref,nsh,nsh,iscmplx=.true.)
              ! Check against pre-computed reference data

              ! Put local DMAT_REF (identical on each MPI rank) into DEM
              do icol = 1, nsh
                 ! This takes the 1-D array containing the refence data and
                 ! extracts columns to put in the Dmat_dem_ref DEM instance
                 call dense_put_col(&
                      DMAT_REF(istart_ref+(icol-1)*nsh:istart_ref+icol*nsh-1),&
                      Dmat_dem_ref,icol)
              end do

              call internal_check_ref(outcome,Dmat_dem,Dmat_dem_ref,"D",&
                   ERR_THRESHOLD)

              ! Deallocate workspace for reference Dmat with total ang mom il
              call dense_destroy(Dmat_dem_ref)

              ! Update start point in reference data array
              ! (Start point in reference data array is the sum over the squares
              ! of angular momentum components for for l = 0 .. il-1)
              istart_ref = istart_ref + nsh*nsh

              ! Deallocate workspace for dense matrix operations (DEM instance)
              call dense_destroy(Dmat_dem)

              ! Deallocate local array for D matrix
              deallocate(Dmat,stat=ierr_alloc)
              call utils_dealloc_check(myself,"Dmat",ierr_alloc)
           end do ! il
        end do ! irot

      end subroutine internal_test_sph_harm_rot_get_Wigner_Dmat

      subroutine internal_test_sph_harm_rot_get_RSH_Rmat(outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Unit test for sph_harm_rot_get_RSH_Rmat routine           !
        !----------------------------------------------------------------------!
        ! ## TESTS                                                             !
        ! * Test unitarity of matrix up to LMAX                                !
        ! * Check matrix agrees with precomputed reference up to LMAX_REF      !
        !                                                                      !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 11/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use constants, only: garbage_real
        use dense, only: DEM, dense_create, dense_copy, &
             dense_destroy, dense_put_col, dense_put_element
        use utils, only: utils_alloc_check, utils_dealloc_check
        implicit none

        ! Parameters
        character(len=*), parameter :: unit_test_routine = &
             "sph_harm_rot_get_RSH_Rmat"
        integer, parameter :: NROT_REF = 5
        integer, parameter :: LMAX_REF = 2
        ! --> max angular momentum in reference data
        real(kind=DP), parameter :: ERR_FACTOR = 1.0e2_DP
        real(kind=DP), parameter :: ERR_THRESHOLD = ERR_FACTOR * EPSILON(1.0_DP)
        ! --> rough error threshold for unit testing based on machine
        !     epsilon (to adjust, change ERR_FACTOR)

        ! Reference data (do not modify)
        ! Store matrices sequentially in 1-D array up to LMAX_REF
        ! We have 5 rotation sets, each containing RSH rotation matrices up
        ! to LMAX_REF (rotation sets are numbered 0 to 4)
        ! Each rotation set has:
        !   l = 0 (1x1 matrix, 1 element)
        !   l = 1 (3x3 matrix, 9 elements)
        !   l = 2 (5x5 matrix, 25 elements)
        ! Euler angles for each rotation set
        !   0: alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        !   1: alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        !   2: alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        !   3: alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        !   4: alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        ! Start indices for each reference matrix:
        ! Array elements per rotation set = 35
        ! Rotation set   l   Start index   End index
        ! +--------------+---+-------------+-----------+
        ! 0              0   1             1
        ! 0              1   2             10
        ! 0              2   11            35
        ! 1              0   36            36
        ! 1              1   37            45
        ! 1              2   46            70
        ! 2              0   71            71
        ! 2              1   72            80
        ! 2              2   81            105
        ! 3              0   106           106
        ! 3              1   107           115
        ! 3              2   116           140
        ! 4              0   141           141
        ! 4              1   142           150
        ! 4              2   151           175

        ! Euler angles for each rotation set
        !   0: alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        !   1: alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        !   2: alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        !   3: alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        !   4: alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        real(kind=DP) :: RMAT_EULER_ANGLES_REF(5*3)

        ! Rotation set 0
        data RMAT_EULER_ANGLES_REF(1:3) &
             / 0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP /

        ! Rotation set 1
        data RMAT_EULER_ANGLES_REF(4:6) &
             / 5.4641481231202804e-01_DP, &
             8.1804498047957241e-01_DP, &
             6.1422738037708546e+00_DP /

        ! Rotation set 2
        data RMAT_EULER_ANGLES_REF(7:9) &
             / 3.2078413032109494e+00_DP, &
             1.9773138676292171e-01_DP, &
             2.2398638613670325e+00_DP /

        ! Rotation set 3
        data RMAT_EULER_ANGLES_REF(10:12) &
             / 1.0649532499842701e-01_DP, &
             1.2464827205876983e+00_DP, &
             2.8973157492962169e+00_DP /

        ! Rotation set 4
        data RMAT_EULER_ANGLES_REF(13:15) &
             / 1.9651774322985169e+00_DP, &
             5.1291187749676292e-01_DP, &
             5.0980521145061619e-01_DP /

        ! Pre-computed matrices from Sage Math implementation of RSH rotation
        ! matrix, Rmat
        ! Rotation set 0
        real(kind=DP) :: RMAT_REF(5*(1+9+25))

        ! Rotation set 0
        ! Rmat for il = 0, alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        data RMAT_REF(1:1) &
             / 1.0000000000000000e+00_DP /

        ! Rmat for il = 1, alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        data RMAT_REF(2:10) &
             / 1.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! Rmat for il = 2, alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        data RMAT_REF(11:35) &
             / 1.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! Rotation set 1
        ! Rmat for il = 0, alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        data RMAT_REF(36:36) &
             / 1.0000000000000000e+00_DP /

        ! Rmat for il = 1, alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        data RMAT_REF(37:45) &
             / 8.9581688197178955e-01_DP, &
             -1.0249872900457570e-01_DP, &
             -4.3244204759341742e-01_DP, &
             3.7922963644088564e-01_DP, &
             6.8364930698757931e-01_DP, &
             6.2354511296322102e-01_DP, &
             2.3172612459378508e-01_DP, &
             -7.2257707935406368e-01_DP, &
             6.5129514628384610e-01_DP /

        ! Rmat for il = 2, alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        data RMAT_REF(46:70) &
             / 4.8323306738713168e-01_DP, &
             -6.7104837945925877e-01_DP, &
             1.2828128121527649e-01_DP, &
             2.4571578703899863e-01_DP, &
             -4.8923158105165160e-01_DP, &
             3.9458739837285617e-01_DP, &
             5.7355403481173006e-01_DP, &
             -1.2137031675665666e-01_DP, &
             -3.5955128770527411e-01_DP, &
             -6.0936743588445352e-01_DP, &
             4.0957248850934164e-01_DP, &
             4.4905162771615753e-01_DP, &
             2.0106456241689621e-01_DP, &
             7.3834932986371904e-01_DP, &
             2.1217050014465635e-01_DP, &
             3.9148211407729738e-01_DP, &
             -1.1560323861450013e-01_DP, &
             -8.5561459984920907e-01_DP, &
             -5.3019311691385740e-03_DP, &
             3.1823449157840200e-01_DP, &
             5.3831098692372925e-01_DP, &
             -7.5619894496062129e-02_DP, &
             4.4306868244679132e-01_DP, &
             -5.1493570484572104e-01_DP, &
             4.9298506612625681e-01_DP /

        ! Rotation set 2
        ! Rmat for il = 0, alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        data RMAT_REF(71:71) &
             / 1.0000000000000000e+00_DP /

        ! Rmat for il = 1, alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        data RMAT_REF(72:80) &
             / 6.6980985817564243e-01_DP, &
             1.5409186813769893e-01_DP, &
             7.2636798529709445e-01_DP, &
             -1.3004726847145001e-02_DP, &
             9.8051475932109300e-01_DP, &
             -1.9601449903803170e-01_DP, &
             -7.4241877062097317e-01_DP, &
             1.2184622656173390e-01_DP, &
             6.5876237453448239e-01_DP /

        ! Rmat for il = 2, alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        data RMAT_REF(81:105) &
             / -9.8023694004310921e-02_DP, &
             -3.2786891572927357e-02_DP, &
             3.2520141893724881e-02_DP, &
             1.9001512305454488e-01_DP, &
             9.7578331023671072e-01_DP, &
             -1.4073866104033353e-01_DP, &
             6.5475452922548827e-01_DP, &
             2.6169443241525436e-01_DP, &
             6.8201028994328217e-01_DP, &
             -1.3366796251017499e-01_DP, &
             4.4151967256333491e-03_DP, &
             -2.2085945560331714e-02_DP, &
             9.4211378986975103e-01_DP, &
             -3.3289169435650751e-01_DP, &
             3.3127689510263264e-02_DP, &
             1.3695781866170365e-01_DP, &
             -7.2953713908467577e-01_DP, &
             2.0693161480474226e-01_DP, &
             6.2204260405729694e-01_DP, &
             -1.3878193014759990e-01_DP, &
             -9.7560598944841737e-01_DP, &
             -1.9367317807209894e-01_DP, &
             -7.7057216162050458e-03_DP, &
             -3.1659690271971561e-02_DP, &
             -9.8091484411666693e-02_DP /

        ! Rotation set 3
        ! Rmat for il = 0, alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        data RMAT_REF(106:106) &
             / 1.0000000000000000e+00_DP /

        ! Rmat for il = 1, alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        data RMAT_REF(107:115) &
             / -9.7300736549449873e-01_DP, &
             2.2924679203237475e-01_DP, &
             2.6506132051283000e-02_DP, &
             1.0075299256179493e-01_DP, &
             3.1865824911272261e-01_DP, &
             9.4249973727437020e-01_DP, &
             2.0761864363129293e-01_DP, &
             9.1972975846999694e-01_DP, &
             -3.3315412379471321e-01_DP /

        ! Rmat for il = 2, alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        data RMAT_REF(116:140) &
             / 3.2966458348152117e-01_DP, &
             -8.4730792123719656e-01_DP, &
             3.6519441995220991e-01_DP, &
             -5.1996035702794291e-02_DP, &
             1.9318384226849336e-01_DP, &
             -9.1438861421918727e-01_DP, &
             -2.8695952312980666e-01_DP, &
             1.2652870408509656e-01_DP, &
             2.2451143889171776e-01_DP, &
             1.2301542635273262e-01_DP, &
             1.6447497141100248e-01_DP, &
             5.5608828671150884e-02_DP, &
             -3.4768538040862090e-01_DP, &
             5.2019602673885534e-01_DP, &
             7.6050418274220744e-01_DP, &
             1.6211424211913489e-01_DP, &
             1.5882491897669521e-01_DP, &
             5.0762854047385420e-01_DP, &
             7.6068274594828555e-01_DP, &
             -3.3491587380586380e-01_DP, &
             -4.3378345585762063e-02_DP, &
             4.1401186212436358e-01_DP, &
             6.8706016026605110e-01_DP, &
             -3.1248820755293805e-01_DP, &
             5.0696346364417499e-01_DP /

        ! Rotation set 4
        ! Rmat for il = 0, alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        data RMAT_REF(141:141) &
             / 1.0000000000000000e+00_DP /

        ! Rmat for il = 1, alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        data RMAT_REF(142:150) &
             / -7.2794582936207863e-01_DP, &
             2.3947320269862019e-01_DP, &
             -6.4245424327622058e-01_DP, &
             4.5304642258917877e-01_DP, &
             8.7131929730390190e-01_DP, &
             -1.8855137529326590e-01_DP, &
             5.1462977808663857e-01_DP, &
             -4.2831678385874672e-01_DP, &
             -7.4276572630365556e-01_DP /

        ! Rmat for il = 2, alpha = 1.9652, beta = 0.5129, gamma = 0.5098
        data RMAT_REF(151:175) &
             / 2.1006712800778060e-01_DP, &
             4.3503145761824108e-01_DP, &
             -1.7765713029664904e-01_DP, &
             9.7301447923772433e-02_DP, &
             8.5181559324762790e-01_DP, &
             -1.5380640932831294e-01_DP, &
             -5.2578077072648943e-01_DP, &
             3.6140560390029852e-01_DP, &
             -6.0493578151605443e-01_DP, &
             4.5092888496392824e-01_DP, &
             -1.4795615523100672e-01_DP, &
             6.8372374906921229e-01_DP, &
             6.3879597678124811e-01_DP, &
             -2.8455594565983239e-01_DP, &
             -1.4696402595781455e-01_DP, &
             -4.3354150754882803e-01_DP, &
             2.5435946995200520e-01_DP, &
             -6.4640253768149858e-01_DP, &
             -5.6642639204656842e-01_DP, &
             -9.3101680704759049e-02_DP, &
             -8.4992124784690626e-01_DP, &
             -4.6101052279569267e-02_DP, &
             1.0921260390280242e-01_DP, &
             4.7198960227556280e-01_DP, &
             2.0200739572995879e-01_DP /

        ! Arguments
        integer,intent(out) :: outcome

        ! Local variables
        real(kind=DP),allocatable :: Rmat(:,:)
        type(DEM) :: Rmat_dem, Rmat_dem_ref
        real(kind=DP) :: alpha_ang, beta_ang, gamma_ang
        integer :: irot ! rotation set counter
        integer :: il ! angular momentum counter
        integer :: nsh ! number of spherical harmonics (2*l + 1)
        integer :: icol ! column counter for putting columns in DEM instance
        integer :: istart_ref ! start index in RMAT_REF array
        integer :: ierr_alloc ! holds status returned by allocate/deallocate

        if (pub_on_root.and.sph_harm_rot_debug) then
           write(stdout,'(a)') module_short_code//&
                " DEBUG: Unit test for "//unit_test_routine
           write(stdout,'(a,es16.8)') module_short_code//&
                " DEBUG: Error threshold = ", ERR_THRESHOLD
        end if

        ! Assume everything is fine, unless it's not
        outcome = TEST_PASSED

        ! Counter for position in 1-D reference data array
        istart_ref = 1

        ! -------------------------------------------------------------------
        ! # Note on use of local allocatable array and DEM object
        ! We are storing the R matrix in a local real array and an
        ! instance of DEM.
        !
        ! This is wasteful, since we have two copies of R, however, it
        ! serves a purpose, as it allows us
        !
        ! (i)  to ensure that Rmat is always evaluated locally on an MPI
        !      process by sph_harm_rot_get_RSH_Rmat; and
        ! (ii) to use the dense module, which provides an agnostic interface
        !      to LAPACK and ScaLAPACK.
        !
        ! In this case, we are not concerned about whether Rmat is distributed
        ! as a DEM, because we expect all MPI ranks to hold the same Rmat,
        ! so it is not important whether operations are performed locally or
        ! in a distributed fashion.
        !
        ! If we had a different Rmat (e.g. different l values) on each MPI
        ! rank, then we would want to avoid using the dense module to do
        ! matrix operations, since when using ScaLAPACK, we would be
        ! attempting to distribute different data into a single distributed
        ! object.
        ! -------------------------------------------------------------------
        do irot = 0, NROT_REF-1
           do il = 0, LMAX_REF
              ! Get alpha, beta and gamma angles for rotation set
              alpha_ang = RMAT_EULER_ANGLES_REF(irot*3+1)
              beta_ang  = RMAT_EULER_ANGLES_REF(irot*3+2)
              gamma_ang = RMAT_EULER_ANGLES_REF(irot*3+3)

              if (sph_harm_rot_debug.and.pub_on_root) then
                 write(stdout,'(a,i0)') module_short_code//" DEBUG: &
                      &Testing R matrix with l = ",il
                 write(stdout,'(a,f10.4,a,f10.4,a,f10.4)') &
                      module_short_code//" DEBUG: &
                      &alpha = ",alpha_ang,", &
                      &beta  = ",beta_ang,", &
                      &gamma = ",gamma_ang
              end if

              nsh = 2*il + 1

              ! Allocate local array for D matrix
              allocate(Rmat(nsh,nsh),stat=ierr_alloc)
              call utils_alloc_check(myself,"Rmat",ierr_alloc)

              ! Initialize to garbage
              Rmat(:,:) = garbage_real

              ! Evaluate D matrix for this angular momentum and set of angles
              call sph_harm_rot_get_RSH_Rmat(Rmat,il,alpha_ang,beta_ang,&
                   gamma_ang)

              ! Allocate workspace for dense matrix operations (DEM instance)
              ! (this may be distributed, if using ScaLAPACK)
              call dense_create(Rmat_dem,nsh,nsh,iscmplx=.false.)

              ! Put local Rmat (identical on each MPI rank) into DEM
              do icol = 1, nsh
                 call dense_put_col(Rmat(:,icol),Rmat_dem,icol)
              end do

              ! [ Check matrix is unitary ]
              call internal_check_unitary(outcome,Rmat_dem,nsh,"R",ERR_THRESHOLD)

              ! [ Check matrix against known results ]
              ! Allocate workspace for reference Rmat with total ang mom il
              call dense_create(Rmat_dem_ref,nsh,nsh,iscmplx=.false.)
              ! Check against pre-computed reference data

              ! Put local RMAT_REF (identical on each MPI rank) into DEM
              do icol = 1, nsh
                 ! This takes the 1-D array containing the refence data and
                 ! extracts columns to put in the Rmat_dem_ref DEM instance
                 call dense_put_col(&
                      RMAT_REF(istart_ref+(icol-1)*nsh:istart_ref+icol*nsh-1),&
                      Rmat_dem_ref,icol)
              end do

              call internal_check_ref(outcome,Rmat_dem,Rmat_dem_ref,"R",&
                   ERR_THRESHOLD)

              ! Deallocate workspace for reference Rmat with total ang mom il
              call dense_destroy(Rmat_dem_ref)

              ! Update start point in reference data array
              ! (Start point in reference data array is the sum over the squares
              ! of angular momentum components for for l = 0 .. il-1)
              istart_ref = istart_ref + nsh*nsh

              ! Deallocate workspace for dense matrix operations (DEM instance)
              call dense_destroy(Rmat_dem)

              ! Deallocate local array for R matrix
              deallocate(Rmat,stat=ierr_alloc)
              call utils_dealloc_check(myself,"Rmat",ierr_alloc)
           end do ! il
        end do ! irot

      end subroutine internal_test_sph_harm_rot_get_RSH_Rmat

      subroutine internal_test_sph_harm_rot_apply_RSH_Rmat_to_atomblock(outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Unit test for sph_harm_rot_apply_RSH_Rmat_to_atomblock    !
        !            routine                                                   !
        !----------------------------------------------------------------------!
        ! ## TESTS                                                             !
        ! * Test unitarity of matrix by passing a unit matrix as Mmat for a    !
        !   range of input arguments                                           !
        ! * Test that a forward rotation followed by a reverse rotation of a   !
        !   randomly generated (precomputed) matrix returns the original       !
        !   matrix                                                             !
        ! * Check that applying a specific RSH rotation matrix to a small      !
        !   example SW electrostatic metric matrix produces a result in        !
        !   agreement with a precomputed result                                !
        !                                                                      !
        ! NOTE: This unit test only tests the case where the num_q_per_l       !
        ! array argument for sph_harm_rot_apply_RSH_Rmat_to_atomblock has all  !
        ! identical values equal to a single value of qmax.                    !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 11/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use constants, only: garbage_real
        use dense, only: DEM, dense_axpy, dense_create, dense_destroy, &
             dense_norm, dense_put_col, dense_put_element
        use utils, only: utils_alloc_check, utils_dealloc_check
        implicit none

        ! Parameters
        character(len=*), parameter :: unit_test_routine = &
             "sph_harm_rot_apply_RSH_Rmat_to_atomblock"
        integer, parameter :: LMAX_REF = 2
        ! --> max angular momentum to test
        integer, parameter :: NROT_REF = 5
        ! --> number of Euler rotations to test
        integer, parameter :: NQ_REF = 6
        ! --> number of qmax values to test
        integer, parameter :: NMAT_MAX_REF = 18*(LMAX_REF*(LMAX_REF+2)+1)
        ! --> maximum dimension of Mmat for tests
        !     (Note: 18 is the maximum qmax value in the RMAT_AB_QMAX_REF
        !     array below)

        integer, parameter :: LMAX_SMALL_REF = 1
        ! --> angular momentum for small example V-matrix
        integer, parameter :: QMAX_SMALL_REF = 2
        ! --> qmax for small example V-matrix
        integer, parameter :: NMAT_SMALL_REF = &
            QMAX_SMALL_REF*(LMAX_SMALL_REF*(LMAX_SMALL_REF+2)+1)
        ! --> dimension of small example V-matrix
        real(kind=DP), parameter :: ABS_ERR_FACTOR = 1.0e2_DP
        real(kind=DP), parameter :: ABS_ERR_THRESHOLD = ABS_ERR_FACTOR * &
             EPSILON(1.0_DP)
        ! --> rough error absolute threshold for unit testing based on machine
        !     epsilon (to adjust, change ERR_FACTOR)
        real(kind=DP), parameter :: REL_ERR_FACTOR = 1.0e2_DP
        real(kind=DP), parameter :: REL_ERR_THRESHOLD = REL_ERR_FACTOR * &
             1.0e-11_DP
        ! --> rough relative error threshold for unit testing

        ! ----------------------------------------------------------------------
        ! Error factors
        ! In this unit test, we use both absolute and relative measures of
        ! error, since the size of elements in the matrices being compared
        ! varies quite significantly (from zero to ~1.0e4).
        !
        ! The value of the ABS_ERR_THRESHOLD is based on machine epsilon, while
        ! the REL_ERR_THRESHOLD is based on the maximum agreement that can be
        ! expected between the largest reference value ~1e4 and a computed
        ! value.
        !
        ! The IEEE 754 double precision floating point format provides ~15
        ! significant digits of agreement in the mantissa. With a number of
        ! the magnitude 1.0e4, the best absolute agreement between the
        ! computed and reference values we can expect is
        !   1.0e4 * 1.0e-15 = 1e-11
        ! We therefore take as out relative error threshold some factor
        ! multiplied by 1.0e-11
        !
        ! In this unit test we check against these thresholds using
        !   absolute(A - Aref) <= (ABS_ERR_THRESHOLD
        !                          + REL_ERR_THRESHOLD * absolute(Aref))
        ! which has the effect that very small reference values Aref are
        ! checked against an error threshold close to ABS_ERR_THRESHOLD
        ! while larger Aref values are checked against a threshold which scales
        ! with the size of Aref.
        !
        ! This is based on the implementation of numpy.isclose():
        !   https://docs.scipy.org/doc/numpy-1.15.1/reference/
        !           generated/numpy.isclose.html
        !
        ! The above error threshold check only applies to checking the result
        ! of rotations of arbitary or precomputed reference matrices. For the
        ! unitarity test it is not necessary, since we expect the result to be
        ! the identity matrix and thus the range of sizes of elements is not
        ! an issue (they should all be either 1 or zero).
        ! ----------------------------------------------------------------------

        ! Reference data (do not modify)
        real(kind=DP) :: RMAT_AB_EULER_ANGLES_REF(NROT_REF*3)

        integer :: RMAT_AB_QMAX_REF(NQ_REF)

        ! Pre-generated random reference vectors (using Sage Math's
        ! RR.random_element() function) to be used to construct an arbitrary
        ! V-matrix, where
        !   V_ij = R1[i]*R2[j]
        ! and R{1,2} are the random vectors.
        ! * Pre-generating the random vectors means the matrix is
        !   identical for each test.
        ! * Using two vectors requires substantially fewer random numbers to
        !   be included in the source compared to an entire pre-generated matrix
        ! * The values in the vectors are random real numbers between
        !   -sqrt(Vmin) and sqrt(Vmax), where Vmin and Vmax are realistic
        !   maximum and minimum values for V-matrix elements obtained from
        !   previous calculations (in this case Vmin = Vmax = 10000)
        real(kind=DP) :: RMAT_AB_RAND_VEC_REF(NMAT_MAX_REF,2)

        ! Pre-generated example V-matrix atomblock before application of
        ! RSH rotation matrices (i.e. matrix in integral frame)
        real(kind=DP) :: RMAT_AB_VMAT_SMALL_BEFORE_REF(&
             NMAT_SMALL_REF,NMAT_SMALL_REF)

        ! Pre-generated example V-matrix atomblock after application of
        ! RSH rotation matrices (i.e. matrix in ONETEP frame)
        real(kind=DP) :: RMAT_AB_VMAT_SMALL_AFTER_REF(&
             NMAT_SMALL_REF,NMAT_SMALL_REF)

        ! Single arbitrary rotation angle for small example
        real(kind=DP) :: RMAT_AB_EULER_ANGLE_SMALL_REF(3)

        ! Euler angles for each rotation set
        !   0: alpha = 0.0000, beta = 0.0000, gamma = 0.0000
        !   1: alpha = 0.5464, beta = 0.8180, gamma = 6.1423
        !   2: alpha = 3.2078, beta = 0.1977, gamma = 2.2399
        !   3: alpha = 0.1065, beta = 1.2465, gamma = 2.8973
        !   4: alpha = 1.9652, beta = 0.5129, gamma = 0.5098

        ! Rotation set 0
        data RMAT_AB_EULER_ANGLES_REF(1:3) &
             / 0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP /

        ! Rotation set 1
        data RMAT_AB_EULER_ANGLES_REF(4:6) &
             / 5.4641481231202804e-01_DP, &
             8.1804498047957241e-01_DP, &
             6.1422738037708546e+00_DP /

        ! Rotation set 2
        data RMAT_AB_EULER_ANGLES_REF(7:9) &
             / 3.2078413032109494e+00_DP, &
             1.9773138676292171e-01_DP, &
             2.2398638613670325e+00_DP /

        ! Rotation set 3
        data RMAT_AB_EULER_ANGLES_REF(10:12) &
             / 1.0649532499842701e-01_DP, &
             1.2464827205876983e+00_DP, &
             2.8973157492962169e+00_DP /

        ! Rotation set 4
        data RMAT_AB_EULER_ANGLES_REF(13:15) &
             / 1.9651774322985169e+00_DP, &
             5.1291187749676292e-01_DP, &
             5.0980521145061619e-01_DP /

        ! Recommendations from Jacek Dziedzic's documentation for SWRI:
        !
        ! Accuracy    qmax
        ! --------    ----
        ! Crude       7
        ! Reasonable  10
        ! Extreme     14
        !
        ! These reference values include the recommended qmax values

        ! qmax values
        data RMAT_AB_QMAX_REF(1:6) &
             / 1, 3, 7, 10, 14, 18 /

        ! Random reference vector 1
        data RMAT_AB_RAND_VEC_REF(1:162,1) &
             / 8.6216971170901804e+00_DP, -3.3625099318433499e+01_DP, &
             -8.0280042543906802e+01_DP, -4.5360236031911398e+00_DP, &
             1.7013138233655500e+01_DP, 2.9715951937536702e+00_DP, &
             -4.6234137422439998e+00_DP, -7.9481058905673706e+01_DP, &
             -7.0975043093662606e+01_DP, 1.3030334609982299e+01_DP, &
             -9.6031138758271993e+01_DP, -8.0611068270235293e+01_DP, &
             -3.3294556824898002e+01_DP, -3.4668569661639403e+01_DP, &
             -3.0995014403311501e+01_DP, -3.0070893020464300e+01_DP, &
             -6.5906423637825199e+01_DP, -1.4630786308382699e+01_DP, &
             -5.1240181017452002e+01_DP, 8.9746022436070803e+01_DP, &
             2.4672466992775600e+01_DP, 6.9605302626563201e+01_DP, &
             3.2140567101369797e+01_DP, 5.0085401090216500e+01_DP, &
             -8.0970057330181405e+01_DP, -8.4188989222868091e+00_DP, &
             -2.5634165158267599e+01_DP, -1.8851361486735801e+01_DP, &
             7.3681815128205898e+01_DP, -9.6493330255627697e+01_DP, &
             3.7078677551298803e+01_DP, -6.5069488921839294e+01_DP, &
             8.9500091611253495e+01_DP, 9.2952155067023099e+01_DP, &
             -8.5412556163818394e+01_DP, 2.2277502276919998e+00_DP, &
             -2.8872181118282601e+01_DP, -4.9646859765006099e+01_DP, &
             -9.6638487911506502e+01_DP, 1.7425329811229901e+01_DP, &
             -2.4539442401747401e+01_DP, 4.3612171941279101e+01_DP, &
             -3.8661653041843998e+01_DP, 5.9488996233333097e+01_DP, &
             2.3045213308759699e+01_DP, -5.6595286044867104e+00_DP, &
             -3.3353332994150300e+01_DP, -6.1644118250514197e+01_DP, &
             9.9326242603193506e+01_DP, -2.0369504334406500e+01_DP, &
             -3.1628747608299101e+01_DP, 6.0028090301505699e+01_DP, &
             7.9188259413964104e+01_DP, 2.6434744769929100e+01_DP, &
             -3.3081820619838098e+01_DP, 4.1143917481945898e+01_DP, &
             9.3681048753244198e+01_DP, 9.4570761556975000e+01_DP, &
             7.4890130366258006e+01_DP, -5.0088519132135502e+01_DP, &
             -2.4233685781180402e+01_DP, -8.1592138041378604e+01_DP, &
             7.9323020779679297e+01_DP, 5.2491022375371699e+01_DP, &
             -4.3045538140553496e+01_DP, -9.9187632474832597e+01_DP, &
             -8.5273872343450293e+01_DP, -9.7287283024706397e+01_DP, &
             8.9438974197288701e+01_DP, -8.9047662718846595e+01_DP, &
             8.5117203783598995e+01_DP, -5.8038402362219102e+01_DP, &
             6.4837181866519202e+01_DP, -4.0676636744110802e+01_DP, &
             -9.3593875015046294e+01_DP, -1.6153408506437901e+00_DP, &
             -2.7213697733422698e+01_DP, -2.0711882017015199e+01_DP, &
             3.4783078136633301e+00_DP, -6.2637688747301297e+01_DP, &
             2.1835779854490500e+01_DP, 2.1572536984528000e+01_DP, &
             4.0190334750500902e+01_DP, 6.1421712931961700e+01_DP, &
             8.3141465357393002e+01_DP, 5.3076748957899902e+01_DP, &
             3.0823148289725999e+01_DP, -5.5065673810024002e+01_DP, &
             -9.8408340642665493e+01_DP, -8.5455390632358402e+01_DP, &
             -9.0068043382945305e+01_DP, 3.2147830376507798e+01_DP, &
             -5.5728332326202597e+01_DP, -1.0296997385104201e+00_DP, &
             6.3731662927797998e+00_DP, 8.3416242319116407e+01_DP, &
             7.9696457401931397e+01_DP, -5.3237510606373597e+01_DP, &
             -4.5699580614125601e+01_DP, -3.2120038842379500e+01_DP, &
             -5.5302871736580101e+01_DP, -5.1339017306654100e+01_DP, &
             -2.2205433601553100e+01_DP, -6.1717429011624297e+01_DP, &
             2.5796472642099399e-01_DP, 8.1477107302883702e+00_DP, &
             5.0021272762937997e+01_DP, 8.0168034891777594e+01_DP, &
             4.6123475718889297e+00_DP, -6.7732006898475106e+01_DP, &
             -1.9888409446181100e+01_DP, -5.4410006745506401e+01_DP, &
             -8.7942937285154699e-01_DP, 3.9326937310901999e+00_DP, &
             -3.5087845797120899e+01_DP, -3.8374993817262002e+01_DP, &
             2.3121100355047801e+01_DP, 7.1099099809235298e+01_DP, &
             -1.3075942021418900e+01_DP, 5.9666770594824897e+01_DP, &
             -5.8406120369069598e+01_DP, -3.4228499674133097e+01_DP, &
             5.7516112759643597e+01_DP, -3.6433480673411800e+01_DP, &
             -3.7162359835507701e+01_DP, -1.0976938487793999e+01_DP, &
             3.8384446006607802e+00_DP, 1.8221385288215100e+01_DP, &
             -9.6934835846791401e+01_DP, -9.7162599240835604e+01_DP, &
             -3.2344111535881197e+01_DP, 2.3691824611358999e+01_DP, &
             8.1464259575998497e+01_DP, 5.6999001351611698e+00_DP, &
             -6.0537229563605301e+01_DP, -2.5965180572576500e+01_DP, &
             5.6905693182789001e+01_DP, -2.0914693031264601e+01_DP, &
             3.9045935186529299e+01_DP, -8.5471264387312303e+01_DP, &
             5.7798911808592500e+01_DP, 6.1896679800646197e+01_DP, &
             -1.1942597288204000e+01_DP, 4.8618170270118398e+01_DP, &
             2.4734982610803900e+01_DP, -7.4241083124253507e+01_DP, &
             -1.9002857170983201e+01_DP, 9.6928475936147905e+01_DP, &
             -8.8661545011650006e+01_DP, -4.0176995885674103e+01_DP, &
             -3.1483016869846502e+01_DP, 8.2099059792567800e+01_DP, &
             -6.4670613489203603e+01_DP, -7.2578076549908701e+01_DP, &
             -6.4934791935166899e+01_DP, 9.4905496246843697e+01_DP, &
             -3.7180863742984897e+01_DP, -7.9621835491829003e+01_DP, &
             2.1810107015411901e+01_DP, -6.6169591943958096e+01_DP, &
             4.9963274404589100e+01_DP, 7.9845837584801103e+01_DP /

        ! Random reference vector 2
        data RMAT_AB_RAND_VEC_REF(1:162,2) &
             / 7.3971911186041794e+01_DP, -4.0415035957346902e+01_DP, &
             -7.8269063255247104e+01_DP, 1.3597744636880700e+01_DP, &
             -6.5740143612787094e+01_DP, 2.0619917773725600e+01_DP, &
             1.6034886888116400e+00_DP, 3.1223840121201000e+01_DP, &
             -4.8340183065077099e+01_DP, 3.1080744969820699e+01_DP, &
             -5.1482218588883502e+00_DP, -1.3755629200934999e-01_DP, &
             3.8078632203480602e+01_DP, 1.7267552162366400e+01_DP, &
             -8.6090103023333199e+01_DP, 3.3918478643896997e+01_DP, &
             2.8939364650910100e+01_DP, -8.8305612732365503e+01_DP, &
             1.1361870128251400e+01_DP, -9.3004942300596099e+01_DP, &
             8.6407323366564000e+01_DP, -1.3073692380243500e+00_DP, &
             4.8883295795623198e+01_DP, 8.4027344334938803e+01_DP, &
             -5.0438841866876501e+01_DP, 8.3553652063159206e-02_DP, &
             -9.3278619016254297e+01_DP, -7.5871146851740804e+01_DP, &
             1.8011542804541700e+00_DP, 6.1425391013155100e+01_DP, &
             2.0180504608536398e+00_DP, -7.1698910879802904e+00_DP, &
             -7.7387782744469405e+01_DP, -2.8658316359273201e+01_DP, &
             1.8102857320498298e+01_DP, 3.1462233299931700e+01_DP, &
             3.2249500874563800e+01_DP, 2.3240081820880999e+01_DP, &
             -7.8534674201559000e+01_DP, -8.9980469553826893e+01_DP, &
             -6.8638667577769796e+00_DP, 8.4271601009807000e+00_DP, &
             -6.0827095386258399e+01_DP, -2.8093162536794100e+01_DP, &
             -1.3842126600065701e+01_DP, -6.4887745184677200e+01_DP, &
             3.1169577460741799e+01_DP, 6.5381799255762104e+01_DP, &
             -6.1238373039226701e+01_DP, -7.8036898808830003e+01_DP, &
             3.8248432080583903e+01_DP, 1.4206794494240301e+01_DP, &
             9.8081513260790203e+01_DP, 5.4789337034934398e+01_DP, &
             -7.4337582552575896e+01_DP, -8.2187255414300793e+01_DP, &
             7.8533884799913906e+01_DP, 3.9330089763389601e+01_DP, &
             -4.0018800159020202e+01_DP, -8.9920004980129804e+01_DP, &
             3.5647004388468801e+01_DP, -8.7951542263618194e+01_DP, &
             -8.7347810001807005e+01_DP, -8.4493287142330004e+01_DP, &
             2.5256973783203598e+01_DP, -1.4584556198765600e+01_DP, &
             -1.3210774276352399e+01_DP, -6.9979437264637895e+01_DP, &
             -3.0439288578826901e+01_DP, -1.4419830673361700e+01_DP, &
             2.7395536368814099e+01_DP, 1.9938817690441098e+01_DP, &
             -6.8047388310696206e+01_DP, -2.3726795002355001e+01_DP, &
             -8.0699877923048106e+01_DP, 3.1588751246709499e+01_DP, &
             -1.7589778574144098e+01_DP, 7.9204723028376605e+01_DP, &
             -8.8727659901788897e+01_DP, -1.4074151076236801e+01_DP, &
             8.7746218030344295e+01_DP, 2.3068340110083799e+01_DP, &
             6.0916270771968698e+01_DP, -1.4595855653759600e+01_DP, &
             3.9053296761345500e+01_DP, 9.7048012395869407e+01_DP, &
             5.6699094526946602e+01_DP, 5.2739748544210300e+01_DP, &
             3.9543264351517699e+01_DP, -1.5939367710816200e+01_DP, &
             -3.6064154263105600e+01_DP, -3.0719286093614400e+01_DP, &
             -4.0967212027042002e+01_DP, -9.9659572904582106e+01_DP, &
             5.9508362660704300e+01_DP, 4.5544489875982501e+00_DP, &
             -3.1429229455471301e+01_DP, -3.9527485303995398e+01_DP, &
             -3.2371881846293697e+01_DP, -8.2739776070906700e+01_DP, &
             5.8675716502001997e+01_DP, -2.0766055644295800e+01_DP, &
             1.7323047009759001e+01_DP, 4.2318973943895003e+01_DP, &
             -4.5982789866013498e+01_DP, -3.0994820613553799e+01_DP, &
             1.5839119376556500e+01_DP, -1.3930246423549301e+01_DP, &
             8.5141713450831293e+00_DP, -2.0164365572172500e+01_DP, &
             6.0568231746513398e+01_DP, 7.5486646496305099e+01_DP, &
             9.6001897763398304e+01_DP, 9.0450666827812000e+01_DP, &
             -6.6799976491790600e+00_DP, -7.1174513136838996e+01_DP, &
             -9.5862181105324296e+01_DP, 8.9531509848305106e+01_DP, &
             -9.6289087210788907e+01_DP, 3.6145395150535599e+01_DP, &
             4.8633022957696902e+01_DP, 2.9248753728318899e+01_DP, &
             8.4227629096738994e+01_DP, 6.1582137078216903e+01_DP, &
             -1.3128519446880600e+00_DP, -5.8925158868156799e+01_DP, &
             8.8745446581858701e+01_DP, -5.2132641085842102e-01_DP, &
             -5.8307737254773102e+01_DP, 2.9397899270521901e+01_DP, &
             -1.9191007564075701e+01_DP, -3.1507238275173801e+01_DP, &
             -5.3642913068600002e+01_DP, 1.1134155861809299e+01_DP, &
             -9.4527749761112698e+00_DP, -4.9308512598376701e+01_DP, &
             5.3872471351058003e+01_DP, 4.8974257854822199e+01_DP, &
             -9.7811373726836806e+01_DP, 7.8935291330783699e+01_DP, &
             1.4615419051743100e+01_DP, -5.1236433628344599e+01_DP, &
             -4.9518611440548497e+01_DP, -8.0548819045553103e+01_DP, &
             -5.2236255061571498e+01_DP, 7.5839119358073404e+01_DP, &
             7.8920840819719103e+01_DP, -7.3349257896073894e+01_DP, &
             2.0278386737485398e+01_DP, 8.4276337309551707e+01_DP, &
             -8.5488637199434706e+01_DP, -2.1278505528534499e+01_DP, &
             7.3943979373552494e+01_DP, -2.6381439180967099e+01_DP, &
             1.5321760003460300e+01_DP, -1.6960626157996298e+01_DP, &
             5.5078106823006202e+01_DP, 3.0734173454331099e+01_DP, &
             2.0363612259719101e+01_DP, 1.7282934958898700e+01_DP, &
             -9.5993779779164100e+01_DP, 5.8315240259467103e+01_DP /

        ! Integral coordinate frame V-matrix evaluated in Sage Math
        ! using 2Dn-1Da approach
        ! lmax = 1
        ! qmax = 2
        ! R = 1.4 bohr
        ! a = 8.0 bohr

        ! Example V-matrix column 1
        data RMAT_AB_VMAT_SMALL_BEFORE_REF(1:8,1) &
             / 6.2359841065007004e+03_DP, &
             -1.0567725951459429e+03_DP, &
             0.0000000000000000e+00_DP, &
             3.8681441064178716e+02_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             -8.1484379495643324e+01_DP, &
             0.0000000000000000e+00_DP /

        ! Example V-matrix column 2
        data RMAT_AB_VMAT_SMALL_BEFORE_REF(1:8,2) &
             / -1.0567725951459429e+03_DP, &
             3.7120632555166486e+02_DP, &
             0.0000000000000000e+00_DP, &
             4.7093586218590666e+01_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             4.5808863708212975e+01_DP, &
             0.0000000000000000e+00_DP /

        ! Example V-matrix column 3
        data RMAT_AB_VMAT_SMALL_BEFORE_REF(1:8,3) &
             / 0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             7.7288979892783289e+02_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             -1.1025318102754959e+02_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP /

        ! Example V-matrix column 4
        data RMAT_AB_VMAT_SMALL_BEFORE_REF(1:8,4) &
             / -3.8681441064178716e+02_DP, &
             -4.7093586218590666e+01_DP, &
             0.0000000000000000e+00_DP, &
             7.1588171041705891e+02_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             -1.1017211365920036e+02_DP, &
             0.0000000000000000e+00_DP /

        ! Example V-matrix column 5
        data RMAT_AB_VMAT_SMALL_BEFORE_REF(1:8,5) &
             / 0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             7.7288979892783289e+02_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             -1.1025318102754959e+02_DP /

        ! Example V-matrix column 6
        data RMAT_AB_VMAT_SMALL_BEFORE_REF(1:8,6) &
             / 0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             -1.1025318102754959e+02_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             8.5005182008246024e+01_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP /

        ! Example V-matrix column 7
        data RMAT_AB_VMAT_SMALL_BEFORE_REF(1:8,7) &
             / 8.1484379495643324e+01_DP, &
             -4.5808863708212975e+01_DP, &
             0.0000000000000000e+00_DP, &
             -1.1017211365920036e+02_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             6.6721613015682905e+01_DP, &
             0.0000000000000000e+00_DP /

        ! Example V-matrix column 8
        data RMAT_AB_VMAT_SMALL_BEFORE_REF(1:8,8) &
             / 0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             -1.1025318102754959e+02_DP, &
             0.0000000000000000e+00_DP, &
             0.0000000000000000e+00_DP, &
             8.5005182008246024e+01_DP /

        ! ONETEP coordinate frame V-matrix evaluated in Sage Math
        ! using 2Dn-1Da approach (with rotation)
        ! lmax = 1
        ! qmax = 2
        ! R = 1.4 bohr
        ! a = 8.0 bohr
        !
        ! Euler rotation angles:
        ! alpha = 0.000000 rad =  0.000000 deg
        ! beta  = 0.818045 rad = 46.870525 deg
        ! gamma = 0.546415 rad = 31.307263 deg

        ! Example V-matrix column 1
        data RMAT_AB_VMAT_SMALL_AFTER_REF(1:8,1) &
             / 6.2359841065007004e+03_DP, &
             -1.0567725951459429e+03_DP, &
             1.4669148831778051e+02_DP, &
             2.6444540376806668e+02_DP, &
             -2.4119623537943502e+02_DP, &
             -3.0901291611744000e+01_DP, &
             -5.5706739572509470e+01_DP, &
             5.0809186617348907e+01_DP /

        ! Example V-matrix column 2
        data RMAT_AB_VMAT_SMALL_AFTER_REF(1:8,2) &
             / -1.0567725951459429e+03_DP, &
             3.7120632555166486e+02_DP, &
             1.7859283580373656e+01_DP, &
             3.2195497581899318e+01_DP, &
             -2.9364975538514315e+01_DP, &
             1.7372078729835703e+01_DP, &
             3.1317197928008266e+01_DP, &
             -2.8563893095654464e+01_DP /

        ! Example V-matrix column 3
        data RMAT_AB_VMAT_SMALL_AFTER_REF(1:8,3) &
             / -1.4669148831778051e+02_DP, &
             -1.7859283580373656e+01_DP, &
             7.6469117399986851e+02_DP, &
             -1.4779921482027220e+01_DP, &
             1.3480519494281651e+01_DP, &
             -1.1024152231447300e+02_DP, &
             2.1017532252979976e-02_DP, &
             -1.9169740083018638e-02_DP /

        ! Example V-matrix column 4
        data RMAT_AB_VMAT_SMALL_AFTER_REF(1:8,4) &
             / -2.6444540376806668e+02_DP, &
             -3.2195497581899318e+01_DP, &
             -1.4779921482027191e+01_DP, &
             7.4624556517714677e+02_DP, &
             2.4301760528504360e+01_DP, &
             2.1017532252976423e-02_DP, &
             -1.1021529205480425e+02_DP, &
             -3.4557899129119107e-02_DP /

        ! Example V-matrix column 5
        data RMAT_AB_VMAT_SMALL_AFTER_REF(1:8,5) &
             / 2.4119623537943502e+02_DP, &
             2.9364975538514315e+01_DP, &
             1.3480519494281680e+01_DP, &
             2.4301760528504474e+01_DP, &
             7.5072456909570985e+02_DP, &
             -1.9169740083025744e-02_DP, &
             -3.4557899129119107e-02_DP, &
             -1.1022166134502234e+02_DP /

        ! Example V-matrix column 6
        data RMAT_AB_VMAT_SMALL_AFTER_REF(1:8,6) &
             / 3.0901291611744000e+01_DP, &
             -1.7372078729835703e+01_DP, &
             -1.1024152231447300e+02_DP, &
             2.1017532252979976e-02_DP, &
             -1.9169740083018638e-02_DP, &
             8.2375728391567463e+01_DP, &
             -4.7401995257258882e+00_DP, &
             4.3234568053042111e+00_DP /

        ! Example V-matrix column 7
        data RMAT_AB_VMAT_SMALL_AFTER_REF(1:8,7) &
             / 5.5706739572509470e+01_DP, &
             -3.1317197928008266e+01_DP, &
             2.1017532252976423e-02_DP, &
             -1.1021529205480425e+02_DP, &
             -3.4557899129119107e-02_DP, &
             -4.7401995257258882e+00_DP, &
             7.6459873811452454e+01_DP, &
             7.7940328621908712e+00_DP /

        ! Example V-matrix column 8
        data RMAT_AB_VMAT_SMALL_AFTER_REF(1:8,8) &
             / -5.0809186617348907e+01_DP, &
             2.8563893095654464e+01_DP, &
             -1.9169740083025744e-02_DP, &
             -3.4557899129119107e-02_DP, &
             -1.1022166134502234e+02_DP, &
             4.3234568053042111e+00_DP, &
             7.7940328621908606e+00_DP, &
             7.7896374829155093e+01_DP /

        ! Euler angles for example V-matrix
        data RMAT_AB_EULER_ANGLE_SMALL_REF(1:3) &
             / 0.0000000000000000e+00_DP, &
             8.1804498047957241e-01_DP, &
             5.4641481231202804e-01_DP /

        ! Arguments
        integer,intent(out) :: outcome

        ! Local variables
        real(kind=DP),allocatable :: Mmat(:,:)
        real(kind=DP),allocatable :: Mmat_ref(:,:)
        type(DEM) :: Mmat_dem, Mmat_dem_ref, Imat_dem
        real(kind=DP) :: alpha_ang, beta_ang, gamma_ang
        real(kind=DP) :: norm
        real(kind=DP) :: max_abs_diff
        integer, allocatable :: num_q_per_l(:)
        integer :: qmax
        integer :: nmat
        integer :: irot
        integer :: il
        integer :: iq
        integer :: irow, icol
        integer :: ierr_alloc ! holds status returned by allocate/deallocate

        if (pub_on_root.and.sph_harm_rot_debug) then
           write(stdout,'(a)') module_short_code//&
                " DEBUG: Unit test for "//unit_test_routine
           write(stdout,'(a,es16.8)') module_short_code//&
                " DEBUG: Absolute error threshold = ", ABS_ERR_THRESHOLD
           write(stdout,'(a,es16.8)') module_short_code//&
                " DEBUG: Relative error threshold = ", REL_ERR_THRESHOLD
        end if

        ! Assume everything is fine, unless it's not
        outcome = TEST_PASSED

        do irot = 0, NROT_REF-1
           ! Get alpha, beta and gamma angles for rotation set
           alpha_ang = RMAT_AB_EULER_ANGLES_REF(irot*3+1)
           beta_ang  = RMAT_AB_EULER_ANGLES_REF(irot*3+2)
           gamma_ang = RMAT_AB_EULER_ANGLES_REF(irot*3+3)
           do il = 0, LMAX_REF
              do iq = 1, NQ_REF
                 ! Get qmax value
                 qmax = RMAT_AB_QMAX_REF(iq)

                 ! Allocate array for number of q-values per l-value
                 ! (all identical in this case)
                 allocate(num_q_per_l(0:il),stat=ierr_alloc)
                 call utils_alloc_check(myself,"num_q_per_l",ierr_alloc)

                 ! Initialize array of number of q-values per l-value
                 num_q_per_l(0:il) = qmax


                 if (sph_harm_rot_debug.and.pub_on_root) then
                    write(stdout,'(a,i0,a,i0)') module_short_code//" DEBUG: &
                         &Testing applying R matrix with lmax = ",il,", qmax = ",qmax
                    write(stdout,'(a,f10.4,a,f10.4,a,f10.4)') &
                         module_short_code//" DEBUG: &
                         &alpha = ",alpha_ang,", &
                         &beta  = ",beta_ang,", &
                         &gamma = ",gamma_ang
                 end if

                 ! Compute dimension of (square) Mmat, i.e.
                 ! nmat = qmax*( \sum_{l=0}^{il} 2*l+1 )
                 !      = qmax*(il*(il+2)+1)
                 nmat = qmax*(il*(il+2)+1)
                 ! <-- qmax is the same for all l-values, we can use
                 !     this the product of qmax with the total number of
                 !     ang mom components, rather than summing over products
                 !     for each entry in num_q_per_l

                 ! Allocate Mmat (matrix in basis of RSHs to be rotated)
                 allocate(Mmat(nmat,nmat),stat=ierr_alloc)
                 call utils_alloc_check(myself,"Mmat",ierr_alloc)

                 ! Allocate workspace for dense matrix operations on Mmat
                 call dense_create(Mmat_dem,nmat,nmat,iscmplx=.false.)

                 ! Allocate workspace for real identity matrix
                 call dense_create(Imat_dem,nmat,nmat,iscmplx=.false.)

                 ! Allocate workspace for reference Mmat
                 call dense_create(Mmat_dem_ref,nmat,nmat,iscmplx=.false.)


                 ! [ Check per-atomblock RSH rotation matrices are unitary ]
                 ! Check that the rotated Mmat returned remains equal to the
                 ! identity. This checks whether the RSH rotation matrices are
                 ! unitary, since M' = R* M R becomes M' = R* R when M == I.

                 ! Set Mmat equal to the identity matrix
                 Mmat(:,:) = 0.0_DP
                 do icol = 1, nmat
                    Mmat(icol,icol) = 1.0_DP
                 end do

                 ! Apply RSH rotation matrices to Mmat on both sides, i.e.
                 ! M' = R* M R
                 call sph_harm_rot_apply_RSH_Rmat_to_atomblock(Mmat, &
                      il,num_q_per_l,alpha_ang,beta_ang,gamma_ang)

                 ! Put local Mmat (identical on each MPI rank) into DEM
                 do icol = 1, nmat
                    call dense_put_col(Mmat(:,icol),Mmat_dem,icol)
                 end do

                 ! Build real identity matrix Imat
                 do icol = 1, nmat
                    ! Contents are always zeroed by dense_create, so just place
                    ! 1.0_DP along diagonal
                    call dense_put_element(1.0_DP,Imat_dem,icol,icol)
                 end do

                 ! ... check for equality of Mmat with identity
                 call dense_axpy(Mmat_dem,Imat_dem,-1.0_DP)

                 ! ... using maximum column sum ("1") norm
                 norm = dense_norm("1",Mmat_dem)

                 ! ... and comparing to err_thr
                 if (sph_harm_rot_debug.and.pub_on_root) then
                    write(stdout,'(a,es16.8)') module_short_code//" DEBUG: &
                         &* 1-norm of R* I R - I                       = ",norm
                 end if
                 if (norm > ABS_ERR_THRESHOLD) then
                    outcome = TEST_FAILED
                    if (sph_harm_rot_debug.and.pub_on_root) then
                       write(stdout,'(a)') module_short_code//" DEBUG: &
                            &Test failure!"
                       write(stdout,'(a)') module_short_code//" DEBUG: &
                            & -->  1-norm of R* I R - I not within absolute &
                            &error threshold."
                    end if
                 end if

                 ! Deallocate workspace for dense matrix operations
                 call dense_destroy(Mmat_dem)

                 ! Deallocate workspace for real identity matrix
                 call dense_destroy(Imat_dem)

                 ! Deallocate workspace for reference Mmat
                 call dense_destroy(Mmat_dem_ref)

                 ! [ Check the reverse rotations restore original Mmat ]

                 ! Allocate Mmat_ref (arbitrary matrix to be rotated)
                 allocate(Mmat_ref(nmat,nmat),stat=ierr_alloc)
                 call utils_alloc_check(myself,"Mmat_ref",ierr_alloc)

                 ! Initialize arbitrary Mmat from pre-generated random vectors
                 do irow = 1, nmat
                   do icol = 1, nmat
                      Mmat(irow,icol) = RMAT_AB_RAND_VEC_REF(irow,1)*&
                           RMAT_AB_RAND_VEC_REF(icol,2)
                   end do
                 end do
                 ! Note that this arbitrary Mmat does not have the same
                 ! symmetry properties as the V-matrix which this routine will
                 ! be used to rotate. However, it should serve to test that the
                 ! per-atomblock RSH rotation matrices are actually behaving as
                 ! rotation matrices.

                 ! Copy initial Mmat to Mmat_ref for later comparison against
                 ! Mmat after forward and backward rotation
                 Mmat_ref(:,:) = Mmat(:,:)

                 ! Apply RSH rotation matrices to Mmat on both sides, i.e.
                 ! M' = R* M R
                 call sph_harm_rot_apply_RSH_Rmat_to_atomblock(Mmat, &
                      il,num_q_per_l,alpha_ang,beta_ang,gamma_ang)

                 ! Apply reverse RSH rotation, using
                 !   [R^{a}(alpha,beta,gamma)]^{-1} = R^{a}(-gamma,-beta,-alpha)
                 ! See e.g. Eq. 32{a,b,c} of Morrison & Parker, 1987
                 call sph_harm_rot_apply_RSH_Rmat_to_atomblock(Mmat, &
                      il,num_q_per_l,-gamma_ang,-beta_ang,-alpha_ang)

                 ! Compare element-by-element using internal_isclose (which uses
                 ! both relative and absolute error thresholds)
                 if (sph_harm_rot_debug.and.pub_on_root) then
                    ! Get maximum absolute difference between matrices
                    max_abs_diff = maxval(abs(Mmat(:,:)-Mmat_ref(:,:)))
                    write(stdout,'(a,es16.8)') module_short_code//" DEBUG: &
                         &* max(abs(R* M R - (R*)^-1 R* M R (R)^{-1})) = ",&
                         max_abs_diff
                 end if
                 if (.not.all(internal_isclose(Mmat(:,:),Mmat_ref(:,:),&
                     REL_ERR_THRESHOLD,ABS_ERR_THRESHOLD))) then
                    outcome = TEST_FAILED
                    if (sph_harm_rot_debug.and.pub_on_root) then
                       write(stdout,'(a)') module_short_code//" DEBUG: &
                            &Test failure!"
                       write(stdout,'(a)') module_short_code//" DEBUG: &
                            & --> Difference between element(s) of R* M R and &
                            &(R*)^-1 R* M R (R)^{-1} fails test against error &
                            &thresholds."
                    end if
                 end if

                 ! Daellocate Mmat_ref
                 deallocate(Mmat_ref,stat=ierr_alloc)
                 call utils_dealloc_check(myself,"Mmat_ref",ierr_alloc)

                 ! Deallocate Mmat
                 deallocate(Mmat,stat=ierr_alloc)
                 call utils_dealloc_check(myself,"Mmat",ierr_alloc)

                 ! Deallocate array for qmax values
                 deallocate(num_q_per_l,stat=ierr_alloc)
                 call utils_dealloc_check(myself,"num_q_per_l",ierr_alloc)

              end do
           end do
        end do

        ! [ Check application of RSH rotation matrix to small example V-matrix ]
        ! Compute dimension of (square) Mmat, i.e.

        ! Get information about example V-matrix and rotation
        nmat = NMAT_SMALL_REF
        il   = LMAX_SMALL_REF
        qmax = QMAX_SMALL_REF
        alpha_ang = RMAT_AB_EULER_ANGLE_SMALL_REF(1)
        beta_ang  = RMAT_AB_EULER_ANGLE_SMALL_REF(2)
        gamma_ang = RMAT_AB_EULER_ANGLE_SMALL_REF(3)

        ! Allocate array for number of q-values per l-value
        ! (all identical in this case)
        allocate(num_q_per_l(0:il),stat=ierr_alloc)
        call utils_alloc_check(myself,"num_q_per_l",ierr_alloc)

        ! Initialize array of number of q-values per l-value
        num_q_per_l(0:il) = qmax

        if (sph_harm_rot_debug.and.pub_on_root) then
           write(stdout,'(a,i0,a,i0)') module_short_code//" DEBUG: &
                &Testing applying R matrix to example V-matrix with &
                &lmax = ",il,", qmax = ",qmax
           write(stdout,'(a,f10.4,a,f10.4,a,f10.4)') &
                module_short_code//" DEBUG: &
                &alpha = ",alpha_ang,", &
                &beta  = ",beta_ang,", &
                &gamma = ",gamma_ang
        end if

        ! Allocate Mmat (matrix in basis of RSHs to be rotated)
        allocate(Mmat(nmat,nmat),stat=ierr_alloc)
        call utils_alloc_check(myself,"Mmat",ierr_alloc)

        ! Allocate Mmat_ref (arbitrary matrix to be rotated)
        allocate(Mmat_ref(nmat,nmat),stat=ierr_alloc)
        call utils_alloc_check(myself,"Mmat_ref",ierr_alloc)

        ! Copy initial example V-matrix atomblock (before rotation) to Mmat
        Mmat(1:nmat,1:nmat) = RMAT_AB_VMAT_SMALL_BEFORE_REF(1:nmat,1:nmat)

        ! Copy reference example V-matrix atomblock (after rotation) to Mmat_ref
        Mmat_ref(1:nmat,1:nmat) = RMAT_AB_VMAT_SMALL_AFTER_REF(1:nmat,1:nmat)

        ! Apply RSH rotation matrices to Mmat on both sides, i.e.
        ! M' = R* M R
        call sph_harm_rot_apply_RSH_Rmat_to_atomblock(Mmat, &
             il,num_q_per_l,alpha_ang,beta_ang,gamma_ang)

        ! Compare element-by-element using internal_isclose (which uses
        ! both relative and absolute error thresholds)
        if (sph_harm_rot_debug.and.pub_on_root) then
           ! Get maximum absolute difference between matrices
           max_abs_diff = maxval(abs(Mmat(:,:)-Mmat_ref(:,:)))
           write(stdout,'(a,es16.8)') module_short_code//" DEBUG: &
                &* max(abs(R* V R - V'ref))                   = ",&
                max_abs_diff
        end if
        if (.not.all(internal_isclose(Mmat(:,:),Mmat_ref(:,:),&
            REL_ERR_THRESHOLD,ABS_ERR_THRESHOLD))) then
           outcome = TEST_FAILED
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   &Test failure!"
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   & --> Difference between element(s) of R* V R and &
                   &V'ref fails test against error thresholds."
           end if
        end if

        ! Daellocate Mmat_ref
        deallocate(Mmat_ref,stat=ierr_alloc)
        call utils_dealloc_check(myself,"Mmat_ref",ierr_alloc)

        ! Deallocate Mmat
        deallocate(Mmat,stat=ierr_alloc)
        call utils_dealloc_check(myself,"Mmat",ierr_alloc)

        ! Deallocate array for qmax values
        deallocate(num_q_per_l,stat=ierr_alloc)
        call utils_dealloc_check(myself,"num_q_per_l",ierr_alloc)

      end subroutine internal_test_sph_harm_rot_apply_RSH_Rmat_to_atomblock

      subroutine internal_test_jacobi_polynomial_realx(outcome)
        !======================================================================!
        ! # INTERNAL PROCEDURE DESCRIPTION                                     !
        ! Host:      sph_harm_rot_run_unit_tests                               !
        ! Function:  Unit test for jacobi_polynomial_real_x function           !
        !----------------------------------------------------------------------!
        ! ## TESTS                                                             !
        ! * Check polynomial satisfies special values                          !
        ! * Check polynomial agrees with precomputed reference values with     !
        !   l up to LMAX_REF and x values in [-1,1] (i.e. [cos(pi),cos(0)])    !
        !   (l, m, mu values are those allowed in definition of Wigner small-d !
        !   matrix in terms of Jacobi polynomials).                            !
        !----------------------------------------------------------------------!
        ! # ARGUMENTS                                                          !
        ! <name>      <in/out/inout> <arg descrption>                          !
        ! outcome     out            integer outcome code for this unit test   !
        !----------------------------------------------------------------------!
        ! # AUTHORS & CHANGELOG                                                !
        ! Author(s): James C. Womack                                           !
        ! Date of creation: 10/2018                                            !
        ! List of major changes:                                               !
        ! <date> <change description> <author>                                 !
        !======================================================================!
        use constants, only: PI
        implicit none

        ! Parameters
        character(len=*), parameter :: unit_test_routine = &
             "jacobi_polynomial_realx"
        integer, parameter :: NX = 100
        ! --> number x values to sample when testing special cases
        integer, parameter :: ABMAX = 3
        ! --> maximum a, b value for testing special cases
        integer, parameter :: LMAX_REF = 2
        integer, parameter :: NX_REF = 5
        real(kind=DP), parameter :: ERR_FACTOR = 1.0e2_DP
        real(kind=DP), parameter :: ERR_THRESHOLD = ERR_FACTOR * EPSILON(1.0_DP)
        ! --> rough error threshold for unit testing based on machine epsilon
        !     (to adjust, change ERR_FACTOR)

        ! Reference data (do not modify)

        ! Arbitary beta angles from which to obtain x argument for
        ! Jacobi polynomial (x = cos(beta)) in evaluating reference
        ! data
        real(kind=DP) :: JACOBI_P_BETA_REF(NX_REF)

        data JACOBI_P_BETA_REF(1:NX_REF) &
             / 1.8757992931910591e-01_DP, &
               4.7098488884217804e-01_DP, &
               7.7217718056222751e-01_DP, &
               2.8273984957925071e+00_DP, &
               1.1056138764903252e+00_DP /

        ! Store matrices sequentially in 1-D array up to LMAX_REF
        !
        ! l = 0 (1x1 matrix, 1 element)
        ! l = 1 (3x3 matrix, 9 elements)
        ! l = 2 (5x5 matrix, 25 elements)
        !
        ! Pre-computed values from Sage Math's built-in jacobi_P function
        ! evaluated as a component of the Wigner small-d matrix, i.e.
        ! The Wigner small-d matrix has the form
        !   d_{mu m}^{l}(beta) = [...] P^(m-mu,m+mu)_{l-m}(cos(beta))
        ! where [...] is an expression in l, mu, m, and beta and P^(a,b)_{n}(x)
        ! is a Jacobi polynomial and we are using as reference data
        !   M_{mu m}^{l}(beta) = P^(m-mu,m+mu)_{l-m}(cos(beta))
        !
        ! For each matrix, row index is mu, column index is m and mu is
        ! the fastest changing index.
        ! The reference data (for each l) are stored in order of reference
        ! beta values
        real(kind=DP) :: JACOBI_P_REF((1+9+25)*NX_REF)

        ! Beta value 0
        ! jP ref matrix for il = 0, beta = 0.1876
        data JACOBI_P_REF(1:1) &
             / 1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 1, beta = 0.1876
        data JACOBI_P_REF(2:10) &
             / 9.8253533771158374e-01_DP, &
             -8.6938677255763690e-03_DP, &
             7.6926837263542163e-05_DP, &
             1.9824584108743202e+00_DP, &
             9.8245841087432018e-01_DP, &
             -1.7541589125679824e-02_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 2, beta = 0.1876
        data JACOBI_P_REF(11:35) &
             / 9.6537568985201583e-01_DP, &
             -8.5420322617690177e-03_DP, &
             7.5583336029818443e-05_DP, &
             -6.6879174771617481e-07_DP, &
             5.9177382913714993e-09_DP, &
             3.8956708884551396e+00_DP, &
             9.4806487532038841e-01_DP, &
             -1.7082726940042599e-02_DP, &
             2.2808167384659629e-04_DP, &
             -2.6988379440301852e-06_DP, &
             5.8952120262695029e+00_DP, &
             2.9215244099580220e+00_DP, &
             9.4783679364654172e-01_DP, &
             -2.5850822664938482e-02_DP, &
             4.6156102358125298e-04_DP, &
             3.9649168217486404e+00_DP, &
             2.9649168217486404e+00_DP, &
             1.9649168217486404e+00_DP, &
             9.6491682174864035e-01_DP, &
             -3.5083178251359648e-02_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! Beta value 1
        ! jP ref matrix for il = 0, beta = 0.4710
        data JACOBI_P_REF(36:36) &
             / 1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 1, beta = 0.4710
        data JACOBI_P_REF(37:45) &
             / 8.9408542809486680e-01_DP, &
             -5.1475478507593563e-02_DP, &
             2.9636148899460403e-03_DP, &
             1.8911218132049208e+00_DP, &
             8.9112181320492079e-01_DP, &
             -1.0887818679507921e-01_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 2, beta = 0.4710
        data JACOBI_P_REF(46:70) &
             / 7.9938875273158128e-01_DP, &
             -4.6023475237849908e-02_DP, &
             2.6497248875857271e-03_DP, &
             -1.5255349457220165e-04_DP, &
             8.7830132159098817e-06_DP, &
             3.3816489118777247e+00_DP, &
             6.9939262759312426e-01_DP, &
             -9.1741843486555411e-02_DP, &
             8.2454986386856742e-03_DP, &
             -6.4534603115244617e-04_DP, &
             5.3645125685692010e+00_DP, &
             2.5278298487618196e+00_DP, &
             6.9114712895443864e-01_DP, &
             -1.4553559085294257e-01_DP, &
             1.7781689339676242e-02_DP, &
             3.7822436264098416e+00_DP, &
             2.7822436264098416e+00_DP, &
             1.7822436264098416e+00_DP, &
             7.8224362640984157e-01_DP, &
             -2.1775637359015843e-01_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! Beta value 2
        ! jP ref matrix for il = 0, beta = 0.7722
        data JACOBI_P_REF(71:71) &
             / 1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 1, beta = 0.7722
        data JACOBI_P_REF(72:80) &
             / 7.3650153916789662e-01_DP, &
             -1.2169513943666230e-01_DP, &
             2.0108181958778750e-02_DP, &
             1.7163933572091179e+00_DP, &
             7.1639335720911790e-01_DP, &
             -2.8360664279088210e-01_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 2, beta = 0.7722
        data JACOBI_P_REF(81:105) &
             / 5.4243451719668079e-01_DP, &
             -8.9628657504353573e-02_DP, &
             1.4809706962508679e-02_DP, &
             -2.4470680072913573e-03_DP, &
             4.0433898168735518e-04_DP, &
             2.5282526988041374e+00_DP, &
             3.1874808130044757e-01_DP, &
             -1.7436317899412446e-01_DP, &
             4.8918917920421398e-02_DP, &
             -1.1405627955914851e-02_DP, &
             4.4190092350073797e+00_DP, &
             1.8444191991937031e+00_DP, &
             2.6982916338002627e-01_DP, &
             -3.0476087243365063e-01_DP, &
             1.2064909175267249e-01_DP, &
             3.4327867144182358e+00_DP, &
             2.4327867144182358e+00_DP, &
             1.4327867144182358e+00_DP, &
             4.3278671441823580e-01_DP, &
             -5.6721328558176420e-01_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! Beta value 3
        ! jP ref matrix for il = 0, beta = 2.8274
        data JACOBI_P_REF(106:106) &
             / 1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 1, beta = 2.8274
        data JACOBI_P_REF(107:115) &
             / 5.9913005557003319e-04_DP, &
             -2.3878003264531673e-02_DP, &
             9.5164486341536658e-01_DP, &
             4.8954266640203414e-02_DP, &
             -9.5104573335979659e-01_DP, &
             -1.9510457333597966e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 2, beta = 2.8274
        data JACOBI_P_REF(116:140) &
             / 3.5895682348735110e-07_DP, &
             -1.4306029422780296e-05_DP, &
             5.7015903990098517e-04_DP, &
             -2.2723379155306922e-02_DP, &
             9.0562794606485175e-01_DP, &
             5.8659944985070591e-05_DP, &
             -1.7387302217250289e-03_DP, &
             4.5418146251768286e-02_DP, &
             -8.5847071063453495e-01_DP, &
             -3.7134053008806349e+00_DP, &
             3.5947803334201991e-03_DP, &
             -6.9836619626884922e-02_DP, &
             8.5673198041280996e-01_DP, &
             2.7833005804525048e+00_DP, &
             5.7098691804921993e+00_DP, &
             9.7908533280406829e-02_DP, &
             -9.0209146671959317e-01_DP, &
             -1.9020914667195932e+00_DP, &
             -2.9020914667195932e+00_DP, &
             -3.9020914667195932e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! Beta value 4
        ! jP ref matrix for il = 0, beta = 1.1056
        data JACOBI_P_REF(141:141) &
             / 1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 1, beta = 1.1056
        data JACOBI_P_REF(142:150) &
             / 5.2460025713807501e-01_DP, &
             -1.9969267886094116e-01_DP, &
             7.6014385140042537e-02_DP, &
             1.4485858719980325e+00_DP, &
             4.4858587199803257e-01_DP, &
             -5.5141412800196743e-01_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! jP ref matrix for il = 2, beta = 1.1056
        data JACOBI_P_REF(151:175) &
             / 2.7520542978933443e-01_DP, &
             -1.0475883067904078e-01_DP, &
             3.9877165990658985e-02_DP, &
             -1.5179516200582411e-02_DP, &
             5.7781867482187190e-03_DP, &
             1.5198570418735009e+00_DP, &
             -5.3943729540724031e-02_DP, &
             -1.7915862895691675e-01_DP, &
             1.4421234362492308e-01_DP, &
             -8.3830811795204518e-02_DP, &
             3.1476015428284501e+00_DP, &
             9.7472273483140170e-01_DP, &
             -1.9815607316564710e-01_DP, &
             -3.7103488116269595e-01_DP, &
             4.5608631084025519e-01_DP, &
             2.8971717439960649e+00_DP, &
             1.8971717439960651e+00_DP, &
             8.9717174399606514e-01_DP, &
             -1.0282825600393486e-01_DP, &
             -1.1028282560039349e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP, &
             1.0000000000000000e+00_DP /

        ! Arguments
        integer,intent(out) :: outcome

        ! Local variables
        integer       :: ix
        integer       :: n, a, b
        integer       :: il, m, mu
        integer       :: istart_ref
        real(kind=DP) :: x
        real(kind=DP) :: jP, jP_ref
        real(kind=DP) :: max_abs_diff, abs_diff

        if (pub_on_root.and.sph_harm_rot_debug) then
           write(stdout,'(a)') module_short_code//&
                " DEBUG: Unit test for "//unit_test_routine
           write(stdout,'(a,es16.8)') module_short_code//&
                " DEBUG: Error threshold = ", ERR_THRESHOLD
        end if

        ! Assume everything is fine, unless it's not
        outcome = TEST_PASSED

        ! [ Check against some special values ]
        ! From Abramowitz and Stegun, "Handbook of Mathematical Functions",
        ! table 22.4, p.777"
        if (sph_harm_rot_debug.and.pub_on_root) then
           write(stdout,'(a,i0)') module_short_code//" DEBUG: &
                &Testing Jacobi polynomial against special values"
        end if

        ! P^(a,b)_{0}(x) = 1
        max_abs_diff = 0.0_DP
        n = 0
        jP_ref = 1.0_DP
        ! n + a, n + b and n + a + b must be >= 0, so for n = 0, a, b >= 0
        do a = 0, ABMAX
           do b = 0, ABMAX
              do ix = 0, NX
                 ! To sample x in range [-1,1], use cos(x) for x = 0, pi
                 x = cos(PI*real(ix,kind=DP)/real(NX,kind=DP))
                 jP = jacobi_polynomial_realx(n,a,b,x)
                 abs_diff = abs(jP-jP_ref)
                 if (abs_diff > max_abs_diff) max_abs_diff = abs_diff
              end do
           end do
        end do

        if (sph_harm_rot_debug.and.pub_on_root) then
           write(stdout,'(a,es16.8)') module_short_code//" DEBUG: &
           &* Max. abs. diff for P^(a,b)_{0}(x) - 1.0          = ", max_abs_diff
        end if
        if (max_abs_diff > ERR_THRESHOLD) then
           outcome = TEST_FAILED
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   &Test failure!"
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   & --> Computed P^(a,b)_{0}(x) not within error &
                   &threshold of expected value."
           end if
        end if

        ! P^(a,b)_{1}(x) = (1/2)(a-b+(a+b+2)x)
        max_abs_diff = 0.0_DP
        n = 1
        ! n + a, n + b and n + a + b must be >= 0, so for n = 1
        ! a >= -n, b >= max(-n-a,-n)
        do a = -n, ABMAX
           do b = max(-n-a,-n), ABMAX
              do ix = 0, NX
                 ! To sample x in range [-1,1], use cos(x) for x = 0, pi
                 x = cos(PI*real(ix,kind=DP)/real(NX,kind=DP))
                 jP_ref = 0.5_DP*(a-b+(a+b+2)*x)
                 jP = jacobi_polynomial_realx(n,a,b,x)
                 abs_diff = abs(jP-jP_ref)
                 if (abs_diff > max_abs_diff) max_abs_diff = abs_diff
              end do
           end do
        end do

        if (sph_harm_rot_debug.and.pub_on_root) then
           write(stdout,'(a,es16.8)') module_short_code//" DEBUG: &
           &* Max. abs. diff for P^(a,b)_{1}(x) - (known form) = ", max_abs_diff
        end if
        if (max_abs_diff > ERR_THRESHOLD) then
           outcome = TEST_FAILED
           if (sph_harm_rot_debug.and.pub_on_root) then
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   &Test failure!"
              write(stdout,'(a)') module_short_code//" DEBUG: &
                   & --> Computed P^(a,b)_{1}(x) not within error &
                   &threshold of expected value."
           end if
        end if

        ! [ Check computed Jacobi polynomial against reference results ]
        ! Use integer arguments possible in Wigner small-d matrix, as defined
        ! in Eq. B9 of Morrison & Parker, Aus. J. Phys. 1987:
        !   d_{mu m}^{l}(beta) = [...] P^(m-mu,m+mu)_{l-m}(cos(beta))
        ! where [...] is an expression in l, mu, m, and beta and P^(a,b)_{n}(x)
        ! is a Jacobi polynomial.
        !
        ! For each value of l:
        !   m = -l, ..., +l
        !  mu = -l, ..., +l
        !
        ! beta is in the range 0 <= beta <= pi.
        !
        ! The arguments for jacobi_polynomial_realx in terms of l, m, mu are
        !   n = l - m
        !   a = m - mu
        !   b = m + mu

        ! Counter for position in 1-D reference data array
        istart_ref = 1

        ! Loop over pre-set random x values (from x = cos(beta))
        do ix = 1, NX_REF
           x = cos(JACOBI_P_BETA_REF(ix))
           ! Loop over total angular momentum
           do il = 0, LMAX_REF ! n is j from d_{mu m}^{j}(beta)
              ! Check abs difference against tolerance per l value
              if (sph_harm_rot_debug.and.pub_on_root) then
                 write(stdout,'(a,i0,a,f10.4)') module_short_code//" DEBUG: &
                      &Testing Jacobi polynomials with l = ",il,", x  = ", x
              end if
              max_abs_diff = 0.0_DP
              ! Loop over m and mu values corresponding to il
              do m = -il, il
                 do mu = -il, il
                    ! Compute integer arguments for Jacobi polynomial
                    n = il - m
                    a = m - mu
                    b = m + mu
                    ! Check against reference result for l, m, mu, x
                    ! where x = cos(beta) and beta is an arbitrary number
                    ! between 0 and pi.
                    jP_ref = JACOBI_P_REF(istart_ref)
                    jP = jacobi_polynomial_realx(n,a,b,x)
                    abs_diff = abs(jP-jP_ref)
                    if (abs_diff > max_abs_diff) max_abs_diff = abs_diff
                    ! We can simply increment a single counter to step through
                    ! reference data, since this is ordered
                    !   mu, m, il, ix
                    ! from fastest to slowest changing indices.
                    istart_ref = istart_ref + 1
                 end do
              end do
              if (sph_harm_rot_debug.and.pub_on_root) then
                 write(stdout,'(a,es16.8)') module_short_code//" DEBUG: &
                 &* Max. abs. diff for P^(m-mu,m+mu)_{l-m}(x) - ref  = ", &
                 max_abs_diff
              end if
              if (max_abs_diff > ERR_THRESHOLD) then
                 outcome = TEST_FAILED
                 if (sph_harm_rot_debug.and.pub_on_root) then
                    write(stdout,'(a)') module_short_code//" DEBUG: &
                         &Test failure!"
                    write(stdout,'(a)') module_short_code//" DEBUG: &
                         & --> Computed P^(m-mu,m+mu)_{l-m}(x) not within &
                         &error threshold of expected value."
                 end if
              end if
           end do
        end do

      end subroutine internal_test_jacobi_polynomial_realx

  end subroutine sph_harm_rot_run_unit_tests

  ! ## PRIVATE PROCEDURES

  elemental function jacobi_polynomial_realx(n,a,b,x) result(P)
    !==========================================================================!
    ! # DESCRIPTION                                                            !
    ! Evaluate the Jacobi polynomial P_^{(a,b)}_{n}(x) for real x in [-1,1]    !
    ! and integers n, a, b satisfying the following inequalities               !
    !   n + a >= 0                                                             !
    !   n + b >= 0                                                             !
    !   n + a + b >= 0                                                         !
    !                                                                          !
    ! In this case, the Jacobi polynomial has the following definition:        !
    !   P^(a,b)_{n}(x)                                                         !
    !     = f1(n,a,b) * \sum_{s} f2(s,n,a,b) * f3(x)^(n-s) * f4(x)^s           !
    ! with                                                                     !
    !   f1 = (n+a)!(n+b)!                                                      !
    !   f2 = 1/(s!(n+a-s)!(b+s)!(n-s)!)                                        !
    !   f3 = (x-1)/2                                                           !
    !   f4 = (x+1)/2                                                           !
    ! and where the summation over s is over all integers for which the        !
    ! arguments of the factorials are >=0.                                     !
    !                                                                          !
    ! TODO It is clear from the definition why n+a and n+b must be positive    !
    !      (avoid factorial for negative number), but why must n+a+b be        !
    !      positive (as stated on Wikipedia)?                                  !
    !--------------------------------------------------------------------------!
    ! # ARGUMENTS AND RESULT                                                   !
    ! <name> <in/out/inout> <arg descrption>                                   !
    ! n      in             integer n satisfying above inequalities            !
    ! a      in             integer a satisfying above inequalities            !
    ! b      in             integer b satisfying above inequalities            !
    ! x      in             real number, -1 <= x <= 1                          !
    ! P      out (result)   value of Jacobi polynomial P^(a,b)_{n}(x)          !
    !--------------------------------------------------------------------------!
    ! # AUTHORS & CHANGELOG                                                    !
    ! Author(s):        James C. Womack                                        !
    ! Date of creation: 10/2018                                                !
    ! List of major changes:                                                   !
    ! <date> <change description> <author>                                     !
    !==========================================================================!

    use services, only: services_factorial

    implicit none

    ! Parameters
    !character(len=*), parameter :: myself = "jacobi_polynomial_realx"

    ! Arguments
    integer, intent(in)           :: n
    integer, intent(in)           :: a
    integer, intent(in)           :: b
    real(kind=DP), intent(in)     :: x

    ! Result
    real(kind=DP) :: P

    ! Local variables
    integer :: s
    integer :: max_s ! maximum value in summation over s
    integer :: min_s ! minimum value in summation over s
    real(kind=DP) :: f1, f2, f3, f4
    ! <-- parts of expression for polynomial (see above)

    !---------------------------------------------------------------------------
    ! @optimize Precomputed factorials
    ! We may be able to optimize the routine by precomputing factorials.
    ! This will require benchmarking to determine whether it is more efficient.
    ! For now, we use naive direct computation of factorials with the
    ! services_factorial function.
    !
    ! We can determine the maximum argument for a factorial, given n, a and b.
    ! This enables the factorials up to a certain argument value to be
    ! precomputed.
    !
    ! Factorials to evaluate in the Jacobi polynomial expression:
    !   (n+a)!, (n+b)!, s!, (n+a-s)!, (b+s)!, (n-s)!
    ! but max_s <= n+a, max_s <= n and min_s >= 0
    ! so s     <= n+a
    !    n+a-s <= n+a
    !    b+s   <= n+b
    ! so we only need to compare
    !   (n+a)!, (n+b)!, (n-max_s)!
    !---------------------------------------------------------------------------

    ! Determine minimum and maximum values of s
    ! From Morrison & Parker:
    ! """
    ! the summation over s extends over all integers for which the arguments of
    ! factorials are non-negative
    ! """
    ! Maximum value: s <= n+a, s <= n
    max_s = min(n+a,n)
    ! Minimum value: s >= -b,  s >= 0
    min_s = max(-b,0)
    ! Provided s is in these ranges, all factorials in f2 will be >= 0:
    ! * max_s is minimum of n+a and n, so (n+a-s) and (n-s) will always be
    !   positive or zero
    ! * min_s is maximum of -b and 0, so (b+s) will always be positive or zero

    ! Compute f1, f3, f4 (independent of s)
    f1 = real(services_factorial(n+a) * services_factorial(n+b),kind=DP)
    f3 = (x-1.0_DP)/2.0_DP
    f4 = (x+1.0_DP)/2.0_DP

    P = 0.0_DP

    do s = min_s, max_s
      f2 = 1.0_DP/real(services_factorial(s) * services_factorial(n+a-s) &
                   * services_factorial(b+s) * services_factorial(n-s),kind=DP)
      P  = P + f2 * f3**(n-s) * f4**s
    end do

    P = f1 * P

  end function jacobi_polynomial_realx

end module sph_harm_rotation
