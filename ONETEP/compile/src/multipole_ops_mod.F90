! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!-----------------------------------------------------------------------!
!                     M U L T I P O L E   O P S                         !
!-----------------------------------------------------------------------!
!                                                                       !
! This module implements a number of simple multipole conversions and   !
! operations (such as translation) to keep the mundane juggling of      !
! conventions in one place.                                             !
!                                                                       !
! Also implemented are: Coulombic and Thole interaction tensors for     !
! orders 0 to 4, and potentials of sets of spherical multipoles.        !
!                                                                       !
! This module is unaware of periodic boundary conditions or minimum     !
! image conventions. OBCs are assumed everywhere.                       !
!                                                                       !
! Written by Jacek Dziedzic in 2014.09.                                 !
! Extended by Jacek Dziedzic in 2015.05 to handle MDQPOLES.             !
! Extended by Jacek Dziedzic in 2016.09 to handle T3 and T4.            !
! Extended by Jacek Dziedzic in 2017.07 with multipole_set_info(),      !
!      multipole_set_merge(), multipole_set_gather(),                   !
!      multipole_set_scale_monopoles().                                 !
!                                                                       !
!-----------------------------------------------------------------------!
!                                                                       !
! [1] "Theory of Molecular Fluids Volume1: Fundamentals"                !
!     by C.G. Gray and K.E. Gubbins.                                    !
! [2] "Jackson notes 2014. 1. Spherical multipole moments.              !
!     www.physics.sfsu.edu/~lea/courses/grad/sphermult.PDF              !
! [3] "Classical electrodynamics", 3rd ed.                              !
!     by John David Jackson.                                            !
! [4] "Distributed multipole expansion of the charge density in the     !
!     ONETEP code" by Chris-Kriton Skylaris.                            !
! [5] JD's notes on multipoles.                                         !
! [6] "Distributed Multipole Analysis of Gaussian wavefunctions"        !
!     GDMA version 2.2.09 by Anthony Stone                              !
!                                                                       !
!-----------------------------------------------------------------------!

!@todo: Handle PBCs correctly

module multipole_ops

  use constants, only: DP
  use geometry, only: POINT

  implicit none

  private

  public :: multipole_0_sph_real_to_cart

  public :: multipole_1_sph_real_to_cart
  public :: multipole_1_cart_to_sph_real
  public :: multipole_1_translate_cartesian
  public :: multipole_1_magnitude_cartesian
  public :: multipole_1_magnitude_sph_real

  public :: multipole_2_sph_real_to_cart_traceless
  public :: multipole_2_cart_traceless_to_sph_real
  public :: multipole_2_cart_traceless_to_primitive
  public :: multipole_2_cart_primitive_to_traceless
  public :: multipole_2_translate_cartesian_primitive
  public :: multipole_2_magnitude_cartesian
  public :: multipole_2_magnitude_sph_real

  public :: multipole_init_spherical_set
  public :: multipole_free_spherical_set
  public :: multipole_pot_of_spherical_set

  public :: multipole_prepare_gamma_lookup
  public :: multipole_destroy_gamma_lookup

  public :: multipole_thole_damping
  public :: multipole_T0
  public :: multipole_T1
  public :: multipole_T2
  public :: multipole_T3
  public :: multipole_T4

  public :: multipole_mdq_cart_to_spherical_set
  public :: multipole_set_info
  public :: multipole_set_merge
  public :: multipole_set_gather
  public :: multipole_set_scale_monopoles
  public :: multipole_set_scale_quadrupoles

  public :: CART_PRIM_MDQPOLE
  public :: SPHERICAL_MULTIPOLE_SET

  ! jd: A monopole, dipole and quadrupole located at a centre.
  !     The dipole and quadrupole are in a Cartesian primitive representation.
  type CART_PRIM_MDQPOLE
     type(POINT)      :: centre
     real(kind=DP)    :: charge
     real(kind=DP)    :: dipole(3)
     real(kind=DP)    :: quadrupole(6)
     real(kind=DP)    :: polarisability
     real(kind=DP)    :: userdata(3)
     character(len=4) :: species
     logical          :: masked ! the mdqpole should be ignored if .true.
  end type CART_PRIM_MDQPOLE

  ! jd: The number of elements in the above (for comms)
  integer, parameter, public :: CART_PRIM_MDQPOLE_NDATA = 18

  integer, parameter, public :: DAMPING_COULOMBIC_MASKED = 1
  integer, parameter, public :: DAMPING_COULOMBIC_SMEARED = 2
  integer, parameter, public :: DAMPING_THOLE = 3
  integer, parameter, public :: DAMPING_ZERO = 4
  integer, parameter, public :: DAMPING_FROM_POTENTIAL = 5

  integer, parameter, public :: set_name_len = 80

  ! jd: A set of multipoles on n_sites sites. Multipoles are in real spherical
  !     representation. There is one of each (l,m) for every site.
  type SPHERICAL_MULTIPOLE_SET
     integer                    :: min_l
     integer                    :: max_l
     integer                    :: n_multipoles_per_site ! == (1+max_l-min_l)*(1+max_l+min_l)
     integer                    :: n_sites
     integer                    :: potential_damping
     integer                    :: energy_damping
     character(len=set_name_len):: set_name
     real(kind=DP)              :: dma_monopole_scaling_factor
     real(kind=DP)              :: dma_nelecs
     type(POINT), allocatable   :: centres(:) ! size: n_sites
     integer, allocatable       :: min_l_mask(:) ! size: n_sites
     integer, allocatable       :: max_l_mask(:) ! size: n_sites
     real(kind=DP), allocatable :: safe_radii(:) ! size: n_sites
     real(kind=DP), allocatable :: polarisabilities(:) ! size: n_sites
     real(kind=DP), allocatable :: mpoles(:)! size: n_sites*n_multipole_per_site
     real(kind=DP), allocatable :: userdata(:,:) ! size: 1:3, n_sites
     character(len=4), allocatable :: species(:) ! size: n_sites
  end type
  ! Note on safe_radii. If the multipole set represents a distributed density,
  ! such as due to NGWFs, the multipole expansion is strictly convergent only
  ! outside of the density, i.e. more than safe_radius away from the centre.
  ! This is checked against when generating the potential of the multipoles.
  ! If the multipole set represents a point density, this is set to zero.

  real(kind=DP), dimension(3), parameter     :: rank1zero = &
       (/0D0,0D0,0D0/)

  real(kind=DP), dimension(3,3), parameter   :: rank2zero = &
       reshape((/0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/),(/3,3/))

  real(kind=DP), dimension(3,3,3), parameter :: rank3zero = &
       reshape((/&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/),(/3,3,3/))

  real(kind=DP), dimension(3,3,3,3), parameter :: rank4zero = &
       reshape((/&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/),(/3,3,3,3/))

  ! jd: A lookup table of \Gamma(2/3,x) needed when calculating Thole-0-damped
  !     potentials. This is allocated and filled with a call to
  !     multipole_prepare_gamma_lookup() and destroyed with a call to
  !     multipole_destroy_gamma_lookup().
  real(kind=DP), allocatable, save :: gamma_two_thirds_lookup(:)

contains

  ! ---------------------------------------------------------------------------
  ! ----- P U B L I C   S U B R O U T I N E S   A N D   F U N C T I O N S -----
  ! ---------------------------------------------------------------------------

  subroutine multipole_prepare_gamma_lookup
    !==========================================================================!
    ! Prepares a module-wide lookup for \Gamma(2/3,x) needed when calculating  !
    ! Thole-0-damped potentials.                                               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.07.                                    !
    !==========================================================================!

    use constants, only: TWO_THIRDS
    use services, only: services_incomplete_gamma_prepare_lookup
    use timer, only: timer_clock

    implicit none

    ! ---------------------------------------------------------------------------

    call timer_clock('multipole_ops_gamma_lookup',1)
    call services_incomplete_gamma_prepare_lookup(gamma_two_thirds_lookup, &
         TWO_THIRDS)
    call timer_clock('multipole_ops_gamma_lookup',2)

  end subroutine multipole_prepare_gamma_lookup

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_destroy_gamma_lookup
    !==========================================================================!
    ! Destroyes a module-wide lookup for \Gamma(2/3,x) needed when calculating !
    ! Thole-0-damped potentials.                                               !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.07.                                    !
    !==========================================================================!

    use services, only: services_incomplete_gamma_destroy_lookup

    implicit none

    ! ---------------------------------------------------------------------------

    call services_incomplete_gamma_destroy_lookup(gamma_two_thirds_lookup)

  end subroutine multipole_destroy_gamma_lookup

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_0_sph_real_to_cart(q, &             ! output
       Q00, is_racah_normalized)                           ! input
    !==========================================================================!
    ! Converts spherical real monopole to Cartesian monopole.                  !
    ! The input monopole may or may not be Racah-normalized -- set the last    !
    ! argument appropriately.                                                  !
    ! Reference: [1] eq. (2.85)                                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    use constants, only: PI

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: q
    real(kind=DP), intent(in)  :: Q00
    logical, intent(in)        :: is_racah_normalized

    ! jd: Local variables
    real(kind=DP), parameter   :: INV_RACAH_FACTOR = SQRT(4.0_DP*PI/1.0_DP)

    ! ---------------------------------------------------------------------------

    q = Q00

    if(.not. is_racah_normalized) then
       q = q * INV_RACAH_FACTOR
    end if

  end subroutine multipole_0_sph_real_to_cart

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_1_sph_real_to_cart(mux, muy, muz, & ! output
       Q1m1, Q10, Q11, is_racah_normalized)                ! input
    !==========================================================================!
    ! Converts spherical real dipoles to Cartesian dipoles.                    !
    ! The input dipoles may or may not be Racah-normalized -- set the last     !
    ! argument appropriately.                                                  !
    ! Reference: [1] eq. (2.86)                                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    use constants, only: PI

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: mux, muy, muz
    real(kind=DP), intent(in)  :: Q1m1, Q10, Q11
    logical, intent(in)        :: is_racah_normalized

    ! jd: Local variables
    real(kind=DP), parameter   :: INV_RACAH_FACTOR = SQRT(4.0_DP*PI/3.0_DP)

    ! ---------------------------------------------------------------------------

    mux = Q11
    muy = Q1m1
    muz = Q10

    if(.not. is_racah_normalized) then
       mux = mux * INV_RACAH_FACTOR
       muy = muy * INV_RACAH_FACTOR
       muz = muz * INV_RACAH_FACTOR
    end if

  end subroutine multipole_1_sph_real_to_cart

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_2_sph_real_to_cart_traceless(&
       Qxx, Qxy, Qxz, Qyy, Qyz, Qzz, &                     ! output
       Q2m2, Q2m1, Q20, Q21, Q22, is_racah_normalized)     ! input
    !==========================================================================!
    ! Converts spherical real quadrupoles to Cartesian traceless quadrupoles.  !
    ! The input quadrupoles may or may not be Racah-normalized -- set the last !
    ! argument appropriately.                                                  !
    ! Reference: [1] eq. (2.87)                                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    use constants, only: PI

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: Qxx, Qxy, Qxz, Qyy, Qyz, Qzz
    real(kind=DP), intent(in)  :: Q2m2, Q2m1, Q20, Q21, Q22
    logical, intent(in)        :: is_racah_normalized

    ! jd: Local variables
    real(kind=DP), parameter   :: SQRT3_OVER_2 = SQRT(3.0_DP)/2.0_DP
    real(kind=DP), parameter   :: INV_RACAH_FACTOR = SQRT(4.0_DP*PI/5.0_DP)

    ! ---------------------------------------------------------------------------

    Qxy = SQRT3_OVER_2 * Q2m2
    Qxz = SQRT3_OVER_2 * Q21
    Qyz = SQRT3_OVER_2 * Q2m1
    Qzz = Q20
    Qxx = SQRT3_OVER_2 * Q22 - 0.5_DP*Q20
    Qyy = -Q20 - Qxx

    if(.not. is_racah_normalized) then
       Qxx = Qxx * INV_RACAH_FACTOR
       Qxy = Qxy * INV_RACAH_FACTOR
       Qxz = Qxz * INV_RACAH_FACTOR
       Qyy = Qyy * INV_RACAH_FACTOR
       Qyz = Qyz * INV_RACAH_FACTOR
       Qzz = Qzz * INV_RACAH_FACTOR
    end if

  end subroutine multipole_2_sph_real_to_cart_traceless

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_2_cart_traceless_to_primitive(&
       Qxx_p, Qxy_p, Qxz_p, Qyy_p, Qyz_p, Qzz_p, &      ! output
       Qxx_t, Qxy_t, Qxz_t, Qyy_t, Qyz_t, Qzz_t, trace) ! input
    !==========================================================================!
    ! Converts Cartesian traceless quadrupoles to Cartesian primitive.         !
    ! The convention of Ref. [1] (GrayGubbins) is used for the traceless       !
    ! quadrupole definition.                                                   !
    ! Racah-normalization has no effect here.                                  !
    ! There is an arbitrariness in the resultant primitive Cartesians, since   !
    ! they are invariant wrt adding \lambda * I. The last parameter, 'trace'   !
    ! specifies the desired trace of the resultant primitive Cartesian that    !
    ! removes the arbitrariness.                                               !
    ! Reference: [1] eq. (2.59)                                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: Qxx_p, Qxy_p, Qxz_p, Qyy_p, Qyz_p, Qzz_p
    real(kind=DP), intent(in)  :: Qxx_t, Qxy_t, Qxz_t, Qyy_t, Qyz_t, Qzz_t
    real(kind=DP), intent(in)  :: trace

    ! ---------------------------------------------------------------------------

    Qxy_p = 2.0_DP/3.0_DP * Qxy_t
    Qxz_p = 2.0_DP/3.0_DP * Qxz_t
    Qyz_p = 2.0_DP/3.0_DP * Qyz_t

    Qxx_p =1.0_DP/3.0_DP * (trace + 2.0_DP * Qxx_t)
    Qyy_p =1.0_DP/3.0_DP * (trace + 2.0_DP * Qyy_t)
    Qzz_p =1.0_DP/3.0_DP * (trace + 2.0_DP * Qzz_t)

  end subroutine multipole_2_cart_traceless_to_primitive

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_2_cart_primitive_to_traceless(&
       Qxx_t, Qxy_t, Qxz_t, Qyy_t, Qyz_t, Qzz_t, &                ! output
       Qxx_p, Qxy_p, Qxz_p, Qyy_p, Qyz_p, Qzz_p, convention)      ! input
    !==========================================================================!
    ! Converts Cartesian primitive quadrupoles to Cartesian traceless in a     !
    ! number of conventions.                                                   !
    ! Racah-normalization has no effect here.                                  !
    ! References:                                                              !
    ! [1] eq. (2.59) for 'GrayGubbins' convention.                             !
    ! [2] eq. (6), [3] eq. (4.9) for 'Jackson' convention.                     !
    ! empirical evidence for 'Gaussian' convention.                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)   :: Qxx_t, Qxy_t, Qxz_t, Qyy_t, Qyz_t, Qzz_t
    real(kind=DP), intent(in)    :: Qxx_p, Qxy_p, Qxz_p, Qyy_p, Qyz_p, Qzz_p
    character(len=*), intent(in) :: convention

    ! jd: Local variables
    real(kind=DP), parameter     :: THIRD = 1.0_DP/3.0_DP

    ! ---------------------------------------------------------------------------

    select case (convention)
    case ('GrayGubbins')
       Qxx_t = 3.0_DP/2.0_DP * Qxx_p - 0.5_DP * (Qxx_p + Qyy_p + Qzz_p)
       Qyy_t = 3.0_DP/2.0_DP * Qyy_p - 0.5_DP * (Qxx_p + Qyy_p + Qzz_p)
       Qzz_t = 3.0_DP/2.0_DP * Qzz_p - 0.5_DP * (Qxx_p + Qyy_p + Qzz_p)
       Qxy_t = 3.0_DP/2.0_DP * Qxy_p
       Qxz_t = 3.0_DP/2.0_DP * Qxz_p
       Qyz_t = 3.0_DP/2.0_DP * Qyz_p
    case ('Jackson')
       Qxx_t = 3.0_DP * Qxx_p - (Qxx_p + Qyy_p + Qzz_p)
       Qyy_t = 3.0_DP * Qyy_p - (Qxx_p + Qyy_p + Qzz_p)
       Qzz_t = 3.0_DP * Qzz_p - (Qxx_p + Qyy_p + Qzz_p)
       Qxy_t = 3.0_DP * Qxy_p
       Qxz_t = 3.0_DP * Qxz_p
       Qyz_t = 3.0_DP * Qyz_p
    case ('Gaussian')
       Qxx_t = Qxx_p - THIRD * (Qxx_p + Qyy_p + Qzz_p)
       Qyy_t = Qyy_p - THIRD * (Qxx_p + Qyy_p + Qzz_p)
       Qzz_t = Qzz_p - THIRD * (Qxx_p + Qyy_p + Qzz_p)
       Qxy_t = Qxy_p
       Qxz_t = Qxz_p
       Qyz_t = Qyz_p
    case default
       call utils_abort('Unrecognized quadrupole convention in &
            &multipole_2_cart_primitive_to_traceless')
    end select
  end subroutine multipole_2_cart_primitive_to_traceless

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_1_translate_cartesian(&
       dipole_at_dest, dipole_at_src, monopole_at_src, src_pt, dest_pt, &
       component)
    !==========================================================================!
    ! Translates a Cartesian dipole from src_pt to dest_pt according to [4],   !
    ! eq. (18). The dipoles are in the GrayGubbins convention, [1], eq. (2.54) !
    ! (this is important, otherwise Racah factor would have to be included.)   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: dipole_at_dest
    real(kind=DP), intent(in)  :: dipole_at_src, monopole_at_src
    real(kind=DP), intent(in)  :: src_pt(3), dest_pt(3)
    integer, intent(in)        :: component ! 1=X, 2=Y, 3=Z

    ! ---------------------------------------------------------------------------

    call utils_assert(component >= 1 .and. component <= 3, &
         'Invalid dipole component in multipole_1_translate_cartesian')

    dipole_at_dest = dipole_at_src - &
         (src_pt(component) - dest_pt(component)) * monopole_at_src

  end subroutine multipole_1_translate_cartesian

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_2_translate_cartesian_primitive(&
       quadrupole_at_dest, quadrupole_at_src, dipole1_at_dest, dipole2_at_dest,&
       monopole_at_src, src_pt, dest_pt, component1, component2)
    !==========================================================================!
    ! Translates a Cartesian primitive quadrupole from src_pt to dest_pt,      !
    ! according to jd's notes, equation (**). The quadrupoles are according    !
    ! to [1] eq (2.54). This subroutine will *NOT* work with traceless         !
    ! quadrupoles, as the normalization factors would be off.                  !
    !--------------------------------------------------------------------------!
    ! NOTE: The specified dipoles (dipole1, dipole2) must be at dest (!), so   !
    !       will likely need to be translated beforehand.                      !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: quadrupole_at_dest
    real(kind=DP), intent(in)  :: quadrupole_at_src
    real(kind=DP), intent(in)  :: dipole1_at_dest, dipole2_at_dest
    real(kind=DP), intent(in)  :: monopole_at_src
    real(kind=DP), intent(in)  :: src_pt(3), dest_pt(3)
    integer, intent(in)        :: component1 ! 1=X, 2=Y, 3=Z
    integer, intent(in)        :: component2 ! 1=X, 2=Y, 3=Z

    ! ---------------------------------------------------------------------------

    call utils_assert(component1 >= 1 .and. component1 <= 3, &
         'Invalid dipole component1 in multipole_2_translate_cartesian_&
         &primitive')
    call utils_assert(component2 >= 1 .and. component2 <= 3, &
         'Invalid dipole component2 in multipole_2_translate_cartesian_&
         &primitive')

    quadrupole_at_dest = quadrupole_at_src - &
         (src_pt(component1) - dest_pt(component1)) * dipole2_at_dest - &
         (src_pt(component2) - dest_pt(component2)) * dipole1_at_dest - &
         (src_pt(component1) - dest_pt(component1)) * &
         (src_pt(component2) - dest_pt(component2)) * monopole_at_src

    ! jd: I honestly don't know where this minus comes from
    quadrupole_at_dest = -quadrupole_at_dest

  end subroutine multipole_2_translate_cartesian_primitive

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_1_cart_to_sph_real(Q1m1, Q10, Q11, & ! output
       mux, muy, muz, want_racah_normalized)                ! input
    !==========================================================================!
    ! Converts Cartesian dipoles to spherical real dipoles.                    !
    ! The output dipoles may or may not be Racah-normalized -- set the last    !
    ! argument appropriately.                                                  !
    ! Reference: [5], [6] Table 2 and cks's observation of the following       !
    ! trivial correspondence:                                                  !
    ! GDMA lingo         real-valued spherical multipoles ('m' denotes a minus)!
    ! Q10           =    Q10                                                   !
    ! Q11c          =    Q11                                                   !
    ! Q11s          =    Q1m1                                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    use constants, only: PI

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: Q1m1, Q10, Q11
    real(kind=DP), intent(in)  :: mux, muy, muz
    logical, intent(in)        :: want_racah_normalized

    ! jd: Local variables
    real(kind=DP), parameter   :: RACAH_FACTOR = SQRT(3.0_DP/(4.0_DP*PI))

    ! ---------------------------------------------------------------------------

    Q11 = mux
    Q1m1 = muy
    Q10 = muz

    if(.not. want_racah_normalized) then
       Q11 = Q11 * RACAH_FACTOR
       Q1m1 = Q1m1 * RACAH_FACTOR
       Q10 = Q10 * RACAH_FACTOR
    end if

  end subroutine multipole_1_cart_to_sph_real

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_2_cart_traceless_to_sph_real(&
       Q2m2, Q2m1, Q20, Q21, Q22, &                             ! output
       Qxx, Qxy, Qxz, Qyy, Qyz, Qzz, want_racah_normalized)     ! input
    !==========================================================================!
    ! Converts Cartesian traceless quadrupoles to spherical real quadrupoles.  !
    ! The output quadrupoles may or may not be Racah-normalized -- set the last!
    ! argument appropriately.                                                  !
    ! Reference: [5], [6] Table 2 and cks's observation of the following       !
    ! trivial correspondence:                                                  !
    ! GDMA lingo         real-valued spherical multipoles ('m' denotes a minus)!
    ! Q20           =    Q20                                                   !
    ! Q21c          =    Q21                                                   !
    ! Q21s          =    Q2m1                                                  !
    ! Q22c          =    Q22                                                   !
    ! Q22s          =    Q2m2                                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    use constants, only: PI

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: Q2m2, Q2m1, Q20, Q21, Q22
    real(kind=DP), intent(in)  :: Qxx, Qxy, Qxz, Qyy, Qyz, Qzz
    logical, intent(in)        :: want_racah_normalized

    ! jd: Local variables
    real(kind=DP), parameter   :: TWO_OVER_SQRT3 = 2.0_DP/SQRT(3.0_DP)
    real(kind=DP), parameter   :: RACAH_FACTOR = SQRT(5.0_DP/(4.0_DP*PI))

    ! ---------------------------------------------------------------------------

    Q20 = Qzz
    Q21 = TWO_OVER_SQRT3 * Qxz
    Q2m1 = TWO_OVER_SQRT3 * Qyz
    Q2m2 = TWO_OVER_SQRT3 * Qxy
    Q22 = (Qxx-Qyy) / sqrt(3.0_DP)

    if(.not. want_racah_normalized) then
       Q20 = Q20 * RACAH_FACTOR
       Q21 = Q21 * RACAH_FACTOR
       Q2m1 = Q2m1 * RACAH_FACTOR
       Q22 = Q22 * RACAH_FACTOR
       Q2m2 = Q2m2 * RACAH_FACTOR
    end if

  end subroutine multipole_2_cart_traceless_to_sph_real

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function multipole_1_magnitude_cartesian(mux, muy, muz)
    !==========================================================================!
    ! Returns the magnitude of a Cartesian dipole.                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: mux, muy, muz

    ! ---------------------------------------------------------------------------

    multipole_1_magnitude_cartesian = sqrt(mux*mux+muy*muy+muz*muz)

  end function multipole_1_magnitude_cartesian

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function multipole_2_magnitude_cartesian(&
       Qxx, Qxy, Qxz, Qyy, Qyz, Qzz)
    !==========================================================================!
    ! Returns the magnitude of a Cartesian quadrupole.                         !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2014.09.                                    !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: Qxx, Qxy, Qxz, Qyy, Qyz, Qzz

    ! ---------------------------------------------------------------------------

    multipole_2_magnitude_cartesian = sqrt(Qxx*Qxx + 2.0_DP*Qxy*Qxy + &
         2.0_DP*Qxz*Qxz + Qyy*Qyy + 2.0_DP*Qyz*Qyz + Qzz*Qzz)

  end function multipole_2_magnitude_cartesian

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function multipole_1_magnitude_sph_real(Q10, Q10s, Q10c)
    !==========================================================================!
    ! Returns the magnitude of a spherical real dipole.                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.03.                                    !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: Q10, Q10s, Q10c

    ! ---------------------------------------------------------------------------

    multipole_1_magnitude_sph_real = sqrt(Q10*Q10+Q10s*Q10s+Q10c*Q10c)

  end function multipole_1_magnitude_sph_real

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function multipole_2_magnitude_sph_real(&
       Q20, Q21c, Q21s, Q22c, Q22s)
    !==========================================================================!
    ! Returns the magnitude of a spherical real quadrupole.                    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2015.03.                                    !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: Q20, Q21c, Q21s, Q22c, Q22s

    ! ---------------------------------------------------------------------------

    multipole_2_magnitude_sph_real = sqrt(Q20*Q20 + Q21c*Q21c + Q21s*Q21s + &
         Q22c*Q22c + Q22s*Q22s)

  end function multipole_2_magnitude_sph_real

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_init_spherical_set(sph_multipoles, &         ! out
       min_l, max_l, n_sites, set_name, potential_damping, energy_damping) ! in
    !==========================================================================!
    ! Sets up the data structure describing a set of multipoles in spherical   !
    ! real representation.                                                     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2015
    !==========================================================================!

    use constants, only: garbage_real
    use geometry, only: POINT
    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! jd: Argument
    type(SPHERICAL_MULTIPOLE_SET), intent(out) :: sph_multipoles
    integer, intent(in)                        :: min_l
    integer, intent(in)                        :: max_l
    integer, intent(in)                        :: n_sites
    character(len=*), intent(in)               :: set_name
    integer, intent(in), optional              :: potential_damping
    integer, intent(in), optional              :: energy_damping

    ! jd: Local variables
    integer :: ierr

    type(POINT), parameter      :: garbage_point = &
         POINT(garbage_real, garbage_real, garbage_real)
    character(len=*), parameter :: myself = 'multipole_init_spherical_set'

    ! -------------------------------------------------------------------------

    call utils_assert(max_l>=0, myself//': max_l must be non-negative')
    call utils_assert(min_l<=max_l, myself//': min_l must be <= max_l')
    call utils_assert(n_sites>=0, myself//': n_sites must be non-negative')

    sph_multipoles%min_l = min_l
    sph_multipoles%max_l = max_l
    sph_multipoles%n_multipoles_per_site = (1+max_l-min_l)*(1+max_l+min_l) ! == sum_min_l^max_l (2l+1)
    sph_multipoles%n_sites = n_sites
    sph_multipoles%set_name = set_name
    sph_multipoles%dma_monopole_scaling_factor = 1D0
    sph_multipoles%dma_nelecs = garbage_real

    if(present(potential_damping)) then
       sph_multipoles%potential_damping = potential_damping
    else
       sph_multipoles%potential_damping = DAMPING_COULOMBIC_SMEARED
    end if
    if(present(energy_damping)) then
       sph_multipoles%energy_damping = energy_damping
    else
       sph_multipoles%energy_damping = DAMPING_FROM_POTENTIAL
    end if

    allocate(sph_multipoles%centres(n_sites),stat=ierr)
    call utils_alloc_check(myself,'sph_multipoles%centres',ierr)
    sph_multipoles%centres(:) = garbage_point

    allocate(sph_multipoles%min_l_mask(n_sites),stat=ierr)
    call utils_alloc_check(myself,'sph_multipoles%min_l_mask',ierr)
    sph_multipoles%min_l_mask(:) = 0

    allocate(sph_multipoles%max_l_mask(n_sites),stat=ierr)
    call utils_alloc_check(myself,'sph_multipoles%max_l_mask',ierr)
    sph_multipoles%max_l_mask(:) = max_l

    allocate(sph_multipoles%safe_radii(n_sites),stat=ierr)
    call utils_alloc_check(myself,'sph_multipoles%safe_radii',ierr)
    sph_multipoles%safe_radii(:) = 0D0

    allocate(sph_multipoles%polarisabilities(n_sites),stat=ierr)
    call utils_alloc_check(myself,'sph_multipoles%polarisabilities',ierr)
    sph_multipoles%polarisabilities(:) = 0D0

    allocate(sph_multipoles%userdata(3,n_sites),stat=ierr)
    call utils_alloc_check(myself,'sph_multipoles%userdata',ierr)
    sph_multipoles%userdata(:,:) = 0D0

    allocate(sph_multipoles%species(n_sites),stat=ierr)
    call utils_alloc_check(myself,'sph_multipoles%species',ierr)
    sph_multipoles%species(:) = '????'

    allocate(sph_multipoles%mpoles(sph_multipoles%n_multipoles_per_site * &
         n_sites),stat=ierr)
    call utils_alloc_check(myself,'sph_multipoles%mpoles',ierr)
    sph_multipoles%mpoles(:) = garbage_real

  end subroutine multipole_init_spherical_set

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_free_spherical_set(sph_multipoles) ! in/out
    !==========================================================================!
    ! Frees the data structure describing a set of multipoles in spherical     !
    ! real representation.                                                     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2015
    !==========================================================================!

    use constants, only: garbage_int, garbage_real
    use utils, only: utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(inout) :: sph_multipoles

    ! jd: Local variables
    integer :: ierr
    character(len=*), parameter :: myself = 'multipole_free_spherical_set'

    ! -------------------------------------------------------------------------

    deallocate(sph_multipoles%mpoles,stat=ierr)
    call utils_dealloc_check(myself,'sph_multipoles%mpoles',ierr)
    deallocate(sph_multipoles%species,stat=ierr)
    call utils_dealloc_check(myself,'sph_multipoles%species',ierr)
    deallocate(sph_multipoles%userdata,stat=ierr)
    call utils_dealloc_check(myself,'sph_multipoles%userdata',ierr)
    deallocate(sph_multipoles%polarisabilities,stat=ierr)
    call utils_dealloc_check(myself,'sph_multipoles%polarisabilities',ierr)
    deallocate(sph_multipoles%safe_radii,stat=ierr)
    call utils_dealloc_check(myself,'sph_multipoles%safe_radii',ierr)
    deallocate(sph_multipoles%min_l_mask,stat=ierr)
    call utils_dealloc_check(myself,'sph_multipoles%min_l_mask',ierr)
    deallocate(sph_multipoles%max_l_mask,stat=ierr)
    call utils_dealloc_check(myself,'sph_multipoles%max_l_mask',ierr)
    deallocate(sph_multipoles%centres,stat=ierr)
    call utils_dealloc_check(myself,'sph_multipoles%centres',ierr)

    sph_multipoles%max_l = garbage_int
    sph_multipoles%n_multipoles_per_site = garbage_int
    sph_multipoles%n_sites = garbage_int
    sph_multipoles%set_name = '!DELETED!'
    sph_multipoles%potential_damping = garbage_int
    sph_multipoles%energy_damping = garbage_int
    sph_multipoles%dma_monopole_scaling_factor = garbage_real
    sph_multipoles%dma_nelecs = garbage_real

  end subroutine multipole_free_spherical_set

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function multipole_pot_of_spherical_set(sph_multipoles, r_point, &
       alpha_dest)
    !==========================================================================!
    ! Computes the electrostatic potential at a single point, due to a set of  !
    ! spherical multipoles. The potential can be Coulombic (in which case      !
    ! there is no damping, and singularities are avoided through the use of an !
    ! exclusion radius -- with masking or smearing), or Thole-damped (in which !
    ! case damping is determined either by the polarisabilities of the source  !
    ! multipole and the target site (hence alpha_dest); or from an averaged    !
    ! polarisation (in the absence of alpha_dest)).                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   a_multipoles (in): A spherical multipole set.                          !
    !   r_point (in): The point where to compute the potential.                !
    !   alpha_dest (in, opt): If specified, it is used for the polarisability  !
    !                         of the target site (and the source used the      !
    !                         polarisability stored in sph_multipoles). If     !
    !                         omitted, pol_emb_pairwise_polarisability is used !
    !                         for A.                                           !
    !--------------------------------------------------------------------------!
    ! Re-written from scratch by Jacek Dziedzic in May 2015 to use             !
    ! SPHERICAL_MULTIPOLE_SET.                                                 !
    ! Support for damping and polarisabilities by Jacek Dziedzic, Oct-Dec 2015.!
    !==========================================================================!

    use constants, only: SAFE_DIV_EPS, ONE_SIXTH
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+), OPERATOR(*)
    use rundat, only: pub_pol_emb_mpole_exclusion_radius, &
         pub_pol_emb_pairwise_polarisability, pub_pol_emb_thole_a, &
         pub_pol_emb_smearing_a
    use spherical_wave, only: sw_real_sph_harm_unit
    use utils, only: utils_abort, utils_assert, utils_int_to_str

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: sph_multipoles
    type(POINT), intent(in)                   :: r_point
    real(kind=DP), intent(in), optional       :: alpha_dest

    ! jd: Local variables
    integer                     :: site
    integer                     :: i_l, i_m
    integer                     :: mpole_idx
    real(kind=DP)               :: result, result_for_this_site
    real(kind=DP)               :: factor
    real(kind=DP)               :: r_norm
    real(kind=DP)               :: r_inv
    type(POINT)                 :: dr_point
    type(POINT)                 :: dr_unit
    real(kind=DP)               :: big_A
    real(kind=DP)               :: damping
    ! jd: Racah coefficients sqrt(4pi/(2l+1)) for l=0..4
    real(kind=DP), parameter, dimension(0:4) :: racah = &
         (/ 3.5449077018110320546_DP, 2.0466534158929769770_DP, &
         1.5853309190424044053_DP, 1.3398491713813575449_DP, &
         1.1816359006036773515_DP /)
    real(kind=DP), parameter, dimension(0:4) :: minus_one_to_l = &
         (/ 1.0_DP, -1.0_DP, 1.0_DP, -1.0_DP, 1.0_DP /)
    character(len=*), parameter :: myself = 'multipole_pot_of_spherical_set'

    !------------------------------------------------------------------------

    result = 0D0
    select case(sph_multipoles%potential_damping)
    case(DAMPING_ZERO)
       goto 999 ! return zero
    case(DAMPING_COULOMBIC_MASKED)
    case(DAMPING_COULOMBIC_SMEARED)
    case(DAMPING_THOLE)
       continue ! Execute the usual below
    case default
       call utils_abort(myself//': Unrecognized potential_damping')
    end select

    if(sph_multipoles%max_l>4) call utils_abort(myself//': max_l too big')

    if(present(alpha_dest)) then
       call utils_assert(alpha_dest > 0D0, myself//': Illegal alpha_dest', &
            alpha_dest)
    end if

    mpole_idx = 1
    loop_sites:                                                               &
    do site = 1, sph_multipoles%n_sites

       ! jd: Skip dummy QM multipoles masquerading as MM multipoles
       if(sph_multipoles%max_l_mask(site) == -1) then
          ! If multipole site is outside the mask, make sure
          ! to fast-forward mpole_idx by the correct number.
          mpole_idx = mpole_idx + sph_multipoles%n_multipoles_per_site
          cycle
       end if

       dr_point = r_point - sph_multipoles%centres(site)
       r_norm = magnitude(dr_point) ! @MIC must be supported later on
       if (r_norm <= sph_multipoles%safe_radii(site)) then
          call utils_abort(trim(myself)//': Attempt to calculate the multipole&
               & potential inside a region where the expansion is not&
               & necessarily convergent.')
       end if

       ! jd: If Thole-damping, figure out big_A
       if(sph_multipoles%potential_damping == DAMPING_THOLE) then
          if(present(alpha_dest)) then
             big_A = &
                  (alpha_dest * sph_multipoles%polarisabilities(site))**ONE_SIXTH
          else
             big_A = pub_pol_emb_pairwise_polarisability
          end if

          if(big_A <= 0D0) call utils_assert(.false., &
               myself//': Illegal Thole big_A. Multipole site: '//&
               trim(utils_int_to_str(site))//', set: "'//&
               trim(sph_multipoles%set_name)//'"', big_A)
       end if

       ! jd: For Coulombic-masked potential crudely exclude immediate vicinity
       !     of multipole to avoid singularities, leaving the potential at 0.
       if(sph_multipoles%potential_damping == DAMPING_COULOMBIC_MASKED) then
          if(r_norm <= pub_pol_emb_mpole_exclusion_radius) then
             ! If multipole site is to be skipped, make sure
             ! to fast-forward mpole_idx by the correct number.
             mpole_idx = mpole_idx + sph_multipoles%n_multipoles_per_site
             cycle
          end if
       end if

       ! jd: Corner case for multipoles exactly on grid points.
       !     If ZERO, we're never in this loop at all.
       !     If THOLE, then 'damping' [*] contains the result for the damped
       !        potential taken in the limit (rather than just the damping
       !        factor, like for usual points).
       !     If COULOMBIC_MASKED or COULOMBIC_SMEARED, the only thing to do
       !        for now is *not* to divide by zero and proceed.
       if(r_norm < SAFE_DIV_EPS) then
          ! No-op
       else
          r_inv = 1.0_DP/r_norm
          dr_unit = r_inv * dr_point
       end if

       result_for_this_site = 0D0
       loop_l:                                                                &
       do i_l = sph_multipoles%min_l, sph_multipoles%max_l

          ! If this l is within the mask for the multipole, include it
          if(i_l <= sph_multipoles%max_l_mask(site) .and. &
               i_l >= sph_multipoles%min_l_mask(site)) then

             ! By default there is no damping
             damping = 1D0

             ! ... unless we're doing Thole damping
             if(sph_multipoles%potential_damping == DAMPING_THOLE) then
                damping = multipole_thole_damping(r_norm, pub_pol_emb_thole_a, &
                     big_A, i_l)
             end if

             ! ... or Thole-like smearing, and we're close to a site
             if(sph_multipoles%potential_damping == DAMPING_COULOMBIC_SMEARED) then
                if(r_norm <= pub_pol_emb_mpole_exclusion_radius) then
                   damping = multipole_thole_damping(r_norm, &
                        pub_pol_emb_thole_a, pub_pol_emb_smearing_a, i_l)
                end if
             end if

             if(r_norm < SAFE_DIV_EPS) then
                ! jd: Corner case for multipoles exactly on grid points.
                !     Avoid calculating 0.0^-1 below.
                factor = 0D0
             else
                ! @ This is the (-1)^l version. The (-1)^l will be moved to the
                ! definition of the multipole in the future to make things more
                ! consistent.
                factor = r_norm**(-(i_l+1)) * racah(i_l) * minus_one_to_l(i_l)
             end if

             loop_m:                                                           &
             do i_m = -i_l, i_l
                if(r_norm < SAFE_DIV_EPS) then
                   result_for_this_site = 1.0_DP ! cf. [*] above
                else
                   ! jd: Standard case
                   result_for_this_site = result_for_this_site + &
                        sph_multipoles%mpoles(mpole_idx) * factor * &
                        sw_real_sph_harm_unit(dr_unit%X, dr_unit%Y, dr_unit%Z, &
                        i_l, i_m)
                end if
                mpole_idx = mpole_idx + 1
             end do loop_m
          else
             ! If this l is outside the mask for the multipole site, make sure
             ! to fast-forward mpole_idx by the correct number.
             mpole_idx = mpole_idx + 2*i_l+1
          end if

       end do loop_l

       result = result + result_for_this_site * damping

    end do loop_sites

999 multipole_pot_of_spherical_set = result

  end function multipole_pot_of_spherical_set

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_mdq_cart_to_spherical_set(sph_multipoles, &         ! out
       mdqpoles_cart, n_sites, set_name, potential_damping, energy_damping)! in
    !==========================================================================!
    ! Converts an array (n_sites elements) of Cartesian primitive              !
    ! mono-di-quadru-poles: (CART_PRIM_MDQPOLE(n_sites)) to a spherical        !
    ! multipole set (SPHERICAL_MULTIPOLE_SET).                                 !
    !                                                                          !
    ! safe_radii(:) are set to zero. sph_multipoles are allocated and initia-  !
    ! lised here, through a call to multipole_init_spherical_set().            !
    ! The resultant multipoles are Racah-normalized.                           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   The caller is responsible for calling multipole_free_spherical_set()   !
    !   on the returned set.                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2015.                                   !
    !==========================================================================!

    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(out) :: sph_multipoles
    type(CART_PRIM_MDQPOLE), intent(in)        :: mdqpoles_cart(:)
    integer, intent(in)                        :: n_sites
    character(len=set_name_len), intent(in)    :: set_name
    integer, intent(in)                        :: potential_damping
    integer, intent(in)                        :: energy_damping

    ! jd: Local variables
    integer :: i
    integer :: offs
    real(kind=DP) :: Qxx_t, Qxy_t, Qxz_t, Qyy_t, Qyz_t, Qzz_t ! traceless qpole
    character(len=*), parameter :: myself='multipole_mdq_cart_to_spherical_set'

    ! -------------------------------------------------------------------------

    call utils_assert(n_sites>=0, myself//': n_sites must be non-negative')
    call utils_assert(ubound(mdqpoles_cart,1)==n_sites, myself//&
         ': Length of mdqpoles_cart(:) does not match n_sites: ', &
         ubound(mdqpoles_cart,1),n_sites)

    ! jd: Initialise a set of up-to-quadrupoles
    call multipole_init_spherical_set(sph_multipoles, 0, 2, n_sites, set_name, &
         potential_damping, energy_damping)

    ! jd: Populate
    do i = 1, n_sites
       offs = (i-1) * sph_multipoles%n_multipoles_per_site

       ! jd: Store centre
       sph_multipoles%centres(i) = mdqpoles_cart(i)%centre
       ! jd: Store polarisability
       sph_multipoles%polarisabilities(i) = mdqpoles_cart(i)%polarisability
       ! jd: Store userdata
       sph_multipoles%userdata(1:3,i) = mdqpoles_cart(i)%userdata(1:3)
       ! jd: Store species
       sph_multipoles%species(i) = mdqpoles_cart(i)%species
       ! jd: If original masked, mask this the sph multipole too
       if(mdqpoles_cart(i)%masked) then
          sph_multipoles%min_l_mask(i) = huge(1)
          sph_multipoles%max_l_mask(i) = -1
       else
          sph_multipoles%min_l_mask(i) = 0
          sph_multipoles%max_l_mask(i) = 2
       end if
       ! jd: Store charge
       sph_multipoles%mpoles(offs+1) = mdqpoles_cart(i)%charge
       ! jd: Store dipole
       call multipole_1_cart_to_sph_real(&
            sph_multipoles%mpoles(offs+2), &
            sph_multipoles%mpoles(offs+3), &
            sph_multipoles%mpoles(offs+4), &
            mdqpoles_cart(i)%dipole(1), &
            mdqpoles_cart(i)%dipole(2), &
            mdqpoles_cart(i)%dipole(3), want_racah_normalized = .true.)
       ! jd: Convert quadrupole to traceless
       call multipole_2_cart_primitive_to_traceless(&
            Qxx_t, Qxy_t, Qxz_t, Qyy_t, Qyz_t, Qzz_t, &
            mdqpoles_cart(i)%quadrupole(1), &
            mdqpoles_cart(i)%quadrupole(2), &
            mdqpoles_cart(i)%quadrupole(3), &
            mdqpoles_cart(i)%quadrupole(4), &
            mdqpoles_cart(i)%quadrupole(5), &
            mdqpoles_cart(i)%quadrupole(6), 'GrayGubbins')
       ! jd: Convert traceless to spherical and store
       call multipole_2_cart_traceless_to_sph_real(&
            sph_multipoles%mpoles(offs+5), &
            sph_multipoles%mpoles(offs+6), &
            sph_multipoles%mpoles(offs+7), &
            sph_multipoles%mpoles(offs+8), &
            sph_multipoles%mpoles(offs+9), &
            Qxx_t, Qxy_t, Qxz_t, Qyy_t, Qyz_t, Qzz_t, &
            want_racah_normalized = .true.)
       ! jd: If all quadrupoles are zero, apply a mask so that quadrupole
       !     potentials will be instantly evaluated as zero for efficiency.
       !     Similarly, if dipoles too are zero, and charges.
       !     Do not do this if it has been masked already.
       !     Do the analogous thing for min_l_mask
       if(sph_multipoles%max_l_mask(i) == 2) then ! only if not masked yet
          if(all(mdqpoles_cart(i)%quadrupole(:) == 0D0)) then
             sph_multipoles%max_l_mask(i) = 1
             if(all(mdqpoles_cart(i)%dipole(:) == 0D0)) then
                sph_multipoles%max_l_mask(i) = 0
                if(mdqpoles_cart(i)%charge == 0D0) then
                   sph_multipoles%max_l_mask(i) = -1
                end if
             end if
          end if

          if(mdqpoles_cart(i)%charge == 0D0) then
             sph_multipoles%min_l_mask(i) = 1
             if(all(mdqpoles_cart(i)%dipole(:) == 0D0)) then
                sph_multipoles%min_l_mask(i) = 2
                if(all(mdqpoles_cart(i)%quadrupole(:) == 0D0)) then
                   sph_multipoles%min_l_mask(i) = huge(1)
                end if
             end if
          end if
       end if

    end do

  end subroutine multipole_mdq_cart_to_spherical_set

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function multipole_thole_damping(r_norm, small_a, big_A, order)
    !==========================================================================!
    ! Calculates Thole damping factors (lambda).                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r_norm (in): |r| for which the damping factor is to be calculated.     !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    !   big_A (in): The value 'A' (distance scale) in the Thole formulas.      !
    !               Usually (alpha_i * alpha_j)^(1/6).                         !
    !   order (in): The order of the damping (of the Thole tensor):            !
    !               order=0 -> lambda_1                                        !
    !               order=1 -> lambda_3                                        !
    !               order=2 -> lambda_5                                        !
    !               order=3 -> lambda_7                                        !
    !               order=4 -> lambda_9                                        !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   As r_norm -> 0, the limiting behaviour has to be considered.           !
    !   Even if a particular lambda tends to 0 as r_norm->0, the potential     !
    !   itself may tend to a non-zero constant. For r_norm==0 or closer than   !
    !   SAFE_DIV_EPS this function returns *the potential* itself, rather than !
    !   the damping factor.                                                    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2015.                               !
    ! Lambda_9 added by Jacek Dziedzic in September 2016.                      !
    !==========================================================================!

    use constants, only: SAFE_DIV_EPS, ONE_THIRD, TWO_THIRDS
    use services, only: services_incomplete_gamma_from_lookup
    use utils, only: utils_abort, utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: r_norm
    real(kind=DP), intent(in) :: small_a
    real(kind=DP), intent(in) :: big_A
    integer, intent(in)       :: order

    ! jd: Local variables
    real(kind=DP) :: r_over_big_A_cubed, third_root_of_small_a
    character(len=*), parameter :: myself = 'multipole_thole_damping'

    ! -------------------------------------------------------------------------

    call utils_assert(big_A > SAFE_DIV_EPS, &
         myself//': Illegal big_A (check polarisations at sites or &
         &pol_emb_pairwise_polarisability)', big_A)
    call utils_assert(small_a > SAFE_DIV_EPS, &
         myself//': Illegal small_A (check pol_emb_thole_a)', small_A)

    r_over_big_A_cubed = (r_norm/big_A)**3.0_DP

    select case(order)
    ! ===================
    case(0) ! aka lambda_1
    ! ===================
       third_root_of_small_a = small_a**ONE_THIRD

       if(r_norm < SAFE_DIV_EPS) then
          ! jd: Corner case for potential exactly on the multipole site.
          !     Here we return not the damping factor (because that would lead
          !     to 0/0), but the limit of the damped potential at zero.
          !     This limit is (a^(1/3) Gamma[2/3]) / A
          multipole_thole_damping = &
               (third_root_of_small_a * &
               services_incomplete_gamma_from_lookup(&
               TWO_THIRDS, 0D0, gamma_two_thirds_lookup)) / big_A
       else
          ! jd: Usual case
          ! jd: cf. Burnham et al, J. Chem. Phys., Vol. 110, No. 9, 1 March 1999
          multipole_thole_damping = 1.0_DP - exp(-small_a*r_over_big_A_cubed) +&
               third_root_of_small_a * (r_norm/big_A) * &
               services_incomplete_gamma_from_lookup(TWO_THIRDS, &
               small_a * r_over_big_A_cubed, gamma_two_thirds_lookup)
       end if

    ! ===================
    case(1) ! aka lambda_3
    ! ===================
       if(r_norm < SAFE_DIV_EPS) then
          ! jd: Corner case for potential exactly on the multipole site.
          !     Here we return not the damping factor (because that would lead
          !     to 0/0), but the limit of the damped potential at zero.
          !     This limit is 0.0.
          multipole_thole_damping = 0D0
       else
          ! jd: cf. slides on "modified T matrix elements", source unknown
          ! jd: cf. Simmonett et al, JCP 143 07415 (2015), eq. (8), lambda_3
          multipole_thole_damping = 1.0_DP - exp(-small_a*r_over_big_A_cubed)
       end if

    ! ===================
    case(2) ! aka lambda_5
    ! ===================
       if(r_norm < SAFE_DIV_EPS) then
          ! jd: Corner case for potential exactly on the multipole site.
          !     Here we return not the damping factor (because that would lead
          !     to 0/0), but the limit of the damped potential at zero.
          !     This limit is 0.0.
          multipole_thole_damping = 0D0
       else
          ! jd: cf. slides on "modified T matrix elements", source unknown
          ! jd: cf. Simmonett et al, JCP 143 07415 (2015), eq. (8), lambda_5
          multipole_thole_damping = 1.0_DP - &
               (1.0_DP + small_a*r_over_big_A_cubed) * &
               exp(-small_a*r_over_big_A_cubed)
       end if
    ! ===================
    case(3) ! aka lambda_7
    ! ===================
       if(r_norm < SAFE_DIV_EPS) then
          ! jd: Corner case for potential exactly on the multipole site.
          !     Here we return not the damping factor (because that would lead
          !     to 0/0), but the limit of the damped potential at zero.
          !     This limit is 0.0.
          multipole_thole_damping = 0D0
       else
          ! jd: cf. Derivation 'thole_field_expressions_3.nb'.
          !     Consistent with esolv1.f in TINKER.
          multipole_thole_damping = 1.0_DP - &
               (1.0_DP + small_a*r_over_big_A_cubed + &
               0.6_DP * (small_a*r_over_big_A_cubed) ** 2.0_DP) * &
               exp(-small_a*r_over_big_A_cubed)
       end if
    ! ===================
    case(4) ! aka lambda_9
    ! ===================
       if(r_norm < SAFE_DIV_EPS) then
          ! jd: Corner case for potential exactly on the multipole site.
          !     Here we return not the damping factor (because that would lead
          !     to 0/0), but the limit of the damped potential at zero.
          !     This limit is 0.0.
          multipole_thole_damping = 0D0
       else
          ! jd: cf. Derivation 'thole_field_expressions_3.nb'
          !     Consistent with esolv1.f in TINKER.
          multipole_thole_damping = 1.0_DP - &
               (1.0_DP + small_a*r_over_big_A_cubed + &
               18.0_DP/35.0_DP * (small_a*r_over_big_A_cubed) ** 2.0_DP &
               - 9.0_DP/35.0_DP * (small_a*r_over_big_A_cubed) ** 3.0_DP) * &
               exp(-small_a*r_over_big_A_cubed)
       end if
    case default
       call utils_abort(myself//': Unrecognized Thole damping order:',order)
    end select

  end function multipole_thole_damping

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_T0(r1vec, r2vec, potential_damping, r0, &
       alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the multipole-multipole interaction tensor of order 0.           !
    ! (charge-charge interaction). Supports four electrostatics schemes:       !
    ! - DAMPING_ZERO (no interactions),                                        !
    ! - DAMPING_COULOMBIC_MASKED (Coulombic interactions that vanish sharply   !
    !   below r0),                                                             !
    ! - DAMPING_COULOMBIC_SMEARED (Coulombic interactions that taper smoothly  !
    !   to a finite value below r0),                                           !
    ! - DAMPING_THOLE (Thole-damped interactions).                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   potential_damping (in): Selects the damping scheme.                    !
    !   r0 (in, opt): Cut-off radius below which masking or smearing occurs.   !
    !                 Needed only for DAMPING_COULOMBIC_{MASKED,SMEARED}.      !
    !   alpha1 (in, opt): The polarisability of the first site.                !
    !                     Needed only for DAMPING_THOLE.                       !
    !   alpha2 (in, opt): The polarisability of the second site.               !
    !                     Needed only for DAMPING_THOLE.                       !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    !                 Needed only for DAMPING_COULOMBIC_SMEARED, DAMPING_THOLE.!
    ! Return value:                                                            !
    !   T0(r1vec-r2vec), a scalar.                                             !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Some, for Thole-damped interactions. See comment in multipole_thole_T0.!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)
    use utils, only: utils_abort

    implicit none

    ! jd: Return value
    real(kind=DP) :: multipole_T0

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    integer, intent(in)                 :: potential_damping
    real(kind=DP), intent(in), optional :: r0
    real(kind=DP), intent(in), optional :: alpha1, alpha2
    real(kind=DP), intent(in), optional :: small_a

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: r_norm
    character(len=*), parameter :: myself = 'multipole_T0'

    ! -------------------------------------------------------------------------

    ! jd: Short-circuit for zero, to avoid calculating r_norm unnecessarily
    if(potential_damping == DAMPING_ZERO) then
       multipole_T0 = 0D0
       return
    end if

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    select case(potential_damping)
    case(DAMPING_COULOMBIC_MASKED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_MASKED interactions')
       multipole_T0 = multipole_coul_masked_T0(r1vec, r2vec, r0)

    case(DAMPING_COULOMBIC_SMEARED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_SMEARED interactions')
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for COULOMBIC_SMEARED interactions')
       multipole_T0 = multipole_coul_smeared_T0(r1vec, r2vec, r0, small_a)

    case(DAMPING_THOLE)
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for THOLE interactions')
       if(.not. present(alpha1)) call utils_abort(myself//&
            ': alpha1 must be present for THOLE interactions')
       if(.not. present(alpha2)) call utils_abort(myself//&
            ': alpha2 must be present for THOLE interactions')
       multipole_T0 = multipole_thole_T0(r1vec, r2vec, alpha1, alpha2, small_a)

    case default
       call utils_abort(myself//': Unrecognized potential_damping')
    end select

  end function multipole_T0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_T1(r1vec, r2vec, potential_damping, r0, &
       alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the multipole-multipole interaction tensor of order 1.           !
    ! (charge-dipole interaction). Supports four electrostatics schemes:       !
    ! - DAMPING_ZERO (no interactions),                                        !
    ! - DAMPING_COULOMBIC_MASKED (Coulombic interactions that vanish sharply   !
    !   below r0),                                                             !
    ! - DAMPING_COULOMBIC_SMEARED (Coulombic interactions that taper smoothly  !
    !   to a finite value below r0),                                           !
    ! - DAMPING_THOLE (Thole-damped interactions).                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   potential_damping (in): Selects the damping scheme.                    !
    !   r0 (in, opt): Cut-off radius below which masking or smearing occurs.   !
    !                 Needed only for DAMPING_COULOMBIC_{MASKED,SMEARED}.      !
    !   alpha1 (in, opt): The polarisability of the first site.                !
    !                     Needed only for DAMPING_THOLE.                       !
    !   alpha2 (in, opt): The polarisability of the second site.               !
    !                     Needed only for DAMPING_THOLE.                       !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    !                 Needed only for DAMPING_COULOMBIC_SMEARED, DAMPING_THOLE.!
    ! Return value:                                                            !
    !   T1(r1vec-r2vec), a 3-vector.                                           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Some, for Thole-damped interactions. See comment in multipole_thole_T1.!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3) :: multipole_T1

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    integer, intent(in)                 :: potential_damping
    real(kind=DP), intent(in), optional :: r0
    real(kind=DP), intent(in), optional :: alpha1, alpha2
    real(kind=DP), intent(in), optional :: small_a

    ! jd: Local variables
    character(len=*), parameter :: myself = 'multipole_T1'

    ! -------------------------------------------------------------------------

    ! jd: Short-circuit for zero, to avoid calculating r_norm unnecessarily
    if(potential_damping == DAMPING_ZERO) then
       multipole_T1 = rank1zero
       return
    end if

    select case(potential_damping)
    case(DAMPING_COULOMBIC_MASKED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_MASKED interactions')
       multipole_T1 = multipole_coul_masked_T1(r1vec, r2vec, r0)

    case(DAMPING_COULOMBIC_SMEARED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_SMEARED interactions')
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for COULOMBIC_SMEARED interactions')
       multipole_T1 = multipole_coul_smeared_T1(r1vec, r2vec, r0, small_a)

    case(DAMPING_THOLE)
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for THOLE interactions')
       if(.not. present(alpha1)) call utils_abort(myself//&
            ': alpha1 must be present for THOLE interactions')
       if(.not. present(alpha2)) call utils_abort(myself//&
            ': alpha2 must be present for THOLE interactions')
       multipole_T1 = multipole_thole_T1(r1vec, r2vec, alpha1, alpha2, small_a)

    case default
       call utils_abort(myself//': Unrecognized potential_damping')
    end select

  end function multipole_T1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_T2(r1vec, r2vec, potential_damping, r0, &
       alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the multipole-multipole interaction tensor of order 2.           !
    ! (charge-quadrupole or dipole-dipole interaction).                        !
    ! Supports four electrostatics schemes:                                    !
    ! - DAMPING_ZERO (no interactions),                                        !
    ! - DAMPING_COULOMBIC_MASKED (Coulombic interactions that vanish sharply   !
    !   below r0),                                                             !
    ! - DAMPING_COULOMBIC_SMEARED (Coulombic interactions that taper smoothly  !
    !   to a finite value below r0),                                           !
    ! - DAMPING_THOLE (Thole-damped interactions).                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   potential_damping (in): Selects the damping scheme.                    !
    !   r0 (in, opt): Cut-off radius below which masking or smearing occurs.   !
    !                 Needed only for DAMPING_COULOMBIC_{MASKED,SMEARED}.      !
    !   alpha1 (in, opt): The polarisability of the first site.                !
    !                     Needed only for DAMPING_THOLE.                       !
    !   alpha2 (in, opt): The polarisability of the second site.               !
    !                     Needed only for DAMPING_THOLE.                       !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    !                 Needed only for DAMPING_COULOMBIC_SMEARED, DAMPING_THOLE.!
    ! Return value:                                                            !
    !   T2(r1vec-r2vec), a 3-matrix.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3) :: multipole_T2

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    integer, intent(in)                 :: potential_damping
    real(kind=DP), intent(in), optional :: r0
    real(kind=DP), intent(in), optional :: alpha1, alpha2
    real(kind=DP), intent(in), optional :: small_a

    ! jd: Local variables
    character(len=*), parameter :: myself = 'multipole_T2'

    ! -------------------------------------------------------------------------

    ! jd: Short-circuit for zero, to avoid calculating r_norm unnecessarily
    if(potential_damping == DAMPING_ZERO) then
       multipole_T2 = rank2zero
       return
    end if

    select case(potential_damping)
    case(DAMPING_COULOMBIC_MASKED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_MASKED interactions')
       multipole_T2 = multipole_coul_masked_T2(r1vec, r2vec, r0)

    case(DAMPING_COULOMBIC_SMEARED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_SMEARED interactions')
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for COULOMBIC_SMEARED interactions')
       multipole_T2 = multipole_coul_smeared_T2(r1vec, r2vec, r0, small_a)

    case(DAMPING_THOLE)
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for THOLE interactions')
       if(.not. present(alpha1)) call utils_abort(myself//&
            ': alpha1 must be present for THOLE interactions')
       if(.not. present(alpha2)) call utils_abort(myself//&
            ': alpha2 must be present for THOLE interactions')
       multipole_T2 = multipole_thole_T2(r1vec, r2vec, alpha1, alpha2, small_a)

    case default
       call utils_abort(myself//': Unrecognized potential_damping')
    end select

  end function multipole_T2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_T3(r1vec, r2vec, potential_damping, r0, &
       alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the multipole-multipole interaction tensor of order 3.           !
    ! (charge-octupole, dipole-quadrupole interaction).                        !
    ! Supports four electrostatics schemes:                                    !
    ! - DAMPING_ZERO (no interactions),                                        !
    ! - DAMPING_COULOMBIC_MASKED (Coulombic interactions that vanish sharply   !
    !   below r0),                                                             !
    ! - DAMPING_COULOMBIC_SMEARED (Coulombic interactions that taper smoothly  !
    !   to a finite value below r0),                                           !
    ! - DAMPING_THOLE (Thole-damped interactions).                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   potential_damping (in): Selects the damping scheme.                    !
    !   r0 (in, opt): Cut-off radius below which masking or smearing occurs.   !
    !                 Needed only for DAMPING_COULOMBIC_{MASKED,SMEARED}.      !
    !   alpha1 (in, opt): The polarisability of the first site.                !
    !                     Needed only for DAMPING_THOLE.                       !
    !   alpha2 (in, opt): The polarisability of the second site.               !
    !                     Needed only for DAMPING_THOLE.                       !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    !                 Needed only for DAMPING_COULOMBIC_SMEARED, DAMPING_THOLE.!
    ! Return value:                                                            !
    !   T3(r1vec-r2vec), a 3-(rank-3 tensor).                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3,3) :: multipole_T3

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    integer, intent(in)                 :: potential_damping
    real(kind=DP), intent(in), optional :: r0
    real(kind=DP), intent(in), optional :: alpha1, alpha2
    real(kind=DP), intent(in), optional :: small_a

    ! jd: Local variables
    character(len=*), parameter :: myself = 'multipole_T3'

    ! -------------------------------------------------------------------------

    ! jd: Short-circuit for zero, to avoid calculating r_norm unnecessarily
    if(potential_damping == DAMPING_ZERO) then
       multipole_T3 = rank3zero
       return
    end if

    select case(potential_damping)
    case(DAMPING_COULOMBIC_MASKED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_MASKED interactions')
       multipole_T3 = multipole_coul_masked_T3(r1vec, r2vec, r0)

    case(DAMPING_COULOMBIC_SMEARED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_SMEARED interactions')
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for COULOMBIC_SMEARED interactions')
       multipole_T3 = multipole_coul_smeared_T3(r1vec, r2vec, r0, small_a)

    case(DAMPING_THOLE)
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for THOLE interactions')
       if(.not. present(alpha1)) call utils_abort(myself//&
            ': alpha1 must be present for THOLE interactions')
       if(.not. present(alpha2)) call utils_abort(myself//&
            ': alpha2 must be present for THOLE interactions')
       multipole_T3 = multipole_thole_T3(r1vec, r2vec, alpha1, alpha2, small_a)

    case default
       call utils_abort(myself//': Unrecognized potential_damping')
    end select

  end function multipole_T3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_T4(r1vec, r2vec, potential_damping, r0, &
       alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the multipole-multipole interaction tensor of order 4.           !
    ! (charge-hexadecapole, dipole-octupole, quadrupole-quadrupole interaxn).  !
    ! Supports four electrostatics schemes:                                    !
    ! - DAMPING_ZERO (no interactions),                                        !
    ! - DAMPING_COULOMBIC_MASKED (Coulombic interactions that vanish sharply   !
    !   below r0),                                                             !
    ! - DAMPING_COULOMBIC_SMEARED (Coulombic interactions that taper smoothly  !
    !   to a finite value below r0),                                           !
    ! - DAMPING_THOLE (Thole-damped interactions).                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   potential_damping (in): Selects the damping scheme.                    !
    !   r0 (in, opt): Cut-off radius below which masking or smearing occurs.   !
    !                 Needed only for DAMPING_COULOMBIC_{MASKED,SMEARED}.      !
    !   alpha1 (in, opt): The polarisability of the first site.                !
    !                     Needed only for DAMPING_THOLE.                       !
    !   alpha2 (in, opt): The polarisability of the second site.               !
    !                     Needed only for DAMPING_THOLE.                       !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    !                 Needed only for DAMPING_COULOMBIC_SMEARED, DAMPING_THOLE.!
    ! Return value:                                                            !
    !   T4(r1vec-r2vec), a 3-(rank-4 tensor).                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3,3,3) :: multipole_T4

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    integer, intent(in)                 :: potential_damping
    real(kind=DP), intent(in), optional :: r0
    real(kind=DP), intent(in), optional :: alpha1, alpha2
    real(kind=DP), intent(in), optional :: small_a

    ! jd: Local variables
    character(len=*), parameter :: myself = 'multipole_T4'

    ! -------------------------------------------------------------------------

    ! jd: Short-circuit for zero, to avoid calculating r_norm unnecessarily
    if(potential_damping == DAMPING_ZERO) then
       multipole_T4 = rank4zero
       return
    end if

    select case(potential_damping)
    case(DAMPING_COULOMBIC_MASKED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_MASKED interactions')
       multipole_T4 = multipole_coul_masked_T4(r1vec, r2vec, r0)

    case(DAMPING_COULOMBIC_SMEARED)
       if(.not. present(r0)) call utils_abort(myself//&
            ': r0 must be present for COULOMBIC_SMEARED interactions')
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for COULOMBIC_SMEARED interactions')
       multipole_T4 = multipole_coul_smeared_T4(r1vec, r2vec, r0, small_a)

    case(DAMPING_THOLE)
       if(.not. present(small_a)) call utils_abort(myself//&
            ': small_a must be present for THOLE interactions')
       if(.not. present(alpha1)) call utils_abort(myself//&
            ': alpha1 must be present for THOLE interactions')
       if(.not. present(alpha2)) call utils_abort(myself//&
            ': alpha2 must be present for THOLE interactions')
       multipole_T4 = multipole_thole_T4(r1vec, r2vec, alpha1, alpha2, small_a)

    case default
       call utils_abort(myself//': Unrecognized potential_damping')
    end select

  end function multipole_T4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_set_info(mpole_set, fileunit, on_root_only)
    !==========================================================================!
    ! Writes human-readable information on a multipole set to a file unit.     !
    ! If called on multiple procs, attempts to order outputs.                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   mpole_set (in): The set for which information is desired.              !
    !   fileunit (in): The file unit to which information is to be written.    !
    !   on_root_only (in): Pass .true. to only print on root. Pass .false. to  !
    !                      print on all procs.                                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2017.                                  !
    !==========================================================================!

    use comms, only: pub_total_num_procs, pub_my_proc_id, pub_root_proc_id, &
         comms_barrier

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: mpole_set
    integer, intent(in) :: fileunit
    logical, intent(in) :: on_root_only

    ! jd: Local variables
    integer :: rank, min_rank, max_rank

    ! -------------------------------------------------------------------------

    if(on_root_only) then
       min_rank = pub_root_proc_id
       max_rank = pub_root_proc_id
    else
       min_rank = 0
       max_rank = pub_total_num_procs
    end if

    do rank = min_rank, max_rank
       if(.not. on_root_only) call comms_barrier
       if(pub_my_proc_id == rank) then
          write(fileunit, '(a,i0,a)') '================ MPI rank #', rank, &
               ' ================'
           write(fileunit, *) "Multipole set '"//trim(mpole_set%set_name)//"'"
          write(fileunit, *) 'min_l: ', mpole_set%min_l
          write(fileunit, *) 'max_l: ', mpole_set%max_l
          write(fileunit, *) 'n_multipoles_per_site: ', &
               mpole_set%n_multipoles_per_site
          write(fileunit, *) 'n_sites: ', mpole_set%n_sites
          write(fileunit, *) 'potential_damping: ', mpole_set%potential_damping
          write(fileunit, *) 'energy_damping: ', mpole_set%energy_damping
          write(fileunit, *) 'dma_monopole_scaling_factor: ', &
               mpole_set%dma_monopole_scaling_factor
          write(fileunit, *) 'dma_nelecs: ', mpole_set%dma_nelecs
          write(fileunit, *) 'mpoles: '
          write(fileunit, *) mpole_set%mpoles(1:mpole_set%n_sites * &
               mpole_set%n_multipoles_per_site)
          write(fileunit, *) 'polarisabilities: '
          write(fileunit, *) mpole_set%polarisabilities(1:mpole_set%n_sites)
          write(fileunit, *) 'userdata: '
          write(fileunit, *) mpole_set%userdata(1:3,1:mpole_set%n_sites)
          write(fileunit, *) 'species: '
          write(fileunit, *) mpole_set%species(1:mpole_set%n_sites)
          write(fileunit, *) 'min_l_mask: '
          write(fileunit, *) mpole_set%min_l_mask(1:mpole_set%n_sites)
          write(fileunit, *) 'max_l_mask: '
          write(fileunit, *) mpole_set%max_l_mask(1:mpole_set%n_sites)
          write(fileunit, *) 'safe_radii: '
          write(fileunit, *) mpole_set%safe_radii(1:mpole_set%n_sites)
          write(fileunit, *) 'centres: '
          write(fileunit, *) mpole_set%centres(1:mpole_set%n_sites)
       end if
    end do

  end subroutine multipole_set_info

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_set_merge(set_out, &
       name_for_set_out, set_in1, set_in2, alpha1, alpha2)
    !==========================================================================!
    ! Merges two multipole sets into one set. Resulant set spans an l range    !
    ! that is an OR of the two input set ranges. Multipoles absent from either !
    ! of the sets are taken as zero.                                           !
    ! Optional arguments alpha1 and alpha2 can be used to weight multipoles    !
    ! from either or both sets.                                                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! @docme. Caller responsible for freeing set_out.
    ! @nsites has to be identical, dampings too
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2017.                                  !
    !==========================================================================!

    use geometry, only: operator(-), geometry_magnitude
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(out) :: set_out
    character(len=*), intent(in)               :: name_for_set_out
    type(SPHERICAL_MULTIPOLE_SET), intent(in)  :: set_in1
    type(SPHERICAL_MULTIPOLE_SET), intent(in)  :: set_in2
    real(kind=DP), intent(in), optional        :: alpha1
    real(kind=DP), intent(in), optional        :: alpha2

    ! jd: Local variables
    integer :: min_l1, max_l1
    integer :: min_l2, max_l2
    integer :: min_l, max_l
    integer :: n_sites1, n_sites2
    integer :: potential_damping1, potential_damping2
    integer :: energy_damping1, energy_damping2
    integer :: site
    integer :: offs_in1, offs_in2, offs_out
    integer :: l_out, m_out
    integer :: i_mpole_out, i_mpole_in1, i_mpole_in2
    real(kind=DP) :: loc_alpha1, loc_alpha2
    type(POINT) :: centre1, centre2
    real(kind=DP) :: pol1, pol2
    real(kind=DP) :: mpole_in1, mpole_in2
    character(len=*), parameter :: myself = 'multipole_set_merge'

    ! -------------------------------------------------------------------------

    ! jd: Span of l for the out set (encompasses both in sets)
    min_l1 = set_in1%min_l
    max_l1 = set_in1%max_l
    min_l2 = set_in2%min_l
    max_l2 = set_in2%max_l
    min_l = min(min_l1,min_l2)
    max_l = max(max_l1,max_l2)

    ! jd: Sanity checks on set compatibility
    n_sites1 = set_in1%n_sites
    n_sites2 = set_in2%n_sites
    call utils_assert(n_sites1 == n_sites2, myself//&
         ': Sets have different number of sites: ', n_sites1, n_sites2)
    potential_damping1 = set_in1%potential_damping
    potential_damping2 = set_in2%potential_damping
    call utils_assert(potential_damping1 == potential_damping2, &
         myself//': Sets have different potential_damping: ', &
         potential_damping1, potential_damping2)
    energy_damping1 = set_in1%energy_damping
    energy_damping2 = set_in2%energy_damping
    call utils_assert(energy_damping1 == energy_damping2, &
         myself//': Sets have different energy_damping: ', &
         energy_damping1, energy_damping2)
    do site = 1, n_sites1
       centre1 = set_in1%centres(site)
       centre2 = set_in2%centres(site)
       call utils_assert(abs(geometry_magnitude(centre1-centre2)) < 1D-20, &
            myself//': Sets have at least one different centre (pos.) at index ', &
            site)
       pol1 = set_in1%polarisabilities(site)
       pol2 = set_in2%polarisabilities(site)
       call utils_assert(abs(pol1 - pol2) < 1D-20, &
            myself//': Sets have at least one different polarisability at index ', &
            real(site,kind=DP), pol1, pol2)
    end do

    ! jd: Local copies of optional multiplicative factors
    if(present(alpha1)) then
       loc_alpha1 = alpha1
    else
       loc_alpha1 = 1D0
    end if
    if(present(alpha2)) then
       loc_alpha2 = alpha2
    else
       loc_alpha2 = 1D0
    end if

    ! jd: Create the out set
    call multipole_init_spherical_set(set_out, min_l, max_l, n_sites1, &
         name_for_set_out, potential_damping1, energy_damping1)

    ! jd: Copy over centres
    set_out%centres(:) = set_in1%centres(:)

    ! jd: Copy over polarisabilities
    set_out%polarisabilities(:) = set_in1%polarisabilities(:)

    ! jd: Drop min_l and max_l masks -- they will be used below when merging
    !     multipoles, but are not copied over to out set.
    set_out%min_l_mask(:) = min_l
    set_out%max_l_mask(:) = max_l

    ! jd: Merge safe_radii
    set_out%safe_radii(:) = min(set_in1%safe_radii(:), set_in2%safe_radii(:))

    ! jd: Merge multipoles themselves
    do site = 1, n_sites1
       l_out = min_l
       m_out = -min_l
       offs_out = (site-1) * set_out%n_multipoles_per_site
       offs_in1 = (site-1) * set_in1%n_multipoles_per_site
       offs_in2 = (site-1) * set_in2%n_multipoles_per_site

       do i_mpole_out = 1, set_out%n_multipoles_per_site

          ! jd: Extract the multipole from both input sets, or fake as zero
          !     if absent.
          if(l_out >= set_in1%min_l .and. l_out <= set_in1%max_l .and. &
               l_out >= set_in1%min_l_mask(site) .and. &
               l_out <= set_in1%max_l_mask(site)) then
             i_mpole_in1 = multipole_lm_to_offset(l_out, m_out, set_in1)
             mpole_in1 = set_in1%mpoles(offs_in1 + i_mpole_in1)
          else
             mpole_in1 = 0D0
          end if
          if(l_out >= set_in2%min_l .and. l_out <= set_in2%max_l .and. &
               l_out >= set_in2%min_l_mask(site) .and. &
               l_out <= set_in2%max_l_mask(site)) then
             i_mpole_in2 = multipole_lm_to_offset(l_out, m_out, set_in2)
             mpole_in2 = set_in2%mpoles(offs_in2 + i_mpole_in2)
          else
             mpole_in2 = 0D0
          end if

          ! jd: Combine the two multipole components
          set_out%mpoles(offs_out + i_mpole_out) = &
               loc_alpha1 * mpole_in1 + loc_alpha2 * mpole_in2

          ! jd: Next multipole in output
          m_out = m_out + 1
          if(m_out > l_out) then
             l_out = l_out + 1
             m_out = -l_out
          end if

       end do
    end do

  end subroutine multipole_set_merge

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_set_gather(all_multipoles, my_multipoles)
    !==========================================================================!
    ! Gathers an MPI-distributed multipole set 'my_multipoles' onto root proc. !
    ! The result is returned in 'all_multipoles', on root.                     !
    ! Caller is responsible for freeing 'all_multipoles' later on.             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   my_multipoles (in): A multipole set for atoms local to this proc.      !
    !   all_multipoles (out): A multipole set gathered over all MPI ranks.     !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Since all_multipoles is intent(out), we must return something on every !
    !   rank. For ranks that are not root we return a dummy set with 0 sites.  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2017.                                  !
    !==========================================================================!

    use comms, only: pub_on_root, comms_send, comms_reduce, comms_allgather
    use utils, only: utils_alloc_check

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(out) :: all_multipoles
    type(SPHERICAL_MULTIPOLE_SET), intent(in)  :: my_multipoles

    ! jd: Local variables
    integer :: n_sites
    character(len=*), parameter :: myself = 'multipole_set_gather'

    ! -------------------------------------------------------------------------

    ! jd: Create a set encompassing all multipoles (on root), or a dummy,
    !     empty set (on non-root ranks)
    n_sites = my_multipoles%n_sites
    call comms_reduce('SUM', n_sites)
    if(.not. pub_on_root) then
       n_sites = 0
    end if
    call multipole_init_spherical_set(all_multipoles, &                  ! out
         my_multipoles%min_l, my_multipoles%max_l, &                     ! in
         n_sites, my_multipoles%set_name, &                              ! in
         my_multipoles%potential_damping, my_multipoles%energy_damping)  ! in

    ! jd: Gather min_l_mask
    call comms_allgather(all_multipoles%min_l_mask(:), &
         my_multipoles%min_l_mask(:), my_multipoles%n_sites, &
         (/-1/), (/-1/), gather_not_allgather = .true.)

    ! jd: Gather max_l_mask
    call comms_allgather(all_multipoles%max_l_mask(:), &
         my_multipoles%max_l_mask(:), my_multipoles%n_sites, &
         (/-1/), (/-1/), gather_not_allgather = .true.)

    ! jd: Gather safe_radii
    call comms_allgather(all_multipoles%safe_radii(:), &
         my_multipoles%safe_radii(:), my_multipoles%n_sites, &
         (/-1/), (/-1/), gather_not_allgather = .true.)

    ! jd: Gather polarisabilities
    call comms_allgather(all_multipoles%polarisabilities(:), &
         my_multipoles%polarisabilities(:), my_multipoles%n_sites, &
         (/-1/), (/-1/), gather_not_allgather = .true.)

    ! jd: Gather centres. Ugh, there is no comms_allgather for type(POINT)
    call comms_allgather(all_multipoles%centres(:)%x, &
         my_multipoles%centres(:)%x, my_multipoles%n_sites, &
         (/-1/), (/-1/), gather_not_allgather = .true.)
    call comms_allgather(all_multipoles%centres(:)%y, &
         my_multipoles%centres(:)%y, my_multipoles%n_sites, &
         (/-1/), (/-1/), gather_not_allgather = .true.)
    call comms_allgather(all_multipoles%centres(:)%z, &
         my_multipoles%centres(:)%z, my_multipoles%n_sites, &
         (/-1/), (/-1/), gather_not_allgather = .true.)

    ! jd: Gather mpoles
    call comms_allgather(all_multipoles%mpoles(:), &
         my_multipoles%mpoles(:), &
         my_multipoles%n_sites * my_multipoles%n_multipoles_per_site, &
         (/-1/), (/-1/), gather_not_allgather = .true.)

    ! jd: Copy dma_nelecs and dma_monopole_scaling_factor from root
    !     We assume all MPI ranks to have the same value anyway
    all_multipoles%dma_nelecs = my_multipoles%dma_nelecs
    all_multipoles%dma_monopole_scaling_factor = &
         my_multipoles%dma_monopole_scaling_factor

  end subroutine multipole_set_gather

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_set_scale_monopoles(multipoles, scaling_factor)
    !==========================================================================!
    ! Scales all monopoles, if any, in a multipole set by a scaling factor.    !
    ! Also stores 'scaling_factor' as 'dma_monopole_scaling_factor'.           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   multipoles (in/out): A multipole set whose monopoles are to be scaled. !
    !   scaling_factor (in): A multiplicative factor to apply to monopoles.    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2017.                                  !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(inout) :: multipoles
    real(kind=DP), intent(in)                    :: scaling_factor

    ! jd: Local variables
    integer :: i_site
    integer :: monopole_offset
    integer :: site_offset

    ! -------------------------------------------------------------------------

    if(multipoles%min_l <= 0 .and. multipoles%max_l >= 0) then
       do i_site = 1, multipoles%n_sites
          if(multipoles%min_l_mask(i_site) <= 0 .and. &
               multipoles%max_l_mask(i_site) >= 0) then
             monopole_offset = multipole_lm_to_offset(0,0,multipoles)
             site_offset = (i_site-1) * multipoles%n_multipoles_per_site
             multipoles%mpoles(site_offset + monopole_offset) = &
                  multipoles%mpoles(site_offset + monopole_offset) * &
                  scaling_factor
          end if
       end do
    end if

    multipoles%dma_monopole_scaling_factor = scaling_factor

  end subroutine multipole_set_scale_monopoles

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine multipole_set_scale_quadrupoles(multipoles, scaling_factor)
    !==========================================================================!
    ! Scales all quadrupoles, if any, in a multipole set by a scaling factor.  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   multipoles (in/out): A multipole set whose q'poles are to be scaled.   !
    !   scaling_factor (in): A multiplicative factor to apply to quadrupoles.  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2017.                                  !
    !==========================================================================!

    implicit none

    ! jd: Arguments
    type(SPHERICAL_MULTIPOLE_SET), intent(inout) :: multipoles
    real(kind=DP), intent(in)                    :: scaling_factor

    ! jd: Local variables
    integer :: i_site
    integer :: i_m
    integer :: qpole_offset
    integer :: site_offset

    ! -------------------------------------------------------------------------

    if(multipoles%min_l <= 0 .and. multipoles%max_l >= 0) then
       do i_site = 1, multipoles%n_sites
          if(multipoles%min_l_mask(i_site) <= 0 .and. &
               multipoles%max_l_mask(i_site) >= 0) then
             site_offset = (i_site-1) * multipoles%n_multipoles_per_site
             do i_m = -2, 2
                qpole_offset = multipole_lm_to_offset(2,i_m,multipoles)
                multipoles%mpoles(site_offset + qpole_offset) = &
                     multipoles%mpoles(site_offset + qpole_offset) * &
                     scaling_factor
             end do
          end if
       end do
    end if

    multipoles%dma_monopole_scaling_factor = scaling_factor

  end subroutine multipole_set_scale_quadrupoles

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ---------------------------------------------------------------------------
  ! ---- P R I V A T E   S U B R O U T I N E S   A N D   F U N C T I O N S ----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function multipole_lm_to_offset(l, m, mpole_set)
    !==========================================================================!
    ! Returns an offset (multipole index) for a given l, m.                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   l (in): Multipole angular momentum number.                             !
    !   m (in): Multipole m number.                                            !
    !   mpole_set (in): Multipole set on which we're operating.                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in July 2017.                                  !
    !==========================================================================!

    use rundat, only: pub_debug
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    integer, intent(in) :: l
    integer, intent(in) :: m
    type(SPHERICAL_MULTIPOLE_SET), intent(in) :: mpole_set

    ! jd: Local variables
    integer :: min_l, max_l
    character(len=*), parameter :: myself = 'multipole_lm_to_offset'

    ! -------------------------------------------------------------------------

    min_l = mpole_set%min_l
    max_l = mpole_set%max_l

    if(pub_debug) then
       call utils_assert(l >= min_l, myself//': l must be >= min_l', l, min_l)
       call utils_assert(l <= max_l, myself//': l must be <= max_l', l, max_l)
       call utils_assert(m >= -l .and. m<= l, myself//': m out of range', m, l)
    end if

    multipole_lm_to_offset = 1 + l + l*l - min_l*min_l + m

    if(pub_debug) then
       call utils_assert(multipole_lm_to_offset >= 1 .and. &
            multipole_lm_to_offset <= mpole_set%n_multipoles_per_site, &
            myself//': offset out of range', &
            multipole_lm_to_offset, mpole_set%n_multipoles_per_site)
    end if

  end function multipole_lm_to_offset

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_masked_T0(r1vec, r2vec, r0)
    !==========================================================================!
    ! Returns the Coulombic-masked multipole-multipole interaction tensor      !
    ! of order 0 (charge-charge interaction).                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    ! Return value:                                                            !
    !   T0(r1vec-r2vec), a scalar.                                             !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   If you pass r0 = 0 (or anything < SAFE_DIV_EPS), expect infinity when  !
    !   r1vec coincides with r2vec.                                            !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)

    implicit none

    ! jd: Return value
    real(kind=DP) :: multipole_coul_masked_T0

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: r_norm

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       multipole_coul_masked_T0 = 0D0
    else
       multipole_coul_masked_T0 = 1D0/r_norm
    end if

  end function multipole_coul_masked_T0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_masked_T1(r1vec, r2vec, r0)
    !==========================================================================!
    ! Returns the Coulombic-masked multipole-multipole interaction tensor      !
    ! of order 1 (charge-dipole interaction).                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    ! Return value:                                                            !
    !   T1(r1vec-r2vec), a 3-vector.                                           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   If you pass r0 = 0 (or anything < SAFE_DIV_EPS), expect infinity or    !
    !   0/0 when r1vec coincides with r2vec.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3) :: multipole_coul_masked_T1

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: r_norm

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       multipole_coul_masked_T1 = rank1zero
    else
       multipole_coul_masked_T1(1) = -dr_point%x / r_norm**3
       multipole_coul_masked_T1(2) = -dr_point%y / r_norm**3
       multipole_coul_masked_T1(3) = -dr_point%z / r_norm**3
    end if

  end function multipole_coul_masked_T1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_masked_T2(r1vec, r2vec, r0)
    !==========================================================================!
    ! Returns the Coulombic-masked multipole-multipole interaction tensor      !
    ! of order 2 (charge-quadrupole or dipole-dipole interaction).             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    ! Return value:                                                            !
    !   T2(r1vec-r2vec), a 3-matrix.                                           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   If you pass r0 = 0 (or anything < SAFE_DIV_EPS), expect infinity or    !
    !   0/0 when r1vec coincides with r2vec.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3) :: multipole_coul_masked_T2

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: r_norm
    real(kind=DP) :: r2, rm5

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       multipole_coul_masked_T2 = rank2zero
    else
       rm5 = r_norm ** (-5)
       r2 = r_norm ** 2
       multipole_coul_masked_T2(1,1) = &
            (3D0 * dr_point%x * dr_point%x - r2) * rm5
       multipole_coul_masked_T2(2,2) = &
            (3D0 * dr_point%y * dr_point%y - r2) * rm5
       multipole_coul_masked_T2(3,3) = &
            (3D0 * dr_point%z * dr_point%z - r2) * rm5
       multipole_coul_masked_T2(2,1) = (3D0 * dr_point%x * dr_point%y) * rm5
       multipole_coul_masked_T2(3,1) = (3D0 * dr_point%x * dr_point%z) * rm5
       multipole_coul_masked_T2(3,2) = (3D0 * dr_point%y * dr_point%z) * rm5
       multipole_coul_masked_T2(2,3) = multipole_coul_masked_T2(3,2)
       multipole_coul_masked_T2(1,2) = multipole_coul_masked_T2(2,1)
       multipole_coul_masked_T2(1,3) = multipole_coul_masked_T2(3,1)
    end if

  end function multipole_coul_masked_T2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_masked_T3(r1vec, r2vec, r0)
    !==========================================================================!
    ! Returns the Coulombic-masked multipole-multipole interaction tensor      !
    ! of order 3 (charge-octupole or dipole-quadrupole interaction).           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    ! Return value:                                                            !
    !   T3(r1vec-r2vec), a 3-(rank-3 tensor).                                  !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   If you pass r0 = 0 (or anything < SAFE_DIV_EPS), expect infinity or    !
    !   0/0 when r1vec coincides with r2vec.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3,3) :: multipole_coul_masked_T3

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: r_norm
    real(kind=DP) :: rm5, rm7, l5, l7

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       multipole_coul_masked_T3 = rank3zero
    else
       rm5 = r_norm ** (-5.0_DP)
       rm7 = r_norm ** (-7.0_DP)
       l5 = 3.0_DP * rm5
       l7 = 15.0_DP * rm7

       multipole_coul_masked_T3(1,1,1) = &
            (3.0_DP * l5 * dr_point%x) - (l7 * dr_point%x * dr_point%x * dr_point%x)
       multipole_coul_masked_T3(2,1,1) = &
            (l5 * dr_point%y) - (l7 * dr_point%x * dr_point%x * dr_point%y)
       multipole_coul_masked_T3(3,1,1) = &
            (l5 * dr_point%z) - (l7 * dr_point%x * dr_point%x * dr_point%z)
       multipole_coul_masked_T3(2,2,1) = &
            (l5 * dr_point%x) - (l7 * dr_point%x * dr_point%y * dr_point%y)
       multipole_coul_masked_T3(3,2,1) = &
            - (l7 * dr_point%x * dr_point%y * dr_point%z)
       multipole_coul_masked_T3(3,3,1) = &
            (l5 * dr_point%x) - (l7 * dr_point%x * dr_point%z * dr_point%z)
       multipole_coul_masked_T3(2,2,2) = &
            (3.0_DP * l5 * dr_point%y) - (l7 * dr_point%y * dr_point%y * dr_point%y)
       multipole_coul_masked_T3(3,2,2) = &
            (l5 * dr_point%z) - (l7 * dr_point%y * dr_point%y * dr_point%z)
       multipole_coul_masked_T3(3,3,2) = &
            (l5 * dr_point%y) - (l7 * dr_point%y * dr_point%z * dr_point%z)
       multipole_coul_masked_T3(3,3,3) = &
            (3.0_DP * l5 * dr_point%z) - (l7 * dr_point%z * dr_point%z * dr_point%z)

       multipole_coul_masked_T3(1,1,2) = multipole_coul_masked_T3(2,1,1)
       multipole_coul_masked_T3(1,2,1) = multipole_coul_masked_T3(2,1,1)
       multipole_coul_masked_T3(1,1,3) = multipole_coul_masked_T3(3,1,1)
       multipole_coul_masked_T3(1,3,1) = multipole_coul_masked_T3(3,1,1)

       multipole_coul_masked_T3(1,2,2) = multipole_coul_masked_T3(2,2,1)
       multipole_coul_masked_T3(2,1,2) = multipole_coul_masked_T3(2,2,1)

       multipole_coul_masked_T3(1,2,3) = multipole_coul_masked_T3(3,2,1)
       multipole_coul_masked_T3(1,3,2) = multipole_coul_masked_T3(3,2,1)
       multipole_coul_masked_T3(2,1,3) = multipole_coul_masked_T3(3,2,1)
       multipole_coul_masked_T3(2,3,1) = multipole_coul_masked_T3(3,2,1)
       multipole_coul_masked_T3(3,1,2) = multipole_coul_masked_T3(3,2,1)

       multipole_coul_masked_T3(1,3,3) = multipole_coul_masked_T3(3,3,1)
       multipole_coul_masked_T3(3,1,3) = multipole_coul_masked_T3(3,3,1)

       multipole_coul_masked_T3(2,2,3) = multipole_coul_masked_T3(3,2,2)
       multipole_coul_masked_T3(2,3,2) = multipole_coul_masked_T3(3,2,2)

       multipole_coul_masked_T3(2,3,3) = multipole_coul_masked_T3(3,3,2)
       multipole_coul_masked_T3(3,2,3) = multipole_coul_masked_T3(3,3,2)
    end if

  end function multipole_coul_masked_T3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_masked_T4(r1vec, r2vec, r0)
    !==========================================================================!
    ! Returns the Coulombic-masked multipole-multipole interaction tensor      !
    ! of order 4 (charge-hexadecapole, dipole-octupole or quadrupole-quadrupole!
    ! interaction).                                                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    ! Return value:                                                            !
    !   T4(r1vec-r2vec), a 3-(rank-4 tensor).                                  !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   If you pass r0 = 0 (or anything < SAFE_DIV_EPS), expect infinity or    !
    !   0/0 when r1vec coincides with r2vec.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3,3,3) :: multipole_coul_masked_T4

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: rm5, rm7, rm9
    real(kind=DP) :: x, y, z, x2, y2, z2, x3, y3, z3
    real(kind=DP) :: r_norm
    real(kind=DP) :: a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s
    real(kind=DP) :: u, v, w, aa, cc

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       multipole_coul_masked_T4 = rank4zero
    else
       rm5 = r_norm ** (-5.0_DP)
       rm7 = r_norm ** (-7.0_DP)
       rm9 = r_norm ** (-9.0_DP)
       x = dr_point%x
       y = dr_point%y
       z = dr_point%z
       x2 = x*x
       y2 = y*y
       z2 = z*z
       x3 = x * x2
       y3 = y * y2
       z3 = z * z2

       a = -((45.0_DP * x * y)*rm7)
       b = (105.0_DP * x3 * y)*rm9
       c = -((45.0_DP * x * z)*rm7)
       d = (105.0_DP * x3 * z)*rm9
       e = -((15.0_DP * y * z)*rm7)
       f = (105.0_DP * x2 * y * z)*rm9
       g = 3.0_DP*rm5
       h = -((15.0_DP * x2)*rm7)
       i = -((15.0_DP * y2)*rm7)
       j = (105.0_DP * x2 * y2)*rm9
       k = (105.0_DP * x * y3)*rm9
       l = -((15.0_DP * x * z)*rm7)
       m = ( 105.0_DP * x * y2 * z)*rm9
       n = -((15.0_DP * x * y)*rm7)
       o = (105.0_DP * x * y * z2)*rm9
       p = -((15.0_DP * z2)*rm7)
       q = (105.0_DP * x2 * z2)*rm9
       r = (105.0_DP * x * z3)*rm9
       s = -((45.0_DP * y * z)*rm7)
       u = (105.0_DP * y3 * z)*rm9
       v = (105.0_DP * y2 * z2)*rm9
       w = (105.0_DP * y * z3)*rm9
       aa = 9.0_DP*rm5
       cc = 105.0_DP *rm9

       ! --- block 1,1 ---
       multipole_coul_masked_T4(1,1,1,1) = aa - 90.0_DP * x2*rm7 + cc*x2*x2
       multipole_coul_masked_T4(2,1,1,1) = a + b
       multipole_coul_masked_T4(3,1,1,1) = c + d
       multipole_coul_masked_T4(1,2,1,1) = multipole_coul_masked_T4(2,1,1,1)
       multipole_coul_masked_T4(2,2,1,1) = g + h + i + j
       multipole_coul_masked_T4(3,2,1,1) = e + f
       multipole_coul_masked_T4(1,3,1,1) = multipole_coul_masked_T4(3,1,1,1)
       multipole_coul_masked_T4(2,3,1,1) = multipole_coul_masked_T4(3,2,1,1)
       multipole_coul_masked_T4(3,3,1,1) = g + h + p + q

       ! --- block 2,1 ---
       multipole_coul_masked_T4(1,1,2,1) = a + b
       multipole_coul_masked_T4(2,1,2,1) = g + h + i + j
       multipole_coul_masked_T4(3,1,2,1) = e + f
       multipole_coul_masked_T4(1,2,2,1) = multipole_coul_masked_T4(2,1,2,1)
       multipole_coul_masked_T4(2,2,2,1) = a + k
       multipole_coul_masked_T4(3,2,2,1) = l + m
       multipole_coul_masked_T4(1,3,2,1) = multipole_coul_masked_T4(3,1,2,1)
       multipole_coul_masked_T4(2,3,2,1) = multipole_coul_masked_T4(3,2,2,1)
       multipole_coul_masked_T4(3,3,2,1) = n + o

       ! --- block 3,1 ---
       multipole_coul_masked_T4(1,1,3,1) = c + d
       multipole_coul_masked_T4(2,1,3,1) = e + f
       multipole_coul_masked_T4(3,1,3,1) = g + h + p + q
       multipole_coul_masked_T4(1,2,3,1) = multipole_coul_masked_T4(2,1,3,1)
       multipole_coul_masked_T4(2,2,3,1) = l + m
       multipole_coul_masked_T4(3,2,3,1) = n + o
       multipole_coul_masked_T4(1,3,3,1) = multipole_coul_masked_T4(3,1,3,1)
       multipole_coul_masked_T4(2,3,3,1) = multipole_coul_masked_T4(3,2,3,1)
       multipole_coul_masked_T4(3,3,3,1) = c + r

       ! --- block 1,2 ---
       multipole_coul_masked_T4(1,1,1,2) = multipole_coul_masked_T4(1,1,2,1)
       multipole_coul_masked_T4(2,1,1,2) = multipole_coul_masked_T4(2,1,2,1)
       multipole_coul_masked_T4(3,1,1,2) = multipole_coul_masked_T4(3,1,2,1)
       multipole_coul_masked_T4(1,2,1,2) = multipole_coul_masked_T4(1,2,2,1)
       multipole_coul_masked_T4(2,2,1,2) = multipole_coul_masked_T4(2,2,2,1)
       multipole_coul_masked_T4(3,2,1,2) = multipole_coul_masked_T4(3,2,2,1)
       multipole_coul_masked_T4(1,3,1,2) = multipole_coul_masked_T4(1,3,2,1)
       multipole_coul_masked_T4(2,3,1,2) = multipole_coul_masked_T4(2,3,2,1)
       multipole_coul_masked_T4(3,3,1,2) = multipole_coul_masked_T4(3,3,2,1)

       ! --- block 2,2 ---
       multipole_coul_masked_T4(1,1,2,2) = g + h + i + j
       multipole_coul_masked_T4(2,1,2,2) = a + k
       multipole_coul_masked_T4(3,1,2,2) = l + m
       multipole_coul_masked_T4(1,2,2,2) = multipole_coul_masked_T4(2,1,2,2)
       multipole_coul_masked_T4(2,2,2,2) = aa - 90.0_DP * y2*rm7 + cc*y2*y2
       multipole_coul_masked_T4(3,2,2,2) = s + u
       multipole_coul_masked_T4(1,3,2,2) = multipole_coul_masked_T4(3,1,2,2)
       multipole_coul_masked_T4(2,3,2,2) = multipole_coul_masked_T4(3,2,2,2)
       multipole_coul_masked_T4(3,3,2,2) = g + i + p + v

       ! --- block 3,2 ---
       multipole_coul_masked_T4(1,1,3,2) = e + f
       multipole_coul_masked_T4(2,1,3,2) = l + m
       multipole_coul_masked_T4(3,1,3,2) = n + o
       multipole_coul_masked_T4(1,2,3,2) = multipole_coul_masked_T4(2,1,3,2)
       multipole_coul_masked_T4(2,2,3,2) = s + u
       multipole_coul_masked_T4(3,2,3,2) = g + i + p + v
       multipole_coul_masked_T4(1,3,3,2) = multipole_coul_masked_T4(3,1,3,2)
       multipole_coul_masked_T4(2,3,3,2) = multipole_coul_masked_T4(3,2,3,2)
       multipole_coul_masked_T4(3,3,3,2) = s + w

       ! --- block 1,3 ---
       multipole_coul_masked_T4(1,1,1,3) = multipole_coul_masked_T4(1,1,3,1)
       multipole_coul_masked_T4(2,1,1,3) = multipole_coul_masked_T4(2,1,3,1)
       multipole_coul_masked_T4(3,1,1,3) = multipole_coul_masked_T4(3,1,3,1)
       multipole_coul_masked_T4(1,2,1,3) = multipole_coul_masked_T4(1,2,3,1)
       multipole_coul_masked_T4(2,2,1,3) = multipole_coul_masked_T4(2,2,3,1)
       multipole_coul_masked_T4(3,2,1,3) = multipole_coul_masked_T4(3,2,3,1)
       multipole_coul_masked_T4(1,3,1,3) = multipole_coul_masked_T4(1,3,3,1)
       multipole_coul_masked_T4(2,3,1,3) = multipole_coul_masked_T4(2,3,3,1)
       multipole_coul_masked_T4(3,3,1,3) = multipole_coul_masked_T4(3,3,3,1)

       ! --- block 2,3 ---
       multipole_coul_masked_T4(1,1,2,3) = multipole_coul_masked_T4(1,1,3,2)
       multipole_coul_masked_T4(2,1,2,3) = multipole_coul_masked_T4(2,1,3,2)
       multipole_coul_masked_T4(3,1,2,3) = multipole_coul_masked_T4(3,1,3,2)
       multipole_coul_masked_T4(1,2,2,3) = multipole_coul_masked_T4(1,2,3,2)
       multipole_coul_masked_T4(2,2,2,3) = multipole_coul_masked_T4(2,2,3,2)
       multipole_coul_masked_T4(3,2,2,3) = multipole_coul_masked_T4(3,2,3,2)
       multipole_coul_masked_T4(1,3,2,3) = multipole_coul_masked_T4(1,3,3,2)
       multipole_coul_masked_T4(2,3,2,3) = multipole_coul_masked_T4(2,3,3,2)
       multipole_coul_masked_T4(3,3,2,3) = multipole_coul_masked_T4(3,3,3,2)

       ! --- block 3,3 ---
       multipole_coul_masked_T4(1,1,3,3) = g + h + p + q
       multipole_coul_masked_T4(2,1,3,3) = n + o
       multipole_coul_masked_T4(3,1,3,3) = c + r
       multipole_coul_masked_T4(1,2,3,3) = multipole_coul_masked_T4(2,1,3,3)
       multipole_coul_masked_T4(2,2,3,3) = g + i + p + v
       multipole_coul_masked_T4(3,2,3,3) = s + w
       multipole_coul_masked_T4(1,3,3,3) = multipole_coul_masked_T4(3,1,3,3)
       multipole_coul_masked_T4(2,3,3,3) = multipole_coul_masked_T4(3,2,3,3)
       multipole_coul_masked_T4(3,3,3,3) = aa - 90.0_DP * z2*rm7 + cc*z2*z2

    end if

  end function multipole_coul_masked_T4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_smeared_T0(r1vec, r2vec, r0, small_a)
    !==========================================================================!
    ! Returns the Coulombic-smeared multipole-multipole interaction tensor     !
    ! of order 0 (charge-charge interaction).                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T0(r1vec-r2vec), a scalar.                                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)
    use rundat, only: pub_pol_emb_smearing_a

    implicit none

    ! jd: Return value
    real(kind=DP) :: multipole_coul_smeared_T0

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0
    real(kind=DP), intent(in)           :: small_a

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: r_norm

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       ! jd: Emulate smearing from damping
       multipole_coul_smeared_T0 = multipole_thole_T0(r1vec, r2vec, 1D0, &
            pub_pol_emb_smearing_a**6, small_a)
    else
       multipole_coul_smeared_T0 = 1D0/r_norm
    end if

  end function multipole_coul_smeared_T0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_smeared_T1(r1vec, r2vec, r0, small_a)
    !==========================================================================!
    ! Returns the Coulombic-smeared multipole-multipole interaction tensor     !
    ! of order 1 (charge-dipole interaction).                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T1(r1vec-r2vec), a 3-vector.                                           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   If you pass r0 = 0 (or anything < SAFE_DIV_EPS), expect infinity or    !
    !   0/0 when r1vec coincides with r2vec.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)
    use rundat, only: pub_pol_emb_smearing_a

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3) :: multipole_coul_smeared_T1

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0
    real(kind=DP), intent(in)           :: small_a

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: r_norm

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       ! jd: Emulate smearing from damping
       multipole_coul_smeared_T1 = multipole_thole_T1(r1vec, r2vec, 1D0, &
            pub_pol_emb_smearing_a**6, small_a)
    else
       multipole_coul_smeared_T1(1) = -dr_point%x / r_norm**3
       multipole_coul_smeared_T1(2) = -dr_point%y / r_norm**3
       multipole_coul_smeared_T1(3) = -dr_point%z / r_norm**3
    end if

  end function multipole_coul_smeared_T1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_smeared_T2(r1vec, r2vec, r0, small_a)
    !==========================================================================!
    ! Returns the Coulombic-smeared multipole-multipole interaction tensor      !
    ! of order 2 (charge-quadrupole or dipole-dipole interaction).             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T2(r1vec-r2vec), a 3-matrix.                                           !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)
    use rundat, only: pub_pol_emb_smearing_a

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3) :: multipole_coul_smeared_T2

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0
    real(kind=DP), intent(in)           :: small_a

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: r_norm
    real(kind=DP) :: r2, rm5

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       ! jd: Emulate smearing from damping
       multipole_coul_smeared_T2 = multipole_thole_T2(r1vec, r2vec, 1D0, &
            pub_pol_emb_smearing_a**6, small_a)
    else
       rm5 = r_norm ** (-5)
       r2 = r_norm ** 2
       multipole_coul_smeared_T2(1,1) = &
            (3D0 * dr_point%x * dr_point%x - r2) * rm5
       multipole_coul_smeared_T2(2,2) = &
            (3D0 * dr_point%y * dr_point%y - r2) * rm5
       multipole_coul_smeared_T2(3,3) = &
            (3D0 * dr_point%z * dr_point%z - r2) * rm5
       multipole_coul_smeared_T2(2,1) = (3D0 * dr_point%x * dr_point%y) * rm5
       multipole_coul_smeared_T2(3,1) = (3D0 * dr_point%x * dr_point%z) * rm5
       multipole_coul_smeared_T2(3,2) = (3D0 * dr_point%y * dr_point%z) * rm5
       multipole_coul_smeared_T2(2,3) = multipole_coul_smeared_T2(3,2)
       multipole_coul_smeared_T2(1,2) = multipole_coul_smeared_T2(2,1)
       multipole_coul_smeared_T2(1,3) = multipole_coul_smeared_T2(3,1)
    end if

  end function multipole_coul_smeared_T2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_smeared_T3(r1vec, r2vec, r0, small_a)
    !==========================================================================!
    ! Returns the Coulombic-smeared multipole-multipole interaction tensor      !
    ! of order 3 (charge-octupole or dipole-quadrupole interaction).           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T3(r1vec-r2vec), a 3-(rank-3 tensor).                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)
    use rundat, only: pub_pol_emb_smearing_a

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3,3) :: multipole_coul_smeared_T3

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0
    real(kind=DP), intent(in)           :: small_a

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: r_norm
    real(kind=DP) :: rm5, rm7, l5, l7

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       ! jd: Emulate smearing from damping
       multipole_coul_smeared_T3 = multipole_thole_T3(r1vec, r2vec, 1D0, &
            pub_pol_emb_smearing_a**6, small_a)
    else
       rm5 = r_norm ** (-5.0_DP)
       rm7 = r_norm ** (-7.0_DP)
       l5 = 3.0_DP * rm5
       l7 = 15.0_DP * rm7

       multipole_coul_smeared_T3(1,1,1) = &
            (3.0_DP * l5 * dr_point%x) - (l7 * dr_point%x * dr_point%x * dr_point%x)
       multipole_coul_smeared_T3(2,1,1) = &
            (l5 * dr_point%y) - (l7 * dr_point%x * dr_point%x * dr_point%y)
       multipole_coul_smeared_T3(3,1,1) = &
            (l5 * dr_point%z) - (l7 * dr_point%x * dr_point%x * dr_point%z)
       multipole_coul_smeared_T3(2,2,1) = &
            (l5 * dr_point%x) - (l7 * dr_point%x * dr_point%y * dr_point%y)
       multipole_coul_smeared_T3(3,2,1) = &
            - (l7 * dr_point%x * dr_point%y * dr_point%z)
       multipole_coul_smeared_T3(3,3,1) = &
            (l5 * dr_point%x) - (l7 * dr_point%x * dr_point%z * dr_point%z)
       multipole_coul_smeared_T3(2,2,2) = &
            (3.0_DP * l5 * dr_point%y) - (l7 * dr_point%y * dr_point%y * dr_point%y)
       multipole_coul_smeared_T3(3,2,2) = &
            (l5 * dr_point%z) - (l7 * dr_point%y * dr_point%y * dr_point%z)
       multipole_coul_smeared_T3(3,3,2) = &
            (l5 * dr_point%y) - (l7 * dr_point%y * dr_point%z * dr_point%z)
       multipole_coul_smeared_T3(3,3,3) = &
            (3.0_DP * l5 * dr_point%z) - (l7 * dr_point%z * dr_point%z * dr_point%z)

       multipole_coul_smeared_T3(1,1,2) = multipole_coul_smeared_T3(2,1,1)
       multipole_coul_smeared_T3(1,2,1) = multipole_coul_smeared_T3(2,1,1)
       multipole_coul_smeared_T3(1,1,3) = multipole_coul_smeared_T3(3,1,1)
       multipole_coul_smeared_T3(1,3,1) = multipole_coul_smeared_T3(3,1,1)

       multipole_coul_smeared_T3(1,2,2) = multipole_coul_smeared_T3(2,2,1)
       multipole_coul_smeared_T3(2,1,2) = multipole_coul_smeared_T3(2,2,1)

       multipole_coul_smeared_T3(1,2,3) = multipole_coul_smeared_T3(3,2,1)
       multipole_coul_smeared_T3(1,3,2) = multipole_coul_smeared_T3(3,2,1)
       multipole_coul_smeared_T3(2,1,3) = multipole_coul_smeared_T3(3,2,1)
       multipole_coul_smeared_T3(2,3,1) = multipole_coul_smeared_T3(3,2,1)
       multipole_coul_smeared_T3(3,1,2) = multipole_coul_smeared_T3(3,2,1)

       multipole_coul_smeared_T3(1,3,3) = multipole_coul_smeared_T3(3,3,1)
       multipole_coul_smeared_T3(3,1,3) = multipole_coul_smeared_T3(3,3,1)

       multipole_coul_smeared_T3(2,2,3) = multipole_coul_smeared_T3(3,2,2)
       multipole_coul_smeared_T3(2,3,2) = multipole_coul_smeared_T3(3,2,2)

       multipole_coul_smeared_T3(2,3,3) = multipole_coul_smeared_T3(3,3,2)
       multipole_coul_smeared_T3(3,2,3) = multipole_coul_smeared_T3(3,3,2)
    end if

  end function multipole_coul_smeared_T3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_coul_smeared_T4(r1vec, r2vec, r0, small_a)
    !==========================================================================!
    ! Returns the Coulombic-smeared multipole-multipole interaction tensor      !
    ! of order 4 (charge-hexadecapole, dipole-octupole or quadrupole-quadrupole!
    ! interaction).                                                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   r0 (in):    Cut-off radius below which masking to zero occurs.         !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T4(r1vec-r2vec), a 3-(rank-4 tensor).                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use geometry, only: POINT, magnitude, operator(-)
    use rundat, only: pub_pol_emb_smearing_a

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3,3,3) :: multipole_coul_smeared_T4

    ! jd: Arguments
    type(POINT), intent(in)             :: r1vec, r2vec
    real(kind=DP), intent(in)           :: r0
    real(kind=DP), intent(in)           :: small_a

    ! jd: Local variables
    type(POINT)   :: dr_point
    real(kind=DP) :: rm5, rm7, rm9
    real(kind=DP) :: x, y, z, x2, y2, z2, x3, y3, z3
    real(kind=DP) :: r_norm
    real(kind=DP) :: a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s
    real(kind=DP) :: u, v, w, aa, cc

    ! -------------------------------------------------------------------------

    dr_point = r2vec - r1vec
    r_norm = magnitude(dr_point) ! @MIC must be supported later on

    if(r_norm <= r0) then
       ! jd: Emulate smearing from damping
       multipole_coul_smeared_T4 = multipole_thole_T4(r1vec, r2vec, 1D0, &
            pub_pol_emb_smearing_a**6, small_a)
    else
       rm5 = r_norm ** (-5.0_DP)
       rm7 = r_norm ** (-7.0_DP)
       rm9 = r_norm ** (-9.0_DP)
       x = dr_point%x
       y = dr_point%y
       z = dr_point%z
       x2 = x*x
       y2 = y*y
       z2 = z*z
       x3 = x * x2
       y3 = y * y2
       z3 = z * z2

       a = -((45.0_DP * x * y)*rm7)
       b = (105.0_DP * x3 * y)*rm9
       c = -((45.0_DP * x * z)*rm7)
       d = (105.0_DP * x3 * z)*rm9
       e = -((15.0_DP * y * z)*rm7)
       f = (105.0_DP * x2 * y * z)*rm9
       g = 3.0_DP*rm5
       h = -((15.0_DP * x2)*rm7)
       i = -((15.0_DP * y2)*rm7)
       j = (105.0_DP * x2 * y2)*rm9
       k = (105.0_DP * x * y3)*rm9
       l = -((15.0_DP * x * z)*rm7)
       m = ( 105.0_DP * x * y2 * z)*rm9
       n = -((15.0_DP * x * y)*rm7)
       o = (105.0_DP * x * y * z2)*rm9
       p = -((15.0_DP * z2)*rm7)
       q = (105.0_DP * x2 * z2)*rm9
       r = (105.0_DP * x * z3)*rm9
       s = -((45.0_DP * y * z)*rm7)
       u = (105.0_DP * y3 * z)*rm9
       v = (105.0_DP * y2 * z2)*rm9
       w = (105.0_DP * y * z3)*rm9
       aa = 9.0_DP*rm5
       cc = 105.0_DP *rm9

       ! --- block 1,1 ---
       multipole_coul_smeared_T4(1,1,1,1) = aa - 90.0_DP * x2*rm7 + cc*x2*x2
       multipole_coul_smeared_T4(2,1,1,1) = a + b
       multipole_coul_smeared_T4(3,1,1,1) = c + d
       multipole_coul_smeared_T4(1,2,1,1) = multipole_coul_smeared_T4(2,1,1,1)
       multipole_coul_smeared_T4(2,2,1,1) = g + h + i + j
       multipole_coul_smeared_T4(3,2,1,1) = e + f
       multipole_coul_smeared_T4(1,3,1,1) = multipole_coul_smeared_T4(3,1,1,1)
       multipole_coul_smeared_T4(2,3,1,1) = multipole_coul_smeared_T4(3,2,1,1)
       multipole_coul_smeared_T4(3,3,1,1) = g + h + p + q

       ! --- block 2,1 ---
       multipole_coul_smeared_T4(1,1,2,1) = a + b
       multipole_coul_smeared_T4(2,1,2,1) = g + h + i + j
       multipole_coul_smeared_T4(3,1,2,1) = e + f
       multipole_coul_smeared_T4(1,2,2,1) = multipole_coul_smeared_T4(2,1,2,1)
       multipole_coul_smeared_T4(2,2,2,1) = a + k
       multipole_coul_smeared_T4(3,2,2,1) = l + m
       multipole_coul_smeared_T4(1,3,2,1) = multipole_coul_smeared_T4(3,1,2,1)
       multipole_coul_smeared_T4(2,3,2,1) = multipole_coul_smeared_T4(3,2,2,1)
       multipole_coul_smeared_T4(3,3,2,1) = n + o

       ! --- block 3,1 ---
       multipole_coul_smeared_T4(1,1,3,1) = c + d
       multipole_coul_smeared_T4(2,1,3,1) = e + f
       multipole_coul_smeared_T4(3,1,3,1) = g + h + p + q
       multipole_coul_smeared_T4(1,2,3,1) = multipole_coul_smeared_T4(2,1,3,1)
       multipole_coul_smeared_T4(2,2,3,1) = l + m
       multipole_coul_smeared_T4(3,2,3,1) = n + o
       multipole_coul_smeared_T4(1,3,3,1) = multipole_coul_smeared_T4(3,1,3,1)
       multipole_coul_smeared_T4(2,3,3,1) = multipole_coul_smeared_T4(3,2,3,1)
       multipole_coul_smeared_T4(3,3,3,1) = c + r

       ! --- block 1,2 ---
       multipole_coul_smeared_T4(1,1,1,2) = multipole_coul_smeared_T4(1,1,2,1)
       multipole_coul_smeared_T4(2,1,1,2) = multipole_coul_smeared_T4(2,1,2,1)
       multipole_coul_smeared_T4(3,1,1,2) = multipole_coul_smeared_T4(3,1,2,1)
       multipole_coul_smeared_T4(1,2,1,2) = multipole_coul_smeared_T4(1,2,2,1)
       multipole_coul_smeared_T4(2,2,1,2) = multipole_coul_smeared_T4(2,2,2,1)
       multipole_coul_smeared_T4(3,2,1,2) = multipole_coul_smeared_T4(3,2,2,1)
       multipole_coul_smeared_T4(1,3,1,2) = multipole_coul_smeared_T4(1,3,2,1)
       multipole_coul_smeared_T4(2,3,1,2) = multipole_coul_smeared_T4(2,3,2,1)
       multipole_coul_smeared_T4(3,3,1,2) = multipole_coul_smeared_T4(3,3,2,1)

       ! --- block 2,2 ---
       multipole_coul_smeared_T4(1,1,2,2) = g + h + i + j
       multipole_coul_smeared_T4(2,1,2,2) = a + k
       multipole_coul_smeared_T4(3,1,2,2) = l + m
       multipole_coul_smeared_T4(1,2,2,2) = multipole_coul_smeared_T4(2,1,2,2)
       multipole_coul_smeared_T4(2,2,2,2) = aa - 90.0_DP * y2*rm7 + cc*y2*y2
       multipole_coul_smeared_T4(3,2,2,2) = s + u
       multipole_coul_smeared_T4(1,3,2,2) = multipole_coul_smeared_T4(3,1,2,2)
       multipole_coul_smeared_T4(2,3,2,2) = multipole_coul_smeared_T4(3,2,2,2)
       multipole_coul_smeared_T4(3,3,2,2) = g + i + p + v

       ! --- block 3,2 ---
       multipole_coul_smeared_T4(1,1,3,2) = e + f
       multipole_coul_smeared_T4(2,1,3,2) = l + m
       multipole_coul_smeared_T4(3,1,3,2) = n + o
       multipole_coul_smeared_T4(1,2,3,2) = multipole_coul_smeared_T4(2,1,3,2)
       multipole_coul_smeared_T4(2,2,3,2) = s + u
       multipole_coul_smeared_T4(3,2,3,2) = g + i + p + v
       multipole_coul_smeared_T4(1,3,3,2) = multipole_coul_smeared_T4(3,1,3,2)
       multipole_coul_smeared_T4(2,3,3,2) = multipole_coul_smeared_T4(3,2,3,2)
       multipole_coul_smeared_T4(3,3,3,2) = s + w

       ! --- block 1,3 ---
       multipole_coul_smeared_T4(1,1,1,3) = multipole_coul_smeared_T4(1,1,3,1)
       multipole_coul_smeared_T4(2,1,1,3) = multipole_coul_smeared_T4(2,1,3,1)
       multipole_coul_smeared_T4(3,1,1,3) = multipole_coul_smeared_T4(3,1,3,1)
       multipole_coul_smeared_T4(1,2,1,3) = multipole_coul_smeared_T4(1,2,3,1)
       multipole_coul_smeared_T4(2,2,1,3) = multipole_coul_smeared_T4(2,2,3,1)
       multipole_coul_smeared_T4(3,2,1,3) = multipole_coul_smeared_T4(3,2,3,1)
       multipole_coul_smeared_T4(1,3,1,3) = multipole_coul_smeared_T4(1,3,3,1)
       multipole_coul_smeared_T4(2,3,1,3) = multipole_coul_smeared_T4(2,3,3,1)
       multipole_coul_smeared_T4(3,3,1,3) = multipole_coul_smeared_T4(3,3,3,1)

       ! --- block 2,3 ---
       multipole_coul_smeared_T4(1,1,2,3) = multipole_coul_smeared_T4(1,1,3,2)
       multipole_coul_smeared_T4(2,1,2,3) = multipole_coul_smeared_T4(2,1,3,2)
       multipole_coul_smeared_T4(3,1,2,3) = multipole_coul_smeared_T4(3,1,3,2)
       multipole_coul_smeared_T4(1,2,2,3) = multipole_coul_smeared_T4(1,2,3,2)
       multipole_coul_smeared_T4(2,2,2,3) = multipole_coul_smeared_T4(2,2,3,2)
       multipole_coul_smeared_T4(3,2,2,3) = multipole_coul_smeared_T4(3,2,3,2)
       multipole_coul_smeared_T4(1,3,2,3) = multipole_coul_smeared_T4(1,3,3,2)
       multipole_coul_smeared_T4(2,3,2,3) = multipole_coul_smeared_T4(2,3,3,2)
       multipole_coul_smeared_T4(3,3,2,3) = multipole_coul_smeared_T4(3,3,3,2)

       ! --- block 3,3 ---
       multipole_coul_smeared_T4(1,1,3,3) = g + h + p + q
       multipole_coul_smeared_T4(2,1,3,3) = n + o
       multipole_coul_smeared_T4(3,1,3,3) = c + r
       multipole_coul_smeared_T4(1,2,3,3) = multipole_coul_smeared_T4(2,1,3,3)
       multipole_coul_smeared_T4(2,2,3,3) = g + i + p + v
       multipole_coul_smeared_T4(3,2,3,3) = s + w
       multipole_coul_smeared_T4(1,3,3,3) = multipole_coul_smeared_T4(3,1,3,3)
       multipole_coul_smeared_T4(2,3,3,3) = multipole_coul_smeared_T4(3,2,3,3)
       multipole_coul_smeared_T4(3,3,3,3) = aa - 90.0_DP * z2*rm7 + cc * z2*z2

    end if

  end function multipole_coul_smeared_T4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_thole_T0(r1vec, r2vec, alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the Thole-damped interaction tensor of order 0.                  !
    ! (charge-charge interaction).                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   alpha1 (in): The polarisability of the first site.                     !
    !   alpha2 (in): The polarisability of the second site.                    !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T0(r1vec-r2vec), a scalar.                                             !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   For very short distances the approximation in \Gamma(a,x) breaks down  !
    !   leading to somewhat inaccurate results. To avoid this, whenever        !
    !   |r1vec-r2vec| is below tiny_r (currently 0.05a0), we return the value  !
    !   for r=0. This makes the hump in T0(r) flat immediately around r=0, with!
    !   no discontiuities in T0, but with a discontinuity in the derivative.   !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2015.                              !
    !==========================================================================!

    use constants, only: ONE_SIXTH, ONE_THIRD, TWO_THIRDS
    use geometry, only: POINT, magnitude, operator(-)
    use services, only: services_incomplete_gamma_from_lookup
    use utils, only: utils_assert

    implicit none

    ! jd: Return value
    real(kind=DP) :: multipole_thole_T0

    ! jd: Arguments
    type(POINT), intent(in)   :: r1vec, r2vec
    real(kind=DP), intent(in) :: alpha1, alpha2
    real(kind=DP), intent(in) :: small_a

    ! jd: Local variables
    real(kind=DP) :: r_norm
    real(kind=DP) :: big_A
    type(POINT)   :: rvec
    real(kind=DP), parameter :: tiny_r = 0.05

    ! -------------------------------------------------------------------------

    rvec = r2vec-r1vec ! @MIC
    r_norm = magnitude(rvec)
    call utils_assert(alpha1 > 0D0, 'multipole_thole_T0: Illegal QM-side &
         &polarisability', alpha1)
    call utils_assert(alpha2 > 0D0, 'multipole_thole_T0: Illegal MM-side &
         &polarisability', alpha2)

    big_A = (alpha1*alpha2)**ONE_SIXTH

    if(r_norm > tiny_r) then
       multipole_thole_T0 = &
            multipole_thole_damping(r_norm, small_a, big_A, 0) / r_norm
    else
       multipole_thole_T0 = &
            (small_a**ONE_THIRD * &
            services_incomplete_gamma_from_lookup(TWO_THIRDS, 0D0, &
            gamma_two_thirds_lookup)) / big_A
    end if

  end function multipole_thole_T0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_thole_T1(r1vec, r2vec, alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the Thole-damped interaction tensor of order 1.                  !
    ! (charge-dipole interaction).                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   alpha1 (in): The polarisability of the first site.                     !
    !   alpha2 (in): The polarisability of the second site.                    !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T1(r1vec-r2vec), a 3-vector.                                           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   For very short distances (much shorter than for T0), the exp() factor  !
    !   in the Thole expression becomes exactly 1.0, making the result         !
    !   slightly inaccurate in the vicinity of 0 (closer than 1E-4 a0).        !
    !   To avoid this, for these short distances we return the limit value,    !
    !   which is exactly zero.                                                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in November 2015.                              !
    !==========================================================================!

    use constants, only: ONE_SIXTH
    use geometry, only: POINT, magnitude, operator(-)
    use utils, only: utils_assert

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3) :: multipole_thole_T1

    ! jd: Arguments
    type(POINT), intent(in)   :: r1vec, r2vec
    real(kind=DP), intent(in) :: alpha1, alpha2
    real(kind=DP), intent(in) :: small_a

    ! jd: Local variables
    real(kind=DP) :: r_norm
    real(kind=DP) :: big_A
    real(kind=DP) :: lambda3_r3
    type(POINT)   :: rvec
    real(kind=DP), parameter :: tiny_r = 1D-4

    ! -------------------------------------------------------------------------

    rvec = r2vec-r1vec ! @MIC
    r_norm = magnitude(rvec)
    call utils_assert(alpha1 > 0D0, 'multipole_thole_T1: Illegal QM-side &
         &polarisability', alpha1)
    call utils_assert(alpha2 > 0D0, 'multipole_thole_T1: Illegal MM-side &
         &polarisability', alpha2)

    big_A = (alpha1*alpha2)**ONE_SIXTH

    if(r_norm > tiny_r) then
       lambda3_r3 = multipole_thole_damping(r_norm,small_a,big_A, 1) / r_norm**3
       multipole_thole_T1(1) = - rvec%x * lambda3_r3
       multipole_thole_T1(2) = - rvec%y * lambda3_r3
       multipole_thole_T1(3) = - rvec%z * lambda3_r3
    else
       multipole_thole_T1(:) = &
          multipole_thole_damping(0D0,small_a,big_A,1) ! which is 0.0
    end if

  end function multipole_thole_T1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_thole_T2(r1vec, r2vec, alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the Thole-damped interaction tensor of order 2.                  !
    ! (charge-quadrupole or dipole-dipole interaction).                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   alpha1 (in): The polarisability of the first site.                     !
    !   alpha2 (in): The polarisability of the second site.                    !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T2(r1vec-r2vec), a 3-matrix.                                           !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Note that at the limit of r->0 the diagonal terms are non-zero.        !
    !   cf. derivation in thole_field_expressions_3.nb.                        !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in December 2015.                              !
    !==========================================================================!

    use constants, only: ONE_SIXTH
    use geometry, only: POINT, magnitude, operator(-)
    use utils, only: utils_assert

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3) :: multipole_thole_T2

    ! jd: Arguments
    type(POINT), intent(in)   :: r1vec, r2vec
    real(kind=DP), intent(in) :: alpha1, alpha2
    real(kind=DP), intent(in) :: small_a

    ! jd: Local variables
    real(kind=DP) :: r_norm
    real(kind=DP) :: big_A
    real(kind=DP) :: lambda3, lambda5, three_lambda5, lambda3_r2, rm5
    type(POINT)   :: rvec
    real(kind=DP), parameter :: tiny_r = 1D-3

    ! -------------------------------------------------------------------------

    rvec = r2vec-r1vec ! @MIC
    r_norm = magnitude(rvec)
    call utils_assert(alpha1 > 0D0, 'multipole_thole_T2: Illegal QM-side &
         &polarisability', alpha1)
    call utils_assert(alpha2 > 0D0, 'multipole_thole_T2: Illegal MM-side &
         &polarisability', alpha2)

    big_A = (alpha1*alpha2)**ONE_SIXTH

    if(r_norm > tiny_r) then
       lambda3 = multipole_thole_damping(r_norm, small_a, big_A, 1)
       lambda5 = multipole_thole_damping(r_norm, small_a, big_A, 2)
       three_lambda5 = 3.0_DP * lambda5
       lambda3_r2 = lambda3 * r_norm*r_norm
       rm5 = r_norm ** (-5.0_DP)
       multipole_thole_T2(1,1) = &
            (three_lambda5 * rvec%x * rvec%x - lambda3_r2) * rm5
       multipole_thole_T2(2,2) = &
            (three_lambda5 * rvec%y * rvec%y - lambda3_r2) * rm5
       multipole_thole_T2(3,3) = &
            (three_lambda5 * rvec%z * rvec%z - lambda3_r2) * rm5
       multipole_thole_T2(2,1) = (three_lambda5 * rvec%x * rvec%y) * rm5
       multipole_thole_T2(3,1) = (three_lambda5 * rvec%x * rvec%z) * rm5
       multipole_thole_T2(3,2) = (three_lambda5 * rvec%y * rvec%z) * rm5
       multipole_thole_T2(2,3) = multipole_thole_T2(3,2)
       multipole_thole_T2(1,2) = multipole_thole_T2(2,1)
       multipole_thole_T2(1,3) = multipole_thole_T2(3,1)
    else
       ! jd: Limit case
       lambda3 = -small_a / big_A**3.0_DP
       multipole_thole_T2(1,1) = lambda3
       multipole_thole_T2(2,1) = 0D0
       multipole_thole_T2(3,1) = 0D0
       multipole_thole_T2(1,2) = 0D0
       multipole_thole_T2(2,2) = lambda3
       multipole_thole_T2(3,2) = 0D0
       multipole_thole_T2(1,3) = 0D0
       multipole_thole_T2(2,3) = 0D0
       multipole_thole_T2(3,3) = lambda3
    end if

  end function multipole_thole_T2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_thole_T3(r1vec, r2vec, alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the Thole-damped interaction tensor of order 3.                  !
    ! (charge-octupole or dipole-quadrupole interaction).                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   alpha1 (in): The polarisability of the first site.                     !
    !   alpha2 (in): The polarisability of the second site.                    !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T3(r1vec-r2vec), a 3-(rank-3-tensor).                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in December 2015.                              !
    !==========================================================================!

    use constants, only: ONE_SIXTH
    use geometry, only: POINT, magnitude, operator(-)
    use utils, only: utils_assert

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3,3) :: multipole_thole_T3

    ! jd: Arguments
    type(POINT), intent(in)   :: r1vec, r2vec
    real(kind=DP), intent(in) :: alpha1, alpha2
    real(kind=DP), intent(in) :: small_a

    ! jd: Local variables
    real(kind=DP) :: r_norm
    real(kind=DP) :: big_A
    real(kind=DP) :: lambda5, lambda7, rm5, rm7, l5, l7
    type(POINT)   :: rvec
    real(kind=DP), parameter :: tiny_r = 1D-3

    ! -------------------------------------------------------------------------

    rvec = r2vec-r1vec ! @MIC
    r_norm = magnitude(rvec)
    call utils_assert(alpha1 > 0D0, 'multipole_thole_T3: Illegal QM-side &
         &polarisability', alpha1)
    call utils_assert(alpha2 > 0D0, 'multipole_thole_T3: Illegal MM-side &
         &polarisability', alpha2)

    big_A = (alpha1*alpha2)**ONE_SIXTH

    if(r_norm > tiny_r) then
       lambda5 = multipole_thole_damping(r_norm, small_a, big_A, 2)
       lambda7 = multipole_thole_damping(r_norm, small_a, big_A, 3)
       rm5 = r_norm ** (-5.0_DP)
       rm7 = r_norm ** (-7.0_DP)
       l5 = 3.0_DP * lambda5 * rm5
       l7 = 15.0_DP * lambda7 * rm7

       multipole_thole_T3(1,1,1) = &
            (3.0_DP * l5 * rvec%x) - (l7 * rvec%x * rvec%x * rvec%x)
       multipole_thole_T3(2,1,1) = &
            (l5 * rvec%y) - (l7 * rvec%x * rvec%x * rvec%y)
       multipole_thole_T3(3,1,1) = &
            (l5 * rvec%z) - (l7 * rvec%x * rvec%x * rvec%z)
       multipole_thole_T3(2,2,1) = &
            (l5 * rvec%x) - (l7 * rvec%x * rvec%y * rvec%y)
       multipole_thole_T3(3,2,1) = &
            - (l7 * rvec%x * rvec%y * rvec%z)
       multipole_thole_T3(3,3,1) = &
            (l5 * rvec%x) - (l7 * rvec%x * rvec%z * rvec%z)
       multipole_thole_T3(2,2,2) = &
            (3.0_DP * l5 * rvec%y) - (l7 * rvec%y * rvec%y * rvec%y)
       multipole_thole_T3(3,2,2) = &
            (l5 * rvec%z) - (l7 * rvec%y * rvec%y * rvec%z)
       multipole_thole_T3(3,3,2) = &
            (l5 * rvec%y) - (l7 * rvec%y * rvec%z * rvec%z)
       multipole_thole_T3(3,3,3) = &
            (3.0_DP * l5 * rvec%z) - (l7 * rvec%z * rvec%z * rvec%z)

       multipole_thole_T3(1,1,2) = multipole_thole_T3(2,1,1)
       multipole_thole_T3(1,2,1) = multipole_thole_T3(2,1,1)
       multipole_thole_T3(1,1,3) = multipole_thole_T3(3,1,1)
       multipole_thole_T3(1,3,1) = multipole_thole_T3(3,1,1)

       multipole_thole_T3(1,2,2) = multipole_thole_T3(2,2,1)
       multipole_thole_T3(2,1,2) = multipole_thole_T3(2,2,1)

       multipole_thole_T3(1,2,3) = multipole_thole_T3(3,2,1)
       multipole_thole_T3(1,3,2) = multipole_thole_T3(3,2,1)
       multipole_thole_T3(2,1,3) = multipole_thole_T3(3,2,1)
       multipole_thole_T3(2,3,1) = multipole_thole_T3(3,2,1)
       multipole_thole_T3(3,1,2) = multipole_thole_T3(3,2,1)

       multipole_thole_T3(1,3,3) = multipole_thole_T3(3,3,1)
       multipole_thole_T3(3,1,3) = multipole_thole_T3(3,3,1)

       multipole_thole_T3(2,2,3) = multipole_thole_T3(3,2,2)
       multipole_thole_T3(2,3,2) = multipole_thole_T3(3,2,2)

       multipole_thole_T3(2,3,3) = multipole_thole_T3(3,3,2)
       multipole_thole_T3(3,2,3) = multipole_thole_T3(3,3,2)

    else
       ! jd: Limit case
       multipole_thole_T3(:,:,:) = rank3zero
    end if

  end function multipole_thole_T3


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function multipole_thole_T4(r1vec, r2vec, alpha1, alpha2, small_a)
    !==========================================================================!
    ! Returns the Thole-damped interaction tensor of order 4.                  !
    ! (charge-hexadecapole, dipole-octupole or quadrupole-quadrupole interaxn).!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   r1vec (in): The position of the first site.                            !
    !   r2vec (in): The position of the second site.                           !
    !   alpha1 (in): The polarisability of the first site.                     !
    !   alpha2 (in): The polarisability of the second site.                    !
    !   small_a (in): Thole coefficient (usually 0.39).                        !
    ! Return value:                                                            !
    !   T4(r1vec-r2vec), a 3-(rank-4-tensor).                                  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in September 2016.                             !
    !==========================================================================!

    use constants, only: ONE_SIXTH
    use geometry, only: POINT, magnitude, operator(-)
    use utils, only: utils_assert

    implicit none

    ! jd: Return value
    real(kind=DP), dimension(3,3,3,3) :: multipole_thole_T4

    ! jd: Arguments
    type(POINT), intent(in)   :: r1vec, r2vec
    real(kind=DP), intent(in) :: alpha1, alpha2
    real(kind=DP), intent(in) :: small_a

    ! jd: Local variables
    real(kind=DP) :: r_norm
    real(kind=DP) :: big_A
    real(kind=DP) :: lambda5, lambda7, lambda9, rm5, rm7, rm9, l5, l7, l9
    real(kind=DP) :: x, y, z
    type(POINT)   :: rvec
    real(kind=DP), parameter :: tiny_r = 1D-3

    ! -------------------------------------------------------------------------

    rvec = r2vec-r1vec ! @MIC
    r_norm = magnitude(rvec)
    call utils_assert(alpha1 > 0D0, 'multipole_thole_T4: Illegal QM-side &
         &polarisability', alpha1)
    call utils_assert(alpha2 > 0D0, 'multipole_thole_T4: Illegal MM-side &
         &polarisability', alpha2)

    big_A = (alpha1*alpha2)**ONE_SIXTH

    if(r_norm > tiny_r) then
       lambda5 = multipole_thole_damping(r_norm, small_a, big_A, 2)
       lambda7 = multipole_thole_damping(r_norm, small_a, big_A, 3)
       lambda9 = multipole_thole_damping(r_norm, small_a, big_A, 4)
       rm5 = r_norm ** (-5.0_DP)
       rm7 = r_norm ** (-7.0_DP)
       rm9 = r_norm ** (-9.0_DP)
       l5 = 3.0_DP * lambda5 * rm5
       l7 = 15.0_DP * lambda7 * rm7
       l9 = 105.0_DP * lambda9 * rm9
       x = rvec%x
       y = rvec%y
       z = rvec%z

       ! --- block 1,1 ---
       multipole_thole_T4(1,1,1,1) = &
            3.0_DP*l5 - 6.0_DP*l7*x**2.0_DP + l9*x**4.0_DP
       multipole_thole_T4(2,1,1,1) = -3.0_DP*l7*x*y + l9*x**3.0_DP*y
       multipole_thole_T4(3,1,1,1) = -3.0_DP*l7*x*z + l9*x**3.0_DP*z
       multipole_thole_T4(1,2,1,1) = multipole_thole_T4(2,1,1,1)
       multipole_thole_T4(2,2,1,1) = &
            l5 + l9*x**2.0_DP*y**2.0_DP - l7*(x**2.0_DP + y**2.0_DP)
       multipole_thole_T4(3,2,1,1) = -l7*y*z + l9*x**2.0_DP*y*z
       multipole_thole_T4(1,3,1,1) = multipole_thole_T4(3,1,1,1)
       multipole_thole_T4(2,3,1,1) = multipole_thole_T4(3,2,1,1)
       multipole_thole_T4(3,3,1,1) = &
            l5 + l9*x**2.0_DP*z**2.0_DP - l7*(x**2.0_DP + z**2.0_DP)

       ! --- block 2,1 ---
       multipole_thole_T4(1,1,2,1) = -3.0_DP*l7*x*y + l9*x**3.0_DP*y
       multipole_thole_T4(2,1,2,1) = &
            l5 + l9*x**2.0_DP*y**2.0_DP - l7*(x**2.0_DP + y**2.0_DP)
       multipole_thole_T4(3,1,2,1) = -l7*y*z + l9*x**2.0_DP*y*z
       multipole_thole_T4(1,2,2,1) = multipole_thole_T4(2,1,2,1)
       multipole_thole_T4(2,2,2,1) = -3.0_DP*l7*x*y + l9*x*y**3.0_DP
       multipole_thole_T4(3,2,2,1) = -l7*x*z + l9*x*y**2.0_DP*z
       multipole_thole_T4(1,3,2,1) = multipole_thole_T4(3,1,2,1)
       multipole_thole_T4(2,3,2,1) = multipole_thole_T4(3,2,2,1)
       multipole_thole_T4(3,3,2,1) = -l7*x*y + l9*x*y*z**2.0_DP

       ! --- block 3,1 ---
       multipole_thole_T4(1,1,3,1) = -3.0_DP*l7*x*z + l9*x**3.0_DP*z
       multipole_thole_T4(2,1,3,1) = -l7*y*z + l9*x**2.0_DP*y*z
       multipole_thole_T4(3,1,3,1) = &
            l5 + l9*x**2.0_DP*z**2.0_DP - l7*(x**2.0_DP + z**2.0_DP)
       multipole_thole_T4(1,2,3,1) = multipole_thole_T4(2,1,3,1)
       multipole_thole_T4(2,2,3,1) = -l7*x*z + l9*x*y**2.0_DP*z
       multipole_thole_T4(3,2,3,1) = -l7*x*y + l9*x*y*z**2.0_DP
       multipole_thole_T4(1,3,3,1) = multipole_thole_T4(3,1,3,1)
       multipole_thole_T4(2,3,3,1) = multipole_thole_T4(3,2,3,1)
       multipole_thole_T4(3,3,3,1) = -3.0_DP*l7*x*z + l9*x*z**3.0_DP

       ! --- block 1,2 ---
       multipole_thole_T4(1,1,1,2) = multipole_thole_T4(1,1,2,1)
       multipole_thole_T4(2,1,1,2) = multipole_thole_T4(2,1,2,1)
       multipole_thole_T4(3,1,1,2) = multipole_thole_T4(3,1,2,1)
       multipole_thole_T4(1,2,1,2) = multipole_thole_T4(1,2,2,1)
       multipole_thole_T4(2,2,1,2) = multipole_thole_T4(2,2,2,1)
       multipole_thole_T4(3,2,1,2) = multipole_thole_T4(3,2,2,1)
       multipole_thole_T4(1,3,1,2) = multipole_thole_T4(1,3,2,1)
       multipole_thole_T4(2,3,1,2) = multipole_thole_T4(2,3,2,1)
       multipole_thole_T4(3,3,1,2) = multipole_thole_T4(3,3,2,1)

       ! --- block 2,2 ---
       multipole_thole_T4(1,1,2,2) = &
            l5 + l9*x**2.0_DP*y**2.0_DP - l7*(x**2.0_DP + y**2.0_DP)
       multipole_thole_T4(2,1,2,2) = -3.0_DP*l7*x*y + l9*x*y**3.0_DP
       multipole_thole_T4(3,1,2,2) = -l7*x*z + l9*x*y**2.0_DP*z
       multipole_thole_T4(1,2,2,2) = multipole_thole_T4(2,1,2,2)
       multipole_thole_T4(2,2,2,2) = &
            3.0_DP*l5 - 6.0_DP*l7*y**2.0_DP + l9*y**4.0_DP
       multipole_thole_T4(3,2,2,2) = -3.0_DP*l7*y*z + l9*y**3.0_DP*z
       multipole_thole_T4(1,3,2,2) = multipole_thole_T4(3,1,2,2)
       multipole_thole_T4(2,3,2,2) = multipole_thole_T4(3,2,2,2)
       multipole_thole_T4(3,3,2,2) = &
            l5 + l9*y**2.0_DP*z**2.0_DP - l7*(y**2.0_DP + z**2.0_DP)

       ! --- block 3,2 ---
       multipole_thole_T4(1,1,3,2) = -l7*y*z + l9*x**2.0_DP*y*z
       multipole_thole_T4(2,1,3,2) = -l7*x*z + l9*x*y**2.0_DP*z
       multipole_thole_T4(3,1,3,2) = -l7*x*y + l9*x*y*z**2.0_DP
       multipole_thole_T4(1,2,3,2) = multipole_thole_T4(2,1,3,2)
       multipole_thole_T4(2,2,3,2) = -3.0_DP*l7*y*z + l9*y**3.0_DP*z
       multipole_thole_T4(3,2,3,2) = &
            l5 + l9*y**2.0_DP*z**2.0_DP - l7*(y**2.0_DP + z**2.0_DP)
       multipole_thole_T4(1,3,3,2) = multipole_thole_T4(3,1,3,2)
       multipole_thole_T4(2,3,3,2) = multipole_thole_T4(3,2,3,2)
       multipole_thole_T4(3,3,3,2) = -3.0_DP*l7*y*z + l9*y*z**3.0_DP

       ! --- block 1,3 ---
       multipole_thole_T4(1,1,1,3) = multipole_thole_T4(1,1,3,1)
       multipole_thole_T4(2,1,1,3) = multipole_thole_T4(2,1,3,1)
       multipole_thole_T4(3,1,1,3) = multipole_thole_T4(3,1,3,1)
       multipole_thole_T4(1,2,1,3) = multipole_thole_T4(1,2,3,1)
       multipole_thole_T4(2,2,1,3) = multipole_thole_T4(2,2,3,1)
       multipole_thole_T4(3,2,1,3) = multipole_thole_T4(3,2,3,1)
       multipole_thole_T4(1,3,1,3) = multipole_thole_T4(1,3,3,1)
       multipole_thole_T4(2,3,1,3) = multipole_thole_T4(2,3,3,1)
       multipole_thole_T4(3,3,1,3) = multipole_thole_T4(3,3,3,1)

       ! --- block 2,3 ---
       multipole_thole_T4(1,1,2,3) = multipole_thole_T4(1,1,3,2)
       multipole_thole_T4(2,1,2,3) = multipole_thole_T4(2,1,3,2)
       multipole_thole_T4(3,1,2,3) = multipole_thole_T4(3,1,3,2)
       multipole_thole_T4(1,2,2,3) = multipole_thole_T4(1,2,3,2)
       multipole_thole_T4(2,2,2,3) = multipole_thole_T4(2,2,3,2)
       multipole_thole_T4(3,2,2,3) = multipole_thole_T4(3,2,3,2)
       multipole_thole_T4(1,3,2,3) = multipole_thole_T4(1,3,3,2)
       multipole_thole_T4(2,3,2,3) = multipole_thole_T4(2,3,3,2)
       multipole_thole_T4(3,3,2,3) = multipole_thole_T4(3,3,3,2)

       ! --- block 3,3 ---
       multipole_thole_T4(1,1,3,3) = &
            l5 + l9*x**2.0_DP*z**2.0_DP - l7*(x**2.0_DP + z**2.0_DP)
       multipole_thole_T4(2,1,3,3) = -l7*x*y + l9*x*y*z**2.0_DP
       multipole_thole_T4(3,1,3,3) = -3.0_DP*l7*x*z + l9*x*z**3.0_DP
       multipole_thole_T4(1,2,3,3) = multipole_thole_T4(2,1,3,3)
       multipole_thole_T4(2,2,3,3) = &
            l5 + l9*y**2.0_DP*z**2.0_DP - l7*(y**2.0_DP + z**2.0_DP)
       multipole_thole_T4(3,2,3,3) = -3.0_DP*l7*y*z + l9*y*z**3.0_DP
       multipole_thole_T4(1,3,3,3) = multipole_thole_T4(3,1,3,3)
       multipole_thole_T4(2,3,3,3) = multipole_thole_T4(3,2,3,3)
       multipole_thole_T4(3,3,3,3) = &
            3.0_DP*l5 - 6.0_DP*l7*z**2.0_DP + l9*z**4.0_DP

    else
       ! jd: Limit case
       multipole_thole_T4 = rank4zero
    end if

  end function multipole_thole_T4

end module multipole_ops
