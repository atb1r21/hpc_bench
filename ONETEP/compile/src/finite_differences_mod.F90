!=============================================================================!
!                        FINITE DIFFERENCES MODULE                            !
!=============================================================================!
! This module defines the subroutines for finite-difference calculations of   !
! the gradient, squared-modulus-of-gradient and the laplacian for a quantity  !
! on a GRID_INFO grid.                                                        !
!                                                                             !
! Finite differences of orders {2,4,6,8,10,12} can be used. Where 2nd order   !
! formulas are unavailable, 12th order is used instead. Slab-parallelism is   !
! automatically taken care of by internally using halos. Grid boundaries are  !
! taken care of by using progressively more forward- or backward-FD's.        !
!-----------------------------------------------------------------------------!
! Written by Jacek Dziedzic in 08/2010.
!-----------------------------------------------------------------------------!

module finite_differences

  use constants, only: DP

  implicit none
  private

  ! --------------------------
  ! --- Public subroutines ---
  ! --------------------------

  public :: finite_difference_initialise
  public :: finite_difference_set_geometry
  public :: finite_difference_gradient
  public :: finite_difference_mod_grad_sqd
  public :: finite_difference_laplacian

  ! ---------------------------------------------------------------------------
  ! ------------------- P U B L I C   V A R I A B L E S  ----------------------
  ! ---------------------------------------------------------------------------

  ! Data structure defining a 3D grid, parallelised over '12' slabs in real
  ! space. This is similar to the usual GRID_INFO, but...
  !  a) it is always orthorhombic,
  !  b) it stores the size of the grid's relevant portion (due to the
  !     granularity). All FD functions operate only on the relevant portion.
  !  c) only the attributes relevant for FDs are stored (no comms group data,
  !     etc).
  !
  ! An instance of this structure is returned by finite_difference_set_geometry,
  ! and can then be passed to other FD routines.
  type FD_GRID_INFO
     ! jd: Shortand for grid%n{1,2,3}
     integer :: pt1f, pt2f, pt3f
     ! jd: pt{1,2,3}f rounded down to the nearest multigrid magic number
     integer :: pq1f, pq2f, pq3f
     ! jd: Shortand for grid%ld{1,2,3}
     integer :: ld1f, ld2f, ld3f
     ! jd: Shortand for grid%{num,max}_slabs12
     integer :: num_slabs12, max_slabs12
     ! jd: Shortand for grid%{num,max}_slabs23
     integer :: num_slabs23, max_slabs23
     ! jd: Equivalent to num_slabs12 on all procs but last
     !     On the last proc, it's less by the required multigrid padding, pt3f-pq3f
     integer :: num_slabs12_pq
     ! jd: Shorthand for |grid%da{1,2,3}|
     real(kind=DP) :: d1f, d2f, d3f
     ! jd: Grid weight
     real(kind=DP) :: weight
     ! jd: Grid granularity (1, if multigrid not used)
     integer :: granularity
     ! jd: Number of first slab on this proc
     integer :: my_first_slab12
  end type FD_GRID_INFO

  public :: FD_GRID_INFO

  ! ---------------------------------------------------------------------------
  ! ------------------------------ P R I V A T E ------------------------------
  ! ---------------------------------------------------------------------------

  ! jd: Initialisation flags
  logical :: finite_difference_initialised = .false.

  ! jd: Halos for finite differences in parallel mode
  real(kind=DP), allocatable :: halo_top(:,:,:)
  real(kind=DP), allocatable :: halo_bot(:,:,:)

  ! jd: Maximum finite difference order that is supported
  integer, parameter :: max_order = 12
  ! jd: Half of the above
  integer, parameter :: max_order_half = max_order/2
  ! jd: Number of coefficients for every formula in one order
  integer, parameter :: ndata1 = (2*max_order+1)
  ! jd: Number of coefficients for all formulas in one order
  integer, parameter :: ndata = (max_order+1)*ndata1

  ! jd: Array of coefficients for finite difference derivative formulae
  !     4th index picks the derivative (1: 1st, 2: 2nd)
  !     3rd index picks the order (one of 4, 6, 8, 10, 12; other values are
  !     not supported but are guaranteed to contain 0.0 everywhere)
  !     2nd index picks the formula (0: central, +ve: progressively more forward
  !     formulae, -ve: progressively more backward formulae), only values from
  !     [-order/2,order/2] are ok, other entries will contain 0.0 everywhere
  !     1st index picks the point from x_{i-12} to x_{i+12}.
  real(kind=DP), dimension(-max_order:max_order, &
       -max_order_half:max_order_half, max_order, 2) :: coeff

  ! jd: To generate a table of these coefficients, run the following command
  !     in Mathematica:
  !     --- cut here ---
  !     minorder = 2; maxorder = 12; Table[Table[{If[m == 1, "1st", "2nd"] <>
  !     " derivative, order " <> IntegerString[j] <> ":", Table[Simplify[
  !     D[InterpolatingPolynomial[Table[{Subscript[x, i] + k h,
  !     f[Subscript[x, i + k]]}, {k, -j + p, p}], z], {z, m}] /. z ->
  !     Subscript[x, i]], {p, 0, j}] // MatrixForm}, {m, 1, 2}], {j, minorder,
  !     maxorder, 2}] // MatrixForm
  !     --- cut here ---

  ! ##########################################################################
  ! # 1st derivatives                                                        #
  ! ##########################################################################

  ! -- Order 12 --------------------------------------------------------------
  DATA coeff(:,-6,12,1) /&
       +   2310D0, -  30240D0, + 182952D0, - 677600D0, +1715175D0, -3136320D0, &
       +4268880D0, -4390848D0, +3430350D0, -2032800D0, + 914760D0, - 332640D0, &
       +  86021D0,                                                             &
              0D0,        0D0,        0D0,        0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0,        0D0  /

  DATA coeff(:,-5,12,1) /&
              0D0, -    210D0, +   2772D0, -  16940D0, +  63525D0, - 163350D0, &
       + 304920D0, - 426888D0, + 457380D0, - 381150D0, + 254100D0, - 152460D0, &
       +  55991D0,                                                             &
       +   2310D0,        0D0,        0D0,        0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0,        0D0  /

  DATA coeff(:,-4,12,1) /&
              0D0,        0D0, +     42D0, -    560D0, +   3465D0, -  13200D0, &
       +  34650D0, -  66528D0, +  97020D0, - 110880D0, + 103950D0, -  92400D0, &
       +  39611D0,                                                             &
       +   5040D0, -    210D0,        0D0,        0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0,        0D0  /

  DATA coeff(:,-3,12,1) /&
              0D0,        0D0,        0D0, -     14D0, +    189D0, -   1188D0, &
       +   4620D0, -  12474D0, +  24948D0, -  38808D0, +  49896D0, -  62370D0, &
       +  27599D0,                                                             &
       +   8316D0, -    756D0, +     42D0,        0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0,        0D0  /

  DATA coeff(:,-2,12,1) /&
              0D0,        0D0,        0D0,        0D0, +      7D0, -     96D0, &
       +    616D0, -   2464D0, +   6930D0, -  14784D0, +  25872D0, -  44352D0, &
       +  17589D0,                                                             &
       +  12320D0, -   1848D0, +    224D0, -     14D0, +      0D0, +      0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0,        0D0  /

  DATA coeff(:,-1,12,1) /&
              0D0,        0D0,        0D0,        0D0,        0D0, -      5D0, &
       +     70D0, -    462D0, +   1925D0, -   5775D0, +  13860D0, -  32340D0, &
       +   8580D0,                                                             &
       +  17325D0, -   3850D0, +    770D0, -    105D0, +      7D0, +      0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0,        0D0  /

  DATA coeff(:,0,12,1) /&
              0D0,        0D0,        0D0,        0D0,        0D0,        0D0, &
       +      5D0, -     72D0, +    495D0, -   2200D0, +   7425D0, -  23760D0, &
              0D0,                                                             &
       +  23760D0, -   7425D0, +   2200D0, -    495D0, +     72D0, -      5D0, &
              0D0,        0D0,        0D0,        0D0,        0D0,        0D0  /

  ! coeff(:,1:6,12,1) will be initialised in finite_difference_initialise

  ! -- Order 10 --------------------------------------------------------------
  DATA coeff(:,-6,10,1) / ndata1*0D0 /
  DATA coeff(:, 6,10,1) / ndata1*0D0 /

  DATA coeff(:,-5,10,1) /                                  0D0, 0D0, &
       +    252D0, -   2800D0, +  14175D0, -  43200D0, +  88200D0, &
       - 127008D0, + 132300D0, - 100800D0, +  56700D0, -  25200D0, &
       +   7381D0,                                                 &
              0D0,        0D0,        0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0, &
                                                         0D0, 0D0  /

  DATA coeff(:,-4,10,1) /                                  0D0, 0D0, &
              0D0, -     28D0, +    315D0, -   1620D0, +   5040D0, &
       -  10584D0, +  15876D0, -  17640D0, +  15120D0, -  11340D0, &
       +   4609D0,                                                 &
       +    252D0,        0D0,        0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0, &
                                                         0D0, 0D0  /

  DATA coeff(:,-3,10,1) /                                  0D0, 0D0, &
              0D0,        0D0, +      7D0, -     80D0, +    420D0, &
       -   1344D0, +   2940D0, -   4704D0, +   5880D0, -   6720D0, &
       +   3069D0,                                                 &
       +    560D0, -     28D0,        0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0, &
                                                         0D0, 0D0  /

  DATA coeff(:,-2,10,1) /                                  0D0, 0D0, &
              0D0,        0D0,        0D0, -      3D0, +     35D0, &
       -    189D0, +    630D0, -   1470D0, +   2646D0, -   4410D0, &
       +   1914D0,                                                 &
       +    945D0, -    105D0, +      7D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0, &
                                                         0D0, 0D0  /

  DATA coeff(:,-1,10,1) /                                  0D0, 0D0, &
              0D0,        0D0,        0D0,        0D0, +      2D0, &
       -     24D0, +    135D0, -    480D0, +   1260D0, -   3024D0, &
       +    924D0,                                                 &
       +   1440D0, -    270D0, +     40D0, -      3D0,        0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0, &
                                                         0D0, 0D0  /

  DATA coeff(:,0,10,1) /                                   0D0, 0D0, &
              0D0,        0D0,        0D0,        0D0,        0D0, &
       -      2D0, +     25D0, -    150D0, +    600D0, -   2100D0, &
              0D0,                                                 &
       +   2100D0, -    600D0, +    150D0, -     25D0, +      2D0, &
              0D0,        0D0,        0D0,        0D0,        0D0, &
                                                         0D0, 0D0  /

  ! coeff(:,1:5,10,1) will be initialised in finite_difference_initialise

  ! -- Order 8 ---------------------------------------------------------------
  DATA coeff(:,-6,8,1) / ndata1*0D0 /
  DATA coeff(:,-5,8,1) / ndata1*0D0 /
  DATA coeff(:, 5,8,1) / ndata1*0D0 /
  DATA coeff(:, 6,8,1) / ndata1*0D0 /

  DATA coeff(:,-4,8,1) /             0D0, 0D0, 0D0, 0D0, &
       +    105D0, -    960D0, +   3920D0, -   9408D0, &
       +  14700D0, -  15680D0, +  11760D0, -   6720D0, &
       +   2283D0,                                     &
              0D0,        0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0, &
                                   0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,-3,8,1) /             0D0, 0D0, 0D0, 0D0, &
              0D0, -     15D0, +    140D0, -    588D0, &
       +   1470D0, -   2450D0, +   2940D0, -   2940D0, &
       +   1338D0,                                     &
       +    105D0,        0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0, &
                                   0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,-2,8,1) /             0D0, 0D0, 0D0, 0D0, &
              0D0,        0D0, +      5D0, -     48D0, &
       +    210D0, -    560D0, +   1050D0, -   1680D0, &
       +    798D0,                                     &
       +    240D0, -     15D0,        0D0,        0D0, &
              0D0,        0D0,        0D0,        0D0, &
                                   0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,-1,8,1) /             0D0, 0D0, 0D0, 0D0, &
              0D0,        0D0,        0D0, -      3D0, &
       +     30D0, -    140D0, +    420D0, -   1050D0, &
       +    378D0,                                     &
       +    420D0, -     60D0, +      5D0,        0D0, &
              0D0,        0D0,        0D0,        0D0, &
                                   0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,0,8,1) /              0D0, 0D0, 0D0, 0D0, &
              0D0,        0D0,        0D0,        0D0, &
       +      3D0, -     32D0, +    168D0, -    672D0, &
              0D0,                                     &
       +    672D0, -    168D0, +     32D0, -      3D0, &
              0D0,        0D0,        0D0,        0D0, &
                                   0D0, 0D0, 0D0, 0D0  /
  ! coeff(:,1:4,8,1) will be initialised in finite_difference_initialise

  ! -- Order 6 ---------------------------------------------------------------
  DATA coeff(:,-6,6,1) / ndata1*0D0 /
  DATA coeff(:,-5,6,1) / ndata1*0D0 /
  DATA coeff(:,-4,6,1) / ndata1*0D0 /
  DATA coeff(:, 4,6,1) / ndata1*0D0 /
  DATA coeff(:, 5,6,1) / ndata1*0D0 /
  DATA coeff(:, 6,6,1) / ndata1*0D0 /

  DATA coeff(:,-3,6,1) /   0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
       +     10D0, -     72D0, +    225D0, &
       -    400D0, +    450D0, -    360D0, &
       +    147D0,                         &
              0D0,        0D0,        0D0, &
              0D0,        0D0,        0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,-2,6,1) /   0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
              0D0, -      2D0, +     15D0, &
       -     50D0, +    100D0, -    150D0, &
       +     77D0,                         &
       +     10D0,        0D0,        0D0, &
              0D0,        0D0,        0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,-1,6,1) /   0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
              0D0,        0D0, +      1D0, &
       -      8D0, +     30D0, -     80D0, &
       +     35D0,                         &
       +     24D0, -      2D0,        0D0, &
              0D0,        0D0,        0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,0,6,1) /   0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
              0D0,        0D0,        0D0, &
       -      1D0, +      9D0, -     45D0, &
              0D0,                         &
       +     45D0, -      9D0, +      1D0, &
              0D0,        0D0,        0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0  /
  ! coeff(:,1:3,6,1) will be initialised in finite_difference_initialise

  ! -- Order 4 ---------------------------------------------------------------
  DATA coeff(:,-6,4,1) / ndata1*0D0 /
  DATA coeff(:,-5,4,1) / ndata1*0D0 /
  DATA coeff(:,-4,4,1) / ndata1*0D0 /
  DATA coeff(:,-3,4,1) / ndata1*0D0 /
  DATA coeff(:, 3,4,1) / ndata1*0D0 /
  DATA coeff(:, 4,4,1) / ndata1*0D0 /
  DATA coeff(:, 5,4,1) / ndata1*0D0 /
  DATA coeff(:, 6,4,1) / ndata1*0D0 /

  DATA coeff(:,-2,4,1) / 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
       +      3D0, -     16D0, &
       +     36D0, -     48D0, &
       +     25D0,             &
              0D0,        0D0, &
              0D0,        0D0, &
                       0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,-1,4,1) / 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
              0D0, -      1D0, &
       +      6D0, -     18D0, &
       +     10D0,             &
       +      3D0,        0D0, &
              0D0,        0D0, &
                       0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,0,4,1) /  0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
              0D0,        0D0, &
       +      1D0, -      8D0, &
              0D0,             &
       +      8D0, -      1D0, &
              0D0,        0D0, &
                       0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  /
  ! coeff(:,1:2,4,1) will be initialised in finite_difference_initialise

  ! -- Order 2 ---------------------------------------------------------------
  DATA coeff(:,-6,2,1) / ndata1*0D0 /
  DATA coeff(:,-5,2,1) / ndata1*0D0 /
  DATA coeff(:,-4,2,1) / ndata1*0D0 /
  DATA coeff(:,-3,2,1) / ndata1*0D0 /
  DATA coeff(:,-2,2,1) / ndata1*0D0 /
  DATA coeff(:, 2,2,1) / ndata1*0D0 /
  DATA coeff(:, 3,2,1) / ndata1*0D0 /
  DATA coeff(:, 4,2,1) / ndata1*0D0 /
  DATA coeff(:, 5,2,1) / ndata1*0D0 /
  DATA coeff(:, 6,2,1) / ndata1*0D0 /

  DATA coeff(:,-1,2,1) / 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
              0D0,        0D0, &
       +      1D0, -      4D0, &
       +      3D0,             &
              0D0,        0D0, &
              0D0,        0D0, &
                       0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  /

  DATA coeff(:,0,2,1) /  0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
              0D0,        0D0, &
              0D0, -      1D0, &
              0D0,             &
       +      1D0,        0D0, &
              0D0,        0D0, &
                       0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  /
  ! coeff(:,1,2,1) will be initialised in finite_difference_initialise

  ! -- Orders that we don't care about ---------------------------------------
  DATA coeff(:,:,1,1) / ndata*0D0/
  DATA coeff(:,:,3,1) / ndata*0D0/
  DATA coeff(:,:,5,1) / ndata*0D0/
  DATA coeff(:,:,7,1) / ndata*0D0/
  DATA coeff(:,:,9,1) / ndata*0D0/
  DATA coeff(:,:,11,1) / ndata*0D0/

  ! ##########################################################################
  ! # 2nd derivatives                                                        #
  ! ##########################################################################

  ! -- Order 12 --------------------------------------------------------------
  DATA coeff(:,-6,12,2) /&
       418555D0,-5465520D0,32966604D0,-121646800D0,306489150D0,-557076960D0, &
       752145240D0,-764853408D0,587250675D0,-337836400D0,142878780D0,&
       -41976720D0,6706804D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 /

  DATA coeff(:,-5,12,2) /&
       0D0,-24305D0,319314D0,-1940070D0,7222325D0,-18396675D0,33904860D0,&
       -46613028D0,48570390D0,-38569575D0,23172050D0,-9329430D0,1265589D0,&
       418555D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 /

  DATA coeff(:,-4,12,2) /&
       0D0,0D0,3349D0,-44280D0,271095D0,-1018600D0,2624325D0,-4905648D0,&
       6863010D0,-7289040D0,5793975D0,-2378200D0,-630201D0,734520D0,-24305D0, &
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 /

  DATA coeff(:,-3,12,2) /&
       0D0,0D0,0D0,-743D0,9873D0,-60786D0,229790D0,-595485D0,1116126D0,&
       -1542156D0,1483812D0,16335D0,-1588015D0,995742D0,-67842D0,3349D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 /

  DATA coeff(:,-2,12,2) /&
       0D0,0D0,0D0,0D0,214D0,-2832D0,17292D0,-64240D0,159885D0,-267168D0,&
       208824D0,972576D0,-2119260D0,1208240D0,-125796D0,13008D0,-743D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0 /

  DATA coeff(:,-1,12,2) /&
       0D0,0D0,0D0,0D0,0D0,-50D0,600D0,-3036D0,6875D0,8250D0,-158400D0,&
       1339800D0,-2394678D0,1361250D0,-187000D0,29700D0,-3525D0,214D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0 /

  DATA coeff(:,0,12,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,-50D0,864D0,-7425D0,44000D0,-222750D0,1425600D0,&
       -2480478D0,1425600D0,-222750D0,44000D0,-7425D0,864D0,-50D0,&
       0D0,0D0,0D0,0D0,0D0,0D0 /

  ! coeff(:,1:6,12,2) will be initialised in finite_difference_initialise

  ! -- Order 10 --------------------------------------------------------------
  DATA coeff(:,-6,10,2) / ndata1*0D0 /
  DATA coeff(:, 6,10,2) / ndata1*0D0 /

  DATA coeff(:,-5,10,2) /&
       0D0,0D0,14258D0,-157800D0,794925D0,-2407200D0,4872700D0,-6932016D0,&
       7088550D0,-5232800D0,2754450D0,-972200D0,177133D0,0D0,0D0,0D0,0D0,0D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0 /

  DATA coeff(:,-4,10,2) /&
       0D0,0D0,0D0, -962D0, 10735D0,-54630D0,167560D0,-344820D0,501354D0,&
       -527660D0,401880D0,-188010D0,20295D0,14258D0,0D0,0D0,0D0,0D0,0D0,0D0,&
       0D0,0D0,0D0,0D0,0D0 /

  DATA coeff(:,-3,10,2) /&
       0D0,0D0,0D0,0D0,153D0,-1720D0,8830D0,-27360D0,56910D0,-83216D0,84420D0,&
       -29280D0,-32615D0,24840D0,-962D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,-2,10,2) /&
       0D0,0D0,0D0,0D0,0D0,-37D0,415D0,-2115D0,6420D0,-12530D0,13734D0,21210D0,&
       -57860D0,33255D0,-2645D0,153D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,-1,10,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,8D0,-80D0,315D0,-320D0,-3360D0,38304D0,-70070D0,&
       39360D0,-4680D0,560D0,-37D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,0,10,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,8D0,-125D0,1000D0,-6000D0,42000D0,-73766D0,&
       42000D0,-6000D0,1000D0,-125D0,8D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  ! coeff(:,1:5,10,2) will be initialised in finite_difference_initialise

  ! -- Order 8 ---------------------------------------------------------------
  DATA coeff(:,-6,8,2) / ndata1*0D0 /
  DATA coeff(:,-5,8,2) / ndata1*0D0 /
  DATA coeff(:, 5,8,2) / ndata1*0D0 /
  DATA coeff(:, 6,8,2) / ndata1*0D0 /

  DATA coeff(:,-4,8,2) /&
       0D0,0D0,0D0,0D0,3267D0,-29664D0,120008D0,-284256D0,435330D0,-448672D0,&
       312984D0,-138528D0,29531D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,&
       0D0/

  DATA coeff(:,-3,8,2) /&
       0D0,0D0,0D0,0D0,0D0,-261D0,2396D0,-9828D0,23688D0,-37030D0,38556D0,&
       -20916D0,128D0,3267D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,-2,8,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,47D0,-432D0,1764D0,-4144D0,5670D0,1008D0,&
       -9268D0,5616D0,-261D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,-1,8,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,-9D0,72D0,-196D0,-252D0,6930D0,-13216D0,&
       7308D0,-684D0,47D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,0,8,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,-9D0,128D0,-1008D0,8064D0,-14350D0,&
       8064D0,-1008D0,128D0,-9D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  ! coeff(:,1:4,8,2) will be initialised in finite_difference_initialise

  ! -- Order 6 ---------------------------------------------------------------
  DATA coeff(:,-6,6,2) / ndata1*0D0 /
  DATA coeff(:,-5,6,2) / ndata1*0D0 /
  DATA coeff(:,-4,6,2) / ndata1*0D0 /
  DATA coeff(:, 4,6,2) / ndata1*0D0 /
  DATA coeff(:, 5,6,2) / ndata1*0D0 /
  DATA coeff(:, 6,6,2) / ndata1*0D0 /

  DATA coeff(:,-3,6,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,137D0,-972D0,2970D0,-5080D0,5265D0,-3132D0,&
       812D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,-2,6,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,-13D0,93D0,-285D0,470D0,-255D0,-147D0,&
       137D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,-1,6,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,2D0,-12D0,15D0,200D0,-420D0,228D0,-13D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,0,6,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,2D0,-27D0,270D0,-490D0,270D0,-27D0,&
       2D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  ! coeff(:,1:3,6,2) will be initialised in finite_difference_initialise

  ! -- Order 4 ---------------------------------------------------------------
  DATA coeff(:,-6,4,2) / ndata1*0D0 /
  DATA coeff(:,-5,4,2) / ndata1*0D0 /
  DATA coeff(:,-4,4,2) / ndata1*0D0 /
  DATA coeff(:,-3,4,2) / ndata1*0D0 /
  DATA coeff(:, 3,4,2) / ndata1*0D0 /
  DATA coeff(:, 4,4,2) / ndata1*0D0 /
  DATA coeff(:, 5,4,2) / ndata1*0D0 /
  DATA coeff(:, 6,4,2) / ndata1*0D0 /

  DATA coeff(:,-2,4,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,11D0,-56D0,114D0,-104D0,35D0,0D0,0D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,-1,4,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,-1D0,4D0,6D0,-20D0,11D0,0D0,0D0,0D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,0,4,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,-1D0,16D0,-30D0,16D0,-1D0,0D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  ! coeff(:,1:2,4,2) will be initialised in finite_difference_initialise

  ! -- Order 2 ---------------------------------------------------------------
  DATA coeff(:,-6,2,2) / ndata1*0D0 /
  DATA coeff(:,-5,2,2) / ndata1*0D0 /
  DATA coeff(:,-4,2,2) / ndata1*0D0 /
  DATA coeff(:,-3,2,2) / ndata1*0D0 /
  DATA coeff(:,-2,2,2) / ndata1*0D0 /
  DATA coeff(:, 2,2,2) / ndata1*0D0 /
  DATA coeff(:, 3,2,2) / ndata1*0D0 /
  DATA coeff(:, 4,2,2) / ndata1*0D0 /
  DATA coeff(:, 5,2,2) / ndata1*0D0 /
  DATA coeff(:, 6,2,2) / ndata1*0D0 /

  DATA coeff(:,-1,2,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,1D0,-2D0,1D0,0D0,0D0,0D0,0D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  DATA coeff(:,0,2,2) /&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,1D0,-2D0,1D0,0D0,0D0,&
       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0/

  ! coeff(:,1,2,2) will be initialised in finite_difference_initialise

  ! -- Orders that we don't care about ---------------------------------------
  DATA coeff(:,:,1,2) / ndata*0D0/
  DATA coeff(:,:,3,2) / ndata*0D0/
  DATA coeff(:,:,5,2) / ndata*0D0/
  DATA coeff(:,:,7,2) / ndata*0D0/
  DATA coeff(:,:,9,2) / ndata*0D0/
  DATA coeff(:,:,11,2) / ndata*0D0/

  ! ##########################################################################

  ! jd: Norms for the finite differences of orders 1..12 for 1st and 2nd deriv.
  real(kind=DP), dimension(max_order,2), parameter :: norms= &
       reshape( source = &
       (/0D0, 1D0/2D0, 0D0, 1D0/12D0, 0D0, 1D0/60D0, 0D0, 1D0/840D0, &
       0D0, 1D0/2520D0, 0D0, 1D0/27720D0, &
       0D0, 1D0, 0D0, 1D0/12D0,0D0,1D0/180D0,0D0,1D0/5040D0,0D0,1D0/25200D0,&
       0D0,1D0/831600D0 /), shape = (/max_order, 2/) )

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine finite_difference_initialise()
    !=========================================================================!
    ! Initialises the finite differences by completing the computation of the !
    ! finite difference coefficients that have been half-prepared by the DATA !
    ! statements.                                                             !
    ! There is no need to call this subroutine explicitly -- it will be called!
    ! upon first use of any computational FD subroutine. You can call it      !
    ! explicitly if you like, and any number of times.                        !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 08/2010                                      !
    !=========================================================================!

    use utils, only: utils_trace_in, utils_trace_out, utils_abort, &
         utils_sanity_check

    implicit none

    ! jd: Internal variables
    integer :: i, j, k, l ! jd: Indices
    character(len=256) :: problem_at

    !------------------------------------------------------------------------

    if(finite_difference_initialised) return

    call utils_trace_in('finite_difference_initialise')

    ! jd: Exploit the antysymmetry of the formulae for 1st derivative and the
    !     symmetry of the formulae for the 2nd derivative to fill the +ve half
    do l=1, 2
       do k=1, max_order
          do j=1, max_order_half
             do i=-max_order, max_order
                coeff(i,j,k,l) = (-1)**l * coeff(-i,-j,k,l)
             end do
          end do
       end do
    end do

    ! jd: Perform a sanity-check -- all the coefficients for any formula (j) for
    !     any order (k) for any derivative (l) should add up to zero
    do l=1, 2
       do k=1, max_order
          do j=-max_order_half, max_order_half
             if(abs(sum(coeff(:,j,k,l)))>1D-12) then
                write(problem_at,'(a,i0,a,i0,a,i0,a)') &
                     'coeff(:,',j,',',k,',',l,')'
                call utils_abort('Inconsistency in the FD coefficient array &
                     &detected in finite_differences_mod, check '//&
                     trim(problem_at))
             end if
          end do
       end do
    end do

    ! jd: Complete the calculation by multiplying the coefficients by the norm
    do l=1, 2
       do k=1, max_order
          coeff(:,:,k,l) = coeff(:,:,k,l) * norms(k,l)
       end do
    end do

    ! jd: This will spot any uninitialised elements if you can force the
    !     compiler to default-initialise everything with a NaN, e.g.
    !     by passing -finit-real=nan to gfortran
    call utils_sanity_check(coeff(:,:,:,1),'FD coeff., 1st derivative')
    call utils_sanity_check(coeff(:,:,:,2),'FD coeff., 2nd derivative')

    finite_difference_initialised = .true.

    call utils_trace_out('finite_difference_initialise')

  end subroutine finite_difference_initialise
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  type(FD_GRID_INFO) function finite_difference_set_geometry(grid, granularity, is_periodic)
    !=========================================================================!
    ! Sets up the geometry of the underlying FD grid, as determined by the    !
    ! input arguments 'grid' and 'multigrid_granularity'.                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   grid (input)        : The grid defining the geometry to be set.       !
    !   granularity (input, optional) : The grid sizes will be rounded down to!
    !                                   an integer multiple of this value --  !
    !                                   this is useful for the multigrid,     !
    !                                   where there are constraints on the    !
    !                                   grid dimensions. If you do not intend !
    !                                   to use the multigrid solver, omit this!
    !                                   argument of pass 1.                   !
    !   is_periodic (input, optional) : Logical array indicating whether the  !
    !                                   grid is periodic (T) or not (F) along !
    !                                   each lattice vector. This is required !
    !                                   to determine whether grid dimensions  !
    !                                   should be rounded to even (periodic)  !
    !                                   or odd (non-periodic) numbers.        !
    !                                   If granularity is not present or 1,   !
    !                                   this can be omitted, otherwise it is  !
    !                                   required.                             !
    ! Return value:                                                           !
    !   An FD_GRID_INFO : This is similar to the usual GRID_INFO, but...      !
    !                     a) is always orthorhombic,                          !
    !                     b) stores the size of the grid's relevant portion   !
    !                        (due to the granularity).                        !
    !                     Store the obtained value and pass it as an argument !
    !                     to other FD functions.                              !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 08/2010                                      !
    ! Updated by James C. Womack, 2017, to add support for periodic BCs       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_round_for_multigrid
    use comms, only: pub_my_proc_id, pub_total_num_procs
    use geometry, only: magnitude
    use utils, only: utils_trace_in, utils_trace_out, utils_assert

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in)      :: grid
    integer, intent(in), optional    :: granularity
    logical, intent(in), optional    :: is_periodic(3)

    ! jd: Return value
    type(FD_GRID_INFO)               :: fd_grid

    ! jd: Local variables
    integer :: local_granularity
    logical :: local_is_periodic(3)
    character(len=*), parameter :: myself = 'finite_difference_set_geometry'
    integer :: i
    integer :: d(3)

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    if (present(granularity)) then
       local_granularity = granularity
    else
       local_granularity = 1
    end if

    if (present(is_periodic)) then
       local_is_periodic(:) = is_periodic(:)
    else
       ! JCW: If is_periodic is not present, local_granularity must be 1
       call utils_assert(local_granularity==1,"Error in "//myself//": &
            &When granularity /= 1, is_periodic must be present.")
       ! JCW: local_is_periodic is only used when local_granularity > 1, so
       ! JCW: this default value should never be
       ! JCW: used for anything.
       local_is_periodic(:) = .false.
    end if

    if (local_granularity /= 1) then
       ! JCW: Set the value of d to determine whether grid dimensions are rounded down
       ! JCW: to odd (non-periodic) or even (periodic) values
       do i = 1, 3
          if (local_is_periodic(i)) then
             d(i) = 0
          else
             d(i) = 1
          end if
       end do

    else
       ! JCW: If granularity == 1, then set d = 1 to prevent grid truncation
       d(1:3) = 1
    end if

    ! jd: Initialise the geometry of the FD grid
    fd_grid%granularity = local_granularity
    fd_grid%pt1f = grid%n1
    fd_grid%pt2f = grid%n2
    fd_grid%pt3f = grid%n3
    ! JCW: Only round down if either (or both) the grid dimension is not a multiple of
    ! JCW: granularity, or the grid along that direction is not periodic
    ! JCW: --> If the dimension is periodic, it is possible that the grid dimensions
    ! JCW:     have been previously determined to satisfy the size constraints of the
    ! JCW:     multigrid solver in cell_grid_distribute. In this case the FD grid can
    ! JCW:     have identical dimension (fd_grid%pq{1,2,3}f = fd_grid%pt{1,2,3}).
    if ((mod(fd_grid%pt1f,local_granularity) /= 0).or.(.not.local_is_periodic(1))) then
       fd_grid%pq1f = cell_grid_round_for_multigrid(fd_grid%pt1f,local_granularity,d(1),"DOWN")
    else
       fd_grid%pq1f = fd_grid%pt1f
    end if
    if ((mod(fd_grid%pt2f,local_granularity) /= 0).or.(.not.local_is_periodic(2))) then
       fd_grid%pq2f = cell_grid_round_for_multigrid(fd_grid%pt2f,local_granularity,d(2),"DOWN")
    else
       fd_grid%pq2f = fd_grid%pt2f
    end if
    if ((mod(fd_grid%pt3f,local_granularity) /= 0).or.(.not.local_is_periodic(3))) then
       fd_grid%pq3f = cell_grid_round_for_multigrid(fd_grid%pt3f,local_granularity,d(3),"DOWN")
    else
       fd_grid%pq3f = fd_grid%pt3f
    end if
    fd_grid%ld1f = grid%ld1
    fd_grid%ld2f = grid%ld2
    fd_grid%ld3f = grid%ld3
    fd_grid%d1f = magnitude(grid%da1)
    fd_grid%d2f = magnitude(grid%da2)
    fd_grid%d3f = magnitude(grid%da3)
    fd_grid%weight = fd_grid%d1f * fd_grid%d2f * fd_grid%d3f
    fd_grid%max_slabs12 = grid%max_slabs12
    fd_grid%num_slabs12 = grid%num_my_slabs12
    fd_grid%max_slabs23 = grid%max_slabs23
    fd_grid%num_slabs23 = grid%num_slabs23
    fd_grid%my_first_slab12 = grid%first_slab12(pub_my_proc_id)

    ! jd: Determine the num_slabs12 equivalent for the FD grid (pq)
    if(pub_my_proc_id < pub_total_num_procs-1) then
       fd_grid%num_slabs12_pq = fd_grid%num_slabs12
    else ! jd: it only differs from num_slabs12 on the last proc
       fd_grid%num_slabs12_pq = &
            fd_grid%num_slabs12 - (fd_grid%pt3f - fd_grid%pq3f)
       call utils_assert(fd_grid%num_slabs12_pq > 0, &
            trim(myself)//&
            ': Cannot distribute the requested FD grid across processors -- &
            &either (a) the no-man''s-land between the FD grid and the host &
            &grid is so thick, it extends to more than one processor; or (b)&
            &you have more MPI ranks than 12-slabs. Neither is allowed. It &
            &is likely (b). Either decrease the number of MPI ranks or reori&
            &ent your system so that the largest dimensions is along Z. &
            &If (a), try to make the FD grid closer in size to the host &
            &grid or, if this is not possible, decrease the number of &
            &processors.', fd_grid%num_slabs12_pq)
    end if

    ! jd: Sanity check -- FD grid cannot extend beyond the original grid
    call utils_assert(fd_grid%pq1f <= fd_grid%pt1f .and. &
         fd_grid%pq2f <= fd_grid%pt2f .and. fd_grid%pq3f <= fd_grid%pt3f, &
         'Inconsistency detected in '//trim(myself))

    call utils_trace_out(myself)

    finite_difference_set_geometry = fd_grid

    !------------------------------------------------------------------------

  end function finite_difference_set_geometry
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine finite_difference_derivative(fcn_der,fcn,n,d,order,fd_grid, &
       cell_grid_info, square_result)
    !=========================================================================!
    ! This subroutine takes a function 'fcn' on a real space grid and uses    !
    ! finite difference methods to find a first or second derivative of the   !
    ! function on the grid.                                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   fcn_der  (inout) the derivative of fcn found by finite differences.   !
    !   fcn      (in)    the scalar function whose gradient we're after.      !
    !   n        (in)    selects between 1st and 2nd derivative.              !
    !   d        (in)    Cartesian direction wrt which we differentiate       !
    !                    (1=x, 2=y, 3=z)!                                     !
    !   order    (in)    the finite difference order to use.                  !
    !   fd_grid  (in)    the grid on which to operate                         !
    !   cell_grid_info (in) a GRID_INFO structure holding info on the parallel!
    !                       distribution of slabs that fd_grid does not hold. !
    !   square_result (in, optional) If true, squares of derivatives will be  !
    !                                returned (added). Not very elegant, but  !
    !                                allows efficient approach for |grad f|^2.!
    !-------------------------------------------------------------------------!
    ! All the arrays are assumed to be in the distributed slab representation.!
    !                                                                         !
    ! Note that this subroutine operates on the passed FD grid (obtained from !
    ! finite_difference_set_geometry, which might be slightly smaller than    !
    ! the original ('host') grid. Values outside the FD grid will be untouched!
    ! in the output. Values close to the boundary of the FD grid will be      !
    ! computed with progressively less central and more  backward or forward  !
    ! differences, but always of the desired order. Details of parallel slab  !
    ! distribution do not affect results.                                     !
    !-------------------------------------------------------------------------!
    ! Notes, caveats:                                                         !
    ! - The value of the computed derivative is actually *added* to fcn_der,  !
    !   this simplifies the calculation of the laplacian and |grad|^2.        !
    !   Thus it is the CALLER'S RESPONSIBILITY TO INITIALIZE FCN_DER prior to !
    !   calling, even if to zeros.                                            !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal 30/9/2008                                      !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    ! Rewritten by Jacek Dziedzic, 08/08/2010.                                !
    ! Optimized by Jacek Dziedzic, 19/04/2011.                                !
    ! Generalized by Jacek Dziedzic, 23/04/2013 to support extra thin slabs.  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id
!$  use rundat, only: pub_threads_max
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, utils_sanity_check, &
         utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(FD_GRID_INFO), intent(in) :: fd_grid
    type(GRID_INFO), intent(in)    :: cell_grid_info
    real(kind=DP), dimension(fd_grid%ld1f,fd_grid%ld2f,fd_grid%max_slabs12), &
         intent(inout)             :: fcn_der
    real(kind=DP), dimension(fd_grid%ld1f,fd_grid%ld2f,fd_grid%max_slabs12), &
         intent(in)                :: fcn
    integer, intent(in)            :: order
    integer, intent(in)            :: n
    integer, intent(in)            :: d
    logical, intent(in), optional  :: square_result

    ! jd: Internal variables
    character(len=*), parameter :: myself = 'finite_difference_derivative'
    logical :: must_square_result  ! jd: .true. if square_result present and T
    integer :: i,j,k               ! jd: Indices
    real(kind=DP) :: ddf = 0.0_DP  ! jd: Temporary

    integer :: order_half          ! jd: order/2
    integer :: m, m_min, m_max, how_far_from_my_slice ! jd: Indices
    real(kind=DP) :: value_here      ! jd: Temporary
    real(kind=DP) :: der             ! jd: Accumulator for results
    integer :: formula               ! jd: Number of relevant FD formula
    integer :: i1, i2, i3            ! jd: Indices
    integer :: i3_global
    integer :: i1max, i2max, i3max   ! jd: Max values for indices
    integer :: i3min_with_halo, i3max_with_halo
    real(kind=DP) :: coeff_times_inv_norm(-max_order:max_order, &
         -max_order_half:max_order_half) ! jd: Look-up
    real(kind=DP) :: fcn_slice(1-max_order:max(fd_grid%pq1f, &
         fd_grid%pq2f,fd_grid%num_slabs12_pq)+max_order) ! jd: Lookup

    !------------------------------------------------------------------------

    call timer_clock(myself,1)
    call utils_trace_in(myself)

    call finite_difference_initialise

    must_square_result = .false.
    if(present(square_result)) then
       if(square_result) then
          must_square_result = .true.
       end if
    end if

    ! jd: Sanity check on the direction, n, order and input and grid size
    if(.not. (d>0 .and. d<=3)) &
         call utils_assert(.false., &
         'Unrecognized direction in '//myself,d)
    if(.not. (n==1 .or. n==2)) &
         call utils_assert(.false., &
         myself//': can only compute 1st and 2nd derivatives')
    call utils_assert(order > 0 .and. order <=12 .and. mod(order,2)==0, &
         'Unsupported FD order in '//myself,order)
    call utils_sanity_check(fcn,'fcn in '//myself)
    call utils_assert(fd_grid%pq1f >= order + 1, &
         'The X dimension of the fine grid is smaller than the FD stencil')
    call utils_assert(fd_grid%pq2f >= order + 1, &
         'The Y dimension of the fine grid is smaller than the FD stencil')
    call utils_assert(fd_grid%pq3f >= order + 1, &
         'The Z dimension of the fine grid is smaller than the FD stencil')

    order_half = order / 2

    ! jd: Prepare parallel halos for fcn
    call parallel_prepare_fd_halos(fcn,order,fd_grid,cell_grid_info)

    ! -------------------------------------------------------------
    ! --- Loop over the array and fill it by finite differences ---
    ! -------------------------------------------------------------

    ! The proposed strategy is a bit more complex than a straightforward scan
    ! through the array and filling it in. The problem of the straightforward
    ! approach is that for derivatives wrt y or z the array is scanned in a
    ! non-cache-efficient fashion. Ideally, we want the repeated reads from
    ! fcn to be performed along the minor direction. The writes to grad_fcn are
    ! less important, because there is only one write for the whole stencil
    ! (and there are 11 reads for an order-10 stencil). Transposing the whole
    ! fcn needs a lot of storage. I propose to transpose only the current row of
    ! fcn -- thus we perform only 1 non-cache-efficient set of reads from fcn,
    ! then subsequent 11 reads are cache-aligned. This is accomplished by
    ! copying to 'fcn_slice' the current row (or 'slice') of fcn. Additionally,
    ! this takes care of the parallel halos, so that we don't have to worry
    ! about them when doing the 11 reads later on. The complication is that
    ! the loop ordering changes depending on the derivative direction. There
    ! are generic indices: i1, i2, i3 which correspond to permutations of
    ! i, j and k. The slice is filled in the innermost loop (over i3).
    ! Timings for a sample system to roughly indicate performance gains:
    ! 70.3s - using derivative_at_point (simple, non-efficient)
    ! 54.4s - +manually inlining derivative_at_point into this subroutine
    ! 47.8s - +clever loop reordering
    ! 32.3s - +current approach with fcn_slice.

    ! jd: Determine the grid spacing and index bounds for this direction
    if(d==1) then
       ddf=fd_grid%d1f
       i1max = fd_grid%num_slabs12_pq  ! z
       i2max = fd_grid%pq2f            ! y
       i3max = fd_grid%pq1f            ! x
    else if(d==2) then
       ddf=fd_grid%d2f
       i1max = fd_grid%num_slabs12_pq  ! z
       i2max = fd_grid%pq1f            ! x
       i3max = fd_grid%pq2f            ! y
    else
       ddf=fd_grid%d3f
       i1max = fd_grid%pq2f            ! y
       i2max = fd_grid%pq1f            ! x
       i3max = fd_grid%num_slabs12_pq  ! z
    end if

    ! jd: Calculate the product of relevant coefficients and 1/ddf,
    !     store these in a lookup to save time in the inner loops
    coeff_times_inv_norm(:,:) = coeff(:,:,order,n) * 1.0_DP/(ddf**n)

    i=0; j=0; k=0 ! jd: Silence compiler warnings

    ! jd: The bounds of the slice are trivially 1..max for slices along X or Y
    if(d<3) then
       i3min_with_halo = 1
       i3max_with_halo = i3max
    ! jd: But can extend by as many as 'order' points for slices along Z, as
    !     these can extend into halos. We must have a slice big enough to fit
    !     the stencil in even the worst cases.
    else
       i3min_with_halo = 1-order
       i3max_with_halo = i3max+order
    end if

    ! jme: initialisation for OpenMP compilations
!$  fcn_slice = 0.0_DP

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) &
!$OMP PRIVATE(i1,i2,i3,i,j,k,value_here,i3_global,how_far_from_my_slice, &
!$OMP      formula,m_max,m_min,der,m) &
!$OMP FIRSTPRIVATE(fcn_slice) &
!$OMP SHARED(fd_grid,d,i1max,i2max,i3max,i3max_with_halo,i3min_with_halo,fcn, &
!$OMP      cell_grid_info,halo_top,halo_bot,order_half,must_square_result, &
!$OMP      order,coeff_times_inv_norm,pub_my_proc_id,fcn_der,pub_threads_max)
    do i1=1, i1max                             ! -- Loop over major index --
       ! jd: Determine index permutation
       if(d==1) then
          k=i1
       else if(d==2) then
          k=i1
       else
          j=i1
       end if

       do i2=1, i2max                   ! -- Loop over intermediate index --
          ! jd: Determine index permutation
          if(d==1) then
             j=i2
          else if(d==2) then
             i=i2
          else
             i=i2
          end if

          ! jd: First fill the 1D slice from fcn to fcn_slice.
          !     This is the only time we scan fcn non-cache-efficiently.
          !     We also take care of slab halos here.
          do i3=i3min_with_halo, i3max_with_halo ! -- Loop over minor index --
             ! jd: Determine index permutation
             if(d==1) then
                i=i3
             else if(d==2) then
                j=i3
             else
                k=i3
             end if

             ! jd: Copy the right point from fcn (or halo) to fcn_slice
             if(d==1 .or. d==2) then ! jd: Slice along x or y
                if(i>=1 .and. i<=fd_grid%pq1f .and. &
                     j>=1 .and. j<=fd_grid%pq2f) then
                   ! This is the most costly line, except for anything in
                   ! the m-loop below. No wonder, thats a lot of cache misses
                   value_here = fcn(i,j,k)
                else
                   ! This should never happen. Once we're sure it doesn't, the
                   ! if above should be eliminated for efficiency
                   call utils_abort('The impossible happened in '//myself)
                end if
             else ! Slice along z
                if(k>=1 .and. k<=fd_grid%num_slabs12_pq) then
                   value_here = fcn(i,j,k)          ! inside my slab of grid
                else                                ! outside my slab of grid
                   ! jd: If out of my slab, but in the FD grid, then
                   !     in someone else's slab. Look in the halos
                   !     - above my slab
                   i3_global = i3 + cell_grid_info%first_slab12(pub_my_proc_id)-1
                   if(k<1 .and. i3_global>=1) then ! not mine, but not outside the grid
                      how_far_from_my_slice = k-1
                      value_here = halo_top(i,j,how_far_from_my_slice)
                   !     - below my slab
                   else if(i3_global <= fd_grid%pq3f .and. &
                        k>fd_grid%num_slabs12_pq) then ! not mine, but not outside the grid
                      how_far_from_my_slice = k-fd_grid%num_slabs12_pq
                      value_here = halo_bot(i,j,how_far_from_my_slice)
                   !     - beyond the grid
                   else
                      value_here = 0D0
                   end if
                end if
             end if

             fcn_slice(i3) = value_here

          end do

          ! jd: The 1D slice is filled, just iterate over it, there are no
          !     more halo or out of grid issues.
          do i3=1, i3max
             if(d==1) then
                i=i3
             else if(d==2) then
                j=i3
             else
                k=i3
             end if

             ! jd: Pick the right FD formula, which is 0 (central diff)
             !     in the bulk and an appropriate mix of forward/backward
             !     when close to the boundaries
             formula = 0
             if(d==1 .or. d==2) then
                if(i3 <= order_half) formula = order_half-i3+1
                if(i3 > i3max-order_half) &
                     formula = i3max-i3-order_half
             else
                ! jd: Along z it's a bit more tricky -- with thin chunks per
                !     proc it might happen that we're not in the bulk on a
                !     proc that's not the first or last.
                i3_global = i3 + cell_grid_info%first_slab12(pub_my_proc_id)-1
                if(i3_global <= order_half) formula = order_half-i3_global+1
                if(i3_global > fd_grid%pq3f-order_half) &
                     formula = fd_grid%pq3f - i3_global - order_half
             end if

             m_min=formula-order_half
             m_max=m_min+order

             ! jd: Go over all the points in the stencil and sum up the
             !     contributions. Accumulate to a local (likely a register),
             !     only then spill to the array, for efficiency.
             der = 0D0
             do m=m_min, m_max
                der = der + coeff_times_inv_norm (m,formula) * &
                     fcn_slice(i3+m)
             end do

             if(must_square_result) der = der*der
             fcn_der(i,j,k) = fcn_der(i,j,k) + der

          end do ! over i3
       end do ! over i2
    end do ! over i1
!$OMP END PARALLEL DO

    ! jd: Clean up by destroying the halo arrays
    call parallel_destroy_fd_halos

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine finite_difference_derivative
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine finite_difference_gradient(grad_fcn,fcn,order,fd_grid, &
       cell_grid_info)
    !=========================================================================!
    ! This subroutine takes a function 'fcn' on the real space grid and uses  !
    ! finite difference methods to find the gradient of the function on the   !
    ! grid.  The gradient is returned in grad_fcn.                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   grad_fcn (out) the gradient of fcn found by finite differences.       !
    !   fcn      (in)  the scalar function whose gradient we're after.        !
    !   order    (in)  the finite difference order to use.                    !
    !   fd_grid  (in)  the grid on which to operate                           !
    !   cell_grid_info (in) a GRID_INFO structure holding info on the parallel!
    !                       distribution of slabs that fd_grid does not hold. !
    !-------------------------------------------------------------------------!
    ! All the arrays are assumed to be in the distributed slab representation.!
    !                                                                         !
    ! Note that this subroutine operates on the passed FD grid (obtained from !
    ! finite_difference_set_geometry, which might be slightly smaller than    !
    ! the original ('host') grid. Values outside the FD grid will be zeroed   !
    ! in the output. Values close to the boundary of the FD grid will be      !
    ! computed with progressively less central and more  backward or forward  !
    ! differences, but always of the desired order. Details of parallel slab  !
    ! distribution do not affect results.                                     !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal 30/9/2008                                      !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    ! Rewritten by Jacek Dziedzic, 08/08/2010.                                !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(FD_GRID_INFO), intent(in) :: fd_grid
    type(GRID_INFO), intent(in)    :: cell_grid_info
    real(kind=DP), dimension(fd_grid%ld1f,fd_grid%ld2f,fd_grid%max_slabs12,3), &
         intent(out)               :: grad_fcn
    real(kind=DP), dimension(fd_grid%ld1f,fd_grid%ld2f,fd_grid%max_slabs12), &
         intent(in)                :: fcn
    integer, intent(in)            :: order

    ! jd: Local variables
    integer :: d
    character(len=*), parameter :: myself = 'finite_difference_gradient'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    call finite_difference_initialise

    grad_fcn = 0D0
    do d=1, 3
       call finite_difference_derivative(grad_fcn(:,:,:,d), fcn, 1, d, order, &
            fd_grid, cell_grid_info)
    end do

    call utils_trace_out(myself)

  end subroutine finite_difference_gradient
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine finite_difference_laplacian(lap_fcn,fcn,order,fd_grid, &
       cell_grid_info)
    !=========================================================================!
    ! This subroutine takes a function 'fcn' on the real space grid and uses  !
    ! finite difference methods to find the laplacian of the function on the  !
    ! grid.  The laplacian is returned in lap_fcn.                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   lap_fcn (out)  the laplacian of fcn found by finite differences.      !
    !   fcn      (in)  the scalar function whose gradient we're after.        !
    !   order    (in)  the finite difference order to use.                    !
    !   fd_grid  (in)  the grid on which to operate                           !
    !   cell_grid_info (in) a GRID_INFO structure holding info on the parallel!
    !                       distribution of slabs that fd_grid does not hold. !
    !-------------------------------------------------------------------------!
    ! All the arrays are assumed to be in the distributed slab representation.!
    !                                                                         !
    ! Note that this subroutine operates on the passed FD grid (obtained from !
    ! finite_difference_set_geometry, which might be slightly smaller than    !
    ! the original ('host') grid. Values outside the FD grid will be zeroed   !
    ! in the output. Values close to the boundary of the FD grid will be      !
    ! computed with progressively less central and more  backward or forward  !
    ! differences, but always of the desired order. Details of parallel slab  !
    ! distribution do not affect results.                                     !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal 30/9/2008                                      !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    ! Rewritten by Jacek Dziedzic, 08/08/2010.                                !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(FD_GRID_INFO), intent(in) :: fd_grid
    type(GRID_INFO), intent(in)    :: cell_grid_info
    real(kind=DP), dimension(fd_grid%ld1f,fd_grid%ld2f,fd_grid%max_slabs12), &
         intent(out)               :: lap_fcn
    real(kind=DP), dimension(fd_grid%ld1f,fd_grid%ld2f,fd_grid%max_slabs12), &
         intent(in)                :: fcn
    integer, intent(in)            :: order

    ! jd: Local variables
    integer :: d
    character(len=*), parameter :: myself = 'finite_difference_laplacian'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    call finite_difference_initialise

    lap_fcn = 0D0
    do d=1, 3
       ! NB: finite_difference_derivative *adds* to lap_fcn
       call finite_difference_derivative(lap_fcn(:,:,:), fcn, 2, d, order, &
            fd_grid, cell_grid_info)
    end do

    call utils_trace_out(myself)

  end subroutine finite_difference_laplacian
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine finite_difference_mod_grad_sqd(mod_grad_sqd_fcn,fcn,order, &
       fd_grid,cell_grid_info)
    !=========================================================================!
    ! This subroutine takes a function 'fcn' on the real space grid and uses  !
    ! finite difference methods to find the |grad fcn|^2 of the function on   !
    ! the grid, which is returned in grad_mod_sqd_fcn.                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mod_grad_sqd_fcn (out)  the result found by finite differences.       !
    !   fcn      (in)  the scalar function whose gradient we're after.        !
    !   order    (in)  the finite difference order to use.                    !
    !   fd_grid  (in)  the grid on which to operate                           !
    !   cell_grid_info (in) a GRID_INFO structure holding info on the parallel!
    !                       distribution of slabs that fd_grid does not hold. !
    !-------------------------------------------------------------------------!
    ! All the arrays are assumed to be in the distributed slab representation.!
    !                                                                         !
    ! Note that this subroutine operates on the passed FD grid (obtained from !
    ! finite_difference_set_geometry, which might be slightly smaller than    !
    ! the original ('host') grid. Values outside the FD grid will be zeroed   !
    ! in the output. Values close to the boundary of the FD grid will be      !
    ! computed with progressively less central and more  backward or forward  !
    ! differences, but always of the desired order. Details of parallel slab  !
    ! distribution do not affect results.                                     !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal 30/9/2008                                      !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    ! Rewritten by Jacek Dziedzic, 08/08/2010.                                !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use utils, only: utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    type(FD_GRID_INFO), intent(in) :: fd_grid
    type(GRID_INFO), intent(in)    :: cell_grid_info
    real(kind=DP), dimension(fd_grid%ld1f,fd_grid%ld2f,fd_grid%max_slabs12), &
         intent(out)               :: mod_grad_sqd_fcn
    real(kind=DP), dimension(fd_grid%ld1f,fd_grid%ld2f,fd_grid%max_slabs12), &
         intent(in)                :: fcn
    integer, intent(in)            :: order

    ! jd: Local variables
    integer :: d
    character(len=*), parameter :: myself = 'finite_difference_mod_grad_sqd'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    call finite_difference_initialise

    mod_grad_sqd_fcn = 0D0
    do d=1, 3
       call finite_difference_derivative(mod_grad_sqd_fcn(:,:,:), fcn, 1, d, &
            order, fd_grid, cell_grid_info, .true.)
                                    ! jd: 'true' asks to square the result
                                    !     otherwise we'd need a temporary here
    end do

    call utils_trace_out(myself)

  end subroutine finite_difference_mod_grad_sqd
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!                 Helper subroutines for parallelization               !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine parallel_prepare_fd_halos(fcn,num_pad_rows,fd_grid,cell_grid_info)
    !=========================================================================!
    ! Prepares halos in the array fcn on the current grid.                    !
    ! Halos are required for parallel operation.                              !
    ! In serial mode dummy halos are allocated and filled with zeroes.        !
    !-------------------------------------------------------------------------!
    ! The halo buffers are allocated here, and it's up to the caller to free  !
    ! them up later by calling parallel_destroy_fd_halos (this is what        !
    ! finite_difference_derivative() does). Another buffer ('halo_buffer') is !
    ! only allocated on root. Its size can be controlled through a constant   !
    ! 'bufsize_in_megs', which is 100MB by default. This buffer is automati-  !
    ! cally destroyed on exit from this subroutine.                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   fcn (in)            Array containing the function in usual distributed!
    !                       representation.                                   !
    !   num_pad_rows (in)   Thickness of the halo (in slabs).                 !
    !   fd_grid      (in)   The grid on which to operate.                     !
    !   cell_grid_info (in) a GRID_INFO structure holding info on the parallel!
    !                       distribution of slabs that fd_grid does not hold. !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2010.                                      !
    ! Rewritten by Jacek Dziedzic in 2013.04 to account for a situation where !
    ! a halo may need to be filled with data from more than one proc.         !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_proc_id, pub_total_num_procs, &
         comms_barrier, comms_send, comms_irecv, comms_wait
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_trace_in, utils_trace_out, &
         utils_assert

    implicit none

    ! jd: Arguments
    type(FD_GRID_INFO), intent(in) :: fd_grid
    type(GRID_INFO), intent(in)    :: cell_grid_info
    real(kind=DP), intent(in), &
         dimension(fd_grid%ld1f,fd_grid%ld2f,fd_grid%max_slabs12) :: fcn
    integer, intent(in) :: num_pad_rows

    ! jd: Local variables
    integer :: ierr ! jd: Error flag
#ifdef MPI
    integer :: elements_in_slab               ! ld1f x ld2f
    integer :: owner
    integer :: halo_slab, my_slab
    integer :: my_first_slab_global
    integer :: this_slab_global
    integer :: this_slab_local
    integer :: proc
    integer :: first_slab_on_proc, last_slab_on_proc
    !integer :: delta_top, delta_bot ! how far is this slab f/m slabs of cur proc
    integer :: request_procs_bot(1:num_pad_rows)
    integer :: request_slabs_bot(1:num_pad_rows)
    integer :: request_handles_bot(1:num_pad_rows)
    integer :: request_procs_top(-num_pad_rows:-1)
    integer :: request_slabs_top(-num_pad_rows:-1)
    integer :: request_handles_top(-num_pad_rows:-1)
    integer, parameter :: max_request_ids = 50000
    integer :: satisfy_procs(1:max_request_ids)
    integer :: satisfy_slabs(1:max_request_ids)
    integer :: satisfy_handles(1:max_request_ids)
    integer :: request_id, n_request_ids

#endif

    character(len=*), parameter :: myself = 'parallel_prepare_fd_halos'

    !------------------------------------------------------------------------

#ifdef MPI
    call comms_barrier ! jd: Start timers only once we're all in sync
#endif

    call timer_clock(myself,1)
    call utils_trace_in(myself)

    call timer_clock(myself//"_alloc",1)

    ! jd: Allocate halo arrays, num_pad_rows thick
    allocate(halo_top(fd_grid%ld1f,fd_grid%ld2f,-num_pad_rows:-1),stat=ierr)
    call utils_alloc_check(myself,'halo_top',ierr)
    allocate(halo_bot(fd_grid%ld1f,fd_grid%ld2f,1:num_pad_rows),stat=ierr)
    call utils_alloc_check(myself,'halo_bot',ierr)

    call timer_clock(myself//"_alloc",2)

    call timer_clock(myself//"_clear",1)

    halo_top = 0.0_DP
    halo_bot = 0.0_DP

    call timer_clock(myself//"_clear",2)

#ifdef MPI

    call timer_clock(myself//"_request",1)

    elements_in_slab = fd_grid%ld1f * fd_grid%ld2f
    my_first_slab_global = cell_grid_info%first_slab12(pub_my_proc_id)

    ! ----------------------------------------------
    ! --- Figuring out requests this proc issues ---
    ! ----------------------------------------------

    ! jd: Come up with a list of procs that have subsequent slabs that we
    !     need in our BOTTOM halo
    do halo_slab = 1, num_pad_rows
       ! Determine the global index of this halo slab and then the owner
       ! (or -1 if slab outside the grid)
       this_slab_global = &
            my_first_slab_global + fd_grid%num_slabs12_pq + halo_slab - 1
       ! Determine who the owner is
       if(this_slab_global <= fd_grid%pq3f) then
          owner = cell_grid_info%proc_slab12(this_slab_global)
       else
          owner = -1
       end if
       request_procs_bot(halo_slab) = owner
       request_slabs_bot(halo_slab) = this_slab_global
    end do

    ! jd: Come up with a list of procs that have subsequent slabs that we
    !     need in our TOP halo
    do halo_slab = -num_pad_rows, -1
       ! Determine the global index of this halo slab and then the owner
       ! (or -1 if slab outside the grid)
       this_slab_global = my_first_slab_global + halo_slab
       ! Determine who the owner is
       if(this_slab_global > 0) then
          owner = cell_grid_info%proc_slab12(this_slab_global)
       else
          owner = -1
       end if
       request_procs_top(halo_slab) = owner
       request_slabs_top(halo_slab) = this_slab_global
    end do

    ! -----------------------------------
    ! --- Issue the requests (irecvs) ---
    ! -----------------------------------

    ! jd: Issue irecvs for these, straight into the BOTTOM halos
    do halo_slab = 1, num_pad_rows
       if(request_procs_bot(halo_slab) == -1) cycle
       ! write(*,'(a,i0,a,i0,a,i0,a,i0)') 'Proc ',pub_my_proc_id,' asks hs ', &
       !     halo_slab,' (',request_slabs_bot(halo_slab),') from #', &
       !     request_procs_bot(halo_slab)
       call comms_irecv(request_procs_bot(halo_slab), &
            halo_bot(:,:,halo_slab), elements_in_slab, &
            tag=request_slabs_bot(halo_slab), &
            handle=request_handles_bot(halo_slab))
    end do

    ! jd: Issue irecvs for these, straight into the TOP halos
    do halo_slab = -num_pad_rows, -1
       if(request_procs_top(halo_slab) == -1) cycle
       ! write(*,'(a,i0,a,i0,a,i0,a,i0)') 'Proc ',pub_my_proc_id,' asks hs ', &
       !     halo_slab,' (',request_slabs_top(halo_slab),') from #', &
       !     request_procs_top(halo_slab)
       call comms_irecv(request_procs_top(halo_slab), &
            halo_top(:,:,halo_slab), elements_in_slab, &
            tag=request_slabs_top(halo_slab), &
            handle=request_handles_top(halo_slab))
    end do

    call timer_clock(myself//"_request",2)

    ! ----------------------------------------------------------
    ! --- Figuring out what requests this proc is to satisfy ---
    ! ----------------------------------------------------------

    call timer_clock(myself//"_satisfy",1)

    ! jd: Come up with a list of procs that are going to need the slabs that
    !     we have for their halos. Each slab may be needed by more than one
    !     proc (or none).
    request_id = 1
    do my_slab = 1, fd_grid%num_slabs12_pq
       this_slab_global = my_first_slab_global + my_slab - 1
       do proc = 0, pub_total_num_procs-1
          first_slab_on_proc = cell_grid_info%first_slab12(proc)
          last_slab_on_proc = cell_grid_info%last_slab12(proc)

          ! Does this guy need our slab for his TOP halo?
          if(first_slab_on_proc > this_slab_global .and. &
               first_slab_on_proc - num_pad_rows <= this_slab_global) then
             satisfy_procs(request_id) = proc
             satisfy_slabs(request_id) = this_slab_global
             request_id = request_id + 1
          end if

          ! Does this guy need our slab for his BOTTOM halo?
          if(last_slab_on_proc < this_slab_global .and. &
               last_slab_on_proc + num_pad_rows >= this_slab_global) then
             ! This guy needs our slab
             satisfy_procs(request_id) = proc
             satisfy_slabs(request_id) = this_slab_global
             request_id = request_id + 1
          end if

          call utils_assert(request_id <= max_request_ids, &
               'max_request_ids insufficient in '//myself)
       end do
    end do

    ! ----------------------------------------------------
    ! --- Satisfy requests issued to this proc (sends) ---
    ! ----------------------------------------------------

    n_request_ids = request_id-1
    do request_id = 1, n_request_ids
       this_slab_local = satisfy_slabs(request_id) - my_first_slab_global+1
       !write(*,'(a,i0,a,i0,a,i0)') 'Proc ',pub_my_proc_id,' sends (', &
       !satisfy_slabs(request_id),') to #', satisfy_procs(request_id)

       call comms_send(satisfy_procs(request_id), fcn(:,:,this_slab_local),&
            elements_in_slab,tag=satisfy_slabs(request_id), &
            return_handle=satisfy_handles(request_id), add_to_stack=.false.)
    end do

    call timer_clock(myself//"_satisfy",2)

    ! -----------------------------------------------------
    ! --- Wait for completion of the receives and sends ---
    ! -----------------------------------------------------

    call timer_clock(myself//"_comms_wait",1)

    ! - receives for the BOTTOM halos
    do halo_slab = 1, num_pad_rows
       if(request_procs_bot(halo_slab) == -1) cycle
       call comms_wait(request_handles_bot(halo_slab))
    end do

    ! - receives for the TOP halos
    do halo_slab = -num_pad_rows, -1
       if(request_procs_top(halo_slab) == -1) cycle
       call comms_wait(request_handles_top(halo_slab))
    end do

    ! - sends (all)
    do request_id = 1, n_request_ids
       call comms_wait(satisfy_handles(request_id))
    end do

    call timer_clock(myself//"_comms_wait",2)

    call timer_clock(myself//"_comms_imbalance",1)
    call comms_barrier
    call timer_clock(myself//"_comms_imbalance",2)

#else
    ! jd: Pointless, but kills the 'fcn unused' warning
    ierr = int(fcn(1,1,1))
#endif

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine parallel_prepare_fd_halos
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine parallel_destroy_fd_halos

    use utils, only: utils_dealloc_check, utils_trace_in, utils_trace_out

    implicit none

    ! jd: Internal variables
    integer :: ierr ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('parallel_destroy_fd_halos')

    ! jd: Deallocate the halo arrays
    !     NB: This takes place regardless of whether we're in MPI or serial mode
    deallocate(halo_bot,stat=ierr)
    call utils_dealloc_check('parallel_destroy_fd_halos','halo_bot',ierr)
    deallocate(halo_top,stat=ierr)
    call utils_dealloc_check('parallel_destroy_fd_halos','halo_top',ierr)

    call utils_trace_out('parallel_destroy_fd_halos')

  end subroutine parallel_destroy_fd_halos
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module finite_differences

