!> \brief computes higher order finite differences in parallel
!!
!! Adapted from ONETEP's
!! finite_differences module.
!!
!> Written by Jacek Dziedzic in 08/2010.
!!
!> Modified for use within DL_MG by J. C. Womack, 2016
!!
!> Optimised by Lucian Anton, 2017
!!
!! This module defines the subroutines for finite-difference calculations of
!! the gradient, squared-modulus-of-gradient and the laplacian for a quantity
!! on a GRID_INFO grid.
!!
!! Finite differences of orders {2,4,6,8,10,12} can be used. Where 2nd order
!! formulas are unavailable, 12th order is used instead. Slab-parallelism is
!! automatically taken care of by internally using halos. Grid boundaries are
!! taken care of by using progressively more forward- or backward-FD's.


module dl_mg_defco_fd
! Original module documentation:


  use dl_mg_params, only: wp
  use dl_mg_common_data, only: &
        nxb, nyb, nzb, &                ! Global grid dimensions
        npx, npy, npz, &                ! Number of MPI ranks in each direction
        dx, dy, dz                      ! Per-grid-point step in each Cartesian direction
  use dl_mg_grids, only: fd


  implicit none
  private

  ! --------------------------
  ! --- Public subroutines ---
  ! --------------------------

  public :: dl_mg_defco_fd_initialise
  !public :: dl_mg_defco_fd_set_geometry
  public :: dl_mg_defco_fd_gradient
  public :: dl_mg_defco_fd_mod_grad_sqd
  public :: dl_mg_defco_fd_laplacian

  ! ---------------------------------------------------------------------------
  ! ------------------- P U B L I C   V A R I A B L E S  ----------------------
  ! ---------------------------------------------------------------------------

  ! jd: Maximum finite difference order that is supported
  integer, parameter, public :: max_order = 12

  ! ---------------------------------------------------------------------------
  ! ------------------------------ P R I V A T E ------------------------------
  ! ---------------------------------------------------------------------------

  ! jd: Initialisation flags
  logical :: dl_mg_defco_fd_initialised = .false.

  ! jd: Halos for finite differences in parallel mode
  real(kind=wp), allocatable :: halo_left(:,:,:)
  real(kind=wp), allocatable :: halo_right(:,:,:)

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
  real(kind=wp), dimension(-max_order:max_order, &
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

  ! coeff(:,1:6,12,1) will be initialised in dl_mg_defco_fd_initialise

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

  ! coeff(:,1:5,10,1) will be initialised in dl_mg_defco_fd_initialise

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
  ! coeff(:,1:4,8,1) will be initialised in dl_mg_defco_fd_initialise

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
  ! coeff(:,1:3,6,1) will be initialised in dl_mg_defco_fd_initialise

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
  ! coeff(:,1:2,4,1) will be initialised in dl_mg_defco_fd_initialise

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
  ! coeff(:,1,2,1) will be initialised in dl_mg_defco_fd_initialise

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

  ! coeff(:,1:6,12,2) will be initialised in dl_mg_defco_fd_initialise

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

  ! coeff(:,1:5,10,2) will be initialised in dl_mg_defco_fd_initialise

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

  ! coeff(:,1:4,8,2) will be initialised in dl_mg_defco_fd_initialise

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

  ! coeff(:,1:3,6,2) will be initialised in dl_mg_defco_fd_initialise

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

  ! coeff(:,1:2,4,2) will be initialised in dl_mg_defco_fd_initialise

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

  ! coeff(:,1,2,2) will be initialised in dl_mg_defco_fd_initialise

  ! -- Orders that we don't care about ---------------------------------------
  DATA coeff(:,:,1,2) / ndata*0D0/
  DATA coeff(:,:,3,2) / ndata*0D0/
  DATA coeff(:,:,5,2) / ndata*0D0/
  DATA coeff(:,:,7,2) / ndata*0D0/
  DATA coeff(:,:,9,2) / ndata*0D0/
  DATA coeff(:,:,11,2) / ndata*0D0/

  ! ##########################################################################

  ! jd: Norms for the finite differences of orders 1..12 for 1st and 2nd deriv.
  real(kind=wp), dimension(max_order,2), parameter :: norms= &
       reshape( source = &
       (/0D0, 1D0/2D0, 0D0, 1D0/12D0, 0D0, 1D0/60D0, 0D0, 1D0/840D0, &
       0D0, 1D0/2520D0, 0D0, 1D0/27720D0, &
       0D0, 1D0, 0D0, 1D0/12D0,0D0,1D0/180D0,0D0,1D0/5040D0,0D0,1D0/25200D0,&
       0D0,1D0/831600D0 /), shape = (/max_order, 2/) )

contains

  !>
  !! Initialises the finite differences by completing the computation of the
  !! finite difference coefficients that have been half-prepared by the DATA
  !! statements.
  !! There is no need to call this subroutine explicitly -- it will be called
  !! upon first use of any computational FD subroutine. You can call it
  !! explicitly if you like, and any number of times.
  !!
  !! Written by Jacek Dziedzic, 08/2010
  !!
  subroutine dl_mg_defco_fd_initialise()


    use dl_mg_defco_utils, only: dl_mg_defco_utils_abort, dl_mg_defco_utils_assert

    implicit none

    ! jd: Internal variables
    integer :: i, j, k, l ! jd: Indices
    integer :: ierror
    character(len=256) :: problem_at

    !------------------------------------------------------------------------

    if(dl_mg_defco_fd_initialised) return

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
                call dl_mg_defco_utils_abort('Inconsistency in the FD coefficient array &
                     &detected in dl_mg_defco_fd_mod, check '//&
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

    ! JCW: Removed sanity checks for moving code into DL_MG (these scanned
    ! JCW: arrays for occurences of NaN)
    ! TODO Consider replacing sanity checks with DL_MG native routines
    !      Santity checks performed for:
    !         coeff(:,:,:,1:2)

    dl_mg_defco_fd_initialised = .true.

  end subroutine dl_mg_defco_fd_initialise


  subroutine dl_mg_defco_fd_derivative(fcn_der,fcn,n,d,order, &
       square_result)
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
    !   square_result (in, optional) If true, squares of derivatives will be  !
    !                                returned (added). Not very elegant, but  !
    !                                allows efficient approach for |grad f|^2.!
    !-------------------------------------------------------------------------!
    ! All the arrays are assumed to be in the distributed slab representation.!
    !                                                                         !
    ! Note that this subroutine operates on the passed FD grid (obtained from !
    ! dl_mg_defco_fd_set_geometry, which might be slightly smaller than    !
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

    use dl_mg_common_data, only: nthreads, isx, isy, isz, bc
    use dl_mg_grids, only:       fd
    use dl_mg_timer
    ! JCW: "use [...], only" syntax is problematic with mpi headers, so
    ! JCW: just import entire dl_mg_mpi_header module
    use dl_mg_mpi_header
    use dl_mg_mpi, only : get_myid
    !use dl_mg_mpi_header, only: MPI_ABORT, MPI_COMM_RANK, MPI_SUCCESS
    use dl_mg_defco_utils, only: dl_mg_defco_utils_assert, dl_mg_defco_utils_abort
    use dl_mg_params, only : dl_mg_bc_periodic
    implicit none

    ! jd: Arguments
    real(kind=wp), dimension(isx:,isy:,isz:), intent(inout) :: fcn_der !< output array
    real(kind=wp), dimension(isx:,isy:,isz:), intent(in)    :: fcn  !< input array
    integer, intent(in)            :: order !< finite difference order
    integer, intent(in)            :: n !< derivative order
    integer, intent(in)            :: d !< dimension (1=x, 2=y, 3=z)
    logical, intent(in), optional  :: square_result

    ! jd: Internal variables
    character(len=*), parameter :: myself = 'dl_mg_defco_fd_derivative'
    logical :: must_square_result  ! jd: .true. if square_result present and T
    integer :: i,j,k               ! jd: Indices
    real(kind=wp) :: ddf = 0.0_wp  ! jd: Temporary

    integer :: order_half          ! jd: order/2
    integer :: m, m_min, m_max, p
    real(kind=wp) :: value_here      ! jd: Temporary
    real(kind=wp) :: der             ! jd: Accumulator for results
    integer :: formula               ! jd: Number of relevant FD formula
    integer :: isxyz(3), iexyz(3)
    integer :: halo_s, halo_e
    real(kind=wp) :: coeff_times_inv_norm(-max_order:max_order, &
         -max_order_half:max_order_half) ! jd: Look-up
    real(kind=wp), allocatable :: fcn_slice(:) ! jd: Lookup
    !
    character(len=256) :: char_buffer
    integer :: myid
    integer :: ierr  ! error flag
    integer :: ierr1 ! to satisfy MPI_ABORT

    !------------------------------------------------------------------------

    ! TODO Insert DL_MG timer start call
    call mg_timer(start, tdefco, tignore, tignore)

    ! TODO Replace per-call determination of MPI process id with single
    !      initialization.
    ! JCW: Determine MPI process id

    myid = get_myid(fd%comm)

    call dl_mg_defco_fd_initialise

    must_square_result = .false.
    if(present(square_result)) then
       must_square_result = square_result
    end if

    ! jd: Sanity check on the direction, n, order and input and grid size
    if (.not. (d>0 .and. d<=3)) then
       char_buffer = ""
       write(char_buffer,'(i32)') d
       call dl_mg_defco_utils_assert(.false., &
            'Unrecognized direction in '//myself//", d="//trim(char_buffer))
    end if
    if (.not. (n==1 .or. n==2)) then
       char_buffer = ""
       write(char_buffer,'(i32)') n
       call dl_mg_defco_utils_assert(.false., &
            myself//': can only compute 1st and 2nd derivatives, n='//trim(char_buffer))
    end if

    order_half = order / 2

    ! la: skectch for computation communication overlap with threads

    ! allocate halos

    ! fill the send haloes

    ! post send recv (master only)

    ! add contribution from the local domain (rest of the threads)

    !   x direction

    !   y direction

    !   z direction

    ! master thread waits for buffers

    ! add contributions from the halos where needed

! ==================================================================================
    ! jd: Prepare parallel halos for fcn
    call parallel_prepare_fd_halos(fcn, d)

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
    select case(d)
    case(1)
       ddf=dx
    case(2)
       ddf=dy
    case(3)
       ddf=dz
    end select

    ! jd: Calculate the product of relevant coefficients and 1/ddf,
    !     store these in a lookup to save time in the inner loops
    coeff_times_inv_norm(:,:) = coeff(:,:,order,n) * 1.0_wp/(ddf**n)

    i=0; j=0; k=0 ! jd: Silence compiler warnings


!$OMP PARALLEL NUM_THREADS(nthreads) DEFAULT(NONE) &
!$OMP PRIVATE(i,j,k,halo_s, halo_e, isxyz, iexyz, &
!$OMP      formula, m, p, m_max,m_min) &
!$OMP SHARED(d,fcn, halo_left,halo_right,order_half,must_square_result, &
!$OMP      order,coeff_times_inv_norm,myid,fcn_der,nthreads, &
!$OMP      fd,bc)

    isxyz(1) = fd%sx
    isxyz(2) = fd%sy
    isxyz(3) = fd%sz
    iexyz(1) = fd%ex
    iexyz(2) = fd%ey
    iexyz(3) = fd%ez

    if (fd%num_neighbours_recv_l(d) > 0) then
       halo_s = fd%halo_indices(2*d-1,1,d) !fd%halo_map_recv_l(d)%m(1, fd%num_neighbours_recv_l(d))
    else
       halo_s = isxyz(d)
    endif
    if (fd%num_neighbours_recv_r(d) > 0) then
       halo_e = fd%halo_indices(2*d,2,d) !fd%halo_map_recv_r(d)%m(2, fd%num_neighbours_recv_r(d))
    else
       halo_e = iexyz(d)
    endif

    ! la: the cases from below probable can be aggregated in 1 blocvk with some index logic
    ! la: however this more explict form could help the compiter to treat the loop nest
    ! la: also it may help with computation comm overlap (to be done)
    ! la: highly improbable to need this at arbitrary dimension

    select case(d)
    case(1)
       !$OMP DO COLLAPSE(3)
       do k = fd%sz, fd%ez
          do j = fd%sy, fd%ey
             do i = fd%sx, fd%ex
                formula = 0
                if (bc(d) /= DL_MG_BC_PERIODIC)  then
                   if(i - halo_s < order_half) formula = order_half -(i - halo_s)
                   if(halo_e -i  < order_half) formula = (halo_e -i) - order_half
                end if

                m_min=formula-order_half
                m_max=m_min+order
                do m = max(m_min, fd%sx -i), min(m_max, fd%ex-i)
                   fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                        fcn(i+m,j,k)
                end do
#ifdef MPI
                do m = m_min, fd%sx - i -1
                   fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                        halo_left(i+m,j,k)
                end do
                do m = fd%ex -i +1, m_max
                   fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                        halo_right(i+m,j,k)
                end do
#else
                ! la: without MPI and PBC the halo data are taken from the other end of fcn grid
                if (bc(d) == DL_MG_BC_PERIODIC) then
                   if ( i - fd%sx < order_half) then
                      p = - (fd%sx - i -1 - m_min +1 ) +1
                      do m = m_min, fd%sx - i -1
                         fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                              fcn(fd%ex+p,j,k)
                         p = p + 1
                      end do
                   endif
                   if (fd%ex -i < order_half) then
                      p = 0
                      do m = fd%ex -i +1, m_max
                         fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                              fcn(fd%sx + p,j,k)
                         p = p + 1
                      end do
                   endif
                end if
#endif
             end do
          enddo
       enddo
       !$OMP ENDDO

    case(2)
       !$OMP DO COLLAPSE(2)
       do k = fd%sz, fd%ez
          do j = fd%sy, fd%ey
             formula = 0
             if (bc(d) /= DL_MG_BC_PERIODIC)  then
                if(j - halo_s < order_half) formula = order_half -(j - halo_s)
                if(halo_e -j  < order_half) formula = (halo_e -j) - order_half
             end if
             m_min=formula-order_half
             m_max=m_min+order
             do i = fd%sx, fd%ex
                do m= max(m_min, fd%sy-j), min(m_max, fd%ey-j)
                   fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                        fcn(i,j+m,k)
                end do
#ifdef MPI
                do m = m_min, fd%sy-j-1
                   fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                        halo_left(i,j+m,k)
                end do
                do m=fd%ey-j+1, m_max
                   fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                        halo_right(i,j+m,k)
                end do
#else
                ! la: without MPI and PBC the halo data are taken from the other end of fcn grid
                if (bc(d) == DL_MG_BC_PERIODIC) then
                   if ( j - fd%sy < order_half) then
                      p = - (fd%sy - j -1 - m_min)
                      do m = m_min, fd%sy - j -1
                         fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                              fcn(i, fd%ey + p, k)
                         p = p + 1
                      end do
                   endif
                   if(fd%ey -j < order_half)then
                      p = 0
                      do m = fd%ey -j +1, m_max
                         fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                              fcn(i,fd%sy+p,k)
                         p = p + 1
                      end do
                   endif
                end if
#endif
             end do
          enddo
       enddo
       !$OMP ENDDO
    case(3)
       ! la: shold thew openmp loop be around i_j loops?
       !$OMP DO
       do k = fd%sz, fd%ez
          formula = 0
          if (bc(d) /= DL_MG_BC_PERIODIC)  then
             if(k - halo_s < order_half) formula = order_half -(k - halo_s)
             if(halo_e -k  < order_half) formula = (halo_e -k) - order_half
          end if
          m_min=formula-order_half
          m_max=m_min+order
          do j = fd%sy, fd%ey
             do i = fd%sx, fd%ex
                do m = max(m_min, fd%sz-k), min(m_max, fd%ez-k)
                   fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                        fcn(i,j,k+m)
                end do
#ifdef MPI
                do m = m_min, fd%sz-k-1
                   fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                        halo_left(i,j,k+m)
                end do
                do m = fd%ez-k+1,  m_max
                   fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                        halo_right(i,j,k+m)
                end do
#else
                ! la: without MPI and PBC the halo data are taken from the other end of fcn grid
                if (bc(d) == DL_MG_BC_PERIODIC) then
                   if ( k - fd%sz < order_half) then
                      p = - (fd%sz - k -1 - m_min)
                      do m = m_min, fd%sz - k -1
                         fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                              fcn(i, j, fd%ez + p)
                         p = p + 1
                      end do
                   endif
                   if(fd%ez -k < order_half)then
                      p = 0
                      do m = fd%ez -k +1, m_max
                         fcn_der(i,j,k) = fcn_der(i,j,k) + coeff_times_inv_norm (m,formula) * &
                              fcn(i,j,fd%sz+p)
                         p = p + 1
                      end do
                   endif
                end if
#endif
             end do
          enddo
       enddo
       !$OMP ENDDO
    end select

    if(must_square_result) then
       !$OMP DO COLLAPSE(3)
       do k = fd%sz, fd%ez
          do j = fd%sy, fd%ey
             do i = fd%sx, fd%ex
                fcn_der(i,j,k) = fcn_der(i,j,k) * fcn_der(i,j,k)
             enddo
          enddo
       enddo
       !$OMP ENDDO NOWAIT
    end if
    !$OMP END PARALLEL

    ! jd: Clean up by destroying the halo arrays
    call parallel_destroy_fd_halos

    ! TODO Insert DL_MG timer end call
    call mg_timer(stop, tdefco, tignore, tignore)

  end subroutine dl_mg_defco_fd_derivative


  subroutine dl_mg_defco_fd_gradient(grad_fcn,fcn,order)
    !=========================================================================!
    ! This subroutine takes a function 'fcn' on the real space grid and uses  !
    ! finite difference methods to find the gradient of the function on the   !
    ! grid.  The gradient is returned in grad_fcn.                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   grad_fcn (out) the gradient of fcn found by finite differences.       !
    !   fcn      (in)  the scalar function whose gradient we're after.        !
    !   order    (in)  the finite difference order to use.                    !
    !-------------------------------------------------------------------------!
    ! All the arrays are assumed to be in the distributed slab representation.!
    !                                                                         !
    ! Note that this subroutine operates on the passed FD grid (obtained from !
    ! dl_mg_defco_fd_set_geometry, which might be slightly smaller than    !
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

    use dl_mg_defco_utils, only: dl_mg_defco_utils_assert, dl_mg_defco_utils_abort

    implicit none

    ! jd: Arguments
    real(kind=wp), dimension(:,:,:,:), intent(out) :: grad_fcn
    ! (Last dimension should be Cartesian coordinate, 1:3)
    real(kind=wp), dimension(:,:,:), intent(in)    :: fcn
    integer, intent(in)            :: order

    ! jd: Local variables
    integer :: d
    character(len=*), parameter :: myself = 'dl_mg_defco_fd_gradient'

    !------------------------------------------------------------------------

    call dl_mg_defco_fd_initialise

    grad_fcn = 0D0
    do d=1, 3
       call dl_mg_defco_fd_derivative(grad_fcn(:,:,:,d), fcn, 1, d, order)
    end do

  end subroutine dl_mg_defco_fd_gradient


  subroutine dl_mg_defco_fd_laplacian(lap_fcn,fcn,order)
    !=========================================================================!
    ! This subroutine takes a function 'fcn' on the real space grid and uses  !
    ! finite difference methods to find the laplacian of the function on the  !
    ! grid.  The laplacian is returned in lap_fcn.                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   lap_fcn (out)  the laplacian of fcn found by finite differences.      !
    !   fcn      (in)  the scalar function whose gradient we're after.        !
    !   order    (in)  the finite difference order to use.                    !
    !-------------------------------------------------------------------------!
    ! All the arrays are assumed to be in the distributed slab representation.!
    !                                                                         !
    ! Note that this subroutine operates on the passed FD grid (obtained from !
    ! dl_mg_defco_fd_set_geometry, which might be slightly smaller than    !
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

    use dl_mg_defco_utils, only: dl_mg_defco_utils_assert, dl_mg_defco_utils_abort

    implicit none

    ! jd: Arguments
    real(kind=wp), dimension(:,:,:), intent(out) :: lap_fcn
    real(kind=wp), dimension(:,:,:), intent(in)  :: fcn
    integer, intent(in)            :: order

    ! jd: Local variables
    integer :: d
    character(len=*), parameter :: myself = 'dl_mg_defco_fd_laplacian'

    !------------------------------------------------------------------------

    call dl_mg_defco_fd_initialise

    lap_fcn = 0D0
    do d=1, 3
       ! NB: dl_mg_defco_fd_derivative *adds* to lap_fcn
       call dl_mg_defco_fd_derivative(lap_fcn(:,:,:), fcn, 2, d, order)

    end do

  end subroutine dl_mg_defco_fd_laplacian


  subroutine dl_mg_defco_fd_mod_grad_sqd(mod_grad_sqd_fcn,fcn,order)
    !=========================================================================!
    ! This subroutine takes a function 'fcn' on the real space grid and uses  !
    ! finite difference methods to find the |grad fcn|^2 of the function on   !
    ! the grid, which is returned in grad_mod_sqd_fcn.                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mod_grad_sqd_fcn (out)  the result found by finite differences.       !
    !   fcn      (in)  the scalar function whose gradient we're after.        !
    !   order    (in)  the finite difference order to use.                    !
    !-------------------------------------------------------------------------!
    ! All the arrays are assumed to be in the distributed slab representation.!
    !                                                                         !
    ! Note that this subroutine operates on the passed FD grid (obtained from !
    ! dl_mg_defco_fd_set_geometry, which might be slightly smaller than    !
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

    use dl_mg_defco_utils, only: dl_mg_defco_utils_assert, dl_mg_defco_utils_abort

    implicit none

    ! jd: Arguments
    real(kind=wp), dimension(:,:,:), intent(out) :: mod_grad_sqd_fcn
    real(kind=wp), dimension(:,:,:), intent(in)  :: fcn
    integer, intent(in)            :: order

    ! jd: Local variables
    integer :: d
    character(len=*), parameter :: myself = 'dl_mg_defco_fd_mod_grad_sqd'

    !------------------------------------------------------------------------

    call dl_mg_defco_fd_initialise

    mod_grad_sqd_fcn = 0D0
    do d=1, 3
       call dl_mg_defco_fd_derivative(mod_grad_sqd_fcn(:,:,:), fcn, 1, d, &
            order, .true.)
                                    ! jd: 'true' asks to square the result
                                    !     otherwise we'd need a temporary here
    end do

  end subroutine dl_mg_defco_fd_mod_grad_sqd


  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!                 Helper subroutines for parallelization               !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

  subroutine parallel_prepare_fd_halos(fcn, axis)
    !=========================================================================!
    ! Prepares halos in the array fcn on the current grid.                    !
    ! Halos are required for parallel operation.                              !
    ! In serial mode dummy halos are allocated and filled with zeroes.        !
    !-------------------------------------------------------------------------!
    ! The halo buffers are allocated here, and it's up to the caller to free  !
    ! them up later by calling parallel_destroy_fd_halos (this is what        !
    ! dl_mg_defco_fd_derivative() does). Another buffer ('halo_buffer') is !
    ! only allocated on root. Its size can be controlled through a constant   !
    ! 'bufsize_in_megs', which is 100MB by default. This buffer is automati-  !
    ! cally destroyed on exit from this subroutine.                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   fcn (in)            Array containing the function in usual distributed!
    !                       representation.                                   !
    !   num_pad_rows (in)   Thickness of the halo (in slabs).                 !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 2010.                                      !
    ! Rewritten by Jacek Dziedzic in 2013.04 to account for a situation where !
    ! a halo may need to be filled with data from more than one node.         !
    !=========================================================================!

    use dl_mg_grids, only : fd
    use dl_mg_mpi_header
    use dl_mg_timer
    use dl_mg_defco_utils, only: dl_mg_defco_utils_alloc_check, dl_mg_defco_utils_assert, &
         dl_mg_defco_tag_ub
    use dl_mg_common_data, only: isx, isy, isz
    use dl_mg_errors, only: handle_error, DL_MG_ERR_UNSPECIFIED
    implicit none

    ! jd: Arguments
    real(kind=wp), intent(in), dimension(isx:,isy:,isz:) :: fcn
    integer, intent(in) :: axis !< xyz

    ! local params
    integer, parameter :: sendK=0, recvK=1, leftK=-1, rightK=1
    type buff_t
       real(wp), allocatable :: b(:,:,:)
    end type buff_t


    ! jd: Local variables
    integer :: ierr ! jd: Error flag
    integer :: ierr1 ! JCW: to satisfy MPI_ABORT
#ifdef MPI
    integer :: hil(6), hir(6), i1s, i1e, i2s, i2e, i3s, i3e
    integer :: isend, irecv, total_sends, total_recvs, idx
    integer :: k, count, ntot, dest, src, tag
    integer, allocatable :: ireq(:)
    type(buff_t), allocatable :: buffsend(:), buffrecv(:)
#endif

    character(len=*), parameter :: myself = 'parallel_prepare_fd_halos'


    !------------------------------------------------------------------------

    call mg_timer(start, tdefco_mpi, tignore, tignore)

#ifdef MPI

    hil = fd%halo_indices(:,1,axis)
    hir = fd%halo_indices(:,2,axis)

    allocate(halo_left(hil(1):hil(2),  &
                       hil(3):hil(4),  &
                       hil(5):hil(6)), &
            halo_right(hir(1):hir(2),  &
                       hir(3):hir(4),  &
                       hir(5):hir(6)))

    total_sends = fd%num_neighbours_send_l(axis) + fd%num_neighbours_send_r(axis)
    total_recvs = fd%num_neighbours_recv_l(axis) + fd%num_neighbours_recv_r(axis)
    allocate(ireq(total_sends+total_recvs))

    allocate(buffsend(total_sends))
    ! send left-right
    do isend = 1, total_sends
       if ( isend <= fd%num_neighbours_send_l(axis)) then
          k = isend
          dest = fd%halo_map_send_l(axis)%m(3,k)
          call get_idx(axis, sendK, leftK, k, i1s, i1e, i2s, i2e, i3s, i3e)
          tag = 0
       else
          k= isend - fd%num_neighbours_send_l(axis)
          dest = fd%halo_map_send_r(axis)%m(3,k)
          call get_idx(axis, sendK, rightK, k, i1s, i1e, i2s, i2e, i3s, i3e)
          tag = 1
       endif

       allocate(buffsend(isend)%b(i1s:i1e, i2s:i2e, i3s:i3e))
       ! no need to do this copy in z axis
       ! a pointer can be used here
       buffsend(isend)%b = fcn(i1s:i1e, i2s:i2e, i3s:i3e)
       ntot = size(buffsend(isend)%b)
       call mpi_isend(buffsend(isend)%b, ntot, MPI_DOUBLE_PRECISION, &
            dest, tag, fd%comm, ireq(isend), ierr)
    end do

    ! receive left right
    allocate(buffrecv(total_recvs))
    do irecv =1, total_recvs
       if ( irecv <= fd%num_neighbours_recv_l(axis)) then
          k = irecv
          src = fd%halo_map_recv_l(axis)%m(3,k)
          call get_idx(axis, recvK, leftK, k, i1s, i1e, i2s, i2e, i3s, i3e)
          tag = 1
       else
          k= irecv - fd%num_neighbours_recv_l(axis)
          src = fd%halo_map_recv_r(axis)%m(3,k)
          call get_idx(axis, recvK, rightK, k, i1s, i1e, i2s, i2e, i3s, i3e)
          tag = 0
       endif
       allocate(buffrecv(irecv)%b(i1s:i1e, i2s:i2e, i3s:i3e))
       ntot = size(buffrecv(irecv)%b)
       call mpi_irecv(buffrecv(irecv)%b, ntot, MPI_DOUBLE_PRECISION, &
            src, tag, fd%comm, ireq(total_sends + irecv), ierr) ! not tag needed as receiveing only one message from each src
    end do

    count = total_sends + total_recvs
    do while( count > 0)
       call mpi_waitany(total_sends + total_recvs, ireq, idx, MPI_STATUS_IGNORE, ierr)
       if (idx == MPI_UNDEFINED) then
          ! this shouldn't happen
          call handle_error(DL_MG_ERR_UNSPECIFIED, &
               msg="error in requests processing")
       endif
       count = count -1
       if (idx <= total_sends) then
          deallocate(buffsend(idx)%b)
       else
          if ( idx - total_sends <= fd%num_neighbours_recv_l(axis)) then
             k = idx  - total_sends
             call get_idx(axis, recvK, leftK, k, i1s, i1e, i2s, i2e, i3s, i3e)
             halo_left(i1s:i1e, i2s:i2e, i3s:i3e)= buffrecv(idx- total_sends)%b
             deallocate(buffrecv(idx - total_sends)%b)
          else
             k = idx - fd%num_neighbours_recv_l(axis) - total_sends
             call get_idx(axis, recvK, rightK, k, i1s, i1e, i2s, i2e, i3s, i3e)
             halo_right(i1s:i1e, i2s:i2e, i3s:i3e)= buffrecv(idx - total_sends)%b
             deallocate(buffrecv(idx - total_sends)%b)
          endif
       end if
    end do

#endif

    call mg_timer(stop, tdefco_mpi, tignore, tignore)

#ifdef MPI
  contains

    !> get indecies from fd data structure
    subroutine get_idx(d,  sendrecv, lr, k, i1s, i1e, i2s, i2e, i3s, i3e)
      use dl_mg_common_data, only : nx, ny, nz
      implicit none

      integer, intent(in ) :: d, sendrecv, lr, k
      integer, intent(out) :: i1s, i1e, i2s, i2e, i3s, i3e

      select case(sendrecv)
      case (sendK)
         select case(lr)
         case(leftK)
            select case(d)
            case(1)
               i1s = fd%halo_map_send_l(d)%m(1,k)
               i1e = fd%halo_map_send_l(d)%m(2,k)
               i2s = fd%sy
               i2e = fd%ey
               i3s = fd%sz
               i3e = fd%ez
            case(2)
               i1s = fd%sx
               i1e = fd%ex
               i2s = fd%halo_map_send_l(d)%m(1,k)
               i2e = fd%halo_map_send_l(d)%m(2,k)
               i3s = fd%sz
               i3e = fd%ez
            case(3)
               i1s = fd%sx
               i1e = fd%ex
               i2s = fd%sy
               i2e = fd%ey
               i3s = fd%halo_map_send_l(d)%m(1,k)
               i3e = fd%halo_map_send_l(d)%m(2,k)
            end select
         case(rightK)
            select case(d)
            case(1)
               i1s = fd%halo_map_send_r(d)%m(1,k)
               i1e = fd%halo_map_send_r(d)%m(2,k)
               i2s = fd%sy
               i2e = fd%ey
               i3s = fd%sz
               i3e = fd%ez
            case(2)
               i1s = fd%sx
               i1e = fd%ex
               i2s = fd%halo_map_send_r(d)%m(1,k)
               i2e = fd%halo_map_send_r(d)%m(2,k)
               i3s = fd%sz
               i3e = fd%ez
            case(3)
               i1s = fd%sx
               i1e = fd%ex
               i2s = fd%sy
               i2e = fd%ey
               i3s = fd%halo_map_send_r(d)%m(1,k)
               i3e = fd%halo_map_send_r(d)%m(2,k)
            end select
         end select
      case (recvK)
         select case(lr)
         case(leftK)
            select case(d)
            case(1)
               i1s = fd%halo_map_recv_l(d)%m(1,k)
               if (i1s > fd%sx) i1s = i1s - nx
               i1e = fd%halo_map_recv_l(d)%m(2,k)
               if (i1e > fd%sx) i1e = i1e - nx
               i2s = fd%sy
               i2e = fd%ey
               i3s = fd%sz
               i3e = fd%ez
            case(2)
               i1s = fd%sx
               i1e = fd%ex
               i2s = fd%halo_map_recv_l(d)%m(1,k)
               if (i2s > fd%sy) i2s = i2s - ny
               i2e = fd%halo_map_recv_l(d)%m(2,k)
               if (i2e > fd%sy) i2e = i2e - ny
               i3s = fd%sz
               i3e = fd%ez
            case(3)
               i1s = fd%sx
               i1e = fd%ex
               i2s = fd%sy
               i2e = fd%ey
               i3s = fd%halo_map_recv_l(d)%m(1,k)
               if (i3s > fd%sz) i3s = i3s - nz
               i3e = fd%halo_map_recv_l(d)%m(2,k)
               if (i3e > fd%sz) i3e = i3e - nz
            end select
         case(rightK)
            select case(d)
            case(1)
               i1s = fd%halo_map_recv_r(d)%m(1,k)
               if (i1s < fd%ex) i1s = i1s + nx
               i1e = fd%halo_map_recv_r(d)%m(2,k)
               if (i1e < fd%ex) i1e = i1e + nx
               i2s = fd%sy
               i2e = fd%ey
               i3s = fd%sz
               i3e = fd%ez
            case(2)
               i1s = fd%sx
               i1e = fd%ex
               i2s = fd%halo_map_recv_r(d)%m(1,k)
               if (i2s < fd%ey) i2s = i2s + ny
               i2e = fd%halo_map_recv_r(d)%m(2,k)
               if (i2e < fd%ey) i2e = i2e + ny
               i3s = fd%sz
               i3e = fd%ez
            case(3)
               i1s = fd%sx
               i1e = fd%ex
               i2s = fd%sy
               i2e = fd%ey
               i3s = fd%halo_map_recv_r(d)%m(1,k)
               if (i3s < fd%ez) i3s = i3s + nz
               i3e = fd%halo_map_recv_r(d)%m(2,k)
               if (i3e < fd%ez) i3e = i3e + nz
            end select
         end select
      end select
    end subroutine get_idx
#endif
!!$#ifdef MPI
!!$    ! JCW: DL_MG routine which invokes mpi_barrier (was comms_barrier)
!!$    call barrier(fd%comm)
!!$#endif
!!$
!!$    ! JCW: Get required extents of dimensions of assumed shape fcn array
!!$    l1 = size(fcn,1); l2 = size(fcn,2)
!!$    ! JCW: These should be equal to are larger than the extents of the
!!$    ! JCW: portion of the array that we are operating on
!!$    call dl_mg_defco_utils_assert(l1>=(fd%ex-fd%sx+1).and.l2>=(fd%ey-fd%sy+1),&
!!$         "Error in "//myself//": &
!!$         &extents of input arrays are inconsistent with start and end indexes &
!!$         &passed to dl_mg_init.")
!!$
!!$    ! TODO Insert DL_MG timer start call
!!$
!!$    ! TODO Insert DL_MG timer start call (for allocation)
!!$
!!$#ifdef MPI
!!$    ! TODO Replace per-call determination of number of MPI ranks with
!!$    !      single initialization
!!$    ! JCW: Determine number of MPI processes
!!$    call MPI_COMM_SIZE(fd%comm, num_ranks , ierr)
!!$    if (ierr /= MPI_SUCCESS) then
!!$       write(*,*) "Error in "//myself//": MPI_COMM_SIZE call failed"
!!$       call MPI_ABORT(fd%comm,ierr,ierr1)
!!$    end if
!!$
!!$    ! TODO Replace per-call determination of MPI process id with single
!!$    !      initialization.
!!$    ! JCW: Determine MPI process id
!!$    call MPI_COMM_RANK(fd%comm, myid, ierr)
!!$    if (ierr /= MPI_SUCCESS) then
!!$       write(*,*) "Error in "//myself//": MPI_COMM_RANK call failed"
!!$       call MPI_ABORT(fd%comm,ierr,ierr1)
!!$    end if
!!$#else
!!$    num_ranks = 1
!!$    myid      = 0
!!$#endif
!!$
!!$    ! jd: Allocate halo arrays, num_pad_rows thick
!!$    allocate(halo_top(l1,l2,-num_pad_rows:-1),stat=ierr)
!!$    call dl_mg_defco_utils_alloc_check(ierr,myself,'halo_top')
!!$    allocate(halo_bot(l1,l2,1:num_pad_rows),stat=ierr)
!!$    call dl_mg_defco_utils_alloc_check(ierr,myself,'halo_bot')
!!$
!!$    ! TODO Insert DL_MG timer end call (for allocation)
!!$
!!$    ! TODO Insert DL_MG timer start call (for zeroing arrays)
!!$
!!$    halo_top = 0.0_wp
!!$    halo_bot = 0.0_wp
!!$
!!$    ! TODO Insert DL_MG timer end call (for zeroing arrays)
!!$
!!$#ifdef MPI
!!$
!!$    ! TODO Insert DL_MG timer start call (requests)
!!$
!!$    elements_in_slab = l1 * l2
!!$    !my_first_slab_global = cell_grid_info%first_slab12(myid)
!!$    my_first_slab_global = fd%sz
!!$
!!$    ! ----------------------------------------------
!!$    ! --- Figuring out requests this node issues ---
!!$    ! ----------------------------------------------
!!$
!!$    ! jd: Come up with a list of nodes that have subsequent slabs that we
!!$    !     need in our BOTTOM halo
!!$    ! Initialize all request_nodes_bot values to -1.
!!$    ! Any requests remaining with a value of -1 after the following if..end if
!!$    ! block will be skipped later (out of grid).
!!$    request_nodes_bot(:) = -1
!!$    if (fd%num_neighbours_r(3) == 0) then
!!$      ! No neighbours, so all slabs must be out of grid (set owner = -1 for
!!$      ! all slabs)
!!$    else
!!$      ineighbour = 1
!!$      neighbour_end_global = fd%ez + fd%halo_map_r(3)%m(ineighbour)
!!$      do halo_slab = 1, num_pad_rows
!!$        this_slab_global = fd%ez + halo_slab
!!$        if (this_slab_global > neighbour_end_global) then
!!$           ! Halo slab outside of this neighbour, move to next
!!$           ineighbour = ineighbour + 1
!!$           if (ineighbour > fd%num_neighbours_r(3)) then
!!$             ! Exceeded number of neighbours on this rank in the
!!$             ! right direction.
!!$             if (this_slab_global <= nzb) then
!!$               ! Remaining slabs should be outside of the grid.
!!$               write(*,'(a)') "Error in "//myself//": &
!!$                    &exhausted number of right neighbours but still inside grid."
!!$               call MPI_ABORT(fd%comm,ierr,ierr1)
!!$             else
!!$               ! End of grid, exit loop (all remaining request_nodes_bot
!!$               ! values have been initialized to -1).
!!$               exit
!!$             end if
!!$           end if
!!$           neighbour_end_global = neighbour_end_global + &
!!$               fd%halo_map_r(3)%m(ineighbour)
!!$        end if
!!$
!!$        ! Determine the global index of this halo slab and then the owner
!!$        ! (or -1 if slab outside the grid)
!!$        if (this_slab_global <= nzb) then
!!$          ! It is possible to have ineighbour > fd%num_neighbours_r, but
!!$          ! only when we are outside of the grid
!!$          if (ineighbour > fd%num_neighbours_r(3)) then
!!$            write(*,'(a)') "Error in "//myself//": ineighbour > fd%num_neighbours_r"
!!$            call MPI_ABORT(fd%comm,ierr,ierr1)
!!$          end if
!!$          owner = fd%neighbour_map_r(3)%m(ineighbour)
!!$        else
!!$          !owner = -1
!!$          write(*,'(a)') "Error in "//myself//": &
!!$               &Requesting slab outside of grid is not allowed (and should not &
!!$               &occur under normal circumstances)."
!!$          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$        end if
!!$        request_nodes_bot(halo_slab) = owner
!!$        request_slabs_bot(halo_slab) = this_slab_global
!!$        !write(*,*) "DL_MG DEBUG: request R ", myid, halo_slab, this_slab_global, owner, fd%num_neighbours_r(3), ineighbour
!!$      end do
!!$    end if
!!$
!!$
!!$    ! jd: Come up with a list of nodes that have subsequent slabs that we
!!$    !     need in our TOP halo
!!$    ! Initialize all request_nodes_top values to -1.
!!$    ! Any requests remaining with a value of -1 after the following if..end if
!!$    ! block will be skipped later (out of grid).
!!$    request_nodes_top(:) = -1
!!$    if (fd%num_neighbours_l(3) == 0) then
!!$      ! No neighbours, so all slabs must be out of grid (set owner = -1 for
!!$      ! all slabs)
!!$    else
!!$      ineighbour = 1
!!$      neighbour_start_global = fd%sz - fd%halo_map_l(3)%m(ineighbour)
!!$      do halo_slab = -1, -num_pad_rows, -1
!!$        this_slab_global = fd%sz + halo_slab
!!$        if (this_slab_global < neighbour_start_global) then
!!$           ! Halo slab outside of this neighbour, move to next
!!$           ineighbour = ineighbour + 1
!!$           if (ineighbour > fd%num_neighbours_l(3)) then
!!$             ! Exceeded number of neighbours on this rank in the
!!$             ! left direction.
!!$             if (this_slab_global > 0) then
!!$               ! Remaining slabs should be outside of the grid.
!!$               write(*,'(a)') "Error in "//myself//": &
!!$                    &exhausted number of left neighbours but still inside grid."
!!$               call MPI_ABORT(fd%comm,ierr,ierr1)
!!$             else
!!$               ! End of grid, exit loop (all remaining request_nodes_top
!!$               ! values have been initialized to -1)
!!$               exit
!!$             end if
!!$           end if
!!$           neighbour_start_global = neighbour_start_global - &
!!$               fd%halo_map_l(3)%m(ineighbour)
!!$        end if
!!$        ! Determine the global index of this halo slab and then the owner
!!$        ! (or -1 if slab outside the grid)
!!$        if (this_slab_global > 0) then
!!$          ! It is possible to have ineighbour > fd%num_neighbours_r, but
!!$          ! only when we are outside of the grid
!!$          call dl_mg_defco_utils_assert(&
!!$               ineighbour <= fd%num_neighbours_l(3),&
!!$               "Error in "//myself//": ineighbour > fd%num_neighbours_l")
!!$          owner = fd%neighbour_map_l(3)%m(ineighbour)
!!$        else
!!$          owner = -1
!!$        end if
!!$        request_nodes_top(halo_slab) = owner
!!$        request_slabs_top(halo_slab) = this_slab_global
!!$        !write(*,*) "DL_MG DEBUG: request L ", halo_slab, this_slab_global, owner
!!$      end do
!!$    end if
!!$
!!$
!!$    ! -----------------------------------
!!$    ! --- Issue the requests (irecvs) ---
!!$    ! -----------------------------------
!!$
!!$    ! jd: Issue irecvs for these, straight into the BOTTOM halos
!!$    do halo_slab = 1, num_pad_rows
!!$       if(request_nodes_bot(halo_slab) == -1) cycle
!!$       ! write(*,'(a,i0,a,i0,a,i0,a,i0)') 'Node ',myid,' asks hs ', &
!!$       !     halo_slab,' (',request_slabs_bot(halo_slab),') from #', &
!!$       !     request_nodes_bot(halo_slab)
!!$
!!$       ! TODO Is checking against the tag upper bound necessary?
!!$       ! JCW: Check if tag is greater than upper bound and, if so, abort
!!$!       if (request_slabs_bot(halo_slab) > dl_mg_defco_tag_ub) then
!!$!          write(*,'(a,i10)') "Error in "//myself//&
!!$!               ": tag > upper bound = ",dl_mg_defco_tag_ub
!!$!          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$!       end if
!!$
!!$       ! JCW: MPI_IRECV replaces ONETEP's comms_irecv routine.
!!$       ! comms_irecv is a wrapper of MPI_IRECV and in DL_MG we call MPI_IRECV
!!$       ! directly.
!!$       call MPI_IRECV(halo_bot(:,:,halo_slab),elements_in_slab,&
!!$            MPI_DOUBLE_PRECISION,request_nodes_bot(halo_slab),&
!!$            request_slabs_bot(halo_slab),fd%comm,&
!!$            request_handles_bot(halo_slab),ierr)
!!$       if (ierr /= MPI_SUCCESS) then
!!$          write(*,*) "Error in "//myself//": MPI_IRECV call failed"
!!$          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$       end if
!!$       !call comms_irecv(request_nodes_bot(halo_slab), &
!!$       !     halo_bot(:,:,halo_slab), elements_in_slab, &
!!$       !     tag=request_slabs_bot(halo_slab), &
!!$       !     handle=request_handles_bot(halo_slab))
!!$    end do
!!$
!!$    ! jd: Issue irecvs for these, straight into the TOP halos
!!$    !do halo_slab = -num_pad_rows, -1
!!$    do halo_slab = -1, -num_pad_rows, -1
!!$       if(request_nodes_top(halo_slab) == -1) cycle
!!$       ! write(*,'(a,i0,a,i0,a,i0,a,i0)') 'Node ',myid,' asks hs ', &
!!$       !     halo_slab,' (',request_slabs_top(halo_slab),') from #', &
!!$       !     request_nodes_top(halo_slab)
!!$
!!$       ! TODO Is checking against the tag upper bound necessary?
!!$       ! JCW: Check if tag is greater than upper bound and, if so, abort
!!$!       if (request_slabs_top(halo_slab) > dl_mg_defco_tag_ub) then
!!$!          write(*,'(a,i10)') "Error in "//myself//&
!!$!               ": tag > upper bound = ",dl_mg_defco_tag_ub
!!$!          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$!       end if
!!$
!!$       ! JCW: MPI_IRECV replaces ONETEP's comms_irecv routine.
!!$       ! comms_irecv is a wrapper of MPI_IRECV and in DL_MG we call MPI_IRECV
!!$       ! directly.
!!$       call MPI_IRECV(halo_top(:,:,halo_slab),elements_in_slab,&
!!$             MPI_DOUBLE_PRECISION,request_nodes_top(halo_slab),&
!!$             request_slabs_top(halo_slab),fd%comm,&
!!$             request_handles_top(halo_slab),ierr)
!!$       if (ierr /= MPI_SUCCESS) then
!!$          write(*,*) "Error in "//myself//": MPI_IRECV call failed"
!!$          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$       end if
!!$       !call comms_irecv(request_nodes_top(halo_slab), &
!!$       !     halo_top(:,:,halo_slab), elements_in_slab, &
!!$       !     tag=request_slabs_top(halo_slab), &
!!$       !     handle=request_handles_top(halo_slab))
!!$    end do
!!$
!!$    ! TODO Insert DL_MG timer end call (requests)
!!$
!!$    ! ----------------------------------------------------------
!!$    ! --- Figuring out what requests this node is to satisfy ---
!!$    ! ----------------------------------------------------------
!!$
!!$    ! TODO Insert DL_MG timer start call (requests to satisfy)
!!$
!!$    ! jd: Come up with a list of nodes that are going to need the slabs that
!!$    !     we have for their halos. Each slab may be needed by more than one
!!$    !     node (or none).
!!$    !do my_slab = 1, fd%mz
!!$    !   this_slab_global = my_first_slab_global + my_slab - 1
!!$    !   do node = 0, num_ranks-1
!!$    !      first_slab_on_node = fd_grid%first_slab12(node)
!!$    !      last_slab_on_node = fd_grid%last_slab12(node)
!!$
!!$    !      ! Does this guy need our slab for his TOP halo?
!!$    !      if(first_slab_on_node > this_slab_global .and. &
!!$    !           first_slab_on_node - num_pad_rows <= this_slab_global) then
!!$    !         satisfy_nodes(request_id) = node
!!$    !         satisfy_slabs(request_id) = this_slab_global
!!$    !         request_id = request_id + 1
!!$    !      end if
!!$
!!$    !      ! Does this guy need our slab for his BOTTOM halo?
!!$    !      if(last_slab_on_node < this_slab_global .and. &
!!$    !           last_slab_on_node + num_pad_rows >= this_slab_global) then
!!$    !         ! This guy needs our slab
!!$    !         satisfy_nodes(request_id) = node
!!$    !         satisfy_slabs(request_id) = this_slab_global
!!$    !         write(*,*) "DL_MG DEBUG: ", my_slab, this_slab_global, myid, node, &
!!$    !           last_slab_on_node, request_id
!!$    !         request_id = request_id + 1
!!$    !      end if
!!$
!!$    !   end do
!!$    !end do
!!$
!!$    request_id = 1
!!$    if (fd%num_neighbours_r(3) == 0) then
!!$      ! No left (TOP) slabs will be requested from this rank
!!$    else
!!$      ! We have >0 neighbours in right (BOTTOM) direction which will request halo
!!$      ! slabs for /their/ left (TOP) halos
!!$      first_slab_on_node = fd%ez + 1
!!$      do ineighbour = 1, fd%num_neighbours_r(3)
!!$        do my_slab = 1 + fd%mz - min(fd%mz,num_pad_rows), fd%mz
!!$          this_slab_global = fd%sz + my_slab - 1
!!$          if (first_slab_on_node - num_pad_rows <= this_slab_global) then
!!$             satisfy_nodes(request_id) = fd%neighbour_map_r(3)%m(ineighbour)
!!$             satisfy_slabs(request_id) = this_slab_global
!!$             !write(*,*) "DL_MG DEBUG: satisfy R", my_slab, this_slab_global, myid, fd%neighbour_map_r(3)%m(ineighbour), &
!!$             !  first_slab_on_node, request_id
!!$             request_id = request_id + 1
!!$          end if
!!$
!!$          call dl_mg_defco_utils_assert(request_id <= max_request_ids, &
!!$               'max_request_ids insufficient in '//myself)
!!$
!!$        end do
!!$        first_slab_on_node = first_slab_on_node + fd%halo_map_r(3)%m(ineighbour)
!!$        if ( first_slab_on_node - num_pad_rows > fd%ez ) then
!!$           ! Halo satisfied
!!$           !write(*,*) "DL_MG DEBUG: Halo satisfied", this_slab_global, request_id
!!$           exit
!!$        end if
!!$      end do
!!$    end if
!!$
!!$    if (fd%num_neighbours_l(3) == 0) then
!!$      ! No right (BOTTOM) slabs will be requested from this rank
!!$    else
!!$      ! We have >0 neighbours in left (TOP) direction which will request halo
!!$      ! slabs for /their/ right (BOTTOM) halos
!!$      last_slab_on_node = fd%sz - 1
!!$      do ineighbour = 1, fd%num_neighbours_l(3)
!!$        do my_slab = 1, min(fd%mz,num_pad_rows)
!!$          this_slab_global = fd%sz + my_slab - 1
!!$          if (last_slab_on_node + num_pad_rows >= this_slab_global) then
!!$             satisfy_nodes(request_id) = fd%neighbour_map_l(3)%m(ineighbour)
!!$             satisfy_slabs(request_id) = this_slab_global
!!$             !write(*,*) "DL_MG DEBUG: satisfy L", my_slab, this_slab_global, myid, fd%neighbour_map_l(3)%m(ineighbour), &
!!$             !  last_slab_on_node, request_id
!!$             request_id = request_id + 1
!!$          end if
!!$
!!$          call dl_mg_defco_utils_assert(request_id <= max_request_ids, &
!!$               'max_request_ids insufficient in '//myself)
!!$
!!$        end do
!!$        last_slab_on_node = last_slab_on_node - fd%halo_map_l(3)%m(ineighbour)
!!$        if ( last_slab_on_node + num_pad_rows < fd%sz) then
!!$           ! Halo satisfied
!!$           !write(*,*) "DL_MG DEBUG: Halo satisfied", this_slab_global, request_id
!!$           exit
!!$        end if
!!$      end do
!!$    end if
!!$
!!$    !if (fd%num_neighbours_l(3) == 0) then
!!$    !  ! No right (BOTTOM) slabs will be requested from this rank
!!$    !else
!!$    !  ! We have >0 neighbours in left (TOP) direction which will request halo
!!$    !  ! slabs for /their/ right (BOTTOM) halos
!!$    !  last_slab_on_node = fd%sz - 1
!!$    !  do ineighbour = 1, fd%num_neighbours_l(3)
!!$    !    do my_slab = 1, fd%mz
!!$    !      this_slab_global = fd%sz + my_slab - 1
!!$    !      if (last_slab_on_node + num_pad_rows >= this_slab_global) then
!!$    !    end do
!!$    !    last_slab_on_node = last_slab_on_node - fd%halo_map_l(3)%m(ineighbour)
!!$    !    if (fd%sz - last_slab_on_node >= num_pad_rows) then
!!$    !       ! Halo satisfied
!!$    !       write(*,*) "DL_MG DEBUG: Halo satisfied"
!!$    !       exit
!!$    !    end if
!!$    !  end do
!!$    !end if
!!$
!!$    !if (fd%num_neighbours_r(3) == 0) then
!!$    !  ! No top slabs will be requested from this rank
!!$    !else
!!$
!!$    !end if
!!$
!!$
!!$    ! ----------------------------------------------------
!!$    ! --- Satisfy requests issued to this node (sends) ---
!!$    ! ----------------------------------------------------
!!$
!!$    n_request_ids = request_id-1
!!$    do request_id = 1, n_request_ids
!!$       this_slab_local = satisfy_slabs(request_id) - my_first_slab_global+1
!!$       !write(*,'(a,i0,a,i0,a,i0)') 'Node ',myid,' sends (', &
!!$       !satisfy_slabs(request_id),') to #', satisfy_nodes(request_id)
!!$
!!$       ! TODO Is checking against the tag upper bound necessary?
!!$       ! JCW: Check if tag is greater than upper bound and, if so, abort
!!$!       if (satisfy_slabs(request_id) > dl_mg_defco_tag_ub) then
!!$!          write(*,'(a,i10)') "Error in "//myself//&
!!$!               ": tag > upper bound = ",dl_mg_defco_tag_ub
!!$!          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$!       end if
!!$
!!$       ! JCW: MPI_ISEND replaces ONETEP's comms_send routine.
!!$       ! comms_send is a wrapper of MPI_ISEND, but also makes checks ONETEP's
!!$       ! send_stack via comms_free_send_stack. DL_MG does not have a concept of
!!$       ! a send_stack, so we use MPI_ISEND directly.
!!$       !
!!$       ! IMPORTANT NOTE: It is important to avoid the creation of temporary array
!!$       !                 copies in this call. If fcn is non-contiguous, then the
!!$       !                 non-blocking MPI call may attempt to write to the dellocated
!!$       !                 address of the temporary array after returning. In practice
!!$       !                 this means ensuring that fcn is contiguous.
!!$       call MPI_ISEND(fcn(:,:,this_slab_local),elements_in_slab,&
!!$            MPI_DOUBLE_PRECISION,satisfy_nodes(request_id),satisfy_slabs(request_id),&
!!$            fd%comm,satisfy_handles(request_id),ierr)
!!$       if (ierr /= MPI_SUCCESS) then
!!$          write(*,*) "Error in "//myself//": MPI_ISEND call failed"
!!$          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$       end if
!!$       !call comms_send(satisfy_nodes(request_id), fcn(:,:,this_slab_local),&
!!$       !     elements_in_slab,tag=satisfy_slabs(request_id), &
!!$       !     return_handle=satisfy_handles(request_id), add_to_stack=.false.)
!!$    end do
!!$
!!$    ! TODO Insert DL_MG timer end call (requests to satisfy)
!!$
!!$    ! -----------------------------------------------------
!!$    ! --- Wait for completion of the receives and sends ---
!!$    ! -----------------------------------------------------
!!$
!!$    ! TODO Insert DL_MG timer start call (comms wait)
!!$
!!$    ! - receives for the BOTTOM halos
!!$    do halo_slab = 1, num_pad_rows
!!$       if(request_nodes_bot(halo_slab) == -1) cycle
!!$       ! JCW: MPI_WAIT replaces comms_wait from ONETEP's comms module
!!$       call MPI_WAIT(request_handles_bot(halo_slab),MPI_STATUS_IGNORE,ierr)
!!$       if (ierr /= MPI_SUCCESS) then
!!$          write(*,*) "Error in "//myself//": MPI_WAIT call failed"
!!$          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$       end if
!!$    end do
!!$
!!$    ! - receives for the TOP halos
!!$    do halo_slab = -num_pad_rows, -1
!!$       if(request_nodes_top(halo_slab) == -1) cycle
!!$       ! JCW: MPI_WAIT replaces comms_wait from ONETEP's comms module
!!$       call MPI_WAIT(request_handles_top(halo_slab),MPI_STATUS_IGNORE,ierr)
!!$       if (ierr /= MPI_SUCCESS) then
!!$          write(*,*) "Error in "//myself//": MPI_WAIT call failed"
!!$          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$       end if
!!$    end do
!!$
!!$    ! - sends (all)
!!$    do request_id = 1, n_request_ids
!!$       ! JCW: MPI_WAIT replaces comms_wait from ONETEP's comms module
!!$       call MPI_WAIT(satisfy_handles(request_id),MPI_STATUS_IGNORE,ierr)
!!$       if (ierr /= MPI_SUCCESS) then
!!$          write(*,*) "Error in "//myself//": MPI_WAIT call failed"
!!$          call MPI_ABORT(fd%comm,ierr,ierr1)
!!$       end if
!!$    end do
!!$
!!$    ! TODO Insert DL_MG timer end call (comms wait)
!!$
!!$    ! TODO Insert DL_MG timer start call (comms imbalance)
!!$    ! JCW: DL_MG routine which invokes mpi_barrier (was comms_barrier)
!!$    call barrier(fd%comm)
!!$    ! TODO Insert DL_MG timer end call (comms imbalance)
!!$
!!$#else
!!$    ! jd: Pointless, but kills the 'fcn unused' warning
!!$    ierr = int(fcn(1,1,1))
!!$#endif

    ! TODO Insert DL_MG timer end call

  end subroutine parallel_prepare_fd_halos


  subroutine parallel_destroy_fd_halos

    use dl_mg_defco_utils, only: dl_mg_defco_utils_dealloc_check

    implicit none

    ! jd: Internal variables
    integer :: ierr ! jd: Error flag

    !------------------------------------------------------------------------

    ! jd: Deallocate the halo arrays
    !     NB: This takes place regardless of whether we're in MPI or serial mode
    if(allocated(halo_left)) then
      deallocate(halo_left,stat=ierr)
      call dl_mg_defco_utils_dealloc_check(ierr,'parallel_destroy_fd_halos','halo_left')
    endif
    if(allocated(halo_right)) then
       deallocate(halo_right,stat=ierr)
       call dl_mg_defco_utils_dealloc_check(ierr,'parallel_destroy_fd_halos','halo_right')
     endif
  end subroutine parallel_destroy_fd_halos

end module dl_mg_defco_fd

