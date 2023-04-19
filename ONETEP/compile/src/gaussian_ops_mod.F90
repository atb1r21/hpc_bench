!================================================================!
!                                                                !
!                    Gaussian operations module                  !
!                                                                !
! This module contains subroutines for calculating analytical    !
! expressions involving Gaussians.                               !
!                                                                !
!----------------------------------------------------------------!
! Written by Jacek Dziedzic in October 2021.                     !
!================================================================!
! Reference:                                                     !
! [1] "The Computational Modelling of Heavy Atom Chemistry",     !
!     Chris-Kriton Skylaris (PhD thesis), Cambridge, 1999.       !
! url: www.southampton.ac.uk/assets/centresresearch/documents/   !
!      compchem/skylaris_phd.pdf                                 !
!================================================================!


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

module gaussian_ops

  use constants, only: DP, PI

  implicit none

  private

  public :: gaussian_overlap_integral

contains

  ! ---------------------------------------------------------------------------
  ! ----- P U B L I C   S U B R O U T I N E S   A N D   F U N C T I O N S -----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function gaussian_overlap_integral(alpha, beta, A, B, &
       ax, bx, ay, by, az, bz)
    !==========================================================================!
    ! Calculates (analytically) the overlap integral between two primitive     !
    ! Cartesian Gaussians, according to [1:(2.14)].                            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   alpha (in): Exponent of Gaussian A.                                    !
    !   beta (in):  Exponent of Gaussian B.                                    !
    !   A (in):     Position of Gaussian A.                                    !
    !   B (in):     Position of Gaussian B.                                    !
    !   ax, ay, az (in): Cartesian exponents of Gaussian A.                    !
    !   bx, by, bz (in): Cartesian exponents of Gaussian B.                    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2021.                               !
    !==========================================================================!

    use constants, only: SAFE_DIV_EPS
    use geometry, only: POINT
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: alpha
    real(kind=DP), intent(in) :: beta
    type(POINT), intent(in) :: A
    type(POINT), intent(in) :: B
    integer, intent(in) :: ax, bx
    integer, intent(in) :: ay, by
    integer, intent(in) :: az, bz

    ! jd: Local variables
    real(kind=DP) :: Lx, Ly, Lz
    real(kind=DP) :: kabx, kaby, kabz, kab
    real(kind=DP) :: Lambda
    character(len=*), parameter :: myself = 'gaussian_overlap_integral'

    ! -------------------------------------------------------------------------

    call utils_assert(alpha > 0.0_DP, myself//': alpha is not positive.')
    call utils_assert(beta > 0.0_DP, myself//': beta is not positive.')

    ! [1:(2.12)]
    Lambda = alpha + beta
    call utils_assert(Lambda > SAFE_DIV_EPS, myself//': Lambda too close to 0.')
    Lx = (alpha * A%x + beta * B%x) / Lambda
    Ly = (alpha * A%y + beta * B%y) / Lambda
    Lz = (alpha * A%z + beta * B%z) / Lambda
    kabx = exp(-(alpha*beta) * (A%x - B%x) * (A%x - B%x) / Lambda)
    kaby = exp(-(alpha*beta) * (A%y - B%y) * (A%y - B%y) / Lambda)
    kabz = exp(-(alpha*beta) * (A%z - B%z) * (A%z - B%z) / Lambda)
    kab = kabx * kaby * kabz

    ! [1:(2.14)]
    gaussian_overlap_integral = kab * &
         gaussian_1d_overlap(Lambda, Lx, A%x, B%x, ax, bx) * &
         gaussian_1d_overlap(Lambda, Ly, A%y, B%y, ay, by) * &
         gaussian_1d_overlap(Lambda, Lz, A%z, B%z, az, bz)

  end function gaussian_overlap_integral

  ! ---------------------------------------------------------------------------
  ! ---- P R I V A T E   S U B R O U T I N E S   A N D   F U N C T I O N S ----
  ! ---------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function gaussian_1d_overlap(Lambda, Lx, A, B, ax, bx)
    !==========================================================================!
    ! Calculates (analytically) the 1D overlap integral between two primitive  !
    ! Cartesian Gaussians, that is, the integral [1:(2.16)].                   !
    !                                                                          !
    ! Note that [1] does this integral via quadratures, while here we derive   !
    ! analytical expressions for ax in [0,6] and bx in [0,6] using Mathematica.!
    ! This can be done with the following Mathematica snippet:                 !
    !                                                                          !
    ! integrand[t_, \[CapitalLambda]_, Lx_, Ax_, Bx_, ax_, bx_] :=             !
    ! 1/Sqrt[\[CapitalLambda]] (t/Sqrt[\[CapitalLambda]] + Lx - Ax)^           !
    ! ax (t/Sqrt[\[CapitalLambda]] + Lx - Bx)^bx Exp[-t^2]                     !
    !                                                                          !
    ! integral[\[CapitalLambda]_, Lx_, Ax_, Bx_, ax_, bx_] := \!\(             !
    ! \*SubsuperscriptBox[\(\[Integral]\), \(-\[Infinity]\), \                 !
    ! \(\[Infinity]\)]\(integrand[t, \[CapitalLambda], \ Lx, Ax, Bx, ax,       !
    ! bx] \[DifferentialD]t\)\)                                                !
    !                                                                          !
    ! Table[integral[\[CapitalLambda], Lx, A, B, ax, bx], {ax, 0, 6}, {bx,     !
    !  0, 6}]//FortranForm                                                     !
    !                                                                          !
    ! The output can be pasted from Mathematica to a text file, and the only   !
    ! change required is the addition of "_DP" to constants. I added it to all !
    ! multiplicative constants, but not to integer exponents, as there is an   !
    ! expectation of a performance gain when they remain integers.             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   Lambda (in): Sum of Gaussian exponents, cf. [1:(2.12)].                !
    !   Lx (in): Gaussian midpoint, cf. [1:(2.12)].                            !
    !   A, B (in): Correspond to Ax and Bx in [1:(2.16)].                      !
    !   ax, bx (in): Cartesian exponents of Gaussians A and B.                 !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in October 2021.                               !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in) :: Lambda
    real(kind=DP), intent(in) :: Lx
    real(kind=DP), intent(in) :: A
    real(kind=DP), intent(in) :: B
    integer, intent(in) :: ax
    integer, intent(in) :: bx

    ! jd: Local variables
    character(len=*), parameter :: myself = 'gaussian_1d_overlap'

    ! -------------------------------------------------------------------------

    select case(ax)

    case(0)
       select case(bx)
       case(0)
          gaussian_1d_overlap = sqrt(PI)/sqrt(Lambda)
       case(1)
          gaussian_1d_overlap = ((-B + Lx)*sqrt(PI))/sqrt(Lambda)
       case(2)
          gaussian_1d_overlap = (sqrt(PI)*(1_DP + 2_DP*(B - Lx)**2*Lambda))/&
               (2_DP*Lambda**1.5_DP)
       case(3)
          gaussian_1d_overlap = ((B - Lx)*sqrt(PI)*(-3_DP - 2_DP*(B - Lx)**2&
              *Lambda))/(2_DP*Lambda**1.5_DP)
       case(4)
          gaussian_1d_overlap = (sqrt(PI)*(3_DP + 4_DP*(B - Lx)**2*Lambda*&
              (3_DP + (B - Lx)**2*Lambda)))/(4_DP*Lambda**2.5_DP)
       case(5)
          gaussian_1d_overlap = ((B - Lx)*sqrt(PI)*(-15_DP + 4_DP*(B - Lx)**2*&
              Lambda*(-5_DP - (B - Lx)**2*Lambda)))/(4_DP*Lambda**2.5_DP)
       case(6)
          gaussian_1d_overlap = (sqrt(PI)*(15_DP + 90_DP*(B - Lx)**2*Lambda + &
              60_DP*(B - Lx)**4*Lambda**2 + 8_DP*(B - Lx)**6*Lambda**3))/&
              (-(8._DP*Lambda**3.5_DP))
       case default
          call utils_abort(myself//"Unsupported bx: ", bx)
       end select

    case(1)
       select case(bx)
       case(0)
          gaussian_1d_overlap = ((-A + Lx)*sqrt(PI))/sqrt(Lambda)
       case(1)
          gaussian_1d_overlap = (sqrt(PI)*(1_DP + 2_DP*(A - Lx)*(B - Lx)*&
               Lambda))/(2_DP*Lambda**1.5_DP)
       case(2)
          gaussian_1d_overlap = (sqrt(PI)*(-A - 2_DP*B + 3_DP*Lx - 2_DP*&
               (A - Lx)*(B - Lx)**2*Lambda))/(2_DP*Lambda**1.5_DP)
       case(3)
          gaussian_1d_overlap = (sqrt(PI)*(3_DP + 6_DP*(A + B - 2_DP*Lx)*&
              (B - Lx)*Lambda + 4_DP*(A - Lx)*(B - Lx)**3*Lambda**2))/&
              (4_DP*Lambda**2.5_DP)
       case(4)
          gaussian_1d_overlap = (sqrt(PI)*(-3_DP*(A + 4_DP*B - 5_DP*Lx) - &
              4_DP*(3_DP*A + 2_DP*B - 5_DP*Lx)*(B - Lx)**2*Lambda - &
              4_DP*(A - Lx)*(B - Lx)**4*Lambda**2))/(4_DP*Lambda**2.5_DP)
       case(5)
          gaussian_1d_overlap = (sqrt(PI)*(15_DP + 30_DP*(A + 2_DP*B - 3_DP*Lx)*&
              (B - Lx)*Lambda + 20_DP*(2_DP*A + B - 3_DP*Lx)*(B - Lx)**3*&
              Lambda**2 + 8_DP*(A - Lx)*(B - Lx)**5*Lambda**3))/&
              (8_DP*Lambda**3.5_DP)
       case(6)
          gaussian_1d_overlap = (sqrt(PI)*(-15_DP*(A + 6_DP*B - 7_DP*Lx) - &
              30_DP*(3_DP*A + 4_DP*B - 7_DP*Lx)*(B - Lx)**2*Lambda - &
              12_DP*(5_DP*A + 2_DP*B - 7_DP*Lx)*(B - Lx)**4*Lambda**2 - 8_DP*&
              (A - Lx)*(B - Lx)**6*Lambda**3))/(8_DP*Lambda**3.5_DP)
       case default
          call utils_abort(myself//"Unsupported bx: ", bx)
       end select

    case(2)
       select case(bx)
       case(0)
          gaussian_1d_overlap = (sqrt(PI)*(1_DP + 2_DP*(A - Lx)**2*Lambda))/&
              (2_DP*Lambda**1.5_DP)
       case(1)
          gaussian_1d_overlap = (sqrt(PI)*(-2_DP*A - B + 3_DP*Lx + 2_DP*&
              (A - Lx)**2*(-B + Lx)*Lambda))/(2_DP*Lambda**1.5_DP)
       case(2)
          gaussian_1d_overlap = (sqrt(PI)*(3_DP + 2_DP*(A**2 + 4_DP*A*B + B**2&
              - 6_DP*(A + B)*Lx + 6_DP*Lx**2)*Lambda + 4_DP*(A - Lx)**2*&
              (B - Lx)**2*Lambda**2))/(4_DP*Lambda**2.5_DP)
       case(3)
          gaussian_1d_overlap = (sqrt(PI)*(-6_DP*A - 9_DP*B + 15_DP*Lx - &
              2_DP*(B - Lx)*(3_DP*A**2 + B**2 + 6_DP*A*(B - 2_DP*Lx) - &
              8_DP*B*Lx + 10_DP*Lx**2)*Lambda + 4_DP*(A - Lx)**2*(-B + Lx)**3*&
              Lambda**2))/(4_DP*Lambda**2.5_DP)
       case(4)
          gaussian_1d_overlap = (sqrt(PI)*(15_DP + 6_DP*(A**2 + 8_DP*A*B + &
              6_DP*B**2 - 10_DP*(A + 2_DP*B)*Lx + 15_DP*Lx**2)*Lambda + &
              4_DP*(B - Lx)**2*(6_DP*A**2 + 8_DP*A*B + B**2 - 10_DP*&
              (2_DP*A + B)*Lx + 15_DP*Lx**2)*Lambda**2 + 8_DP*(A - Lx)**2*&
              (B - Lx)**4*Lambda**3))/(8_DP*Lambda**3.5_DP)
       case(5)
          gaussian_1d_overlap = (sqrt(PI)*(-30_DP*A - 75_DP*B + 105_DP*Lx - &
              30_DP*(B - Lx)*(A**2 + 4_DP*A*B + 2_DP*B**2 - 6_DP*A*Lx - &
              8_DP*B*Lx + 7_DP*Lx**2)*Lambda - 4_DP*(B - Lx)**3*(10_DP*A**2 + &
              B**2 + 10_DP*A*(B - 3_DP*Lx) - 12_DP*B*Lx + 21_DP*Lx**2)*&
              Lambda**2 + 8_DP*(A - Lx)**2*(-B + Lx)**5*Lambda**3))/&
              (8_DP*Lambda**3.5_DP)
       case(6)
          gaussian_1d_overlap = (sqrt(PI)*(105_DP + 30_DP*(A**2 + 12_DP*A*B + &
               15_DP*B**2 - 14_DP*(A + 3_DP*B)*Lx + 28_DP*Lx**2)*Lambda + &
               60_DP*(B - Lx)**2*(3_DP*A**2 + 8_DP*A*B + 3_DP*B**2 - 14_DP*&
               (A + B)*Lx + 14_DP*Lx**2)*Lambda**2 + 8_DP*(B - Lx)**4*&
               (15_DP*A**2 + 12_DP*A*B + B**2 - 14_DP*(3_DP*A + B)*Lx + &
               28_DP*Lx**2)*Lambda**3 + 16_DP*(A - Lx)**2*(B - Lx)**6*&
               Lambda**4))/(16_DP*Lambda**4.5_DP)
       case default
          call utils_abort(myself//"Unsupported bx: ", bx)
       end select

    case(3)
       select case(bx)
       case(0)
          gaussian_1d_overlap = ((A - Lx)*sqrt(PI)*(-3_DP - 2_DP*(A - Lx)**2*&
               Lambda))/(2_DP*Lambda**1.5_DP)
       case(1)
          gaussian_1d_overlap = (sqrt(PI)*(3_DP + 6_DP*(A + B - 2_DP*Lx)*&
               (A - Lx)*Lambda + 4_DP*(A - Lx)**3*(B - Lx)*Lambda**2))/(4_DP*&
               Lambda**2.5_DP)
       case(2)
          gaussian_1d_overlap = (sqrt(PI)*(-9_DP*A - 6_DP*B + 15_DP*Lx - &
               2_DP*(A - Lx)*(A**2 + 6_DP*A*B + 3_DP*B**2 - 8_DP*A*Lx &
               - 12_DP*B*Lx + 10_DP*Lx**2)*Lambda - 4_DP*(A - Lx)**3* &
               (B - Lx)**2*Lambda**2))/(4_DP*Lambda**2.5_DP)
       case(3)
          gaussian_1d_overlap = (sqrt(PI)*(15_DP + 2_DP*Lambda*(9_DP*(A**2 + &
               3_DP*A*B + B**2 - 5_DP*(A + B)*Lx + 5_DP*Lx**2) + &
               6_DP*(A - Lx)*(B - Lx)*(A**2 + 3_DP*A*B + B**2 - 5_DP*(A + B)*&
               Lx + 5_DP*Lx**2)*Lambda + 4_DP*(A - Lx)**3*(B - Lx)**3*&
               Lambda**2)))/(8_DP*Lambda**3.5_DP)
       case(4)
          gaussian_1d_overlap = (sqrt(PI)*(-15_DP*(3_DP*A + 4_DP*B - 7_DP*Lx) -&
               6_DP*(A**3 + 12_DP*A**2*B + 18_DP*A*B**2 + 4_DP*B**3 - 15_DP*&
               (A**2 + 4_DP*A*B + 2_DP*B**2)*Lx + 15_DP*(3_DP*A + 4_DP*B)*Lx**2-&
               35_DP*Lx**3)*Lambda - 12_DP*(A - Lx)*(B - Lx)**2*(2_DP*A**2 + &
               B**2 + 4_DP*A*(B - 2_DP*Lx) - 6_DP*B*Lx + 7_DP*Lx**2)* &
               Lambda**2 - 8_DP*(A - Lx)**3*(B - Lx)**4*Lambda**3))/&
               (8_DP*Lambda**3.5_DP)
       case(5)
          gaussian_1d_overlap = (sqrt(PI)*(105_DP + 30_DP*(3_DP*A**2 + &
               15_DP*A*B + 10_DP*B**2 - 7_DP*(3_DP*A + 5_DP*B)*Lx + &
               28_DP*Lx**2)*Lambda + 60_DP*(A + B - 2_DP*Lx)*(B - Lx)*(A**2 + &
               5_DP*A*B + B**2 - 7_DP*(A + B)*Lx + 7_DP*Lx**2)*Lambda**2 + &
               8_DP*(A - Lx)*(B - Lx)**3*(10_DP*A**2 + 15_DP*A*B + 3_DP*B**2 - &
               7_DP*(5_DP*A + 3_DP*B)*Lx + 28_DP*Lx**2)*Lambda**3 - &
               16_DP*(A - Lx)**3*(-B + Lx)**5*Lambda**4))/(16_DP*Lambda**4.5_DP)
       case(6)
          gaussian_1d_overlap = (sqrt(PI)*(-315_DP*(A + 2_DP*B - 3_DP*Lx) - &
               30_DP*(A**3 + 18_DP*A**2*B + 45_DP*A*B**2 + 20_DP*B**3 - &
               21_DP*(A + B)*(A + 5_DP*B)*Lx + 84_DP*(A + 2_DP*B)*Lx**2 - &
               84_DP*Lx**3)*Lambda - 36_DP*(B - Lx)**2*(5_DP*A**3 + &
               20_DP*A**2*B + 15_DP*A*B**2 + 2_DP*B**3 - 7_DP*(5_DP*A**2 + &
               10_DP*A*B + 3_DP*B**2)*Lx + 14_DP*(5_DP*A + 4_DP*B)*Lx**2 - &
               42_DP*Lx**3)*Lambda**2 - 24_DP*(5_DP*A + B - 6_DP*Lx)*(A + B - &
               2_DP*Lx)*(A - Lx)*(B - Lx)**4*Lambda**3 - 16_DP*(A - Lx)**3*&
               (B - Lx)**6*Lambda**4))/(16_DP*Lambda**4.5_DP)
       case default
          call utils_abort(myself//"Unsupported bx: ", bx)
       end select

    case(4)
       select case(bx)
       case(0)
          gaussian_1d_overlap = (sqrt(PI)*(3_DP + 4_DP*(A - Lx)**2*Lambda* &
               (3_DP + (A - Lx)**2*Lambda)))/(4_DP*Lambda**2.5_DP)
       case(1)
          gaussian_1d_overlap = (sqrt(PI)*(-3_DP*(4_DP*A + B - 5_DP*Lx) - &
               4_DP*(2_DP*A + 3_DP*B - 5_DP*Lx)*(A - Lx)**2*Lambda + &
               4_DP*(A - Lx)**4*(-B + Lx)*Lambda**2))/(4_DP*Lambda**2.5_DP)
       case(2)
          gaussian_1d_overlap = (sqrt(PI)*(15_DP + 6_DP*(6_DP*A**2 + 8_DP*A*B +&
               B**2 - 10_DP*(2_DP*A + B)*Lx + 15_DP*Lx**2)*Lambda + &
               4_DP*(A - Lx)**2*(A**2 + 8_DP*A*B + 6_DP*B**2 - &
               10_DP*(A + 2_DP*B)*Lx + 15_DP*Lx**2)*Lambda**2 + &
               8_DP*(A - Lx)**4*(B - Lx)**2*Lambda**3))/(8_DP*Lambda**3.5_DP)
       case(3)
          gaussian_1d_overlap = (sqrt(PI)*(-15_DP*(4_DP*A + 3_DP*B - 7_DP*Lx) -&
               6_DP*(4_DP*A**3 + 18_DP*A**2*B + 12_DP*A*B**2 + B**3 - 15_DP*&
               (2_DP*A**2 + 4_DP*A*B + B**2)*Lx + 15_DP*(4_DP*A + 3_DP*B)*&
               Lx**2 - 35_DP*Lx**3)*Lambda + 12_DP*(A - Lx)**2*(-B + Lx)*&
               (A**2 + 4_DP*A*B + 2_DP*B**2 - 6_DP*A*Lx - 8_DP*B*Lx + &
               7_DP*Lx**2)*Lambda**2 + 8_DP*(A - Lx)**4*(-B + Lx)**3*&
               Lambda**3))/(8_DP*Lambda**3.5_DP)
       case(4)
          gaussian_1d_overlap = (sqrt(PI)*(105_DP + 4_DP*Lambda*(15_DP*(3_DP*&
               A**2 + 8_DP*A*B + 3_DP*B**2 - 14_DP*(A + B)*Lx + 14_DP*Lx**2) + &
               3_DP*(A**4 + 16_DP*A**3*B + 36_DP*A**2*B**2 + 16_DP*A*B**3 + &
               B**4 - 20_DP*(A + B)*(A**2 + 5_DP*A*B + B**2)*Lx + 30_DP*&
               (3_DP*A**2 + 8_DP*A*B + 3_DP*B**2)*Lx**2 - 140_DP*(A + B)*&
               Lx**3 + 70_DP*Lx**4)*Lambda + 4_DP*(A - Lx)**2*(B - Lx)**2*&
               (3_DP*A**2 + 8_DP*A*B + 3_DP*B**2 - 14_DP*(A + B)*Lx + &
               14_DP*Lx**2)*Lambda**2 + 4_DP*(A - Lx)**4*(B - Lx)**4*&
               Lambda**3)))/(16_DP*Lambda**4.5_DP)
       case(5)
          gaussian_1d_overlap = (sqrt(PI)*(-105_DP*(4_DP*A + 5_DP*B - 9_DP*Lx)-&
               60_DP*(2_DP*A**3 + 15_DP*A**2*B + 20_DP*A*B**2 + 5_DP*B**3 - &
               7_DP*(3_DP*A**2 + 10_DP*A*B + 5_DP*B**2)*Lx + 14_DP*(4_DP*A + &
               5_DP*B)*Lx**2 - 42_DP*Lx**3)*Lambda - 12_DP*(B - Lx)*(5_DP*A**4+&
               B**4 + 20_DP*A**3*(2_DP*B - 3_DP*Lx) - 24_DP*B**3*Lx + &
               126_DP*B**2*Lx**2 - 224_DP*B*Lx**3 + 126_DP*Lx**4 + &
               30_DP*A**2*(2_DP*B**2 - 8_DP*B*Lx + 7_DP*Lx**2) + 20_DP*A*&
               (B - 2_DP*Lx)*(B**2 - 7_DP*B*Lx + 7_DP*Lx**2))*Lambda**2 + &
               16_DP*(A - Lx)**2*(-B + Lx)**3*(5_DP*A**2 + 3_DP*B**2 + &
               10_DP*A*(B - 2_DP*Lx) - 16_DP*B*Lx + 18_DP*Lx**2)*Lambda**3 + &
               16_DP*(A - Lx)**4*(-B + Lx)**5*Lambda**4))/(16_DP*Lambda**4.5_DP)
       case(6)
          gaussian_1d_overlap = (sqrt(PI)*(945_DP + 630_DP*(2_DP*A**2 + &
               8_DP*A*B + 5_DP*B**2 - 6_DP*(2_DP*A + 3_DP*B)*Lx + &
               15_DP*Lx**2)*Lambda + 60_DP*(A**4 + 24_DP*A**3*B + &
               90_DP*A**2*B**2 + 80_DP*A*B**3 + 15_DP*B**4 - 28_DP*(A**3 + &
               9_DP*A**2*B + 15_DP*A*B**2 + 5_DP*B**3)*Lx + 84_DP*(2_DP*A**2 + &
               8_DP*A*B + 5_DP*B**2)*Lx**2 - 168_DP*(2_DP*A + 3_DP*B)*Lx**3 + &
               210_DP*Lx**4)*Lambda**2 + 24_DP*(B - Lx)**2*(15_DP*A**4 + &
               80_DP*A**3*B + 90_DP*A**2*B**2 + 24_DP*A*B**3 + B**4 - &
               28_DP*(5_DP*A**3 + 15_DP*A**2*B + 9_DP*A*B**2 + B**3)*Lx + &
               84_DP*(5_DP*A**2 + 8_DP*A*B + 2_DP*B**2)*Lx**2 - 168_DP*&
               (3_DP*A + 2_DP*B)*Lx**3 + 210_DP*Lx**4)*Lambda**3 + 48_DP*&
               (A - Lx)**2*(B - Lx)**4*(5_DP*A**2 + 8_DP*A*B + 2_DP*B**2 - &
               6_DP*(3_DP*A + 2_DP*B)*Lx + 15_DP*Lx**2)*Lambda**4 + &
               32_DP*(A - Lx)**4*(B - Lx)**6*Lambda**5))/(32_DP*Lambda**5.5_DP)
       case default
          call utils_abort(myself//"Unsupported bx: ", bx)
       end select

    case(5)
       select case(bx)
       case(0)
          gaussian_1d_overlap = ((A - Lx)*sqrt(PI)*(-15_DP + 4_DP*(A - Lx)**2*&
               Lambda*(-5_DP - (A - Lx)**2*Lambda)))/(4_DP*Lambda**2.5_DP)
       case(1)
          gaussian_1d_overlap = (sqrt(PI)*(15_DP + 30_DP*&
               (2_DP*A + B - 3_DP*Lx)*(A - Lx)*Lambda + 20_DP*&
               (A + 2_DP*B - 3_DP*Lx)*(A - Lx)**3*Lambda**2 + &
               8_DP*(A - Lx)**5*(B - Lx)*Lambda**3))/(8_DP*Lambda**3.5_DP)
       case(2)
          gaussian_1d_overlap = (sqrt(PI)*(-15_DP*(5_DP*A + 2_DP*B - 7_DP*Lx)- &
               30_DP*(A - Lx)*(2_DP*A**2 + B**2 + 4_DP*A*(B - 2_DP*Lx) - &
               6_DP*B*Lx + 7_DP*Lx**2)*Lambda - 4_DP*(A - Lx)**3*(A**2 + &
               10_DP*A*B + 10_DP*B**2 - 6_DP*(2_DP*A + 5_DP*B)*Lx + &
               21_DP*Lx**2)*Lambda**2 - 8_DP*(A - Lx)**5* &
               (B - Lx)**2*Lambda**3))/(8_DP*Lambda**3.5_DP)
       case(3)
          gaussian_1d_overlap = (sqrt(PI)*(105_DP + 30_DP*(10_DP*A**2 + &
               15_DP*A*B + 3_DP*B**2 - 7_DP*(5_DP*A + 3_DP*B)*Lx + &
               28_DP*Lx**2)*Lambda + 60_DP*(A + B - 2_DP*Lx)*(A - Lx)*(A**2 + &
               5_DP*A*B + B**2 - 7_DP*(A + B)*Lx + 7_DP*Lx**2)*Lambda**2 + &
               8_DP*(A - Lx)**3*(B - Lx)*(3_DP*A**2 + 15_DP*A*B + 10_DP*B**2 - &
               7_DP*(3_DP*A + 5_DP*B)*Lx + 28_DP*Lx**2)* Lambda**3 - &
               16_DP*(A - Lx)**5*(-B + Lx)**3*Lambda**4))/(16_DP*Lambda**4.5_DP)
       case(4)
          gaussian_1d_overlap = (sqrt(PI)*(-105_DP*(5_DP*A + 4_DP*B - 9_DP*Lx) - &
               60_DP*(5_DP*A**3 + 20_DP*A**2*B + 15_DP*A*B**2 + 2_DP*B**3 - &
               7_DP*(5_DP*A**2 + 10_DP*A*B + 3_DP*B**2)*Lx + 14_DP*(5_DP*A + &
               4_DP*B)*Lx**2 - 42_DP*Lx**3)*Lambda -12_DP*(A - Lx)*(A**4 + &
               20_DP*A**3*B + 60_DP*A**2*B**2 + 40_DP*A*B**3 + 5_DP*B**4 - &
               12_DP*(2_DP*A**3 + 15_DP*A**2*B + 20_DP*A*B**2 + 5_DP*B**3)*Lx +&
               42_DP*(3_DP*A**2 + 10_DP*A*B + 5_DP*B**2)*Lx**2 - &
               56_DP*(4_DP*A + 5_DP*B)*Lx**3 + 126_DP*Lx**4)*Lambda**2 - &
               16_DP*(A - Lx)**3*(B - Lx)**2*(3_DP*A**2 + 10_DP*A*B + &
               5_DP*B**2 - 4_DP*(4_DP*A + 5_DP*B)*Lx + 18_DP*Lx**2)*Lambda**3 -&
               16_DP*(A - Lx)**5*(B - Lx)**4*Lambda**4))/(16_DP*Lambda**4.5_DP)
       case(5)
          gaussian_1d_overlap = (sqrt(PI)*(945_DP + 9450_DP*Lx**2*Lambda + &
               12600_DP*Lx**4*Lambda**2 + 5040_DP*Lx**6*Lambda**3 + &
               720_DP*Lx**8*Lambda**4 + 32_DP*Lx**10*Lambda**5 - 8_DP*B**5*Lx*&
               Lambda**3*(15_DP + 20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) + &
               20_DP*B**4*Lambda**2*(15_DP + 90_DP*Lx**2*Lambda + 60_DP*Lx**4*&
               Lambda**2 + 8_DP*Lx**6*Lambda**3) - 40_DP*B**3*Lx*Lambda**2*&
               (105_DP + 210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) + 20_DP*B**2*Lambda*(105_DP + &
               840_DP*Lx**2*Lambda + 840_DP*Lx**4*Lambda**2 + &
               224_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4) - &
               10_DP*B*Lx*Lambda*(945_DP + 2520_DP*Lx**2*Lambda + &
               1512_DP*Lx**4*Lambda**2 + 288_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4) + 8_DP*A**5*(B - Lx)*Lambda**3*&
               (15_DP + 20_DP*Lx**2*Lambda + 4_DP*B**4*Lambda**2 - &
               16_DP*B**3*Lx*Lambda**2 + 4_DP*Lx**4*Lambda**2 - &
               8_DP*B*Lx*Lambda*(5_DP + 2_DP*Lx**2*Lambda) + 4_DP*B**2*Lambda*&
               (5_DP + 6_DP*Lx**2*Lambda)) - 20_DP*A**4*Lambda**2*(-15_DP - &
               90_DP*Lx**2*Lambda - 60_DP*Lx**4*Lambda**2 + &
               8_DP*B**5*Lx*Lambda**3 - 8_DP*Lx**6*Lambda**3 - &
               20_DP*B**4*Lambda**2*(1_DP + 2_DP*Lx**2*Lambda) + &
               40_DP*B**3*Lx*Lambda**2*(3_DP + 2_DP*Lx**2*Lambda) - &
               20_DP*B**2*Lambda*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) + 10_DP*B*Lx*Lambda*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2)) + &
               40_DP*A**3*Lambda**2*(4_DP*B**5*Lambda**2*(1_DP + &
               2_DP*Lx**2*Lambda) - 20_DP*B**4*Lx*Lambda**2*(3_DP + &
               2_DP*Lx**2*Lambda) + 20_DP*B**3*Lambda*(3_DP + &
               12_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               20_DP*B**2*Lx*Lambda*(15_DP + 20_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) + 5_DP*B*(15_DP + 90_DP*Lx**2*Lambda + &
               60_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) - &
               Lx*(105_DP + 210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3)) - 20_DP*A**2*Lambda*(-105_DP - &
               840_DP*Lx**2*Lambda - 840_DP*Lx**4*Lambda**2 - &
               224_DP*Lx**6*Lambda**3 - 16_DP*Lx**8*Lambda**4 + &
               8_DP*B**5*Lx*Lambda**3*(3_DP + 2_DP*Lx**2*Lambda) - &
               20_DP*B**4*Lambda**2*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) + 40_DP*B**3*Lx*Lambda**2*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - 20_DP*B**2*Lambda*&
               (15_DP + 90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) + 10_DP*B*Lx*Lambda*(105_DP + &
               210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3)) + 10_DP*A*Lambda*(4_DP*B**5*Lambda**2*&
               (3_DP + 12_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               20_DP*B**4*Lx*Lambda**2*(15_DP + 20_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) + 20_DP*B**3*Lambda*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) - 20_DP*B**2*Lx*Lambda*(105_DP + &
               210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) + 5_DP*B*(105_DP + 840_DP*Lx**2*Lambda + &
               840_DP*Lx**4*Lambda**2 + 224_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4) - Lx*(945_DP + 2520_DP*Lx**2*Lambda + &
               1512_DP*Lx**4*Lambda**2 + 288_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4))))/(32_DP*Lambda**5.5_DP)
       case(6)
          gaussian_1d_overlap = (sqrt(PI)*(8_DP*B**6*Lx*Lambda**3*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               24_DP*B**5*Lambda**2*(15_DP + 90_DP*Lx**2*Lambda + &
               60_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) + &
               60_DP*B**4*Lx*Lambda**2*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) - &
               40_DP*B**3*Lambda*(105_DP + 840_DP*Lx**2*Lambda + &
               840_DP*Lx**4*Lambda**2 + 224_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4) + 30_DP*B**2*Lx*Lambda*(945_DP + &
               2520_DP*Lx**2*Lambda + 1512_DP*Lx**4*Lambda**2 + &
               288_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4) - &
               6_DP*B*(945_DP + 9450_DP*Lx**2*Lambda + &
               12600_DP*Lx**4*Lambda**2 + 5040_DP*Lx**6*Lambda**3 + &
               720_DP*Lx**8*Lambda**4 + 32_DP*Lx**10*Lambda**5) + &
               Lx*(10395_DP + 34650_DP*Lx**2*Lambda + &
               27720_DP*Lx**4*Lambda**2 + 7920_DP*Lx**6*Lambda**3 + &
               880_DP*Lx**8*Lambda**4 + 32_DP*Lx**10*Lambda**5) - &
               4_DP*A**5*Lambda**2*(15_DP + 90_DP*Lx**2*Lambda + &
               60_DP*Lx**4*Lambda**2 + 8_DP*B**6*Lambda**3 - &
               48_DP*B**5*Lx*Lambda**3 + 8_DP*Lx**6*Lambda**3 + &
               60_DP*B**4*Lambda**2*(1_DP + 2_DP*Lx**2*Lambda) - &
               80_DP*B**3*Lx*Lambda**2*(3_DP + 2_DP*Lx**2*Lambda) + &
               30_DP*B**2*Lambda*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) - 12_DP*B*Lx*Lambda*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2)) + &
               20_DP*A**4*Lambda**2*(8_DP*B**6*Lx*Lambda**3 - &
               24_DP*B**5*Lambda**2*(1_DP + 2_DP*Lx**2*Lambda) + &
               60_DP*B**4*Lx*Lambda**2*(3_DP + 2_DP*Lx**2*Lambda) - &
               40_DP*B**3*Lambda*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) + 30_DP*B**2*Lx*Lambda*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - 6_DP*B*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) + Lx*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8*Lx**6*Lambda**3)) - &
               20*A**3*Lambda*(105_DP + 840_DP*Lx**2*Lambda + &
               840_DP*Lx**4*Lambda**2 + 224_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4 + 8_DP*B**6*Lambda**3*(1_DP + &
               2_DP*Lx**2*Lambda) - 48_DP*B**5*Lx*Lambda**3*(3_DP + &
               2_DP*Lx**2*Lambda) + 60_DP*B**4*Lambda**2*(3_DP + &
               12_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               80_DP*B**3*Lx*Lambda**2*(15_DP + 20_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) + 30_DP*B**2*Lambda*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) - 12_DP*B*Lx*Lambda*(105_DP + &
               210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3)) + 20_DP*A**2*Lambda*(8_DP*B**6*Lx*&
               Lambda**3*(3_DP + 2_DP*Lx**2*Lambda) - 24_DP*B**5*Lambda**2*&
               (3_DP + 12_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) + &
               60_DP*B**4*Lx*Lambda**2*(15_DP + 20_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) - 40_DP*B**3*Lambda*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) + 30_DP*B**2*Lx*Lambda*(105_DP + &
               210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) - 6_DP*B*(105_DP + 840_DP*Lx**2*Lambda + &
               840_DP*Lx**4*Lambda**2 + 224_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4) + Lx*(945_DP + 2520_DP*Lx**2*Lambda + &
               1512_DP*Lx**4*Lambda**2 + 288_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4)) - 5_DP*A*(945_DP + 9450_DP*Lx**2*Lambda+&
               12600_DP*Lx**4*Lambda**2 + 5040_DP*Lx**6*Lambda**3 + &
               720_DP*Lx**8*Lambda**4 + 32_DP*Lx**10*Lambda**5 + &
               8_DP*B**6*Lambda**3*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) - 48_DP*B**5*Lx*Lambda**3*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) + &
               60_DP*B**4*Lambda**2*(15_DP + 90_DP*Lx**2*Lambda + &
               60_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) - &
               80_DP*B**3*Lx*Lambda**2*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) + &
               30_DP*B**2*Lambda*(105_DP + 840_DP*Lx**2*Lambda + &
               840_DP*Lx**4*Lambda**2 + 224_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4) - 12_DP*B*Lx*Lambda*&
               (945_DP + 2520_DP*Lx**2*Lambda + 1512_DP*Lx**4*Lambda**2 + &
               288_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4))))/&
               (32_DP*Lambda**5.5_DP)
       case default
          call utils_abort(myself//"Unsupported bx: ", bx)
       end select

    case(6)
       select case(bx)
       case(0)
          gaussian_1d_overlap = (sqrt(PI)*(15_DP + 90_DP*(A - Lx)**2*Lambda + &
               60_DP*(A - Lx)**4*Lambda**2 + 8_DP*(A - Lx)**6*Lambda**3))/&
               (8_DP*Lambda**3.5_DP)
       case(1)
          gaussian_1d_overlap = (sqrt(PI)*(-15_DP*(6_DP*A + B - 7_DP*Lx) - &
               30_DP*(4_DP*A + 3_DP*B - 7_DP*Lx)*(A - Lx)**2*Lambda - &
               12_DP*(2_DP*A + 5_DP*B - 7_DP*Lx)*(A - Lx)**4*Lambda**2 + &
               8_DP*(A - Lx)**6*(-B + Lx)*Lambda**3))/(8_DP*Lambda**3.5_DP)
       case(2)
          gaussian_1d_overlap = (sqrt(PI)*(105_DP + 30_DP*(15_DP*A**2 + &
               12_DP*A*B + B**2 - 14_DP*(3_DP*A + B)*Lx + 28_DP*Lx**2)*Lambda +&
               60_DP*(A - Lx)**2*(3_DP*A**2 + 8_DP*A*B + 3_DP*B**2 - &
               14_DP*(A + B)*Lx + 14_DP*Lx**2)*Lambda**2 + &
               8_DP*(A - Lx)**4*(A**2 + 12_DP*A*B + 15_DP*B**2 - &
               14_DP*(A + 3_DP*B)*Lx + 28_DP*Lx**2)*Lambda**3 + &
               16_DP*(A - Lx)**6*(B - Lx)**2*Lambda**4))/(16_DP*Lambda**4.5_DP)
       case(3)
          gaussian_1d_overlap = (sqrt(PI)*(-315_DP*(2_DP*A + B - 3_DP*Lx) - &
               30_DP*(20_DP*A**3 + 45_DP*A**2*B + 18_DP*A*B**2 + B**3 - &
               21_DP*(A + B)*(5_DP*A + B)*Lx + 84_DP*(2_DP*A + B)*Lx**2 - &
               84_DP*Lx**3)*Lambda - 36_DP*(A - Lx)**2*(2_DP*A**3 + &
               15_DP*A**2*B + 20_DP*A*B**2 + 5_DP*B**3 - 7_DP*(3_DP*A**2 + &
               10_DP*A*B + 5_DP*B**2)*Lx + 14_DP*(4_DP*A + 5_DP*B)*Lx**2 - &
               42_DP*Lx**3)*Lambda**2 - 24_DP*(A + 5_DP*B - 6_DP*Lx)*(A + B - &
               2_DP*Lx)*(A - Lx)**4*(B - Lx)*Lambda**3 + &
               16_DP*(A - Lx)**6*(-B + Lx)**3*Lambda**4))/(16_DP*Lambda**4.5_DP)
       case(4)
          gaussian_1d_overlap = (sqrt(PI)*(945_DP + 630_DP*(5_DP*A**2 + &
               8_DP*A*B + 2_DP*B**2 - 6_DP*(3_DP*A + 2_DP*B)*Lx + &
               15_DP*Lx**2)*Lambda + 60_DP*(15_DP*A**4 + 80_DP*A**3*B + &
               90_DP*A**2*B**2 + 24_DP*A*B**3 + B**4 - 28_DP*(5_DP*A**3 + &
               15_DP*A**2*B + 9_DP*A*B**2 + B**3)*Lx + 84_DP*(5_DP*A**2 + &
               8_DP*A*B + 2_DP*B**2)*Lx**2 - 168_DP*(3_DP*A + 2_DP*B)*Lx**3 + &
               210_DP*Lx**4)*Lambda**2 + 24_DP*(A - Lx)**2*(A**4 + &
               24_DP*A**3*B + 90_DP*A**2*B**2 + 80_DP*A*B**3 + 15_DP*B**4 - &
               28_DP*(A**3 + 9_DP*A**2*B + 15_DP*A*B**2 + 5_DP*B**3)*Lx + &
               84_DP*(2_DP*A**2 + 8_DP*A*B + 5_DP*B**2)*Lx**2 - &
               168_DP*(2_DP*A + 3_DP*B)*Lx**3 + 210_DP*Lx**4)*Lambda**3 + &
               48_DP*(A - Lx)**4*(B - Lx)**2* (2_DP*A**2 + 8_DP*A*B + &
               5_DP*B**2 - 6_DP*(2_DP*A + 3_DP*B)*Lx + 15_DP*Lx**2)*Lambda**4 +&
               32_DP*(A - Lx)**6*(B - Lx)**4*Lambda**5))/(32_DP*Lambda**5.5_DP)
       case(5)
          gaussian_1d_overlap = (sqrt(PI)*(-4_DP*B**5*Lambda**2*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) + 20_DP*B**4*Lx*Lambda**2*(105_DP + &
               210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) - 20_DP*B**3*Lambda*(105_DP + &
               840_DP*Lx**2*Lambda + 840_DP*Lx**4*Lambda**2 + &
               224_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4) + &
               20_DP*B**2*Lx*Lambda*(945_DP + 2520_DP*Lx**2*Lambda + &
               1512_DP*Lx**4*Lambda**2 + 288_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4) - 5_DP*B*(945_DP + 9450_DP*Lx**2*Lambda +&
               12600_DP*Lx**4*Lambda**2 + 5040_DP*Lx**6*Lambda**3 + &
               720_DP*Lx**8*Lambda**4 + 32_DP*Lx**10*Lambda**5) + &
               Lx*(10395_DP + 34650_DP*Lx**2*Lambda + 27720_DP*Lx**4*Lambda**2+&
               7920_DP*Lx**6*Lambda**3 + 880_DP*Lx**8*Lambda**4 + &
               32_DP*Lx**10*Lambda**5) - 8_DP*A**6*(B - Lx)*Lambda**3* &
               (15_DP + 20_DP*Lx**2*Lambda + 4_DP*B**4*Lambda**2 - &
               16_DP*B**3*Lx*Lambda**2 + 4_DP*Lx**4*Lambda**2 - &
               8_DP*B*Lx*Lambda*(5_DP + 2_DP*Lx**2*Lambda) + &
               4_DP*B**2*Lambda*(5_DP + 6_DP*Lx**2*Lambda)) + &
               24_DP*A**5*Lambda**2*(-15_DP - 90_DP*Lx**2*Lambda - &
               60_DP*Lx**4*Lambda**2 + 8_DP*B**5*Lx*Lambda**3 - &
               8_DP*Lx**6*Lambda**3 - 20_DP*B**4*Lambda**2*(1_DP + &
               2_DP*Lx**2*Lambda) + 40_DP*B**3*Lx*Lambda**2*(3_DP + &
               2_DP*Lx**2*Lambda) - 20_DP*B**2*Lambda*(3_DP + &
               12_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) + &
               10_DP*B*Lx*Lambda*(15_DP + 20_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2)) - 60_DP*A**4*Lambda**2*(4_DP*B**5*&
               Lambda**2*(1_DP + 2_DP*Lx**2*Lambda) - &
               20_DP*B**4*Lx*Lambda**2*(3_DP + 2_DP*Lx**2*Lambda) + &
               20_DP*B**3*Lambda*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) - 20_DP*B**2*Lx*Lambda*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) + 5_DP*B*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) - Lx*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3)) + &
               40_DP*A**3*Lambda*(-105_DP - 840_DP*Lx**2*Lambda - &
               840_DP*Lx**4*Lambda**2 - 224_DP*Lx**6*Lambda**3 - &
               16_DP*Lx**8*Lambda**4 + 8_DP*B**5*Lx*Lambda**3*(3_DP + &
               2_DP*Lx**2*Lambda) - 20_DP*B**4*Lambda**2*(3_DP + &
               12_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) + &
               40_DP*B**3*Lx*Lambda**2*(15_DP + 20_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) - 20_DP*B**2*Lambda*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) + 10_DP*B*Lx*Lambda*(105_DP + &
               210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3)) - 30_DP*A**2*Lambda*&
               (4_DP*B**5*Lambda**2*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) - 20_DP*B**4*Lx*Lambda**2*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) + &
               20_DP*B**3*Lambda*(15_DP + 90_DP*Lx**2*Lambda + &
               60_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) - &
               20_DP*B**2*Lx*Lambda*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) + &
               5_DP*B*(105_DP + 840_DP*Lx**2*Lambda + 840_DP*Lx**4*Lambda**2 + &
               224_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4) - &
               Lx*(945_DP + 2520_DP*Lx**2*Lambda + 1512_DP*Lx**4*Lambda**2 + &
               288_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4)) + &
               6_DP*A*(-945_DP - 9450_DP*Lx**2*Lambda - &
               12600_DP*Lx**4*Lambda**2 - 5040_DP*Lx**6*Lambda**3 - &
               720_DP*Lx**8*Lambda**4 - 32_DP*Lx**10*Lambda**5 + 8_DP*B**5*Lx*&
               Lambda**3*(15_DP + 20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               20_DP*B**4*Lambda**2*(15_DP + 90_DP*Lx**2*Lambda + &
               60_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) + &
               40_DP*B**3*Lx*Lambda**2*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) - &
               20_DP*B**2*Lambda*(105_DP + 840_DP*Lx**2*Lambda + &
               840_DP*Lx**4*Lambda**2 + 224_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4) + 10_DP*B*Lx*Lambda*(945_DP + &
               2520_DP*Lx**2*Lambda + 1512_DP*Lx**4*Lambda**2 + &
               288_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4))))/&
               (32_DP*Lambda**5.5_DP)
       case(6)
          gaussian_1d_overlap = (sqrt(PI)*(10395_DP + 124740_DP*Lx**2*Lambda + &
               207900_DP*Lx**4*Lambda**2 + 110880_DP*Lx**6*Lambda**3 + &
               23760_DP*Lx**8*Lambda**4 + 2112_DP*Lx**10*Lambda**5 + &
               64_DP*Lx**12*Lambda**6 + 8_DP*B**6*Lambda**3*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) - 48_DP*B**5*Lx*Lambda**3*(105_DP + &
               210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) + 60_DP*B**4*Lambda**2*(105_DP + &
               840_DP*Lx**2*Lambda + 840_DP*Lx**4*Lambda**2 + &
               224_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4) - &
               80_DP*B**3*Lx*Lambda**2* (945_DP + 2520_DP*Lx**2*Lambda + &
               1512_DP*Lx**4*Lambda**2 + 288_DP*Lx**6*Lambda**3 +&
               16_DP*Lx**8*Lambda**4) + 30_DP*B**2*Lambda*(945_DP + &
               9450_DP*Lx**2*Lambda + 12600_DP*Lx**4*Lambda**2 + &
               5040_DP*Lx**6*Lambda**3 + 720_DP*Lx**8*Lambda**4 + &
               32_DP*Lx**10*Lambda**5) - 12_DP*B*Lx*Lambda*(10395_DP + &
               34650_DP*Lx**2*Lambda + 27720_DP*Lx**4*Lambda**2 + &
               7920_DP*Lx**6*Lambda**3 + 880_DP*Lx**8*Lambda**4 + &
               32_DP*Lx**10*Lambda**5) + 8_DP*A**6*Lambda**3*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*B**6*Lambda**3 - 48_DP*B**5*Lx*Lambda**3 + &
               8_DP*Lx**6*Lambda**3 + 60_DP*B**4*Lambda**2*(1_DP + &
               2_DP*Lx**2*Lambda) - 80_DP*B**3*Lx*Lambda**2*(3_DP + &
               2_DP*Lx**2*Lambda) + 30_DP*B**2*Lambda*(3_DP + &
               12_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               12_DP*B*Lx*Lambda*(15_DP + 20_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2)) - 48_DP*A**5*Lambda**3*(8_DP*B**6*Lx*&
               Lambda**3 - 24_DP*B**5*Lambda**2*(1_DP + 2_DP*Lx**2*Lambda) + &
               60_DP*B**4*Lx*Lambda**2*(3_DP + 2_DP*Lx**2*Lambda) - &
               40_DP*B**3*Lambda*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) + 30_DP*B**2*Lx*Lambda*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               6_DP*B*(15_DP + 90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) + Lx*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3)) + &
               60_DP*A**4*Lambda**2*(105_DP + 840_DP*Lx**2*Lambda + &
               840_DP*Lx**4*Lambda**2 + 224_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4 + 8_DP*B**6*Lambda**3*(1_DP + &
               2_DP*Lx**2*Lambda) - 48_DP*B**5*Lx*Lambda**3*(3_DP + &
               2_DP*Lx**2*Lambda) + 60_DP*B**4*Lambda**2*(3_DP + &
               12_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               80_DP*B**3*Lx*Lambda**2*(15_DP + 20_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) + 30_DP*B**2*Lambda*(15_DP + &
               90_DP*Lx**2*Lambda + 60_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3) - 12_DP*B*Lx*Lambda*(105_DP + &
               210_DP*Lx**2*Lambda + 84_DP*Lx**4*Lambda**2 + &
               8_DP*Lx**6*Lambda**3)) - 80_DP*A**3*Lambda**2*&
               (8_DP*B**6*Lx*Lambda**3*(3_DP + 2_DP*Lx**2*Lambda) - &
               24_DP*B**5*Lambda**2*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) + 60_DP*B**4*Lx*Lambda**2*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               40_DP*B**3*Lambda*(15_DP + 90_DP*Lx**2*Lambda + &
               60_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) + &
               30_DP*B**2*Lx*Lambda*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) - &
               6_DP*B*(105_DP + 840_DP*Lx**2*Lambda + 840_DP*Lx**4*Lambda**2 + &
               224_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4) + &
               Lx*(945_DP + 2520_DP*Lx**2*Lambda + 1512_DP*Lx**4*Lambda**2 + &
               288_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4)) + &
               30_DP*A**2*Lambda*(945_DP + 9450_DP*Lx**2*Lambda + &
               12600_DP*Lx**4*Lambda**2 + 5040_DP*Lx**6*Lambda**3 + &
               720_DP*Lx**8*Lambda**4 + 32_DP*Lx**10*Lambda**5 + &
               8_DP*B**6*Lambda**3*(3_DP + 12_DP*Lx**2*Lambda + &
               4_DP*Lx**4*Lambda**2) - 48_DP*B**5*Lx*Lambda**3*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) + &
               60_DP*B**4*Lambda**2*(15_DP + 90_DP*Lx**2*Lambda + &
               60_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) - &
               80_DP*B**3*Lx*Lambda**2*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) + &
               30_DP*B**2*Lambda*(105_DP + 840_DP*Lx**2*Lambda + &
               840_DP*Lx**4*Lambda**2 + 224_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4) - 12_DP*B*Lx*Lambda*(945_DP + &
               2520_DP*Lx**2*Lambda + 1512_DP*Lx**4*Lambda**2 + &
               288_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4)) - &
               12_DP*A*Lambda*(8_DP*B**6*Lx*Lambda**3*(15_DP + &
               20_DP*Lx**2*Lambda + 4_DP*Lx**4*Lambda**2) - &
               24_DP*B**5*Lambda**2*(15_DP + 90_DP*Lx**2*Lambda + &
               60_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) + &
               60_DP*B**4*Lx*Lambda**2*(105_DP + 210_DP*Lx**2*Lambda + &
               84_DP*Lx**4*Lambda**2 + 8_DP*Lx**6*Lambda**3) - &
               40_DP*B**3*Lambda*(105_DP + 840_DP*Lx**2*Lambda + &
               840_DP*Lx**4*Lambda**2 + 224_DP*Lx**6*Lambda**3 + &
               16_DP*Lx**8*Lambda**4) + 30_DP*B**2*Lx*Lambda* &
               (945_DP + 2520_DP*Lx**2*Lambda + 1512_DP*Lx**4*Lambda**2 + &
               288_DP*Lx**6*Lambda**3 + 16_DP*Lx**8*Lambda**4) - &
               6_DP*B*(945_DP + 9450_DP*Lx**2*Lambda + &
               12600_DP*Lx**4*Lambda**2 + 5040_DP*Lx**6*Lambda**3 + &
               720_DP*Lx**8*Lambda**4 + 32_DP*Lx**10*Lambda**5) + &
               Lx*(10395_DP + 34650_DP*Lx**2*Lambda + &
               27720_DP*Lx**4*Lambda**2 + 7920_DP*Lx**6*Lambda**3 + &
               880_DP*Lx**8*Lambda**4 + 32_DP*Lx**10*Lambda**5))))/&
               (64_DP*Lambda**6.5_DP)
       case default
          call utils_abort(myself//"Unsupported bx: ", bx)
       end select

    case default
       call utils_abort(myself//"Unsupported ax: ", ax)
    end select

end function gaussian_1d_overlap

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module gaussian_ops
