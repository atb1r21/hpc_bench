! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!==============================================================================!
!                                                                              !
!  ONETEP exchange-correlation functionals module: xc_funcs_mod.F90            !
!                                                                              !
!  The subroutines in this file were written by Peter Haynes, Quintin Hill,    !
!  Nicholas Hine, Jacek Dziedzic, David O'Regan, Lampros Andrinopoulos and     !
!  Gabriel Constantinescu between February 2004 and April 2014.                !
!                                                                              !
!  Addition of hybrid functionals, BLYP and XLYP, by Quintin Hill in           !
!  February/March 2009 under the supervision of Chris-Kriton Skylaris.         !
!                                                                              !
!==============================================================================!

module xc_funcs

  use constants, only: DP

  implicit none

  ! Local tolerance values for testing
  real(kind=DP), parameter :: dentol = 1.0e-15_DP
  ! JCW: tautol is set in rundat, as pub_xc_min_tau, allowing users to
  ! JCW: change the tolerance, if necessary
  !real(kind=DP), parameter :: tautol = 1.0e-10_DP
  real(kind=DP), parameter :: THIRD = 1.0_DP / 3.0_DP
  real(kind=DP), parameter :: FTHRD = 4.0_DP / 3.0_DP

  private

  ! Public subroutines
  public :: xc_lda_x_point
  public :: xc_lsda_x_point
  public :: xc_capz_c_point
  public :: xc_capz_c_point_sp
  public :: xc_vwn_c_point
  public :: xc_vwn_c_point_sp
  public :: xc_pw92_c_point
  public :: xc_pw92_c_point_sp
  public :: xc_pbe_x_point
  public :: xc_rpbe_x_point
  public :: xc_revpbe_x_point
  public :: xc_pbesol_x_point
  public :: xc_wc_x_point
  public :: xc_pbe0_x_point
  public :: xc_pbe_c_point
  public :: xc_pbesol_c_point
  public :: xc_pbe_x_point_sp
  public :: xc_rpbe_x_point_sp
  public :: xc_revpbe_x_point_sp
  public :: xc_pbesol_x_point_sp
  public :: xc_wc_x_point_sp
  public :: xc_pbe0_x_point_sp
  public :: xc_pbe_c_point_sp
  public :: xc_pbesol_c_point_sp
  public :: xc_pw91_x_point
  public :: xc_pw91_c_point
  public :: xc_pw91_x_point_sp
  public :: xc_pw91_c_point_sp
  public :: xc_b88_x_point
  public :: xc_b88_x_point_sp
  public :: xc_b1_x_point
  public :: xc_b3_x_point
  public :: xc_x_x_point
  public :: xc_x3_x_point
  public :: xc_b1_x_point_sp
  public :: xc_b3_x_point_sp
  public :: xc_x_x_point_sp
  public :: xc_x3_x_point_sp
  public :: xc_lyp_c_point
  public :: xc_lyp_c_point_sp
  public :: xc_b3lyp_c_point
  public :: xc_b3lyp_c_point_sp
  public :: xc_b3pw91_c_point
  public :: xc_b3pw91_c_point_sp
  public :: xc_x3lyp_c_point
  public :: xc_x3lyp_c_point_sp
  public :: xc_none_c_point
  public :: xc_none_c_point_sp
  public :: xc_none_x_point
  public :: xc_none_x_point_sp
  public :: xc_pw91_eq10
  public :: xc_fxc_capz_point
  public :: xc_rpw86_x_point
  public :: xc_rpw86_x_point_sp
  public :: xc_c09_x_point
  public :: xc_c09_x_point_sp
  public :: xc_vdwoptb88_x_point
  public :: xc_vdwoptb88_x_point_sp
  public :: xc_vdwpbek1_x_point
  public :: xc_vdwpbek1_x_point_sp
  public :: xc_vdwoptpbe_x_point
  public :: xc_vdwoptpbe_x_point_sp
  public :: xc_am05_x_point ! gcc32
  public :: xc_am05_c_point ! gcc32
  public :: xc_am05_x_point_sp ! JA
  public :: xc_am05_c_point_sp ! JA
  !
  public :: xc_none_mgga_x_point
  public :: xc_none_mgga_c_point
  public :: xc_none_mgga_x_point_sp
  public :: xc_none_mgga_c_point_sp
  public :: xc_lta_x_point
  public :: xc_lta_x_point_sp
  public :: xc_pkzb_x_point
  public :: xc_pkzb_c_point
  public :: xc_pkzb_x_point_sp
  public :: xc_pkzb_c_point_sp
  !
  public :: xc_b97mv_x_point
  public :: xc_b97mv_c_point
  public :: xc_b97mv_x_point_sp
  public :: xc_b97mv_c_point_sp

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines for evaluating the energy and potential contibution of each
! point.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lda_x_point(den,x_energy,x_pot)

    !==========================================================================!
    ! This subroutine calculates the local density appoximation (LDA) exchange !
    ! energy and potential at a point given the density at that point.         !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in March 2009 based on existing xc_capz          !
    ! subroutine written by Peter Haynes in February 2004.                     !
    !==========================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den      ! Density
    real(kind=DP), intent(out) :: x_energy ! Exchange energy
    real(kind=DP), intent(out) :: x_pot    ! Exchange potential

    real(kind=DP) :: rs
    real(kind=DP) :: epsx
    real(kind=DP), parameter :: xf = -0.458165293283142893475554_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    if (den > 0.0_DP) then   ! Positive charge density

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       epsx = xf / rs
       x_energy = den * epsx
       x_pot = FTHRD * epsx
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
    end if

  end subroutine xc_lda_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lsda_x_point(den1,den2,x_energy,x_pot1,x_pot2)

    !==========================================================================!
    ! This subroutine calculates the local spin density appoximation (LSDA)    !
    ! exchange energy and potential at a point given the density at that point !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in March 2009 based on existing xc_capz          !
    ! subroutine written by Peter Haynes in February 2004.                     !
    !==========================================================================!

    use constants, only: DP
    implicit none

    real(kind=DP), intent(in)  :: den1
    real(kind=DP), intent(in)  :: den2
    real(kind=DP), intent(out) :: x_energy
    real(kind=DP), intent(out) :: x_pot1
    real(kind=DP), intent(out) :: x_pot2

    real(kind=DP) :: epsx1
    real(kind=DP) :: epsx2
    real(kind=DP), parameter :: xfp = -0.930525736349100025002010_DP

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    ! Check for Positive charge densities
    if (den1 >= 0.0_DP) then
       epsx1 = xfp * den1**THIRD
       x_energy = den1 * epsx1
       x_pot1 = FTHRD * epsx1
    end if
    if (den2 >= 0.0_DP) then
       epsx2 = xfp * den2**THIRD
       x_energy = x_energy + den2 * epsx2
       x_pot2 = FTHRD * epsx2
    end if

  end subroutine xc_lsda_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_capz_c_point(den,c_energy,c_pot)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Zunger.                                                          !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew & A Zunger, Phys. Rev. B 23, 5048 (1981)                      !
    !---------------------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.                        !
    ! Split into point contribution by Quintin Hill on 09/03/2009.              !
    !===========================================================================!


    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot

    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: logrs        ! Natural logarithm of rs
    real(kind=DP) :: epsc         ! Correlation energy per particle

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! Constants for spin-unpolarised case:
    real(kind=DP), parameter :: g = -0.1423_DP
    real(kind=DP), parameter :: b1 = 1.0529_DP
    real(kind=DP), parameter :: b2 = 0.3334_DP
    real(kind=DP), parameter :: aa = 0.0311_DP
    real(kind=DP), parameter :: bb = -0.048_DP
    real(kind=DP), parameter :: cc = 0.0020_DP
    real(kind=DP), parameter :: dd = -0.0116_DP

    if (den > 0.0_DP) then   ! Positive charge density

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius

       if (rs >= 1.0_DP) then  ! Low charge density

          sqrtrs = sqrt(rs)
          epsc = g / (1.0_DP + b1 * sqrtrs + b2 * rs)
          c_energy = den * epsc
          c_pot = epsc * epsc * (1.0_DP + SSXTH * b1 * sqrtrs + &
               FTHRD * b2 * rs) / g

       else  ! High charge density

          logrs = log(rs)
          epsc = aa * logrs + bb + &
               (cc * logrs + dd) * rs
          c_energy = den * epsc
          c_pot = epsc - THIRD * (aa + &
               (cc * logrs + cc + dd) * rs)

       end if

    else   ! Negative charge density
       c_energy = 0.0_DP
       c_pot = 0.0_DP

    end if

  end subroutine xc_capz_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_capz_c_point_sp(den,den1,den2,c_energy,c_pot1,c_pot2)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Zunger.                                                          !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew & A Zunger, Phys. Rev. B 23, 5048 (1981)                      !
    !---------------------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.                        !
    ! Split into point contribution by Quintin Hill on 09/03/2009.              !
    !===========================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(in)  :: den1
    real(kind=DP), intent(in)  :: den2
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot1
    real(kind=DP), intent(out) :: c_pot2

    real(kind=DP) :: zeta         ! Spin polarisation
    real(kind=DP) :: fzeta, dfdz  ! Function of zeta and its derivative
    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: logrs        ! Natural logarithm of rs
    real(kind=DP) :: epsc         ! Correlation energy per particle
    real(kind=DP) :: epscp        ! Polarised correlation energy per particle
    real(kind=DP) :: vcu,vcp   ! Unpolarised and polarised correlation potentials

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! Constants for spin-unpolarised case:
    real(kind=DP), parameter :: g = -0.1423_DP
    real(kind=DP), parameter :: b1 = 1.0529_DP
    real(kind=DP), parameter :: b2 = 0.3334_DP
    real(kind=DP), parameter :: aca = 0.0311_DP
    real(kind=DP), parameter :: bca = -0.048_DP
    real(kind=DP), parameter :: cca = 0.0020_DP
    real(kind=DP), parameter :: dca = -0.0116_DP

    ! Constants for spin-polarised case:
    real(kind=DP), parameter :: gp = -0.0843_DP
    real(kind=DP), parameter :: b1p = 1.3981_DP
    real(kind=DP), parameter :: b2p = 0.2611_DP
    real(kind=DP), parameter :: acap = 0.01555_DP
    real(kind=DP), parameter :: bcap = -0.0269_DP
    real(kind=DP), parameter :: ccap = 0.0007_DP
    real(kind=DP), parameter :: dcap = -0.0048_DP

    ! Positive charge densities
    if (den > 0.0_DP) then

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius

       if (rs >= 1.0_DP) then    ! Low charge density

          sqrtrs = sqrt(rs)
          epsc = g / (1.0_DP + b1 * sqrtrs + b2 * rs)
          vcu = epsc * epsc * (1.0_DP + SSXTH * b1 * sqrtrs + &
               FTHRD * b2 * rs) / g
          epscp = gp / (1.0_DP + b1p * sqrtrs + b2p * rs)
          vcp = epscp * epscp * (1.0_DP + SSXTH * b1p * sqrtrs + &
               FTHRD * b2p * rs) / gp

       else   ! High charge density

          logrs = log(rs)
          epsc = aca * logrs + bca + (cca * logrs + dca) * rs
          vcu = epsc - &
               (aca + (cca * logrs + cca + dca) * rs) * THIRD
          epscp = acap * logrs + bcap + (ccap * logrs + dcap) * rs
          vcp = epscp - &
               (acap + (ccap * logrs + ccap + dcap) * rs) * THIRD

       end if

       zeta = (den1 - den2) / den   ! Polarisation
       !pdh: check that -1 <= zeta <= 1
       zeta = min(zeta,1.0_DP)
       zeta = max(zeta,-1.0_DP)

       fzeta = ((1.0_DP + zeta)**FTHRD + &
            (1.0_DP - zeta)**FTHRD - 2.0_DP) / &
            (2.0_DP**FTHRD - 2.0_DP)
       dfdz = FTHRD * ((1.0_DP + zeta)**THIRD - &
            (1.0_DP - zeta)**THIRD) / &
            (2.0_DP**FTHRD - 2.0_DP)

       c_pot1 = vcu + fzeta * (vcp - vcu) + &
            (epscp - epsc) * (1.0_DP - zeta) * dfdz
       c_pot2 = vcu + fzeta * (vcp - vcu) - &
            (epscp - epsc) * (1.0_DP + zeta) * dfdz

       c_energy = den * (epsc + fzeta * (epscp - epsc))

    else  ! Negative charge density
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
    end if

  end subroutine xc_capz_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_fxc_capz_point(den, pot1, pot2)
    !===========================================================================!
    ! Subroutine calculating fxc in the local density approximation from the    !
    ! Monte-Carlo data by Ceperley and Alder, Perdew Zunger parameterisation at !
    ! a single density value den.                                               !
    !===========================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in) :: den
    real(kind=DP), intent(inout) :: pot1
    real(kind=DP), intent(inout) :: pot2

    ! internal variables
    real(kind=DP) :: rs  !Wigner-Seitz radius
    real(kind=DP) :: sqrtrs  ! sqrt(rs)
    real(kind=DP) :: logrs  ! log(rs)
    real(kind=DP) :: drs  ! d(rs)/d(rho)
    real(kind=DP) :: f_x ! exchange part
    real(kind=DP) :: f_c ! correlation part
    real(kind=DP) :: num
    real(kind=DP) :: denom

    ! constants
    real(kind=DP), parameter :: TONFPI=3.0_DP/(4.0_DP*PI)
    real(kind=DP), parameter :: MIN_DENS = 5.0D-13

    ! Constants for spin-unpolarised case:
    real(kind=DP), parameter :: g = -0.1423_DP
    real(kind=DP), parameter :: b1 = 1.0529_DP
    real(kind=DP), parameter :: b2 = 0.3334_DP
    real(kind=DP), parameter :: aa = 0.0311_DP
!   real(kind=DP), parameter :: bb = -0.048_DP  ! fxc doesn't depend on bb
    real(kind=DP), parameter :: cc = 0.0020_DP
    real(kind=DP), parameter :: dd = -0.0116_DP

    ! exchange factor
    real(kind=DP), parameter :: xf = -0.458165293283142893475554_DP

    if (den>MIN_DENS) then
       ! Define parameters
       rs=(TONFPI/den)**THIRD
       drs=-THIRD*rs/den

       ! calculate exchange part. Diverges for small rho...
       f_x=THIRD*FTHRD*xf/(rs*den)

       ! calculate correlation part. distinguish between high and low density
       if(rs>=1.0_DP) then ! low charge density.

         sqrtrs=sqrt(rs)
         denom=12.0_DP*sqrtrs*(b1*sqrtrs+b2*rs+1.0_DP)**3.0_DP
         num=-7.0_DP*b1*b1*sqrtrs-b1*(21.0_DP*b2*rs+5.0_DP)
         num=num-8.0_DP*b2*sqrtrs*(2.0_DP*b2*rs+1.0_DP)
         f_c=g*drs*num/denom

       else ! high charge density

         logrs=log(rs)
         f_c=(aa/rs+THIRD*(2.0_DP*cc*logrs+2.0_DP*dd+cc))*drs

       endif

    else
       f_x=0.0_DP  ! For negative charge density, set to 0
       f_c=0.0_DP  ! For negative charge density, set to 0
    endif

    pot1=f_x+f_c
    pot2=f_c

  end subroutine xc_fxc_capz_point



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vwn_c_point(den,c_energy,c_pot)
    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the correlation energy and potential, according to the local density      !
    ! approximation (LDA), from the parameterisation by Vosko, Wilk and Nusair. !
    !                                                                           !
    ! Spin-unpolarized version.                                                 !
    !                                                                           !
    ! Reference:                                                                !
    !  [1] S H Vosko, L Wilk & M Nusair, Can. J. Phys. 58, 1200 (1980)          !
    !  [2] S H Vosko and L Wilk, Phys. Rev. B 22, 3812 (1980)                   !
    !                                                                           !
    !---------------------------------------------------------------------------!
    ! Written in 2011/05/31 by Jacek Dziedzic, basing on a template by Peter    !
    ! Haynes, modified by Quintin Hill.                                         !
    !===========================================================================!
    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot

    real(kind=DP) :: rs ! Wigner-Seitz radius corresponding to den

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: minus_fourpi_over_nine = -4.0_DP/9.0_DP * PI
    real(kind=DP), parameter :: MIN_DENS = 5.0D-13

    ! Local variables
    real(kind=DP) :: c_eps     ! Correlation energy density
    real(kind=DP) :: drs_drho  ! Derivative of rs wrt electronic density

    ! -------------------------------------------------------------------------

    ! We're only concerned about positive densities.
    ! Furthermore, to get perfect agreement with LIBXC, we ignore densities
    ! below MIN_DENS just like they do
    if (den > MIN_DENS) then

       rs = (TONFPI / den)**THIRD ! the Wigner-Seitz radius

       ! Calculate energy and potential (eq. (4.4) in [1] and its derivative)
       call xc_vwn_eps_c_helper(c_eps,c_pot,rs,1)

       ! Energy density -> energy
       c_energy = c_eps * den

       ! drs/drho
       drs_drho = minus_fourpi_over_nine * rs**4.0_DP

       ! Resulting potential (functional derivative of E_c (not eps_c) wrt rho)
       ! (v_c = eps_c + rho * de_c/drho), where de_c/drho = de_c/drs * drs/drho
       c_pot = c_eps + den * c_pot * drs_drho

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
    end if

  end subroutine xc_vwn_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vwn_c_point_sp(den,den1,den2,c_energy,c_pot1,c_pot2)
    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the correlation energy and potential, according to the local density      !
    ! approximation (LDA), from the parameterisation by Vosko, Wilk and Nusair. !
    !                                                                           !
    ! Spin-polarized version.                                                   !
    !                                                                           !
    ! Reference:                                                                !
    !  [1] S H Vosko, L Wilk & M Nusair, Can. J. Phys. 58, 1200 (1980)          !
    !  [2] S H Vosko and L Wilk, Phys. Rev. B 22, 3812 (1980)                   !
    !                                                                           !
    !---------------------------------------------------------------------------!
    ! Written in 2011/05/31 by Jacek Dziedzic, basing on a template by Peter    !
    ! Haynes, modified by Quintin Hill.                                         !
    !===========================================================================!

    use constants, only: DP, PI, SAFE_DIV_EPS
    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(in)  :: den, den1, den2
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot1, c_pot2

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: d2f_at_zero = &
         1.7099209341613656176_DP ! f''(0) = 4/(9 * (-1 + 2**(1/3))), cf. [1]
    real(kind=DP), parameter :: minus_fourpi_over_nine = -4.0_DP/9.0_DP * PI
    real(kind=DP), parameter :: MIN_DENS = 5.0D-13

    ! Local variables
    real(kind=DP) :: c_eps     ! Correlation energy density
    real(kind=DP) :: rs        ! Wigner-Seitz radius
    real(kind=DP) :: alpha_c   ! See discussion below (9) in [2]
    real(kind=DP) :: eps_f     ! Ferromagnetic contribution to eps
    real(kind=DP) :: eps_p     ! Paramagnetic contribution to eps
    real(kind=DP) :: pot_c, pot_f, pot_p ! Corresponding 1st derivatives
    real(kind=DP) :: eps_delta ! (7) in [2]
    real(kind=DP) :: pot_delta ! (7) in [2] with epsilons replaced by derivs
    real(kind=DP) :: beta      ! (8) in [2]
    real(kind=DP) :: beta_pot  ! (8) in [2] with epsilons replaced by derivs
    real(kind=DP) :: f_of_zeta ! (4) in [2]
    real(kind=DP) :: fprime_of_zeta ! d/dzeta of (4) in [2]
    real(kind=DP) :: zeta      ! Polarization
    real(kind=DP) :: zeta3     ! zeta**3
    real(kind=DP) :: zeta4     ! zeta**4
    real(kind=DP) :: de_drs    ! deps/drs
    real(kind=DP) :: de_dzeta, de_dzeta1, de_dzeta2 ! deps/dzeta and its terms
    real(kind=DP) :: c_pot     ! Resulting potential, before zeta correction
    real(kind=DP) :: drs_drho  ! Derivative of rs wrt electronic density

    ! -------------------------------------------------------------------------

    ! We're only concerned about positive densities.
    ! Furthermore, to get perfect agreement with LIBXC, we ignore densities
    ! below MIN_DENS just like they do
    if (den > MIN_DENS) then

       rs = (TONFPI / den)**THIRD ! the Wigner-Seitz radius

       ! Calculate spin-stiffness contribution, alpha_c(r_s)
       ! see discussion below (9) in [2].
       call xc_vwn_eps_c_helper(alpha_c,pot_c,rs,3)

       ! Calculate paramagnetic contribution eps_p
       ! see discussion below (9) in [2].
       call xc_vwn_eps_c_helper(eps_p,pot_p,rs,1)

       ! Calculate ferromagnetic contribution eps_f
       ! see discussion below (9) in [2].
       call xc_vwn_eps_c_helper(eps_f,pot_f,rs,2)

       ! This should never happen. Since this is executed many times, an
       ! explicit check is better than an assert.
       if(abs(alpha_c) < SAFE_DIV_EPS) then
          call utils_abort("Division by (almost) zero in xc_vwn_c_point_sp", &
               opt_real_to_print1=alpha_c)
       end if

       ! (8) in [2]
       beta = d2f_at_zero * (eps_f - eps_p) / alpha_c - 1.0_DP
       beta_pot = d2f_at_zero * (pot_f - pot_p) / pot_c - 1.0_DP

       ! Polarization
       zeta = (den1 - den2) / den
       zeta = min(zeta,1.0_DP)
       zeta = max(zeta,-1.0_DP)

       ! zeta**3 and zeta**4
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta*zeta3

       ! (4) in [2]
       f_of_zeta = 0.5_DP * ((1+zeta)**(4.0_DP/3.0_DP) + &
            (1-zeta)**(4.0_DP/3.0_DP) - 2.0_DP) / (2.0_DP**THIRD - 1.0_DP)

       ! 1st derivative of the above
       fprime_of_zeta = ( -4.0_DP/3.0_DP*(1.0_DP-zeta)**THIRD + &
            4.0_DP/3.0_DP*(1.0_DP+zeta)**THIRD ) / &
            (2.0_DP * (-1.0_DP + 2.0_DP**THIRD))

       ! Corrections to energy density and potential due to polarization:
       ! - delta_eps_c
       eps_delta = alpha_c * f_of_zeta/d2f_at_zero * (1.0_DP + beta * zeta4)

       ! - delta_v_c
       pot_delta = pot_c * f_of_zeta/d2f_at_zero * (1.0_DP + beta_pot * zeta4)

       ! Corrected energy density: eps_paramagnetic + delta_eps
       c_eps = eps_p + eps_delta

       ! Corrected derivative of c_eps wrt rs
       de_drs = pot_p + pot_delta

       ! derivative of c_eps wrt zeta (= derivative of eps_delta wrt zeta)
       ! term1 = d/dzeta [ alpha_c * f(zeta)/f''(0) * (1-zeta**4) ]
       ! term2 = d/dzeta [ f(zeta) * zeta**4 * (eps_f - eps_p) ]
       de_dzeta1 = alpha_c / d2f_at_zero * (fprime_of_zeta * (1.0_DP - zeta4) +&
            f_of_zeta * (-4.0_DP*zeta3))
       de_dzeta2 = (fprime_of_zeta * zeta4 + 4.0_DP * zeta3 * f_of_zeta) * &
            (eps_f - eps_p)

       de_dzeta = de_dzeta1 + de_dzeta2

       ! Energy density -> energy
       c_energy = c_eps * den

       ! drs/drho
       drs_drho = minus_fourpi_over_nine * rs**4.0_DP

       ! Resulting potential (functional derivative of E_c (not eps_c) wrt rho)
       ! (v_c = eps_c + rho * de_c/drho), where de_c/drho = de_c/drs * drs/drho
       c_pot = c_eps + den * de_drs * drs_drho

       ! This is the spin-polarized case, take zeta into account
       c_pot1 = c_pot - (zeta - 1.0_DP) * de_dzeta
       c_pot2 = c_pot - (zeta + 1.0_DP) * de_dzeta

    else
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
    end if

  end subroutine xc_vwn_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vwn_eps_c_helper(c_eps, c_pot, rs, paramset)
    !===========================================================================!
    ! Given the Wigner-Seitz radius, rs, calculates the correlation energy      !
    ! density, according to expression (4.4) in [1]. The first derivative       !
    ! wrt rs is also calculated.                                                !
    !---------------------------------------------------------------------------!
    ! Arguments:                                                                !
    !   c_eps    (out): Returned energy density.                                !
    !   c_pot    (out): Returned potential.                                     !
    !   rs       (in):  Wigner-Seitz radius (a simple function of density)      !
    !   paramset (in):  Selects the VWN parameter set:                          !
    !                   1 - paramagnetic case                                   !
    !                   2 - ferromagnetic case                                  !
    !                   3 - calculation of alpha_c (cf. text below (9) in [2].  !
    !---------------------------------------------------------------------------!
    ! Reference:                                                                !
    !  [1] S H Vosko, L Wilk & M Nusair, Can. J. Phys. 58, 1200 (1980)          !
    !  [2] S H Vosko and L Wilk, Phys. Rev. B 22, 3812 (1980)                   !
    !---------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 01/06/2011.                                    !
    !===========================================================================!
    use constants, only: DP, PI
    use utils, only: utils_abort
    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: c_eps
    real(kind=DP), intent(out) :: c_pot
    real(kind=DP), intent(in)  :: rs
    integer                    :: paramset

    ! jd: Parameters of the method, following [1], [2]
    !     Note A(:) are divided by 2 because of Ry->Ha conversion.
    !     This also applies to A(3), which becomes 1/(6 pi^2) rather than
    !     1/(3 pi^2)
    real(kind=DP), parameter :: As(3) = &
         (/ 0.0310907_DP, 0.01554535_DP, -0.01688686394038963_DP /)
    real(kind=DP), parameter :: bs(3) = (/ 3.72744_DP, 7.06042_DP, 1.13107_DP /)
    real(kind=DP), parameter :: cs(3) = (/ 12.9352_DP, 18.0578_DP, 13.0045_DP /)
    real(kind=DP), parameter :: x0s(3) = &
         (/ -0.10498_DP, -0.325_DP, -0.0047584_DP /)
    real(kind=DP), parameter :: Qs(3) = &
         (/ 6.15199082_DP, 4.73092691_DP, 7.123108918_DP /) ! sqrt(4c-b2)
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! jd: Internal variables
    real(kind=DP) :: x                 ! jd: =rs^(1/2), cf. [1], below [4.3]
    real(kind=DP) :: X_of_x, X_of_x0   ! jd:            cf. [1], below [4.4]
    real(kind=DP) :: den               ! jd: the density corresponding to rs
    real(kind=DP) :: ln1, ln2, arctan1 ! jd: Temporaries
    real(kind=DP) :: A,b,c,x0,Q        ! jd: Parameters for current paramset

    ! -------------------------------------------------------------------------

    if(paramset < 1 .or. paramset > 3) then
       call utils_abort("Bad paramset in xc_vwn_eps_c_helper")
    endif

    A = As(paramset)
    b = bs(paramset)
    c = cs(paramset)
    x0 = x0s(paramset)
    Q = Qs(paramset)

    x = sqrt(rs) ! rs>0 is guaranteed by the caller
    den = TONFPI / (rs*rs*rs)
    X_of_x = x*x + b*x + c
    X_of_x0 = x0*x0 + b*x0 + c

    ln1 = log(x*x/X_of_x)
    ln2 = log((x-x0)*(x-x0)/X_of_x)
    arctan1 = atan(Q/(2.0_DP*x+b))

    ! Calculate energy density according to (4.4)
    c_eps = A* &
         (ln1 + 2.0_DP*b/Q * arctan1 - b*x0/X_of_x0 * &
         (ln2 + 2.0_DP*(b+2.0_DP*x0)/Q * arctan1))

    ! Derivative of c_eps wrt rs
    c_pot = (A*c*x-A * (c+b*x)*x0) / (rs*(c+b*x+rs)*(x-x0))

  end subroutine xc_vwn_eps_c_helper


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw92_c_point(den,c_energy,c_pot)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Wang.                                                            !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)                     !
    !---------------------------------------------------------------------------!
    ! Originally written by Nicholas Hine, March 2011.                          !
    !===========================================================================!


    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot

    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: ecunif,eurs

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! Constants for spin-unpolarised case:

    tol: if (den > dentol) then   ! Positive charge density
       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,ecunif,eurs)
       c_energy = den * ecunif
       c_pot = ecunif-THIRD*rs*eurs

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
    end if tol

  end subroutine xc_pw92_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw92_c_point_sp(den,den1,den2,c_energy,c_pot1,c_pot2)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Zunger.                                                          !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)                     !
    !---------------------------------------------------------------------------!
    ! Originally written by Nicholas Hine, March 2011.                          !
    !===========================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(in)  :: den1
    real(kind=DP), intent(in)  :: den2
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot1
    real(kind=DP), intent(out) :: c_pot2

    real(kind=DP) :: zeta         ! Spin polarisation
    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: ecunif               ! LDA correlation energy per particle
    real(kind=DP) :: ecd,ecz,eu,eurs,ecrs
    real(kind=DP) :: eprs,ep,alfm,alfrsm
    real(kind=DP) :: eczet
    real(kind=DP) :: zeta3,zeta4          ! zeta**3 and zeta**4
    real(kind=DP) :: fzeta,dfdz           ! Function of zeta and derivative

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: zetamin = -1.0_DP + 1.0e-10_DP
    real(kind=DP), parameter :: zetamax = 1.0_DP - 1.0e-10_DP
    real(kind=DP), parameter :: gaminv = 1.92366105093153631981_DP
    real(kind=DP), parameter :: fzzinv = 1.0_DP / (8.0_DP * gaminv / 9.0_DP)

    tol: if (den > dentol) then
    ! .and. den1 > 0.0_DP .and. den2 > 0.0_DP) then !! ddor

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vectors

       zeta = (den1 - den2) / den   ! Polarisation
       zeta = max(zeta,zetamin)
       zeta = min(zeta,zetamax)

       ! PW91 correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,eu,eurs)
       call xc_pw91_eq10(0.01554535_DP,0.20548_DP,14.1189_DP, &
            6.1977_DP,3.3662_DP,0.62517_DP,rs,sqrtrs,ep,eprs)
       call xc_pw91_eq10(0.0168869_DP,0.11125_DP,10.357_DP, &
            3.6231_DP,0.88026_DP,0.49671_DP,rs,sqrtrs,alfm,alfrsm)

       fzeta = gaminv * &
            ((1.0_DP+zeta)**FTHRD+(1.0_DP-zeta)**FTHRD-2.0_DP)
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta3*zeta

       ecunif = eu*(1.0_DP-fzeta*zeta4) + ep*fzeta*zeta4 - &
            alfm*fzeta*(1.0_DP-zeta4)*fzzinv
       ecrs = eurs*(1.0_DP-fzeta*zeta4) + eprs*fzeta*zeta4 - &
            alfrsm*fzeta*(1.0_DP-zeta4)*fzzinv
       dfdz = gaminv * FTHRD * &
            ((1.0_DP+zeta)**THIRD-(1.0_DP-zeta)**THIRD)
       eczet = (4.0_DP*zeta3*fzeta + zeta4*dfdz) * &
            (ep-eu+alfm*fzzinv)-dfdz*alfm*fzzinv

       ecd = ecunif-THIRD*rs*ecrs
       ecz = eczet * den

       c_energy = den * ecunif
       c_pot1 = ecd+(1.0_DP-zeta)*ecz/den
       c_pot2 = ecd-(1.0_DP+zeta)*ecz/den

    else

       ! qoh: Set all values to zero for very small densities
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP

    end if tol

  end subroutine xc_pw92_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PBE exchange at a point       !
    ! in the spin unpolarised case.  (As in original paper.)       !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fxpbe,fs             ! PBE exchange variables
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       exunif = ax*den**THIRD

       ! qoh: PBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s*s
       p1 = 1.0_DP / p0
       fxpbe = 1.0_DP + uk - uk * p1
       fs = 2.0_DP * um * s * p1 * p1

       x_energy =  den * exunif * fxpbe
       x_pot = FTHRD * exunif*(fxpbe - fs*s)
       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pbe_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_rpbe_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the RPBE exchange at a point      !
    ! in the spin unpolarised case.                                !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Uses  the RPBE exchange enhancement  factor and derivative   !
    !   Hammer, Hansen & Norskov, Phys. Rev. B 59, 7413 (1999)     !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fxpbe,fs             ! PBE exchange variables
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: p0

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       exunif = ax*den**THIRD
       ! qoh: RPBE exchange enhancement factor and derivative
       p0 = exp(-ul*s*s)
       fxpbe = 1.0_DP + uk - uk * p0
       fs = 2.0_DP * um * s * p0

       x_energy =  den * exunif * fxpbe
       x_pot = FTHRD * exunif*(fxpbe - fs*s)
       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_rpbe_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_revpbe_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the revPBE exchange at a point    !
    ! in the spin unpolarised case.                                !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Uses the revised PBE exchange enhancement factor and         !
    ! derivative                                                   !
    ! Zhang & Yang, Phys. Rev. Lett. 80, 890 (1998)                !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fxpbe,fs             ! PBE exchange variables
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 1.245_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       exunif = ax*den**THIRD
       ! qoh: revPBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s*s
       p1 = 1.0_DP / p0
       fxpbe = 1.0_DP + uk - uk * p1
       fs = 2.0_DP * um * s * p1 * p1

       x_energy =  den * exunif * fxpbe
       x_pot = FTHRD * exunif*(fxpbe - fs*s)
       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_revpbe_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbesol_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PBEsol exchange at a point    !
    ! in the spin unpolarised case.                                !
    ! J P Perdew et al. Phys. Rev. Lett. 100, 136406 (2008).       !
    !--------------------------------------------------------------!
    ! Written David O'Regan in 12/2013 based on xc_pbe_x_point.    !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fxpbe,fs             ! PBE exchange variables
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.123456790123456_DP
    real(kind=DP), parameter :: ul = um / uk

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       exunif = ax*den**THIRD

       ! qoh: PBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s*s
       p1 = 1.0_DP / p0
       fxpbe = 1.0_DP + uk - uk * p1
       fs = 2.0_DP * um * s * p1 * p1

       x_energy =  den * exunif * fxpbe
       x_pot = FTHRD * exunif*(fxpbe - fs*s)
       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pbesol_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_wc_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the WC exchange at a point        !
    ! in the spin unpolarised case.                                !
    ! Z Wu and R E Cohen, Phys. Rev. B 73, 235116 (2006)           !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    ! Adapted for Wu-Cohen functional by Nicholas Hine, 03/05/2010.!
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fxwc,fs              ! WC exchange variables
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: p1
    real(kind=DP) :: x

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: c  = 0.007937469335162_DP

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       exunif = ax*den**THIRD

       ! ndmh: WC x function
       x = 10.0_DP/81.0_DP*s*s+(um-10.0_DP/81.0_DP)*s*s*exp(-s*s)+log(1.0_DP+c*s**4)
       ! ndmh: WC exchange enhancement factor and derivative
       ! ndmh: fxwc = 1 + uk - uk/(1+x/uk) (Eqn 4 of Ref)
       p1 = 1.0_DP / (1.0_DP + x/uk)
       fxwc = 1.0_DP + uk - uk * p1
       fs = (2.0_DP*s*(10.0_dp/81.0_dp+(um-10.0_dp/81.0_dp)*(1.0_DP-s*s)* &
            exp(-s*s)+2.0_DP*c*s*s/(1.0_DP+c*s**4)))*p1*p1

       x_energy =  den * exunif * fxwc
       x_pot = FTHRD * exunif*(fxwc - fs*s)
       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_wc_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe0_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PBE0 exchange at a point      !
    ! in the spin unpolarised case.                                !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Hybrid functional form from:                                 !
    !  Adamo, C., Cossi, Maurizio, and Barone, Vincenzo:           !
    !    J.Mol. Struc. (Theochem) 493, 145 (1999)                  !
    ! PBE0 = PBE + 1/4 (HF_X - PBE_X)                              !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP), parameter   :: pbe_fac = 0.75_DP ! PBE0 factor (1 - 1/4)

    call xc_pbe_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)
    x_energy = pbe_fac * x_energy
    x_pot    = pbe_fac * x_pot
    x_dfdmgd = pbe_fac * x_dfdmgd

  end subroutine xc_pbe0_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PBE correlation at a point    !
    ! in the spin unpolarised case.                                !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation per particle
    real(kind=DP) :: q1,q4                ! LDA correlation variables
    real(kind=DP) :: t2,bt2,q5,q6,q7      ! PBE correlation variables
    real(kind=DP) :: h,h1,ha,ht,pon,b
    real(kind=DP) :: eurs,aec,ecn,ect

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI

    ! PBE correlation constants
    real(kind=dp), parameter :: gamma = 0.03109069086965489503494086371273_DP
    real(kind=dp), parameter :: beta = 0.06672455060314922_DP
    real(kind=dp), parameter :: delta = beta / gamma

    tol: if (den > dentol) then

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector
       t = mgd / (2.0_DP*den*ks)    ! Dimensionless gradient

       ! PBE correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,ecunif,eurs)
       pon = -ecunif/gamma
       q1 = exp(pon)
       b = delta/(q1-1.0_DP)
       t2 = t*t
       bt2 = b*t2
       q4 = 1.0_DP+bt2
       q5 = 1.0_DP+bt2+bt2*bt2
       h = gamma*log(1.0_DP+delta*q4*t2/q5)
       ect = ecunif + h
       c_energy = den * ect
       q6 = 1.0_DP/((gamma*q5+beta*t2*q4)*q5)
       q7 = gamma*q6
       ha = -beta*t2*t2*t2*b*(1.0_DP+q4)*q7
       aec = b*b*q1/beta
       ht = 2.0_DP*beta*t*(1.0_DP+2.0_DP*bt2)*q7
       h1 = (rs*ha*aec*eurs+3.5_DP*t*ht)*THIRD
       ecn = THIRD*rs*eurs
       c_pot = ect - (ecn + h1)
       c_dfdmgd = 0.5_DP * ht / ks

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
       c_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pbe_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbesol_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PBEsol correlation at a point !
    ! in the spin unpolarised case.                                !
    ! J P Perdew et al. Phys. Rev. Lett. 100, 136406 (2008).       !
    !--------------------------------------------------------------!
    ! Written David O'Regan in 12/2013 based on xc_pbe_c_point.    !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation per particle
    real(kind=DP) :: q1,q4                ! LDA correlation variables
    real(kind=DP) :: t2,bt2,q5,q6,q7      ! PBE correlation variables
    real(kind=DP) :: h,h1,ha,ht,pon,b
    real(kind=DP) :: eurs,aec,ecn,ect

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI

    ! PBEsol correlation constants
    real(kind=dp), parameter :: gamma = 0.03109069086965489503494086371273_DP
    real(kind=dp), parameter :: beta = 0.046_DP
    real(kind=dp), parameter :: delta = beta / gamma

    tol: if (den > dentol) then

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector
       t = mgd / (2.0_DP*den*ks)    ! Dimensionless gradient

       ! PBEsol correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,ecunif,eurs)
       pon = -ecunif/gamma
       q1 = exp(pon)
       b = delta/(q1-1.0_DP)
       t2 = t*t
       bt2 = b*t2
       q4 = 1.0_DP+bt2
       q5 = 1.0_DP+bt2+bt2*bt2
       h = gamma*log(1.0_DP+delta*q4*t2/q5)
       ect = ecunif + h
       c_energy = den * ect
       q6 = 1.0_DP/((gamma*q5+beta*t2*q4)*q5)
       q7 = gamma*q6
       ha = -beta*t2*t2*t2*b*(1.0_DP+q4)*q7
       aec = b*b*q1/beta
       ht = 2.0_DP*beta*t*(1.0_DP+2.0_DP*bt2)*q7
       h1 = (rs*ha*aec*eurs+3.5_DP*t*ht)*THIRD
       ecn = THIRD*rs*eurs
       c_pot = ect - (ecn + h1)
       c_dfdmgd = 0.5_DP * ht / ks

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
       c_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pbesol_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the PBE exchange at a point       !
    ! in the spin polarised case.  (As in original paper.)         !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fs1,fs2              ! PBE exchange variables
    real(kind=DP) :: fxpbe1,fxpbe2
    real(kind=DP) :: exunif1,exunif2      ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    !qoh: Spin 1:
    tol1: if (den1 > dentol*0.5_DP) then
       kf1 = (SXPISQ * den1)**THIRD

       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradient

       exunif1 = ax*(2.0_DP*den1)**THIRD
       ! qoh: PBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s1*s1
       p1 = 1.0_DP / p0
       fxpbe1 = 1.0_DP + uk - uk * p1
       fs1 = 2.0_DP * um * s1 * p1 * p1

       x_energy = den1 * exunif1 * fxpbe1
       x_pot1 = FTHRD*exunif1*(fxpbe1 - fs1*s1)
       x_dfdmgd1 = 0.5_DP * ax * fs1 * TTPI23

    end if tol1

    !qoh:  Spin 2:
    tol2: if (den2 > 0.5_DP*dentol) then

       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) !Dimensionless gradient
       exunif2 = ax*(2.0_DP*den2)**THIRD
       ! qoh: PBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s2*s2
       p1 = 1.0_DP / p0
       fxpbe2 = 1.0_DP + uk - uk * p1
       fs2 = 2.0_DP * um * s2 * p1 * p1

       x_energy = x_energy + den2 * exunif2 * fxpbe2
       x_pot2 = FTHRD*exunif2*(fxpbe2 - fs2*s2)
       x_dfdmgd2 = 0.5_DP * ax * fs2 * TTPI23

    end if tol2

  end subroutine xc_pbe_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_rpbe_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the RPBE exchange at a point      !
    ! in the spin polarised case.                                  !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Uses  the RPBE exchange enhancement  factor and derivative   !
    !   Hammer, Hansen & Norskov, Phys. Rev. B 59, 7413 (1999)     !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fs1,fs2              ! PBE exchange variables
    real(kind=DP) :: fxpbe1,fxpbe2
    real(kind=DP) :: exunif1,exunif2      ! LDA exchange per particle
    real(kind=DP) :: p0

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    !qoh: Spin 1:
    tol1: if (den1 > dentol*0.5_DP) then
       kf1 = (SXPISQ * den1)**THIRD

       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradient

       exunif1 = ax*(2.0_DP*den1)**THIRD
       ! qoh: RPBE exchange enhancement factor and derivative
       p0 = exp(-ul*s1*s1)
       fxpbe1 = 1.0_DP + uk - uk * p0
       fs1 = 2.0_DP * um * s1 * p0

       x_energy = den1 * exunif1 * fxpbe1
       x_pot1 = FTHRD*exunif1*(fxpbe1 - fs1*s1)
       x_dfdmgd1 = 0.5_DP * ax * fs1 * TTPI23

    end if tol1

    !qoh:  Spin 2:
    tol2: if (den2 > 0.5_DP*dentol) then

       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) !Dimensionless gradient
       exunif2 = ax*(2.0_DP*den2)**THIRD
       ! qoh: RPBE exchange enhancement factor and derivative
       p0 = exp(-ul*s2*s2)
       fxpbe2 = 1.0_DP + uk - uk * p0
       fs2 = 2.0_DP * um * s2 * p0

       x_energy = x_energy + den2 * exunif2 * fxpbe2
       x_pot2 = FTHRD*exunif2*(fxpbe2 - fs2*s2)
       x_dfdmgd2 = 0.5_DP * ax * fs2 * TTPI23

    end if tol2

  end subroutine xc_rpbe_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_revpbe_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the revPBE exchange at a point    !
    ! in the spin polarised case.                                  !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Uses the revised PBE exchange enhancement factor and         !
    ! derivative                                                   !
    ! Zhang & Yang, Phys. Rev. Lett. 80, 890 (1998)                !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fs1,fs2              ! PBE exchange variables
    real(kind=DP) :: fxpbe1,fxpbe2
    real(kind=DP) :: exunif1,exunif2      ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 1.245_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    !qoh: Spin 1:
    tol1: if (den1 > dentol*0.5_DP) then
       kf1 = (SXPISQ * den1)**THIRD

       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradient

       exunif1 = ax*(2.0_DP*den1)**THIRD
       ! qoh: revPBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s1*s1
       p1 = 1.0_DP / p0
       fxpbe1 = 1.0_DP + uk - uk * p1
       fs1 = 2.0_DP * um * s1 * p1 * p1

       x_energy = den1 * exunif1 * fxpbe1
       x_pot1 = FTHRD*exunif1*(fxpbe1 - fs1*s1)
       x_dfdmgd1 = 0.5_DP * ax * fs1 * TTPI23

    end if tol1

    !qoh:  Spin 2:
    tol2: if (den2 > 0.5_DP*dentol) then

       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) !Dimensionless gradient
       exunif2 = ax*(2.0_DP*den2)**THIRD
       ! qoh: revPBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s2*s2
       p1 = 1.0_DP / p0
       fxpbe2 = 1.0_DP + uk - uk * p1
       fs2 = 2.0_DP * um * s2 * p1 * p1

       x_energy = x_energy + den2 * exunif2 * fxpbe2
       x_pot2 = FTHRD*exunif2*(fxpbe2 - fs2*s2)
       x_dfdmgd2 = 0.5_DP * ax * fs2 * TTPI23

    end if tol2

  end subroutine xc_revpbe_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbesol_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the PBEsol exchange at a point    !
    ! in the spin polarised case.                                  !
    ! J P Perdew et al. Phys. Rev. Lett. 100, 136406 (2008).       !
    !--------------------------------------------------------------!
    ! Written David O'Regan in 12/2013 based on xc_pbe_x_point_sp. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fs1,fs2              ! PBE exchange variables
    real(kind=DP) :: fxpbe1,fxpbe2
    real(kind=DP) :: exunif1,exunif2      ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.123456790123456_DP
    real(kind=DP), parameter :: ul = um / uk

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    !qoh: Spin 1:
    tol1: if (den1 > dentol*0.5_DP) then
       kf1 = (SXPISQ * den1)**THIRD

       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradient

       exunif1 = ax*(2.0_DP*den1)**THIRD
       ! qoh: PBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s1*s1
       p1 = 1.0_DP / p0
       fxpbe1 = 1.0_DP + uk - uk * p1
       fs1 = 2.0_DP * um * s1 * p1 * p1

       x_energy = den1 * exunif1 * fxpbe1
       x_pot1 = FTHRD*exunif1*(fxpbe1 - fs1*s1)
       x_dfdmgd1 = 0.5_DP * ax * fs1 * TTPI23

    end if tol1

    !qoh:  Spin 2:
    tol2: if (den2 > 0.5_DP*dentol) then

       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) !Dimensionless gradient
       exunif2 = ax*(2.0_DP*den2)**THIRD
       ! qoh: PBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s2*s2
       p1 = 1.0_DP / p0
       fxpbe2 = 1.0_DP + uk - uk * p1
       fs2 = 2.0_DP * um * s2 * p1 * p1

       x_energy = x_energy + den2 * exunif2 * fxpbe2
       x_pot2 = FTHRD*exunif2*(fxpbe2 - fs2*s2)
       x_dfdmgd2 = 0.5_DP * ax * fs2 * TTPI23

    end if tol2

  end subroutine xc_pbesol_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_wc_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the WC exchange at a point        !
    ! in the spin polarised case.                                  !
    ! Z Wu and R E Cohen, Phys. Rev. B 73, 235116 (2006)           !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    ! Adapted for Wu-Cohen functional by Nicholas Hine, 03/05/2010.!
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fs1,fs2              ! PBE exchange variables
    real(kind=DP) :: fxwc1,fxwc2
    real(kind=DP) :: exunif1,exunif2      ! LDA exchange per particle
    real(kind=DP) :: p1
    real(kind=DP) :: x1,x2

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: c  = 0.007937469335162_DP

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    ! ndmh: Spin 1:
    tol1: if (den1 > dentol*0.5_DP) then
       kf1 = (SXPISQ * den1)**THIRD

       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradient
       exunif1 = ax*(2.0_DP*den1)**THIRD

       ! ndmh: WC x function
       x1 = 10.0_DP/81.0_DP*s1*s1+(um-10.0_DP/81.0_DP)*s1*s1*exp(-s1*s1)+log(1+c*s1**4)
       ! ndmh: WC exchange enhancement factor and derivative
       ! ndmh: fxwc = 1 + uk - uk/(1+x/uk) (Eqn 4 of Ref)
       p1 = 1.0_DP / (1.0_DP + x1/uk)
       fxwc1 = 1.0_DP + uk - uk * p1
       fs1 = (2.0_DP*s1*(10.0_dp/81.0_dp+(um-10.0_dp/81.0_dp)*(1.0_DP-s1*s1)* &
            exp(-s1*s1)+2.0_DP*c*s1*s1/(1.0_DP+c*s1**4)))*p1*p1

       x_energy = den1 * exunif1 * fxwc1
       x_pot1 = FTHRD*exunif1*(fxwc1 - fs1*s1)
       x_dfdmgd1 = 0.5_DP * ax * fs1 * TTPI23

    end if tol1

    ! ndmh:  Spin 2:
    tol2: if (den2 > 0.5_DP*dentol) then

       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) !Dimensionless gradient
       exunif2 = ax*(2.0_DP*den2)**THIRD

       ! ndmh: WC x function
       x2 = 10.0_DP/81.0_DP*s2*s2+(um-10.0_DP/81.0_DP)*s2*s2*exp(-s2*s2)+log(1+c*s2**4)
       ! ndmh: WC exchange enhancement factor and derivative
       ! ndmh: fxwc = 1 + uk - uk/(1+x/uk) (Eqn 4 of Ref)
       p1 = 1.0_DP / (1.0_DP + x2/uk)
       fxwc2 = 1.0_DP + uk - uk * p1
       fs2 = (2.0_DP*s2*(10.0_dp/81.0_dp+(um-10.0_dp/81.0_dp)*(1.0_DP-s2*s2)* &
            exp(-s2*s2)+2.0_DP*c*s2*s2/(1.0_DP+c*s2**4)))*p1*p1

       x_energy = x_energy + den2 * exunif2 * fxwc2
       x_pot2 = FTHRD*exunif2*(fxwc2 - fs2*s2)
       x_dfdmgd2 = 0.5_DP * ax * fs2 * TTPI23

    end if tol2

  end subroutine xc_wc_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe0_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the PBE0 exchange at a point      !
    ! in the spin polarised case.                                  !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Hybrid functional form from:                                 !
    !  Adamo, C., Cossi, Maurizio, and Barone, Vincenzo:           !
    !    J.Mol. Struc. (Theochem) 493, 145 (1999)                  !
    ! PBE0 = PBE + 1/4 (HF_X - PBE_X)                              !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP), parameter   :: xfac = 0.75_DP ! PBE0 factor (1 - 1/4)

    call xc_pbe_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)
    x_energy  = xfac * x_energy
    x_pot1    = xfac * x_pot1
    x_pot2    = xfac * x_pot2
    x_dfdmgd1 = xfac * x_dfdmgd1
    x_dfdmgd2 = xfac * x_dfdmgd2

  end subroutine xc_pbe0_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the PBE correlation at a point    !
    ! in the spin polarised case.                                  !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation per particle
    real(kind=DP) :: q1,q2,q3,q4          ! LDA correlation variables
    real(kind=DP) :: t2,q5,q6,q7,q8       ! PBE correlation variables
    real(kind=DP) :: h,ht,pon,b,hb
    real(kind=DP) :: eurs
    real(kind=DP) :: eprs,ep,alfm,alfrsm
    real(kind=DP) :: ecd,ecz,eu,ecrs
    real(kind=DP) :: ec1d,ec1z,eczet
    real(kind=DP) :: g2,g3,g4,be,bg,gz
    real(kind=DP) :: b2,t3,t4,g
    real(kind=DP) :: zeta                 ! Spin polarisation
    real(kind=DP) :: zeta3,zeta4          ! zeta**3 and zeta**4
    real(kind=DP) :: fzeta,dfdz           ! Function of zeta and derivative
    real(kind=DP) :: mgd12                ! mgd1 + mgd2 to avoid compiler warning

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI
    real(kind=DP), parameter :: TWTHRD = 2.0_DP / 3.0_DP

    ! PBE correlation constants
    real(kind=dp), parameter :: gamma = 0.03109069086965489503494086371273_DP
    real(kind=dp), parameter :: beta = 0.06672455060314922_DP
    real(kind=dp), parameter :: delta = beta / gamma
    real(kind=DP), parameter :: zetamin = -1.0_DP + 1.0e-10_DP
    real(kind=DP), parameter :: zetamax = 1.0_DP - 1.0e-10_DP

    ! PW91 correlation parameters
    real(kind=DP), parameter :: gaminv = 1.92366105093153631981_DP
    real(kind=DP), parameter :: fzzinv = 1.0_DP / (8.0_DP * gaminv / 9.0_DP)

    ! qoh: Avoid compiler warning by making use of mgd1 and mgd2
    mgd12 = mgd1 + mgd2

    tol: if (den > dentol) then
    ! .and. den1 > 0.0_DP .and. den2 > 0.0_DP) then !! ddor

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vectors
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector

       zeta = (den1 - den2) / den   ! Polarisation
       zeta = max(zeta,zetamin)
       zeta = min(zeta,zetamax)
       g = 0.5_DP*((1.0_DP-zeta)**TWTHRD+(1.0_DP+zeta)**TWTHRD)

       t = mgd / (den*ks*g*2.0_DP)! Dimensionless gradient

       ! PBE correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,eu,eurs)
       call xc_pw91_eq10(0.01554535_DP,0.20548_DP, &
            14.1189_DP,6.1977_DP,3.3662_DP,0.62517_DP,rs,sqrtrs, &
            ep,eprs)
       call xc_pw91_eq10(0.0168869_DP,0.11125_DP,10.357_DP, &
            3.6231_DP,0.88026_DP,0.49671_DP,rs,sqrtrs, &
            alfm,alfrsm)

       fzeta = gaminv * &
            ((1.0_DP+zeta)**FTHRD+(1.0_DP-zeta)**FTHRD-2.0_DP)
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta3*zeta

       ecunif = eu*(1.0_DP-fzeta*zeta4) + ep*fzeta*zeta4 - &
            alfm*fzeta*(1.0_DP-zeta4)*fzzinv
       ecrs = eurs*(1.0_DP-fzeta*zeta4) + eprs*fzeta*zeta4 - &
            alfrsm*fzeta*(1.0_DP-zeta4)*fzzinv
       dfdz = gaminv * FTHRD * &
            ((1.0_DP+zeta)**THIRD-(1.0_DP-zeta)**THIRD)
       eczet = (4.0_DP*zeta3*fzeta+zeta4*dfdz) * &
            (ep-eu+alfm*fzzinv)-dfdz*alfm*fzzinv
       ecd = ecunif - THIRD*rs*ecrs
       ecz = eczet * den

       g2 = g*g
       g3 = g2*g
       g4 = g2*g2
       pon = -ecunif/(g3*gamma)
       ep = exp(pon)
       b = delta/(ep-1.0_DP)
       b2 = b*b
       t2 = t*t
       t3 = t2*t
       t4 = t2*t2
       q1 = t2 + b*t4
       q2 = 1.0_DP+b*t2+b2*t4
       q3 = q1/q2
       q4 = q2*q2 + delta*q1*q2
       q5 = beta*g3/q4
       q6 = 2.0_DP*t+4.0_DP*b*t3
       q7 = -b*t4*t2*(2.0_DP+b*t2)
       q8 = ep*b2/beta
       h = gamma*g3*log(1.0_DP + delta*q3)
       ht = q5*q6
       hb = q5*q7
       be = q8/g3
       bg = -3.0_DP*ecunif*q8/g4
       gz = ((1.0_DP+zeta)**(-THIRD)-(1.0_DP-zeta)**(-THIRD))*THIRD
       ec1d = h - THIRD*rs*ecrs*be*hb - SSXTH*t*ht
       ec1z = den * (3.0_DP*h*gz/g + hb*(bg*gz + be*eczet) - ht*t/g*gz)

       c_energy = den * (ecunif + h)
       c_pot1 = ecd+ec1d+(1.0_DP-zeta)*(ecz+ec1z)/den
       c_pot2 = ecd+ec1d-(1.0_DP+zeta)*(ecz+ec1z)/den

       ! qoh: Correlation gradient dependent part
       ! ndmh: separated c_dfdmgd from c_dfdmgd1
       c_dfdmgd = 0.5_DP * ht / (g*ks)
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    else
       c_energy = 0.0_DP ! ddor
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
       c_dfdmgd = 0.0_DP
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    end if tol

  end subroutine xc_pbe_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbesol_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the PBEsol correlationat a point  !
    ! in the spin polarised case.                                  !
    ! J P Perdew et al. Phys. Rev. Lett. 100, 136406 (2008).       !
    !--------------------------------------------------------------!
    ! Written David O'Regan in 12/2013 based on xc_pbe_c_point_sp. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation per particle
    real(kind=DP) :: q1,q2,q3,q4          ! LDA correlation variables
    real(kind=DP) :: t2,q5,q6,q7,q8       ! PBE correlation variables
    real(kind=DP) :: h,ht,pon,b,hb
    real(kind=DP) :: eurs
    real(kind=DP) :: eprs,ep,alfm,alfrsm
    real(kind=DP) :: ecd,ecz,eu,ecrs
    real(kind=DP) :: ec1d,ec1z,eczet
    real(kind=DP) :: g2,g3,g4,be,bg,gz
    real(kind=DP) :: b2,t3,t4,g
    real(kind=DP) :: zeta                 ! Spin polarisation
    real(kind=DP) :: zeta3,zeta4          ! zeta**3 and zeta**4
    real(kind=DP) :: fzeta,dfdz           ! Function of zeta and derivative
    real(kind=DP) :: mgd12                ! mgd1 + mgd2 to avoid compiler warning

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI
    real(kind=DP), parameter :: TWTHRD = 2.0_DP / 3.0_DP

    ! PBEsol correlation constants
    real(kind=dp), parameter :: gamma = 0.03109069086965489503494086371273_DP
    real(kind=dp), parameter :: beta = 0.046_DP
    real(kind=dp), parameter :: delta = beta / gamma
    real(kind=DP), parameter :: zetamin = -1.0_DP + 1.0e-10_DP
    real(kind=DP), parameter :: zetamax = 1.0_DP - 1.0e-10_DP

    ! PW91 correlation parameters
    real(kind=DP), parameter :: gaminv = 1.92366105093153631981_DP
    real(kind=DP), parameter :: fzzinv = 1.0_DP / (8.0_DP * gaminv / 9.0_DP)

    ! qoh: Avoid compiler warning by making use of mgd1 and mgd2
    mgd12 = mgd1 + mgd2

    tol: if (den > dentol) then
    ! .and. den1 > 0.0_DP .and. den2 > 0.0_DP) then !! ddor

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vectors
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector

       zeta = (den1 - den2) / den   ! Polarisation
       zeta = max(zeta,zetamin)
       zeta = min(zeta,zetamax)
       g = 0.5_DP*((1.0_DP-zeta)**TWTHRD+(1.0_DP+zeta)**TWTHRD)

       t = mgd / (den*ks*g*2.0_DP)! Dimensionless gradient

       ! PBEsol correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,eu,eurs)
       call xc_pw91_eq10(0.01554535_DP,0.20548_DP, &
            14.1189_DP,6.1977_DP,3.3662_DP,0.62517_DP,rs,sqrtrs, &
            ep,eprs)
       call xc_pw91_eq10(0.0168869_DP,0.11125_DP,10.357_DP, &
            3.6231_DP,0.88026_DP,0.49671_DP,rs,sqrtrs, &
            alfm,alfrsm)

       fzeta = gaminv * &
            ((1.0_DP+zeta)**FTHRD+(1.0_DP-zeta)**FTHRD-2.0_DP)
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta3*zeta

       ecunif = eu*(1.0_DP-fzeta*zeta4) + ep*fzeta*zeta4 - &
            alfm*fzeta*(1.0_DP-zeta4)*fzzinv
       ecrs = eurs*(1.0_DP-fzeta*zeta4) + eprs*fzeta*zeta4 - &
            alfrsm*fzeta*(1.0_DP-zeta4)*fzzinv
       dfdz = gaminv * FTHRD * &
            ((1.0_DP+zeta)**THIRD-(1.0_DP-zeta)**THIRD)
       eczet = (4.0_DP*zeta3*fzeta+zeta4*dfdz) * &
            (ep-eu+alfm*fzzinv)-dfdz*alfm*fzzinv
       ecd = ecunif - THIRD*rs*ecrs
       ecz = eczet * den

       g2 = g*g
       g3 = g2*g
       g4 = g2*g2
       pon = -ecunif/(g3*gamma)
       ep = exp(pon)
       b = delta/(ep-1.0_DP)
       b2 = b*b
       t2 = t*t
       t3 = t2*t
       t4 = t2*t2
       q1 = t2 + b*t4
       q2 = 1.0_DP+b*t2+b2*t4
       q3 = q1/q2
       q4 = q2*q2 + delta*q1*q2
       q5 = beta*g3/q4
       q6 = 2.0_DP*t+4.0_DP*b*t3
       q7 = -b*t4*t2*(2.0_DP+b*t2)
       q8 = ep*b2/beta
       h = gamma*g3*log(1.0_DP + delta*q3)
       ht = q5*q6
       hb = q5*q7
       be = q8/g3
       bg = -3.0_DP*ecunif*q8/g4
       gz = ((1.0_DP+zeta)**(-THIRD)-(1.0_DP-zeta)**(-THIRD))*THIRD
       ec1d = h - THIRD*rs*ecrs*be*hb - SSXTH*t*ht
       ec1z = den * (3.0_DP*h*gz/g + hb*(bg*gz + be*eczet) - ht*t/g*gz)

       c_energy = den * (ecunif + h)
       c_pot1 = ecd+ec1d+(1.0_DP-zeta)*(ecz+ec1z)/den
       c_pot2 = ecd+ec1d-(1.0_DP+zeta)*(ecz+ec1z)/den

       ! qoh: Correlation gradient dependent part
       ! ndmh: separated c_dfdmgd from c_dfdmgd1
       c_dfdmgd = 0.5_DP * ht / (g*ks)
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    else
       c_energy = 0.0_DP ! ddor
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
       c_dfdmgd = 0.0_DP
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    end if tol

  end subroutine xc_pbesol_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PW91 exchange at a point      !
    ! in the spin unpolarised case.                                !
    !    J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)        !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution  by Quintin Hill on 10/03/2009.!
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: s                    ! Dimensionless density gradients
    real(kind=DP) :: fac,ss,s3,s4,f,fs    ! PW91 exchange variables
    real(kind=DP) :: p0,p1,p2,p3,p4,p5,p6

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP

    ! PW91 exchange constants
    real(kind=DP), parameter :: a1 = 0.19645_DP
    real(kind=DP), parameter :: a2 = 0.27430_DP
    real(kind=DP), parameter :: a3 = 0.15084_DP
    real(kind=DP), parameter :: a4 = 100.0_DP
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: a = 7.7956_DP
    real(kind=DP), parameter :: b1 = 0.004_DP

    tol: if (den > dentol) then   ! Positive charge density

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless

       fac = ax*den**THIRD
       ss = s*s
       s3 = ss*s
       s4 = ss*ss
       p0 = 1.0_DP/sqrt(1.0_DP+a*a*ss)
       p1 = log(a*s+1.0_DP/p0)
       p2 = exp(-a4*ss)
       p3 = 1.0_DP/(1.0_DP+a1*s*p1+b1*s4)
       p4 = 1.0_DP+a1*s*p1+(a2-a3*p2)*ss
       f = p3*p4
       x_energy = den * fac * f

       p5 = 2.0_DP*(s*(a2-a3*p2)+a3*a4*s3*p2-2.0_DP*b1*s3)
       p6 = (a1*(p1+a*s*p0)+4.0_DP*b1*s3)*((a2-a3*p2)*ss-b1*s4)
       fs = (p5*p3-p6*p3*p3)
       x_pot = FTHRD*fac*(f-s*fs)

       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only

       fac = ax*den**THIRD
       x_energy = den * fac
       x_pot = FTHRD * fac
       x_dfdmgd = 0.0_DP

    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pw91_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PW91 correlation at a point   !
    ! in the spin unpolarised case.                                !
    !    J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)        !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 10/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation energy per particle
    real(kind=DP) :: ecd,eurs,ec1d        ! PW91 correlation variables
    real(kind=DP) :: bet,delt,pon,b,b2,t2,t4,t6,rs2
    real(kind=DP) :: rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3
    real(kind=DP) :: h0,h1,h,q8,h0t,h0b,h0rs,h1t,ccrs,r1rs
    real(kind=DP) :: h1rs,ht,hrs!,g,h0z,h1z,hz,gz,g3,g4

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI
    real(kind=DP), parameter :: THSVTH = 3.0_DP / 7.0_DP

    ! PW91 correlation parameters
    real(kind=DP), parameter :: a4 = 100.0_DP
    real(kind=DP), parameter :: xnu = 15.75592_dp
    real(kind=DP), parameter :: cc0 = 0.004235_dp
    real(kind=DP), parameter :: cx = -0.001667212_dp
    real(kind=DP), parameter :: alf = 0.09_dp
    real(kind=DP), parameter :: c1 = 0.002568_dp
    real(kind=DP), parameter :: c2 = 0.023266_dp
    real(kind=DP), parameter :: c3 = 0.000007389_dp
    real(kind=DP), parameter :: c4 = 8.723_dp
    real(kind=DP), parameter :: c5 = 0.472_dp
    real(kind=DP), parameter :: c6 = 0.07389_dp

    tol: if (den > dentol) then   ! Positive charge density
       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector

       t = mgd / (2.0_DP*den*ks)    ! gradients

       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,ecunif,eurs)
       ecd = ecunif-THIRD*rs*eurs
       bet = xnu*cc0
       delt = 2.0_DP*alf/bet
       pon = -delt*ecunif/bet
       b = delt/(exp(pon)-1.0_DP)
       b2 = b*b
       t2 = t*t
       t4 = t2*t2
       t6 = t4*t2
       rs2 = rs*rs
       rs3 = rs2*rs
       q4 = 1.0_DP+b*t2
       q5 = 1.0_DP+b*t2+b2*t4
       q6 = c1+c2*rs+c3*rs2
       q7 = 1.0_DP+c4*rs+c5*rs2+c6*rs3
       cc = -cx + q6/q7
       r0 = ks*ks/(kf*kf)
       r1 = a4*r0
       coeff = cc-cc0-THSVTH*cx
       r2 = xnu*coeff
       r3 = exp(-r1*t2)
       h0 = (bet/delt)*log(1.0_DP+delt*q4*t2/q5)
       h1 = r3*r2*t2
       h = h0+h1
       q8 = q5*q5+delt*q4*q5*t2
       h0t = 2.0_DP*bet*t*(1.0_DP+2.0_DP*b*t2)/q8
       h0b = -bet*t6*(2.0_DP*b+b2*t2)/q8
       h0rs = h0b*b*eurs*(b+delt)/bet
       h1t = 2.0_DP*r3*r2*t*(1.0_DP-r1*t2)
       ccrs = (c2+2.0_DP*c3*rs)/q7 - &
            q6*(c4+2.0_DP*c5*rs+3.0_DP*c6*rs2)/(q7*q7)
       r1rs = 100.0_DP*r0/rs
       h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs)
       ht = h0t + h1t
       hrs = h0rs + h1rs
       ec1d = h-THIRD*rs*hrs-SSXTH*t*ht

       c_energy = den * (ecunif + h)
       c_pot = ecd + ec1d
       c_dfdmgd = 0.5_DP * ht / ks

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
       c_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pw91_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the PW91 exchange at a point      !
    ! in the spin polarised case.                                  !
    !    J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)        !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 10/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fac,ss,s3,s4,f,fs    ! PW91 exchange variables
    real(kind=DP) :: p0,p1,p2,p3,p4,p5,p6

    ! Constants
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP

    ! PW91 exchange constants
    real(kind=DP), parameter :: a1 = 0.19645_DP
    real(kind=DP), parameter :: a2 = 0.27430_DP
    real(kind=DP), parameter :: a3 = 0.15084_DP
    real(kind=DP), parameter :: a4 = 100.0_DP
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: a = 7.7956_DP
    real(kind=DP), parameter :: b1 = 0.004_DP

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    tol1: if (den1 > dentol*0.5_DP) then
       ! Spin 1:
       kf1 = (SXPISQ * den1)**THIRD
       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradients
       fac = ax*(2.0_DP*den1)**THIRD
       ss = s1*s1
       s3 = ss*s1
       s4 = ss*ss
       p0 = 1.0_DP/sqrt(1.0_DP+a*a*ss)
       p1 = log(a*s1+1.0_DP/p0)
       p2 = exp(-a4*ss)
       p3 = 1.0_DP/(1.0_DP+a1*s1*p1+b1*s4)
       p4 = 1.0_DP+a1*s1*p1+(a2-a3*p2)*ss
       f = p3*p4
       x_energy = den1 * fac * f
       p5 = 2.0_DP*(s1*(a2-a3*p2)+a3*a4*s3*p2-2.0_DP*b1*s3)
       p6 = (a1*(p1+a*s1*p0)+4.0_DP*b1*s3)*((a2-a3*p2)*ss-b1*s4)
       fs = (p5*p3-p6*p3*p3)
       x_pot1 = FTHRD * fac * (f - s1*fs)
       x_dfdmgd1 = 0.5_DP * ax * fs * TTPI23

    end if tol1

    tol2: if (den2 > 0.5_DP*dentol) then

       ! Spin 2:
       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) ! Dimensionless gradients
       fac = ax*(2.0_DP*den2)**THIRD
       ss = s2*s2
       s3 = ss*s2
       s4 = ss*ss
       p0 = 1.0_DP/sqrt(1.0_DP+a*a*ss)
       p1 = log(a*s2+1.0_DP/p0)
       p2 = exp(-a4*ss)
       p3 = 1.0_DP/(1.0_DP+a1*s2*p1+b1*s4)
       p4 = 1.0_DP+a1*s2*p1+(a2-a3*p2)*ss
       f = p3*p4
       x_energy = x_energy + den2 * fac * f
       p5 = 2.0_DP*(s2*(a2-a3*p2)+a3*a4*s3*p2-2.0_DP*b1*s3)
       p6 = (a1*(p1+a*s2*p0)+4.0_DP*b1*s3)*((a2-a3*p2)*ss-b1*s4)
       fs = (p5*p3-p6*p3*p3)
       x_pot2 = FTHRD * fac * (f - s2*fs)
       x_dfdmgd2 = 0.5_DP * ax * fs * TTPI23

    end if tol2

  end subroutine xc_pw91_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the PW91 correlation at a point   !
    ! in the spin polarised case.                                  !
    !    J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)        !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 10/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad_n_2|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation energy per particle
    real(kind=DP) :: ecd,eu,eurs,ec1d     ! PW91 correlation variables
    real(kind=DP) :: bet,delt,pon,b,b2,t2,t4,t6,rs2
    real(kind=DP) :: rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3
    real(kind=DP) :: h0,h1,h,q8,h0t,h0b,h0rs,h1t,ccrs,r1rs
    real(kind=DP) :: h1rs,ht,hrs,g,h0z,h1z,hz,gz,g3,g4
    real(kind=DP) :: ecz,ecrs,eczet,ec1z
    real(kind=DP) :: ep,eprs,alfm,alfrsm
    real(kind=DP) :: zeta                 ! Spin polarisation
    real(kind=DP) :: zeta3,zeta4          ! zeta**3 and zeta**4
    real(kind=DP) :: fzeta,dfdz           ! Function of zeta and derivative
    real(kind=DP) :: mgd12                ! mgd1 + mgd 2 avoids compiler warning

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI
    real(kind=DP), parameter :: THSVTH = 3.0_DP / 7.0_DP
    real(kind=DP), parameter :: TWTHRD = 2.0_DP / 3.0_DP

    ! PW91 correlation parameters
    real(kind=DP), parameter :: a4 = 100.0_DP
    real(kind=DP), parameter :: xnu = 15.75592_dp
    real(kind=DP), parameter :: cc0 = 0.004235_dp
    real(kind=DP), parameter :: cx = -0.001667212_dp
    real(kind=DP), parameter :: alf = 0.09_dp
    real(kind=DP), parameter :: c1 = 0.002568_dp
    real(kind=DP), parameter :: c2 = 0.023266_dp
    real(kind=DP), parameter :: c3 = 0.000007389_dp
    real(kind=DP), parameter :: c4 = 8.723_dp
    real(kind=DP), parameter :: c5 = 0.472_dp
    real(kind=DP), parameter :: c6 = 0.07389_dp
    real(kind=DP), parameter :: gaminv = 1.92366105093153631981_DP
    real(kind=DP), parameter :: fzzinv = 1.0_DP / (8.0_DP * gaminv / 9.0_DP)
    real(kind=DP), parameter :: zetamin = -1.0_DP + 1.0e-10_DP
    real(kind=DP), parameter :: zetamax = 1.0_DP - 1.0e-10_DP

    ! qoh: Avoid compiler warning by making use of mgd1 and mgd2
    mgd12 = mgd1 + mgd2

    tol: if (den > dentol .and. den1 > 0.0_DP .and. &
         den2 > 0.0_DP) then

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vectors
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector

       zeta = (den1 - den2) / den   ! Polarisation
       zeta = max(zeta,zetamin)
       zeta = min(zeta,zetamax)
       g = 0.5_DP*((1.0_dp-zeta)**TWTHRD+(1.0_DP+zeta)**TWTHRD)

       t = mgd / (den*ks*g*2.0_DP) ! Dimensionless gradient

       ! PW91 correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,eu,eurs)
       call xc_pw91_eq10(0.01554535_DP,0.20548_DP,14.1189_DP, &
            6.1977_DP,3.3662_DP,0.62517_DP,rs,sqrtrs,ep,eprs)
       call xc_pw91_eq10(0.0168869_DP,0.11125_DP,10.357_DP, &
            3.6231_DP,0.88026_DP,0.49671_DP,rs,sqrtrs,alfm,alfrsm)

       fzeta = gaminv * &
            ((1.0_DP+zeta)**FTHRD+(1.0_DP-zeta)**FTHRD-2.0_DP)
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta3*zeta

       ecunif = eu*(1.0_DP-fzeta*zeta4) + ep*fzeta*zeta4 - &
            alfm*fzeta*(1.0_DP-zeta4)*fzzinv
       !JA: why is term 3 negative?
       !JA: because alfm is minus the spin stiffness alfc
       ecrs = eurs*(1.0_DP-fzeta*zeta4) + eprs*fzeta*zeta4 - &
            alfrsm*fzeta*(1.0_DP-zeta4)*fzzinv
       dfdz = gaminv * FTHRD * &
            ((1.0_DP+zeta)**THIRD-(1.0_DP-zeta)**THIRD)
       eczet = (4.0_DP*zeta3*fzeta+zeta4*dfdz) * &
            (ep-eu+alfm*fzzinv)-dfdz*alfm*fzzinv

       ecd = ecunif-THIRD*rs*ecrs
       ecz = eczet * den

       bet = xnu*cc0
       delt = 2.0_DP*alf/bet
       gz = ((1.0_DP+zeta)**(-THIRD)-(1.0_DP-zeta)**(-THIRD))*THIRD
       g3 = g*g*g
       g4 = g3*g
       pon = -delt*ecunif/(g3*bet)
       b = delt/(exp(pon)-1.0_DP)
       b2 = b*b
       t2 = t*t
       t4 = t2*t2
       t6 = t4*t2
       rs2 = rs*rs
       rs3 = rs2*rs
       q4 = 1.0_DP+b*t2
       q5 = 1.0_DP+b*t2+b2*t4
       q6 = c1+c2*rs+c3*rs2
       q7 = 1.0_DP+c4*rs+c5*rs2+c6*rs3
       cc = -cx + q6/q7
       r0 = ks*ks/(kf*kf)
       r1 = a4*r0*g4
       coeff = cc-cc0-THSVTH*cx
       r2 = xnu*coeff*g3
       r3 = exp(-r1*t2)
       h0 = g3*(bet/delt)*log(1.0_DP+delt*q4*t2/q5)
       h1 = r3*r2*t2
       h = h0+h1
       q8 = q5*q5+delt*q4*q5*t2
       h0t = 2.0_DP*bet*t*g3*(1.0_DP+2.0_DP*b*t2)/q8
       h0b = -bet*t6*g3*(2.0_DP*b+b2*t2)/q8
       h0rs = h0b*b*ecrs*(b+delt)/(bet*g3)
       h0z = 3.0_DP*h0*gz/g + &
            h0b*b*(b+delt)*(eczet-3.0_DP*ecunif*gz/g)/(g3*bet)
       h1t = 2.0_DP*r3*r2*t*(1.0_DP-r1*t2)
       ccrs = (c2+2.0_DP*c3*rs)/q7 - &
            q6*(c4+2.0_DP*c5*rs+3.0_DP*c6*rs2)/(q7*q7)
       r1rs = r1/rs
       h1rs = xnu*t2*r3*g3*(ccrs - coeff*t2*r1rs)
       h1z = h1*(3.0_DP-4.0_DP*r1*t2)*gz/g
       ht = h0t + h1t
       hrs = h0rs + h1rs
       hz = h0z + h1z
       ec1d = h-THIRD*rs*hrs-SSXTH*t*ht
       ec1z = den*(hz-ht*t*gz/g)

       c_energy = den * (ecunif + h)
       c_pot1 = ecd+ec1d+(1.0_DP-zeta)*(ecz+ec1z)/den
       c_pot2=  ecd+ec1d-(1.0_DP+zeta)*(ecz+ec1z)/den

       ! qoh: Correlation gradient dependent par
       ! ndmh: separated c_dfdmgd from c_dfdmgd1t
       c_dfdmgd = 0.5_DP * ht / (g*ks)
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    else

       ! qoh: Set all values to zero for very small densities
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
       c_dfdmgd = 0.0_DP
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP

    end if tol

  end subroutine xc_pw91_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_eq10(A,alpha1,beta1,beta2,beta3,beta4,rs,sqrtrs,G,dGdrs)

    !---------------------------------------------------------!
    ! This subroutine evaluates Eq. 10 of:                    !
    !  J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)     !
    ! for G(rs,A,alpha1,beta1,beta2,beta3,beta4,p=1)          !
    ! and its derivative with respect to rs.                  !
    !---------------------------------------------------------!

    implicit none

    !------------------------
    ! INPUT/OUTPUT variables:
    !------------------------

    real(kind=DP), intent(in)  :: A,alpha1,beta1,beta2,beta3,beta4,rs,sqrtrs
    real(kind=DP), intent(out) :: G,dGdrs

    !------------------------
    ! Local variables:
    !------------------------

    real(kind=DP) :: q0,rs32,rs2,q1,q2,q3

    q0 = -2.0_DP*A*(1.0_DP+alpha1*rs)
    rs32 = sqrtrs*rs
    rs2 = rs*rs
    q1 = 2.0_DP*A*(beta1*sqrtrs+beta2*rs+beta3*rs32+beta4*rs2)
    q2 = log(1.0_DP+1.0_DP/q1)
    G = q0*q2
    q3 = A*(beta1/sqrtrs+2.0_DP*beta2+3.0_DP*beta3*sqrtrs+4.0_DP*beta4*rs)
    dGdrs = -2.0_DP*A*alpha1*q2-q0*q3/(q1*q1+q1)

  end subroutine xc_pw91_eq10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b88_x_point(denin,mgdin,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the Becke 88 exchange at a point  !
    ! in the spin unpolarised case.                                !
    !  A.D. Becke                                                  !
    !  Density-functional exchange-energy approximation with       !
    !  correct asymptotic behaviour                                !
    !  Phys. Rev. A38 (1988) 3098-3100                             !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 19/02/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: denin    ! The charge density at a point
    real(kind=dp), intent(in)  :: mgdin    ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: den ! Density
    real(kind=DP) :: xsigma ! x_{sigma} in paper
    real(kind=DP) :: mgd, mgdsq ! Modulus of the the density gradient and square
    real(kind=DP) :: crden, densq ! Cube root of den and den squared
    real(kind=DP) :: crdenbq ! den ^ (4/3)
    real(kind=DP) :: asinhx ! Inverse hyperbolic sine of x
    real(kind=DP) :: dasinhx  ! Derivative wrt x of asinh(x)
    real(kind=DP) :: edenom, edenomsq ! Energy denominator and its square
    real(kind=DP) :: dnedenom ! Energy denominator density derivative
    real(kind=DP) :: dgedenom ! Energy denominator density gradient derivative

    ! qoh: Parameters:
    real(kind=DP), parameter :: beta = 0.0042_DP ! As in paper
    real(kind=DP), parameter :: sixbeta = 6.0_DP * beta
    real(kind=DP), parameter :: ldafac =  -0.73855876638202234_DP
    real(kind=DP), parameter :: ftldafac = fthrd*ldafac
    real(kind=DP), parameter :: crtwo = 1.2599210498948732_DP ! spin factor
    real(kind=DP), parameter :: rcrtwo = 1.0_DP / crtwo ! spin factor

    den = max(0.0_dp,denin)
    tol: if(den.gt.dentol) then
       mgd = max(0.0_dp,mgdin)
       crden = den**THIRD
       crdenbq = crden*den
       mgdsq = mgd**2
       xsigma = crtwo * mgd / crdenbq
       asinhx = log(xsigma+sqrt(1.0_DP+xsigma**2))
       edenom = rcrtwo*crdenbq+sixbeta*mgd*asinhx
       x_energy = ldafac * crdenbq - beta * mgdsq /edenom

       edenomsq = edenom**2
       densq = den**2
       dasinhx = 1.0_DP / sqrt(1.0_dp+xsigma**2)
       dnedenom = fthrd *(rcrtwo * crden - sixbeta* mgdsq * dasinhx * crtwo / &
            (crden * densq))
       dgedenom = sixbeta* mgd * dasinhx * crtwo / crdenbq + sixbeta*asinhx

       x_pot = ftldafac * crden + (beta * mgdsq * dnedenom)/ edenomsq
       x_dfdmgd = - 2.0_DP * beta * mgd / edenom &
            + (beta * mgdsq * dgedenom) / edenomsq
    else
       ! qoh: Set all values to zero for very small densities
       x_energy = 0.0_dp
       x_pot = 0.0_dp
       x_dfdmgd = 0.0_dp
    endif tol

  end subroutine xc_b88_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b88_x_point_sp(den1in,den2in,mgd1in,mgd2in,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the Becke 88 exchange at a point  !
    ! in the spin polarised case.                                  !
    !  A.D. Becke                                                  !
    !  Density-functional exchange-energy approximation with       !
    !  correct asymptotic behaviour                                !
    !  Phys. Rev. A38 (1988) 3098-3100                             !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 23/02/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: den ! Density
    real(kind=DP) :: xsigma ! x_{sigma} in paper
    real(kind=DP) :: mgd, mgdsq ! Modulus of the the density gradient and square
    real(kind=DP) :: crden, densq ! Cube root of den and den squared
    real(kind=DP) :: crdenbq ! den ^ (4/3)
    real(kind=DP) :: asinhx ! Inverse hyperbolic sine of x
    real(kind=DP) :: dasinhx  ! Derivative wrt x of asinh(x)
    real(kind=DP) :: edenom, edenomsq ! Energy denominator and its square
    real(kind=DP) :: dnedenom ! Energy denominator density derivative
    real(kind=DP) :: dgedenom ! Energy denominator density gradient derivative
    ! qoh: Parameters:
    real(kind=DP), parameter :: beta = 0.0042_DP
    real(kind=DP), parameter :: sixbeta = 6.0_DP * beta
    real(kind=DP), parameter :: ldafac = -0.930525736349100025002010_DP
    real(kind=DP), parameter :: ftldafac = fthrd*ldafac

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    den = max(0.0_dp,den1in)
    spin1: if (den.gt.0.5_DP*dentol) then
       crden = den**THIRD
       crdenbq = crden*den
       densq = den**2
       mgd = max(0.0_dp,mgd1in)
       mgdsq = mgd**2
       xsigma = mgd / crdenbq
       asinhx = log(xsigma+sqrt(1.0_DP+xsigma**2))
       edenom = crdenbq+sixbeta*mgd*asinhx
       edenomsq = edenom**2

       x_energy = ldafac * crdenbq - beta * mgdsq /edenom

       dasinhx = 1.0_DP / sqrt(1.0_dp+xsigma**2)
       dnedenom = fthrd *(crden - sixbeta* mgdsq * dasinhx / (crden*densq))
       dgedenom = sixbeta* mgd * dasinhx / crdenbq + sixbeta*asinhx

       x_pot1 = ftldafac * crden + (beta * mgdsq * dnedenom)/ edenomsq
       x_dfdmgd1 = - 2.0_DP * beta * mgd / edenom &
            + (beta * mgdsq * dgedenom) / edenomsq
    end if spin1

    den = max(0.0_dp,den2in)
    spin2: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd2in)
       crden = den**THIRD
       crdenbq = crden*den
       mgdsq = mgd**2
       xsigma = mgd / crdenbq
       asinhx = log(xsigma+sqrt(1.0_DP+xsigma**2))
       edenom = crdenbq+sixbeta*mgd*asinhx

       !qoh: Add to x_energy from spin 1
       x_energy = x_energy + ldafac * crdenbq - beta * mgdsq /edenom

       edenomsq = edenom**2
       densq = den**2
       dasinhx = 1.0_DP / sqrt(1.0_dp+xsigma**2)
       dnedenom = fthrd *(crden - sixbeta* mgdsq * dasinhx / (crden*densq))
       dgedenom = sixbeta* mgd * dasinhx / crdenbq + sixbeta*asinhx

       x_pot2 = ftldafac * crden + (beta * mgdsq * dnedenom)/ edenomsq
       x_dfdmgd2 = - 2.0_DP * beta * mgd / edenom &
            + (beta * mgdsq * dgedenom) / edenomsq
    end if spin2

  end subroutine xc_b88_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b1_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !===================================================================!
    ! This subroutine calculates the B1{LYP,PW91} exchange at a         !
    ! point in the spin unpolarised case.                               !
    ! B1{LYP,PW91}_X = 1/4 HF_X + 3/4 B88_X                             !
    !  C. Adamo and V. Barone, Chemical Physics Letters 274, 242 (1997) !
    !-------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                            !
    !===================================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP), parameter :: b88_fac = 0.75_DP ! B1{LYP,PW91} a_x

    call xc_b88_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)
    x_energy = x_energy*b88_fac
    x_pot = x_pot * b88_fac
    x_dfdmgd = x_dfdmgd * b88_fac

  end subroutine xc_b1_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==================================================================!
    ! This subroutine calculates the B3{LYP,PW91} exchange at a point  !
    ! in the spin unpolarised case.                                    !
    ! B3{LYP,PW91} = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X     !
    !                + (1 - a_c) VWN_C + a_c {LYP,PW91}_C              !
    ! a_0 = 0.2, a_x = 0.72, a_c = 0.81
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    ! P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and             !
    !   M. J. Frisch, J. Phys. Chem. 98, 11623 (1994)                  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot

    real(kind=DP), parameter :: b88_fac = 0.72_DP ! B3{LYP,PW91} a_x
    real(kind=DP), parameter :: lda_fac = 0.08_DP ! B3{LYP,PW91} 1 - a_0 - a_x

    call xc_lda_x_point(den,lda_energy,lda_pot)
    call xc_b88_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)
    x_energy = x_energy*b88_fac + lda_energy*lda_fac
    x_pot = x_pot * b88_fac + lda_pot * lda_fac
    x_dfdmgd = x_dfdmgd * b88_fac

  end subroutine xc_b3_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the X exchange at a point         !
    ! in the spin unpolarised case.                                !
    ! X = 0.722 B88_X + 0.347 PW91_X - 0.069 LDA_X                 !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot
    real(kind=DP) :: pw91_energy
    real(kind=DP) :: pw91_pot
    real(kind=DP) :: pw91_dfdmgd

    real(kind=DP), parameter :: b88_fac  =  0.722_DP
    real(kind=DP), parameter :: pw91_fac =  0.347_DP
    real(kind=DP), parameter :: lda_fac  = -0.069_DP

    call xc_lda_x_point(den,lda_energy,lda_pot)
    call xc_pw91_x_point(den,mgd,pw91_energy,pw91_pot,pw91_dfdmgd)
    call xc_b88_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    x_energy = x_energy*b88_fac + lda_energy*lda_fac + pw91_energy * pw91_fac
    x_pot = x_pot * b88_fac + lda_pot * lda_fac + pw91_pot * pw91_fac
    x_dfdmgd = x_dfdmgd * b88_fac + pw91_dfdmgd * pw91_fac

   end subroutine xc_x_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x3_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the X3LYP exchange at a point     !
    ! in the spin unpolarised case.                                !
    ! X3LYP = a_0 HF_X + (1 - a_0 - a_x) LDA_X                     !
    !         + a_x ( 0.765 B88_X + 0.235 PW91_X) + a_c VWN_C      !
    !         + (1 - a_C) LYP_C                                    !
    ! a_0 = 0.218, a_x = 0.709, a_c = 0.129.                       !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot
    real(kind=DP) :: pw91_energy
    real(kind=DP) :: pw91_pot
    real(kind=DP) :: pw91_dfdmgd

    real(kind=DP), parameter :: b88_fac = 0.765_DP * 0.709_DP
    real(kind=DP), parameter :: pw91_fac = 0.235_DP * 0.709_DP
    real(kind=DP), parameter :: lda_fac = 0.073_DP

    call xc_lda_x_point(den,lda_energy,lda_pot)
    call xc_pw91_x_point(den,mgd,pw91_energy,pw91_pot,pw91_dfdmgd)
    call xc_b88_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    x_energy = x_energy*b88_fac + lda_energy*lda_fac + pw91_energy * pw91_fac
    x_pot = x_pot * b88_fac + lda_pot * lda_fac + pw91_pot * pw91_fac
    x_dfdmgd = x_dfdmgd * b88_fac + pw91_dfdmgd * pw91_fac

  end subroutine xc_x3_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b1_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !===================================================================!
    ! This subroutine calculates the B1{LYP,PW91} exchange at a         !
    ! point in the spin polarised case.                                 !
    ! B1{LYP,PW91}_X = 1/4 HF_X + 3/4 B88_X                             !
    !  C. Adamo and V. Barone, Chemical Physics Letters 274, 242 (1997) !
    !-------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                            !
    !===================================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP), parameter :: b88_fac = 0.75_DP ! B1{LYP,PW91} a_x

    call xc_b88_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)
    x_energy = x_energy*b88_fac
    x_pot1 = x_pot1 * b88_fac
    x_pot2 = x_pot2 * b88_fac
    x_dfdmgd1 = x_dfdmgd1 * b88_fac
    x_dfdmgd2 = x_dfdmgd2 * b88_fac

  end subroutine xc_b1_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==================================================================!
    ! This subroutine calculates the B3{LYP,PW91} exchange at a point  !
    ! in the spin polarised case.                                      !
    ! B3{LYP,PW91} = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X     !
    !                + (1 - a_c) VWN_C + a_c {LYP,PW91}_C              !
    ! a_0 = 0.2, a_x = 0.72, a_c = 0.81                                !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    ! P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and             !
    !   M. J. Frisch, J. Phys. Chem. 98, 11623 (1994)                  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2

    real(kind=DP), parameter :: b88_fac = 0.72_DP ! B3{LYP,PW91} a_x
    real(kind=DP), parameter :: lda_fac = 0.08_DP ! B3{LYP,PW91} 1 - a_0 - a_x

    call xc_lsda_x_point(den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_b88_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)
    x_energy = x_energy*b88_fac + lda_energy*lda_fac
    x_pot1 = x_pot1 * b88_fac + lda_pot1 * lda_fac
    x_pot2 = x_pot2 * b88_fac + lda_pot2 * lda_fac
    x_dfdmgd1 = x_dfdmgd1 * b88_fac
    x_dfdmgd2 = x_dfdmgd2 * b88_fac

  end subroutine xc_b3_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the X exchange at a point         !
    ! in the spin polarised case.                                  !
    ! X = 0.722 B88_X + 0.347 PW91_X - 0.069 * LDA_X               !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2
    real(kind=DP) :: pw91_energy
    real(kind=DP) :: pw91_pot1
    real(kind=DP) :: pw91_pot2
    real(kind=DP) :: pw91_dfdmgd1
    real(kind=DP) :: pw91_dfdmgd2

    real(kind=DP), parameter :: b88_fac  =  0.722_DP
    real(kind=DP), parameter :: pw91_fac =  0.347_DP
    real(kind=DP), parameter :: lda_fac  = -0.069_DP

    call xc_lsda_x_point(den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_pw91_x_point_sp(den1,den2,mgd1,mgd2,&
         pw91_energy,pw91_pot1,pw91_pot2,pw91_dfdmgd1,pw91_dfdmgd2)
    call xc_b88_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    x_energy = x_energy*b88_fac + lda_energy*lda_fac + pw91_energy * pw91_fac
    x_pot1 = x_pot1 * b88_fac + lda_pot1 * lda_fac + pw91_pot1 * pw91_fac
    x_pot2 = x_pot2 * b88_fac + lda_pot2 * lda_fac + pw91_pot2 * pw91_fac
    x_dfdmgd1 = x_dfdmgd1 * b88_fac + pw91_dfdmgd1 * pw91_fac
    x_dfdmgd2 = x_dfdmgd2 * b88_fac + pw91_dfdmgd2 * pw91_fac

  end subroutine xc_x_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x3_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the X3LYP exchange at a point     !
    ! in the spin polarised case.                                  !
    ! X3LYP = a_0 HF_X + (1 - a_0 - a_x) LDA_X                     !
    !         + a_x ( 0.765 B88_X + 0.235 PW91_X) + a_c VWN_C      !
    !         + (1 - a_C) LYP_C                                    !
    ! a_0 = 0.218, a_x = 0.709, a_c = 0.129.                       !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2
    real(kind=DP) :: pw91_energy
    real(kind=DP) :: pw91_pot1
    real(kind=DP) :: pw91_pot2
    real(kind=DP) :: pw91_dfdmgd1
    real(kind=DP) :: pw91_dfdmgd2

    real(kind=DP), parameter :: b88_fac = 0.765_DP * 0.709_DP
    real(kind=DP), parameter :: pw91_fac = 0.235_DP * 0.709_DP
    real(kind=DP), parameter :: lda_fac = 0.073_DP

    call xc_lsda_x_point(den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_pw91_x_point_sp(den1,den2,mgd1,mgd2,&
         pw91_energy,pw91_pot1,pw91_pot2,pw91_dfdmgd1,pw91_dfdmgd2)
    call xc_b88_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    x_energy = x_energy*b88_fac + lda_energy*lda_fac + pw91_energy * pw91_fac
    x_pot1 = x_pot1 * b88_fac + lda_pot1 * lda_fac + pw91_pot1 * pw91_fac
    x_pot2 = x_pot2 * b88_fac + lda_pot2 * lda_fac + pw91_pot2 * pw91_fac
    x_dfdmgd1 = x_dfdmgd1 * b88_fac + pw91_dfdmgd1 * pw91_fac
    x_dfdmgd2 = x_dfdmgd2 * b88_fac + pw91_dfdmgd2 * pw91_fac

  end subroutine xc_x3_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lyp_c_point(denin,mgdin,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the LYP correlation at a point    !
    ! in the spin unpolarised case.                                !
    !     C. Lee, W. Yang, and R.G. Parr                           !
    !     Development of the Colle-Salvetti correlation-energy     !
    !     formula into a functional of the electron density        !
    !     Phys. Rev. B37 (1988) 785-789                            !
    !  Using reformulation (to avoid Laplacian) given by:          !
    !     B. Miehlich, A. Savin, H. Stoll and H. Preuss            !
    !     Results obtained with the correlation energy density     !
    !     functionals of becke and Lee, Yang and Parr              !
    !     Chem. Phys. Lett. 157 (1989) 200-206                     !
    !  Note: There is a mistake in Handy's formula in ESQC book.   !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 24/02/2009                        !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: denin    ! The charge density at a point
    real(kind=dp), intent(in)  :: mgdin    ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: den, rden ! Density and its reciprocal
    real(kind=DP) :: mgd, mgdsq ! Modulus of density gradient and it's square
    real(kind=DP) :: crden ! Cube root of den
    real(kind=DP) :: rcrden ! Reciprocal cube root of den
    real(kind=DP) :: densq, rdensq ! den squared and its reciprocal
    !real(kind=DP) :: rcrdensq ! Reciprocal cube root of den squared
    real(kind=DP) :: redenom, redenomsq ! Reciprocal of E denominator and square
    real(kind=DP) :: expo ! Exponential term
    real(kind=DP) :: omega ! As in paper
    real(kind=DP) :: delta ! As in paper
    real(kind=DP) :: bracket
    real(kind=DP) :: dndelta ! Density derivative of density
    real(kind=DP) :: dnomega ! Density derivative of omega

    ! qoh: Parameters:
    real(kind=DP), parameter :: aa = 0.04918_DP ! As in paper
    real(kind=DP), parameter :: bb = 0.132_DP   ! As in paper
    real(kind=DP), parameter :: ab = aa * bb
    real(kind=DP), parameter :: cc = 0.25330_DP ! As in paper
    real(kind=DP), parameter :: dd = 0.3490_DP  ! As in paper
    real(kind=DP), parameter :: ETHRD = 8.0_DP / 3.0_DP
    real(kind=DP), parameter :: cf = 2.8712340001881911_DP ! 0.3*(3*PI^2)^(2/3)
    real(kind=DP), parameter :: se72 = 7.0_DP/72.0_DP
    real(kind=DP), parameter :: o24 = 1.0_DP/24.0_DP

    den = max(0.0_dp,denin)
    tol: if(den.gt.dentol) then
       ! qoh: Calculate useful powers of den
       densq = den**2
       rdensq = 1.0_DP / densq
       rden = 1.0_DP / den
       crden = den**THIRD
       rcrden = 1/crden
       !rcrdensq = rcrden**2
       ! qoh: Initialise density gradient and square
       mgd = max(0.0_dp,mgdin)
       mgdsq = mgd**2

       redenom = 1/( 1.0_dp+dd*rcrden )
       redenomsq = redenom**2
       expo = exp(-cc*rcrden)

       omega = expo * crden * redenom * rdensq * rdensq
       delta = cc*rcrden + dd*rcrden * redenom

       bracket = CF*densq*den*rcrden - mgdsq*(o24+se72*delta)
       dnomega = ( THIRD * rdensq * rdensq * redenom * expo * rden ) * &
            (cc - 11.0_DP*crden + dd * redenom )
       dndelta = - THIRD * rden * rcrden * ( cc + dd * redenomsq )

       c_energy = - aa*den*redenom - ab * omega*densq * bracket

       c_pot = - aa * ( redenom + THIRD*rcrden*dd*redenomsq ) &
            - ab * ( dnomega*densq + 2.0_DP*omega*den) * bracket &
            - ab * omega*densq  * (ETHRD*CF*densq*rcrden &
            - mgdsq * (se72*dndelta))

       !qoh: Note: No minus sign since we have merged brackets
       c_dfdmgd = ab * omega * densq * mgd * 2.0_DP * ( o24 + se72 * delta)
    else
       ! qoh: Set all values to zero for very small densities
       c_energy = 0.00_DP
       c_pot = 0.0_DP
       c_dfdmgd = 0.0_DP
    endif tol

  end subroutine xc_lyp_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_rpw86_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the refitted PW86 exchange        !
    ! at a point in the spin unpolarised case.                     !
    !  Refitted version of PW86 by Perdew et al                    !
    !  E. Murray et al.                                            !
    !  J. Chem. Theory Comput. 5, 2754 (2009)                      !
    !--------------------------------------------------------------!
    ! Written by Lampros Andrinopoulos in June 2013                !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FOURTEEN = 14.0_DP
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP

    ! rPW86 enhancement factor parameters
    real(kind=DP), parameter :: a = 0.1234_DP
    real(kind=DP), parameter :: b = 17.3300_DP
    real(kind=DP), parameter :: c = 0.1630_DP

    real(kind=DP) :: s2,s3,s4,s5,s6       ! Powers of s
    real(kind=DP) :: dfx_ds               ! Derivatives for chain rule

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient
       s2 = s*s
       s3 = s2*s
       s4 = s2*s2
       s5 = s4*s
       s6 = s4*s2

       exunif = ax*den**THIRD
       fx = (1.0_DP + 15.0_DP*a*s2 + b*s4 + c*s6)**(1.0_DP/15.0_DP)

       x_energy = den * exunif * fx

       dfx_ds = (30.0_DP*a*s+4.0_DP*b*s3+6.0_DP*c*s5)*fx**(-FOURTEEN)/15.0_DP

       x_pot = FTHRD*ax*(fx*den**THIRD - dfx_ds*mgd/(2.0_DP*den*(THPISQ**THIRD)))
       x_dfdmgd = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_rpw86_x_point

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_rpw86_x_point_sp(den1in,den2in,mgd1in,mgd2in,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the refitted PW86 exchange        !
    ! at a point in the spin unpolarised case. Spin polarized.     !
    !  Refitted version of PW86 by Perdew et al                    !
    !  E. Murray et al.                                            !
    !  J. Chem. Theory Comput. 5, 2754 (2009)                      !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons on 21/05/2021.                      !
    !==============================================================!

    use constants, only: DP, Pi
    implicit none

    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: den ! Density
    real(kind=DP) :: mgd
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: FOURTEEN = 14.0_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP

    ! rPW86 enhancement factor parameters
    real(kind=DP), parameter :: a = 0.1234_DP
    real(kind=DP), parameter :: b = 17.3300_DP
    real(kind=DP), parameter :: c = 0.1630_DP

    real(kind=DP) :: s2,s3,s4,s5,s6       ! Powers of s
    real(kind=DP) :: dfx_ds               ! Derivatives for chain rule

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    den = max(0.0_dp,den1in)
    spin1: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd1in)

       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       s3 = s2*s
       s4 = s2*s2
       s5 = s4*s
       s6 = s4*s2

       exunif = ax*(2.0_dp*den)**THIRD
       fx = (1.0_DP + 15.0_DP*a*s2 + b*s4 + c*s6)**(1.0_DP/15.0_DP)

       x_energy = x_energy + den * exunif * fx

       dfx_ds = (30.0_DP*a*s+4.0_DP*b*s3+6.0_DP*c*s5)*fx**(-FOURTEEN)/15.0_DP

       x_pot1 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*(THPISQ**THIRD)))
       x_dfdmgd1 = ax * dfx_ds/(2.0*THPISQ**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy = x_energy + den * exunif
       x_pot1 = FTHRD * exunif
       x_dfdmgd1 = 0.0_DP
    end if spin1

    den = max(0.0_dp,den2in)
    spin2: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd2in)

       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       s3 = s2*s
       s4 = s2*s2
       s5 = s4*s
       s6 = s4*s2

       exunif = ax*(2.0_dp*den)**THIRD
       fx = (1.0_DP + 15.0_DP*a*s2 + b*s4 + c*s6)**(1.0_DP/15.0_DP)

       x_energy = x_energy + den * exunif * fx

       dfx_ds = (30.0_DP*a*s+4.0_DP*b*s3+6.0_DP*c*s5)*fx**(-FOURTEEN)/15.0_DP

       x_pot2 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*(THPISQ**THIRD)))
       x_dfdmgd2 = ax * dfx_ds/(2.0*THPISQ**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy = x_energy + den * exunif
       x_pot2 = FTHRD * exunif
       x_dfdmgd2 = 0.0_DP
    end if spin2

  end subroutine xc_rpw86_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_c09_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the C09x exchange at a point      !
    ! in the spin unpolarised case.                                !
    !  Exchange to be used with vdw-DF1 or vdW-DF2                 !
    !  Valentino R. Cooper                                         !
    !  Phys. Rev. B 81 161104(R) (2010)                            !
    !--------------------------------------------------------------!
    ! Written by Lampros Andrinopoulos in June 2013                !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP

    ! C09x enhancement factor parameters
    real(kind=DP), parameter :: mu = 0.0617_DP
    real(kind=DP), parameter :: kappa = 1.245_DP
    real(kind=DP), parameter :: alpha = 0.0483_DP

    real(kind=DP) :: s2,s3,s4,s5,s6       ! Powers of s
    real(kind=DP) :: exp1, exp2, dfx_ds

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient
       s2 = s*s
       s3 = s2*s
       s4 = s2*s2
       s5 = s4*s
       s6 = s4*s2

       exp1 = exp(-alpha*s2)
       exp2 = exp(-0.5_DP*alpha*s2)

       exunif = ax*den**THIRD
       fx = 1.0_DP + mu*s2*exp1 + kappa*(1.0_DP - exp2)

       x_energy = den * exunif * fx

       dfx_ds = 2.0_DP*mu*exp1*(s-alpha*s3)+alpha*s*kappa*exp2

       x_pot = FTHRD*ax*(fx*den**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_c09_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_c09_x_point_sp(den1in,den2in,mgd1in,mgd2in,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the C09x exchange at a point      !
    ! in the spin polarised case.                                  !
    !  Exchange to be used with vdw-DF1 or vdW-DF2                 !
    !  Valentino R. Cooper                                         !
    !  Phys. Rev. B 81 161104(R) (2010)                            !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons on 25/05/2021.                      !
    !==============================================================!

    use constants, only: DP, Pi
    implicit none

    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: den ! Density
    real(kind=DP) :: mgd
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: FOURTEEN = 14.0_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP

    ! C09x enhancement factor parameters
    real(kind=DP), parameter :: mu = 0.0617_DP
    real(kind=DP), parameter :: kappa = 1.245_DP
    real(kind=DP), parameter :: alpha = 0.0483_DP

    real(kind=DP) :: s2,s3,s4,s5,s6       ! Powers of s
    real(kind=DP) :: exp1, exp2
    real(kind=DP) :: dfx_ds               ! Derivatives for chain rule

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    den = max(0.0_dp,den1in)
    spin1: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd1in)

       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       s3 = s2*s
       s4 = s2*s2
       s5 = s4*s
       s6 = s4*s2

       exunif = ax*(2.0_dp*den)**THIRD

       exp1 = exp(-alpha*s2)
       exp2 = exp(-0.5_DP*alpha*s2)

       fx = 1.0_DP + mu*s2*exp1 + kappa*(1.0_DP - exp2)

       x_energy = x_energy + den * exunif * fx

       dfx_ds = 2.0_DP*mu*exp1*(s-alpha*s3)+alpha*s*kappa*exp2

       x_pot1 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd1 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    end if spin1

    den = max(0.0_dp,den2in)
    spin2: if (den.gt.0.5_DP*dentol) then
              mgd = max(0.0_dp,mgd2in)

       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       s3 = s2*s
       s4 = s2*s2
       s5 = s4*s
       s6 = s4*s2

       exunif = ax*(2.0_dp*den)**THIRD

       exp1 = exp(-alpha*s2)
       exp2 = exp(-0.5_DP*alpha*s2)

       fx = 1.0_DP + mu*s2*exp1 + kappa*(1.0_DP - exp2)

       x_energy = x_energy + den * exunif * fx

       dfx_ds = 2.0_DP*mu*exp1*(s-alpha*s3)+alpha*s*kappa*exp2

       x_pot2 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd2 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    end if spin2

  end subroutine xc_c09_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vdwoptb88_x_point(denin,mgdin,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the vdW-optB88                    !
    ! exchange at a point in the spin unpolarised case.            !
    !  Optimized B88 exchange to be used with vdw-DF1              !
    !  Klimes et al                                                !
    !  J. Phys.: Condens. Matter 22 022201 (2010)                  !
    !--------------------------------------------------------------!
    ! Written by Lampros Andrinopoulos in June 2013                !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: denin      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgdin      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: den, mgd

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP

    ! vdW-optB88 enhancement factor parameters
    real(kind=DP), parameter :: mu = 0.22_DP
    real(kind=DP), parameter :: beta = mu/1.2_DP

    real(kind=DP) :: s2                   ! Powers of s
    real(kind=DP) :: x, x2, c, dfx_ds, denom, temp, root

    den = max(0.0_DP,denin)
    tol: if (den > dentol) then

       mgd = max(0.0_DP,mgdin)
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient
       s2 = s*s

       c = (2.0_DP)**(FTHRD)*(THPISQ)**THIRD
       !c = (2.0_DP)*(THPISQ)**THIRD
       x = c*s
       x2 = x*x
       root = sqrt(1.0_DP+x2)
       temp = x+root
       denom = 1.0_DP + beta*s*log(temp)

       exunif = ax*den**THIRD
       fx = 1.0_DP + mu*s2/denom

       x_energy = den * exunif * fx

       dfx_ds = 2.0_DP*mu*s/denom - &
            mu*s2*(beta*log(temp)+beta*s*c/root)/(denom*denom)

       x_pot = FTHRD*ax*(fx*den**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_vdwoptb88_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vdwoptb88_x_point_sp(den1in,den2in,mgd1in,mgd2in,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the vdW-optB88                    !
    ! exchange at a point in the spin polarised case.              !
    !  Optimized B88 exchange to be used with vdw-DF1              !
    !  Klimes et al                                                !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons on 25/05/2021.                      !
    !==============================================================!

    use constants, only: DP, Pi
    implicit none

    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: den ! Density
    real(kind=DP) :: mgd
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: FOURTEEN = 14.0_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP

    ! vdW-optB88 enhancement factor parameters
    real(kind=DP), parameter :: mu = 0.22_DP
    real(kind=DP), parameter :: beta = mu/1.2_DP

    real(kind=DP) :: s2                   ! Powers of s
    real(kind=DP) :: x, x2, c, denom, temp, root
    real(kind=DP) :: dfx_ds               ! Derivatives for chain rule

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    den = max(0.0_dp,den1in)
    spin1: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd1in)

       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s

       c = (2.0_DP)**(FTHRD)*(THPISQ)**THIRD
       !c = (2.0_DP)*(THPISQ)**THIRD
       x = c*s
       x2 = x*x
       root = sqrt(1.0_DP+x2)
       temp = x+root
       denom = 1.0_DP + beta*s*log(temp)

       exunif = ax*(2.0_dp*den)**THIRD
       fx = 1.0_DP + mu*s2/denom

       x_energy = x_energy + den * exunif * fx

       dfx_ds = 2.0_DP*mu*s/denom - &
            mu*s2*(beta*log(temp)+beta*s*c/root)/(denom*denom)

       x_pot1 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd1 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy =  x_energy + den * exunif
       x_pot1 = FTHRD * exunif
       x_dfdmgd1 = 0.0_DP
    end if spin1

    den = max(0.0_dp,den2in)
    spin2: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd2in)

       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s

       c = (2.0_DP)**(FTHRD)*(THPISQ)**THIRD
       !c = (2.0_DP)*(THPISQ)**THIRD
       x = c*s
       x2 = x*x
       root = sqrt(1.0_DP+x2)
       temp = x+root
       denom = 1.0_DP + beta*s*log(temp)

       exunif = ax*(2.0_dp*den)**THIRD
       fx = 1.0_DP + mu*s2/denom

       x_energy = x_energy + den * exunif * fx

       dfx_ds = 2.0_DP*mu*s/denom - &
            mu*s2*(beta*log(temp)+beta*s*c/root)/(denom*denom)

       x_pot2 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd2 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy =  x_energy + den * exunif
       x_pot2 = FTHRD * exunif
       x_dfdmgd2 = 0.0_DP
    end if spin2

  end subroutine xc_vdwoptb88_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vdwoptb86b_x_point(denin,mgdin,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the vdW-optB86b                   !
    ! exchange at a point in the spin unpolarised case.            !
    !  Optimized B88 exchange to be used with vdw-DF1              !
    !  Klimes et al                                                !
    !  Phys. Rev. B 83 195131 (2011)                               !
    !--------------------------------------------------------------!
    ! Written by Lampros Andrinopoulos in June 2013                !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: denin      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgdin      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: den, mgd

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: FIFTH = 0.2_DP

    ! vdW-optB88 enhancement factor parameters
    real(kind=DP), parameter :: mu = 0.1234_DP

    real(kind=DP) :: s2                   ! Powers of s
    real(kind=DP) ::  dfx_ds, denom, temp

    den = max(0.0_DP,denin)
    tol: if (den > dentol) then

       mgd = max(0.0_DP,mgdin)
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient
       s2 = s*s


       exunif = ax*den**THIRD
       temp = (1.0_DP+mu*s2)
       denom = temp**(4.0_DP*FIFTH)
       fx = 1.0_DP + mu*s2/denom

       x_energy = den * exunif * fx

       dfx_ds = 2.0_DP*mu*s*(1.0_DP/denom &
            -4.0_DP*FIFTH*mu*s2*temp**(-9.0_DP*FIFTH))

       x_pot = FTHRD*ax*(fx*den**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_vdwoptb86b_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vdwoptb86b_x_point_sp(den1in,den2in,mgd1in,mgd2in,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the vdW-optB86b                   !
    ! exchange at a point in the spin polarised case.              !
    !  Optimized B88 exchange to be used with vdw-DF1              !
    !  Klimes et al                                                !
    !  Phys. Rev. B 83 195131 (2011)                               !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons on 25/05/2021.                      !
    !==============================================================!

    use constants, only: DP, Pi
    implicit none

    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: den ! Density
    real(kind=DP) :: mgd
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: FIFTH = 0.2_DP
    real(kind=DP), parameter :: FOURTEEN = 14.0_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP

    ! vdW-optB88 enhancement factor parameters
    real(kind=DP), parameter :: mu = 0.1234_DP
    real(kind=DP) :: s2                   ! Powers of s
    real(kind=DP) :: x, x2, c, denom, temp, root
    real(kind=DP) :: dfx_ds               ! Derivatives for chain rule

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    den = max(0.0_dp,den1in)
    spin1: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd1in)

       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s

       exunif = ax*(2.0_dp*den)**THIRD
       temp = (1.0_DP+mu*s2)
       denom = temp**(4.0_DP*FIFTH)
       fx = 1.0_DP + mu*s2/denom

       x_energy = x_energy + den * exunif * fx

       dfx_ds = 2.0_DP*mu*s*(1.0_DP/denom &
            -4.0_DP*FIFTH*mu*s2*temp**(-9.0_DP*FIFTH))

       x_pot1 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd1 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy =  x_energy + den * exunif
       x_pot1 = FTHRD * exunif
       x_dfdmgd1 = 0.0_DP
    end if spin1

    den = max(0.0_dp,den2in)
    spin2: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd2in)

       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s

       exunif = ax*(2.0_dp*den)**THIRD

       temp = (1.0_DP+mu*s2)
       denom = temp**(4.0_DP*FIFTH)
       fx = 1.0_DP + mu*s2/denom

       x_energy = x_energy + den * exunif * fx

       dfx_ds = 2.0_DP*mu*s*(1.0_DP/denom &
            -4.0_DP*FIFTH*mu*s2*temp**(-9.0_DP*FIFTH))

       x_pot2 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd2 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy =  x_energy + den * exunif
       x_pot2 = FTHRD * exunif
       x_dfdmgd2 = 0.0_DP
    end if spin2

  end subroutine xc_vdwoptb86b_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vdwpbek1_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the vdW-PBE_kappa=1               !
    ! exchange at a point in the spin unpolarised case.            !
    !  Optimized PBE_kappa=1 exchange to be used with vdw-DF1      !
    !  Klimes et al                                                !
    !  J. Phys.: Condens. Matter 22 022201 (2010)                  !
    !--------------------------------------------------------------!
    ! Written by Lampros Andrinopoulos in June 2013                !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! PBE_kappa=1 enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: mu = 0.2195149727645171_DP
    real(kind=DP), parameter :: kappa = 1.0_DP

    ! Local variables
    real(kind=DP) :: s2, dfx_ds

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       exunif = ax*den**THIRD

       fx = 1.0_DP + kappa - kappa/(1.0_DP + mu*s2/kappa)
       dfx_ds = 2.0_DP*mu*s/((1.0_DP + mu*s2/kappa)**2.0_DP)

       x_energy =  den * exunif * fx
       x_pot = FTHRD*ax*(fx*den**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_vdwpbek1_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vdwpbek1_x_point_sp(den1in,den2in,mgd1in,mgd2in,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the vdW-PBE_kappa=1               !
    ! exchange at a point in the spin polarised case.              !
    !  Optimized PBE_kappa=1 exchange to be used with vdw-DF1      !
    !  Klimes et al                                                !
    !  J. Phys.: Condens. Matter 22 022201 (2010)                  !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons on 25/05/2021.                      !
    !==============================================================!

    use constants, only: DP, Pi
    implicit none

    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: den ! Density
    real(kind=DP) :: mgd
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: FIFTH = 0.2_DP
    real(kind=DP), parameter :: FOURTEEN = 14.0_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI

    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: mu = 0.2195149727645171_DP
    real(kind=DP), parameter :: kappa = 1.0_DP
    ! Local variables
    real(kind=DP) :: s2, dfx_ds

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    den = max(0.0_dp,den1in)
    spin1: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd1in)
       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       exunif = ax*(2.0_dp*den)**THIRD

       fx = 1.0_DP + kappa - kappa/(1.0_DP + mu*s2/kappa)
       dfx_ds = 2.0_DP*mu*s/((1.0_DP + mu*s2/kappa)**2.0_DP)

       x_energy =  x_energy + den * exunif * fx
       x_pot1 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd1 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy =  den * exunif
       x_pot1 = FTHRD * exunif
       x_dfdmgd1 = 0.0_DP
    end if spin1

    den = max(0.0_dp,den2in)
    spin2: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd1in)
       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       exunif = ax*(2.0_dp*den)**THIRD

       fx = 1.0_DP + kappa - kappa/(1.0_DP + mu*s2/kappa)
       dfx_ds = 2.0_DP*mu*s/((1.0_DP + mu*s2/kappa)**2.0_DP)

       x_energy =  x_energy + den * exunif * fx
       x_pot2 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd2 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy =  x_energy + den * exunif
       x_pot2 = FTHRD * exunif
       x_dfdmgd2 = 0.0_DP
    end if spin2

  end subroutine xc_vdwpbek1_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vdwoptpbe_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the vdW-optPBE                    !
    ! exchange at a point in the spin unpolarised case.            !
    !  Optmized mixing of PBE exchange with RPBE exchange          !
    !  to be used with vdw-DF1                                     !
    !  Klimes et al                                                !
    !  J. Phys.: Condens. Matter 22 022201 (2010)                  !
    !--------------------------------------------------------------!
    ! Written by Lampros Andrinopoulos in June 2013                !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! optPBE enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: mu = 0.175519_DP
    real(kind=DP), parameter :: kappa = 1.04804_DP
    ! Mixing ratio
    real(kind=DP), parameter :: x = 0.945268_DP

    ! Local variables
    real(kind=DP) :: s2, dfx1_ds, dfx2_ds, dfx_ds, fx1,fx2

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       exunif = ax*den**THIRD

       ! Mixing of PBE and RPBE with ratio x, with same mu and kappa

       ! PBE enhancement factor
       fx1 = 1.0_DP + kappa - kappa/(1.0_DP + mu*s2/kappa)
       ! RPBE enhancement factor
       fx2 = 1.0_DP + kappa*(1.0_DP-exp(-mu*s2/kappa))

       ! Derivatives
       dfx1_ds = 2.0_DP*mu*s/((1.0_DP + mu*s2/kappa)**2.0_DP)
       dfx2_ds = 2.0_DP*mu*s*exp(-mu*s2/kappa)

       fx = fx1*x + (1.0_DP - x)*fx2
       dfx_ds = dfx1_ds*x + (1.0_DP -x)*dfx2_ds

       x_energy =  den * exunif * fx
       x_pot = FTHRD*ax*(fx*den**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_vdwoptpbe_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vdwoptpbe_x_point_sp(den1in,den2in,mgd1in,mgd2in,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the vdW-optPBE                    !
    ! exchange at a point in the spin polarised case.              !
    !  Optmized mixing of PBE exchange with RPBE exchange          !
    !  to be used with vdw-DF1                                     !
    !  Klimes et al                                                !
    !  J. Phys.: Condens. Matter 22 022201 (2010)                  !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons on 25/05/2021.                      !
    !==============================================================!

    use constants, only: DP, Pi
    implicit none

    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: den ! Density
    real(kind=DP) :: mgd
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! Enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: FIFTH = 0.2_DP
    real(kind=DP), parameter :: FOURTEEN = 14.0_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: mu = 0.175519_DP
    real(kind=DP), parameter :: kappa = 1.04804_DP
    ! Mixing ratio
    real(kind=DP), parameter :: x = 0.945268_DP

    ! Local variables
    real(kind=DP) :: s2, dfx1_ds, dfx2_ds, dfx_ds, fx1,fx2


    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    den = max(0.0_dp,den1in)
    spin1: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd1in)
       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       exunif = ax*(2.0_dp*den)**THIRD

       ! Mixing of PBE and RPBE with ratio x, with same mu and kappa

       ! PBE enhancement factor
       fx1 = 1.0_DP + kappa - kappa/(1.0_DP + mu*s2/kappa)
       ! RPBE enhancement factor
       fx2 = 1.0_DP + kappa*(1.0_DP-exp(-mu*s2/kappa))

       ! Derivatives
       dfx1_ds = 2.0_DP*mu*s/((1.0_DP + mu*s2/kappa)**2.0_DP)
       dfx2_ds = 2.0_DP*mu*s*exp(-mu*s2/kappa)

       fx = fx1*x + (1.0_DP - x)*fx2
       dfx_ds = dfx1_ds*x + (1.0_DP -x)*dfx2_ds

       x_energy =  x_energy + den * exunif * fx
       x_pot1 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd1 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy =  x_energy + den * exunif
       x_pot1 = FTHRD * exunif
       x_dfdmgd1 = 0.0_DP
    end if spin1

    den = max(0.0_dp,den2in)
    spin2: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd2in)
       kf = (THPISQ * 2.0_dp*den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       exunif = ax*(2.0_dp*den)**THIRD

       ! Mixing of PBE and RPBE with ratio x, with same mu and kappa

       ! PBE enhancement factor
       fx1 = 1.0_DP + kappa - kappa/(1.0_DP + mu*s2/kappa)
       ! RPBE enhancement factor
       fx2 = 1.0_DP + kappa*(1.0_DP-exp(-mu*s2/kappa))

       ! Derivatives
       dfx1_ds = 2.0_DP*mu*s/((1.0_DP + mu*s2/kappa)**2.0_DP)
       dfx2_ds = 2.0_DP*mu*s*exp(-mu*s2/kappa)

       fx = fx1*x + (1.0_DP - x)*fx2
       dfx_ds = dfx1_ds*x + (1.0_DP -x)*dfx2_ds

       x_energy =  x_energy + den * exunif * fx
       x_pot2 = FTHRD*ax*(fx*(2.0_dp*den)**THIRD - dfx_ds*mgd/(2.0_DP*den*THPISQ**THIRD))
       x_dfdmgd2 = ax * dfx_ds/(2.0_DP*(THPISQ)**THIRD)

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*(2.0_dp*den)**THIRD
       x_energy =  x_energy + den * exunif
       x_pot2 = FTHRD * exunif
       x_dfdmgd2 = 0.0_DP
    end if spin2

  end subroutine xc_vdwoptpbe_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_am05_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the AM05 exchange at a point in   !
    ! the spin unpolarised case. Is used in conjunction with a     !
    ! version of the VV10 non-local correlation                    !
    ! as part of the AM05-VV10sol functional                       !
    !  Armiento and Mattsson                                       !
    !  PHYSICAL REVIEW B 72, 085108 2005                           !
    !--------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in July 2014               !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! AM05 enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    ! Exchange constants
    real(kind=DP), parameter :: c_val = 0.7168_DP
    real(kind=DP), parameter ::  alpha_val = 2.804_DP

    ! Local variables
    real(kind=DP) :: chi_tilde
    real(kind=DP) :: dchitilde_ds
    real(kind=DP) :: chi_dtilde
    real(kind=DP) :: dchidtilde_ds
    real(kind=DP) :: lamb_val
    real(kind=DP) :: big_x
    real(kind=DP) :: s2
    real(kind=DP) :: dbigx_ds
    real(kind=DP) :: dfx_ds
    real(kind=DP) :: ds_drho

    if(den>dentol) then

       kf = (THPISQ * den)**(1.0_DP/3.0_DP)  ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       exunif = ax*(den**(1.0_DP/3.0_DP))

       ds_drho = -(4.0_DP/3.0_DP) * s/den

       if(mgd<1.0D-20) then
          x_energy =  den * exunif
          x_pot = 4.0_DP * exunif / 3.0_DP
          x_dfdmgd = 1.0D-10
       else

         big_x = 1.0_DP / (1.0_DP+alpha_val*s2)
         dbigx_ds = -2.0_DP* alpha_val*s / &
             ((1.0_DP+alpha_val*s2)**2.0_DP)

         lamb_val = lamb_w((s**(1.5_DP))/(2.0_DP*sqrt(6.0_DP)))

         chi_tilde = (1.5_DP*lamb_val)**(2.0_DP/3.0_DP)
         dchitilde_ds = lamb_val/(1.0_DP+lamb_val)*(3.0_DP/(2.0_DP*s))/sqrt(chi_tilde)

         chi_dtilde = ( ((4.0_DP/3.0_DP)**(4.0_DP/3.0_DP)) * &
             ((2.0_DP*PI/3.0_DP)**4.0_DP) * (chi_tilde**2.0_DP) + &
             (chi_tilde**4.0_DP) )**0.25_DP
         dchidtilde_ds = 0.5_DP * (chi_dtilde**(-3.0_DP))*dchitilde_ds * &
             (((4.0_DP/3.0_DP)**(4.0_DP/3.0_DP))*((2.0_DP*PI/3.0_DP)**4.0_DP)* &
             chi_tilde+2.0_DP*(chi_tilde**3.0_DP))

         fx = (c_val*s2 + 1.0_DP)/(1.0_DP + c_val*s*(3.0_DP/PI)*sqrt(chi_tilde)*chi_dtilde)

         x_energy = den*exunif*(big_x+(1.0_DP-big_x)*fx) ! energy

         dfx_ds =  (2.0_DP*c_val*s*(3.0_DP*c_val*s*chi_dtilde* &
            sqrt(chi_tilde)/PI+1.0_DP) - ((c_val*s2+1.0_DP)*c_val*3.0_DP/PI) &
            * (chi_dtilde*sqrt(chi_tilde)+dchidtilde_ds *sqrt(chi_tilde)*s+ &
            0.5_DP*s*chi_dtilde*dchitilde_ds/sqrt(chi_tilde)) ) / ((c_val*s* &
            3.0_DP*chi_dtilde*sqrt(chi_tilde)/PI+1.0_DP)**2.0_DP)

         ! potential:
         x_pot=(ax*4.0_DP/3.0_DP)*(den**(1.0_DP/3.0_DP))*(big_x+(1.0_DP-big_x) &
             * fx)+den*exunif*(dbigx_ds+dfx_ds-big_x*dfx_ds-dbigx_ds*fx)*ds_drho

         x_dfdmgd = den*exunif * (s/mgd) * (dbigx_ds + dfx_ds - big_x*dfx_ds - &
             dbigx_ds*fx )
       endif

    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 1.0D-10
    endif

  contains

    real(kind=DP) function lamb_w(z) ! computes Lambert W function

      implicit none

      real(kind=DP), intent(in) :: z

      ! local variables
      real(kind=DP) :: e_val,t_val,p_val
      integer :: iter

      !  If z too low, go with the first term of the power expansion, z
      if( z .lt. 1.0D-20) then
         lamb_w = z
      else

         e_val = exp(1.0_DP)
         !Inital guess
         if( abs(z + 1.0_DP/e_val) .gt. 1.45_DP ) then
            !Asymptotic expansion at 0 and Inf
            lamb_w = log(z)
            lamb_w = lamb_w - log(lamb_w)
         else
            ! Series expansion about -1/e to first order
            lamb_w = 1.0_DP*sqrt(2.0_DP*e_val*z + 2.0_DP) - 1.0_DP
         endif

         ! Find result through iteration
         do iter=1,10
            p_val = exp(lamb_w)
            t_val = lamb_w*p_val - z
            if(lamb_w .ne. -1.0_DP ) then
               t_val = t_val/(p_val*(lamb_w + 1.0_DP) - &
                   0.5_DP*(lamb_w + 2.0_DP)*t_val/(lamb_w + 1.0_DP))
            else
               t_val = 0.0_DP
            endif
            lamb_w = lamb_w - t_val
            if(abs(t_val) < (2.48_DP*1.0D-14)*(1.0_DP + abs(lamb_w))) exit
         enddo

      endif
    end function lamb_w


  end subroutine xc_am05_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_am05_x_point_sp(den1in,den2in,mgd1in,mgd2in,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the AM05 exchange at a point in   !
    ! the spin polarised case. Is used in conjunction with a       !
    ! version of the VV10 non-local correlation                    !
    ! as part of the AM05-VV10sol functional                       !
    !  Armiento and Mattsson                                       !
    !  PHYSICAL REVIEW B 72, 085108 2005                           !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons on 25/05/2021.                      !
    !==============================================================!

    use constants, only: DP, Pi
    implicit none

    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fx                   ! AM05 enhancement factor
    real(kind=DP) :: exunif               ! LDA exchange per particle

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    ! Exchange constants
    real(kind=DP), parameter :: c_val = 0.7168_DP
    real(kind=DP), parameter ::  alpha_val = 2.804_DP

    ! Local variables
    real(kind=DP) :: chi_tilde
    real(kind=DP) :: dchitilde_ds
    real(kind=DP) :: chi_dtilde
    real(kind=DP) :: dchidtilde_ds
    real(kind=DP) :: lamb_val
    real(kind=DP) :: big_x
    real(kind=DP) :: s2
    real(kind=DP) :: dbigx_ds
    real(kind=DP) :: dfx_ds
    real(kind=DP) :: ds_drho

    real(kind=dp) :: den, mgd

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    den = max(0.0_dp,den1in)
    spin1: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd1in)

       kf = (THPISQ * (2.0_dp*den))**(1.0_DP/3.0_DP)  ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       exunif = ax*((2.0_dp*den)**(1.0_DP/3.0_DP))

       ds_drho = -(4.0_DP/3.0_DP) * s/den

       if(mgd<1.0D-20) then
          x_energy =  (2.0_dp*den) * exunif
          x_pot1 = 4.0_DP * exunif / 3.0_DP
          x_dfdmgd1 = 1.0D-10
       else

          big_x = 1.0_DP / (1.0_DP+alpha_val*s2)
          dbigx_ds = -2.0_DP* alpha_val*s / &
               ((1.0_DP+alpha_val*s2)**2.0_DP)

          lamb_val = lamb_w((s**(1.5_DP))/(2.0_DP*sqrt(6.0_DP)))

          chi_tilde = (1.5_DP*lamb_val)**(2.0_DP/3.0_DP)
          dchitilde_ds = lamb_val/(1.0_DP+lamb_val)*(3.0_DP/(2.0_DP*s))/sqrt(chi_tilde)

          chi_dtilde = ( ((4.0_DP/3.0_DP)**(4.0_DP/3.0_DP)) * &
               ((2.0_DP*PI/3.0_DP)**4.0_DP) * (chi_tilde**2.0_DP) + &
               (chi_tilde**4.0_DP) )**0.25_DP
          dchidtilde_ds = 0.5_DP * (chi_dtilde**(-3.0_DP))*dchitilde_ds * &
               (((4.0_DP/3.0_DP)**(4.0_DP/3.0_DP))*((2.0_DP*PI/3.0_DP)**4.0_DP)* &
               chi_tilde+2.0_DP*(chi_tilde**3.0_DP))

          fx = (c_val*s2 + 1.0_DP)/(1.0_DP + c_val*s*(3.0_DP/PI)*sqrt(chi_tilde)*chi_dtilde)

          x_energy = x_energy + den*exunif*(big_x+(1.0_DP-big_x)*fx) ! energy

          dfx_ds =  (2.0_DP*c_val*s*(3.0_DP*c_val*s*chi_dtilde* &
               sqrt(chi_tilde)/PI+1.0_DP) - ((c_val*s2+1.0_DP)*c_val*3.0_DP/PI) &
               * (chi_dtilde*sqrt(chi_tilde)+dchidtilde_ds *sqrt(chi_tilde)*s+ &
               0.5_DP*s*chi_dtilde*dchitilde_ds/sqrt(chi_tilde)) ) / ((c_val*s* &
               3.0_DP*chi_dtilde*sqrt(chi_tilde)/PI+1.0_DP)**2.0_DP)

          ! potential:
          x_pot1=(ax*4.0_DP/3.0_DP)*((2.0_dp*den)**(1.0_DP/3.0_DP))*(big_x+(1.0_DP-big_x) &
               * fx)+den*exunif*(dbigx_ds+dfx_ds-big_x*dfx_ds-dbigx_ds*fx)*ds_drho

          x_dfdmgd1 = den*exunif * (s/mgd) * (dbigx_ds + dfx_ds - big_x*dfx_ds - &
               dbigx_ds*fx )
       endif

    else
       x_energy = 0.0_DP
       x_pot1 = 0.0_DP
       x_dfdmgd1 = 1.0D-10
    endif spin1

    den = max(0.0_dp,den2in)
    spin2: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd2in)

       kf = (THPISQ * (2.0_dp*den))**(1.0_DP/3.0_DP)  ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       s2 = s*s
       exunif = ax*((2.0_dp*den)**(1.0_DP/3.0_DP))

       ds_drho = -(4.0_DP/3.0_DP) * s/den

       if(mgd<1.0D-20) then
          x_energy =  (2.0_dp*den) * exunif
          x_pot2 = 4.0_DP * exunif / 3.0_DP
          x_dfdmgd2 = 1.0D-10
       else

          big_x = 1.0_DP / (1.0_DP+alpha_val*s2)
          dbigx_ds = -2.0_DP* alpha_val*s / &
               ((1.0_DP+alpha_val*s2)**2.0_DP)

          lamb_val = lamb_w((s**(1.5_DP))/(2.0_DP*sqrt(6.0_DP)))

          chi_tilde = (1.5_DP*lamb_val)**(2.0_DP/3.0_DP)
          dchitilde_ds = lamb_val/(1.0_DP+lamb_val)*(3.0_DP/(2.0_DP*s))/sqrt(chi_tilde)

          chi_dtilde = ( ((4.0_DP/3.0_DP)**(4.0_DP/3.0_DP)) * &
               ((2.0_DP*PI/3.0_DP)**4.0_DP) * (chi_tilde**2.0_DP) + &
               (chi_tilde**4.0_DP) )**0.25_DP
          dchidtilde_ds = 0.5_DP * (chi_dtilde**(-3.0_DP))*dchitilde_ds * &
               (((4.0_DP/3.0_DP)**(4.0_DP/3.0_DP))*((2.0_DP*PI/3.0_DP)**4.0_DP)* &
               chi_tilde+2.0_DP*(chi_tilde**3.0_DP))

          fx = (c_val*s2 + 1.0_DP)/(1.0_DP + c_val*s*(3.0_DP/PI)*sqrt(chi_tilde)*chi_dtilde)

          x_energy = x_energy + den*exunif*(big_x+(1.0_DP-big_x)*fx) ! energy

          dfx_ds =  (2.0_DP*c_val*s*(3.0_DP*c_val*s*chi_dtilde* &
               sqrt(chi_tilde)/PI+1.0_DP) - ((c_val*s2+1.0_DP)*c_val*3.0_DP/PI) &
               * (chi_dtilde*sqrt(chi_tilde)+dchidtilde_ds *sqrt(chi_tilde)*s+ &
               0.5_DP*s*chi_dtilde*dchitilde_ds/sqrt(chi_tilde)) ) / ((c_val*s* &
               3.0_DP*chi_dtilde*sqrt(chi_tilde)/PI+1.0_DP)**2.0_DP)

          ! potential:
          x_pot2=(ax*4.0_DP/3.0_DP)*((2.0_dp*den)**(1.0_DP/3.0_DP))*(big_x+(1.0_DP-big_x) &
               * fx)+den*exunif*(dbigx_ds+dfx_ds-big_x*dfx_ds-dbigx_ds*fx)*ds_drho

          x_dfdmgd2 = den*exunif * (s/mgd) * (dbigx_ds + dfx_ds - big_x*dfx_ds - &
               dbigx_ds*fx )
       endif

    else
       x_energy = 0.0_DP
       x_pot2 = 0.0_DP
       x_dfdmgd2 = 1.0D-10
    endif spin2

  contains

    real(kind=DP) function lamb_w(z) ! computes Lambert W function

      implicit none

      real(kind=DP), intent(in) :: z

      ! local variables
      real(kind=DP) :: e_val,t_val,p_val
      integer :: iter

      !  If z too low, go with the first term of the power expansion, z
      if( z .lt. 1.0D-20) then
         lamb_w = z
      else

         e_val = exp(1.0_DP)
         !Inital guess
         if( abs(z + 1.0_DP/e_val) .gt. 1.45_DP ) then
            !Asymptotic expansion at 0 and Inf
            lamb_w = log(z)
            lamb_w = lamb_w - log(lamb_w)
         else
            ! Series expansion about -1/e to first order
            lamb_w = 1.0_DP*sqrt(2.0_DP*e_val*z + 2.0_DP) - 1.0_DP
         endif

         ! Find result through iteration
         do iter=1,10
            p_val = exp(lamb_w)
            t_val = lamb_w*p_val - z
            if(lamb_w .ne. -1.0_DP ) then
               t_val = t_val/(p_val*(lamb_w + 1.0_DP) - &
                    0.5_DP*(lamb_w + 2.0_DP)*t_val/(lamb_w + 1.0_DP))
            else
               t_val = 0.0_DP
            endif
            lamb_w = lamb_w - t_val
            if(abs(t_val) < (2.48_DP*1.0D-14)*(1.0_DP + abs(lamb_w))) exit
         enddo

      endif
    end function lamb_w


  end subroutine xc_am05_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_am05_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the AM05 correlation at a point in!
    ! the spin unpolarised case. Is used in conjunction with a     !
    ! version of the VV10 non-local correlation                    !
    ! as part of the AM05-VV10sol functional                       !
    !  Armiento and Mattsson                                       !
    !  PHYSICAL REVIEW B 72, 085108 2005                           !
    !--------------------------------------------------------------!
    ! Written by Gabriel Constantinescu in July 2014               !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! correlation energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution df_{c}/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative df_{c}/d|grad n|

    real(kind=DP) :: s                    ! Dimensionless density gradient
    !real(kind=DP) :: fc                   ! AM05 correlation enhancement factor

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP*PI)
    ! Exchange constants
    real(kind=DP), parameter :: gamma_val = 0.8098_DP
    real(kind=DP), parameter ::  alpha_val = 2.804_DP


    ! Local variables
    real(kind=DP) :: s2
    real(kind=DP) :: big_x
    real(kind=DP) :: kf

    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: ecunif,eurs


    if(den>dentol) then

       kf = (THPISQ * den)**(1.0_DP/3.0_DP)   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient
       s2 = s*s
       big_x = 1.0_DP / (1.0_DP+alpha_val*s2)

       rs = (TONFPI / den)**(1.0_DP/3.0_DP)  ! Wigner-Seitz radius
       sqrtrs=sqrt(rs)

       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,ecunif,eurs)

       c_energy = den * ecunif * (big_x + (1.0_DP-big_x)*gamma_val)
       c_pot = (ecunif-1.0_DP*rs*eurs/3.0_DP) * (big_x+(1.0_DP-big_x)* &
               gamma_val) + ecunif*(1.0_DP-gamma_val)*(8.0_DP/3.0_DP)* &
               alpha_val*s2 / ((1.0_DP+alpha_val*s2)*(1.0_DP+alpha_val*s2))

       c_dfdmgd = - den*ecunif * (1.0_DP-gamma_val) * 2.0_DP*alpha_val*s/ &
               ((2.0_DP*den*kf)*(1.0_DP+alpha_val*s2)*(1.0_DP+alpha_val*s2))

       if(abs(c_dfdmgd)<1.0D-10) then
          c_dfdmgd = 1.0D-10
       endif

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
       c_dfdmgd = 0.0_DP
    endif

  end subroutine xc_am05_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_am05_c_point_sp(denin,den1in,den2in,mgdin,mgd1in,mgd2in,&
       c_energy,c_pot1,c_pot2,c_dfdmgd, c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the AM05 correlation at a point in!
    ! the spin polarised case. Is used in conjunction with a       !
    ! version of the VV10 non-local correlation                    !
    ! as part of the AM05-VV10sol functional                       !
    !  Armiento and Mattsson                                       !
    !  PHYSICAL REVIEW B 72, 085108 2005                           !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons on 25/05/2021.                      !
    !==============================================================!

    use constants, only: DP, Pi
    implicit none
    real(kind=dp), intent(in)  :: denin    ! Density
    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgdin
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{x}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    real(kind=DP) :: s                    ! Dimensionless density gradient
    !real(kind=DP) :: fc                   ! AM05 correlation enhancement factor

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP*PI)
    ! Exchange constants
    real(kind=DP), parameter :: gamma_val = 0.8098_DP
    real(kind=DP), parameter ::  alpha_val = 2.804_DP

    real(kind=DP), parameter :: zetamin = -1.0_DP + 1.0e-10_DP
    real(kind=DP), parameter :: zetamax = 1.0_DP - 1.0e-10_DP
    real(kind=DP), parameter :: gaminv = 1.92366105093153631981_DP
    real(kind=DP), parameter :: fzzinv = 1.0_DP / (8.0_DP * gaminv / 9.0_DP)


    ! Local variables
    real(kind=DP) :: s2
    real(kind=DP) :: big_x

    real(kind=DP) :: zeta         ! Spin polarisation
    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: ecunif               ! LDA correlation energy per particle
    real(kind=DP) :: ecd,ecz,eu,eurs,ecrs
    real(kind=DP) :: eprs,ep,alfm,alfrsm
    real(kind=DP) :: eczet
    real(kind=DP) :: zeta3,zeta4          ! zeta**3 and zeta**4
    real(kind=DP) :: fzeta,dfdz           ! Function of zeta and derivative
    real(kind=dp) :: den
    real(kind=dp) :: fac

    den = denin
    if(den>dentol) then

       zeta = (den1in - den2in) / den   ! Polarisation
       zeta = max(zeta,zetamin)
       zeta = min(zeta,zetamax)

       kf = (THPISQ * den)**(1.0_DP/3.0_DP)   ! Fermi wave-vector
       s = mgdin / (2.0_DP*den*kf)    ! Dimensionless gradient
       s2 = s*s
       big_x = 1.0_DP / (1.0_DP+alpha_val*s2)

       rs = (TONFPI / den)**(1.0_DP/3.0_DP)  ! Wigner-Seitz radius
       sqrtrs=sqrt(rs)

       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,eu,eurs)
       call xc_pw91_eq10(0.01554535_DP,0.20548_DP,14.1189_DP, &
            6.1977_DP,3.3662_DP,0.62517_DP,rs,sqrtrs,ep,eprs)
       call xc_pw91_eq10(0.0168869_DP,0.11125_DP,10.357_DP, &
            3.6231_DP,0.88026_DP,0.49671_DP,rs,sqrtrs,alfm,alfrsm)

       fac = (big_x + (1.0_DP-big_x)*gamma_val)

       fzeta = gaminv * &
            ((1.0_DP+zeta)**FTHRD+(1.0_DP-zeta)**FTHRD-2.0_DP)
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta3*zeta

       ecunif = eu*(1.0_DP-fzeta*zeta4) + ep*fzeta*zeta4 - &
            alfm*fzeta*(1.0_DP-zeta4)*fzzinv

       ecrs = eurs*(1.0_DP-fzeta*zeta4) + eprs*fzeta*zeta4 - &
            alfrsm*fzeta*(1.0_DP-zeta4)*fzzinv
       dfdz = gaminv * FTHRD * &
            ((1.0_DP+zeta)**THIRD-(1.0_DP-zeta)**THIRD)
       eczet = (4.0_DP*zeta3*fzeta+zeta4*dfdz) * &
            (ep-eu+alfm*fzzinv)-dfdz*alfm*fzzinv
       ecd = ecunif - THIRD*rs*ecrs
       ecz = eczet * den

       c_energy = den * ecunif * fac


       c_pot1 = (ecunif-1.0_DP*rs*ecrs/3.0_DP) * (big_x+(1.0_DP-big_x)* &
               gamma_val) + ecunif*(1.0_DP-gamma_val)*(8.0_DP/3.0_DP)* &
               alpha_val*s2 / ((1.0_DP+alpha_val*s2)*(1.0_DP+alpha_val*s2)) +&
               (1.0_DP-zeta)*ecz/den

       c_pot2 = (ecunif-1.0_DP*rs*ecrs/3.0_DP) * (big_x+(1.0_DP-big_x)* &
               gamma_val) + ecunif*(1.0_DP-gamma_val)*(8.0_DP/3.0_DP)* &
               alpha_val*s2 / ((1.0_DP+alpha_val*s2)*(1.0_DP+alpha_val*s2)) -&
               (1.0_DP+zeta)*ecz/den

       c_dfdmgd = - den*ecunif * (1.0_DP-gamma_val) * 2.0_DP*alpha_val*s/ &
               ((2.0_DP*den*kf)*(1.0_DP+alpha_val*s2)*(1.0_DP+alpha_val*s2))

       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP

       if(abs(c_dfdmgd)<1.0D-10) then
          c_dfdmgd = 1.0D-10
       endif

    else
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
       c_dfdmgd = 0.0_DP
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
    endif

  end subroutine xc_am05_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lyp_c_point_sp(denin,den1in,den2in,mgdin,mgd1in,mgd2in,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the LYP correlation at a point    !
    ! in the spin polarised case.                                  !
    !     C. Lee, W. Yang, and R.G. Parr                           !
    !     Development of the Colle-Salvetti correlation-energy     !
    !     formula into a functional of the electron density        !
    !     Phys. Rev. B37 (1988) 785-789                            !
    !  Using reformulation (to avoid Laplacian) given by:          !
    !     B. Miehlich, A. Savin, H. Stoll and H. Preuss            !
    !     Results obtained with the correlation energy density     !
    !     functionals of becke and Lee, Yang and Parr              !
    !     Chem. Phys. Lett. 157 (1989) 200-206                     !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 04/03/2009                        !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: denin    ! Density at point
    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgdin    ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    ! qoh: Local variables
    real(kind=DP) :: den, rden   ! Total density and its reciprocal
    real(kind=DP) :: den1!, rden1 ! Spin 1 density and its reciprocal
    real(kind=DP) :: den2!, rden2 ! Spin 2 density and its reciprocal
    real(kind=DP) :: mgd, mgdsq  ! Mod density gradient and square
    real(kind=DP) :: mgd1, mgd1sq ! Mod spin 1 density gradient and square
    real(kind=DP) :: mgd2, mgd2sq ! Mod spin 2 density gradient and square
    real(kind=DP) :: crden ! Cube root of den
    real(kind=DP) :: rcrden ! Reciprocal cube root of den
    real(kind=DP) :: densq, rdensq ! den squared and its reciprocal
    real(kind=DP) :: crden1sq ! Cube root of den spin 1 squared
    real(kind=DP) :: den1sq ! den spin 1 squared
    real(kind=DP) :: crden2sq ! Cube root of den spin 2 squared
    real(kind=DP) :: den2sq ! den spin 2 squared
    real(kind=DP) :: redenom, redenomsq ! Reciprocal of E denominator and square
    real(kind=DP) :: expo ! Exponential term
    !real(kind=DP) :: rcrdensq ! Reciprocal of the cube root of den squared
    real(kind=DP) :: omega ! As in paper
    real(kind=DP) :: delta ! As in paper
    real(kind=DP) :: dndelta ! delta density derivative
    real(kind=DP) :: dnomega ! omega density derivative
    real(kind=DP) :: inbracket ! Inner bracket [] from paper
    real(kind=DP) :: outbracket ! Outer bracket {} from paper
    real(kind=DP) :: dn1inbracket ! Spin 1 density derivative of inner bracket
    real(kind=DP) :: dn2inbracket ! Spin 2 density derivative of inner bracket

    ! qoh: Parameters:
    real(kind=DP), parameter :: aa = 0.04918_DP ! As in paper
    real(kind=DP), parameter :: bb = 0.132_DP   ! As in paper
    real(kind=DP), parameter :: ab = aa * bb
    real(kind=DP), parameter :: cc = 0.25330_DP ! As in paper
    real(kind=DP), parameter :: dd = 0.3490_DP  ! As in paper
    real(kind=DP), parameter :: TTHRD = 2.0_DP / 3.0_DP
    real(kind=DP), parameter :: ETHRD = 8.0_DP / 3.0_DP
    real(kind=DP), parameter :: cfp = 36.462398978764767_DP ! C_F * 2^{11/3}
    real(kind=DP), parameter :: on9 = 1.0_DP/9.0_DP
    real(kind=DP), parameter :: tw9 = 2.0_DP/9.0_DP
    real(kind=DP), parameter :: el9 = 11.0_DP/9.0_DP
    real(kind=DP), parameter :: tt9 = 22.0_DP/9.0_DP
    real(kind=DP), parameter :: on18 = 1.0_DP/18.0_DP
    real(kind=DP), parameter :: se18 = 7.0_DP/18.0_DP
    real(kind=DP), parameter :: fs18 = 47.0_DP/18.0_DP

    den = max(0.0_dp,denin)
    den1 = max(0.0_dp,den1in)
    den2 = max(0.0_dp,den2in)
    tol: if(den1.gt.dentol .and. den2.gt.dentol) then

       ! qoh: Calculate useful powers of den[ 12]
       densq = den**2
       rdensq = 1.0_DP / densq
       rden = 1.0_DP / den
       crden = den**THIRD
       rcrden = 1/crden
       !rcrdensq = rcrden**2
       den1sq = den1**2
       crden1sq = den1**TTHRD
       den2sq = den2**2
       crden2sq = den2**TTHRD
       ! qoh: Modulus of density gradient and square initialised
       mgd = max(0.0_dp,mgdin)
       mgdsq = mgd**2
       mgd1 = max(0.0_dp,mgd1in)
       mgd1sq = mgd1**2
       mgd2 = max(0.0_dp,mgd2in)
       mgd2sq = mgd2**2

       redenom = 1/( 1.0_dp+dd*rcrden )
       redenomsq = redenom**2
       expo = exp(-cc*rcrden)

       omega = expo * crden * redenom * rdensq * rdensq
       delta = cc*rcrden + dd*rcrden * redenom

       inbracket = CFP*(den1sq*crden1sq + den2sq*crden2sq)  &
            + mgdsq*(fs18-se18*delta) - (2.5_DP - on18*delta)*(mgd1sq+mgd2sq) &
            - (on9*delta -el9)*(den1*rden*mgd1sq+den2*rden*mgd2sq)

       outbracket = den1*den2 * inbracket - TTHRD*densq*mgdsq&
            + (TTHRD*densq-den1sq)*mgd2sq +  (TTHRD*densq-den2sq)*mgd1sq

       dnomega = ( THIRD * rdensq * rdensq * redenom * expo * rden ) * &
            (cc - 11.0_DP*crden + dd * redenom )
       dndelta = - THIRD * rden * rcrden * ( cc + dd * redenomsq )

       dn1inbracket = CFP*ETHRD*den1*crden1sq &
            - dndelta*(se18*mgdsq - on18*(mgd1sq + mgd2sq)) &
            - dndelta*on9*rden*(den1*mgd1sq + den2*mgd2sq) &
            - (on9*delta-el9)*den2*rdensq*(mgd1sq - mgd2sq)

       dn2inbracket =  CFP*ETHRD*den2*crden2sq &
            - dndelta*(se18*mgdsq - on18*(mgd1sq + mgd2sq)) &
            - dndelta*on9*rden*(den1*mgd1sq + den2*mgd2sq) &
            - (on9*delta-el9)*den1*rdensq*(mgd2sq - mgd1sq)

       c_energy = - aa*4.0_DP*den1*den2*rden*redenom - ab * omega * outbracket

       ! qoh: Calculate potential
       c_pot1 = -aa * 4.0_DP*(redenom *(den2sq*rdensq)&
            + THIRD*redenomsq*dd*rcrden*den1*den2 * rdensq)&
            - ab*dnomega*outbracket -ab*omega*den2*inbracket &
            - ab * omega*den1*den2*dn1inbracket &
            + ab*omega*FTHRD*den*(mgdsq-mgd1sq-mgd2sq)&
            + ab*omega*2.0_DP*den1*mgd2sq

       c_pot2 = -aa * 4.0_DP*(redenom *(den1sq*rdensq)&
            + THIRD*redenomsq*dd*rcrden*den1*den2 * rdensq)&
            - ab*dnomega*outbracket -ab*omega*den1*inbracket &
            - ab*omega*den1*den2* dn2inbracket &
            + ab*omega*FTHRD*den*(mgdsq-mgd1sq-mgd2sq)&
            + ab*omega*2.0_DP*den2*mgd1sq

       c_dfdmgd = 0.0_DP

       c_dfdmgd12 = - ab * omega * ( den1 * den2 *(fs18 - se18 * delta) &
            - TTHRD*densq) * (mgdsq - mgd1sq - mgd2sq)

       c_dfdmgd1 = -ab * omega * mgd1 * (den1 * den2 *((tw9 - TTHRD*delta) &
            - (tw9*delta - tt9) * den1 * rden) - 2.0_DP * den2sq )

       c_dfdmgd2 = -ab * omega * mgd2 * ( den1 * den2 *((tw9 - TTHRD*delta) &
            - (tw9*delta - tt9) * den2 * rden) - 2.0_DP * den1sq )

     else
       ! qoh: Set all values to zero for very small densities
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
       c_dfdmgd = 0.0_DP
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    endif tol

  end subroutine xc_lyp_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3lyp_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==================================================================!
    ! This subroutine calculates the B3LYP correlation at a point      !
    ! in the spin unpolarised case.                                    !
    ! B3LYP = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X            !
    !                + (1 - a_c) VWN_C + a_c LYP_C                     !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    ! P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and             !
    !   M. J. Frisch, J. Phys. Chem. 98, 11623 (1994)                  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot

    real(kind=DP), parameter :: lda_fac = 0.19_DP ! 1- a_c
    real(kind=DP), parameter :: lyp_fac = 0.81_DP ! a_c

    call xc_vwn_c_point(den,lda_energy,lda_pot)
    call xc_lyp_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)
    c_energy = c_energy*lyp_fac + lda_fac * lda_energy
    c_pot = c_pot * lyp_fac + lda_fac * lda_pot
    c_dfdmgd = c_dfdmgd * lyp_fac

   end subroutine xc_b3lyp_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3lyp_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==================================================================!
    ! This subroutine calculates the B3LYP correlation at a point      !
    ! in the spin polarised case.                                      !
    ! B3LYP = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X            !
    !                + (1 - a_c) VWN_C + a_c LYP_C                     !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    ! P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and             !
    !   M. J. Frisch, J. Phys. Chem. 98, 11623 (1994)                  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den   ! Density at point
    real(kind=dp), intent(in)  :: den1   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd    ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd1   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2

    real(kind=DP), parameter :: lda_fac = 0.19_DP ! 1- a_c
    real(kind=DP), parameter :: lyp_fac = 0.81_DP ! a_c

    call xc_vwn_c_point_sp(den,den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_lyp_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
         c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
    c_energy = c_energy*lyp_fac + lda_fac * lda_energy
    c_pot1 = c_pot1 * lyp_fac + lda_fac * lda_pot1
    c_pot2 = c_pot2 * lyp_fac + lda_fac * lda_pot2
    c_dfdmgd = c_dfdmgd * lyp_fac
    c_dfdmgd1 = c_dfdmgd1 * lyp_fac
    c_dfdmgd2 = c_dfdmgd2 * lyp_fac
    c_dfdmgd12 = c_dfdmgd12 * lyp_fac

  end subroutine xc_b3lyp_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3pw91_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==================================================================!
    ! This subroutine calculates the B3PW91 correlation at a point     !
    ! in the spin unpolarised case.                                    !
    ! B3PW91 = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X           !
    !                + (1 - a_c) VWN_C + a_c PW91_C                    !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot

    real(kind=DP), parameter :: lda_fac  = 0.19_DP ! 1- a_c
    real(kind=DP), parameter :: pw91_fac = 0.81_DP ! a_c

    call xc_vwn_c_point(den,lda_energy,lda_pot)
    call xc_pw91_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)
    c_energy = c_energy*pw91_fac + lda_fac * lda_energy
    c_pot = c_pot * pw91_fac + lda_fac * lda_pot
    c_dfdmgd = c_dfdmgd * pw91_fac

   end subroutine xc_b3pw91_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3pw91_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==================================================================!
    ! This subroutine calculates the B3PW91 correlation at a point     !
    ! in the spin polarised case.                                      !
    ! B3PW91 = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X           !
    !                + (1 - a_c) VWN_C + a_c PW91_C                    !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den   ! Density at point
    real(kind=dp), intent(in)  :: den1   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd    ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd1   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad_n_2|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2

    real(kind=DP), parameter :: lda_fac  = 0.19_DP ! 1- a_c
    real(kind=DP), parameter :: pw91_fac = 0.81_DP ! a_c

    call xc_vwn_c_point_sp(den,den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_pw91_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
         c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
    c_energy = c_energy*pw91_fac + lda_fac * lda_energy
    c_pot1 = c_pot1 * pw91_fac + lda_fac * lda_pot1
    c_pot2 = c_pot2 * pw91_fac + lda_fac * lda_pot2
    c_dfdmgd = c_dfdmgd * pw91_fac
    c_dfdmgd1 = c_dfdmgd1 * pw91_fac
    c_dfdmgd2 = c_dfdmgd2 * pw91_fac
    c_dfdmgd12 = c_dfdmgd12 * pw91_fac

  end subroutine xc_b3pw91_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x3lyp_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the X3LYP correlation at a point  !
    ! in the spin unpolarised case.                                !
    ! X3LYP = a_0 HF_X + (1 - a_0 - a_x) LDA_X                     !
    !         + a_x ( 0.765 B88_X + 0.235 PW91_X) + a_c VWN_C      !
    !         + (1 - a_C) LYP_C                                    !
    ! a_0 = 0.218, a_x = 0.709, a_c = 0.129.                       !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot

    real(kind=DP), parameter :: lda_fac = 0.129_DP ! a_c
    real(kind=DP), parameter :: lyp_fac = 0.871_DP ! 1 - a_c

    call xc_vwn_c_point(den,lda_energy,lda_pot)
    call xc_lyp_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)
    c_energy = c_energy*lyp_fac + lda_fac * lda_energy
    c_pot = c_pot * lyp_fac + lda_fac * lda_pot
    c_dfdmgd = c_dfdmgd * lyp_fac

  end subroutine xc_x3lyp_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x3lyp_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the X3LYP correlation at a point  !
    ! in the spin polarised case.                                  !
    ! X3LYP = a_0 HF_X + (1 - a_0 - a_x) LDA_X                     !
    !         + a_x ( 0.765 B88_X + 0.235 PW91_X) + a_c VWN_C      !
    !         + (1 - a_C) LYP_C                                    !
    ! a_0 = 0.218, a_x = 0.709, a_c = 0.129.                       !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den   ! Density at point
    real(kind=dp), intent(in)  :: den1   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd    ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd1   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2

    real(kind=DP), parameter :: lda_fac = 0.129_DP ! a_c
    real(kind=DP), parameter :: lyp_fac = 0.871_DP ! 1 - a_c

    call xc_vwn_c_point_sp(den,den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_lyp_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
         c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
    c_energy = c_energy*lyp_fac + lda_fac * lda_energy
    c_pot1 = c_pot1 * lyp_fac + lda_fac * lda_pot1
    c_pot2 = c_pot2 * lyp_fac + lda_fac * lda_pot2
    c_dfdmgd = c_dfdmgd * lyp_fac
    c_dfdmgd1 = c_dfdmgd1 * lyp_fac
    c_dfdmgd2 = c_dfdmgd2 * lyp_fac
    c_dfdmgd12 = c_dfdmgd12 * lyp_fac

  end subroutine xc_x3lyp_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine is a dummy for adding no exchange.           !
    !==============================================================!

    use constants, only: DP
    use utils, only: utils_use_var

    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: x_energy ! The X energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution dfxc/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative  dfxc/d|grad n|

    call utils_use_var(den+mgd)
    x_energy = 0.0_DP
    x_pot = 0.0_DP
    x_dfdmgd = 0.0_DP

  end subroutine xc_none_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine is a dummy for adding no correlation.        !
    !==============================================================!

    use constants, only: DP
    use utils, only: utils_use_var

    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    call utils_use_var(den+mgd)
    c_energy = 0.0_DP
    c_pot = 0.0_DP
    c_dfdmgd = 0.0_DP

  end subroutine xc_none_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine is a dummy for adding no exchange.           !
    !==============================================================!

    use constants, only: DP
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The correlation energy
    real(kind=dp), intent(out) :: x_pot1   ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: x_pot2   ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{c}/d|grad n_2|

    x_energy = 0.0_DP*(den1+den2+mgd1+mgd2)
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

  end subroutine xc_none_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine is a dummy for adding no correlation.        !
    !==============================================================!

    use constants, only: DP
    use utils, only: utils_use_var

    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1   ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2   ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|


    call utils_use_var(den + den1 + den2 + mgd + mgd1 + mgd2)

    c_energy = 0.0_DP
    c_pot1 = 0.0_DP
    c_pot2 = 0.0_DP
    c_dfdmgd = 0.0_DP
    c_dfdmgd1 = 0.0_DP
    c_dfdmgd2 = 0.0_DP
    c_dfdmgd12 = 0.0_DP

  end subroutine xc_none_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_mgga_x_point(den,grad,mgd,tau,x_energy,x_pot,x_dfdmgd,&
       x_dfdgrad,x_dfdtau,use_dfdgrad)

    !==============================================================!
    ! This subroutine is a dummy for adding no exchange with a     !
    ! meta-GGA function interface.                                 !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert, utils_use_var
    implicit none

    ! Arguments
    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: grad(3)  ! The charge density gradient
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau      ! The KE density at a point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|
    real(kind=dp), intent(out) :: x_dfdgrad(3) ! Derivative df_{x}/d(grad n)
    real(kind=DP), intent(out) :: x_dfdtau ! Derivative df_{x}/dtau
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! This routine outputs dfdmgd, rather than dfdgrad
    use_dfdgrad = .false.

    x_energy = 0.0_DP
    x_pot = 0.0_DP
    x_dfdmgd = 0.0_DP
    x_dfdgrad = 0.0_DP
    x_dfdtau = 0.0_DP

    call utils_use_var(den)
    call utils_use_var(sum(grad))
    call utils_use_var(mgd)
    call utils_use_var(tau)

  end subroutine xc_none_mgga_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_mgga_c_point(den,grad,mgd,tau,c_energy,c_pot,c_dfdmgd,&
       c_dfdgrad,c_dfdtau,use_dfdgrad)

    !==============================================================!
    ! This subroutine is a dummy for adding no correlation with a  !
    ! meta-GGA function interface.                                 !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert, utils_use_var
    implicit none

    ! Arguments
    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: grad(3)  ! The charge density gradient
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau      ! The KE density at a point
    real(kind=dp), intent(out) :: c_energy ! Exchange energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution df_{c}/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdgrad(3) ! Derivative df_{c}/d(grad n)
    real(kind=DP), intent(out) :: c_dfdtau ! Derivative df_{c}/dtau
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! This routine outputs dfdmgd, rather than dfdgrad
    use_dfdgrad = .false.

    c_energy = 0.0_DP
    c_pot = 0.0_DP
    c_dfdmgd = 0.0_DP
    c_dfdgrad = 0.0_DP
    c_dfdtau = 0.0_DP

    call utils_use_var(den)
    call utils_use_var(sum(grad))
    call utils_use_var(mgd)
    call utils_use_var(tau)

  end subroutine xc_none_mgga_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_mgga_x_point_sp(den1,den2,grad1,grad2,mgd1,mgd2,tau1,tau2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2,&
       x_dfdgrad1,x_dfdgrad2,x_dfdtau1,x_dfdtau2,use_dfdgrad)

    !==============================================================!
    ! This subroutine is a dummy for adding no exchange with a     !
    ! meta-GGA function interface.                                 !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert, utils_use_var
    implicit none

    real(kind=dp), intent(in)  :: den1      ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2      ! Density at point for spin 2
    real(kind=dp), intent(in)  :: grad1(3)  ! The charge density gradient for spin 1
                                            ! Cartesian components at a point
    real(kind=dp), intent(in)  :: grad2(3)  ! The charge density gradient for spin 2
                                            ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd1      ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2      ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau1      ! The KE density at point for spin 1
    real(kind=dp), intent(in)  :: tau2      ! The KE density at point for spin 2
    real(kind=dp), intent(out) :: x_energy  ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1    ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2    ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2
    real(kind=dp), intent(out) :: x_dfdgrad1(3) ! Derivative df_{x}/d(grad n_1)
    real(kind=dp), intent(out) :: x_dfdgrad2(3) ! Derivative df_{x}/d(grad n_2)
    real(kind=DP), intent(out) :: x_dfdtau1 ! Derivative df_{x}/dtau spin 1
    real(kind=DP), intent(out) :: x_dfdtau2 ! Derivative df_{x}/dtau spin 2
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! This routine outputs dfdmgd, rather than dfdgrad
    use_dfdgrad = .false.

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP
    x_dfdgrad1 = 0.0_DP
    x_dfdgrad2 = 0.0_DP
    x_dfdtau1 = 0.0_DP
    x_dfdtau2 = 0.0_DP

    call utils_use_var(den1+den2)
    call utils_use_var(sum(grad1)+sum(grad2))
    call utils_use_var(mgd1+mgd2)
    call utils_use_var(tau1+tau2)

  end subroutine xc_none_mgga_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_mgga_c_point_sp(den,den1,den2,grad1,grad2,mgd,mgd1,mgd2,tau1,tau2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12,&
       c_dfdgrad1,c_dfdgrad2,c_dfdtau1,c_dfdtau2,use_dfdgrad)

    !==============================================================!
    ! This subroutine is a dummy for adding no correlation with a  !
    ! meta-GGA function interface.                                 !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert, utils_use_var
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: grad1(3) ! The charge density gradient for spin 1
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: grad2(3) ! The charge density gradient for spin 2
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau1     ! The KE density at point for spin 1
    real(kind=dp), intent(in)  :: tau2     ! The KE density at point for spin 2
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1   ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2   ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|
    real(kind=dp), intent(out) :: c_dfdgrad1(3) ! Derivative df_{c}/d(grad n_1)
    real(kind=dp), intent(out) :: c_dfdgrad2(3) ! Derivative df_{c}/d(grad n_2)
    real(kind=DP), intent(out) :: c_dfdtau1 ! Derivative df_{c}/dtau spin 1
    real(kind=DP), intent(out) :: c_dfdtau2 ! Derivative df_{c}/dtau spin 2
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! This routine outputs dfdmgd, rather than dfdgrad
    use_dfdgrad = .false.

    c_energy   = 0.0_DP
    c_pot1     = 0.0_DP
    c_pot2     = 0.0_DP
    c_dfdmgd   = 0.0_DP
    c_dfdmgd1  = 0.0_DP
    c_dfdmgd2  = 0.0_DP
    c_dfdmgd12 = 0.0_DP
    c_dfdgrad1 = 0.0_DP
    c_dfdgrad2 = 0.0_DP
    c_dfdtau1  = 0.0_DP
    c_dfdtau2  = 0.0_DP

    call utils_use_var(den+den1+den2)
    call utils_use_var(sum(grad1)+sum(grad2))
    call utils_use_var(mgd+mgd1+mgd2)
    call utils_use_var(tau1+tau2)

  end subroutine xc_none_mgga_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lta_x_point(den,grad,mgd,tau,x_energy,x_pot,x_dfdmgd,&
       x_dfdgrad,x_dfdtau,use_dfdgrad)

    !==============================================================!
    ! This subroutine calculates the LTA (local-tau approximation) !
    ! exchange at a point in the spin unpolarized case.            !
    !  M Ernzerhof, GE Scuseria, J. Chem. Phys. 111, 911 (1999)    !
    ! Note that LTA is an unusual XC functional in that the        !
    ! the dependence of the total energy on the charge density is  !
    ! entirely through the implicit dependence of the kinetic      !
    ! energy density on the charge density. As a consequence,      !
    ! the derivatives dfdn and dfdmgd are always zero, and the     !
    ! XC potential is expressed entirely by dfdtau.                !
    !                                                              !
    ! This code should be considered experimental and should not   !
    ! be used in production calculations. It has not been tested   !
    ! rigorously.                                                  !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI, ONE_THIRD, TWO_THIRDS
    use rundat, only: tautol => pub_xc_min_tau
    use utils, only: utils_abort, utils_assert, utils_use_var
    implicit none

    ! Arguments
    real(kind=DP), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: grad(3)  ! The charge density gradient
                                           ! Cartesian components at a point
    real(kind=DP), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=DP), intent(in)  :: tau      ! The KE density at a point
    real(kind=DP), intent(out) :: x_energy ! Exchange energy
    real(kind=DP), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=DP), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|
    real(kind=dp), intent(out) :: x_dfdgrad(3) ! Derivative df_{x}/d(grad n)
    real(kind=DP), intent(out) :: x_dfdtau ! Derivative df_{x}/dtau
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! Parameters
    real(kind=DP), parameter   :: CxCtau45 = &
         -0.75_DP * (3.0_DP / PI)**ONE_THIRD * &
         ( 10.0_DP / (3.0_DP * (3.0_DP * PI*PI)**TWO_THIRDS ) )**0.8_DP
        ! Cx * Ctau^(4/5), defined as in Eqs. 16 and 17 of Ernzerhof and Scuseria's paper

    ! This routine outputs dfdmgd, rather than dfdgrad
    use_dfdgrad = .false.
    x_dfdgrad = 0.0_DP

    x_energy = 0.0_DP
    x_pot    = 0.0_DP ! LTA does not depend on density
    x_dfdmgd = 0.0_DP ! LTA does not depend on gradient of density
    x_dfdtau = 0.0_DP

    if ( (den < dentol).or.(tau < tautol) ) return

    x_energy = CxCtau45 * tau**(0.8_DP)

    ! WARNING: This potential is formally correct, based on the form of the
    ! energy expression given in Eq. 17 of Ernzerhof and Scuseria's 1999 paper.
    ! However, the form of the potential is problematic, with the exchange
    ! potential in regions with tau = 0.0 being formally infinite.
    x_dfdtau = 0.8_DP * CxCtau45 * tau**(-0.2_DP)

    call utils_use_var(sum(grad))
    call utils_use_var(mgd)

  end subroutine xc_lta_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lta_x_point_sp(den1,den2,grad1,grad2,mgd1,mgd2,tau1,tau2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2,&
       x_dfdgrad1,x_dfdgrad2,x_dfdtau1,x_dfdtau2,use_dfdgrad)

    ! WARNING: Not implemented
    !==============================================================!
    ! This subroutine calculates the LTA (local-tau approximation) !
    ! exchange at a point in the spin polarized case.              !
    !  M Ernzerhof, GE Scuseria, J. Chem. Phys. 111, 911 (1999)    !
    ! Note that LTA is an unusual XC functional in that the        !
    ! the dependence of the total energy on the charge density is  !
    ! entirely through the implicit dependence of the kinetic      !
    ! energy density on the charge density. As a consequence,      !
    ! the derivatives dfdn and dfdmgd are always zero, and the     !
    ! XC potential is expressed entirely by dfdtau.                !
    !                                                              !
    ! TODO Implement spin-polarized case                           !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI
    !use rundat, only: tautol => pub_xc_min_tau
    use utils, only: utils_abort, utils_assert, utils_use_var
    implicit none

    ! Arguments
    real(kind=dp), intent(in)  :: den1      ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2      ! Density at point for spin 2
    real(kind=dp), intent(in)  :: grad1(3)  ! The charge density gradient for spin 1
                                            ! Cartesian components at a point
    real(kind=dp), intent(in)  :: grad2(3)  ! The charge density gradient for spin 2
                                            ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd1      ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2      ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau1      ! The KE density at point for spin 1
    real(kind=dp), intent(in)  :: tau2      ! The KE density at point for spin 2
    real(kind=dp), intent(out) :: x_energy  ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1    ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2    ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2
    real(kind=dp), intent(out) :: x_dfdgrad1(3) ! Derivative df_{x}/d(grad n_1)
    real(kind=dp), intent(out) :: x_dfdgrad2(3) ! Derivative df_{x}/d(grad n_2)
    real(kind=DP), intent(out) :: x_dfdtau1 ! Derivative df_{x}/dtau spin 1
    real(kind=DP), intent(out) :: x_dfdtau2 ! Derivative df_{x}/dtau spin 2
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! This routine outputs dfdmgd, rather than dfdgrad
    use_dfdgrad = .false.

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP
    x_dfdgrad1 = 0.0_DP
    x_dfdgrad2 = 0.0_DP
    x_dfdtau1 = 0.0_DP
    x_dfdtau2 = 0.0_DP

    call utils_use_var(den1+den2)
    call utils_use_var(sum(grad1)+sum(grad2))
    call utils_use_var(mgd1+mgd2)
    call utils_use_var(tau1+tau2)

    ! JCW: Not implemented: Abort.
    call utils_abort("Error in xc_lta_x_point_sp: Spin polarized LTA exchange &
         &not implemented.")

  end subroutine xc_lta_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pkzb_x_point(den,grad,mgd,tau,x_energy,x_pot,x_dfdmgd,&
       x_dfdgrad,x_dfdtau,use_dfdgrad)

    !==============================================================!
    ! This subroutine calculates the PKZB exchange at a point in   !
    ! the spin unpolarized case.                                   !
    !    JP Perdew, S Kurth, A Zupan & P Blaha,                    !
    !    Phys. Rev. Lett. 82, 2544 (1999)                          !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI, ONE_THIRD, TWO_THIRDS, ONE_SIXTH
    !use rundat, only: tautol => pub_xc_min_tau <-- unused
    use utils, only: utils_abort, utils_assert, utils_use_var
    implicit none

    ! Arguments
    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: grad(3)  ! The charge density gradient
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau      ! The KE density at a point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|
    real(kind=dp), intent(out) :: x_dfdgrad(3) ! Derivative df_{x}/d(grad n)
    real(kind=DP), intent(out) :: x_dfdtau ! Derivative df_{x}/dtau
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! Local variables
    ! (References to Eqs are from original PKZB paper, cited in description of
    ! subroutine above)
    real(kind=DP) :: p        ! Defined in Eq. 10
    real(kind=DP) :: q        ! Defined in Eq. 12
    real(kind=DP) :: q_term_1 ! First term of Eq. 12
    real(kind=DP) :: q_term_1_o_tau ! First term of Eq. 12, divided by tau
    real(kind=DP) :: x        ! Defined in Eq. 14
    real(kind=DP) :: o_o_opxkappa ! 1/(1+x/kappa)
    real(kind=DP) :: Fx       ! Defined in Eq. 13
    real(kind=DP) :: exunif   ! Exchange energy of uniform electron (per unit V)
    real(kind=DP) :: o_o_den  ! 1.0/den
    real(kind=DP) :: o_o_mgd  ! 1.0/mgd
    real(kind=DP) :: dexunifdn! d(exunif)/dn
    real(kind=DP) :: dFxdn    ! d(Fx)/dn (Fx from Eq. 13)
    real(kind=DP) :: dFx_coeff! Common coefficient in derivatives of Fx (Eq. 13)
    real(kind=DP) :: dpdn     ! dp/dn (p from Eq. 10)
    real(kind=DP) :: dqdn     ! dq/dn (\widetilde{q} from Eq. 12)
    real(kind=DP) :: dxdn     ! dx/dn (x from Eq. 14)
    real(kind=DP) :: dFxdmgd  ! d(Fx)/d|grad(n)| (Fx from Eq. 13)
    real(kind=DP) :: dpdmgd   ! dp/|grad(n)| (p from Eq. 10)
    real(kind=DP) :: dqdmgd   ! dq/|grad(n)| (\widetilde{q} from Eq. 12)
    real(kind=DP) :: dxdmgd   ! dx/|grad(n)| (x from Eq. 14)
    real(kind=DP) :: dFxdtau  ! d(Fx)/dtau (Fx from Eq. 13)
    real(kind=DP) :: dqdtau   ! dq/tau (\widetilde{q} from Eq. 12)
    real(kind=DP) :: dxdtau   ! dx/tau (x from Eq. 14)


    ! Parameters
    real(kind=DP), parameter :: FOUR_THIRDS = 4.0_DP * ONE_THIRD
    real(kind=DP), parameter :: FIVE_THIRDS = 5.0_DP * ONE_THIRD
    real(kind=DP), parameter :: EIGHT_THIRDS = 2.0_DP * FOUR_THIRDS
    real(kind=DP), parameter :: NINE_TWENTIETHS = 9.0_DP / 20.0_DP
    real(kind=DP), parameter :: ONE_TWELFTH = 0.5_DP * ONE_SIXTH
    ! kappa = 0.804, D = 0.113
    real(kind=DP), parameter :: kappa = 0.804_DP
    real(kind=DP), parameter :: o_o_kappa = 1.0_DP/kappa
    real(kind=DP), parameter :: D = 0.113_DP
    ! exunif_coeff = -3/(4*PI) * (3*PI^2)^{1/3}
    real(kind=DP), parameter :: exunif_coeff = -3.0_DP**FOUR_THIRDS/(4.0_DP*PI**ONE_THIRD)
    ! p_coeff      = 1/( 4*(3*PI^2)^{2/3} )
    real(kind=DP), parameter :: p_coeff = 1.0_DP / (4.0_DP * (3.0_DP * PI**2.0_DP)**TWO_THIRDS )
    ! q_tau_coeff  = 3/( 2*(3*PI^2)^{2/3} ) = 6 * p_coeff
    real(kind=DP), parameter :: q_tau_coeff = 6.0_DP * p_coeff
    ! x_p_coeff    = 10 / 81
    real(kind=DP), parameter :: x_p_coeff = 10.0_DP / 81.0_DP
    ! x_q2_coeff   = 146 / 2025
    real(kind=DP), parameter :: x_q2_coeff = 146.0_DP / 2025.0_DP
    ! x_qp_coeff   = -73/405
    real(kind=DP), parameter :: x_qp_coeff = -73.0_DP / 405.0_DP
    ! x_p2_coeff   = D + (1/kappa)*(10/81)^{2} = D + (1/kappa)*x_p_coeff^{2}
    real(kind=DP), parameter :: x_p2_coeff = D + o_o_kappa*x_p_coeff*x_p_coeff

    ! This routine outputs dfdmgd, rather than dfdgrad
    use_dfdgrad = .false.

    x_energy = 0.0_DP
    x_pot = 0.0_DP
    x_dfdmgd = 0.0_DP
    x_dfdgrad = 0.0_DP
    x_dfdtau = 0.0_DP

    ! To agree with Libxc, we need to set the energy and potential to
    ! zero for the case where EITHER den < dentol or tau < tautol.
    ! The equations for PKZB do not seem to imply this behaviour, however.
    !tol: if ((den > dentol).and.(tau > tautol)) then
    tol: if (den > dentol) then

       ! Uniform electron gas exchange per unit volume
       exunif = exunif_coeff*den**FOUR_THIRDS

       ! EXCHANGE ENERGY
       ! [ Exchange energy per unit volume ]
       ! Intermediate quantities
       o_o_den = 1.0_DP/den
       p = p_coeff*mgd*mgd*o_o_den**EIGHT_THIRDS
       q_term_1_o_tau = q_tau_coeff*o_o_den**FIVE_THIRDS
       q_term_1 = q_term_1_o_tau*tau
       q = q_term_1 - NINE_TWENTIETHS - ONE_TWELFTH*p
       x = x_p_coeff*p + x_q2_coeff*q*q + x_qp_coeff*q*p + x_p2_coeff*p*p
       o_o_opxkappa = 1.0_DP/(1.0_DP + x*o_o_kappa)

       ! PKZB exchange enhancement factor
       Fx = 1.0_DP + kappa - kappa*o_o_opxkappa

       ! Energy expression
       x_energy = exunif*Fx

       ! DERIVATIVE WRT DENSITY
       ! [ Derivative of exchange energy per unit volume wrt density ]
       ! Intermediate quantities
       dpdn = -EIGHT_THIRDS*p*o_o_den
       dqdn = -FIVE_THIRDS*q_term_1*o_o_den - ONE_TWELFTH*dpdn
       dxdn = x_p_coeff*dpdn + 2.0_DP*x_q2_coeff*q*dqdn + &
            x_qp_coeff*(dqdn*p + q*dpdn) + 2.0_DP*x_p2_coeff*p*dpdn

       ! Derivative of uniform electron gas exchange per unit volume
       dexunifdn = FOUR_THIRDS*exunif*o_o_den

       ! Derivative of PKZB exchange enhancement factor
       dFx_coeff = o_o_opxkappa*o_o_opxkappa
       dFxdn = dFx_coeff*dxdn

       ! Derivative expression
       x_pot = dexunifdn*Fx + exunif*dFxdn

       ! DERIVATIVE WRT MODULUS OF DENSITY GRADIENT
       ! [ Derivative of exchange energy per unit volume wrt modulus of density
       ! gradient ]
       ! Intermediate quantities
       if (mgd.gt.0.0_DP) then
          o_o_mgd = 1.0_DP/mgd
       else
          o_o_mgd = 0.0_DP
       end if
       dpdmgd = 2.0_DP*p*o_o_mgd
       dqdmgd = -ONE_TWELFTH*dpdmgd
       dxdmgd = x_p_coeff*dpdmgd + 2.0_DP*x_q2_coeff*q*dqdmgd + &
            x_qp_coeff*(dqdmgd*p + q*dpdmgd) + 2.0_DP*x_p2_coeff*p*dpdmgd

       ! Derivative of PKZB exchange enhancement factor
       dFxdmgd = dFx_coeff*dxdmgd

       ! Derivative expression
       x_dfdmgd = exunif*dFxdmgd

       ! DERIVATIVE WRT KINETIC ENERGY DENSITY (TAU)
       ! [ Derivative of exchange energy per unit volume wrt KE density (tau) ]
       ! Intermediate quantities
       dqdtau = q_term_1_o_tau
       dxdtau = 2.0_DP*x_q2_coeff*q*dqdtau + x_qp_coeff*p*dqdtau

       ! Derivative of PKZB exchange enhancement factor
       dFxdtau = dFx_coeff*dxdtau

       ! Derivative expression
       x_dfdtau = exunif*dFxdtau
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
       x_dfdtau = 0.0_DP
    end if tol

    call utils_use_var(sum(grad))

  end subroutine xc_pkzb_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pkzb_c_point(den,grad,mgd,tau,c_energy,c_pot,c_dfdmgd,&
        c_dfdgrad,c_dfdtau,use_dfdgrad)

    !==============================================================!
    ! This subroutine calculates the PKZB correlation at a point   !
    ! in the spin unpolarized case.                                !
    !    JP Perdew, S Kurth, A Zupan & P Blaha,                    !
    !    Phys. Rev. Lett. 82, 2544 (1999)                          !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert
    implicit none

    ! Arguments
    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: grad(3)  ! The charge density gradient
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau      ! The KE density at a point
    real(kind=dp), intent(out) :: c_energy ! Exchange energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution df_{c}/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdgrad(3) ! Derivative df_{c}/d(grad n)
    real(kind=DP), intent(out) :: c_dfdtau ! Derivative df_{c}/dtau
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! Local variables
    real(kind=DP) :: c_pot1, c_pot2, c_dfdmgd1, c_dfdmgd2, c_dfdmgd12
    real(kind=DP) :: c_dfdgrad1(3), c_dfdgrad2(3), c_dfdtau1, c_dfdtau2
    logical       :: c_use_dfdgrad
    real(kind=DP) :: half_grad(3)

    ! This routine outputs dfdgrad, rather than dfdmgd
    use_dfdgrad = .true.

    ! Do not use dfdmgd
    c_dfdmgd = 0.0_DP

    !subroutine xc_pkzb_c_point_sp(den,den1,den2,grad1,grad2,mgd,mgd1,mgd2,tau1,tau2,&
    !     c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12,&
    !     c_dfdgrad1,c_dfdgrad2,c_dfdtau1,c_dfdtau2,use_dfdgrad)

    ! jme: Prevent creation of array temporary
    half_grad(:) = 0.5_DP * grad(:)

    ! Naive implementation using spin-polarized version of function
    ! Temporary: May be more efficient to write a specific spin-unpolarized
    !            routine.
    call xc_pkzb_c_point_sp(den,0.5_DP*den,0.5_DP*den,half_grad,half_grad,&
                            mgd,0.5_DP*mgd,0.5_DP*mgd,&
                            0.5_DP*tau,0.5_DP*tau,c_energy,c_pot1,c_pot2,&
                            c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12,&
                            c_dfdgrad1,c_dfdgrad2,c_dfdtau1,c_dfdtau2,&
                            c_use_dfdgrad)

    call utils_assert(c_use_dfdgrad,"Error in xc_pkzb_c_point: spin-polarized &
         &routine should return c_use_dfdgrad == .true.")

    c_pot     = 0.5_DP * c_pot1 + 0.5_DP * c_pot2
    c_dfdgrad(:) = 0.5_DP * c_dfdgrad1(:) + 0.5_DP * c_dfdgrad2(:)
    c_dfdtau  = 0.5_DP * c_dfdtau1 + 0.5_DP * c_dfdtau2

  end subroutine xc_pkzb_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pkzb_x_point_sp(den1,den2,grad1,grad2,mgd1,mgd2,tau1,tau2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2,&
       x_dfdgrad1,x_dfdgrad2,x_dfdtau1,x_dfdtau2,use_dfdgrad)

    ! WARNING: Not implemented
    !==============================================================!
    ! This subroutine calculates the PKZB exchange at a point in   !
    ! the spin polarized case.                                     !
    !    JP Perdew, S Kurth, A Zupan & P Blaha,                    !
    !    Phys. Rev. Lett. 82, 2544 (1999)                          !
    !                                                              !
    ! TODO Implement spin-polarized case                           !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert, utils_use_var
    implicit none

    ! Arguments
    real(kind=dp), intent(in)  :: den1      ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2      ! Density at point for spin 2
    real(kind=dp), intent(in)  :: grad1(3)  ! The charge density gradient for spin 1
                                            ! Cartesian components at a point
    real(kind=dp), intent(in)  :: grad2(3)  ! The charge density gradient for spin 2
                                            ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd1      ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2      ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau1      ! The KE density at point for spin 1
    real(kind=dp), intent(in)  :: tau2      ! The KE density at point for spin 2
    real(kind=dp), intent(out) :: x_energy  ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1    ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2    ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2
    real(kind=dp), intent(out) :: x_dfdgrad1(3) ! Derivative df_{x}/d(grad n_1)
    real(kind=dp), intent(out) :: x_dfdgrad2(3) ! Derivative df_{x}/d(grad n_2)
    real(kind=DP), intent(out) :: x_dfdtau1 ! Derivative df_{x}/dtau spin 1
    real(kind=DP), intent(out) :: x_dfdtau2 ! Derivative df_{x}/dtau spin 2
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! This routine outputs dfdgrad, rather than dfdmgd
    use_dfdgrad = .true.

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP
    x_dfdgrad1 = 0.0_DP
    x_dfdgrad2 = 0.0_DP
    x_dfdtau1 = 0.0_DP
    x_dfdtau2 = 0.0_DP

    call utils_use_var(den1+den2)
    call utils_use_var(sum(grad1)+sum(grad2))
    call utils_use_var(mgd1+mgd2)
    call utils_use_var(tau1+tau2)

    ! JCW: Not implemented: Abort.
    call utils_abort("Error in xc_pkzb_x_point_sp: Spin polarized PKZB exchange &
         &not implemented.")

  end subroutine xc_pkzb_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pkzb_c_point_sp(den,den1,den2,grad1,grad2,mgd,mgd1,mgd2,tau1,tau2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12,&
       c_dfdgrad1,c_dfdgrad2,c_dfdtau1,c_dfdtau2,use_dfdgrad)

    ! WARNING: Only tested when called from xc_pkzb_c_point (spin
    !          unpolarized case)
    !==============================================================!
    ! This subroutine calculates the PKZB correlation at a point   !
    ! in the spin polarized case.                                  !
    !    JP Perdew, S Kurth, A Zupan & P Blaha,                    !
    !    Phys. Rev. Lett. 82, 2544 (1999)                          !
    !                                                              !
    ! TODO Test spin-polarized case                                !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI, ONE_THIRD
    use rundat, only: tautol => pub_xc_min_tau
    use utils, only: utils_abort, utils_assert, utils_sanity_check
    implicit none

    ! Arguments
    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: grad1(3) ! The charge density gradient for spin 1
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: grad2(3) ! The charge density gradient for spin 2
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau1     ! The KE density at point for spin 1
    real(kind=dp), intent(in)  :: tau2     ! The KE density at point for spin 2
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1   ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2   ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|
    real(kind=dp), intent(out) :: c_dfdgrad1(3) ! Derivative df_{c}/d(grad n_1)
    real(kind=dp), intent(out) :: c_dfdgrad2(3) ! Derivative df_{c}/d(grad n_2)
    real(kind=DP), intent(out) :: c_dfdtau1 ! Derivative df_{c}/dtau spin 1
    real(kind=DP), intent(out) :: c_dfdtau2 ! Derivative df_{c}/dtau spin 2
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! Parameters
    real(kind=DP), parameter :: C = 0.53_DP

    ! Local variables
    real(kind=DP) :: tw(2)      ! Weizacker KE density spin 1 and 2
    real(kind=DP) :: dtwdgrad(3,2) ! dtw(is)/(grad n(is)) for spin 1 and spin 2
    real(kind=DP) :: tw_o_t12(2)    ! [ (tau_w) / (tau1 + tau2) ], tau_w spin 1 and 2
    real(kind=DP) :: tw12_o_t12     ! [ (tau_w1 + tau_w2) / (tau1 + tau2) ]
    real(kind=DP) :: tw12_o_t12_2   ! [ (tau_w1 + tau_w2) / (tau1 + tau2) ]**2
    real(kind=DP) :: tw12_o_t12_t12 ! (tau_w1 + tau_w2) / [ (tau1 + tau2) ]**2
    real(kind=DP) :: tw12_2_o_t12_3 ! [ (tau_w1 + tau_w2) ]**2 / [ (tau1 + tau2) ]**3
    real(kind=DP) :: tw_o_t(2)   ! [ tau_w / tau ] for tau_w, tau spin 1 and 2
    real(kind=DP) :: tw_o_t_2(2) ! [ tau_w / tau ]**2 for tau_w, tau spin 1 and 2
    real(kind=DP) :: tw_o_t_t(2) ! tau_w / [ tau ]**2 for tau_w, tau spin 1 and 2
    real(kind=DP) :: tw_2_o_t_3(2) ! [ tau_w ]**2 / [ tau ]**3 for tau_w, tau spin 1 and 2
    real(kind=DP) :: c_energy_pbe_os
    real(kind=DP) :: c_pot1_pbe_os
    real(kind=DP) :: c_pot2_pbe_os
    real(kind=DP) :: c_dfdmgd_pbe_os
    real(kind=DP) :: c_dfdmgd1_pbe_os
    real(kind=DP) :: c_dfdmgd2_pbe_os
    real(kind=DP) :: c_dfdmgd12_pbe_os
    real(kind=DP) :: c_dfdgrad_pbe_os(3)
    real(kind=DP) :: c_energy_pbe_ss(2)
    real(kind=DP) :: c_pot1_pbe_ss(2)
    real(kind=DP) :: c_pot2_pbe_ss(2)
    real(kind=DP) :: c_dfdmgd_pbe_ss(2)
    real(kind=DP) :: c_dfdmgd1_pbe_ss(2)
    real(kind=DP) :: c_dfdmgd2_pbe_ss(2)
    real(kind=DP) :: c_dfdmgd12_pbe_ss(2)
    real(kind=DP) :: c_dfdgrad_pbe_ss(3,2)
    real(kind=DP) :: f1
    real(kind=DP) :: df1dden(2) ! df1/dn
    real(kind=DP) :: df1dgrad(3,2) ! df1/d(grad n(is))
    real(kind=DP) :: df1dtau ! df1/dtau_1 = df1/dtau_2
    real(kind=DP) :: f2(2)
    real(kind=DP) :: df2dden(2) ! df2(is)/dn(is)
    real(kind=DP) :: df2dgrad(3,2) ! df2/d(grad n(is))
    real(kind=DP) :: df2dtau(2) ! df2(is)/dtau(is)

    integer :: is

    real(kind=DP) :: o_o_tau1, o_o_tau2, o_o_tau12
    real(kind=DP) :: o_o_den1, o_o_den2
    real(kind=DP) :: o_o_mgd1, o_o_mgd2, o_o_mgd
    ! JCW: For division by mgd, simply avoid divide-by-zero errors
    real(kind=DP) :: mgdtol = 0.0_DP

    ! This routine outputs dfdgrad, rather than dfdmgd
    use_dfdgrad = .true.

    ! Do not use c_dfdmgd
    c_dfdmgd   = 0.0_DP
    c_dfdmgd1  = 0.0_DP
    c_dfdmgd2  = 0.0_DP
    c_dfdmgd12 = 0.0_DP

    ! Initialize energy to zero
    c_energy   = 0.0_DP

    if ( (den1+den2).gt.dentol ) then

       if ( den1.gt.dentol ) then
          o_o_den1 = 1.0_DP / den1
       else
          o_o_den1 = 0.0_DP
       end if
       if ( den2.gt.dentol ) then
          o_o_den2 = 1.0_DP / den2
       else
          o_o_den2 = 0.0_DP
       end if

       if (mgd1.gt.mgdtol ) then
          o_o_mgd1 = 1.0_DP / mgd1
       else
          o_o_mgd1 = 0.0_DP
       end if
       if (mgd2.gt.mgdtol ) then
          o_o_mgd2 = 1.0_DP / mgd2
       else
          o_o_mgd2 = 0.0_DP
       end if
       if (mgd.gt.mgdtol) then
          o_o_mgd = 1.0_DP / mgd
       else
          o_o_mgd = 0.0_DP
       end if

       if (tau1.gt.tautol ) then
          o_o_tau1 = 1.0_DP / tau1
       else
          o_o_tau1 = 0.0_DP
       end if
       if (tau2.gt.tautol ) then
          o_o_tau2 = 1.0_DP / tau2
       else
          o_o_tau2 = 0.0_DP
       end if
       if ((tau1+tau2).gt.tautol) then
          o_o_tau12 = 1.0_DP / (tau1+tau2)
       else
          o_o_tau12 = 0.0_DP
       end if

       tw(1)         = 0.125_DP * mgd1*mgd1 * o_o_den1
       dtwdgrad(:,1) = 0.25_DP  * grad1(:) * o_o_den1 ! dtw/d(grad n_1)
       tw(2)         = 0.125_DP * mgd2*mgd2 * o_o_den2
       dtwdgrad(:,2) = 0.25_DP  * grad2(:) * o_o_den2 ! dtw/d(grad n_2)

       tw_o_t(1)   = tw(1) * o_o_tau1
       tw_o_t_t(1) = tw_o_t(1) * o_o_tau1
       tw_o_t_2(1) = tw_o_t(1) * tw_o_t(1)
       tw_2_o_t_3(1) = tw_o_t_2(1) * o_o_tau1

       tw_o_t(2)   = tw(2) * o_o_tau2
       tw_o_t_t(2) = tw_o_t(2) * o_o_tau2
       tw_o_t_2(2) = tw_o_t(2) * tw_o_t(2)
       tw_2_o_t_3(2) = tw_o_t_2(2) *o_o_tau2

       tw_o_t12(1) = tw(1) * o_o_tau12
       tw_o_t12(2) = tw(2) * o_o_tau12
       tw12_o_t12 = tw_o_t12(1) + tw_o_t12(2)
       tw12_o_t12_t12 = tw12_o_t12 * o_o_tau12
       tw12_o_t12_2   = tw12_o_t12 * tw12_o_t12
       tw12_2_o_t12_3 = tw12_o_t12_2 * o_o_tau12

       ! Opposite spin
       call xc_pbe_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
            c_energy_pbe_os,c_pot1_pbe_os,c_pot2_pbe_os,&
            c_dfdmgd_pbe_os,c_dfdmgd1_pbe_os,c_dfdmgd2_pbe_os,c_dfdmgd12_pbe_os)

       ! Same spin
       call xc_pbe_c_point_sp(den1,den1,0.0_DP,mgd1,mgd1,0.0_DP,&
            c_energy_pbe_ss(1),c_pot1_pbe_ss(1),c_pot2_pbe_ss(1),&
            c_dfdmgd_pbe_ss(1),c_dfdmgd1_pbe_ss(1),&
            c_dfdmgd2_pbe_ss(1),c_dfdmgd12_pbe_ss(1))
       call xc_pbe_c_point_sp(den2,den2,0.0_DP,mgd2,mgd2,0.0_DP,&
            c_energy_pbe_ss(2),c_pot1_pbe_ss(2),c_pot2_pbe_ss(2),&
            c_dfdmgd_pbe_ss(2),c_dfdmgd1_pbe_ss(2),&
            c_dfdmgd2_pbe_ss(2),c_dfdmgd12_pbe_ss(2))

       ! Evaluate spin-polarized energy expression
       ! Opposite spin part
       f1 = 1.0_DP + C * tw12_o_t12_2
       c_energy = c_energy_pbe_os * f1
       ! Same spin part
       do is = 1, 2
         f2(is) = (1.0_DP + C) * tw_o_t_2(is)
         c_energy = c_energy - f2(is) * c_energy_pbe_ss(is)
       end do

       df1dden(1) = -( 2.0_DP * C * tw12_o_t12 * tw_o_t12(1) ) * o_o_den1
       df2dden(1) = -( 2.0_DP * (1.0_DP + C) * tw_o_t_2(1) ) * o_o_den1
       df1dden(2) = -( 2.0_DP * C * tw12_o_t12 * tw_o_t12(2) ) * o_o_den2
       df2dden(2) = -( 2.0_DP * (1.0_DP + C) * tw_o_t_2(2) ) * o_o_den2

       df1dgrad(:,1) = 2.0_DP * C * tw12_o_t12_t12 * dtwdgrad(:,1)
       df1dgrad(:,2) = 2.0_DP * C * tw12_o_t12_t12 * dtwdgrad(:,2)
       df2dgrad(:,1) = 2.0_DP * (1.0_DP + C) * tw_o_t_t(1) * dtwdgrad(:,1)
       df2dgrad(:,2) = 2.0_DP * (1.0_DP + C) * tw_o_t_t(2) * dtwdgrad(:,2)

       ! NB. df1dtau is the same for both tau_1, tau_2, since the expression for
       ! f1 is symmetrical wrt both kinetic energy densities.
       df1dtau = -2.0_DP * C * tw12_2_o_t12_3

       df2dtau(1) = -2.0_DP * (1.0_DP + C) * tw_2_o_t_3(1)
       df2dtau(2) = -2.0_DP * (1.0_DP + C) * tw_2_o_t_3(2)

       !! NB. c_dfdgrad1 and c_dfdgrad2 are identical for PBE, since PBE correlation
       !! depends on grad(n_1) and grad(n_2) via |grad(n)| and the is therefore
       !! symmetrical wrt both densities.
       c_dfdgrad_pbe_os(:) = c_dfdmgd_pbe_os * (grad1(:) + grad2(:)) * o_o_mgd
       c_dfdgrad_pbe_ss(:,1) = c_dfdmgd_pbe_ss(1) * grad1(:) * o_o_mgd1
       c_dfdgrad_pbe_ss(:,2) = c_dfdmgd_pbe_ss(2) * grad2(:) * o_o_mgd2

       ! Evaluate gradient of energy per unit volume wrt density
       c_pot1 = c_pot1_pbe_os * f1 + df1dden(1) * c_energy_pbe_os &
            - ( c_pot1_pbe_ss(1) * f2(1) + df2dden(1) * c_energy_pbe_ss(1) )
       c_pot2 = c_pot2_pbe_os * f1 + df1dden(2) * c_energy_pbe_os &
            - ( c_pot1_pbe_ss(2) * f2(2) + df2dden(2) * c_energy_pbe_ss(2) )

       ! Evaluate gradient of energy per unit volume wrt density gradient
       c_dfdgrad1(:) = c_dfdgrad_pbe_os(:) * f1 + c_energy_pbe_os * df1dgrad(:,1) &
            - ( c_dfdgrad_pbe_ss(:,1) * f2(1) + c_energy_pbe_ss(1) * df2dgrad(:,1) )
       c_dfdgrad2(:) = c_dfdgrad_pbe_os(:) * f1 + c_energy_pbe_os * df1dgrad(:,2) &
            - ( c_dfdgrad_pbe_ss(:,2) * f2(2) + c_energy_pbe_ss(2) * df2dgrad(:,2) )

       ! Evaluate gradient of energy per unit volume wrt kinetic energy density
       c_dfdtau1 = c_energy_pbe_os * df1dtau - c_energy_pbe_ss(1) * df2dtau(1)
       c_dfdtau2 = c_energy_pbe_os * df1dtau - c_energy_pbe_ss(2) * df2dtau(2)

    else
       c_energy = 0.0_DP
       c_pot1   = 0.0_DP
       c_pot2   = 0.0_DP
       c_dfdgrad1 = 0.0_DP
       c_dfdgrad2 = 0.0_DP
       c_dfdtau1 = 0.0_DP
       c_dfdtau2 = 0.0_DP
    end if

  end subroutine xc_pkzb_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b97mv_x_point(den,grad,mgd,tau,x_energy,x_pot,x_dfdmgd,&
       x_dfdgrad,x_dfdtau,use_dfdgrad)

    !==============================================================!
    ! This subroutine calculates the the local part of B97M-V      !
    ! exchange in the spin unpolarized case. It forms a wrapper    !
    ! around the corresponding routine in the xc_b97mv module,     !
    ! which is a modified version of code originally provided by   !
    ! N. Mardirossian.                                             !
    !                                                              !
    ! Note that global dentol and tautol values do not apply to    !
    ! B97M-V, since these are set in xc_b97mv on a per-routine     !
    ! basis.                                                       !
    !                                                              !
    ! B97M-V citation:                                             !
    !       N. Mardirossian & M. Head-Gordon,                      !
    !            J. Chem. Phys. 142, 074111 (2015)                 !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert, utils_use_var
    use xc_b97mv, only: N_D1_META, b97m_v_x
    implicit none

    ! Arguments
    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: grad(3)  ! The charge density gradient
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau      ! The KE density at a point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|
    real(kind=dp), intent(out) :: x_dfdgrad(3) ! Derivative df_{x}/d(grad n)
    real(kind=DP), intent(out) :: x_dfdtau ! Derivative df_{x}/dtau
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! Local variables
    real(kind=DP) :: D1F(N_D1_META)
    real(kind=DP) :: half_grad(3)

    !subroutine b97m_v_x(F,D1F,RhoA,RhoB,D1RhoA,D1RhoB,TauA,TauB)
    !  real(kind=DP),intent(out) :: F
    !  real(kind=DP),intent(out) :: D1F(N_D1_META)
    !  real(kind=DP),intent(in) :: RhoA,RhoB
    !  real(kind=DP),intent(in) :: D1RhoA(3),D1RhoB(3)
    !  real(kind=DP),intent(in) :: TauA,TauB
    !end subroutine b97m_v_x

    ! This routine outputs dfdgrad, rather than dfdmgd
    use_dfdgrad = .true.

    x_energy = 0.0_DP
    x_pot = 0.0_DP
    x_dfdmgd = 0.0_DP
    x_dfdgrad = 0.0_DP
    x_dfdtau = 0.0_DP

    ! Avoid temporary array copy in subroutine call
    half_grad(:) = 0.5_DP * grad(:)



    ! Notes:
    ! * In b97m_v_x tau is defined as 2*tau from ONETEP (i.e. the definition
    !   does not include the usual factor of a 1/2).
    call b97m_v_x(x_energy,D1F,&
         0.5_DP*den,0.5_DP*den,&
         half_grad, half_grad,&
         tau,tau)

    ! Unpack D1F array
    !integer,parameter :: POS_RA  = 1 ! df/dnA
    !integer,parameter :: POS_RB  = 2 ! df/dnB
    !integer,parameter :: POS_GAA = 3 ! df/d(grad(nA).grad(nA))
    !integer,parameter :: POS_GBB = 4 ! df/d(grad(nB).grad(nB))
    !integer,parameter :: POS_GAB = 5 ! df/d(grad(nA).grad(nB)) [ unused ]
    !integer,parameter :: POS_TA  = 6 ! df/dtauA
    !integer,parameter :: POS_TB  = 7 ! df/dtauB
    x_pot     = 0.5_DP * D1F(1) + 0.5_DP * D1F(2)
    ! Use dfdgrad, rather than dfdmgd. We need dfdgrad in the calling routine
    ! (xc_mgga), so getting this directly from df/d(grad(nA).grad(nA)) avoids
    ! the transformation to dfdmgd altogether.
    x_dfdgrad(:) = (0.5_DP * D1F(3) + 0.5_DP * D1F(4))*grad(:)
    x_dfdtau  = D1F(6) + D1F(7)

    call utils_use_var(mgd)

  end subroutine xc_b97mv_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b97mv_c_point(den,grad,mgd,tau,c_energy,c_pot,c_dfdmgd,&
       c_dfdgrad,c_dfdtau,use_dfdgrad)

    !==============================================================!
    ! This subroutine calculates the the local part of B97M-V      !
    ! correlation in the spin unpolarized case. It forms a wrapper !
    ! around the corresponding routine in the xc_b97mv module,     !
    ! which is a modified version of code originally provided by   !
    ! N. Mardirossian.                                             !
    !                                                              !
    ! Note that global dentol and tautol values do not apply to    !
    ! B97M-V, since these are set in xc_b97mv on a per-routine     !
    ! basis.                                                       !
    !                                                              !
    ! B97M-V citation:                                             !
    !       N. Mardirossian & M. Head-Gordon,                      !
    !            J. Chem. Phys. 142, 074111 (2015)                 !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by James C. Womack, 2016                             !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert, utils_use_var
    use xc_b97mv, only: N_D1_META, b97m_v_css, b97m_v_cos
    implicit none

    ! Arguments
    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: grad(3)  ! The charge density gradient
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau      ! The KE density at a point
    real(kind=dp), intent(out) :: c_energy ! Exchange energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution df_{c}/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdgrad(3) ! Derivative df_{c}/d(grad n)
    real(kind=DP), intent(out) :: c_dfdtau ! Derivative df_{c}/dtau
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! Local variables
    real(kind=DP) :: Fss, Fos
    real(kind=DP) :: D1Fss(N_D1_META), D1Fos(N_D1_META)
    real(kind=DP) :: half_grad(3)

    ! This routine outputs dfdmgd, rather than dfdgrad
    use_dfdgrad = .true.

    c_energy = 0.0_DP
    c_pot = 0.0_DP
    c_dfdmgd = 0.0_DP
    c_dfdgrad = 0.0_DP
    c_dfdtau = 0.0_DP

    ! Avoid temporary array copy in subroutine call
    half_grad(:) = 0.5_DP * grad(:)

    ! Notes:
    ! * In b97m_v_css, b97m_v_cos, tau is defined as 2*tau from ONETEP
    !   (i.e. the definition does not include the usual factor of a 1/2).

    ! Same-spin contribution
    !subroutine b97m_v_css(F,D1F,RhoA,RhoB,D1RhoA,D1RhoB,TauA,TauB)
    !  real(kind=DP),intent(out) :: F
    !  real(kind=DP),intent(out) :: D1F(N_D1_META)
    !  real(kind=DP),intent(in) :: RhoA,RhoB
    !  real(kind=DP),intent(in) :: D1RhoA(3),D1RhoB(3)
    !  real(kind=DP),intent(in) :: TauA,TauB
    !end subroutine b97m_v_css
    call b97m_v_css(Fss,D1Fss,&
         0.5_DP*den,0.5_DP*den,&
         half_grad,half_grad,&
         tau,tau)

    ! Opposite-spin contribution
    !subroutine b97m_v_cos(F,D1F,RhoA,RhoB,D1RhoA,D1RhoB,TauA,TauB)
    !  real(kind=DP),intent(out) :: F
    !  real(kind=DP),intent(out) :: D1F(N_D1_META)
    !  real(kind=DP),intent(in) :: RhoA,RhoB
    !  real(kind=DP),intent(in) :: D1RhoA(3),D1RhoB(3)
    !  real(kind=DP),intent(in) :: TauA,TauB
    !end subroutine b97m_v_cos
    call b97m_v_cos(Fos,D1Fos,&
         0.5_DP*den,0.5_DP*den,&
         half_grad,half_grad,&
         tau,tau)

    ! Evaluate correlation energy as sum of same-spin and opposite spin
    ! contributions
    c_energy = Fss + Fos

    ! Unpack D1F array
    !integer,parameter :: POS_RA  = 1 ! df/dnA
    !integer,parameter :: POS_RB  = 2 ! df/dnB
    !integer,parameter :: POS_GAA = 3 ! df/d(grad(nA).grad(nA))
    !integer,parameter :: POS_GBB = 4 ! df/d(grad(nB).grad(nB))
    !integer,parameter :: POS_GAB = 5 ! df/d(grad(nA).grad(nB)) [ unused ]
    !integer,parameter :: POS_TA  = 6 ! df/dtauA
    !integer,parameter :: POS_TB  = 7 ! df/dtauB
    c_pot = 0.5_DP * ( D1Fss(1) + D1Fos(1) ) +&
         0.5_DP * ( D1Fss(2) + D1Fos(2) )
    ! Use dfdgrad, rather than dfdmgd. We need dfdgrad in the calling routine
    ! (xc_mgga), so getting this directly from df/d(grad(nA).grad(nA)) avoids
    ! the transformation to dfdmgd altogether.
    c_dfdgrad(:) = grad(:) * (0.5_DP * ( D1Fss(3) + D1Fos(3) ) +&
         0.5_DP * ( D1Fss(4) + D1Fos(4) ) )
    c_dfdtau  = D1Fss(6) + D1Fos(6) + D1Fss(7) + D1Fos(7)

    call utils_use_var(mgd)

  end subroutine xc_b97mv_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b97mv_x_point_sp(den1,den2,grad1,grad2,mgd1,mgd2,tau1,tau2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2,&
       x_dfdgrad1,x_dfdgrad2,x_dfdtau1,x_dfdtau2,use_dfdgrad)
    !==============================================================!
    ! This subroutine calculates the the local part of B97M-V      !
    ! exchange in the spin polarized case. It forms a wrapper      !
    ! around the corresponding routine in the xc_b97mv module,     !
    ! which is a modified version of code originally provided by   !
    ! N. Mardirossian.                                             !
    !                                                              !
    ! Note that global dentol and tautol values do not apply to    !
    ! B97M-V, since these are set in xc_b97mv on a per-routine     !
    ! basis.                                                       !
    !                                                              !
    ! B97M-V citation:                                             !
    !       N. Mardirossian & M. Head-Gordon,                      !
    !            J. Chem. Phys. 142, 074111 (2015)                 !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons, 2021, based on stub by             !
    ! James C. Womack, 2016                                        !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert, utils_use_var
    use xc_b97mv, only: N_D1_META, b97m_v_x
    implicit none

    real(kind=dp), intent(in)  :: den1      ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2      ! Density at point for spin 2
    real(kind=dp), intent(in)  :: grad1(3)  ! The charge density gradient for spin 1
                                            ! Cartesian components at a point
    real(kind=dp), intent(in)  :: grad2(3)  ! The charge density gradient for spin 2
                                            ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd1      ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2      ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau1      ! The KE density at point for spin 1
    real(kind=dp), intent(in)  :: tau2      ! The KE density at point for spin 2
    real(kind=dp), intent(out) :: x_energy  ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1    ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2    ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2
    real(kind=dp), intent(out) :: x_dfdgrad1(3) ! Derivative df_{x}/d(grad n_1)
    real(kind=dp), intent(out) :: x_dfdgrad2(3) ! Derivative df_{x}/d(grad n_2)
    real(kind=DP), intent(out) :: x_dfdtau1 ! Derivative df_{x}/dtau spin 1
    real(kind=DP), intent(out) :: x_dfdtau2 ! Derivative df_{x}/dtau spin 2
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential

    ! Local variables
    real(kind=DP) :: D1F(N_D1_META)
    real(kind=dp) :: twotau1, twotau2

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP
    x_dfdgrad1 = 0.0_DP
    x_dfdgrad2 = 0.0_DP
    x_dfdtau1 = 0.0_DP
    x_dfdtau2 = 0.0_DP

    !subroutine b97m_v_x(F,D1F,RhoA,RhoB,D1RhoA,D1RhoB,TauA,TauB)
    !  real(kind=DP),intent(out) :: F
    !  real(kind=DP),intent(out) :: D1F(N_D1_META)
    !  real(kind=DP),intent(in) :: RhoA,RhoB
    !  real(kind=DP),intent(in) :: D1RhoA(3),D1RhoB(3)
    !  real(kind=DP),intent(in) :: TauA,TauB
    !end subroutine b97m_v_x

    ! This routine outputs dfdgrad, rather than dfdmgd
    use_dfdgrad = .true.

    twotau1 = tau1*2.0_dp
    twotau2 = tau2*2.0_dp

    ! Notes:
    ! * In b97m_v_x tau is defined as 2*tau from ONETEP (i.e. the definition
    !   does not include the usual factor of a 1/2).
    call b97m_v_x(x_energy,D1F,&
         den1,den2,&
         grad1, grad2,&
         twotau1,twotau2)

    ! Unpack D1F array
    !integer,parameter :: POS_RA  = 1 ! df/dnA
    !integer,parameter :: POS_RB  = 2 ! df/dnB
    !integer,parameter :: POS_GAA = 3 ! df/d(grad(nA).grad(nA))
    !integer,parameter :: POS_GBB = 4 ! df/d(grad(nB).grad(nB))
    !integer,parameter :: POS_GAB = 5 ! df/d(grad(nA).grad(nB)) [ unused ]
    !integer,parameter :: POS_TA  = 6 ! df/dtauA
    !integer,parameter :: POS_TB  = 7 ! df/dtauB
    x_pot1 = D1F(1)
    x_pot2 = D1F(2)
    ! Use dfdgrad, rather than dfdmgd. We need dfdgrad in the calling routine
    ! (xc_mgga), so getting this directly from df/d(grad(nA).grad(nA)) avoids
    ! the transformation to dfdmgd altogether.
    x_dfdgrad1(:) = D1F(3) * grad1(:) * 2.0_dp
    x_dfdgrad2(:) = D1F(4) * grad2(:) * 2.0_dp
    x_dfdtau1 = D1F(6) * 2.0_dp
    x_dfdtau2 = D1F(7) * 2.0_dp

    call utils_use_var(mgd1)
    call utils_use_var(mgd2)

  end subroutine xc_b97mv_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b97mv_c_point_sp(den,den1,den2,grad1,grad2,mgd,mgd1,mgd2,tau1,tau2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12,&
       c_dfdgrad1,c_dfdgrad2,c_dfdtau1,c_dfdtau2,use_dfdgrad)

    !==============================================================!
    ! This subroutine calculates the the local part of B97M-V      !
    ! correlation in the spin polarized case. It forms a wrapper   !
    ! around the corresponding routine in the xc_b97mv module,     !
    ! which is a modified version of code originally provided by   !
    ! N. Mardirossian.                                             !
    !                                                              !
    ! Note that global dentol and tautol values do not apply to    !
    ! B97M-V, since these are set in xc_b97mv on a per-routine     !
    ! basis.                                                       !
    !                                                              !
    ! B97M-V citation:                                             !
    !       N. Mardirossian & M. Head-Gordon,                      !
    !            J. Chem. Phys. 142, 074111 (2015)                 !
    !                                                              !
    ! The meta-GGA routines include {x,c}_dfdgrad to allow the     !
    ! gradient of the energy per unit volume with respect to the   !
    ! Cartesian components of grad(n) to be directly output.       !
    !--------------------------------------------------------------!
    ! Written by Jolyon Aarons, 2021, based on stub by             !
    ! James C. Womack, 2016.                                       !
    !==============================================================!

    use constants, only: DP, PI
    use utils, only: utils_abort, utils_assert, utils_use_var
    use xc_b97mv, only: N_D1_META, b97m_v_css, b97m_v_cos
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: grad1(3) ! The charge density gradient for spin 1
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: grad2(3) ! The charge density gradient for spin 2
                                           ! Cartesian components at a point
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: tau1     ! The KE density at point for spin 1
    real(kind=dp), intent(in)  :: tau2     ! The KE density at point for spin 2
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1   ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2   ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|
    real(kind=dp), intent(out) :: c_dfdgrad1(3) ! Derivative df_{c}/d(grad n_1)
    real(kind=dp), intent(out) :: c_dfdgrad2(3) ! Derivative df_{c}/d(grad n_2)
    real(kind=DP), intent(out) :: c_dfdtau1 ! Derivative df_{c}/dtau spin 1
    real(kind=DP), intent(out) :: c_dfdtau2 ! Derivative df_{c}/dtau spin 2
    logical, intent(out)       :: use_dfdgrad ! If .true. then calling unit
         ! should use {x,c}_dfdgrad, rather than {x,c}_dfdmgd to evaluate
         ! gradient part of XC potential


    ! Local variables
    real(kind=DP) :: Fss, Fos
    real(kind=DP) :: D1Fss(N_D1_META), D1Fos(N_D1_META)
    real(kind=dp) :: twotau1, twotau2

    ! This routine outputs dfdmgd, rather than dfdgrad
    use_dfdgrad = .true.

    c_energy = 0.0_DP
    c_pot1 = 0.0_DP
    c_pot2 = 0.0_DP
    c_dfdmgd = 0.0_DP
    c_dfdmgd1 = 0.0_DP
    c_dfdmgd2 = 0.0_DP
    c_dfdmgd12 = 0.0_DP
    c_dfdgrad1 = 0.0_DP
    c_dfdgrad2 = 0.0_DP
    c_dfdtau1 = 0.0_DP
    c_dfdtau2 = 0.0_DP

    twotau1 = tau1*2.0_dp
    twotau2 = tau2*2.0_dp

    ! Notes:
    ! * In b97m_v_css, b97m_v_cos, tau is defined as 2*tau from ONETEP
    !   (i.e. the definition does not include the usual factor of a 1/2).

    ! Same-spin contribution
    !subroutine b97m_v_css(F,D1F,RhoA,RhoB,D1RhoA,D1RhoB,TauA,TauB)
    !  real(kind=DP),intent(out) :: F
    !  real(kind=DP),intent(out) :: D1F(N_D1_META)
    !  real(kind=DP),intent(in) :: RhoA,RhoB
    !  real(kind=DP),intent(in) :: D1RhoA(3),D1RhoB(3)
    !  real(kind=DP),intent(in) :: TauA,TauB
    !end subroutine b97m_v_css
    call b97m_v_css(Fss,D1Fss,&
         den1,den2,&
         grad1,grad2,&
         twotau1,twotau2)

    ! Opposite-spin contribution
    !subroutine b97m_v_cos(F,D1F,RhoA,RhoB,D1RhoA,D1RhoB,TauA,TauB)
    !  real(kind=DP),intent(out) :: F
    !  real(kind=DP),intent(out) :: D1F(N_D1_META)
    !  real(kind=DP),intent(in) :: RhoA,RhoB
    !  real(kind=DP),intent(in) :: D1RhoA(3),D1RhoB(3)
    !  real(kind=DP),intent(in) :: TauA,TauB
    !end subroutine b97m_v_cos
    call b97m_v_cos(Fos,D1Fos,&
         den1,den2,&
         grad1,grad2,&
         twotau1,twotau2)

    ! Evaluate correlation energy as sum of same-spin and opposite spin
    ! contributions
    c_energy = Fss + Fos

    ! Unpack D1F array
    !integer,parameter :: POS_RA  = 1 ! df/dnA
    !integer,parameter :: POS_RB  = 2 ! df/dnB
    !integer,parameter :: POS_GAA = 3 ! df/d(grad(nA).grad(nA))
    !integer,parameter :: POS_GBB = 4 ! df/d(grad(nB).grad(nB))
    !integer,parameter :: POS_GAB = 5 ! df/d(grad(nA).grad(nB)) [ unused ]
    !integer,parameter :: POS_TA  = 6 ! df/dtauA
    !integer,parameter :: POS_TB  = 7 ! df/dtauB
    c_pot1 = (D1Fss(1) + D1Fos(1))
    c_pot2 = (D1Fss(2) + D1Fos(2))

    ! Use dfdgrad, rather than dfdmgd. We need dfdgrad in the calling routine
    ! (xc_mgga), so getting this directly from df/d(grad(nA).grad(nA)) avoids
    ! the transformation to dfdmgd altogether.
    c_dfdgrad1(:) = grad1(:) * ( D1Fss(3) + D1Fos(3) ) * 2.0_dp
    c_dfdtau1  = (D1Fss(6) + D1Fos(6)) * 2.0_dp

    c_dfdgrad2(:) = grad2(:) * ( D1Fss(4) + D1Fos(4) ) * 2.0_dp
    c_dfdtau2  = (D1Fss(7) + D1Fos(7)) * 2.0_dp

    call utils_use_var(mgd1)
    call utils_use_var(mgd2)

  end subroutine xc_b97mv_c_point_sp

end module
