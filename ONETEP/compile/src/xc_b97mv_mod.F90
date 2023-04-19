! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!             B97M-V exchange-correlation functional             !
!                                                                !
!                      James C. Womack, 2016                     !
!----------------------------------------------------------------!
!                         Acknowledgement                        !
! Original B97M-V routines provided by Narbe Mardirossian of the !
! Head-Gordon group at University of California, Berkeley.       !
! Permission to use modified versions of these routines in       !
! ONETEP was kindly provided by the original authors.            !
!----------------------------------------------------------------!
!                             Citation                           !
!                 N. Mardirossian & M. Head-Gordon,              !
!                 J. Chem. Phys. 142, 074111 (2015)              !
!----------------------------------------------------------------!
! This module contains modified versions of source code for      !
! evaluating the B97M-V meta-GGA exchange-correlation functional !
! provided by Narbe Mardirossian.                                !
!                                                                !
! The main modifications are to convert the original fixed-form  !
! F77 code to free-form F90/F2003 adhering to ONETEP style       !
! guidelines.                                                    !
!                                                                !
! The routines in this module are intended to be imported by a   !
! a wrapper routine in the xc_funcs module, which will perform   !
! conversions from the original Q-Chem data representation to    !
! ONETEP data representation.                                    !
!================================================================!
module xc_b97mv

  use constants, only: DP, PI
  implicit none
  private

  ! Parameters for imported code
  real(kind=DP),parameter :: VALUE_OF_PI = PI
  integer,parameter :: N_D1_META = 7
  integer,parameter :: POS_RA  = 1 ! df/dnA
  integer,parameter :: POS_RB  = 2 ! df/dnB
  integer,parameter :: POS_GAA = 3 ! df/d(grad(nA).grad(nA))
  integer,parameter :: POS_GBB = 4 ! df/d(grad(nB).grad(nB))
  integer,parameter :: POS_GAB = 5 ! df/d(grad(nA).grad(nB))
  integer,parameter :: POS_TA  = 6 ! df/dtauA
  integer,parameter :: POS_TB  = 7 ! df/dtauB

  public :: N_D1_META
  public :: b97m_v_x
  public :: b97m_v_cos
  public :: b97m_v_css

contains
  subroutine b97m_v_x(F,D1F,RhoA,RhoB,D1RhoA,D1RhoB,TauA,TauB)
    ! TODO Describe function of routine
    ! TODO Describe arguments
    ! Converted from F77 provided by N. Mardirossian to F2003 conforming to
    ! ONETEP style
    ! J. C. Womack, 09/2016
    implicit none
    ! Arguments
    real(kind=DP),intent(out) :: F
    real(kind=DP),intent(out) :: D1F(N_D1_META)
    real(kind=DP),intent(in) :: RhoA,RhoB
    real(kind=DP),intent(in) :: D1RhoA(3),D1RhoB(3)
    real(kind=DP),intent(in) :: TauA,TauB


    ! Parameters
    real(kind=DP), parameter :: tol=1.0e-14_DP
    real(kind=DP), parameter :: E=2.718281828459045_DP
    real(kind=DP), parameter :: &
         cLDA=-(3.0_DP/2.0_DP)*((3.0_DP/(4.0_DP*Pi))**(1.0_DP/3.0_DP))
    real(kind=DP), parameter :: &
         cGGA=1.0_DP/(4.0_DP*(6.0_DP*(Pi**2.0_DP))**(2.0_DP/3.0_DP))
    real(kind=DP), parameter :: &
         cMGGA=(3.0_DP/5.0_DP)*((6.0_DP*(Pi**2.0_DP))**(2.0_DP/3.0_DP))
    real(kind=DP), parameter :: caa=1.0_DP
    real(kind=DP), parameter :: cab=0.416_DP
    real(kind=DP), parameter :: cac=1.308_DP
    real(kind=DP), parameter :: cad=3.07_DP
    real(kind=DP), parameter :: cae=1.901_DP
    real(kind=DP), parameter :: cgamma=0.004_DP

    ! Local variables
    real(kind=DP) :: ExUEG
    real(kind=DP) :: Ex
    real(kind=DP) :: RA
    real(kind=DP) :: RB
    real(kind=DP) :: GA
    real(kind=DP) :: GB
    real(kind=DP) :: TA
    real(kind=DP) :: TB
    real(kind=DP) :: gx
    real(kind=DP) :: ts
    real(kind=DP) :: u
    real(kind=DP) :: w
    real(kind=DP) :: xs
    real(kind=DP) :: dExUEGdRA
    real(kind=DP) :: dExUEGdRB
    real(kind=DP) :: dExdRA
    real(kind=DP) :: dExdRB
    real(kind=DP) :: dExdGA
    real(kind=DP) :: dExdGB
    real(kind=DP) :: dExdTA
    real(kind=DP) :: dExdTB
    real(kind=DP) :: dgxdu
    real(kind=DP) :: dgxdw
    real(kind=DP) :: dtsdRA
    real(kind=DP) :: dtsdRB
    real(kind=DP) :: dtsdTA
    real(kind=DP) :: dtsdTB
    real(kind=DP) :: dudxs
    real(kind=DP) :: dwdts
    real(kind=DP) :: dxsdRA
    real(kind=DP) :: dxsdRB
    real(kind=DP) :: dxsdGA
    real(kind=DP) :: dxsdGB

    F = 0.0_DP
    D1F = 0.0_DP

    RA = max(RhoA,0.0_DP)
    RB = max(RhoB,0.0_DP)
    GA = D1RhoA(1)**2.0_DP+D1RhoA(2)**2.0_DP+D1RhoA(3)**2.0_DP
    GB = D1RhoB(1)**2.0_DP+D1RhoB(2)**2.0_DP+D1RhoB(3)**2.0_DP
    TA = max(TauA,0.0_DP)
    TB = max(TauB,0.0_DP)
    if ((RA.gt.Tol).or.(RB.gt.Tol).or.(TA.gt.Tol).or.(TB.gt.Tol)) then
       if ((RA.gt.Tol).AND.(TA.gt.Tol)) then
          ExUEG = cLDA*RA**(4.0_DP/3.0_DP)
          dExUEGdRA = (4.0_DP*cLDA*RA**(1.0_DP/3.0_DP))/3.0_DP
          xs = GA/RA**(8.0_DP/3.0_DP)
          dxsdRA = (-8.0_DP*GA)/(3.0_DP*RA**(11.0_DP/3.0_DP))
          dxsdGA = RA**(-8.0_DP/3.0_DP)
          ts = (cMGGA*RA**(5.0_DP/3.0_DP))/TA
          dtsdRA = (5.0_DP*cMGGA*RA**(2.0_DP/3.0_DP))/(3.0_DP*TA)
          dtsdTA = -((cMGGA*RA**(5.0_DP/3.0_DP))/TA**2.0_DP)
          u = (cgamma*xs)/(1.0_DP+cgamma*xs)
          dudxs = cgamma/(1.0_DP+cgamma*xs)**2.0_DP
          w = (-1.0_DP+ts)/(1.0_DP+ts)
          dwdts = 2.0_DP/(1.0_DP+ts)**2.0_DP
          gx = caa+cac*u+cae*u**2.0_DP+cab*w+cad*u*w
          dgxdu = cac+2.0_DP*cae*u+cad*w
          dgxdw = cab+cad*u
          Ex = ExUEG*gx
          dExdRA = dgxdw*dtsdRA*dwdts*ExUEG+dgxdu*dudxs*dxsdRA*ExUEG+dExUEGdRA*gx
          dExdGA = dgxdu*dudxs*dxsdGA*ExUEG
          dExdTA = dgxdw*dtsdTA*dwdts*ExUEG
          F = F+Ex
          D1F(POS_RA) = D1F(POS_RA)+dExdRA
          D1F(POS_GAA) = D1F(POS_GAA)+dExdGA
          D1F(POS_TA) = D1F(POS_TA)+dExdTA
       endif
       if ((RB.gt.Tol).and.(TB.gt.Tol)) then
          ExUEG = cLDA*RB**(4.0_DP/3.0_DP)
          dExUEGdRB = (4.0_DP*cLDA*RB**(1.0_DP/3.0_DP))/3.0_DP
          xs = GB/RB**(8.0_DP/3.0_DP)
          dxsdRB = (-8.0_DP*GB)/(3.0_DP*RB**(11.0_DP/3.0_DP))
          dxsdGB = RB**(-8.0_DP/3.0_DP)
          ts = (cMGGA*RB**(5.0_DP/3.0_DP))/TB
          dtsdRB = (5.0_DP*cMGGA*RB**(2.0_DP/3.0_DP))/(3.0_DP*TB)
          dtsdTB = -((cMGGA*RB**(5.0_DP/3.0_DP))/TB**2.0_DP)
          u = (cgamma*xs)/(1.0_DP+cgamma*xs)
          dudxs = cgamma/(1.0_DP+cgamma*xs)**2.0_DP
          w = (-1.0_DP+ts)/(1.0_DP+ts)
          dwdts = 2.0_DP/(1.0_DP+ts)**2.0_DP
          gx = caa+cac*u+cae*u**2.0_DP+cab*w+cad*u*w
          dgxdu = cac+2.0_DP*cae*u+cad*w
          dgxdw = cab+cad*u
          Ex = ExUEG*gx
          dExdRB = dgxdw*dtsdRB*dwdts*ExUEG+dgxdu*dudxs*dxsdRB*ExUEG+dExUEGdRB*gx
          dExdGB = dgxdu*dudxs*dxsdGB*ExUEG
          dExdTB = dgxdw*dtsdTB*dwdts*ExUEG
          F = F+Ex
          D1F(POS_RB) = D1F(POS_RB)+dExdRB
          D1F(POS_GBB) = D1F(POS_GBB)+dExdGB
          D1F(POS_TB) = D1F(POS_TB)+dExdTB
       endif
    endif

  end subroutine b97m_v_x

  subroutine b97m_v_css(F,D1F,RhoA,RhoB,D1RhoA,D1RhoB,TauA,TauB)
    ! TODO Describe function of routine
    ! TODO Describe arguments
    ! Converted from F77 provided by N. Mardirossian to F2003 conforming to
    ! ONETEP style
    ! J. C. Womack, 09/2016
    implicit none
    ! Arguments
    real(kind=DP),intent(out) :: F
    real(kind=DP),intent(out) :: D1F(N_D1_META)
    real(kind=DP),intent(in) :: RhoA,RhoB
    real(kind=DP),intent(in) :: D1RhoA(3),D1RhoB(3)
    real(kind=DP),intent(in) :: TauA,TauB

    ! Parameters
    real(kind=DP), parameter :: Tol=1.0e-14_DP
    real(kind=DP), parameter :: E=2.718281828459045_DP
    real(kind=DP), parameter :: paa=0.0621814_DP/2.0_DP
    real(kind=DP), parameter :: pab=0.0621814_DP/4.0_DP
    real(kind=DP), parameter :: pac=0.0168869_DP
    real(kind=DP), parameter :: pba=0.2137_DP
    real(kind=DP), parameter :: pbb=0.20548_DP
    real(kind=DP), parameter :: pbc=0.11125_DP
    real(kind=DP), parameter :: pca=7.5957_DP
    real(kind=DP), parameter :: pcb=14.1189_DP
    real(kind=DP), parameter :: pcc=10.357_DP
    real(kind=DP), parameter :: pda=3.5876_DP
    real(kind=DP), parameter :: pdb=6.1977_DP
    real(kind=DP), parameter :: pdc=3.6231_DP
    real(kind=DP), parameter :: pea=1.6382_DP
    real(kind=DP), parameter :: peb=3.3662_DP
    real(kind=DP), parameter :: pec=0.88026_DP
    real(kind=DP), parameter :: pfa=0.49294_DP
    real(kind=DP), parameter :: pfb=0.62517_DP
    real(kind=DP), parameter :: pfc=0.49671_DP
    real(kind=DP), parameter :: &
         fppz=4.0_DP/(9.0_DP*(2.0_DP**(1.0_DP/3.0_DP)-1.0_DP))
    real(kind=DP), parameter :: ppa=1.0_DP
    real(kind=DP), parameter :: ppb=1.0_DP
    real(kind=DP), parameter :: ppc=1.0_DP
    real(kind=DP), parameter :: &
         cGGA=1.0_DP/(4.0_DP*(6.0_DP*(Pi**2.0_DP))**(2.0_DP/3.0_DP))
    real(kind=DP), parameter :: &
         cMGGA=(3.0_DP/5.0_DP)*((6.0_DP*(Pi**2.0_DP))**(2.0_DP/3.0_DP))
    real(kind=DP), parameter :: cgamma=0.2_DP
    real(kind=DP), parameter :: caa=1.0_DP
    real(kind=DP), parameter :: cab=-5.668_DP
    real(kind=DP), parameter :: cac=-1.855_DP
    real(kind=DP), parameter :: cad=-20.497_DP
    real(kind=DP), parameter :: cae=-20.364_DP

    ! Local variables
    real(kind=DP) :: EcUEG
    real(kind=DP) :: Ecss
    real(kind=DP) :: RA
    real(kind=DP) :: RB
    real(kind=DP) :: GA
    real(kind=DP) :: GB
    real(kind=DP) :: TA
    real(kind=DP) :: TB
    real(kind=DP) :: GPWb
    real(kind=DP) :: gx
    real(kind=DP) :: rs
    real(kind=DP) :: ts
    real(kind=DP) :: u
    real(kind=DP) :: w
    real(kind=DP) :: xs
    real(kind=DP) :: dECssdRA
    real(kind=DP) :: dECssdRB
    real(kind=DP) :: dECssdGA
    real(kind=DP) :: dECssdGB
    real(kind=DP) :: dECssdTA
    real(kind=DP) :: dECssdTB
    real(kind=DP) :: dEcUEGdGPWb
    real(kind=DP) :: dGPWbdrs
    real(kind=DP) :: dgxdu
    real(kind=DP) :: dgxdw
    real(kind=DP) :: drsdRA
    real(kind=DP) :: drsdRB
    real(kind=DP) :: dtsdRA
    real(kind=DP) :: dtsdRB
    real(kind=DP) :: dtsdTA
    real(kind=DP) :: dtsdTB
    real(kind=DP) :: dudxs
    real(kind=DP) :: dwdts
    real(kind=DP) :: dxsdGA
    real(kind=DP) :: dxsdGB
    real(kind=DP) :: dxsdRA
    real(kind=DP) :: dxsdRB

    F = 0.0_DP
    D1F = 0.0_DP

    RA = max(RhoA,0.0_DP)
    RB = max(RhoB,0.0_DP)
    GA = D1RhoA(1)**2.0_DP+D1RhoA(2)**2.0_DP+D1RhoA(3)**2.0_DP
    GB = D1RhoB(1)**2.0_DP+D1RhoB(2)**2.0_DP+D1RhoB(3)**2.0_DP
    TA = max(TauA,0.0_DP)
    TB = max(TauB,0.0_DP)
    if ((RA.gt.Tol).OR.(RB.gt.Tol).OR.(TA.gt.Tol).OR.(TB.gt.Tol)) THEN
       if ((RA.gt.Tol).AND.(TA.gt.Tol)) THEN
          rs = ((3.0_DP/Pi)**(1.0_DP/3.0_DP)*(RA**(-1.0_DP))**(1.0_DP/3.0_DP))/2.0_DP**(2.0_DP/3.0_DP)
          drsdRA = -((RA**(-1.0_DP))**(4.0_DP/3.0_DP)/(6.0_DP**(2.0_DP/3.0_DP)*Pi**(1.0_DP/3.0_DP)))
          GPWb = -2.0_DP*pab*(1.0_DP+pbb*rs)*log(1.0_DP+1.0_DP/(2.0_DP*pab*(pcb*sqrt(rs)+rs*&
               (pdb+peb*sqrt(rs)+pfb*rs**ppb))))
          dGPWbdrs = (pab*(1.0_DP+pbb*rs)*(pcb+2.0_DP*pdb*sqrt(rs)+3.0_DP*peb*rs+2.0_DP*pfb*&
               rs**(1.0_DP/2.0_DP+ppb)+2.0_DP*pfb*ppb*rs**(1.0_DP/2.0_DP+ppb)))/(rs*(pcb+pdb*&
               sqrt(rs)+peb*rs+pfb*rs**(1.0_DP/2.0_DP+ppb))*(1.0_DP+2.0_DP*pab*sqrt(rs)*&
               (pcb+pdb*sqrt(rs)+peb*rs+pfb*rs**(1.0_DP/2.0_DP+ppb))))-2.0_DP*pab*pbb*&
               log(1.0_DP+1.0_DP/(2.0_DP*pab*(pcb*sqrt(rs)+rs*(pdb+peb*sqrt(rs)+pfb*rs**ppb))))
          EcUEG = GPWb
          dEcUEGdGPWb = 1.0_DP
          xs = GA/RA**(8.0_DP/3.0_DP)
          dxsdRA = (-8.0_DP*GA)/(3.0_DP*RA**(11.0_DP/3.0_DP))
          dxsdGA = RA**(-8.0_DP/3.0_DP)
          ts = (cMGGA*RA**(5.0_DP/3.0_DP))/TA
          dtsdRA = (5.0_DP*cMGGA*RA**(2.0_DP/3.0_DP))/(3.0_DP*TA)
          dtsdTA = -((cMGGA*RA**(5.0_DP/3.0_DP))/TA**2.0_DP)
          u = (cgamma*xs)/(1.0_DP+cgamma*xs)
          dudxs = cgamma/(1.0_DP+cgamma*xs)**2.0_DP
          w = (-1.0_DP+ts)/(1.0_DP+ts)
          dwdts = 2.0_DP/(1.0_DP+ts)**2.0_DP
          gx = caa+cac*u**2.0_DP+cab*w+cad*u**2.0_DP*w**3.0_DP+cae*u**2.0_DP*w**4.0_DP
          dgxdu = 2.0_DP*u*(cac+w**3.0_DP*(cad+cae*w))
          dgxdw = cab+u**2.0_DP*w**2.0_DP*(3.0_DP*cad+4.0_DP*cae*w)
          Ecss = EcUEG*gx*RA
          dEcssdRA = dEcUEGdGPWb*dGPWbdrs*drsdRA*gx*RA+EcUEG*(gx+dgxdw*dtsdRA*dwdts*&
               RA+dgxdu*dudxs*dxsdRA*RA)
          dEcssdGA = dgxdu*dudxs*dxsdGA*EcUEG*RA
          dEcssdTA = dgxdw*dtsdTA*dwdts*EcUEG*RA
          F = F+Ecss
          D1F(POS_RA) = D1F(POS_RA)+dEcssdRA
          D1F(POS_GAA) = D1F(POS_GAA)+dEcssdGA
          D1F(POS_TA) = D1F(POS_TA)+dEcssdTA
       endif
       if ((RB.gt.Tol).AND.(TB.gt.Tol)) THEN
          rs = ((3.0_DP/Pi)**(1.0_DP/3.0_DP)*(RB**(-1.0_DP))**(1.0_DP/3.0_DP))/2.0_DP**(2.0_DP/3.0_DP)
          drsdRB = -((RB**(-1.0_DP))**(4.0_DP/3.0_DP)/(6.0_DP**(2.0_DP/3.0_DP)*Pi**(1.0_DP/3.0_DP)))
          GPWb = -2.0_DP*pab*(1.0_DP+pbb*rs)*log(1.0_DP+1.0_DP/(2.0_DP*pab*(pcb*sqrt(rs)+rs*(pdb+&
               peb*sqrt(rs)+pfb*rs**ppb))))
          dGPWbdrs = (pab*(1.0_DP+pbb*rs)*(pcb+2.0_DP*pdb*sqrt(rs)+3.0_DP*peb*rs+2.0_DP*pfb*&
               rs**(1.0_DP/2.0_DP+ppb)+2.0_DP*pfb*ppb*rs**(1.0_DP/2.0_DP+ppb)))/(rs*(pcb+pdb*&
               sqrt(rs)+peb*rs+pfb*rs**(1.0_DP/2.0_DP+ppb))*(1.0_DP+2.0_DP*pab*sqrt(rs)*&
               (pcb+pdb*sqrt(rs)+peb*rs+pfb*rs**(1.0_DP/2.0_DP+ppb))))-2.0_DP*pab*pbb*&
               log(1.0_DP+1.0_DP/(2.0_DP*pab*(pcb*sqrt(rs)+rs*(pdb+peb*sqrt(rs)+pfb*rs**ppb))))
          EcUEG = GPWb
          dEcUEGdGPWb = 1.0_DP
          xs = GB/RB**(8.0_DP/3.0_DP)
          dxsdRB = (-8.0_DP*GB)/(3.0_DP*RB**(11.0_DP/3.0_DP))
          dxsdGB = RB**(-8.0_DP/3.0_DP)
          ts = (cMGGA*RB**(5.0_DP/3.0_DP))/TB
          dtsdRB = (5.0_DP*cMGGA*RB**(2.0_DP/3.0_DP))/(3.0_DP*TB)
          dtsdTB = -((cMGGA*RB**(5.0_DP/3.0_DP))/TB**2.0_DP)
          u = (cgamma*xs)/(1.0_DP+cgamma*xs)
          dudxs = cgamma/(1.0_DP+cgamma*xs)**2.0_DP
          w = (-1.0_DP+ts)/(1.0_DP+ts)
          dwdts = 2.0_DP/(1.0_DP+ts)**2.0_DP
          gx = caa+cac*u**2.0_DP+cab*w+cad*u**2.0_DP*w**3.0_DP+cae*u**2.0_DP*w**4.0_DP
          dgxdu = 2.0_DP*u*(cac+w**3.0_DP*(cad+cae*w))
          dgxdw = cab+u**2.0_DP*w**2.0_DP*(3.0_DP*cad+4.0_DP*cae*w)
          Ecss = EcUEG*gx*RB
          dEcssdRB = dEcUEGdGPWb*dGPWbdrs*drsdRB*gx*RB+EcUEG*(gx+dgxdw*dtsdRB*&
               dwdts*RB+dgxdu*dudxs*dxsdRB*RB)
          dEcssdGB = dgxdu*dudxs*dxsdGB*EcUEG*RB
          dEcssdTB = dgxdw*dtsdTB*dwdts*EcUEG*RB
          F = F+Ecss
          D1F(POS_RB) = D1F(POS_RB)+dEcssdRB
          D1F(POS_GBB) = D1F(POS_GBB)+dEcssdGB
          D1F(POS_TB) = D1F(POS_TB)+dEcssdTB
       endif
    endif

  end subroutine b97m_v_css

  subroutine b97m_v_cos(F,D1F,RhoA,RhoB,D1RhoA,D1RhoB,TauA,TauB)
    ! TODO Describe function of routine
    ! TODO Describe arguments
    ! Converted from F77 provided by N. Mardirossian to F2003 conforming to
    ! ONETEP style
    ! J. C. Womack, 09/2016
    implicit none
    ! Arguments
    real(kind=DP),intent(out) :: F
    real(kind=DP),intent(out) :: D1F(N_D1_META)
    real(kind=DP),intent(in) :: RhoA,RhoB
    real(kind=DP),intent(in) :: D1RhoA(3),D1RhoB(3)
    real(kind=DP),intent(in) :: TauA,TauB

    ! Parameters
    real(kind=DP), parameter :: Tol=1D-14
    real(kind=DP), parameter :: E=2.718281828459045_DP
    real(kind=DP), parameter :: cLDA=-(3.0_DP/2.0_DP)*((3.0_DP/(4.0_DP*Pi))**(1.0_DP/3.0_DP))
    real(kind=DP), parameter :: paa=0.0621814_DP/2.0_DP
    real(kind=DP), parameter :: pab=0.0621814_DP/4.0_DP
    real(kind=DP), parameter :: pac=0.0168869_DP
    real(kind=DP), parameter :: pba=0.2137_DP
    real(kind=DP), parameter :: pbb=0.20548_DP
    real(kind=DP), parameter :: pbc=0.11125_DP
    real(kind=DP), parameter :: pca=7.5957_DP
    real(kind=DP), parameter :: pcb=14.1189_DP
    real(kind=DP), parameter :: pcc=10.357_DP
    real(kind=DP), parameter :: pda=3.5876_DP
    real(kind=DP), parameter :: pdb=6.1977_DP
    real(kind=DP), parameter :: pdc=3.6231_DP
    real(kind=DP), parameter :: pea=1.6382_DP
    real(kind=DP), parameter :: peb=3.3662_DP
    real(kind=DP), parameter :: pec=0.88026_DP
    real(kind=DP), parameter :: pfa=0.49294_DP
    real(kind=DP), parameter :: pfb=0.62517_DP
    real(kind=DP), parameter :: pfc=0.49671_DP
    real(kind=DP), parameter :: fppz=4.0_DP/(9.0_DP*(2.0_DP**(1.0_DP/3.0_DP)-1.0_DP))
    real(kind=DP), parameter :: ppa=1.0_DP
    real(kind=DP), parameter :: ppb=1.0_DP
    real(kind=DP), parameter :: ppc=1.0_DP
    real(kind=DP), parameter :: &
         cGGA=1.0_DP/(4.0_DP*(6.0_DP*(Pi**2.0_DP))**(2.0_DP/3.0_DP))
    real(kind=DP), parameter :: &
         cMGGA=(3.0_DP/5.0_DP)*((6.0_DP*(Pi**2.0_DP))**(2.0_DP/3.0_DP))
    real(kind=DP), parameter :: cgamma=0.003_DP
    real(kind=DP), parameter :: caa=1.0_DP
    real(kind=DP), parameter :: cab=2.535_DP
    real(kind=DP), parameter :: cac=1.573_DP
    real(kind=DP), parameter :: cad=-6.427_DP
    real(kind=DP), parameter :: cae=-6.298_DP

    ! Local variables
    real(kind=DP) :: EcUEG
    real(kind=DP) :: EcUEGa
    real(kind=DP) :: EcUEGb
    real(kind=DP) :: Ecos
    real(kind=DP) :: fpol
    real(kind=DP) :: RA
    real(kind=DP) :: RB
    real(kind=DP) :: GA
    real(kind=DP) :: GB
    real(kind=DP) :: TA
    real(kind=DP) :: TB
    real(kind=DP) :: GAB
    real(kind=DP) :: GPWa
    real(kind=DP) :: GPWb
    real(kind=DP) :: GPWba
    real(kind=DP) :: GPWbb
    real(kind=DP) :: GPWc
    real(kind=DP) :: gx
    real(kind=DP) :: rs
    real(kind=DP) :: rsa
    real(kind=DP) :: rsb
    real(kind=DP) :: ts
    real(kind=DP) :: u
    real(kind=DP) :: w
    real(kind=DP) :: xs
    real(kind=DP) :: zeta
    real(kind=DP) :: dEcosdRA
    real(kind=DP) :: dEcosdRB
    real(kind=DP) :: dEcosdGA
    real(kind=DP) :: dEcosdGB
    real(kind=DP) :: dEcosdTA
    real(kind=DP) :: dEcosdTB
    real(kind=DP) :: dEcUEGadGPWba
    !real(kind=DP) :: dEcUEGadGPWbb <-- unused
    real(kind=DP) :: dEcUEGdfpol
    real(kind=DP) :: dEcUEGdGPWa
    real(kind=DP) :: dEcUEGdGPWb
    real(kind=DP) :: dEcUEGdGPWc
    real(kind=DP) :: dEcUEGbdGPWbb
    real(kind=DP) :: dEcUEGdzeta
    real(kind=DP) :: dfpoldzeta
    real(kind=DP) :: dGPWadrs
    real(kind=DP) :: dGPWbadrsa
    !real(kind=DP) :: dGPWbadrsb <-- unused
    real(kind=DP) :: dGPWbdrs
    real(kind=DP) :: dGPWcdrs
    real(kind=DP) :: dGPWbbdrsb
    real(kind=DP) :: dgxdu
    real(kind=DP) :: dgxdw
    real(kind=DP) :: drsadRA
    ! real(kind=DP) :: drsadRB <-- unused
    real(kind=DP) :: drsdRA
    real(kind=DP) :: drsdRB
    real(kind=DP) :: drsbdRB
    real(kind=DP) :: dtsdRA
    real(kind=DP) :: dtsdRB
    real(kind=DP) :: dtsdTA
    real(kind=DP) :: dtsdTB
    real(kind=DP) :: dudxs
    real(kind=DP) :: dwdts
    real(kind=DP) :: dxsdRA
    real(kind=DP) :: dxsdRB
    real(kind=DP) :: dxsdGA
    real(kind=DP) :: dxsdGB
    real(kind=DP) :: dzetadRA
    real(kind=DP) :: dzetadRB

    ! JCW: Zero energy and derivatives, unlike in original code.
    ! JCW: This make the behaviour of the same- and opposite-spin routines
    ! JCW: consistent.
    F = 0.0_DP
    D1F = 0.0_DP

    RA = max(RhoA,0.0_DP)
    RB = max(RhoB,0.0_DP)
    GA = D1RhoA(1)**2.0_DP+D1RhoA(2)**2.0_DP+D1RhoA(3)**2.0_DP
    GB = D1RhoB(1)**2.0_DP+D1RhoB(2)**2.0_DP+D1RhoB(3)**2.0_DP
    GAB = D1RhoA(1)*D1RhoB(1)+D1RhoA(2)*D1RhoB(2)+D1RhoA(3)*D1RhoB(3)
    TA = max(TauA,0.0_DP)
    TB = max(TauB,0.0_DP)
    if ((RA.gt.Tol).or.(RB.gt.Tol).or.(TA.gt.Tol).or.(TB.gt.Tol)) then
       if ((RA.gt.Tol).and.(RB.gt.Tol).and.(TA.gt.Tol).and.(TB.gt.Tol)) then
          rs = ((3.0_DP/Pi)**(1.0_DP/3.0_DP)*((RA+RB)**(-1.0_DP))**(1.0_DP/3.0_DP))/&
               2.0_DP**(2.0_DP/3.0_DP)
          drsdRA = -(((RA+RB)**(-1.0_DP))**(4.0_DP/3.0_DP)/(6.0_DP**(2.0_DP/3.0_DP)*&
               Pi**(1.0_DP/3.0_DP)))
          drsdRB = -(((RA+RB)**(-1.0_DP))**(4.0_DP/3.0_DP)/(6.0_DP**(2.0_DP/3.0_DP)*&
               Pi**(1.0_DP/3.0_DP)))
          zeta = (RA-RB)/(RA+RB)
          dzetadRA = (2.0_DP*RB)/(RA+RB)**2.0_DP
          dzetadRB = (-2.0_DP*RA)/(RA+RB)**2.0_DP
          if (abs(zeta).lt.Tol) then
             zeta = 0.0_DP
          end if
          GPWa = -2.0_DP*paa*(1.0_DP+pba*rs)*log(1.0_DP+1.0_DP/(2.0_DP*paa*(pca*&
               sqrt(rs)+rs*(pda+pea*sqrt(rs)+pfa*rs**ppa))))
          dGPWadrs = &
               (paa*(1.0_DP+pba*rs)*(pca+2.0_DP*pda*sqrt(rs)+3.0_DP*pea*rs+2.0_DP*pfa*&
               rs**(1.0_DP/2.0_DP+ppa)+2.0_DP*pfa*ppa*rs**(1.0_DP/2.0_DP+ppa)))/(rs*&
               (pca+pda*sqrt(rs)+pea*rs+pfa*rs**(1.0_DP/2.0_DP+ppa))*(1.0_DP+2.0_DP*&
               paa*sqrt(rs)*(pca+pda*sqrt(rs)+pea*rs+pfa*rs**(1.0_DP/2.0_DP+ppa))))-&
               2.0_DP*paa*pba*log(1.0_DP+1.0_DP/(2.0_DP*paa*(pca*sqrt(rs)+rs*(pda+&
               pea*sqrt(rs)+pfa*rs**ppa))))
          GPWb = -2.0_DP*pab*(1.0_DP+pbb*rs)*log(1.0_DP+1.0_DP/(2.0_DP*pab*(pcb*&
              sqrt(rs)+rs*(pdb+peb*sqrt(rs)+pfb*rs**ppb))))
          dGPWbdrs = &
               (pab*(1.0_DP+pbb*rs)*(pcb+2.0_DP*pdb*sqrt(rs)+3.0_DP*peb*rs+2.0_DP*pfb*&
               rs**(1.0_DP/2.0_DP+ppb)+2.0_DP*pfb*ppb*rs**(1.0_DP/2.0_DP+ppb)))/(rs*&
               (pcb+pdb*sqrt(rs)+peb*rs+pfb*rs**(1.0_DP/2.0_DP+ppb))*(1.0_DP+2.0_DP*&
               pab*sqrt(rs)*(pcb+pdb*sqrt(rs)+peb*rs+pfb*rs**(1.0_DP/2.0_DP+ppb))))-&
               2.0_DP*pab*pbb*log(1.0_DP+1.0_DP/(2.0_DP*pab*(pcb*sqrt(rs)+rs*(pdb+&
               peb*sqrt(rs)+pfb*rs**ppb))))
          GPWc = -2.0_DP*pac*(1.0_DP+pbc*rs)*log(1.0_DP+1.0_DP/(2.0_DP*pac*(pcc*&
               sqrt(rs)+rs*(pdc+pec*sqrt(rs)+pfc*rs**ppc))))
          dGPWcdrs = &
               (pac*(1.0_DP+pbc*rs)*(pcc+2.0_DP*pdc*sqrt(rs)+3.0_DP*pec*rs+2.0_DP*pfc*&
               rs**(1.0_DP/2.0_DP+ppc)+2.0_DP*pfc*ppc*rs**(1.0_DP/2.0_DP+ppc)))/(rs*&
               (pcc+pdc*sqrt(rs)+pec*rs+pfc*rs**(1.0_DP/2.0_DP+ppc))*(1.0_DP+2.0_DP*&
               pac*sqrt(rs)*(pcc+pdc*sqrt(rs)+pec*rs+pfc*rs**(1.0_DP/2.0_DP+ppc))))-&
               2.0_DP*pac*pbc*log(1.0_DP+1.0_DP/(2.0_DP*pac*(pcc*sqrt(rs)+rs*(pdc+&
               pec*sqrt(rs)+pfc*rs**ppc))))
          fpol = (-2.0_DP+(1.0_DP-zeta)**(4.0_DP/3.0_DP)+(1.0_DP+zeta)**(4.0_DP/3.0_DP))/&
               (-2.0_DP+2.0_DP*2.0_DP**(1.0_DP/3.0_DP))
          dfpoldzeta = (2.0_DP*(1.0_DP-zeta)**(1.0_DP/3.0_DP)-2.0_DP*&
               (1.0_DP+zeta)**(1.0_DP/3.0_DP))/(3.0_DP-3.0_DP*2.0_DP**(1.0_DP/3.0_DP))
          EcUEG = GPWa+fpol*(-GPWa+GPWb)*zeta**4.0_DP+&
               (fpol*GPWc*(-1.0_DP+zeta**4.0_DP))/fppz
          dEcUEGdGPWa = 1.0_DP-fpol*zeta**4.0_DP
          dEcUEGdGPWb = fpol*zeta**4.0_DP
          dEcUEGdGPWc = (fpol*(-1.0_DP+zeta**4.0_DP))/fppz
          dEcUEGdfpol = (fppz*(-GPWa+GPWb)*zeta**4.0_DP+GPWc*(-1.0_DP+zeta**4.0_DP))/fppz
          dEcUEGdzeta = (4.0_DP*fpol*(fppz*(-GPWa+GPWb)+GPWc)*zeta**3.0_DP)/fppz
          rsa = ((3.0_DP/Pi)**(1.0_DP/3.0_DP)*(RA**(-1.0_DP))**(1.0_DP/3.0_DP))/&
               2.0_DP**(2.0_DP/3.0_DP)
          drsadRA = -((RA**(-1.0_DP))**(4.0_DP/3.0_DP)/(6.0_DP**(2.0_DP/3.0_DP)*&
               Pi**(1.0_DP/3.0_DP)))
          GPWba = -2.0_DP*pab*(1.0_DP+pbb*rsa)*log(1.0_DP+1.0_DP/(2.0_DP*pab*(pcb*&
               sqrt(rsa)+rsa*(pdb+peb*sqrt(rsa)+pfb*rsa))))
          dGPWbadrsa = &
               (pab*(1.0_DP+pbb*rsa)*(pcb+2.0_DP*pdb*sqrt(rsa)+3.0_DP*peb*rsa+&
               4.0_DP*pfb*rsa**(3.0_DP/2.0_DP)))/(rsa*(pcb+pdb*sqrt(rsa)+peb*rsa+&
               pfb*rsa**(3.0_DP/2.0_DP))*(1.0_DP+2.0_DP*pab*(pcb*sqrt(rsa)+rsa*&
               (pdb+peb*sqrt(rsa)+pfb*rsa))))-2.0_DP*pab*pbb*log(1.0_DP+1.0_DP/&
               (2.0_DP*pab*(pcb*sqrt(rsa)+rsa*(pdb+peb*sqrt(rsa)+pfb*rsa))))
          EcUEGa = GPWba
          dEcUEGadGPWba = 1.0_DP
          rsb = ((3.0_DP/Pi)**(1.0_DP/3.0_DP)*(RB**(-1.0_DP))**(1.0_DP/3.0_DP))/&
               2.0_DP**(2.0_DP/3.0_DP)
          drsbdRB = -((RB**(-1.0_DP))**(4.0_DP/3.0_DP)/(6.0_DP**(2.0_DP/3.0_DP)*&
               Pi**(1.0_DP/3.0_DP)))
          GPWbb = -2.0_DP*pab*(1.0_DP+pbb*rsb)*log(1.0_DP+1.0_DP/(2.0_DP*pab*(pcb*&
               sqrt(rsb)+rsb*(pdb+peb*sqrt(rsb)+pfb*rsb))))
          dGPWbbdrsb = &
               (pab*(1.0_DP+pbb*rsb)*(pcb+2.0_DP*pdb*sqrt(rsb)+3.0_DP*peb*rsb+&
               4.0_DP*pfb*rsb**(3.0_DP/2.0_DP)))/(rsb*(pcb+pdb*sqrt(rsb)+peb*rsb+&
               pfb*rsb**(3.0_DP/2.0_DP))*(1.0_DP+2.0_DP*pab*(pcb*sqrt(rsb)+rsb*&
               (pdb+peb*sqrt(rsb)+pfb*rsb))))-2.0_DP*pab*pbb*log(1.0_DP+1.0_DP/&
               (2.0_DP*pab*(pcb*sqrt(rsb)+rsb*(pdb+peb*sqrt(rsb)+pfb*rsb))))
          EcUEGb = GPWbb
          dEcUEGbdGPWbb = 1.0_DP
          xs = GA/RA**(8.0_DP/3.0_DP)+GB/RB**(8.0_DP/3.0_DP)
          dxsdRA = (-8.0_DP*GA)/(3.0_DP*RA**(11.0_DP/3.0_DP))
          dxsdRB = (-8.0_DP*GB)/(3.0_DP*RB**(11.0_DP/3.0_DP))
          dxsdGA = RA**(-8.0_DP/3.0_DP)
          dxsdGB = RB**(-8.0_DP/3.0_DP)
          ts = (cMGGA*(RA**(5.0_DP/3.0_DP)/TA+RB**(5.0_DP/3.0_DP)/TB))/2.0_DP
          dtsdRA = (5.0_DP*cMGGA*RA**(2.0_DP/3.0_DP))/(6.0_DP*TA)
          dtsdRB = (5.0_DP*cMGGA*RB**(2.0_DP/3.0_DP))/(6.0_DP*TB)
          dtsdTA = -(cMGGA*RA**(5.0_DP/3.0_DP))/(2.0_DP*TA**2.0_DP)
          dtsdTB = -(cMGGA*RB**(5.0_DP/3.0_DP))/(2.0_DP*TB**2.0_DP)
          u = (cgamma*xs)/(1.0_DP+cgamma*xs)
          dudxs = cgamma/(1.0_DP+cgamma*xs)**2.0_DP
          w = (-1.0_DP+ts)/(1.0_DP+ts)
          dwdts = 2.0_DP/(1.0_DP+ts)**2.0_DP
          gx = caa+cac*u+cae*u**3.0_DP+cab*w+cad*u**2.0_DP*w**3.0_DP
          dgxdu = cac+u*(3.0_DP*cae*u+2.0_DP*cad*w**3.0_DP)
          dgxdw = cab+3.0_DP*cad*u**2.0_DP*w**2.0_DP
          Ecos = gx*(-(EcUEGa*RA)-EcUEGb*RB+EcUEG*(RA+RB))
          dEcosdRA = gx*(EcUEG-EcUEGa-dEcUEGadGPWba*dGPWbadrsa*drsadRA*RA+&
               (dEcUEGdGPWa*dGPWadrs*drsdRA+dEcUEGdGPWb*dGPWbdrs*drsdRA+&
               dEcUEGdGPWc*dGPWcdrs*drsdRA+dEcUEGdzeta*dzetadRA+dEcUEGdfpol*&
               dfpoldzeta*dzetadRA)*(RA+RB))+(dgxdw*dtsdRA*dwdts+dgxdu*dudxs*dxsdRA)*&
               (-(EcUEGa*RA)-EcUEGb*RB+EcUEG*(RA+RB))
          dEcosdRB = gx*(EcUEG-EcUEGb-dEcUEGbdGPWbb*dGPWbbdrsb*drsbdRB*RB+&
               (dEcUEGdGPWa*dGPWadrs*drsdRB+dEcUEGdGPWb*dGPWbdrs*drsdRB+&
               dEcUEGdGPWc*dGPWcdrs*drsdRB+dEcUEGdzeta*dzetadRB+dEcUEGdfpol*&
               dfpoldzeta*dzetadRB)*(RA+RB))+(dgxdw*dtsdRB*dwdts+dgxdu*dudxs*dxsdRB)*&
               (-(EcUEGa*RA)-EcUEGb*RB+EcUEG*(RA+RB))
          dEcosdGA = dgxdu*dudxs*dxsdGA*(-(EcUEGa*RA)-EcUEGb*RB+EcUEG*(RA+RB))
          dEcosdGB = dgxdu*dudxs*dxsdGB*(-(EcUEGa*RA)-EcUEGb*RB+EcUEG*(RA+RB))
          dEcosdTA = dgxdw*dtsdTA*dwdts*(-(EcUEGa*RA)-EcUEGb*RB+EcUEG*(RA+RB))
          dEcosdTB = dgxdw*dtsdTB*dwdts*(-(EcUEGa*RA)-EcUEGb*RB+EcUEG*(RA+RB))
          F = F+Ecos
          D1F(POS_RA) = D1F(POS_RA)+dEcosdRA
          D1F(POS_RB) = D1F(POS_RB)+dEcosdRB
          D1F(POS_GAA) = D1F(POS_GAA)+dEcosdGA
          D1F(POS_GBB) = D1F(POS_GBB)+dEcosdGB
          D1F(POS_TA) = D1F(POS_TA)+dEcosdTA
          D1F(POS_TB) = D1F(POS_TB)+dEcosdTB
       end if
       if ((RA.gt.Tol).and.(RB.lt.Tol).and.(TA.gt.Tol).and.(TB.lt.Tol)) then
          Ecos = 0.0_DP
          dEcosdRA = 0.0_DP
          dEcosdGA = 0.0_DP
          dEcosdTA = 0.0_DP
          F = F+Ecos
          D1F(POS_RA) = D1F(POS_RA)+dEcosdRA
          D1F(POS_GAA) = D1F(POS_GAA)+dEcosdGA
          D1F(POS_TA) = D1F(POS_TA)+dEcosdTA
       end if
       if ((RA.lt.Tol).and.(RB.gt.Tol).and.(TA.lt.Tol).and.(TB.gt.Tol)) then
          Ecos = 0.0_DP
          dEcosdRB = 0.0_DP
          dEcosdGB = 0.0_DP
          dEcosdTB = 0.0_DP
          F = F+Ecos
          D1F(POS_RB) = D1F(POS_RB)+dEcosdRB
          D1F(POS_GBB) = D1F(POS_GBB)+dEcosdGB
          D1F(POS_TB) = D1F(POS_TB)+dEcosdTB
       end if
    end if

  end subroutine b97m_v_cos

end module
