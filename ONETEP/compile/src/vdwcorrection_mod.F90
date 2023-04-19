! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!          Van der Waals empirical correction module             !
!                                                                !
! This module contains subroutines for calculating a Van der     !
! Waals energy correction as a sum of damped London potentials.  !
!----------------------------------------------------------------!
! Written by Quintin Hill in 2007/8 with assistance from         !
! Chris-Kriton Skylaris.                                         !
! Forces added July 2008 by Quintin Hill                         !
! Grimme's -D2 correction added February 2013 by Max Phipps      !
!================================================================!

module vdwcorrection

  use constants, only: DP

  implicit none

  private

  ! qoh: Structure to store parameters in

  TYPE VDWPARAMETERS
     !qoh: C_6 coefficients
     REAL(kind=DP), DIMENSION(109) :: C6COEFF
     !qoh: R_0 values
     REAL(kind=DP), DIMENSION(109) :: RADZERO
     !qoh: damping coefficient
     REAL(kind=DP), DIMENSION(5) :: DCOEFF
     !qoh: Effective number of electrons
     REAL(kind=DP), DIMENSION(109) :: NEFF
     !mjsp: Grimme's s6 scaling coefficient
     REAL(kind=DP) :: S6COEFF
  END TYPE VDWPARAMETERS

  public :: vdwcorrection_calculate_energy
  public :: vdwcorrection_calculate_forces
  public :: vdwcorrection_override_alloc
  public :: vdwcorrection_override_dealloc

  integer :: dispersion_index   ! mjsp: dispersion correction setting index

  ! qoh: Share array of Van der Waals parameters accross module
  type(VDWPARAMETERS) :: vdwparams ! Array of parameters
  real(kind=DP), save, public :: pub_dispersion_energy
  ! The dispersion correction energy
  character(len=80), save, allocatable :: vdwparam_override(:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_initializeparams

    !==================================================================!
    ! This subroutine populates the array of Van der Waals parameters  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Some unoptimised parameters (e.g. for P) were taken from Halgren.!
    ! Journal of the American Chemical Society 114(20), 7827-7843, 1992!
    ! - Modified for Grimme's -D2 correction Feb. 2013 Max Phipps      !
    !   J Comput Chem. 27(15), 1787-99, 2006                           !
    !==================================================================!

#ifdef LIBXC
    ! JCW: Used in F2003 Libxc interface
    use iso_c_binding, only: c_int
#endif
    use bibliography, only: bibliography_cite
    use comms, only: pub_on_root
    use constants, only: VERBOSE, stdout
    use rundat, only: pub_xc_functional, pub_dispersion, pub_output_detail
    use utils, only: utils_banner
#ifdef LIBXC
    use xc_f03_lib_m  !! External dependency
    use rundat, only: pub_libxc_c_func_id, pub_libxc_x_func_id
#endif

    implicit none

    ! ndmh: Local Variables
    character(len=80) :: xc_func

    call bibliography_cite('DISP')

    ! ndmh: set up xc_func, overriding parameter for libxc functionals
    xc_func = pub_xc_functional
#ifdef LIBXC
    if (xc_func=='LIBXC') then
       if ((int(pub_libxc_x_func_id,kind=c_int)==XC_GGA_X_B88) .and. &
            (int(pub_libxc_c_func_id,kind=c_int)==XC_GGA_C_LYP)) then
          xc_func = 'BLYP'
       else if ((int(pub_libxc_x_func_id,kind=c_int)==XC_GGA_X_PBE) .and. &
            (int(pub_libxc_c_func_id,kind=c_int)==XC_GGA_C_PBE)) then
          xc_func = 'PBE'
       else if ((int(pub_libxc_x_func_id,kind=c_int)==XC_GGA_X_PW91) .and. &
            (int(pub_libxc_c_func_id,kind=c_int)==XC_GGA_C_PW91)) then
          xc_func = 'PW91'
       else if ((int(pub_libxc_x_func_id,kind=c_int)==XC_GGA_X_PBE_R) .and. &
            (int(pub_libxc_c_func_id,kind=c_int)==XC_GGA_C_PBE)) then
          xc_func = 'REVPBE'
       else if ((int(pub_libxc_x_func_id,kind=c_int)==XC_GGA_X_RPBE) .and. &
            (int(pub_libxc_c_func_id,kind=c_int)==XC_GGA_C_PBE)) then
          xc_func = 'RPBE'
       else if ((int(pub_libxc_x_func_id,kind=c_int)==XC_GGA_XC_XLYP)) then
          xc_func = 'XLYP'
       end if
    end if
#endif

    ! mjsp: convert pub_dispersion string to dispersion_index integer
    if (pub_output_detail>=VERBOSE .and. pub_on_root) then
       ! mjsp: start of dispersion printout block
       write(stdout,'(a)') utils_banner('=','Empirical dispersion information')
    end if
    select case (pub_dispersion)
    case('ELSTNER','1')
      dispersion_index = 1
      if (pub_on_root) write(stdout,'(2a)') "Model : Elstner"
    case('WUYANG1','2')
      dispersion_index = 2
      if (pub_on_root) write(stdout,'(2a)') "Model : Wu and Yang (First Damping Function)"
    case('WUYANG2','3')
      dispersion_index = 3
      if (pub_on_root) write(stdout,'(2a)') "Model : Wu and Yang (Second Damping Function)"
    case('GRIMMED2','4')
      dispersion_index = 4
      if (pub_on_root) write(stdout,'(2a)') "Model : Grimme D2"
    case default
      dispersion_index = 0
    end select


    ! qoh: Initialise the vdwparams array
    vdwparams%c6coeff = 0.0_DP
    vdwparams%neff    = 0.0_DP
    ! mjsp: note s6coeff is updated only in case(4) D2 correction
    vdwparams%s6coeff = 1.0_DP
    ! qoh: Setting to 1 prevents division by zero
    vdwparams%radzero = 1.0_DP

    ! qoh: Populate the vdwparams array
    select case (dispersion_index)
    ! mjsp: Populate vdwparams with qoh's optimized dispersion values
    case(1:3)
        select case (xc_func)
        case('PBE')
           pbedamp: select case (dispersion_index)

           case(1) pbedamp

              vdwparams%dcoeff(1)=3.2607_DP
              vdwparams%dcoeff(2)=3.5400_DP
              vdwparams%dcoeff(3)=23.0000_DP
              vdwparams%c6coeff(1)=2.9239_DP
              vdwparams%c6coeff(6)=27.3561_DP
              vdwparams%c6coeff(7)=19.5089_DP
              vdwparams%c6coeff(8)=11.7697_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0368_DP
              vdwparams%radzero(1)=2.9635_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.3610_DP
              vdwparams%radzero(7)=3.5136_DP
              vdwparams%radzero(8)=3.7294_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=5.1033_DP

           case(2) pbedamp

              vdwparams%dcoeff(1)=3.0000_DP
              vdwparams%dcoeff(2)=3.5400_DP
              vdwparams%dcoeff(3)=23.0000_DP
              vdwparams%c6coeff(1)=2.8450_DP
              vdwparams%c6coeff(6)=27.3200_DP
              vdwparams%c6coeff(7)=19.4800_DP
              vdwparams%c6coeff(8)=11.7600_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0232_DP
              vdwparams%radzero(1)=3.0000_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.8000_DP
              vdwparams%radzero(7)=3.8000_DP
              vdwparams%radzero(8)=3.8000_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=5.5998_DP

           case(3) pbedamp

              vdwparams%dcoeff(1)=3.0000_DP
              vdwparams%dcoeff(2)=3.5400_DP
              vdwparams%dcoeff(3)=23.0025_DP
              vdwparams%c6coeff(1)=2.9865_DP
              vdwparams%c6coeff(6)=27.3784_DP
              vdwparams%c6coeff(7)=19.5223_DP
              vdwparams%c6coeff(8)=11.7733_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0388_DP
              vdwparams%radzero(1)=2.8996_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.0001_DP
              vdwparams%radzero(7)=3.2659_DP
              vdwparams%radzero(8)=3.6630_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=4.7948_DP

           end select pbedamp

        case('RPBE')
           rpbedamp: select case (dispersion_index)

           case(1) rpbedamp

              vdwparams%dcoeff(1)=3.4967_DP
              vdwparams%dcoeff(2)=3.5400_DP
              vdwparams%dcoeff(3)=23.0000_DP
              vdwparams%c6coeff(1)=3.2054_DP
              vdwparams%c6coeff(6)=27.4608_DP
              vdwparams%c6coeff(7)=19.5765_DP
              vdwparams%c6coeff(8)=11.7926_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0465_DP
              vdwparams%radzero(1)=2.8062_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.0375_DP
              vdwparams%radzero(7)=3.2484_DP
              vdwparams%radzero(8)=3.6109_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=4.0826_DP

           case(2) rpbedamp

              vdwparams%dcoeff(1)=3.0000_DP
              vdwparams%dcoeff(2)=3.9308_DP
              vdwparams%dcoeff(3)=23.0000_DP
              vdwparams%c6coeff(1)=3.2290_DP
              vdwparams%c6coeff(6)=27.4707_DP
              vdwparams%c6coeff(7)=19.5922_DP
              vdwparams%c6coeff(8)=11.8007_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0388_DP
              vdwparams%radzero(1)=2.9406_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.4751_DP
              vdwparams%radzero(7)=3.6097_DP
              vdwparams%radzero(8)=3.7393_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=4.8007_DP

           case(3) rpbedamp

              vdwparams%dcoeff(1)=3.0000_DP
              vdwparams%dcoeff(2)=3.5400_DP
              vdwparams%dcoeff(3)=23.0033_DP
              vdwparams%c6coeff(1)=3.0210_DP
              vdwparams%c6coeff(6)=27.3845_DP
              vdwparams%c6coeff(7)=19.5208_DP
              vdwparams%c6coeff(8)=11.7723_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0437_DP
              vdwparams%radzero(1)=2.8667_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=2.9914_DP
              vdwparams%radzero(7)=3.2931_DP
              vdwparams%radzero(8)=3.6758_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=3.9912_DP

           end select rpbedamp

        case('REVPBE')
           revpbedamp: select case (dispersion_index)

           case(1) revpbedamp
              vdwparams%dcoeff(1)=3.4962_DP
              vdwparams%dcoeff(2)=3.5400_DP
              vdwparams%dcoeff(3)=23.0000_DP
              vdwparams%c6coeff(1)=3.2167_DP
              vdwparams%c6coeff(6)=27.4616_DP
              vdwparams%c6coeff(7)=19.5760_DP
              vdwparams%c6coeff(8)=11.7927_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0475_DP
              vdwparams%radzero(1)=2.7985_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.0398_DP
              vdwparams%radzero(7)=3.2560_DP
              vdwparams%radzero(8)=3.6122_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=3.9811_DP

           case(2) revpbedamp

              vdwparams%dcoeff(1)=3.0000_DP
              vdwparams%dcoeff(2)=3.8282_DP
              vdwparams%dcoeff(3)=23.0000_DP
              vdwparams%c6coeff(1)=3.1160_DP
              vdwparams%c6coeff(6)=27.4184_DP
              vdwparams%c6coeff(7)=19.5513_DP
              vdwparams%c6coeff(8)=11.7867_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0389_DP
              vdwparams%radzero(1)=2.9540_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.5630_DP
              vdwparams%radzero(7)=3.6696_DP
              vdwparams%radzero(8)=3.7581_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=4.7980_DP

           case(3) revpbedamp

              vdwparams%dcoeff(1)=3.0000_DP
              vdwparams%dcoeff(2)=3.5400_DP
              vdwparams%dcoeff(3)=23.0034_DP
              vdwparams%c6coeff(1)=3.0258_DP
              vdwparams%c6coeff(6)=27.3854_DP
              vdwparams%c6coeff(7)=19.5210_DP
              vdwparams%c6coeff(8)=11.7724_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0435_DP
              vdwparams%radzero(1)=2.8623_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=2.9923_DP
              vdwparams%radzero(7)=3.2952_DP
              vdwparams%radzero(8)=3.6750_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=3.9910_DP

           end select revpbedamp

        case('PW91')
           pw91damp: select case (dispersion_index)

           case(1) pw91damp

              vdwparams%dcoeff(1)=3.2106_DP
              vdwparams%dcoeff(2)=3.5400_DP
              vdwparams%dcoeff(3)=23.0000_DP
              vdwparams%c6coeff(1)=2.8701_DP
              vdwparams%c6coeff(6)=27.3422_DP
              vdwparams%c6coeff(7)=19.5030_DP
              vdwparams%c6coeff(8)=11.7670_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0346_DP
              vdwparams%radzero(1)=3.0013_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.4423_DP
              vdwparams%radzero(7)=3.5445_DP
              vdwparams%radzero(8)=3.7444_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=5.4111_DP

           case(2) pw91damp

              vdwparams%dcoeff(1)=3.0000_DP
              vdwparams%dcoeff(2)=3.3102_DP
              vdwparams%dcoeff(3)=23.0000_DP
              vdwparams%c6coeff(1)=2.5834_DP
              vdwparams%c6coeff(6)=27.2743_DP
              vdwparams%c6coeff(7)=19.4633_DP
              vdwparams%c6coeff(8)=11.7472_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0233_DP
              vdwparams%radzero(1)=3.0759_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.9644_DP
              vdwparams%radzero(7)=3.8390_DP
              vdwparams%radzero(8)=3.8330_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=5.6250_DP

           case(3) pw91damp

              vdwparams%dcoeff(1)=3.0000_DP
              vdwparams%dcoeff(2)=3.5400_DP
              vdwparams%dcoeff(3)=23.0004_DP
              vdwparams%c6coeff(1)=2.9252_DP
              vdwparams%c6coeff(6)=27.3608_DP
              vdwparams%c6coeff(7)=19.5150_DP
              vdwparams%c6coeff(8)=11.7706_DP
              vdwparams%c6coeff(15)=190.5033_DP
              vdwparams%c6coeff(16)=161.0381_DP
              vdwparams%radzero(1)=2.9543_DP
              vdwparams%radzero(2)=3.0000_DP
              vdwparams%radzero(3)=3.8000_DP
              vdwparams%radzero(4)=3.8000_DP
              vdwparams%radzero(5)=3.8000_DP
              vdwparams%radzero(6)=3.0785_DP
              vdwparams%radzero(7)=3.3119_DP
              vdwparams%radzero(8)=3.6858_DP
              vdwparams%radzero(9)=3.8000_DP
              vdwparams%radzero(10)=3.8000_DP
              vdwparams%radzero(11)=4.8000_DP
              vdwparams%radzero(12)=4.8000_DP
              vdwparams%radzero(13)=4.8000_DP
              vdwparams%radzero(14)=4.8000_DP
              vdwparams%radzero(15)=4.8000_DP
              vdwparams%radzero(16)=4.9014_DP

           end select pw91damp

        case('BLYP')
           blypdamp: select case (dispersion_index)

           case(1) blypdamp
              vdwparams%dcoeff(1)=3.4942_DP
              vdwparams%c6coeff(1)=3.2260_DP
              vdwparams%c6coeff(6)=27.4628_DP
              vdwparams%c6coeff(7)=19.5741_DP
              vdwparams%c6coeff(8)=11.7910_DP
              vdwparams%c6coeff(16)=161.0455_DP
              vdwparams%radzero(1)=2.7898_DP
              vdwparams%radzero(6)=3.0362_DP
              vdwparams%radzero(7)=3.2664_DP
              vdwparams%radzero(8)=3.6239_DP
              vdwparams%radzero(16)=4.1775_DP

           case(2) blypdamp
              vdwparams%dcoeff(2)=3.9215_DP
              vdwparams%c6coeff(1)=3.2514_DP
              vdwparams%c6coeff(6)=27.4690_DP
              vdwparams%c6coeff(7)=19.5853_DP
              vdwparams%c6coeff(8)=11.7973_DP
              vdwparams%c6coeff(16)=161.0389_DP
              vdwparams%radzero(1)=2.9339_DP
              vdwparams%radzero(6)=3.4803_DP
              vdwparams%radzero(7)=3.6233_DP
              vdwparams%radzero(8)=3.7464_DP
              vdwparams%radzero(16)=4.7968_DP

           case(3) blypdamp

              vdwparams%dcoeff(3)=23.0035_DP
              vdwparams%c6coeff(1)=3.0293_DP
              vdwparams%c6coeff(6)=27.3861_DP
              vdwparams%c6coeff(7)=19.5207_DP
              vdwparams%c6coeff(8)=11.7721_DP
              vdwparams%c6coeff(16)=161.0430_DP
              vdwparams%radzero(1)=2.8570_DP
              vdwparams%radzero(6)=2.9890_DP
              vdwparams%radzero(7)=3.3017_DP
              vdwparams%radzero(8)=3.6800_DP
              vdwparams%radzero(16)=4.0893_DP

           end select blypdamp

        case('XLYP')
           xlypdamp: select case (dispersion_index)

           case(1) xlypdamp
              vdwparams%dcoeff(1)=3.4939_DP
              vdwparams%c6coeff(1)=3.2150_DP
              vdwparams%c6coeff(6)=27.4630_DP
              vdwparams%c6coeff(7)=19.5753_DP
              vdwparams%c6coeff(8)=11.7908_DP
              vdwparams%c6coeff(16)=161.0455_DP
              vdwparams%radzero(1)=2.7993_DP
              vdwparams%radzero(6)=3.0332_DP
              vdwparams%radzero(7)=3.2581_DP
              vdwparams%radzero(8)=3.6244_DP
              vdwparams%radzero(16)=4.1801_DP

           case(2) xlypdamp
              vdwparams%dcoeff(2)=4.0506_DP
              vdwparams%c6coeff(1)=3.4379_DP
              vdwparams%c6coeff(6)=27.5653_DP
              vdwparams%c6coeff(7)=19.6598_DP
              vdwparams%c6coeff(8)=11.8190_DP
              vdwparams%c6coeff(16)=161.0374_DP
              vdwparams%radzero(1)=2.9208_DP
              vdwparams%radzero(6)=3.3614_DP
              vdwparams%radzero(7)=3.5352_DP
              vdwparams%radzero(8)=3.7252_DP
              vdwparams%radzero(16)=4.9034_DP

           case(3) xlypdamp

              vdwparams%dcoeff(3)=23.0034_DP
              vdwparams%c6coeff(1)=3.0238_DP
              vdwparams%c6coeff(6)=27.3853_DP
              vdwparams%c6coeff(7)=19.5206_DP
              vdwparams%c6coeff(8)=11.7719_DP
              vdwparams%c6coeff(16)=161.0439_DP
              vdwparams%radzero(1)=2.8619_DP
              vdwparams%radzero(6)=2.9877_DP
              vdwparams%radzero(7)=3.2994_DP
              vdwparams%radzero(8)=3.6813_DP
              vdwparams%radzero(16)=3.9901_DP

           end select xlypdamp

        case default

           ! qoh: The damping coefficient for the three damping functions
           vdwparams%dcoeff(1)=3.0000_DP
           vdwparams%dcoeff(2)=3.5400_DP
           vdwparams%dcoeff(3)=23.0000_DP

           ! qoh:  Values from Wu and Yang in hartree/bohr ^ 6
           vdwparams%c6coeff(1)=2.8450_DP
           vdwparams%c6coeff(6)=27.3200_DP
           vdwparams%c6coeff(7)=19.4800_DP
           vdwparams%c6coeff(8)=11.7600_DP
           ! qoh: C6 from Halgren
           vdwparams%c6coeff(16)=161.0388_DP

           ! qoh: vdw radii from Elstner in Angstrom
           vdwparams%radzero(1)=3.0000_DP
           vdwparams%radzero(2)=3.0000_DP
           vdwparams%radzero(3)=3.8000_DP
           vdwparams%radzero(4)=3.8000_DP
           vdwparams%radzero(5)=3.8000_DP
           vdwparams%radzero(6)=3.8000_DP
           vdwparams%radzero(7)=3.8000_DP
           vdwparams%radzero(8)=3.8000_DP
           vdwparams%radzero(16)=4.8000_DP

        end select

        ! qoh: Unoptimised parameters from Halgren
        vdwparams%c6coeff(9)=6.2413_DP
        vdwparams%c6coeff(15)=190.5033_DP
        vdwparams%c6coeff(17)=103.5612_DP
        vdwparams%c6coeff(35)=201.8972_DP

        vdwparams%radzero(9)=3.09_DP
        vdwparams%radzero(15)=4.8000_DP
        vdwparams%radzero(17)=4.09_DP
        vdwparams%radzero(35)=4.33_DP

        ! qoh: Array containing the Neff from Wu and Yang
        vdwparams%neff(1)=0.5300_DP
        vdwparams%neff(2)=0.0000_DP
        vdwparams%neff(3)=0.0000_DP
        vdwparams%neff(4)=0.0000_DP
        vdwparams%neff(5)=0.0000_DP
        vdwparams%neff(6)=2.0200_DP
        vdwparams%neff(7)=2.5200_DP
        vdwparams%neff(8)=2.6500_DP
        ! qoh: Neff from Halgren
        vdwparams%neff(9)=3.48_DP
        vdwparams%neff(15)=4.5_DP
        vdwparams%neff(16)=4.8_DP
        vdwparams%neff(17)=5.10_DP
        vdwparams%neff(35)=6.00_DP

    ! mjsp: Populate vdwparams with Grimme's -D2 dispersion values
    case(4)
        select case (xc_func)
        case('PBE')
            vdwparams%s6coeff=0.7500_DP
        case('BLYP')
            vdwparams%s6coeff=1.2000_DP
        case('B3LYP')
            vdwparams%s6coeff=1.050_DP
        case('RPBE')
            vdwparams%s6coeff=1.250_DP
        case('REVPBE')
            vdwparams%s6coeff=1.250_DP
        end select

        vdwparams%dcoeff(4)=20.0000_DP

        ! mjsp: C6 parameter values from Grimme converted in hartree*bohr^6
        vdwparams%c6coeff(1)=2.4283387943_DP
        vdwparams%c6coeff(2)=1.3876221682_DP
        vdwparams%c6coeff(3)=27.9258961350_DP
        vdwparams%c6coeff(4)=27.9258961350_DP
        vdwparams%c6coeff(5)=54.2907173308_DP
        vdwparams%c6coeff(6)=30.3542349293_DP
        vdwparams%c6coeff(7)=21.3346908361_DP
        vdwparams%c6coeff(8)=12.1416939717_DP
        vdwparams%c6coeff(9)=13.0089578269_DP
        vdwparams%c6coeff(10)=10.9275245746_DP
        vdwparams%c6coeff(11)=99.0415322552_DP
        vdwparams%c6coeff(12)=99.0415322552_DP
        vdwparams%c6coeff(13)=187.1555399358_DP
        vdwparams%c6coeff(14)=160.0969076559_DP
        vdwparams%c6coeff(15)=135.9869724834_DP
        vdwparams%c6coeff(16)=96.6131934608_DP
        vdwparams%c6coeff(17)=87.9405549096_DP
        vdwparams%c6coeff(18)=79.9617274424_DP
        vdwparams%c6coeff(19)=187.3289927068_DP
        vdwparams%c6coeff(20)=187.3289927068_DP
        vdwparams%c6coeff(21)=187.3289927068_DP
        vdwparams%c6coeff(22)=187.3289927068_DP
        vdwparams%c6coeff(23)=187.3289927068_DP
        vdwparams%c6coeff(24)=187.3289927068_DP
        vdwparams%c6coeff(25)=187.3289927068_DP
        vdwparams%c6coeff(26)=187.3289927068_DP
        vdwparams%c6coeff(27)=187.3289927068_DP
        vdwparams%c6coeff(28)=187.3289927068_DP
        vdwparams%c6coeff(29)=187.3289927068_DP
        vdwparams%c6coeff(30)=187.3289927068_DP
        vdwparams%c6coeff(31)=294.6962579711_DP
        vdwparams%c6coeff(32)=296.6042384524_DP
        vdwparams%c6coeff(33)=283.9421861676_DP
        vdwparams%c6coeff(34)=219.2443025754_DP
        vdwparams%c6coeff(35)=216.2956054679_DP
        vdwparams%c6coeff(36)=208.3167780008_DP
        vdwparams%c6coeff(37)=427.9079861182_DP
        vdwparams%c6coeff(38)=427.9079861182_DP
        vdwparams%c6coeff(39)=427.9079861182_DP
        vdwparams%c6coeff(40)=427.9079861182_DP
        vdwparams%c6coeff(41)=427.9079861182_DP
        vdwparams%c6coeff(42)=427.9079861182_DP
        vdwparams%c6coeff(43)=427.9079861182_DP
        vdwparams%c6coeff(44)=427.9079861182_DP
        vdwparams%c6coeff(45)=427.9079861182_DP
        vdwparams%c6coeff(46)=427.9079861182_DP
        vdwparams%c6coeff(47)=427.9079861182_DP
        vdwparams%c6coeff(48)=427.9079861182_DP
        vdwparams%c6coeff(49)=647.3257414646_DP
        vdwparams%c6coeff(50)=671.4356766370_DP
        vdwparams%c6coeff(51)=666.7524518194_DP
        vdwparams%c6coeff(52)=550.5390952327_DP
        vdwparams%c6coeff(53)=546.3762287281_DP
        vdwparams%c6coeff(54)=520.1848603034_DP


        ! mjsp: vdw radii from Grimme in Bohr
        vdwparams%radzero(1)=1.8916157150_DP
        vdwparams%radzero(2)=1.9124027009_DP
        vdwparams%radzero(3)=1.5590239409_DP
        vdwparams%radzero(4)=2.6607341925_DP
        vdwparams%radzero(5)=2.8062430937_DP
        vdwparams%radzero(6)=2.7438821360_DP
        vdwparams%radzero(7)=2.6399472066_DP
        vdwparams%radzero(8)=2.5360122772_DP
        vdwparams%radzero(9)=2.4320773478_DP
        vdwparams%radzero(10)=2.3489294043_DP
        vdwparams%radzero(11)=2.1618465314_DP
        vdwparams%radzero(12)=2.5775862490_DP
        vdwparams%radzero(13)=3.0972608960_DP
        vdwparams%radzero(14)=3.2427697971_DP
        vdwparams%radzero(15)=3.2219828112_DP
        vdwparams%radzero(16)=3.1804088395_DP
        vdwparams%radzero(17)=3.0972608960_DP
        vdwparams%radzero(18)=3.0141129525_DP
        vdwparams%radzero(19)=2.8062430937_DP
        vdwparams%radzero(20)=2.7854561078_DP
        vdwparams%radzero(21)=2.9517519948_DP
        vdwparams%radzero(22)=2.9517519948_DP
        vdwparams%radzero(23)=2.9517519948_DP
        vdwparams%radzero(24)=2.9517519948_DP
        vdwparams%radzero(25)=2.9517519948_DP
        vdwparams%radzero(26)=2.9517519948_DP
        vdwparams%radzero(27)=2.9517519948_DP
        vdwparams%radzero(28)=2.9517519948_DP
        vdwparams%radzero(29)=2.9517519948_DP
        vdwparams%radzero(30)=2.9517519948_DP
        vdwparams%radzero(31)=3.1180478818_DP
        vdwparams%radzero(32)=3.2635567830_DP
        vdwparams%radzero(33)=3.3259177406_DP
        vdwparams%radzero(34)=3.3467047265_DP
        vdwparams%radzero(35)=3.3051307548_DP
        vdwparams%radzero(36)=3.2635567830_DP
        vdwparams%radzero(37)=3.0764739101_DP
        vdwparams%radzero(38)=3.0348999383_DP
        vdwparams%radzero(39)=3.0972608960_DP
        vdwparams%radzero(40)=3.0972608960_DP
        vdwparams%radzero(41)=3.0972608960_DP
        vdwparams%radzero(42)=3.0972608960_DP
        vdwparams%radzero(43)=3.0972608960_DP
        vdwparams%radzero(44)=3.0972608960_DP
        vdwparams%radzero(45)=3.0972608960_DP
        vdwparams%radzero(46)=3.0972608960_DP
        vdwparams%radzero(47)=3.0972608960_DP
        vdwparams%radzero(48)=3.0972608960_DP
        vdwparams%radzero(49)=3.1596218536_DP
        vdwparams%radzero(50)=3.4090656842_DP
        vdwparams%radzero(51)=3.5545745853_DP
        vdwparams%radzero(52)=3.5753615712_DP
        vdwparams%radzero(53)=3.5753615712_DP
        vdwparams%radzero(54)=3.5545745853_DP


    end select

    call internal_param_override()

  contains

    subroutine internal_param_override()

    !==================================================================!
    ! This subroutine overrides the parameters set above with those    !
    ! specified in the input file.                                     !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in December 2010.                        !
    !==================================================================!

      use constants, only: DP, VERBOSE, stdout
      use rundat, only: pub_output_detail, pub_vdw_dcoeff

      implicit none

      real(kind=DP) :: c6coeff
      real(kind=DP) :: radzero
      real(kind=DP) :: neff
      integer :: nzatom ! atomic number
      integer :: row ! current row
      integer :: ios ! I/O status flag

      if (allocated(vdwparam_override)) then
         if (pub_on_root) write(stdout, '(a)') "Dispersion Atomic Parameter Override(s):"
         do row=1,size(vdwparam_override)
            ! qoh: Initialise variables
            nzatom=0
            c6coeff=-1.0_DP
            radzero=-1.0_DP
            neff=-1.0_DP

            ! qoh: read a row of the override block
            read (vdwparam_override(row),*,iostat=ios) nzatom, c6coeff, &
                 radzero, neff
            ! qoh: Verify content of the row
            if ((nzatom .gt. 109 .or. nzatom .lt. 1) .or. ios .ne. 0) then
               if (pub_output_detail >= VERBOSE) then
                  if (pub_on_root) write(stdout,*) "WARNING: Invalid line in VDW_PARAMS block &
                       &ignoring. Line is:"
                  if (pub_on_root) write(stdout,*) vdwparam_override(row)
               end if
               cycle
            end if
            if (c6coeff .ge. 0.0_DP) then ! mjsp: override c6coeff
               vdwparams%c6coeff(nzatom) = c6coeff
               if (pub_on_root) write(stdout,'(2x,a, i5,1x,a, 2f15.5)') "Atoms with Z=", nzatom, " overriden with &
                       &c6coeff=", c6coeff
            end if
            if (radzero .gt. epsilon(1.0_DP)) then ! mjsp: override radzero
               vdwparams%radzero(nzatom) = radzero
               if (pub_on_root) write(stdout,'(2x,a, i5,1x,a, 2f15.5)') "Atoms with Z=", nzatom, " overriden with &
                       &radzero=", radzero
            end if
            if (neff .ge. 0.0_DP) then ! mjsp: override neff
               vdwparams%neff(nzatom) = neff
               if (pub_on_root) write(stdout,'(2x,a, i5,1x,a, 2f18.5)') "Atoms with Z=", nzatom, " overriden with &
                       &neff=", neff
            end if
         end do
         if (pub_on_root) write(stdout, '(a)') "Dispersion Atomic Parameter Override(s) Complete"
      end if

      if (pub_vdw_dcoeff .ge. 0.0_DP) then ! mjsp: override dcoeff
         vdwparams%dcoeff(dispersion_index) = pub_vdw_dcoeff
         if (pub_on_root) write(stdout,'(a, 2f16.5)') "Dispersion damping dcoeff parameter override=",&
                 pub_vdw_dcoeff
      end if

    end subroutine internal_param_override

  end subroutine vdwcorrection_initializeparams

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_calculate_energy(elements,cell,par) !input

    !==================================================================!
    ! This subroutine calculates the dispersion correction to the      !
    ! total energy.                                                    !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Modified by Quintin Hill on 29/05/2009 to account for periodic   !
    ! boundary conditions.                                             !
    ! Modified by Max Phipps Feb 2013 to include Grimme's D2 correction!
    ! Modified to remove pub_par by Robert Charlton, 13/06/2018.       !
    ! Parallelized by Arihant Bhandari in July 2022.                   !
    !==================================================================!

    use comms,           only: pub_my_proc_id, comms_reduce, pub_on_root, &
         comms_bcast, pub_root_proc_id
    use constants,       only: DP, stdout, VERBOSE
    use geometry,        only: point, operator(.DOT.), operator(+), &
         operator(*), operator(-)
    use ion,             only: element
    use parallel_strategy, only: PARAL_INFO
    use rundat,          only: pub_output_detail, pub_vdw_bc_is_periodic, &
         pub_vdw_radial_cutoff, pub_print_qc, pub_threads_max
    use simulation_cell, only: CELL_INFO, simulation_cell_num_periodic_images
    use timer, only: timer_clock
    use utils,           only: utils_qc_print

    implicit none

    type(ELEMENT), intent(in)   :: elements(:)
    type(CELL_INFO), intent(in) :: cell
    type(PARAL_INFO), intent(in):: par

    ! Internal variables
    integer       :: i, zi, zj
    integer       :: orig_iat, orig_jat ! atom counters for loops
    integer       :: local_iat, global_iat, global_jat ! atom ids
    integer       :: a1_neighbour,a2_neighbour, a3_neighbour ! Periodic images
    integer       :: periodicity(3) ! number of periodic cells along each direction
    type(POINT)   :: lat(3)       ! simulation cell lattice vectors
    type(POINT)   :: ri, rj
    real(kind=DP) :: cutoff       ! radial cutoff for van der Waals interactions
    real(kind=DP) :: distance     ! Distance between pairs of atoms
    real(kind=DP) :: sqdist       ! Square distance between pairs of atoms
    real(kind=DP) :: c6coeff      ! The c6coefficient of the pair
    real(kind=DP) :: damping      ! The damping for the pair
    real(kind=DP) :: distfromco   ! distance - cutoff
    real(kind=DP) :: smoothco     ! cutoff smoothing
    real(kind=DP) :: energy, de   ! dispersion energy
    type(point)   :: displacement ! Vector between pairs of atoms
    character(len=*), parameter :: myself = 'vdwcorrection_calculate_energy'
    real(kind=DP), parameter :: low_cut = 1D-6 ! ab: For I==J

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    energy = 0.0_DP

    call vdwcorrection_initializeparams

    call vdwcorrection_warnings(elements)

    ! ab: broadcast to non-root procs
    call comms_bcast(pub_root_proc_id, vdwparams%c6coeff(:))
    call comms_bcast(pub_root_proc_id, vdwparams%radzero(:))
    call comms_bcast(pub_root_proc_id, vdwparams%dcoeff(:))
    call comms_bcast(pub_root_proc_id, vdwparams%neff(:))
    call comms_bcast(pub_root_proc_id, vdwparams%s6coeff)

    ! ab: lattice vectors
    lat(1) = cell%a1
    lat(2) = cell%a2
    lat(3) = cell%a3

    ! ab: radial cutoff for van der Waals interactions
    cutoff = pub_vdw_radial_cutoff

    ! ab: count number of periodic cells required along each direction
    call simulation_cell_num_periodic_images(lat, cutoff, periodicity)
    do i=1,3
       if (.not. pub_vdw_bc_is_periodic(i)) then
          periodicity(i) = 0
       end if
    end do

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(STATIC)  &
!$OMP PRIVATE(local_iat,global_iat,orig_iat,ri,zi,displacement,sqdist,distance,&
!$OMP      global_jat,orig_jat,rj,zj,damping,de,                               &
!$OMP      a1_neighbour,a2_neighbour,a3_neighbour,distfromco,smoothco,c6coeff) &
!$OMP SHARED(par,elements,periodicity,lat,pub_my_proc_id,cutoff,vdwparams)     &
!$OMP REDUCTION(+:energy)
    ! ab: loop over local atoms
    loop_atom_I:                                                               &
    do local_iat = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_iat = par%first_atom_on_proc(pub_my_proc_id) + local_iat-1
       orig_iat = par%orig_atom(global_iat)
       ri = elements(orig_iat)%centre
       zi = elements(orig_iat)%atomic_number
       ! ab: loop over atoms J <= I
       loop_atom_J:                                                            &
       do global_jat = 1, global_iat
          orig_jat = par%orig_atom(global_jat)
          rj = elements(orig_jat)%centre
          zj = elements(orig_jat)%atomic_number

          ! qoh: Calculate c6 coefficient
          c6coeff = vdwcorrection_c6(zi, zj)

          ! qoh: Loop over periodic images of atom2 in adjacent cells
          a1: do a1_neighbour = -periodicity(1),periodicity(1)
             a2: do a2_neighbour = -periodicity(2),periodicity(2)
                a3: do a3_neighbour = -periodicity(3),periodicity(3)

                   displacement =  ri - rj &
                        - real(a1_neighbour,kind=DP)*lat(1) &
                        - real(a2_neighbour,kind=DP)*lat(2) &
                        - real(a3_neighbour,kind=DP)*lat(3)

                   sqdist = displacement .DOT. displacement
                   distance = sqrt(sqdist)

                   within_cutoff: if (distance > low_cut .and. distance < cutoff) then

                      distfromco = distance - cutoff
                      smoothco = 1.0_DP - exp(-1.0_DP * distfromco**2)

                      ! qoh : Get damping function
                      damping = vdwcorrection_damping(zi,zj,distance)

                      ! qoh: distance**6 = sqdist**3
                      if (orig_iat == orig_jat) then
                         de = (c6coeff*damping/sqdist**3)*smoothco/2.0_DP
                      else
                         de = (c6coeff*damping/sqdist**3)*smoothco
                      end if

                      energy = energy - de

                   end if within_cutoff

                end do a3
             end do a2
          end do a1

       enddo loop_atom_J
    enddo loop_atom_I
!$OMP END PARALLEL DO

    call comms_reduce('SUM', energy)

    ! msjp: s6coeff=1 for all dispersion functions that are not Grimme's
    pub_dispersion_energy = vdwparams%s6coeff * energy

    ! ab: print information on root only
    if (pub_on_root) then

       write(stdout,'(a, e15.8,1x,a)') &
            'Dispersion Correction Energy: ',pub_dispersion_energy, 'Hartree'

       if(pub_print_qc) then
          call utils_qc_print('Dispersion Energy',pub_dispersion_energy)
       end if

       ! mjsp: end of dispersion printout block
       if (pub_output_detail>=VERBOSE) then
          write(stdout,'(a)') repeat('=',80)
       end if

    end if

    call timer_clock(myself,2)

  end subroutine vdwcorrection_calculate_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_calculate_forces(vdw_forces,elements,cell,par)

    !==================================================================!
    ! This subroutine calculates the dispersion correction to the      !
    ! total forces.                                                    !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Modified by Quintin Hill on 29/05/2009 to account for periodic   !
    ! boundary conditions.                                             !
    ! Modified by Max Phipps Feb 2013 to include Grimme's D2 correction!
    ! Modified to remove pub_par by Joseph Prentice, May 2018          !
    ! PBC bug fixed by Arihant Bhandari, Nov 2020                      !
    ! Parallelized and corrected by Arihant Bhandari in July 2022.     !
    !==================================================================!

    use comms,           only: pub_on_root, comms_reduce, pub_my_proc_id
    use constants,       only: stdout, VERBOSE
    use geometry,        only: point, operator(.DOT.), operator(+), &
         operator(*), operator(-)
    use ion,             only: element
    use parallel_strategy, only: PARAL_INFO
    use rundat,          only: pub_output_detail, pub_vdw_radial_cutoff, &
         pub_vdw_bc_is_periodic, pub_threads_max
    use simulation_cell, only: CELL_INFO, simulation_cell_num_periodic_images
    use timer, only: timer_clock

    implicit none

    ! Arguments

    real(kind=DP), intent(out) :: vdw_forces(:,:)
    type(ELEMENT), intent(in)  :: elements(:)
    type(CELL_INFO), intent(in) :: cell
    type(PARAL_INFO), intent(in):: par

    ! Internal variables
    integer       :: i, zi, zj
    integer       :: orig_iat, orig_jat ! atom counters for loops
    integer       :: local_iat, global_iat, global_jat ! atom ids
    integer       :: a1_neighbour,a2_neighbour, a3_neighbour ! Periodic images
    integer       :: periodicity(3) ! number of periodic cells along each direction
    type(POINT)   :: lat(3)       ! simulation cell lattice vectors
    type(POINT)   :: ri, rj
    real(kind=DP) :: cutoff       ! radial cutoff for van der Waals interactions
    real(kind=DP) :: distance     ! Distance between pairs of atoms
    real(kind=DP) :: sqdist       ! Square distance between pairs of atoms
    real(kind=DP) :: c6coeff      ! The c6coefficient of the pair
    real(kind=DP) :: damping      ! The damping for the pair
    real(kind=DP) :: dampingdrv   ! The damping derivative for the pair
    real(kind=DP) :: drvcommon    ! The common part of the derivative
    real(kind=DP) :: distfromco   ! distance - cutoff
    real(kind=DP) :: smoothco     ! cutoff smoothing
    real(kind=DP) :: smoothcodrv  ! cutoff smoothing derivative
    type(point)   :: displacement ! Vector between pairs of atoms
    character(len=*), parameter :: myself = 'vdwcorrection_calculate_forces'

    ! -------------------------------------------------------------------------

    call timer_clock(myself,1)

    vdw_forces = 0.0_DP

    ! ab: lattice vectors
    lat(1) = cell%a1
    lat(2) = cell%a2
    lat(3) = cell%a3

    ! ab: radial cutoff for van der Waals interactions
    cutoff = pub_vdw_radial_cutoff

    ! ab: count number of periodic cells required along each direction
    call simulation_cell_num_periodic_images(lat, cutoff, periodicity)
    do i=1,3
       if (.not. pub_vdw_bc_is_periodic(i)) then
          periodicity(i) = 0
       end if
    end do

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(NONE) SCHEDULE(STATIC)  &
!$OMP PRIVATE(local_iat,global_iat,orig_iat,ri,zi,displacement,sqdist,distance,&
!$OMP      global_jat,orig_jat,rj,zj,damping,                                  &
!$OMP      smoothcodrv,dampingdrv,drvcommon,c6coeff,                           &
!$OMP      a1_neighbour,a2_neighbour,a3_neighbour,distfromco,smoothco)         &
!$OMP SHARED(par,elements,periodicity,lat,pub_my_proc_id,cutoff,vdwparams)     &
!$OMP REDUCTION(+:vdw_forces)
    ! ab: loop over local atoms
    loop_atom_I:                                                               &
    do local_iat = 1, par%num_atoms_on_proc(pub_my_proc_id)
       global_iat = par%first_atom_on_proc(pub_my_proc_id) + local_iat-1
       orig_iat = par%orig_atom(global_iat)
       ri = elements(orig_iat)%centre
       zi = elements(orig_iat)%atomic_number
       ! ab: loop over J < I, as J == I does not lead to any force
       loop_atom_J:                                                            &
       do global_jat = 1, global_iat - 1
          orig_jat = par%orig_atom(global_jat)
          rj = elements(orig_jat)%centre
          zj = elements(orig_jat)%atomic_number

          ! qoh: Calculate c6 coefficient
          c6coeff = vdwcorrection_c6(zi, zj)

          ! qoh: Loop over periodic images of atom2 in adjacent cells
          a1: do a1_neighbour = -periodicity(1),periodicity(1)
             a2: do a2_neighbour = -periodicity(2),periodicity(2)
                a3: do a3_neighbour = -periodicity(3),periodicity(3)

                   displacement =  ri - rj &
                        - real(a1_neighbour,kind=DP)*lat(1) &
                        - real(a2_neighbour,kind=DP)*lat(2) &
                        - real(a3_neighbour,kind=DP)*lat(3)

                   sqdist = displacement .DOT. displacement
                   distance = sqrt(sqdist)

                   within_cutoff: if (distance < cutoff) then

                      distfromco = distance - cutoff
                      smoothco = 1.0_DP - exp(-1.0_DP * distfromco**2)
                                        ! mjsp: corrected for qoh missing
                                        ! division by distance in smoothcodrv
                      smoothcodrv = 2.0_DP*distfromco*(1.0_DP - smoothco)/distance

                      ! qoh : Get damping function
                      damping = vdwcorrection_damping(zi,zj,distance)

                      dampingdrv = vdwcorrection_drvdamping(zi,zj,distance)

                      ! qoh: distance**6 = sqdist**3
                      ! mjsp: s6coeff=1 for all dispersion functions that are not Grimme's
                      drvcommon = vdwparams%s6coeff*c6coeff/sqdist**3* &
                           ( smoothco * (dampingdrv-6.0_DP*damping/sqdist) &
                           + smoothcodrv*damping )

                      ! cks: This is the *negative* of the derivative
                      ! cks: of the dispersion energy w.r.t. atomic
                      ! cks: coordinates.
                      vdw_forces(1,orig_iat) = vdw_forces(1,orig_iat) &
                           +drvcommon*displacement%x
                      vdw_forces(2,orig_iat) = vdw_forces(2,orig_iat) &
                           +drvcommon*displacement%y
                      vdw_forces(3,orig_iat) = vdw_forces(3,orig_iat) &
                           +drvcommon*displacement%z

                      ! ab: equal and opposite force on atom J.
                      vdw_forces(1,orig_jat) = vdw_forces(1,orig_jat) &
                           -drvcommon*displacement%x
                      vdw_forces(2,orig_jat) = vdw_forces(2,orig_jat) &
                           -drvcommon*displacement%y
                      vdw_forces(3,orig_jat) = vdw_forces(3,orig_jat) &
                           -drvcommon*displacement%z

                   end if within_cutoff
                end do a3
             end do a2
          end do a1
       enddo loop_atom_J
    enddo loop_atom_I
!$OMP END PARALLEL DO

    call comms_reduce('SUM', vdw_forces)

    if (pub_on_root) then
       ! mjsp: end of dispersion printout block
       if (pub_output_detail>=VERBOSE) then
          write(stdout,'(a)') repeat('=',80)
       end if
    end if

    call timer_clock(myself,2)

  end subroutine vdwcorrection_calculate_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_warnings(elements)

    !==================================================================!
    ! This subroutine warns about the use of unoptimised or unavailable!
    ! dispersion parameters.                                           !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 13/02/2009.                           !
    ! Quintin Hill added functional check on 23/02/2009.               !
    ! Modified by Max Phipps for Grimme's correction 18/02/2013        !
    ! Modified to remove pub_par by Robert Charlton, 13/06/2018.       !
    !==================================================================!

    use comms,           only: pub_on_root
    use constants,       only: periodic_table_name, stdout, VERBOSE
    use ion,             only: element
    use rundat,          only: pub_xc_functional, pub_output_detail

    implicit none

    type(element), intent(in)    :: elements(:)

    logical, save                :: already_warned = .false.
    integer                      :: atcnum ! atomic number


    ! mjsp: qoh implemented dispersion functionality warnings
    ! qoh: Atomic numbers of elements with unoptimised dispersion parameters:
    integer, parameter           :: qoh_unoptimised(4) = (/9,15,17,35/)
    ! qoh: Atomic numbers of elements with optimised dispersion parameters:
    integer, parameter           :: qoh_optimised(5) = (/1,6,7,8,16/)

    ! qoh: XC functionals with optimised parameters:
    character(len=6), parameter  :: qoh_xcfoptimised(6) &
         = (/'BLYP  ','PBE   ','PW91  ','REVPBE','RPBE  ','XLYP  '/)


    ! mjsp: Grimme dispersion warnings' arrays follow
    ! mjsp: XC functionals with parameters:
    character(len=6), parameter  :: grimmeD2_xcfoptimised(5) &
         = (/'BLYP  ','PBE   ','B3LYP ','REVPBE', 'RPBE  '/)
    ! mjsp: Atomic numbers of elements with Grimme optimised dispersion parameters (Hydrogen to Xenon):
    integer                      :: grimmeD2_optimised(54)
    integer                      :: i
    do i = 1,54
        grimmeD2_optimised(i) = i
    end do



    if ( pub_output_detail /= VERBOSE .or. already_warned ) return
    if (pub_on_root) then
       select case (dispersion_index)
       case(1:3)

           ! qoh: Loop over all elements and check if each present in this calculation
           do atcnum=1,109
              elpresent: if (any(elements%atomic_number == atcnum) &
                   .and. any(qoh_xcfoptimised == pub_xc_functional)) then
                 if (any(qoh_unoptimised == atcnum)) then
                    write(stdout,'(a,a2)') 'WARNING: Unoptimised dispersion &
                         &parameters used for ', periodic_table_name(atcnum)
                 elseif (.not. any(qoh_optimised == atcnum)) then
                    write(stdout,'(a,a2)') 'WARNING: No dispersion parameters &
                         &available for ', periodic_table_name(atcnum)
                 end if
              end if elpresent
           end do

           if (.not. any(qoh_xcfoptimised == pub_xc_functional)) &
                write(stdout,'(2a)') 'WARNING: No optimised dispersion parameters &
                &available for ', trim(pub_xc_functional)

       ! mjsp: Grimme dispersion functionality warnings
       case (4)

           do atcnum=1,109
               if (any(elements%atomic_number == atcnum) &
                   .and. any(grimmeD2_xcfoptimised == pub_xc_functional) &
                   .and. (.not. any(grimmeD2_optimised == atcnum))) then
                    write(stdout,'(a,a2)') 'WARNING: No dispersion parameters &
                         &available for ', periodic_table_name(atcnum)
               end if
           end do

           if (.not. any(grimmeD2_xcfoptimised == pub_xc_functional)) &
                write(stdout,'(2a)') 'WARNING: No optimised dispersion parameters &
                &available for ', trim(pub_xc_functional)

           ! mjsp/lgv: 01/09/2015: add warning if using RPBE or revPBE functionals:
           if (pub_xc_functional == 'RPBE  ') &
              write(stdout,'(2a)') 'WARNING: s6 parameter for RPBE = 1.25 &
                & (see J. Chem. Phys. 132, 154105, (2010) and J. Chem. Phys. 137, 114105 (2012)).'
           if (pub_xc_functional == 'REVPBE') &
              write(stdout,'(2a)') 'WARNING: s6 parameter for revPBE = 1.25 &
                & (see J. Chem. Phys. 132, 154105, (2010)).'

       end select
    end if

    already_warned = .true.

  end subroutine vdwcorrection_warnings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_override_alloc(override_block,nrows)

    !==================================================================!
    ! This subroutine allocates an array to hold parameter overides.   !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in December 2010.                        !
    !==================================================================!

    use utils, only: utils_alloc_check

    implicit none

    integer, intent(in) :: nrows ! Number of rows of override data
    character(len=80), intent(in) :: override_block(nrows) ! Override data
    integer             :: ierr ! error flag

    allocate(vdwparam_override(nrows), stat=ierr)
    call utils_alloc_check('vdwcorrection_overide_alloc','vdwparam_override',&
         ierr)

    vdwparam_override = override_block

  end subroutine vdwcorrection_override_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_override_dealloc()

    !==================================================================!
    ! This subroutine deallocates the array to hold parameter overides.!
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in December 2010.                        !
    !==================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    integer             :: ierr ! error flag

    if (allocated(vdwparam_override)) then
       deallocate(vdwparam_override, stat=ierr)
       call utils_dealloc_check('vdwcorrection_overide_alloc',&
            'vdwparam_override',ierr)
    end if

  end subroutine vdwcorrection_override_dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function vdwcorrection_c6(nzatom1,nzatom2)

    !==================================================================!
    ! This function calculates a heteroatomic C_6 coefficient from     !
    ! homoatomic C_6 coefficients using the formula given in Elstner's !
    ! paper (J. Chem. Phys. 114(12), 5149-5155). (Corrected error in   !
    ! the formula appearing in this paper.) or by the formula given in !
    ! Grimme's paper (J Comput Chem. 27(15), 1787-99, 2006)            !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Modified by Max Phipps Feb. 2013 for Grimme's C_6 formula        !
    !==================================================================!

    implicit none

    real(kind=DP)       :: vdwcorrection_c6
    integer, intent(in) :: nzatom1 ! Atomic number of atom 1
    integer, intent(in) :: nzatom2 ! Atomic number of atom 2
    real(kind=DP)       :: c61     ! c6 coefficient of atom 1
    real(kind=DP)       :: c62     ! c6 coefficient of atom 2
    real(kind=DP)       :: ne1     ! Effective number of electrons for atom 1
    real(kind=DP)       :: ne2     ! Effective number of electrons for atom 2
    real(kind=DP), parameter :: third = 1.0_DP/3.0_DP

    ! qoh: Set up shorthands

    c61 = vdwparams%c6coeff(nzatom1)
    c62 = vdwparams%c6coeff(nzatom2)
    ne1 = vdwparams%neff(nzatom1)
    ne2 = vdwparams%neff(nzatom2)

    select case (dispersion_index)

    ! mjsp: Elstner's C_6 coefficient formula
    case(1:3)
        if (c61 .lt. epsilon(1.0_DP) .or. c62 .lt.epsilon(1.0_DP) .or. &
             ne1 .lt. epsilon(1.0_DP) .or. ne2 .lt. epsilon(1.0_DP)) then
           vdwcorrection_c6=0.0_DP
        else
           vdwcorrection_c6 = 2.0_DP * (c61**2 * c62**2 * ne1 * ne2)**(third)/&
                (((c61 * ne2**2)**(third)) + ((c62* ne1**2)**(third)))
        end if

    ! mjsp: Grimme's C_6 coefficient formula
    case(4)
        vdwcorrection_c6 = sqrt(c61 * c62)

    end select

  end function vdwcorrection_c6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function vdwcorrection_damping(nzatom1,nzatom2,separation)

    !==================================================================!
    ! This function calculates the damping function specified by       !
    ! dispersion_index:                                                  !
    ! (1) Damping function from Elstner                                !
    !     (J. Chem. Phys. 114(12), 5149-5155).                         !
    ! (2) First damping function from Wu and Yang (I)                  !
    !     (J. Chem. Phys. 116(2), 515-524).                            !
    ! (3) Second damping function from Wu and Yang (II)                !
    !     (J. Chem. Phys. 116(2), 515-524).                            !
    ! (4) D2 damping function from Grimme                              !
    !     (J Comput Chem. 27(15), 1787-99, 2006)                       !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Merged into a single function in July 2008                       !
    ! Modified by Max Phipps Feb. 2013 for Grimme's -D2 correction     !
    !==================================================================!

    implicit none

    real(kind=DP)             :: vdwcorrection_damping
    integer, intent(in)       :: nzatom1 ! Atomic number of atom 1
    integer, intent(in)       :: nzatom2 ! Atomic number of atom 2
    real(kind=DP), intent(in) :: separation
    real(kind=DP)             :: radzero
    real(kind=DP)             :: expo ! Exponent
    integer, parameter        :: mexpo(2) = (/4,2/)
    integer, parameter        :: nexpo(2) = (/7,3/)
    integer                   :: mexp
    integer                   :: nexp

    radzero = vdwcorrection_radzero(nzatom1,nzatom2)

    select case (dispersion_index)
    case(1,2)

       mexp = mexpo(dispersion_index)
       nexp = nexpo(dispersion_index)

       expo = -vdwparams%dcoeff(dispersion_index)*(separation/radzero)**nexp
       vdwcorrection_damping = ( 1.0_DP - exp(expo))**mexp

    case(3,4)

       expo = -vdwparams%dcoeff(dispersion_index)*((separation/radzero)-1.0_DP)
       vdwcorrection_damping = 1.0_DP/( 1.0_DP + exp(expo))

    case default
       vdwcorrection_damping = 1.0_DP
    end select

  end function vdwcorrection_damping


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function vdwcorrection_drvdamping(nzatom1,nzatom2,separation)

    !==================================================================!
    ! This function calculates the derivative with respect to atomic   !
    !cooridinates of the damping function specified by dispersion_index: !
    ! (1) Damping funtion from Elstner                                 !
    !     (J. Chem. Phys. 114(12), 5149-5155).                         !
    ! (2) First damping function from Wu and Yang (I)                  !
    !     (J. Chem. Phys. 116(2), 515-524).                            !
    ! (3) Second damping function from Wu and Yang (II)                !
    !     (J. Chem. Phys. 116(2), 515-524).                            !
    ! (4) D2 damping function from Grimme                              !
    !     (J Comput Chem. 27(15), 1787-99, 2006)                       !
    ! Note: For simplicity the (r_{A,i}-r_{B_i}) (where i is x,y or z) !
    !       term is omitted here and included in the calling routine.  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in July 2008.                            !
    ! Modified for Grimme D2 correction May 2013                       !
    !==================================================================!

    implicit none

    real(kind=DP)             :: vdwcorrection_drvdamping
    integer,       intent(in) :: nzatom1 ! Atomic number of atom 1
    integer,       intent(in) :: nzatom2 ! Atomic number of atom 2
    real(kind=DP), intent(in) :: separation
    real(kind=DP)             :: radzero
    real(kind=DP)             :: expo ! Exponent
    integer, parameter        :: mexpo(2) = (/4,2/)
    integer, parameter        :: nexpo(2) = (/7,3/)
    integer                   :: mexp
    integer                   :: nexp

    radzero = vdwcorrection_radzero(nzatom1,nzatom2)

    select case (dispersion_index)
    case(1,2)

       mexp = mexpo(dispersion_index)
       nexp = nexpo(dispersion_index)

       expo = -vdwparams%dcoeff(dispersion_index)*(separation/radzero)**nexp

       vdwcorrection_drvdamping = real((mexp*nexp),kind=DP)*exp(expo)* &
            vdwparams%dcoeff(dispersion_index)*( 1.0_DP - exp(expo))**(mexp-1)*&
            separation**(nexp-2)/radzero**nexp

    case(3,4)

       expo = -vdwparams%dcoeff(dispersion_index)*((separation/radzero)-1.0_DP)

       vdwcorrection_drvdamping = vdwparams%dcoeff(dispersion_index)*exp(expo)/&
            (separation*radzero*( 1.0_DP + exp(expo))**2)

    case default
       vdwcorrection_drvdamping = 0.0_DP
    end select

  end function vdwcorrection_drvdamping

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function vdwcorrection_radzero(nzatom1,nzatom2)

    !==================================================================!
    ! Function to calculate the R_0 for an atom pair. Uses expression  !
    ! found in Elstner's paper (J. Chem. Phys. 114(12), 5149-5155).    !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Modified by Max Phipps Feb. 2013 for Grimme's correction         !
    !==================================================================!

    use constants, only: DP, ANGSTROM
    implicit none

    real(kind=DP) :: vdwcorrection_radzero
    integer, intent(in) :: nzatom1 ! Atomic number of atom 1
    integer, intent(in) :: nzatom2 ! Atomic number of atom 2

    select case (dispersion_index)

    ! mjsp: qoh implemented Elstner's R_0 for an atomic pair
    case(1:3)

        vdwcorrection_radzero = ANGSTROM*(vdwparams%radzero(nzatom1)**3 + &
             vdwparams%radzero(nzatom2)**3)/&
             (vdwparams%radzero(nzatom1)**2 + vdwparams%radzero(nzatom2)**2)

    ! mjsp: Grimme's -D2 correction R_0 for an atomic pair
    case(4)

        vdwcorrection_radzero = (vdwparams%radzero(nzatom1) + &
                                            vdwparams%radzero(nzatom2))
    end select


  end function vdwcorrection_radzero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module vdwcorrection
