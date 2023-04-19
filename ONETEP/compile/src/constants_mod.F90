! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by Chris-Kriton Skylaris
!
!   February 2000
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module constants

  integer, parameter :: SP = kind(1.0)                 ! Single precision real type
  integer, parameter :: DP = kind(1.0d0)               ! Double precision real type
  integer, parameter :: LONG = selected_int_kind(18)   ! Long integer type
  integer, parameter :: I4B = selected_int_kind(9)     ! 4-byte integer kind

  character, parameter :: CRLF = ACHAR(10)

  real(kind=DP), parameter :: PI=3.141592653589793238462643383279502884197_DP
  real(kind=DP), parameter :: SQRT_PI=1.772453850905516027298167483341145182798_DP
  real(kind=DP), parameter :: TWO_PI=2.0_DP * pi
  real(kind=DP), parameter :: FOUR_PI=4.0_DP * pi
  real(kind=DP), parameter :: SQRT_TWO=1.414213562373095048801688724210_DP
  real(kind=DP), parameter :: ONE_THIRD = 1.0_DP/3.0_DP
  real(kind=DP), parameter :: TWO_THIRDS = 2.0_DP/3.0_DP
  real(kind=DP), parameter :: ONE_SIXTH = 1.0_DP/6.0_DP
  real(kind=DP), parameter :: ONE_FIFTH = 1.0_DP/5.0_DP

  real(kind=DP), parameter :: SAFE_DIV_EPS = 1D-30
  ! jd: A constant to compare 'abs(x)' against if you plan to divide by x.
  !     It's not enough to compare against zero or TINY(1.0_DP) to protect
  !     yourself from overflow. The value for SAFE_DIV_EPS is not foolproof
  !     either, but should be sufficient for our purposes.
  !     cf. http://wiki.seas.harvard.edu/geos-chem/index.php/Floating_point_math_issues

#ifndef OLD_CONSTANTS

  ! ndmh 11/02/2011:
  ! For consistency with CASTEP pseudopotentials, the following conversion
  ! factors have been obtained by calling CASTEP's internal conversion routines
  ! as follows:
  !   ANGSTROM = io_atomic_to_unit(1.0_dp,'ang')
  !   HARTREE_IN_EVS = io_atomic_to_unit(1.0_dp,'eV')
  ! based on the CODATA 2002 values (CASTEP's default choice). This produced
  !   electron_mass_si*speed_light_si*fine_structure_si/hbar_si*1e-10_dp
  !   elementary_charge_si/(fine_structure_si**2*electron_mass_si*speed_light_si**2)
  ! where these were derived from the following definitions (CODATA2002):
  !real(kind=dp), parameter, public :: speed_light_si = 299792458.0_dp
  !real(kind=dp), parameter, public :: planck_si = 6.6260693e-34_dp
  !real(kind=dp), parameter, public :: elementary_charge_si = 1.60217653e-19_dp
  !real(kind=dp), parameter, public :: electron_mass_si = 9.1093826e-31_dp
  !real(kind=dp), parameter, public :: avogadro_si = 6.0221415e23_dp
  !real(kind=dp), parameter, public :: pi=3.141592653589793238462643383279502884197_dp
  !real(kind=dp), parameter, public :: two_pi=2.0_dp*pi
  !real(kind=dp), parameter, public :: hbar_si = planck_si/two_pi
  !real(kind=dp), parameter, public :: mu_0_si = 4.0_dp*pi*1e-7_dp
  !real(kind=dp), parameter, public :: epsilon_0_si = 1.0_dp/(mu_0_si*speed_light_si**2)
  !real(kind=dp), parameter, public :: fine_structure_si = elementary_charge_si**2
  !                                 / (4.0_dp*pi*epsilon_0_si*hbar_si*speed_light_si)


  ! The value of one Angstrom in terms of Bohr radii
  real(kind=DP), parameter :: ANGSTROM=1.889726134583548707935_DP

  ! cks: The value of one Hartree in terms of electron-volts
  real(kind=DP), parameter :: HARTREE_IN_EVS=27.2113846081672_DP

#else

  ! ndmh: equivalent values in terms of the old versions of the CASTEP constants
  real(kind=DP), parameter :: ANGSTROM=1.889726313_DP
  real(kind=DP), parameter :: HARTREE_IN_EVS=27.2116529_DP

#endif

  ! tjz07: Additional constants:
  ! fine structure constant. Also, value of c in Hartree atomic units
  real(kind=DP), parameter :: FINE_STRUCTURE=1.0_DP/137.035999074_DP
  ! vv: The value of the speed of light in terms of m/s
  real(kind=DP), parameter :: SPEED_OF_LIGHT_SI=299792458.0_dp
  ! tjz07: conversion factor from Hartrees to nano-seconds
  real(kind=DP), parameter :: HARTREE_IN_NS=2.418884326505_DP*1D-8
  ! vv: Boltzmann constant in terms of Ha/K
  real(kind=DP), parameter :: k_B= 3.1668115744561575d-06
  real(kind=DP), parameter :: HARTREE_IN_DEBYE=2.5417462310548_DP

  ! jd: Avogadro number, CODATA 2018 cf.
  !     https://physics.nist.gov/cgi-bin/cuu/Value?na
  real(kind=DP), parameter :: AVOGADRO = 6.02214076D23
  ! jd: conversion factor from mol/L to atomic units of concentration
  real(kind=DP), parameter :: CONC_MOL_PER_L_TO_AU = &
       (1.0_DP / ANGSTROM)**3 * 1.0D-27 * AVOGADRO
  ! jd: conversion factor from atomic units of concentration to mol/L
  real(kind=DP), parameter :: CONC_AU_TO_MOL_PER_L=1D0 / CONC_MOL_PER_L_TO_AU
  ! jd: thermodynamic standard reference (particles/a_0^3)
  real(kind=DP), parameter :: CONC_0 = CONC_MOL_PER_L_TO_AU
  ! jd: conversion fractor from Ha to kcal/mol
  real(kind=DP), parameter :: HARTREE_TO_KCAL_PER_MOL = 627.50947406
  ! jd: converstion factor from atomic units of electric field (Ha/(e a0))
  !     to MV/cm. 1 V/Angstrom = 100 MV/cm.
  real(kind=DP), parameter :: HARTREE_PER_BOHR_TO_MV_PER_CM = &
       HARTREE_IN_EVS * ANGSTROM * 1D2
  ! ab: electrochemical potential of standard Hydrogen electrode (SHE)
  real(kind=DP), parameter :: mu_ref_SHE = -4.44_DP/HARTREE_IN_EVS
  ! ab: elementary charge (e), CODATA 2018
  real(kind=DP), parameter :: e_SI=1.602176634*1.0D-19
  ! ab: atmospheric pressure
  real(kind=DP), parameter :: ATM_SI=1.01325D5

  ! pdh: square root of minus one
  complex(kind=DP), parameter :: cmplx_i = (0.0_DP,1.0_DP)

  ! ndmh: complex 0
  complex(kind=DP), parameter :: cmplx_0 = (0.0_DP,0.0_DP)
  complex(kind=DP), parameter :: cmplx_1 = (1.0_DP,0.0_DP)

  ! pdh: units for standard output and standard error
  integer :: stdout = 6
  integer :: stderr = 6!stdout

  !  vv: maximum length of file names
  integer, parameter, public :: file_maxsize = 256
  ! JCW: maximum file name length which is conventionally hard-coded
  ! JCW: (use this to warn users of overlong filepaths)
  integer, parameter, public :: file_maxsize_legacy = 80
  ! JCW: maximum file name length which does not cause ONETEP to abort
  ! JCW: in get_rundat, determined through testing. This is shorter than
  ! JCW: file_maxsize_legacy, since esdf_string can only accept label:value
  ! JCW: combinations which are <= llength (80) characters, and some
  ! JCW: keywords use "pub_rootname" as the default value.
  ! JCW: Any label:value combination passed to esdf_string which includes
  ! JCW: a filename will need to have combined (after esdf_reduce is applied
  ! JCW: to the label) length of 80 characters or less, including the filename
  ! JCW: and colon.
  integer, parameter, public :: file_maxsize_rundat = 64

  ! ndmh: sizes of datatypes, for size estimation
  integer, parameter :: logical_size = 1
  integer, parameter :: char_size = 1
  integer, parameter :: int_size = 4
  integer, parameter :: real_size = 8
  integer, parameter :: cmplx_size = 16

  ! cks: output level integer constants
  integer, parameter :: BRIEF   = 0
  integer, parameter :: NORMAL  = 1
  integer, parameter :: VERBOSE = 2
  integer, parameter :: PROLIX = 3
  integer, parameter :: MAXIMUM = 10

  ! pdh: spin polarisation
  integer, parameter :: max_spins = 2
  integer, parameter :: UP = 1
  integer, parameter :: DN = 2

  ! ndmh: PAW energy component labels
  integer, parameter, public :: paw_en_size     = 9
  integer, parameter, public :: paw_en_dijhat   = 1
  integer, parameter, public :: paw_en_dij0     = 2
  integer, parameter, public :: paw_en_ehart    = 3
  integer, parameter, public :: paw_en_exc      = 4
  integer, parameter, public :: paw_en_exc_dc   = 5
  integer, parameter, public :: paw_en_etxc     = 6
  integer, parameter, public :: paw_en_etxc_dc  = 7
  integer, parameter, public :: paw_en_dijxc    = 8
  integer, parameter, public :: paw_en_exc_core = 9

  ! aam: The symbols of the elements in the periodic table
  character(len=2), parameter, dimension(109) :: periodic_table_name= (/ &
       & 'H ',                                                                                'He', &
       & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
       & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
       & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
       & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
       & 'Cs','Ba', &
       & 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', &
       & 'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
       & 'Fr','Ra', &
       & 'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr', &
       & 'Rf','Db','Sg','Bh','Hs','Mt' /)

  ! aam: The atomic masses of the elements in the periodic table
  real(kind=dp), parameter, dimension(109) :: periodic_table_mass = (/ &
     & 1.00794_dp,                                                                                                    4.00260_dp, &
     & 6.941_dp, 9.012187_dp,                            10.811_dp, 12.0107_dp, 14.00674_dp, 15.9994_dp, 18.99840_dp, 20.1797_dp, &
     & 22.98977_dp, 24.3050_dp,                           26.98154_dp, 28.0855_dp, 30.97376_dp, 32.066_dp, 35.4527_dp, 39.948_dp, &
     & 39.0983_dp, 40.078_dp, &
     &   44.95591_dp, 47.867_dp, 50.9415_dp, 51.9961_dp, 54.93805_dp, 55.845_dp, 58.93320_dp, 58.6934_dp, 63.546_dp, 65.39_dp, &
     &                                                           69.723_dp, 72.61_dp, 74.92160_dp, 78.96_dp, 79.904_dp, 83.80_dp, &
     & 85.4678_dp, 87.62_dp, &
     &   88.90585_dp, 91.224_dp, 92.90638_dp, 95.94_dp, 98.0_dp, 101.07_dp, 102.90550_dp, 106.42_dp, 107.8682_dp, 112.411_dp, &
     &                                                    114.818_dp, 118.710_dp, 121.760_dp, 127.60_dp, 126.90447_dp, 131.29_dp, &
     & 132.90545_dp, 137.327_dp, &
     &   138.9055_dp, 140.116_dp, 140.90765_dp, 144.24_dp, 145.0_dp, 150.36_dp, 151.964_dp, 157.25_dp, 158.92534_dp, &
     &                                    162.50_dp, 164.93032_dp, 167.26_dp, 168.93421_dp, 173.04_dp, 174.967_dp, &
     &   178.49_dp, 180.9479_dp, 183.84_dp, 186.207_dp, 190.23_dp, 192.217_dp, 195.078_dp, 196.96655_dp, 200.59_dp, &
     &                                                      204.3833_dp, 207.2_dp, 208.98038_dp, 209.0_dp, 210.0_dp, 222.0_dp, &
     & 223.0_dp, 226.0_dp, &
     &   227.0_dp, 232.0381_dp, 231.03588_dp, 238.0289_dp, 237.0_dp, 244.0_dp, 243.0_dp, 247.0_dp, 247.0_dp, 251.0_dp, 252.0_dp, &
     &                                                                           257.0_dp, 258.0_dp, 259.0_dp, 262.0_dp, &
     &   261.0_dp, 262.0_dp, 263.0_dp, 264.0_dp, 265.0_dp, 268.0_dp /)

  ! gab: Equilibrium Van Der Waals radii of elements following the
  ! Alvarez model[1]. Only filled for the first 4 rows, I, Pt and Pd.
  !
  ! 1. Alvarez, S. A cartography of the van der Waals territories.
  ! Dalt. Trans. 42, 8617 (2013).

  real(kind=dp), parameter, dimension(109) :: periodic_table_ss_radii = (/ &
     & 2.26_dp, 2.70_dp,                                              & ! H-He
     & 4.01_dp, 3.74_dp, 3.61_dp, 3.34_dp, 3.14_dp, 2.83_dp, 2.76_dp, & ! Li-Ne
     & 2.99_dp,                                                       & ! Li-Ne
     & 4.72_dp, 4.74_dp, 4.25_dp, 4.14_dp, 3.59_dp, 3.57_dp, 3.44_dp, & ! Na-Ar
     & 3.67_dp,                                                       & ! Na-Ar
     & 5.16_dp, 4.95_dp, 4.88_dp, 4.65_dp, 4.57_dp, 4.63_dp, 4.63_dp, & ! Ca-Kr
     & 4.61_dp, 4.54_dp, 4.54_dp, 4.50_dp, 4.52_dp, 4.38_dp, 4.33_dp, & ! Ca-Kr
     & 3.55_dp, 3.44_dp, 3.51_dp, 3.91_dp, 85.1_dp, 87.1_dp,          & ! Ca-Kr
     ! gab1u17: Values from this point are untested - seperate parametrisation should be used.
     & 6.06_dp, 5.36_dp, 5.19_dp, 4.76_dp, 4.84_dp,                   & ! Rb-Xe
     & 4.63_dp, 4.61_dp, 4.64_dp, 4.61_dp, 4.06_dp,                   & ! Rb-Xe
     & 4.78_dp, 4.70_dp, 4.59_dp, 4.57_dp, 4.67_dp,                   & ! Rb-Xe
     & 3.89_dp,                                                       & ! Rb-Xe
     & 6.57_dp, 5.72_dp, 5.63_dp, 5.44_dp,                            & ! Cs-Rn
     & 5.52_dp, 5.57_dp, 5.48_dp, 5.42_dp, 5.34_dp,                   & ! Cs-Rn
     & 5.27_dp, 5.42_dp, 5.30_dp, 5.35_dp, 5.27_dp,                   & ! Cs-Rn
     & 5.29_dp, 5.18_dp, 4.97_dp, 4.78_dp, 4.85_dp,                   & ! Cs-Rn
     & 4.70_dp, 4.68_dp, 4.55_dp, 4.33_dp, 4.38_dp,                   & ! Cs-Rn
     & 4.63_dp, 4.67_dp, 4.91_dp, 4.80_dp, 5.29_dp,                   & ! Cs-Rn
     & 5.53_dp, 5.44_dp, 5.12_dp,                                     & ! Cs-Rn
     & 5.31_dp, 5.35_dp, 5.76_dp, 6.42_dp, 5.76_dp, 5.10_dp,          & ! Cs-Rn
     ! gab1u17: Species from Fr onwards aren't calculated Alverez's model.
     ! These elements are probably beyond the scope of the ISM at this point
     ! anyway.
     & 101.0_dp,101.0_dp,101.0_dp,101.0_dp,101.0_dp,                  & ! Cs-Rn
     & 101.0_dp,101.0_dp,101.0_dp,101.0_dp,101.0_dp,101.0_dp,         & ! Cs-Rn
     & 101.0_dp,101.0_dp,101.0_dp,101.0_dp,101.0_dp,101.0_dp /)         ! Cs-Rn

  ! jd: Size of real(kind=DP). AFAIK only Fortran 2008 has a portable
  !     sort-of-sizeof. For now let's assume we all run on x86_64.
  !     This is only used for memory estimates.
  integer, parameter :: SIZEOF_DOUBLE = 8
  integer, parameter :: SIZEOF_INT = 4

  ! jd: Garbage values for initialisations
  real(kind=DP), parameter :: garbage_real = -9D99
  complex(kind=DP), parameter :: garbage_complex = (-9D99,-9D99)
  integer, parameter       :: garbage_int = -42424242

  ! jd: Generic large double
  real(kind=DP), parameter :: very_big_double = 9D99

  ! jd: Constants for describing coefficient types and metric matrix types in
  !     SW RI and SW expansion. Do not reorder, there are do loops using these.
  integer, parameter :: SW_V = 1
  integer, parameter :: SW_O = 2
  integer, parameter :: SW_W = 3
  integer, parameter :: SW_D = 4 ! mask added to above when using D-coeffs
  character, parameter, dimension(8) :: SW_LETTERS = (/'V','O','W','?','v','o','w','?'/)

  ! jd: Constants for describing the metric in HFx
  integer, parameter :: METRIC_ELECTROSTATIC = 1
  integer, parameter :: METRIC_OVERLAP = 2

  ! jd: Maximum number of SWRI block entries in rundat. Logically this belongs
  !     in sw_resolution_of_identity_mod, but this value is also used to
  !     dimension a member of ELEMENT in ion. Storing it here avoids a
  !     dependency nightmare.
  integer, parameter :: max_rundat_swri_block_entries = 1

  ! jd: Length of words returned by utils_nth_word_of
  integer, parameter :: MAX_WORD_LENGTH = 512

  ! jd: A threshold for kernel cutoff (in bohr) below which a warning is issued
  real(kind=DP), parameter :: kernel_cutoff_sanity_threshold = 20.0_DP

  ! mjsp: Constants for identifying EDA modes and handling
  ! EDA stage ordering (for restarts).
  integer, parameter :: EDA_INIT = 0             ! Initialisation mode
  integer, parameter :: EDA_ISOLATED = 1         ! Isolated fragment calculations
  integer, parameter :: EDA_FROZENNONIDEM = 2    ! 'Pure' frozen state
  integer, parameter :: EDA_FROZENIDEM = 3       ! Normalised frozen state
  integer, parameter :: EDA_POLFRAGLOC_DEVEL = 4 ! DEVEL: iterate the fragments,
                                     ! polarising each fragment in the field
                                     ! of the surrounding frozen fragments.
  integer, parameter :: EDA_POLFRAGLOC_DEVEL_OFF = 5
  integer, parameter :: EDA_POLSIMUL = 6         ! Stoll et. al. SCF MI equations
  integer, parameter :: EDA_POLSIMUL_OFF = 7     ! Stoll et. al. SCF MI equations temporary off switch
  integer, parameter :: EDA_CTFRAGLOC_DEVEL = 8 ! DEVEL: iterate the fragments,
                                     ! delocalising between fragment pairs in the field
                                     ! of the surrounding fragments in their
                                     ! polarised states.
  integer, parameter :: EDA_CTFRAGLOC_DEVEL_OFF = 9
  integer, parameter :: EDA_FULL = 10            ! Fully optimized supermolecule
  integer, parameter :: EDA_SUPER_PREP = 11      ! Tool to join fragment data to prepare
                                                 ! the supermolecule data

  ! jd: SWx, DMA and polemb related constants
  integer, parameter :: REP_N_SWEXES = 8
  integer, parameter :: REP_SWEX_HFX = 1         ! val-val or cond-cond
  integer, parameter :: REP_SWEX_HFX_OTHER = 2   ! cond-val and joint-val. val-cond is not used.
  integer, parameter :: REP_SWEX_PROPERTIES_DMA_1 = 3  ! DMA for properties calc
  integer, parameter :: REP_SWEX_PROPERTIES_DMA_2 = 4  ! same, but with Bessel averaging
  integer, parameter :: REP_SWEX_POL_EMB_DMA_1 = 5     ! DMA for polarisable embedding
  integer, parameter :: REP_SWEX_POL_EMB_DMA_2 = 6     ! same, but with Bessel averaging
  integer, parameter :: REP_SWEX_POL_EMB_AUX_DMA_1 = 7 ! DMA for pol emb, auxiliary density
  integer, parameter :: REP_SWEX_POL_EMB_AUX_DMA_2 = 8 ! same, but with Bessel averaging

  ! jd: Here we acknowledge the fact that HFX_OTHER SW_EX's do not have Bb-Cc
  !     symmetry. All other SW_EX's are symmetrical under Bb-Cc label exchange.
  logical, parameter :: SWEX_BB_CC_SYMMETRIC(1:REP_N_SWEXES) = &
       (/.true., .false., .true., .true., .true., .true., .true., .true./)

  character, parameter :: dma_dummy_atom_symbol = 'J'

  ! jd: Identifiers for DFTB methods
  integer, parameter :: DFTB_GFN0 = 0
  integer, parameter :: DFTB_GFN1 = 1
  integer, parameter :: DFTB_GFN2 = 2

  character(len=*), parameter :: DFTB_DEFAULT_METHOD_PARAM_FILES(0:2) = &
      (/ 'param_gfn0-xtb.txt', 'param_gfn1-xtb.txt', 'param_gfn2-xtb.txt' /)
  character(len=*), parameter :: DFTB_DEFAULT_COMMON_PARAM_FILE = &
      'param_gfn_common.txt'

end module constants


