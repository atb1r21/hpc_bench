! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!
!  E l e c t r o n i c   S t r u c t u r e   D a t a   F o r m a t
!  ---------------------------------------------------------------
!
!                             E S D F
!                             =======
!
!  Author: Chris J. Pickard (c)
!  Email : cp@min.uni-kiel.de
!  Place : Kiel, Germany
!  Date  : 5/6th August 1999
!
!  Summary
!  -------
!
!  This module is designed to simplify and enhance the input of data into
!  electronic structure codes (for example, CASTEP). It works from a
!  highly flexible input file. data is input in a "label <value>"
!  fashion, and input is independent of the ordering of the input
!  file. An important feature is the requirement that most inputs require
!  default settings to be supplied within the main program calling
!  ESDF. This means that rarely used variables will not clutter everyday
!  input files, and, even more usefully, "intelligence" may be built into
!  the main code as the defaults may be dependent of other set
!  variables. Block data may also be read in. Another important feature
!  is the ability to define "physical" values. This means that the input
!  files need not depend on the internal physical units used by the main
!  program.
!
!
!  History
!  -------
!
!  ESDF has been written from scratch in F90, but is heavily based
!  (especially for the concept) on the FDF package developed by Alberto
!  Garcia and Jose Soler. It is not as "flexible" as FDF - there is no
!  provision for branching to further input files. This simplifies the
!  code, and I hope that it is still useful without this feature. Also,
!  the input and defaults are not dumped to a output file currently. I've
!  not found this a hindrance as of now.
!
!
!  Future
!  ------
!
!  My intention is to make this release available to Alberto Garcia and
!  Jose Soler for their comments. It might be a good idea to use this as
!  a base for fully converting the FDF package to F90. Or it may remain
!  as a cut down version of FDF. I certainly hope that a package of the
!  FDF sort becomes widely used in the electronic structure community. My
!  experience has been very positive.
!
!
!  Usage
!  -----
!
!  First, "Use esdf" wherever you wish to make use of its features. In
!  the main program call the initialisation routine: call
!  esdf_init('input.esdf'). "input.esdf" is the name of the input file -
!  it could be anything. This routine opens the input file, and reads
!  into a dynamically allocated storage array. The comments and blank
!  lines are stripped out. You are now ready to use the
!  esdf_functions. For example, if you want to read in the number of
!  atoms in your calculation, you would use: natom =
!  esdf_integer('NumberOfAtoms',1), where 'NumberOfAtoms' is the label to
!  search for in the input file, and '1' is the default. call esdf_close to
!  deallocate the data arrays. You may then open another input file using
!  esdf_init. It is not currently possible to open more that on input
!  file at one time.
!
!
!  Syntax
!  ------
!
!  The input file can contain comments. These are defined as anything to
!  the right of, and including, '#', ';', or '!'. It is straightforward
!  to modify the routine to accept further characters. Blank lines are
!  ignored -- use comments and blank lines to make you input file
!  readable.
!
!  The "labels" are case insensitive (e.g. unitCell is equivalent to
!  UnItceLL) and punctuation insensitive (unit.cell is equivalent to
!  unit_cell is equivalent to unitcell). Punctuation characters are '.',
!  '_', and '-' at the moment. Again - use this feature to improve
!  readability.
!
!  The following are equivalent ways of defining a physical quantity:
!
!  "AgeOfUniverse = 24.d0 s" or "AgeOfUniverse : 24.d0 S" or
!  "AgeOfUniverse 24.d0 S"
!
!  It would be read in by the main program in the following way:
!
!  aou = esdf_physical('ageofuniverse',77.d0,ns)
!
!  "aou" is the double precision variable, 77.d0 is the default number of
!  "ns" or nanoseconds. 24s will be converted automatically to its
!  equivalent number of nanoseconds.
!
!  Block data should be placed in the input file as follows:
!
!  %block cellvectors
!  1.0 1.0 0.0
!  0.0 1.0 1.0
!  1.0 0.0 1.0
!  %endblock cellvectors
!
!  And it may be read:
!
!    if (esdf_block('CellVectors',nlines))
!      if (nlines.ne.3) then (... break out here if the incorrect number
!  of lines)
!      do i=1,nlines
!        read(block_data(i),*) x,y,z
!      end do
!    end if
!
!
!  List of functions
!  -----------------
!
!  Self explanatory:
!
!  esdf_string(label,default)
!  esdf_integer(label,default)
!  esdf_single(label,default)
!  esdf_double(label,default)
!  esdf_physical(label,default,unit)
!
!  A little more explanation:
!
!  esdf_defined(label) is true if "label" found, false otherwise
!
!  esdf_boolean(label,default) is true if "label yes/true/t (case/punct.insens)
!                              is false if"label no/false/f (case/punct.insens)
!
!  The Help feature
!  ----------------
!
!  The routine "esdf_help(helpword,searchword)" can be used to access the
!  information contained within the "esdf_key_mod" module.
!
!  if "helpword" is "search" (case insensitive), then all variables whose
!  description contains "searchword" will be output.
!
!  if "helpword" is "basic", "inter", "expert" or "dummy" the varibles of
!  that type will be displayed.
!
!  if "helpword" is one of the valid labels, then a description of this
!  label will be output.
!
!
!  Finishing off
!  -------------
!
!  Three routines, "esdf_dump", "esdf_warnout" and "esdf_close", can be
!  used to finish the use of ESDF. "esdf_dump" outputs a file "ESDF.esdf"
!  which could be used as an input file for further runs. "esdf_warnout"
!  outputs ESDF warnings to screen, and "esdf_close" deallocates the
!  allocated ESDF arrays.
!
!  Contact the Author
!  ------------------
!
!  This code is under development, and the author would be very happy to
!  receive comments by email. Any use in a commercial software package is
!  forbidden without prior arrangement with the author (Chris J. Pickard).
!
!
!  Code cleanup by Nicholas Hine, 21 October 2011


module esdf

  Use esdf_key
  Use constants, only: DP, TWO_PI, AVOGADRO, ANGSTROM, e_SI, ATM_SI, k_B, &
       HARTREE_IN_DEBYE, HARTREE_IN_NS, HARTREE_IN_EVS, HARTREE_TO_KCAL_PER_MOL, &
       SPEED_OF_LIGHT_SI

  Implicit None

  ! Kind parameters

  !qoh: Remove non standard integer kind
  !integer, private, parameter :: I4B = selected_int_kind(9)
  !qoh: Use ONETEP's standard DP from constants
  !integer, private, parameter :: DP  = Kind(1.d0)
  integer, private, parameter :: SP  = Kind(1.0)

  ! Set the length of the lines

  integer, public, parameter :: llength=80
  integer, public, parameter :: esdf_maxstring = 3000
  !cks 18/4/2000 the size of ndump has been changed to 20000 from its original value of 200
  !ndmh 08/09/2008 the size of ndump has been further increased to 200000
  integer, private, parameter ::  nphys = 98, ndump = 200000
  integer, private :: nrecords,nwarns,nerrors,ndmp
  character(llength), private, dimension(:), allocatable :: llist,warns,errors,dump
  character(llength), private, dimension(:,:), allocatable :: tlist

  ! The public block data array
  character(llength), public, dimension(:), allocatable :: block_data

#ifndef OLD_UNITS
  ! ab: Derived units used in ESDF physical unit datastructure
  real(kind=DP), parameter :: ANG_SI=1.0D-10
  real(kind=DP), parameter :: BOHR_SI=ANG_SI/ANGSTROM
  real(kind=DP), parameter :: BOHR_PER_NS_SI=BOHR_SI*1.D9
  real(kind=DP), parameter :: BOHR_PER_PS_SI=BOHR_SI*1.D12
  real(kind=DP), parameter :: BOHR_PER_FS_SI=BOHR_SI*1.D15
  real(kind=DP), parameter :: eBOHR_SI=e_SI*BOHR_SI
  real(kind=DP), parameter :: INV_BOHR_SI=1.0D0/BOHR_SI
  real(kind=DP), parameter :: eANG_SI=e_SI*ANG_SI
  real(kind=DP), parameter :: EV_SI=e_SI
  real(kind=DP), parameter :: EV_PER_ANG_SI=EV_SI/ANG_SI
  real(kind=DP), parameter :: EV_PER_ANG3_SI=EV_SI/ANG_SI**3
  real(kind=DP), parameter :: meV_SI=1.0D-3*EV_SI
  real(kind=DP), parameter :: HARTREE_SI=HARTREE_IN_EVS*EV_SI
  real(kind=DP), parameter :: SI_HARTREE=1.0D0/HARTREE_SI
  real(kind=DP), parameter :: mHARTREE_SI=HARTREE_SI*1.0D-3
  real(kind=DP), parameter :: HA_BOHR12_SI=HARTREE_SI*BOHR_SI**12
  real(kind=DP), parameter :: HA_PER_BOHR_SI=HARTREE_SI/BOHR_SI
  real(kind=DP), parameter :: HA_PER_BOHR2_SI=HARTREE_SI/BOHR_SI**2
  real(kind=DP), parameter :: HA_PER_BOHR3_SI=HARTREE_SI/BOHR_SI**3
  real(kind=DP), parameter :: RYD_SI=0.5D0*HARTREE_SI
  real(kind=DP), parameter :: mRYD_SI=1.0D-3*RYD_SI
  real(kind=DP), parameter :: RYD_PER_BOHR_SI=RYD_SI/BOHR_SI
  real(kind=DP), parameter :: RYD_PER_BOHR3_SI=RYD_SI/BOHR_SI**3
  real(kind=DP), parameter :: RYD_FS2_SI=RYD_SI*1.0D-30
  real(kind=DP), parameter :: AUT_SI=HARTREE_IN_NS*1.0D-9
  real(kind=DP), parameter :: AUV_SI=BOHR_SI/AUT_SI
  real(kind=DP), parameter :: KJ_PER_MOL_SI=1.0D3/AVOGADRO
  real(kind=DP), parameter :: ATOMIC_MASS_UNIT=1.0D-3/AVOGADRO
  real(kind=DP), parameter :: KCAL_PER_MOL_SI=HARTREE_SI/HARTREE_TO_KCAL_PER_MOL
  real(kind=DP), parameter :: PLANCK_SI=TWO_PI*AUT_SI*HARTREE_SI
  real(kind=DP), parameter :: E_THZ_SI=PLANCK_SI*1.0D12
  real(kind=DP), parameter :: INV_CM_SI=PLANCK_SI*SPEED_OF_LIGHT_SI*1.0D2
  real(kind=DP), parameter :: DEBYE_SI=e_SI*BOHR_SI/HARTREE_IN_DEBYE
  real(kind=DP), parameter :: SI_DEBYE=1.0D0/DEBYE_SI
  real(kind=DP), parameter :: EF_AU_SI=HARTREE_IN_EVS/BOHR_SI
  real(kind=DP), parameter :: EF_RY_SI=0.5_DP*EF_AU_SI
  real(kind=DP), parameter :: kB_SI=HARTREE_SI*k_B
#endif

  ! Set the physical units database
  type phys_unit
     character(10) :: d   ! d - dimension
     character(20) :: n   ! n - name
     real(kind=DP) :: u   ! u - unit
  end type phys_unit

  type(phys_unit), private, dimension(nphys) :: phy

  !
  !     We allow case variations in the units. This could be dangerous
  !     (meV --> MeV !!) in real life, but not in this restricted
  !     field.
  !
  ! m - mass l - length t - time e - energy f - force p - pressure c- charge
  ! d - dipole mom - mom inert ef - efield s - surface tension
  !
  ! jd: 2010/06/08: Added a new category: surface tension, with two
  !                 units: ha/bohr**2 and J/m**2
  ! jd: 2014/07/08: Added a new category: 1/energy
  ! jd: 2014/07/08: Added a new category: energy * length^12 (el12)
  ! jd: 2014/07/08: Increased length of name to 20.
  ! ab: 2020/07/27: Added a new category, potential (pot) with two units:
  !                 volts (v) and hartree/elementary charge (ha/e)
  ! ab: 2022/09/06: Converted hard-coded numbers into derived parameters.
#ifndef OLD_UNITS
  data phy(1)%d /'m'/, phy(1)%n /'kg'/, phy(1)%u /1.d0/
  data phy(2)%d /'m'/, phy(2)%n /'g'/, phy(2)%u /1.d-3/
  data phy(3)%d /'m'/, phy(3)%n /'amu'/, phy(3)%u /ATOMIC_MASS_UNIT/
  data phy(4)%d /'l'/, phy(4)%n /'m'/, phy(4)%u /1.d0/
  data phy(5)%d /'l'/, phy(5)%n /'nm'/, phy(5)%u /1.d-9/
  data phy(6)%d /'l'/, phy(6)%n /'ang'/, phy(6)%u /1.d-10/
  data phy(7)%d /'l'/, phy(7)%n /'bohr'/, phy(7)%u /BOHR_SI/
  data phy(8)%d /'t'/, phy(8)%n /'aut'/, phy(8)%u /AUT_SI/
  data phy(9)%d /'t'/, phy(9)%n /'s'/, phy(9)%u /1.d0/
  data phy(10)%d /'t'/, phy(10)%n /'ns'/, phy(10)%u /1.d-9/
  data phy(11)%d /'t'/, phy(11)%n /'ps'/, phy(11)%u /1.d-12/
  data phy(12)%d /'t'/, phy(12)%n /'fs'/, phy(12)%u /1.d-15/
  data phy(13)%d /'e'/, phy(13)%n /'j'/, phy(13)%u /1.d0/
  data phy(14)%d /'e'/, phy(14)%n /'erg'/, phy(14)%u /1.d-7/
  data phy(15)%d /'e'/, phy(15)%n /'ev'/, phy(15)%u /EV_SI/
  data phy(16)%d /'e'/, phy(16)%n /'mev'/, phy(16)%u /meV_SI/
  data phy(17)%d /'e'/, phy(17)%n /'ry'/, phy(17)%u /RYD_SI/
  data phy(18)%d /'e'/, phy(18)%n /'mry'/, phy(18)%u /mRYD_SI/
  data phy(19)%d /'e'/, phy(19)%n /'hartree'/, phy(19)%u /HARTREE_SI/
  data phy(20)%d /'e'/, phy(20)%n /'kcal/mol'/, phy(20)%u /KCAL_PER_MOL_SI/
  data phy(21)%d /'e'/, phy(21)%n /'mhartree'/, phy(21)%u /mHARTREE_SI/
  data phy(22)%d /'e'/, phy(22)%n /'kj/mol'/, phy(22)%u /KJ_PER_MOL_SI/
  data phy(23)%d /'e'/, phy(23)%n /'hz'/, phy(23)%u /PLANCK_SI/
  data phy(24)%d /'e'/, phy(24)%n /'thz'/, phy(24)%u /E_THZ_SI/
  data phy(25)%d /'e'/, phy(25)%n /'cm-1'/, phy(25)%u /INV_CM_SI/
  data phy(26)%d /'e'/, phy(26)%n /'cm^-1'/, phy(26)%u /INV_CM_SI/
  data phy(27)%d /'e'/, phy(27)%n /'cm**-1'/, phy(27)%u /INV_CM_SI/
  data phy(28)%d /'f'/, phy(28)%n /'n'/, phy(28)%u /1.d0/
  data phy(29)%d /'f'/, phy(29)%n /'ev/ang'/, phy(29)%u /EV_PER_ANG_SI/
  data phy(30)%d /'f'/, phy(30)%n /'ry/bohr'/, phy(30)%u /RYD_PER_BOHR_SI/
  data phy(31)%d /'l'/, phy(31)%n /'cm'/, phy(31)%u /1.d-2/
  data phy(32)%d /'p'/, phy(32)%n /'pa'/, phy(32)%u /1.d0/
  data phy(33)%d /'p'/, phy(33)%n /'mpa'/, phy(33)%u /1.d6/
  data phy(34)%d /'p'/, phy(34)%n /'gpa'/, phy(34)%u /1.d9/
  data phy(35)%d /'p'/, phy(35)%n /'atm'/, phy(35)%u /ATM_SI/
  data phy(36)%d /'p'/, phy(36)%n /'bar'/, phy(36)%u /1.d5/
  data phy(37)%d /'p'/, phy(37)%n /'mbar'/, phy(37)%u /1.d11/
  data phy(38)%d /'p'/, phy(38)%n /'ry/bohr**3'/, phy(38)%u /RYD_PER_BOHR3_SI/
  data phy(39)%d /'p'/, phy(39)%n /'ev/ang**3'/, phy(39)%u /EV_PER_ANG3_SI/
  data phy(40)%d /'c'/, phy(40)%n /'c'/, phy(40)%u /1.d0/
  data phy(41)%d /'c'/, phy(41)%n /'e'/, phy(41)%u /e_SI/
  data phy(42)%d /'d'/, phy(42)%n /'c*m'/, phy(42)%u /1.d0/
  data phy(43)%d /'d'/, phy(43)%n /'d'/, phy(43)%u /DEBYE_SI/
  data phy(44)%d /'d'/, phy(44)%n /'debye'/, phy(44)%u /DEBYE_SI/
  data phy(45)%d /'d'/, phy(45)%n /'e*bohr'/, phy(45)%u /eBOHR_SI/
  data phy(46)%d /'d'/, phy(46)%n /'e*ang'/, phy(46)%u /eANG_SI/
  data phy(47)%d /'mom'/, phy(47)%n /'kg*m**2'/, phy(47)%u /1.d0/
  data phy(48)%d /'mom'/, phy(48)%n /'ry*fs**2'/, phy(48)%u /RYD_FS2_SI/
  data phy(49)%d /'ef'/, phy(49)%n /'v/m'/, phy(49)%u /1.d0/
  data phy(50)%d /'ef'/, phy(50)%n /'v/nm'/, phy(50)%u /1.d9/
  data phy(51)%d /'ef'/, phy(51)%n /'v/ang'/, phy(51)%u /1.d10/
  data phy(52)%d /'ef'/, phy(52)%n /'v/bohr'/, phy(52)%u /INV_BOHR_SI/
  data phy(53)%d /'ef'/, phy(53)%n /'ry/bohr/e'/, phy(53)%u /EF_RY_SI/
  data phy(54)%d /'ef'/, phy(54)%n /'ha/bohr/e'/, phy(54)%u /EF_AU_SI/
  data phy(55)%d /'e'/, phy(55)%n /'k'/, phy(55)%u /kB_SI/
  data phy(56)%d /'f'/, phy(56)%n /'ha/bohr'/, phy(56)%u /HA_PER_BOHR_SI/
  data phy(57)%d /'p'/, phy(57)%n /'ha/bohr**3'/, phy(57)%u /HA_PER_BOHR3_SI/
  data phy(58)%d /'1/l'/, phy(58)%n /'1/ang'/, phy(58)%u /1.d10/
  data phy(59)%d /'1/l'/, phy(59)%n /'1/bohr'/, phy(59)%u /INV_BOHR_SI/
  data phy(60)%d /'1/l'/, phy(60)%n /'1/nm'/, phy(60)%u /1.d9/
  data phy(61)%d /'s'/, phy(61)%n /'ha/bohr**2'/, phy(61)%u /HA_PER_BOHR2_SI/
  data phy(62)%d /'s'/, phy(62)%n /'j/m**2'/, phy(62)%u /1.d0/
  data phy(63)%d /'v'/, phy(63)%n /'m/s'/, phy(63)%u /1.d0/
  data phy(64)%d /'v'/, phy(64)%n /'auv'/, phy(64)%u /AUV_SI/
  data phy(65)%d /'v'/, phy(65)%n /'ang/ns'/, phy(65)%u /0.1d0/
  data phy(66)%d /'v'/, phy(66)%n /'ang/ps'/, phy(66)%u /1.d2/
  data phy(67)%d /'v'/, phy(67)%n /'ang/fs'/, phy(67)%u /1.d5/
  data phy(68)%d /'v'/, phy(68)%n /'bohr/ns'/, phy(68)%u /BOHR_PER_NS_SI/
  data phy(69)%d /'v'/, phy(69)%n /'bohr/ps'/, phy(69)%u /BOHR_PER_PS_SI/
  data phy(70)%d /'v'/, phy(70)%n /'bohr/fs'/, phy(70)%u /BOHR_PER_FS_SI/
  data phy(71)%d /'1/l'/, phy(71)%n /'bohr-1'/, phy(71)%u /INV_BOHR_SI/
  data phy(72)%d /'1/l'/, phy(72)%n /'bohr^-1'/, phy(72)%u /INV_BOHR_SI/
  data phy(73)%d /'1/l'/, phy(73)%n /'bohr**-1'/, phy(73)%u /INV_BOHR_SI/
  data phy(74)%d /'1/l'/, phy(74)%n /'ang-1'/, phy(74)%u /1.d10/
  data phy(75)%d /'1/l'/, phy(75)%n /'ang^-1'/, phy(75)%u /1.d10/
  data phy(76)%d /'1/l'/, phy(76)%n /'ang**-1'/, phy(76)%u /1.d10/
  data phy(77)%d /'1/l'/, phy(77)%n /'nm-1'/, phy(77)%u /1.d9/
  data phy(78)%d /'1/l'/, phy(78)%n /'nm^-1'/, phy(78)%u /1.d9/
  data phy(79)%d /'1/l'/, phy(79)%n /'nm**-1'/, phy(79)%u /1.d9/
  data phy(80)%d /'e'/, phy(80)%n /'ha'/, phy(80)%u /HARTREE_SI/
  data phy(81)%d /'e'/, phy(81)%n /'mha'/, phy(81)%u /mHARTREE_SI/
  data phy(82)%d /'te'/, phy(82)%n /'kelvin'/, phy(82)%u /1.d0/
  data phy(83)%d /'fr'/, phy(83)%n /'1/cm'/, phy(83)%u /1.d2/
  data phy(84)%d /'1/e'/, phy(84)%n /'1/ha'/, phy(84)%u /SI_HARTREE/
  data phy(85)%d /'1/e'/, phy(85)%n /'1/hartree'/, phy(85)%u /SI_HARTREE/
  data phy(86)%d /'1/e'/, phy(86)%n /'ha-1'/, phy(86)%u /SI_HARTREE/
  data phy(87)%d /'1/e'/, phy(87)%n /'ha^-1'/, phy(87)%u /SI_HARTREE/
  data phy(88)%d /'1/e'/, phy(88)%n /'ha**-1'/, phy(88)%u /SI_HARTREE/
  data phy(89)%d /'el12'/, phy(89)%n /'ha*bohr^12'/, phy(89)%u /HA_BOHR12_SI/
  data phy(90)%d /'el12'/, phy(90)%n /'ha*bohr**12'/, phy(90)%u /HA_BOHR12_SI/
  data phy(91)%d /'el12'/, phy(91)%n /'hartree*bohr^12'/, phy(91)%u /HA_BOHR12_SI/
  data phy(92)%d /'el12'/, phy(92)%n /'hartree*bohr**12'/, phy(92)%u /HA_BOHR12_SI/
  data phy(93)%d /'s'/, phy(93)%n /'n/m'/, phy(93)%u /1.d0/
  data phy(94)%d /'invd'/, phy(94)%n /'debye^-1'/, phy(94)%u /SI_DEBYE/
  data phy(95)%d /'invd'/, phy(95)%n /'debye**-1'/, phy(95)%u /SI_DEBYE/
  data phy(96)%d /'invd'/, phy(96)%n /'1/debye'/, phy(96)%u /SI_DEBYE/
  data phy(97)%d /'pot'/, phy(97)%n /'v'/, phy(97)%u /1.d0/
  data phy(98)%d /'pot'/, phy(98)%n /'ha/e'/, phy(98)%u /HARTREE_IN_EVS/
  ! jd: If you add anything to this list, make sure to use lowercase letters
  !     (e.g. c/s and not C/s for coulomb/second).
#else
  data phy(1)%d /'m'/;data phy(1)%n /'kg'/;data phy(1)%u /1.d0/
  data phy(2)%d /'m'/;data phy(2)%n /'g'/;data phy(2)%u /1.d-3/
  data phy(3)%d /'m'/;data phy(3)%n /'amu'/;data phy(3)%u /1.66054d-27/
  data phy(4)%d /'l'/;data phy(4)%n /'m'/;data phy(4)%u /1.d0/
  data phy(5)%d /'l'/;data phy(5)%n /'nm'/;data phy(5)%u /1.d-9/
  data phy(6)%d /'l'/;data phy(6)%n /'ang'/;data phy(6)%u /1.d-10/
  data phy(7)%d /'l'/;data phy(7)%n /'bohr'/;data phy(7)%u /0.529177d-10/
  data phy(8)%d /'t'/;data phy(8)%n /'aut'/;data phy(8)%u /2.41888468d-17/
  data phy(9)%d /'t'/;data phy(9)%n /'s'/;data phy(9)%u /1.d0/
  data phy(10)%d /'t'/;data phy(10)%n /'ns'/;data phy(10)%u /1.d-9/
  data phy(11)%d /'t'/;data phy(11)%n /'ps'/;data phy(11)%u /1.d-12/
  data phy(12)%d /'t'/;data phy(12)%n /'fs'/;data phy(12)%u /1.d-15/
  data phy(13)%d /'e'/;data phy(13)%n /'j'/;data phy(13)%u /1.d0/
  data phy(14)%d /'e'/;data phy(14)%n /'erg'/;data phy(14)%u /1.d-7/
  data phy(15)%d /'e'/;data phy(15)%n /'ev'/;data phy(15)%u /1.60219d-19/
  data phy(16)%d /'e'/;data phy(16)%n /'mev'/;data phy(16)%u /1.60219d-22/
  data phy(17)%d /'e'/;data phy(17)%n /'ry'/;data phy(17)%u /2.17991d-18/
  data phy(18)%d /'e'/;data phy(18)%n /'mry'/;data phy(18)%u /2.17991d-21/
  data phy(19)%d /'e'/;data phy(19)%n /'hartree'/;data phy(19)%u /4.35982d-18/
  data phy(20)%d /'e'/;data phy(20)%n /'kcal/mol'/;data phy(20)%u /6.94780d-21/
  data phy(21)%d /'e'/;data phy(21)%n /'mhartree'/;data phy(21)%u /4.35982d-21/
  data phy(22)%d /'e'/;data phy(22)%n /'kj/mol'/;data phy(22)%u /1.6606d-21/
  data phy(23)%d /'e'/;data phy(23)%n /'hz'/;data phy(23)%u /6.6262d-34/
  data phy(24)%d /'e'/;data phy(24)%n /'thz'/;data phy(24)%u /6.6262d-22/
  data phy(25)%d /'e'/;data phy(25)%n /'cm-1'/;data phy(25)%u /1.986d-23/
  data phy(26)%d /'e'/;data phy(26)%n /'cm^-1'/;data phy(26)%u /1.986d-23/
  data phy(27)%d /'e'/;data phy(27)%n /'cm**-1'/;data phy(27)%u /1.986d-23/
  data phy(28)%d /'f'/;data phy(28)%n /'n'/;data phy(28)%u /1.d0/
  data phy(29)%d /'f'/;data phy(29)%n /'ev/ang'/;data phy(29)%u /1.60219d-9/
  data phy(30)%d /'f'/;data phy(30)%n /'ry/bohr'/;data phy(30)%u /4.11943d-8/
  data phy(31)%d /'l'/;data phy(31)%n /'cm'/;data phy(31)%u /1.d-2/
  data phy(32)%d /'p'/;data phy(32)%n /'pa'/;data phy(32)%u /1.d0/
  data phy(33)%d /'p'/;data phy(33)%n /'mpa'/;data phy(33)%u /1.d6/
  data phy(34)%d /'p'/;data phy(34)%n /'gpa'/;data phy(34)%u /1.d9/
  data phy(35)%d /'p'/;data phy(35)%n /'atm'/;data phy(35)%u /1.01325d5/
  data phy(36)%d /'p'/;data phy(36)%n /'bar'/;data phy(36)%u /1.d5/
  data phy(37)%d /'p'/;data phy(37)%n /'mbar'/;data phy(37)%u /1.d11/
  data phy(38)%d /'p'/;data phy(38)%n /'ry/bohr**3'/;data phy(38)%u /1.47108d13/
  data phy(39)%d /'p'/;data phy(39)%n /'ev/ang**3'/;data phy(39)%u /1.60219d11/
  data phy(40)%d /'c'/;data phy(40)%n /'c'/;data phy(40)%u /1.d0/
  data phy(41)%d /'c'/;data phy(41)%n /'e'/;data phy(41)%u /1.602177d-19/
  data phy(42)%d /'d'/;data phy(42)%n /'c*m'/;data phy(42)%u /1.d0/
  data phy(43)%d /'d'/;data phy(43)%n /'d'/;data phy(43)%u /3.33564d-30/
  data phy(44)%d /'d'/;data phy(44)%n /'debye'/;data phy(44)%u /3.33564d-30/
  data phy(45)%d /'d'/;data phy(45)%n /'e*bohr'/;data phy(45)%u /8.47835d-30/
  data phy(46)%d /'d'/;data phy(46)%n /'e*ang'/;data phy(46)%u /1.602177d-29/
  data phy(47)%d /'mom'/;data phy(47)%n /'kg*m**2'/;data phy(47)%u /1.d0/
  data phy(48)%d /'mom'/;data phy(48)%n /'ry*fs**2'/;data phy(48)%u /2.1799d-48/
  data phy(49)%d /'ef'/;data phy(49)%n /'v/m'/;data phy(49)%u /1.d0/
  data phy(50)%d /'ef'/;data phy(50)%n /'v/nm'/;data phy(50)%u /1.d9/
  data phy(51)%d /'ef'/;data phy(51)%n /'v/ang'/;data phy(51)%u /1.d10/
  data phy(52)%d /'ef'/;data phy(52)%n /'v/bohr'/;data phy(52)%u /1.8897268d10/
  data phy(53)%d /'ef'/;data phy(53)%n /'ry/bohr/e'/;data phy(53)%u /2.5711273d11/
  data phy(54)%d /'ef'/;data phy(54)%n /'ha/bohr/e'/;data phy(54)%u /5.1422546d11/
  data phy(55)%d /'e'/;data phy(55)%n /'k'/;data phy(55)%u /1.38066d-23/
  data phy(56)%d /'f'/;data phy(56)%n /'ha/bohr'/;data phy(56)%u /8.23886d-8/
  data phy(57)%d /'p'/;data phy(57)%n /'ha/bohr**3'/;data phy(57)%u /2.94216d13/
  data phy(58)%d /'1/l'/;data phy(58)%n /'1/ang'/;data phy(58)%u /1.d10/
  data phy(59)%d /'1/l'/;data phy(59)%n /'1/bohr'/;data phy(59)%u /1.889727d10/
  data phy(60)%d /'1/l'/;data phy(60)%n /'1/nm'/;data phy(60)%u /1.d9/
  data phy(61)%d /'s'/;data phy(61)%n /'ha/bohr**2'/;data phy(61)%u /1.556894d3/
  data phy(62)%d /'s'/;data phy(62)%n /'j/m**2'/;data phy(62)%u /1.d0/
  data phy(63)%d /'v'/;data phy(63)%n /'m/s'/;data phy(63)%u /1.d0/
  data phy(64)%d /'v'/;data phy(64)%n /'auv'/;data phy(64)%u /2.1876912633d6/
  data phy(65)%d /'v'/;data phy(65)%n /'ang/ns'/;data phy(65)%u /0.1d0/
  data phy(66)%d /'v'/;data phy(66)%n /'ang/ps'/;data phy(66)%u /1.d2/
  data phy(67)%d /'v'/;data phy(67)%n /'ang/fs'/;data phy(67)%u /1.d5/
  data phy(68)%d /'v'/;data phy(68)%n /'bohr/ns'/;data phy(68)%u /0.529177d-1/
  data phy(69)%d /'v'/;data phy(69)%n /'bohr/ps'/;data phy(69)%u /0.529177d2/
  data phy(70)%d /'v'/;data phy(70)%n /'bohr/fs'/;data phy(70)%u /0.529177d5/
  data phy(71)%d /'1/l'/;data phy(71)%n /'bohr-1'/;data phy(71)%u /1.889727d10/
  data phy(72)%d /'1/l'/;data phy(72)%n /'bohr^-1'/;data phy(72)%u /1.889727d10/
  data phy(73)%d /'1/l'/;data phy(73)%n /'bohr**-1'/;data phy(73)%u /1.889727d10/
  data phy(74)%d /'1/l'/;data phy(74)%n /'ang-1'/;data phy(74)%u /1.d10/
  data phy(75)%d /'1/l'/;data phy(75)%n /'ang^-1'/;data phy(75)%u /1.d10/
  data phy(76)%d /'1/l'/;data phy(76)%n /'ang**-1'/;data phy(76)%u /1.d10/
  data phy(77)%d /'1/l'/;data phy(77)%n /'nm-1'/;data phy(77)%u /1.d9/
  data phy(78)%d /'1/l'/;data phy(78)%n /'nm^-1'/;data phy(78)%u /1.d9/
  data phy(79)%d /'1/l'/;data phy(79)%n /'nm**-1'/;data phy(79)%u /1.d9/
  data phy(80)%d /'e'/;data phy(80)%n /'ha'/;data phy(80)%u /4.35982d-18/
  data phy(81)%d /'e'/;data phy(81)%n /'mha'/;data phy(81)%u /4.35982d-21/
  data phy(82)%d /'te'/;data phy(82)%n /'kelvin'/;data phy(82)%u /1.d0/
  data phy(83)%d /'fr'/;data phy(83)%n /'1/cm'/;data phy(83)%u /1.d2/
  data phy(84)%d /'1/e'/;data phy(84)%n /'1/ha'/;data phy(84)%u /2.2936727d+17/
  data phy(85)%d /'1/e'/;data phy(85)%n /'1/hartree'/;data phy(85)%u /2.2936727d+17/
  data phy(86)%d /'1/e'/;data phy(86)%n /'ha-1'/;data phy(86)%u /2.2936727d+17/
  data phy(87)%d /'1/e'/;data phy(87)%n /'ha^-1'/;data phy(87)%u /2.2936727d+17/
  data phy(88)%d /'1/e'/;data phy(88)%n /'ha**-1'/;data phy(88)%u /2.2936727d+17/
  data phy(89)%d /'el12'/;data phy(89)%n /'ha*bohr^12'/;data phy(89)%u /2.1022293D-141/
  data phy(90)%d /'el12'/;data phy(90)%n /'ha*bohr**12'/;data phy(90)%u /2.1022293D-141/
  data phy(91)%d /'el12'/;data phy(91)%n /'hartree*bohr^12'/;data phy(91)%u /2.1022293D-141/
  data phy(92)%d /'el12'/;data phy(92)%n /'hartree*bohr**12'/;data phy(92)%u /2.1022293D-141/
  data phy(93)%d /'s'/;data phy(93)%n /'n/m'/;data phy(93)%u /1.d0/
  data phy(94)%d /'invd'/;data phy(94)%n /'debye^-1'/;data phy(94)%u /2.99792543D+29/
  data phy(95)%d /'invd'/;data phy(95)%n /'debye**-1'/;data phy(95)%u /2.99792543D+29/
  data phy(96)%d /'invd'/;data phy(96)%n /'1/debye'/;data phy(96)%u /2.99792543D+29/
  data phy(97)%d /'pot'/;data phy(97)%n /'v'/;data phy(97)%u /1.d0/
  data phy(98)%d /'pot'/;data phy(98)%n /'ha/e'/;data phy(98)%u /HARTREE_IN_EVS/
#endif

contains

  subroutine esdf_init(filename,error)

    use utils, only: utils_alloc_check, utils_abort
    implicit none

    character(*), intent(in) :: filename
    integer, intent(out) :: error

    ! Local

    integer, parameter :: ncomm=3,ndiv=4
    integer :: unit,ierr,i,j,ic,nt,ndef
    character(llength) :: cjunk,ctemp,ctemp2
    character(1) :: comment(ncomm),divide(ndiv)
    ! jd: coding up tab as a constant rather than a literal prevents nasty
    !     breakage due to editors replacing tabs by spaces
    character(1), parameter :: tab = achar(9)
    logical :: inblock

    ! smmd : local variable to handle inclusion of external files
    integer :: ixf, xfile_num_loc,xfile_num_tot, xfile_unit(llength)
    character(llength) :: xfile_name(llength)

    ! Define comment characters
    data comment /'#',';','!'/
    data divide /' ','=',':',tab/

    ! Set error flag to success by default
    error = 0
    call esdf_key_init

    ! "reduce" the keyword list for comparison
    do i = 1,numkw
       ctemp = kw(i)%label
       ! qoh: Explicit character trunctation to avoid compiler warning
       ctemp2 = esdf_reduce(ctemp)
       kw(i)%label = ctemp2(1:32)
    end do

    ! open the esdf file
    call esdf_file(unit,filename,ierr)
    cjunk = 'Unable to open the input file "'//trim(filename)//'"'
    if (ierr.eq.1) then
       call utils_abort(cjunk)
    else

       ! Count the number of records (excluding blank lines, commented lines, and included lines)
       nrecords = 0
       xfile_num_tot = 0

       do
          read(unit,'(a)',end=100,err=102,iostat=ierr) cjunk
          do j=1,ncomm
             ic=index(cjunk,comment(j))
             if (ic.gt.0) cjunk(ic:) = ' '
          end do

          ! smmd: Check for inclusion of external files via the "includefile" keyword
          ixf=index(esdf_reduce(cjunk),'includefile')
          if (len_trim(cjunk).gt.0 .and. ixf.eq.0) then
             nrecords = nrecords + 1
          end if

          ! smmd:  Count the number of records coming from external files
          if (ixf.gt.0) then
             xfile_num_loc=0
             ctemp = adjustl(cjunk)
             ixf = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
             ctemp = adjustl(ctemp(ixf:))
             do while(len_trim(ctemp).gt.0)
                ixf = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
                if (ixf.gt.1) then
                   xfile_num_tot=xfile_num_tot+1
                   xfile_num_loc=xfile_num_loc+1
                   xfile_name(xfile_num_tot) = adjustl(ctemp(:ixf-1))
                end if
                ctemp = adjustl(ctemp(ixf+1:))
             end do

             do ixf=xfile_num_tot-xfile_num_loc+1,xfile_num_tot
                call esdf_file(xfile_unit(ixf),xfile_name(ixf),ierr)
                cjunk = 'Unable to open secondary input file "'//trim(xfile_name(ixf))//'"'
                if (ierr.eq.1) then
                   call utils_abort(cjunk)
                else
                   do
                      read(xfile_unit(ixf),'(a)',end=112,err=102,iostat=ierr) cjunk
                      do j=1,ncomm
                         ic=index(cjunk,comment(j))
                         if (ic.gt.0) cjunk(ic:) = ' '
                      end do
                      if (len_trim(cjunk).gt.0) then
                         nrecords = nrecords + 1
                      end if
                   end do
112                rewind(xfile_unit(ixf))
                end if
             end do
          end if

       end do
100    rewind(unit)

    end if

    ! allocate the array to hold the records and tokens
    allocate(llist(nrecords),stat=ierr)
    call utils_alloc_check('esdf_[init]','llist',ierr)
    allocate(block_data(nrecords),stat=ierr)
    call utils_alloc_check('esdf_[init]','block_data',ierr)
    allocate(tlist(llength,nrecords),stat=ierr)
    call utils_alloc_check('esdf_[init]','tlist',ierr)
    allocate(warns(nrecords),stat=ierr)
    call utils_alloc_check('esdf_[init]','warns',ierr)
    allocate(errors(nrecords),stat=ierr)
    call utils_alloc_check('esdf_[init]','errors',ierr)
    allocate(dump(ndump),stat=ierr)
    call utils_alloc_check('esdf_[init]','dump',ierr)

    ! Set the number of warnings and errors to zero
    nwarns = 0 ; warns = ' ' ; nerrors = 0 ; errors = ' ' ; ndmp = 0 ; dump = ' '

    if (ierr.eq.1) then
       nrecords = 0
    else

       ! Read in the records
       nrecords = 0
       xfile_num_tot = 0
       do
          read(unit,'(a)',end=101,err=103,iostat=ierr) cjunk
          do j=1,ncomm
             ic=index(cjunk,comment(j))
             if (ic.gt.0) cjunk(ic:) = ' '
          end do

          ixf=index(esdf_reduce(cjunk),'includefile')

          if (len_trim(cjunk).gt.0 .and. ixf.eq.0) then
             nrecords=nrecords+1
             llist(nrecords) = adjustl(cjunk)
          end if

          ! smmd:  Read in the records coming from external files
          if (ixf.gt.0) then
             xfile_num_loc=0
             ctemp = adjustl(cjunk)
             ixf = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
             ctemp = adjustl(ctemp(ixf:))
             do while(len_trim(ctemp).gt.0)
                ixf = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
                if (ixf.gt.1) then
                   xfile_num_loc=xfile_num_loc+1
                   xfile_num_tot=xfile_num_tot+1
                end if
                ctemp = adjustl(ctemp(ixf+1:))
             end do
             do ixf=xfile_num_tot-xfile_num_loc+1,xfile_num_tot
                do
                   read(xfile_unit(ixf),'(a)',end=113,err=103,iostat=ierr) cjunk
                   do j=1,ncomm
                      ic=index(cjunk,comment(j))
                      if (ic.gt.0) cjunk(ic:) = ' '
                   end do

                   if (len_trim(cjunk).gt.0) then
                      nrecords=nrecords+1
                      llist(nrecords) = adjustl(cjunk)
                   end if
                end do
113             close(xfile_unit(ixf))
             end do
          end if

       end do
101    close(unit)

    end if

    ! Now read in the tokens from llist
    tlist = ' '

    do i=1,nrecords
       ctemp = llist(i)
       nt=0
       do while(len_trim(ctemp).gt.0)
          ic = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
          if (ic.gt.1) then
             nt=nt+1
             tlist(nt,i) = adjustl(ctemp(:ic-1))
          end if
          ctemp = adjustl(ctemp(ic+1:))
       end do
    end do

    ! Check if any of the "labels" in the input file are unrecognised
    inblock=.false.
    do i=1,nrecords
       ! Check if we are in a block
       if (esdf_reduce(tlist(1,i)).eq.'%block') then
          inblock = .true.
          ! Check if block label is recognised
          ! jd: Elide the check if block name contains '-', to allow blocks
          !     whose name only becomes known at runtime
          if(index(esdf_reduce(tlist(2,i)),'-') .eq. 0) then
             if ((count(esdf_reduce(tlist(2,i)).eq.kw%label).eq.0)) then
                ctemp='Unrecognised block label: "'//trim(esdf_reduce(tlist(2,i)))//'"'
                if (count(ctemp.eq.warns).eq.0) call esdf_error(ctemp)
             end if
          end if
          ! Check if "label" is multiply defined in the input file
          ndef=0
          do j=1,nrecords
             if (esdf_reduce(tlist(2,i)).eq.esdf_reduce(tlist(2,j))) ndef=ndef+1
          end do
          ctemp='Label "'//trim(esdf_reduce(tlist(2,i)))//&
               &'" is multiply defined in the input file. '
          if ((ndef.gt.2).and.(count(ctemp.eq.errors).eq.0))&
               call esdf_error(ctemp)
       end if
       ! Check it is in the list of keywords
       if ((count(esdf_reduce(tlist(1,i)).eq.kw%label).eq.0)&
            .and.(.not.inblock)) then
          ctemp='Unrecognised parameter: "'//trim(esdf_reduce(tlist(1,i)))//'"'
          if (count(ctemp.eq.warns).eq.0) call esdf_error(ctemp)
       end if
       if (.not.inblock) then
          ! Check if "label" is multiply defined in the input file
          ndef=0
          do j=1,nrecords
             if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(tlist(1,j))) ndef=ndef+1
          end do
          ctemp='Label "'//trim(esdf_reduce(tlist(1,i)))//&
               &'" is multiply defined in the input file. '
          if ((ndef.gt.1).and.(count(ctemp.eq.errors).eq.0)) &
               call esdf_error(ctemp)
       end if
       ! Check if we have left a block
       if (esdf_reduce(tlist(1,i)).eq.'%endblock') inblock= .false.

    end do

    return

    ! Error on read while counting records
102 write(6,'(3a,i6)') 'Error in esdf_init: counting records in file "', &
         trim(filename),'" failed with code ',ierr
    error = 1
    return

    ! Error on read while reading records
103 write(6,'(3a,i6)') 'Error in esdf_init: reading records in file "', &
         trim(filename),'" failed with code ',ierr
    error = 1
    return

  end subroutine esdf_init

  !
  ! return the string attached to the "label"
  !
  function esdf_string(label,default,force_no_case_change)

    character(*), intent(in) :: label,default
    character(llength) :: esdf_string
    logical, intent(in), optional :: force_no_case_change

    ! Local
    integer :: i
    character(llength) :: ctemp
    ! aam: For conversion of string to upper case
    integer   :: ipos      ! Loop variable over characters in string
    character :: letter    ! character in string
    logical :: loc_no_case_change
    character(len=256) :: loc_buffer ! JCW
    character(len=6)   :: int_str ! JCW

    ! rab207: added to prevent case change if, e.g. reading in a file name
    if (present(force_no_case_change)) then
      loc_no_case_change = force_no_case_change
    else
      loc_no_case_change = .false.
    endif

    ! Check "label" is defined
    call esdf_lblchk(label,'T')

    ! Set to default
    esdf_string = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          !jme: avoid matches of trim(tlist(2,i)) within tlist(1,i)
          ipos = index(llist(i),trim(tlist(1,i)))+len_trim(tlist(1,i))-1
          ctemp = repeat(' ',ipos) // llist(i)(ipos+1:)
          esdf_string = llist(i)(index(ctemp,trim(tlist(2,i))):)
          exit
       end if

    end do

    ! aam: Convert string to uppercase
    ! aam: Loop over characters in string
    ! rab207: now optional (on by default)
    if (.not.loc_no_case_change) then
      do ipos=1,len(esdf_string)

         letter = esdf_string(ipos:ipos)
         if (letter >= 'a' .and. letter <= 'z') &
              esdf_string(ipos:ipos) = achar(iachar(letter)-32)

      end do
    endif

    ! JCW: Check the string length before dumping and raise an error if string is
    ! JCW: longer than llength
    write(loc_buffer,*) trim(esdf_reduce(label)),':',trim(esdf_string)
    if (len(trim(loc_buffer)) > llength) then
       ! JCW: llength is the length of elements of the dump array
       write(int_str,'(i6)') llength
       call esdf_die("Problem in esdf_string: "//&
            "Cannot write "//trim(loc_buffer)//" to dump array, as length is &
            &greater than llength ("//trim(adjustl(int_str))//"). Please &
            &shorten length of string value.")
    end if

    ! Dump the string used
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',trim(esdf_string)
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_string'
    call esdf_warn(ctemp)
  end function esdf_string

  !
  ! return the integer attached to the "label"
  !
  function esdf_integer(label,default)

    integer, intent(in) :: default
    character(*), intent(in) :: label
    integer :: esdf_integer

    ! Local
    integer :: i
    character(llength) :: ctemp

    ! Check "label" is defined
    call esdf_lblchk(label,'I')

    ! Set to default
    esdf_integer = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          read(tlist(2,i),*,err=100) esdf_integer
          exit
       end if

    end do

    ! Dump the value used

    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',esdf_integer
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return


100 ctemp = 'Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_integer'
    call esdf_die(ctemp)

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_integer'
    call esdf_warn(ctemp)

  end function esdf_integer

  !
  ! return the single precisioned value attached to the "label"
  !
  function esdf_single(label,default)

    real(SP), intent(in) :: default
    character(*), intent(in) :: label
    real(SP) :: esdf_single

    ! Local
    integer :: i
    character(llength) :: ctemp

    ! Check "label" is defined
    call esdf_lblchk(label,'S')

    ! Set to default
    esdf_single = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          read(tlist(2,i),*,err=100) esdf_single
          exit
       end if

    end do

    ! Dump the value used
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',esdf_single
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

100 ctemp = 'Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_single'
    call esdf_die(ctemp)
101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_single'
    call esdf_warn(ctemp)

  end function esdf_single

  !
  ! return the double precisioned value attached to the "label"
  !
  function esdf_double(label,default)

    real(kind=DP), intent(in) :: default
    character(*), intent(in) :: label
    real(kind=DP) :: esdf_double

    ! Local
    integer :: i
    character(llength) :: ctemp

    ! Check "label" is defined
    call esdf_lblchk(label,'D')

    ! Set to default
    esdf_double = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          read(tlist(2,i),*,err=100) esdf_double
          exit
       end if

    end do

    ! Dump the value used
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',esdf_double
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

100 esdf_double = default
    ctemp = 'Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_double'
    call esdf_die(ctemp)
101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_double'
    call esdf_warn(ctemp)

  end function esdf_double

  !
  ! return the double precisioned physical value attached to the "label"
  ! Units converted to "dunit"
  !
  function esdf_physical(label,default,dunit)

    real(kind=DP), intent(in) :: default
    character(*), intent(in) :: label,dunit
    real(kind=DP) :: esdf_physical

    ! Local
    integer :: i
    character(llength) :: ctemp,iunit

    ! Check "label" is defined
    call esdf_lblchk(label,'P')

    ! Set to default
    esdf_physical = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          read(tlist(2,i),*,err=100) esdf_physical
          ! aam_2013-07-24: if unit not present, assume default and issue warning.
          ctemp = tlist(3,i)
          if (len(trim(adjustl(ctemp))).gt.0) then
             ! jd: Changed to '(a)', or otherwise parsing stopped at '/',
             !     so 'ha/bohr' would be interpreted as merely 'ha'.
             read(tlist(3,i),'(a)',err=100) iunit
             esdf_physical = esdf_convfac(iunit,dunit)*esdf_physical
          else
             ctemp = '"'//trim(esdf_reduce(label))//'" (esdf_physical) &
                  &has no units. Assume "'//trim(esdf_reduce(dunit))//'".'
             call esdf_warn(ctemp)
          endif
          exit
       end if

    end do

    ! Dump the value used
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',esdf_physical,' ',trim(dunit)
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

100 esdf_physical = default
    ctemp = 'Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_physical'
    call esdf_die(ctemp)
101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_physical'
    call esdf_warn(ctemp)

  end function esdf_physical

  !
  ! Is the "label" defined in the input file
  !
  function esdf_defined(label)

    character(*), intent(in) :: label
    logical :: esdf_defined

    ! Local
    integer :: i
    character(llength) :: ctemp

    ! Check "label" is defined
    call esdf_lblchk(label,'E')

    ! Set to default
    esdf_defined = .false.

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          esdf_defined = .true.
          exit
       end if

    end do

    ! Dump the value used
    if (esdf_defined) then

       ndmp=ndmp+1
       write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':'
       if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    end if

    return

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_defined'
    call esdf_warn(ctemp)

  end function esdf_defined

  !
  ! Is the "label" defined in the input file
  !
  function esdf_boolean(label,default)

    character(*), intent(in) :: label
    logical, intent(in) :: default
    logical :: esdf_boolean

    ! Local
    integer :: int_buffer(3)
    integer :: i
    character(llength) :: ctemp,positive(3),negative(3)

    data positive /'yes','true','t'/
    data negative /'no','false','f'/

    ! Check "label" is defined
    call esdf_lblchk(label,'L')

    ! Set to default
    esdf_boolean = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          if (len_trim(tlist(2,i)).eq.0) then
             esdf_boolean = .true.
             exit
          end if

          int_buffer =index(positive,esdf_reduce(tlist(2,i)) )

          if ( Sum( int_buffer ) .gt.0 ) then

             esdf_boolean = .true.
             exit
          end if

          int_buffer =index(negative,esdf_reduce(tlist(2,i)))

          if (  Sum( int_buffer ) .gt.0 ) then

             esdf_boolean = .false.
             exit
          end if
          call esdf_die('Unable to parse boolean value')

       end if

    end do

    ! Dump the value used
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),': ',esdf_boolean
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_boolean'
    call esdf_warn(ctemp)

  end function esdf_boolean

  function esdf_block(label,nlines,dont_check)
    ! jd: Modified for the possibility to override the label check --
    !     this is useful when the block name only becomes known at runtime.

    use utils, only: utils_abort, utils_int_to_str

    character(*), intent(in) :: label
    integer, intent(out) :: nlines
    logical, intent(in), optional :: dont_check
    logical :: esdf_block

    ! Local
    integer :: i,j
    character(llength) :: ctemp

    ! Check "label" is defined (unless overridden)
    if(present(dont_check)) then
       if(.not. dont_check) then
          call esdf_lblchk(label,'B')
       end if
    end if

    ctemp ='Block "'//trim(esdf_reduce(label))//'" not closed correctly '

    esdf_block=.false.

    nlines = 0

    do i=1,nrecords
       if ((esdf_reduce(tlist(1,i)).eq.esdf_reduce('%block'))&
            .and.(esdf_reduce(tlist(2,i)).eq.esdf_reduce(label))) then
          esdf_block = .true.
          do while(esdf_reduce(tlist(1,i+nlines+1))&
               .Ne.esdf_reduce('%endblock'))
             nlines=nlines+1
             if (nlines+i.gt.nrecords) call esdf_die(ctemp)
             block_data(nlines)=llist(i+nlines)
          end do

          if (esdf_reduce(tlist(2,i+nlines+1)).Ne.esdf_reduce(label))&
               call esdf_die(ctemp)
          exit
       end if
    end do

    if (.not.esdf_block) return

    ! Dump the block
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101) '%block ',trim(esdf_reduce(label)),': '
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) then
       ndmp=ndmp-1
       return
    end if
    do j=1,nlines
       ndmp=ndmp+1
       ! jd: Handle overflows gracefully.
       if(ndmp > ndump) then
          call utils_abort("Maximum %block size ("//trim(&
               utils_int_to_str(ndump))//" lines) exceeded. &
               &Adjust 'ndump' in esdf_mod.F90, recompile and try again.")
       end if
       dump(ndmp)=block_data(j)
    end do
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101) '%endblock ',trim(esdf_reduce(label)),': '

    return

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_block'
    call esdf_warn(ctemp)

  end function esdf_block

  !
  ! Reduce the string to lower case and remove punctuation
  !
  function esdf_reduce(string)

    character(*), intent(in) :: string
    character(llength) :: esdf_reduce

    ! Local

    integer, parameter :: npunct=4
    integer :: iA,iZ,ishift,ic,i,ln
    character(llength) :: ctemp
    character(1) :: punct(npunct)
    logical :: keep_going

    ! Define the punctuation to be removed
    data punct /'.','_','"',''''/
    ! jd: Removed '-' from punct, as otherwise it breaks unit names, e.g.
    !     (cm^-1 -> cm^1)

    ! Initialise system dependent bounds in collating sequence
    iA = ichar('A');iZ = ichar('Z')
    ishift = ichar('a')-iA

    ! Initialise output
    ln = len_trim(string)



    esdf_reduce(1:ln) = string(1:ln)
    if (ln < llength) esdf_reduce(ln+1:)=' '

    ! Drop all upper case characters to lower case
    do i=1,llength
       ic = ichar(esdf_reduce(i:i))
       if ((ic.ge.iA).and.(ic.le.iZ)) esdf_reduce(i:i) = char(ishift+ic)
    end do

    ! Now remove punctuation
    do i=1,npunct

       keep_going=.true.

       do while(keep_going)
          ic = index(esdf_reduce,punct(i))
          if (ic.gt.0) then
             ctemp = esdf_reduce
             esdf_reduce(ic:)=ctemp(ic+1:)//' '
          else
             keep_going=.false.
          end If
       end do

    end do

    esdf_reduce = trim(adjustl(esdf_reduce))

  end function esdf_reduce

  function esdf_reduce_max(string)

    character(*), intent(in) :: string
    character(esdf_maxstring) :: esdf_reduce_max

    ! Local

    integer, parameter :: npunct=4
    integer :: iA,iZ,ishift,ic,i,ln
    character(esdf_maxstring) :: ctemp
    character(1) :: punct(npunct)
    logical :: keep_going

    ! Define the punctuation to be removed
    data punct /'.','_','"',''''/
    ! jd: Removed '-' from punct, as otherwise it breaks unit names, e.g.
    !     (cm^-1 -> cm^1)

    ! Initialise system dependent bounds in collating sequence
    iA = ichar('A');iZ = ichar('Z')
    ishift = ichar('a')-iA

    ! Initialise output
    ln = len(string)

    esdf_reduce_max(1:ln) = string(1:ln)
    if (ln < esdf_maxstring) esdf_reduce_max(ln+1:)=' '

    ! Drop all upper case characters to lower case
    do i=1,esdf_maxstring
       ic = ichar(esdf_reduce_max(i:i))
       if ((ic.ge.iA).and.(ic.le.iZ)) esdf_reduce_max(i:i) = char(ishift+ic)
    end do

    ! Now remove punctuation
    do i=1,npunct

       keep_going=.true.

       do while(keep_going)
          ic = index(esdf_reduce_max,punct(i))
          if (ic.gt.0) then
             ctemp = esdf_reduce_max
             esdf_reduce_max(ic:)=ctemp(ic+1:)//' '
          else
             keep_going=.false.
          end If
       end do

    end do

    esdf_reduce_max = trim(adjustl(esdf_reduce_max))

  end function esdf_reduce_max

  !
  ! Find the conversion factor between physical units
  !
  function esdf_convfac(from,to)

    use constants, only: CRLF
    use utils, only: utils_assert

    character(*), intent(in) :: from,to
    real(kind=DP) :: esdf_convfac

    ! Local
    integer :: i,ifrom,ito
    character(llength) :: ctemp

    ! Find the index numbers of the from and to units
    ifrom = 0 ; ito = 0
    do i=1, nphys
       if (esdf_reduce(from).eq.phy(i)%n) ifrom = i
       if (esdf_reduce(to).eq.phy(i)%n) ito = i
    end do

    ! Check that the units were recognised
    call utils_assert(ifrom /= 0, &
         'Units not recognised in input file: "'//trim(esdf_reduce(from))//&
         '".'//CRLF//'Also it might be worth checking if the atomic species &
         &used are consistent between %blocks.')
    call utils_assert(ito /= 0, &
         'Units not recognised in program: "'//trim(esdf_reduce(to))//&
         '". Double-check capitalisation.')
    ! Check that from and to are of the same dimensions
    call utils_assert(phy(ifrom)%d == phy(ito)%d, &
         'Unit dimensions do not match: "'//trim(esdf_reduce(from))//&
         '" vs "'//trim(esdf_reduce(to))//'".')

    ! Set the conversion factor
    esdf_convfac = phy(ifrom)%u/phy(ito)%u

  end function esdf_convfac

  !
  ! Find an unused i/o unit
  !
  function esdf_unit(ierr)
    integer, intent(out) :: ierr
    integer :: esdf_unit
    ! Local
    logical :: op
    ierr=0
    do esdf_unit=10,99
       inquire(unit=esdf_unit,opened=op,err=100)
       if (.not.op) return
    end do
    call esdf_warn('Unable to find a free i/o unit using esdf_unit')
    ierr = 1
    return
100 call esdf_die('Error opening files by esdf_unit')
  end function esdf_unit

  !
  ! open an old file
  !
  subroutine esdf_file(unit,filename,ierr)
    character(*), intent(in) :: filename
    integer, intent(out) :: unit,ierr
    logical :: ex
    unit = esdf_unit(ierr)
    if (ierr.gt.0) return
    inquire(file=trim(filename),exist=ex,err=100)
    if (.not.ex) goto 100
    open(unit=unit,file=trim(filename),form='formatted',status='old',err=100)
    return
100 ierr=1
    return
  end subroutine esdf_file

  ! open a new file

  subroutine esdf_newfile(unit,filename,ierr)
    character(*), intent(in) :: filename
    integer, intent(out) :: unit,ierr
    unit = esdf_unit(ierr)
    if (ierr.gt.0) return
    open(unit=unit,file=trim(filename),form='formatted',status='replace',err=100)
    return
100 ierr=1
    return
  end subroutine esdf_newfile

  !
  ! Check that the label is known, and used correctly
  !
  subroutine esdf_lblchk(string,typ)
    character(*), intent(in) :: string
    character(1), intent(in) :: typ
    ! Local
    character(llength) :: ctemp
    character(1) :: tp
    integer :: i
    ! Check if label is recognised
    i=count(esdf_reduce(string).eq.kw%label)
    ctemp = 'Label "'//trim(esdf_reduce(string))//'" not recognised in&
         & keyword list'
    if (i.eq.0) call esdf_die(ctemp)
    ctemp = 'Label "'//trim(esdf_reduce(string))//'" is multiply defined'
    if (i.gt.1) call esdf_die(ctemp)

    ! Type E means we just want to know if something is defined.
    if(typ/='E') then
       ctemp = 'Label "'//trim(esdf_reduce(string))//'" has been used with the wrong type'
       tp = ' '
       i=0
       do while(tp.eq.' ')
          i=i+1
          if (esdf_reduce(string).eq.kw(i)%label) tp=kw(i)%typ(1:1)
       end do
       if (typ.Ne.tp) call esdf_die(ctemp)
    end if
  end subroutine esdf_lblchk


  subroutine esdf_help(helpword,searchword)

    use constants, only: file_maxsize
    use utils, only: utils_abort

    Implicit None

    character(len=*), intent(in) :: helpword
    character(len=*), intent(in) :: searchword

    ! Local
    integer  :: i,indx,indx2,ln
    character(20) :: ctyp,clev,cgroup
    character(file_maxsize) :: title,fm
    character(file_maxsize) :: ctemp
    character(file_maxsize) :: ctemp2
    character(1)  :: cl
    logical       :: found

    call esdf_key_init

    if(esdf_reduce(helpword).eq.'help') then
       if(esdf_reduce(searchword).eq.'search') then
          write(6,'(a)') '<searchword> is empty'
          return
       end if
       write(6,'(a)') 'Usage:'
       write(6,'(a)') 'onetep <inputname>                            : Run files <inputname>.dat'
       write(6,'(a)') '  "    [-s|--search|-help search] <text>      : print list of keywords with <text> match in description'
       write(6,'(a)') '  "    [-h|--help] <keyword>                  : describe specific keyword in kewyord list'
       write(6,'(a)') '  "         "      all                        : print list of all keywords'
       write(6,'(a)') '  "         "      basic                      : print list of basic-level keywords'
       write(6,'(a)') '  "         "      inter                      : print list of intermediate-level keywords'
       write(6,'(a)') '  "         "      expert                     : print list of expert-level keywords'
       write(6,'(a)') '  "         "      dummy                      : print list of dummy keywords'
       return
    end if

    if (esdf_reduce(helpword).eq.'search') then
       if (esdf_reduce(searchword)=="") then
          write(6,'(a)') '<searchword> is empty'
          return
       end if

       found=.false.
       ctemp=esdf_reduce_max(searchword)
       ! Search for useful keywords
       do i=1,numkw
          if ((index(esdf_reduce_max(kw(i)%label),trim(ctemp)).gt.0).Or.&
               (index(esdf_reduce_max(kw(i)%dscrpt),trim(ctemp)).gt.0)) then
             found=.true.
             indx=index(kw(i)%dscrpt,'!*')-1
             if (indx.eq.-1) &
                  call esdf_die('help: keyword description incorrectly formatted')
             title = kw(i)%dscrpt(1:indx)
             ln=len(trim(title))
             if (ln.gt.file_maxsize) call esdf_die('help: keyword title too long')
             write(6,'(a,T33,a)') kw(i)%label,trim(title(3:))
          end If
       end do
       if (.not. found) then
          write(6,'(a)') 'None found'
          return
       else
          return
       end if
    end if

    ! All keywords, short description
    if ('all'.eq.esdf_reduce(helpword)) then
       do i=1,numkw
          if (len_trim(kw(i)%label).gt.0) then
             indx=index(kw(i)%dscrpt,'!*')-1
             if (indx.eq.-1) &
                  call esdf_die('help: keyword description incorrectly formatted')
             title = kw(i)%dscrpt(1:indx)
             ln=len(trim(title))
             if (ln.gt.file_maxsize) call esdf_die('help: keyword title too long')
             write(6,'(a,T33,a)') kw(i)%label,trim(title(3:))
          end If
       end do
       return
    end If

    ! All specific levels of keywords
    if (any((/'basic ','inter ','expert','dummy '/).eq.esdf_reduce(helpword))) then
       cl = ' ' ! qoh: Initialise to prevent compiler warning
       select case(esdf_reduce(helpword))
       case('basic')  ; cl = 'B'
       case('inter')  ; cl = 'I'
       case('expert') ; cl = 'E'
       case('dummy')  ; cl = 'D'
       end select

       do i=1,numkw
          if (kw(i)%typ(3:3).eq.cl) then
             indx=index(kw(i)%dscrpt,'!*')-1
             if (indx.eq.-1) &
                  call esdf_die('help: keyword description incorrectly formatted')
             title = kw(i)%dscrpt(1:indx)
             ln=len(trim(title))
             if (ln.gt.file_maxsize) call esdf_die('help: keyword title too long')

             write(6,'(a,T33,a)') kw(i)%label,trim(title(3:))

          end If
       end do
       return
    end If

    ! "reduce" the keyword list for comparison
    do i = 1,numkw
       ctemp = kw(i)%label
       ! qoh: Explicit character trunctation to avoid compiler warning
       ctemp2 = esdf_reduce_max(ctemp)
       kw(i)%label = ctemp2
    end do

    ! More information about a specific keyword
    if (.not.any(kw%label.eq.esdf_reduce(helpword))) then
       write (6,'(a)') 'Keyword not recognised '
       write(6,'(a)') 'It might not exist or you might have spelled it wrong'
       write(6,'(a)') 'Hint: try -help search <keyword> or -help all'
       return
    end if
    if (count(kw%label.eq.esdf_reduce(helpword)).gt.1) &
         call esdf_die('help: keyword entry duplicated')
    found = .false.
    do i=1,numkw
       if (esdf_reduce_max(kw(i)%label).eq.esdf_reduce(helpword)) then
          found = .true.
          indx=index(kw(i)%dscrpt,'!*')+1
          if (indx.eq.1) &
               call esdf_die('help: keyword description incorrectly formatted')
          title = kw(i)%dscrpt(1:indx)
          ln=len(trim(title))
          if (ln.gt.file_maxsize) &
               call esdf_die('help: keyword title too long')
          if (ln.le.9) write(fm,'("(",i2,"x,a",i1,")")') 40-ln/2,ln
          if (ln.gt.9) write(fm,'("(",i0,"x,a",i0,")")') max(40-ln/2,1),ln

          write(6,'(a)')
          write(6,fm) trim(title)
          write(6,'(a)')
          select case(kw(i)%typ(1:1))
          case('I') ; ctyp ='integer'
          case('S') ; ctyp ='Single Precision'
          case('D') ; ctyp ='Double Precision'
          case('P') ; ctyp ='Physical'
          case('T') ; ctyp ='String'
          case('E') ; ctyp ='Defined'
          case('B') ; ctyp ='Block'
          case('L') ; ctyp ='Boolean'
          end select
          select case(kw(i)%typ(3:3))
          case('B') ; clev ='Basic'
          case('I') ; clev ='Intermediate'
          case('E') ; clev ='Expert'
          case('D') ; clev ='Dummy'
          end select
          if(kw(i)%kwgrp.ne."") then
             ln=len(kw(i)%kwgrp)
             if (ln.eq.1) call esdf_die('help: keyword description incorrectly formatted')
             cgroup = kw(i)%kwgrp(1:ln)
             ln=len(trim(title))
             write(6,'(a,a,18x,a,a,18x,a,a)') 'type:',trim(ctyp),'level:',trim(clev),'group:',trim(cgroup)
             indx=indx+1
          else
             write(fm,'(a,i0,a)') '("type: ",a,',&
                  78-(7+len_trim(clev))-(6+len_trim(ctyp)),'x," Level: ",a)'
             write(ctemp,fm) trim(ctyp),trim(clev)
             write(6,'(a)') trim(ctemp)
             write(6,*)
             indx=indx+1
             ln = len(trim(kw(i)%dscrpt))
             do while (indx.lt.ln)
                ctemp = kw(i)%dscrpt(indx:Min(indx+file_maxsize,ln))
                indx2=index(ctemp,' ',back=.true.)
                indx=indx+len(ctemp(:indx2))
             end do
          end if

       end If
    end do
    if( .not. found) then
       write(6,'(a)') 'None found, sorry.'
       return
    end if

  end subroutine esdf_help

  !
  ! Stop execution due to an error cause by esdf
  !
  subroutine esdf_die(string)
    use utils, only: utils_abort
    character(*), intent(in) :: string
    call utils_abort('ESDF error: '//trim(string))
  end subroutine esdf_die

  !
  ! Warning due to an error cause by esdf
  !
  subroutine esdf_warn(string)
    character(*), intent(in) :: string
    nwarns=nwarns+1
    warns(nwarns) = string
  end subroutine esdf_warn

  !
  ! Dump the warnings to screen
  !
  subroutine esdf_warnout

    use constants, only: stdout

    integer :: i
    character(len=80) :: ctemp

    do i=1,nwarns
       ! aam_2013-07-24: some prettification
       ! write(6,*)'i=',i
       ! write(6,*)'warns(i)=',warns(i)
       write(ctemp,*) i
       write(stdout,*) 'ESDF WARNING(',trim(adjustl(ctemp)),'): '//trim(warns(i))
    end do


  end subroutine esdf_warnout

  !
  ! Error by esdf
  !
  subroutine esdf_error(string)
    character(*), intent(in) :: string
    nerrors=nerrors+1
    errors(nerrors) = string
  end subroutine esdf_error

  !
  ! Dump the errors to screen
  !
  subroutine esdf_errorout

    use constants, only: stdout
    use utils, only: utils_assert

    integer :: i
    character(len=80) :: ctemp

    ! return if no errors
    if (nerrors == 0) return

    write(stdout,'(a)') repeat('~',80)
    write(stdout,'(a)') 'Problems found with input file: '
    write(stdout,'(a)') repeat('~',80)

    do i=1,nerrors
       write(ctemp,*) i
       write(stdout,'(1x,3a)') &
            'ESDF ERROR(',trim(adjustl(ctemp)),'): '//trim(errors(i))
    end do

    call utils_assert(nerrors == 0, 'Problems found with input file: &
         &see output file for details.')


  end subroutine esdf_errorout

  !
  ! Deallocate the data arrays --- call this before re-initialising
  !
  subroutine esdf_close

    use constants, only: stdout
    use utils, only: utils_dealloc_check
    implicit none

    integer :: ierr

    deallocate(llist,stat=ierr)
    call utils_dealloc_check('esdf_[close]','llist',ierr)
    deallocate(block_data,stat=ierr)
    call utils_dealloc_check('esdf_[close]','block_data',ierr)
    deallocate(tlist,stat=ierr)
    call utils_dealloc_check('esdf_[close]','tlist',ierr)
    deallocate(warns,stat=ierr)
    call utils_dealloc_check('esdf_[close]','warns',ierr)
    deallocate(errors,stat=ierr)
    call utils_dealloc_check('esdf_[close]','errors',ierr)
    deallocate(dump,stat=ierr)
    call utils_dealloc_check('esdf_[close]','dump',ierr)

  end subroutine esdf_close

  ! smmd : 12/07/2011
  ! Split an input line into pieces according to divide symbols
  subroutine esdf_line_divide(snum,sline,line)

    ! local
    integer, intent(inout) :: snum
    character(llength), intent(in) :: line
    character(llength), intent(inout) :: sline(snum)
    ! jd: coding up tab as a constant rather than a literal prevents nasty
    !     breakage due to editors replacing tabs by spaces
    character(1), parameter :: tab = achar(9)
    character(1) :: divide(4)

    data divide /' ','=',':',tab/

    integer :: ixf, ntmp
    character(llength) :: cjunk

    ntmp = 0
    cjunk = adjustl(line)

    split_loop : do while(len_trim(cjunk).gt.0 .and. ntmp.lt.snum)
       ixf = minval(index(cjunk,divide),mask=index(cjunk,divide)>0)
       if (ixf .gt. 1) then
          ntmp = ntmp + 1
          sline(ntmp) = adjustl(cjunk(:ixf-1))
       end if
       cjunk = adjustl(cjunk(ixf+1:))
    end do split_loop

    snum = ntmp

  end subroutine esdf_line_divide


  ! cks : 11/7/2001
  ! Dump to the standard output an input file which contains all set variables
  ! including defaults
  subroutine esdf_stdout_dump( )

    ! Local Variables
    integer :: i,j,indx, bpos
    character(llength) :: cjunk, keyword
    character(128)     :: keyword_buffer1, keyword_buffer2, keyword_buffer3
    logical :: inblock

    ! Find maximum index at which a ':' separator occurs
    indx = 0
    inblock = .false.
    do i=1,ndmp
       if (.not.inblock.and.(index(esdf_reduce(dump(i)),'%block')>0)) &
            inblock =.true.
       if (inblock.and.(index(esdf_reduce(dump(i)),'%endblock')>0)) &
            inblock =.false.
       if (.not.inblock) indx = max(indx,index(dump(i),':'))
    end do

    ! Re-format the dump lines so that the ':' characters line up
    inblock = .false.
    do i=1,ndmp
       if (inblock) then
          j = 0
       else
          j = index(dump(i),':')
       end if
       if (.not.inblock.and.(index(esdf_reduce(dump(i)),'%block')>0)) &
            inblock =.true.
       if (inblock.and.(index(esdf_reduce(dump(i)),'%endblock')>0)) then
          j = index(dump(i),':')
          inblock =.false.
       end if
       if (j.gt.0) then
          cjunk = dump(i)(1:j-1)
          keyword_buffer1 = trim(adjustl(dump(i)(j+1:)))
          dump(i) = trim(cjunk)//repeat(' ',indx-j+1)//': '// keyword_buffer1
       end if
    end do

    ! Find last incidence of '#' character on any dump line
    indx = maxval(index(dump(1:ndmp),'#',back=.true.))

    ! Re-format the dump lines so that the '#' characters line up
    do i=1,ndmp
       j=index(dump(i),'#',back=.true.)
       if (j.gt.0) then
          dump(i)=dump(i)(1:j-1)//repeat(' ',indx-j+1)//'#'
       end if
    end do

    ! Now add the keyword descriptions for each variable
    do i=1,ndmp
       if (inblock) then
          j = 0
       else
          j = index(dump(i),':')
       end if
       if (.not.inblock.and.(index(esdf_reduce(dump(i)),'%block')>0)) &
            inblock =.true.
       if (inblock.and.(index(esdf_reduce(dump(i)),'%endblock')>0)) then
          j = index(dump(i),':')
          inblock =.false.
       end if
       if (j.gt.0) then
          cjunk = dump(i)(1:j-1)
          ! qoh,jd: handle "%block" and "%endblock" manually
          bpos = index(cjunk,"%block")
          if ( bpos > 0 .and. len(cjunk) > bpos+6 ) cjunk = dump(i)(bpos+6:j-1)
          bpos = index(cjunk,"%endblock")
          if ( bpos > 0 .and. len(cjunk) > bpos+9 ) cjunk = dump(i)(bpos+9:j-1)
          keyword = trim(adjustl(cjunk))

          do j=1,numkw
             if (keyword == trim(adjustl(kw(j)%label))) exit
          ! jd: Originally was:
          !  if (index(cjunk,trim(kw(j)%label)).gt.0) exit
          !     ... but this gave false matches if one keyword was
          !         a substring of another
          end do

          if(j <= numkw) then
             select case(kw(j)%typ(1:1))
             case('I') ; cjunk ='Integer'
             case('S') ; cjunk ='Single Precision'
             case('D') ; cjunk ='Double Precision'
             case('P') ; cjunk ='Physical'
             case('T') ; cjunk ='String'
             case('E') ; cjunk ='Defined'
             case('B') ; cjunk ='Block'
             case('L') ; cjunk ='Boolean'
             end select
             indx=index(kw(j)%dscrpt,'!*')

             ! cks: HACK in order to make it run on Franklin
             keyword_buffer1 =trim(dump(i))
             keyword_buffer2 =trim(kw(j)%dscrpt(1:indx-1))
             keyword_buffer3 =trim(adjustl(cjunk))
             dump(i) = trim(keyword_buffer1)//' '//trim(keyword_buffer2)//' ('//trim(keyword_buffer3)//')'
          else
             dump(i) = trim(dump(i))//' !* [special block]'
          end if
       end if
    end do

    do i=1,ndmp
       write(6,'(a)') adjustl(dump(i))
    end do

  end subroutine esdf_stdout_dump

  !=====================================================

  subroutine esdf_dump_input_file
    !====================================================!
    ! Routine to dump the input file to the output file. !
    ! Position data is suppressed, and instead a hash of !
    ! the missing lines (based directly on characters) is!
    ! generated and printed.                             !
    ! Robert Bell, 23/07/14                              !
    !====================================================!

    use constants, only: stdout, LONG
    use utils, only: utils_fnv1_hash, utils_banner

    implicit none

    ! internal
    integer :: iline
    logical :: do_print
    integer(kind=LONG) :: hash
    character(llength) :: filename

    write(stdout,'(//a)') repeat('-',80)
    write(stdout,'(a)') utils_banner('-', 'INPUT FILE')
    write(stdout,'(a/)') repeat('-',80)

    do_print = .true.

    do iline = 1, nrecords

       ! print the hash of position_* data instead of whole block
#ifndef ACCELRYS
       ! jd: Accelrys asked not to print out the hash in their version
       call internal_hash_block('positions')
#endif
       ! ...put other block exceptions here...


       ! either print line, or generate checksum hash for missing lines
       if (do_print) then
          write(stdout,'(a)') llist(iline)
       else
          hash = utils_fnv1_hash(llist(iline),hash)
       endif


       ! if xyz_file input, print hash of entire file
#ifndef ACCELRYS
       ! jd: Accelrys asked not to print out the hash in their version
       call internal_hash_tag_file('xyz_file')
#endif
       ! ...put other input files here...


    enddo

    write(stdout,'(/a)') repeat('-',80)
    write(stdout,'(a)') utils_banner('-', 'END INPUT FILE')
    write(stdout,'(a//)') repeat('-',80)

    contains

    subroutine internal_hash_block(block_tag)
      !===============================================!
      ! data between lines containing start_tag and   !
      ! end_tag (exclusive) are hashed, not printed.  !
      !===============================================!

      implicit none

      ! arguments
      character(len=*), intent(in) :: block_tag

      ! if block start
      if (index(tlist(1,iline),'%block') > 0 .and. &
          index(tlist(2,iline),block_tag) > 0) then
         do_print = .false.
         write(stdout,'(a)') llist(iline)
         hash = utils_fnv1_hash('') ! init hash
      endif

      ! if block end
      if (index(tlist(1,iline),'%endblock') > 0 .and. &
          index(tlist(2,iline),block_tag) > 0) then
         do_print = .true.
         write(stdout,'(5x,a,i26.24)') 'Hash = ', hash
      endif

    end subroutine internal_hash_block

    !===============================================!

    subroutine internal_hash_tag_file(tag)
    !===============================================!
    ! the second entry on the line containing tag is!
    ! taken to be file input data. The hash of the  !
    ! contents of this file is printed.             !
    !===============================================!

       use utils, only: utils_fnv1_hash_file

       implicit none

       ! arguments
       character(len=*), intent(in) :: tag

       if (index(tlist(1,iline),tag) > 0) then
          filename = tlist(2,iline)
          hash = utils_fnv1_hash_file(trim(filename))
          write(stdout,'(5x,a,i26.24)') 'Hash = ', hash
       endif

    end subroutine internal_hash_tag_file

  end subroutine esdf_dump_input_file



end module esdf

