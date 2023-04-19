!================================================================!
!                                                                !
!  Density Functional Tight Binding parameter definition module  !
!  ------------------------------------------------------------  !
!                                                                !
! This module defines the datastructures holding the parameters  !
! of the DFTB method. So far, we are focusing on xTB methods by  !
! the Grimme group. The datastructures have been isolated from   !
! the main DFTB module so that they can be included in           !
! model_type without causing a circular dependency.              !
!----------------------------------------------------------------!
! Work started in 2021/06. Written by Jacek Dziedzic and Arihant !
! Bhandari under the supervision of Chris-Kriton Skylaris.       !
!================================================================!
! References:                                                    !
! [1] P. Pracht, E. Caldeweyher, S. Ehlert, and S. Grimme,       !
!     A robust non-self-consistent tight-binding quantum         !
!     chemistry method for large molecules.                      !
!     ChemRxiv. 2019: 1-19.                                      !
!     https://doi.org/10.26434/chemrxiv.8326202.                 !
! [2] S. Grimme, J. Antony, S. Ehrlich, and H. Krieg,            !
!     A consistent and accurate ab initio parametrization of     !
!     density functional dispersion correction (DFT-D) for the   !
!     94 elements H-Pu.                                          !
!     J. Chem. Phys. 2010 (132), 154104.                         !
!================================================================!


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

module dftb_parameters

  use constants, only: DP, garbage_int, garbage_real

  implicit none

  private

  public :: DFTB_PARAMS

  ! ===========================================================================
  ! jd: Module-wide parameters
  ! ===========================================================================

  ! ab: Number of species in the gfn parameter file
  integer, public :: dftb_num_species = garbage_int

  ! jd: Currently only Z=1..86 are defined, but we can be future-proof up till
  !     row 8.
  integer, parameter, public :: dftb_max_species = 150

  ! ab: maximum supported angular momentum in DFTB
  integer, parameter, public :: dftb_max_angmom = 2

  ! jd: Maximum string length for strings used in this module
  integer, parameter, public :: dftb_strlen = 200

  ! ===========================================================================
  ! jd: Private datastructures
  ! ===========================================================================

  ! ab: STO-nG datastructure
  type sto_ng
     character(len=2) :: ao = "0z" ! ab: name of the ao
     integer :: nprim    ! ab: number of primitives
     real(kind=DP), allocatable :: alpha(:) ! ab: exponents
     real(kind=DP), allocatable :: coeff(:) ! ab: contraction coefficients
  end type

  ! ============================================
  ! jd: Per-species parameters for GFN{0,1,2}.
  !     Some fields are unused in each method.
  ! ============================================
  type DFTB_PER_SPECIES_PARAMS

     ! jd: Common to GFN0, GFN1, GFN2
     character(len=dftb_strlen) :: ao = "*undefined*"
     real(kind=DP)     :: lev(3) = [ garbage_real, garbage_real, garbage_real ]
     real(kind=DP)     :: exp(3) = [ garbage_real, garbage_real, garbage_real ]
     real(kind=DP)     :: gam = garbage_real
     real(kind=DP)     :: polys = garbage_real
     real(kind=DP)     :: polyp = garbage_real
     real(kind=DP)     :: polyd = garbage_real
     real(kind=DP)     :: repa = garbage_real ! jd: aka 'alpha', e.g. [1] (3)
     real(kind=DP)     :: repb = garbage_real ! jd: aka 'Z_eff', e.g. [1] (3)

     ! jd: GFN0 only
     real(kind=DP)     :: alpg = garbage_real
     real(kind=DP)     :: en = garbage_real
     real(kind=DP)     :: kappa = garbage_real
     real(kind=DP)     :: kqat2 = garbage_real
     real(kind=DP)     :: kqs = garbage_real
     real(kind=DP)     :: kqp = garbage_real
     real(kind=DP)     :: kqd = garbage_real
     real(kind=DP)     :: xi = garbage_real

     ! jd: GFN1 only
     real(kind=DP)     :: cxb = garbage_real

     ! jd: GFN2 only
     real(kind=DP)     :: dpol = garbage_real
     real(kind=DP)     :: qpol = garbage_real

     ! jd: GFN0 and GFN2 only
     real(kind=DP)     :: kcns = garbage_real
     real(kind=DP)     :: kcnp = garbage_real
     real(kind=DP)     :: kcnd = garbage_real

     ! jd: GFN1 and GFN2 only
     real(kind=DP)     :: gam3 = garbage_real

     ! jd: Hardcoded parameters, in ONETEP provided via a separate file
     real(kind=DP)     :: elecneg = garbage_real
     real(kind=DP)     :: r0 = garbage_real
     real(kind=DP)     :: cnfactor = garbage_real
     real(kind=DP)     :: covrad = garbage_real
     real(kind=DP)     :: atomicrad = garbage_real
     real(kind=DP)     :: paulingen = garbage_real

     ! ab: inferred parameters
     integer, allocatable :: ang_mom(:) ! ab: angular momentum for each function
     integer, allocatable :: ngwf_type(:) ! ab: shell type for each function
     integer              :: nshells

     ! ab: STO-nG information
     integer, allocatable :: ang(:) ! ab: angular momentum, l for each shell
     type(sto_ng), allocatable :: sto(:) ! ab: STO-nG parameters for each shell

  end type DFTB_PER_SPECIES_PARAMS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ================================================
  ! jd: Per-species-pair parameters, for GFN1 only.
  ! ================================================
  type DFTB_PER_SPECIES_PAIR_PARAMS

     real(kind=DP) :: x = garbage_real

  end type DFTB_PER_SPECIES_PAIR_PARAMS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ===========================================================================
  ! jd: Public datastructures
  !     DFTB_PARAMS must be public so that an instance can reside in e&f.
  ! ===========================================================================

  ! ============================================
  ! jd: Parameter datastructure for GFN{0,1,2}.
  !     Some fields are unused in each method.
  ! ============================================
  type DFTB_PARAMS

     ! ========================================================================
     ! Global parameters
     ! ========================================================================

     ! jd: Identifiers of the parameter set that are not themselves parameters.
     character(len=dftb_strlen) :: level
     character(len=dftb_strlen) :: name
     character(len=dftb_strlen) :: doi

     ! jd: Common to GFN0, GFN1, GFN2
     real(kind=DP) :: ks = garbage_real
     real(kind=DP) :: kp = garbage_real
     real(kind=DP) :: kd = garbage_real
     real(kind=DP) :: kdiff = garbage_real
     real(kind=DP) :: ipeashift = garbage_real
     real(kind=DP) :: a1 = garbage_real
     real(kind=DP) :: a2 = garbage_real
     real(kind=DP) :: s8 = garbage_real
     real(kind=DP) :: s9 = garbage_real
     real(kind=DP) :: kexp = garbage_real
     real(kind=DP) :: kexplight = garbage_real

     ! jd: GFN0 only
     real(kind=DP) :: ens = garbage_real
     real(kind=DP) :: enp = garbage_real
     real(kind=DP) :: end = garbage_real
     real(kind=DP) :: enscale4 = garbage_real
     real(kind=DP) :: srbshift = garbage_real
     real(kind=DP) :: srbpre = garbage_real
     real(kind=DP) :: srbexp = garbage_real
     real(kind=DP) :: srbken = garbage_real
     real(kind=DP) :: renscale = garbage_real

     ! jd: GFN1 only
     real(kind=DP) :: ksp = garbage_real
     real(kind=DP) :: cns = garbage_real
     real(kind=DP) :: cnp = garbage_real
     real(kind=DP) :: cnd1 = garbage_real
     real(kind=DP) :: cnd2 = garbage_real
     real(kind=DP) :: xbdamp = garbage_real
     real(kind=DP) :: xbrad = garbage_real

     ! jd: GFN2 only
     real(kind=DP) :: ksd = garbage_real
     real(kind=DP) :: kpd = garbage_real
     real(kind=DP) :: gam3s = garbage_real
     real(kind=DP) :: gam3p = garbage_real
     real(kind=DP) :: gam3d1 = garbage_real
     real(kind=DP) :: gam3d2 = garbage_real
     real(kind=DP) :: aesshift = garbage_real
     real(kind=DP) :: aesexp = garbage_real
     real(kind=DP) :: aesrmax = garbage_real
     real(kind=DP) :: aesdmp3 = garbage_real
     real(kind=DP) :: aesdmp5 = garbage_real

     ! jd: GFN1 and GFN2 only
     real(kind=DP) :: enscale = garbage_real
     real(kind=DP) :: alphaj = garbage_real

     ! jd: Hardcoded parameters, in ONETEP provided via a separate file
     real(kind=DP) :: p(4,2) = garbage_real

     ! ========================================================================
     ! Per-species parameters
     ! ========================================================================

     type(DFTB_PER_SPECIES_PARAMS) :: z(0:dftb_max_species)
     ! jd: NB: We need index zero here  ^ to be able to safely pass
     !     dftb_par%z(cur_z) to the parser before the "$Z=" tag is ever hit

     ! ========================================================================
     ! Per-species-pair parameters
     ! ========================================================================

     type(DFTB_PER_SPECIES_PAIR_PARAMS) :: zz(dftb_max_species,dftb_max_species)

  end type DFTB_PARAMS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module dftb_parameters
