! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by Chris-Kriton Skylaris
!
!   January 2001
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module ion

  use constants, only: DP, max_rundat_swri_block_entries
  use geometry, only: POINT

  implicit none

  private

  ! ndmh: Define type for NGWF radial functions
  type RADIAL_NGWF_TYPE

     ! ndmh: number of shells
     integer :: nshells

     ! ndmh: total number of functions
     integer :: nfuncs

     ! ndmh: number of points of radial grid
     integer :: npts

     ! ndmh: angular momentum of each shell
     integer, allocatable, dimension(:) :: angmom

     ! ndmh: cutoff radius for each shell
     real(kind=DP), allocatable, dimension(:) :: rc

     ! ndmh: radial grid positions
     real(kind=DP), allocatable, dimension(:) :: rad

     ! ndmh: radial parts of NGWF on real space radial grid
     real(kind=DP), allocatable, dimension(:,:) :: func_real

  end type RADIAL_NGWF_TYPE

  ! ndmh: Define type for radial densities
  type RADIAL_DENSITY_TYPE

    ! ndmh: flag to show whether each initial guess density exists
    logical :: present

    ! ndmh: number of points of radial grid
    integer :: npts

    ! ndmh: number of partial waves (if PAW)
    integer :: npwtot

    ! ndmh: pspecies number
    integer :: ipsp

    ! ndmh: radial grid positions
    real(kind=DP), allocatable, dimension(:) :: rad

    ! ndmh: density on radial grid
    real(kind=DP), allocatable, dimension(:,:) :: den

    ! ndmh: initial projector density kernel corresponding to aug_den
    real(kind=DP), allocatable, dimension(:,:,:) :: rhoij0

    ! ndmh: augmentation charge on radial grid (if present)
    real(kind=DP), allocatable, dimension(:,:) :: aug_den

  end type RADIAL_DENSITY_TYPE

  ! ndmh: please note I am aware that SPECIE is not the singular of SPECIES!
  type SPECIE
    ! cks: chemical symbol of atom
    character(len=2) :: symbol
    ! ndmh: species id string
    character(len=4) :: species_id
    ! cks: atomic number
    integer          :: atomic_number
    ! cks: charge of the ionic core of a particular element.
    real(kind=DP)          :: ion_charge
    ! cks: ngwf region radius for atom in a.u.
    ! cks: note that when a halo is defined this is the
    ! cks: sum of the ngwf radius and the halo length
    real(kind=dp) :: radius
    ! cks: conduction ngwf region radius for atom in a.u.
    real(kind=dp) :: radius_cond
    ! cks: name of pseudopotential (file) on atom
    character(len=64) :: pseudo_name
    ! ndmh: names of NGWF sets for valence NGWFs, conduction NGWFs, and
    ! ndmh: auxiliary basis
    character(len=128) :: ngwf_set
    character(len=128) :: cond_ngwf_set
    character(len=128) :: aux_ngwf_set
    ! ndmh: name of core wavefunction file for this atom
    character(len=64) :: core_wf_name
    ! ndmh: number of NGWF functions on atom
    integer :: nfunctions
    ! ndmh: number of conduction NGWFs, radius of conduction NGWFs
    integer :: nfunctions_cond
    ! ndmh: number of auxiliary NGWFs
    integer :: nfunctions_aux
    ! rc2013: embedding region iterator
    integer          :: region
    ! rc2013: NGWF constraint information
    ! This may not be necessary
    character(len=10) :: ngwf_constraint
    ! jd: electrolyte cavity radial steric cutoff (if in use, otherwise 0D0)
    real(kind=DP) :: hc_steric_cutoff
    ! ab: Solvent radius
    real(kind=DP) :: solvent_radius

  end type SPECIE

  ! cks: this structure defines quantities which are all associated with
  !      a particular atom.
  type ELEMENT
     ! cks: cartesian coordinates of atom in a.u.
     type(POINT)      :: centre
     ! cks: chemical symbol of atom
     character(len=2) :: symbol
     ! cks: atomic number
     integer          :: atomic_number
     ! cks: number of NGWF functions on atom
     integer          :: nfunctions
     !gibo: number of angular momentum shells for NGWFs of this atom
     integer              :: num_ang_mom_shells
     !gibo: ang. mom. of each shell in norb
     integer              :: orb_ang_mom(10)
     ! ndmh: number of conduction NGWFs, radius of conduction NGWFs
     integer :: nfunctions_cond
     real(kind=DP) :: radius_cond
     ! ndmh: number of auxiliary NGWFs
     integer :: nfunctions_aux
     ! qoh: "species number" to which element belongs (in %block species)
     integer          :: global_species_number
     ! rc2013: "local species number" in this subsystem
     integer          :: species_number
     ! qoh: "pseudopotential species number" of element
     integer          :: pspecies_number
     ! ndmh: species id string
     character(len=4) :: species_id
     ! cks: charge of the ionic core of a particular element.
     real(kind=DP)          :: ion_charge
     ! cks: number of non-local projectors on atom
     integer :: nprojectors
     ! ndmh: number of PAW partial waves on atom
     integer :: npawpws
     ! ndmh: number of PAW core orbitals on atom
     integer :: ncorewfs
     ! cks: NGWF region radius for atom in a.u.
     ! cks: Note that when a halo is defined this is the
     ! cks: sum of the NGWF radius and the halo length
     real(kind=DP) :: radius
     ! cks: max core radius for any projector of this atom in a.u.
     real(kind=DP) :: max_core_radius
     ! ndmh: max core wavefunction radius for any core wf of this atom in a.u.
     real(kind=DP) :: max_core_wf_radius
     ! aam: ionic constraint information
     character(len=5) :: ion_constraint_type
     real(kind=dp)    :: ion_constraint(3)
     ! cks: selective NGWF plotting information
     logical          :: ngwf_plot
     ! vv: selective local dipole calculation information
     logical          :: loc_dipole
     ! jd: .true. if participates in a given swri
     logical          :: in_swri(max_rundat_swri_block_entries)
     ! jd: Thole point-dipole polarisability (only used in pol_emb context)
     real(kind=DP)    :: thole_polarisability
     ! smmd: ionic velocity
     real(kind=DP)    :: ion_velocity(3)
     ! smmd: group of atoms to which element belongs
     integer          :: group_id
     ! ab: solvent_radius
     real(kind=DP) :: solvent_radius
     ! ja531 -> number of SWs in pDOS
     integer :: nfunctions_pdos

     ! rc2013: NGWF constraint information
     character(len=10) :: ngwf_constraint
     ! rc2013: embedding region iterator
     integer          :: region
     ! rc2013: global atom id
     integer          :: global_atom_number

     ! gab: Soft sphere radius for IS cavity
     real(kind=DP) :: soft_sphere_radius

     ! ab: coordination number (used in dftb)
     real(kind=DP) :: dftb_coord

     ! ab: total charge (used in dftb)
     real(kind=DP) :: dftb_charge

  end type ELEMENT

  ! cks: public type definitions
  public :: ELEMENT
  public :: SPECIE
  public :: RADIAL_NGWF_TYPE
  public :: RADIAL_DENSITY_TYPE

end module ion
