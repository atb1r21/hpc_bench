! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=====================================================================!
!                                                                     !
!                   Pseudo-atomic solver module                       !
!                                                                     !
! This module contains routines associated with the solution of       !
! the Kohn-Sham equations for an isolated atom subject to sphere      !
! boundary conditions, using a spherical Bessel basis.                !
!---------------------------------------------------------------------!
! This module was created by Nicholas Hine in September 2010.         !
!=====================================================================!

module atom

  use constants, only: DP, stdout, paw_en_size
  use paw, only: PAW_SPECIES
  use pseudopotentials, only: PSEUDO_SPECIES

  implicit none

  private

  ! Type describing a spherical Bessel basis and an associated real-space grid
  type, public :: BASIS_TYPE

     ! number of points in regular real space grid
     integer :: npts

     ! real space cutoff
     real(kind=DP) :: rmax
     real(kind=DP),allocatable :: rmax_l(:)

     ! regular real space grid spacing
     real(kind=DP) :: dr

     ! grid positions of regular real space grid
     real(kind=DP),allocatable :: rad(:)

     ! maximum energy of spherical wave expansion
     real(kind=DP) :: cutoff_energy

     ! maximum angular momentum of all channels
     integer :: lmax

     ! number of spherical waves
     integer,allocatable :: nsws(:)

     ! number of basis functions
     integer,allocatable :: nbfn(:)

     ! maximum number of spherical waves in any channel
     integer :: nsws_max

     ! number of derivatives to fix to zero at sphere boundary
     ! (currently only supports nfixderiv=0 and nfixderiv=1)
     integer :: nfixderiv

     ! maximum number of basis functions in any channel
     integer :: nbfn_max

     ! spherical Bessel function q values
     real(kind=DP),allocatable :: qb(:,:)

     ! spherical Bessel functions on real grid
     real(kind=DP),allocatable :: bessfn(:,:,:)

     ! norm of Bessel functions
     real(kind=DP),allocatable :: bessnorm(:,:)

     ! coefficients of basis functions in terms of Bessel functions
     real(kind=DP),allocatable :: bcoef(:,:,:)

     ! basis functions on real grid (linear combinations of Bessel fns)
     real(kind=DP),allocatable :: bfn(:,:,:)

  end type

  ! Type describing a pseudo-atom
  type, public :: ATOM_TYPE

     ! atomic number of the atom
     real(kind=DP) :: atomic_number

     ! charge of pseudised core of atom
     real(kind=DP) :: ion_charge

     ! total valence charge
     real(kind=DP) :: val_charge

     ! pspecies number for paw or norm_conserv_pseudo
     integer :: pspecies_number

     ! the valence electron density
     real(kind=DP),allocatable :: den(:)

     ! the compensation density and its integral
     real(kind=DP),allocatable :: comp_den(:)
     real(kind=DP) :: comp_den_int

     ! the pseudo core density
     real(kind=DP),allocatable :: core_den(:)

     ! the local pseudopotential
     real(kind=DP),allocatable :: locpspot(:)

     ! the confining potential
     real(kind=DP),allocatable :: confpot(:,:)

     ! the total effective potential
     real(kind=DP),allocatable :: effpot(:)

     ! the xc potential
     real(kind=DP),allocatable :: xcpot(:)

     ! workspace array
     real(kind=DP),allocatable :: work(:)

     ! parameters for confining potential
     real(kind=DP) :: confpot_scale
     real(kind=DP),allocatable :: confpot_width(:)

     ! the pseudopotential projectors
     integer :: nshells, npwtot, proj_lmax, log_npts_max
     integer, allocatable :: proj_ang_mom(:)
     real(kind=DP), allocatable :: proj_q(:,:)
     real(kind=DP), allocatable :: dij0(:,:)
     real(kind=DP), allocatable :: dij(:,:)
     real(kind=DP), allocatable :: dij_xc(:,:)
     real(kind=DP), allocatable :: dij_hart(:,:)
     real(kind=DP), allocatable :: dijhat(:,:)
     real(kind=DP), allocatable :: dij_hfco(:,:)

     ! maximum angular momentum of any orbital
     integer :: orb_lmax

     ! the augmentation functions and partial wave occupancies
     real(kind=DP), allocatable :: shape_l(:,:)
     real(kind=DP), allocatable :: grad_shape_l(:,:)
     real(kind=DP), allocatable :: qijL(:,:,:)
     real(kind=DP), allocatable :: qfunc(:,:,:)
     real(kind=DP), allocatable :: rhoij0(:,:)
     real(kind=DP), allocatable :: rhoij(:,:)

     ! the single-electron wavefunctions, eigenvalues and occupancies
     integer :: norbs
     integer, allocatable :: orb_ang_mom(:)
     integer, allocatable :: polarise(:)
     real(kind=DP),allocatable :: psi_cin(:,:)
     real(kind=DP),allocatable :: occ(:)
     real(kind=DP),allocatable :: psi_r(:,:)
     real(kind=DP),allocatable :: eigs(:)
     real(kind=DP),allocatable :: splitnorm(:,:)

     ! the total energy and its components
     real(kind=DP) :: total_energy
     real(kind=DP) :: kinetic_energy
     real(kind=DP) :: hartree_energy
     real(kind=DP) :: xc_energy
     real(kind=DP) :: hf_total_energy
     real(kind=DP) :: band_energy
     real(kind=DP) :: xc_dc_energy
     real(kind=DP) :: locpsp_energy
     real(kind=DP) :: nonlocal_energy
     real(kind=DP) :: paw_sphere_energies(paw_en_size)

     ! total energy history
     real(kind=DP) :: old_total_energy
     real(kind=DP) :: older_total_energy

     ! the guess valence electron density (including guess charge and spin)
     real(kind=DP),allocatable :: guess_den(:,:)
     real(kind=DP),allocatable :: guess_comp_den(:,:)
     real(kind=DP),allocatable :: guess_occ(:,:)
     real(kind=DP),allocatable :: guess_rhoij(:,:,:)

     ! lpl: Compensated charged stuff for ddec_mod. Irrelevant if atom_solve
     !      is called without setting 'include_rcomppot' to .true.
     real(kind=DP) :: rcomp      ! Compensation shell radius
     real(kind=DP) :: qcomp      ! Total charge on compensation sphere
     integer :: rcomp_ipt
     logical :: include_rcomppot ! Whether to include comppot in SCF
     logical :: update_rcomppot  ! Whether to update comppot every iter
     real(kind=DP), allocatable :: rcomppot(:) ! Compensation potential
     ! lpl: 'rcomp_pot' gets added to 'locpspot' - original 'locpspot' is
     !      stored in 'locpspot_orig'
     real(kind=DP), allocatable :: locpspot_orig(:)
     ! lpl: Forgot about the 'no compensation charge' case
     logical :: rcomp_extra_scf

     ! pointer to paw_sp in mdl
     type(PAW_SPECIES), pointer :: paw_sp(:)
     type(PSEUDO_SPECIES), pointer :: pseudo_sp(:)

     ! jcap: region label
     integer :: ireg

  end type

  ! Default atom occupations
  ! ndmh: Borrowed from CASTEP - originally from www.webelements.com
  integer, dimension(0:109,28) :: gs_occ
  data gs_occ(0,:)   /10,20,21,30,31,32,40,41,42,43,50,51,52,53,54,60,61,62,63,64,65,70,71,72,73,74,75,76/
  data gs_occ(1,:)   / 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(2,:)   / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(3,:)   / 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(4,:)   / 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(5,:)   / 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(6,:)   / 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(7,:)   / 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(8,:)   / 2, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(9,:)   / 2, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(10,:)  / 2, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(11,:)  / 2, 2, 6, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(12,:)  / 2, 2, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(13,:)  / 2, 2, 6, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(14,:)  / 2, 2, 6, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(15,:)  / 2, 2, 6, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(16,:)  / 2, 2, 6, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(17,:)  / 2, 2, 6, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(18,:)  / 2, 2, 6, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(19,:)  / 2, 2, 6, 2, 6, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(20,:)  / 2, 2, 6, 2, 6, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(21,:)  / 2, 2, 6, 2, 6, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(22,:)  / 2, 2, 6, 2, 6, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(23,:)  / 2, 2, 6, 2, 6, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(24,:)  / 2, 2, 6, 2, 6, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(25,:)  / 2, 2, 6, 2, 6, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(26,:)  / 2, 2, 6, 2, 6, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(27,:)  / 2, 2, 6, 2, 6, 7, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(28,:)  / 2, 2, 6, 2, 6, 8, 2, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(29,:)  / 2, 2, 6, 2, 6,10, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(30,:)  / 2, 2, 6, 2, 6,10, 2, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(31,:)  / 2, 2, 6, 2, 6,10, 2, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(32,:)  / 2, 2, 6, 2, 6,10, 2, 2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(33,:)  / 2, 2, 6, 2, 6,10, 2, 3, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(34,:)  / 2, 2, 6, 2, 6,10, 2, 4, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(35,:)  / 2, 2, 6, 2, 6,10, 2, 5, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(36,:)  / 2, 2, 6, 2, 6,10, 2, 6, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(37,:)  / 2, 2, 6, 2, 6,10, 2, 6, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(38,:)  / 2, 2, 6, 2, 6,10, 2, 6, 0,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(39,:)  / 2, 2, 6, 2, 6,10, 2, 6, 1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(40,:)  / 2, 2, 6, 2, 6,10, 2, 6, 2,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(41,:)  / 2, 2, 6, 2, 6,10, 2, 6, 4,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(42,:)  / 2, 2, 6, 2, 6,10, 2, 6, 5,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(43,:)  / 2, 2, 6, 2, 6,10, 2, 6, 6,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(44,:)  / 2, 2, 6, 2, 6,10, 2, 6, 7,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(45,:)  / 2, 2, 6, 2, 6,10, 2, 6, 8,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(46,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(47,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(48,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(49,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(50,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(51,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(52,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(53,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(54,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(55,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 2, 6, 0,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(56,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,-1, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(57,:)  / 2, 2, 6, 2, 6,10, 2, 6,10, 0, 2, 6, 1,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(58,:)  / 2, 2, 6, 2, 6,10, 2, 6,10, 1, 2, 6, 1,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(59,:)  / 2, 2, 6, 2, 6,10, 2, 6,10, 3, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(60,:)  / 2, 2, 6, 2, 6,10, 2, 6,10, 4, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(61,:)  / 2, 2, 6, 2, 6,10, 2, 6,10, 5, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(62,:)  / 2, 2, 6, 2, 6,10, 2, 6,10, 6, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(63,:)  / 2, 2, 6, 2, 6,10, 2, 6,10, 7, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(64,:)  / 2, 2, 6, 2, 6,10, 2, 6,10, 7, 2, 6, 1,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(65,:)  / 2, 2, 6, 2, 6,10, 2, 6,10, 9, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(66,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,10, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(67,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,11, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(68,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,12, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(69,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,13, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(70,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 0,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(71,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 1,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(72,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 2,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(73,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 3,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(74,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 4,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(75,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 5,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(76,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 6,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(77,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 7,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(78,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6, 9,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(79,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(80,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,-1,-1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(81,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,-1,-1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(82,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,-1,-1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(83,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,-1,-1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(84,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,-1,-1, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(85,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,-1,-1, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(86,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,-1,-1, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  data gs_occ(87,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 0,-1, 2, 6, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0/
  data gs_occ(88,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 0,-1, 2, 6, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(89,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 0,-1, 2, 6, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(90,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 0,-1, 2, 6, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(91,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 2,-1, 2, 6, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(92,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 3,-1, 2, 6, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(93,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 4,-1, 2, 6, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(94,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 6,-1, 2, 6, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(95,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 7,-1, 2, 6, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(96,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 7,-1, 2, 6, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(97,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10, 9,-1, 2, 6, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(98,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,10,-1, 2, 6, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(99,:)  / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,11,-1, 2, 6, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(100,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,12,-1, 2, 6, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(101,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,13,-1, 2, 6, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(102,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,14,-1, 2, 6, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(103,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,14,-1, 2, 6, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(104,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,14,-1, 2, 6, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(105,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,14,-1, 2, 6, 3, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(106,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,14,-1, 2, 6, 4, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(107,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,14,-1, 2, 6, 5, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(108,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,14,-1, 2, 6, 6, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/
  data gs_occ(109,:) / 2, 2, 6, 2, 6,10, 2, 6,10,14, 2, 6,10,14,-1, 2, 6, 7, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0/

  ! Public Subroutines
  public :: atom_create
  public :: atom_destroy
  public :: atom_create_basis
  public :: atom_destroy_basis
  public :: atom_get_lmax
  public :: atom_set_nfunc
  public :: atom_set_nfunc_nl_cdft
  public :: atom_solve
  public :: atom_split_orbital
  public :: atom_polarise_orbital
  public :: atom_write_orbitals
  public :: atom_write_fireball
  public :: atom_set_density_guess

  ! lpl: For accessing gs_occ in ddec_mod
  public :: atom_gs_occ

contains

  subroutine atom_solve(atom,basis,report,num_report_lines,config)

    !==========================================================================!
    ! This subroutine solves for the self-consistent Kohn-Sham orbitals of an  !
    ! isolated atom, expressed in terms of a Spherical Bessel function basis.  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom   (input) : initialised atom, created with atom_create.           !
    !   basis  (input) : initialised basis, created with atom_create_basis.    !
    !   report (input) : array of strings to write SCF report to.              !
    !   num_report_lines (out) : number of lines of report.                    !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use constants, only: periodic_table_name, paw_en_ehart, paw_en_exc, &
         paw_en_exc_dc, paw_en_etxc, paw_en_etxc_dc, &
         paw_en_dij0, paw_en_exc_core, VERBOSE
    use rundat, only: pub_paw, pub_output_detail !, pub_ddec_atomsolve_maxit
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(ATOM_TYPE),intent(inout) :: atom
    type(BASIS_TYPE),intent(in) :: basis
    character(len=128),intent(out) :: report(400)
    character(len=128),intent(in) :: config
    integer, intent(out) :: num_report_lines

    ! Local variables
    character(len=80) :: fmt,tmp,message
    integer :: iscf, il
    integer :: norbs_write
    logical :: converged
    integer :: nscf
    real(kind=DP) :: etol
    real(kind=DP) :: etol_est
    real(kind=DP) :: mix_alpha

    ! lpl: More SCF iterations possibly needed for anionic calculations
    !      in ddec_mod, which often includes compensation charge
    if( atom%rcomp_extra_scf ) then
       nscf = 120 ! pub_ddec_atomsolve_maxit
    else
       nscf = 120
    end if

    mix_alpha = 0.6_DP
    num_report_lines = 0
    etol = 1e-7_DP
    etol_est = 1e-5_DP

    ! ndmh: override to tolerance on Harris-Foulkes estimator
    ! ndmh: in cases where confining potentials for different
    ! ndmh: angular momentum channels are different
    if (any(atom%confpot_width(:)/=atom%confpot_width(0))) then
       etol_est = 1.0_DP
    end if

    ! Write report header
    write(report(num_report_lines+1),*)
    num_report_lines = num_report_lines + 1

    write(report(num_report_lines+1),'(2a,2(a,i3))') &
         'Atom SCF Calculation for ', &
         periodic_table_name(int(atom%atomic_number)), &
         ' : Z (AE atom) = ',int(atom%atomic_number), &
         ' : Z (PS atom) = ',int(atom%ion_charge)
    num_report_lines = num_report_lines + 1

    write(report(num_report_lines+1),'(2a)') 'Config String: ', &
         trim(adjustl(config))
    num_report_lines = num_report_lines + 1

    norbs_write = min(atom%norbs,15)
    write(tmp,'(i10)') norbs_write
    write(fmt,'(5a)') '(a,i3,4x,',trim(adjustl(tmp)),'f5.2)'
    write(report(num_report_lines+1),fmt) 'Orbitals (num,occ):', &
         atom%norbs,atom%occ(1:norbs_write)
    num_report_lines = num_report_lines + 1
    write(fmt,'(5a)') '(a,i3,4x,',trim(adjustl(tmp)),'i5)'
    write(report(num_report_lines+1),fmt) 'Orbitals   (num,l):', &
         atom%norbs,atom%orb_ang_mom(1:norbs_write)
    num_report_lines = num_report_lines + 1

    if (pub_output_detail >= VERBOSE) then

       write(report(num_report_lines+1),'(a,i5,2f16.10)') &
            'Real-space Grid (nr,dr,rmax):',basis%npts,basis%dr,basis%rmax
       num_report_lines = num_report_lines + 1

       write(report(num_report_lines+1),'(a,2i5,f16.10)') &
            'Spherical Bessel Basis (nsws,lmax,cutoff):',basis%nsws_max, &
            basis%lmax,basis%cutoff_energy
       num_report_lines = num_report_lines + 1

       do il=0,atom%orb_lmax
          write(report(num_report_lines+1),'(a,i3,2f16.10)') &
               'Confining potential (l,width,scale):',il, &
               atom%confpot_width(il), atom%confpot_scale
          num_report_lines = num_report_lines + 1
       end do

    end if

    atom%older_total_energy = huge(1.0_DP)
    atom%old_total_energy = huge(1.0_DP)
    atom%paw_sphere_energies(:) = 0.0_DP
    converged = .false.

    ! Initialise Hamiltonian in sw basis for this atom
    call atom_initialise_ham(atom,basis)

    do iscf=1,nscf

       ! Form and then diagonalise the Hamiltonian matrix in the sw basis
       call atom_diagonalise(atom,basis)

       ! Shift old total energies into history
       atom%older_total_energy = atom%old_total_energy
       atom%old_total_energy = atom%total_energy

       ! Harris-Foulkes total energy
       atom%hf_total_energy = atom%band_energy - atom%hartree_energy + &
            atom%xc_dc_energy

       ! Normal total energy
       atom%total_energy = atom%kinetic_energy + atom%locpsp_energy + &
            atom%hartree_energy + atom%xc_energy

       if (pub_paw) then
          ! Add up PAW sphere terms
          atom%hf_total_energy = atom%hf_total_energy &
              - atom%paw_sphere_energies(paw_en_ehart) &
              + atom%paw_sphere_energies(paw_en_exc) &
              - atom%paw_sphere_energies(paw_en_etxc) &
              - atom%paw_sphere_energies(paw_en_exc_dc) &
              + atom%paw_sphere_energies(paw_en_etxc_dc) &
              - atom%paw_sphere_energies(paw_en_exc_core)
          atom%total_energy = atom%total_energy &
               + atom%paw_sphere_energies(paw_en_dij0) &
               + atom%paw_sphere_energies(paw_en_ehart) &
               + atom%paw_sphere_energies(paw_en_exc) &
               - atom%paw_sphere_energies(paw_en_etxc) &
               - atom%paw_sphere_energies(paw_en_exc_core)
       else
          atom%total_energy = atom%total_energy + atom%nonlocal_energy
       end if

       ! Add this cycle to the report
       ! jd: Increased width of total_energy to prevent "output conversion error"
       !     that terminates debug mode runs using extremely crude KE cutoffs
       if (pub_output_detail >= VERBOSE) then
          write(report(num_report_lines+1),'(i3,f22.10,3f16.10)') iscf, &
               atom%total_energy,atom%band_energy, &
               atom%hartree_energy,atom%xc_energy
          num_report_lines = num_report_lines + 1
       end if

       ! Update the Hamiltonian for the new density
       call atom_update_ham(atom,basis,mix_alpha)

       ! Check if we are oscillating in energy, and if so, reduce mixing param
       if ((atom%total_energy>atom%old_total_energy) .and. &
            (atom%old_total_energy<atom%older_total_energy)) then
          if (mix_alpha>0.1_DP) mix_alpha = mix_alpha * 0.8_DP
       end if

       ! Check if we have converged:
       ! Change in energy must be less than etol, and difference of
       ! Harris-Foulkes total energy with normal total energy must
       ! be less than etol_est
       if ((abs(atom%total_energy-atom%old_total_energy)<etol).and. &
            (abs(atom%total_energy-atom%hf_total_energy)<etol_est)) then
          converged = .true.
       end if
       if (converged) exit

    end do

    ! Finalise report
    if (.not.converged) then
       write(message,'(a,i3,a)') 'ERROR in atom_solve: Failed to converge atom &
            &after ',nscf,' SCF iterations'
       write(stdout,'(a)') message
       ! Print report now if convergence failed
       do iscf=1,num_report_lines
          write(stdout,'(a)') trim(report(iscf))
       end do
       call utils_abort(message)
    else
       write(report(num_report_lines+1),'(a,i3,a,f15.8)') &
            'Atom SCF converged after ',iscf,' iterations to a total &
            &energy of ',atom%total_energy
       num_report_lines = num_report_lines + 1
    end if

  end subroutine atom_solve

  subroutine atom_create_basis(basis,rmax_l,lmax,cutoff_energy,fix_deriv)

    !==========================================================================!
    ! This subroutine allocates storage for, and creates the contents of a     !
    ! BASIS_TYPE object describing a spherical wave basis for solving the      !
    ! electronic structure of an atom, plus an associated real-space grid.     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   basis (inout) : BASIS_TYPE object to be initialised.                   !
    !   rmax (input)  : Sphere cutoff radius for the real-space grid and the   !
    !                   Spherical Bessel function basis.                       !
    !   lmax (input)  : Maximum angular momentum value of any wavefunction.    !
    !   cutoff_energy (input) : energy cutoff at which to truncate Spherical   !
    !                           Bessel basis.                                  !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010, in places using ideas from   !
    ! routines by Chris Pickard in ion_atom.F90 of the CASTEP code.            !
    !==========================================================================!

    use services, only: services_regular_integral
    use spherical_wave, only: sw_bessel_zeros_init, sw_bessel_zeros, &
         sw_bessel_accurate
    use utils, only: utils_alloc_check, utils_abort

    implicit none

    ! Arguments
    type(BASIS_TYPE),intent(inout) :: basis
    real(kind=DP),intent(in) :: rmax_l(0:5)
    integer,intent(in) :: lmax
    real(kind=DP),intent(in) :: cutoff_energy
    logical, intent(in) :: fix_deriv

    ! Local variables
    integer :: ierr
    integer :: ipt
    integer :: il
    integer :: isw
    integer :: ibfn
    real(kind=DP) :: qq
    real(kind=DP) :: jp
    real(kind=DP) :: jj,jjp(1:2)

    ! set up basic information about the basis
    allocate(basis%rmax_l(0:lmax),stat=ierr)
    call utils_alloc_check('atom_create_basis','basis%rmax_l',ierr)
    do il=0,lmax
       if (il<=5) then
          basis%rmax_l(il) = rmax_l(il)
       else
          basis%rmax_l(il) = rmax_l(5)
       end if
    end do
    basis%rmax = maxval(rmax_l(:))
    basis%lmax = lmax
    basis%cutoff_energy = cutoff_energy
    basis%dr = 0.125_DP/sqrt(2.0_DP*cutoff_energy)
    basis%npts = int(basis%rmax/basis%dr)
    basis%npts = 2*(basis%npts/2)+1
    basis%dr = basis%rmax / real(basis%npts-1,kind=DP)
    if (fix_deriv) then
       basis%nfixderiv = 1
    else
       basis%nfixderiv = 0
    end if

    ! set up the regular real space grid
    allocate(basis%rad(basis%npts),stat=ierr)
    call utils_alloc_check('atom_create_basis','basis%rad',ierr)
    do ipt=1,basis%npts
       basis%rad(ipt) = basis%dr * real(ipt-1,kind=DP)
    end do

    ! set rmax_l values to exactly lie at previous point on full grid
    do il=0,lmax
       do ipt=1,basis%npts
          if (basis%rad(ipt) > basis%rmax_l(il)) then
             basis%rmax_l(il) = basis%rad(ipt-1)
             exit
          end if
       end do
    end do

    ! set up the spherical bessel functions
    allocate(basis%nsws(0:basis%lmax),stat=ierr)
    call utils_alloc_check('atom_create_basis','basis%nsws',ierr)

    ! guess number of spherical waves required, then keep increasing it until
    ! we have enough
    basis%nsws_max = 36
    basis%nsws(:) = basis%nsws_max + 1
    do
       ! calculate the first nsws_max zeros of the spherical bessel functions
       call sw_bessel_zeros_init(basis%nsws_max,basis%lmax)

       ! count the number of spherical waves with energies up to the cutoff
       do il=0,lmax
          do isw=1,basis%nsws_max
             qq = sw_bessel_zeros(isw,il)/basis%rmax_l(il)
             if (0.5_DP*qq**2 > basis%cutoff_energy) then
                basis%nsws(il) = isw - 1 + basis%nfixderiv
                exit
             end if
          end do
       end do

       if (all(basis%nsws(:)<basis%nsws_max)) exit
       basis%nsws_max = int(basis%nsws_max * 1.4_DP)

    end do

    ! Allocate storage for bessel functions
    basis%nsws_max = maxval(basis%nsws(:))
    allocate(basis%qb(basis%nsws_max,0:lmax),stat=ierr)
    call utils_alloc_check('atom_create_basis','basis%qb',ierr)
    allocate(basis%bessfn(basis%npts,basis%nsws_max,0:lmax),stat=ierr)
    call utils_alloc_check('atom_create_basis','basis%bessfn',ierr)
    allocate(basis%bessnorm(basis%nsws_max,0:lmax),stat=ierr)
    call utils_alloc_check('atom_create_basis','basis%bessnorm',ierr)

    ! Evaluate the bessel functions and their norms
    basis%qb(:,:) = 0.0_DP
    basis%bessfn(:,:,:) = 0.0_DP
    do il=0,lmax
       do isw=1,basis%nsws(il)
          basis%qb(isw,il) = sw_bessel_zeros(isw,il)/basis%rmax_l(il)
          do ipt=1,basis%npts
             if (basis%rad(ipt)<basis%rmax_l(il)) then
                call sw_bessel_accurate(il,basis%qb(isw,il) * basis%rad(ipt), &
                     basis%bessfn(ipt,isw,il),jp)
             else
                basis%bessfn(ipt,isw,il) = 0.0_DP
             end if
          end do
          basis%bessnorm(isw,il) = sqrt(services_regular_integral(basis%npts, &
               basis%dr,(basis%bessfn(:,isw,il)*basis%rad(:))**2))
          basis%bessfn(:,isw,il) = basis%bessfn(:,isw,il)/basis%bessnorm(isw,il)
       end do
    end do

    ! Allocate storage for basis functions
    basis%nbfn_max = basis%nsws_max - basis%nfixderiv
    allocate(basis%nbfn(0:basis%lmax),stat=ierr)
    call utils_alloc_check('atom_create_basis','basis%nbfn',ierr)
    allocate(basis%bcoef(2,1:basis%nsws_max,0:basis%lmax),stat=ierr)
    call utils_alloc_check('atom_create_basis','basis%bcoef',ierr)
    allocate(basis%bfn(basis%npts,basis%nbfn_max,0:lmax),stat=ierr)
    call utils_alloc_check('atom_create_basis','basis%bfn',ierr)

    basis%nbfn(0:lmax) = basis%nsws(0:lmax) - basis%nfixderiv

    if (basis%nfixderiv==0) then
       basis%bfn(1:basis%npts,1:basis%nbfn_max,0:lmax) = &
            basis%bessfn(1:basis%npts,1:basis%nbfn_max,0:lmax)
       basis%bcoef(1,:,:) = 1.0_DP
       basis%bcoef(2,:,:) = 0.0_DP
    else

       ! reduce the basis by fixing the derivatives to zero at the sphere edge
       do il=0,lmax
          do ibfn=1,basis%nbfn(il)
             isw = ibfn
             call sw_bessel_accurate(il,basis%qb(isw,il)*basis%rmax_l(il),jj,jjp(1))
             call sw_bessel_accurate(il,basis%qb(isw+1,il)*basis%rmax_l(il),jj,jjp(2))
             jjp(1) = jjp(1)*basis%qb(isw,il)/basis%bessnorm(isw,il)
             jjp(2) = jjp(2)*basis%qb(isw+1,il)/basis%bessnorm(isw+1,il)
             basis%bcoef(1,ibfn,il) =  jjp(2)/sqrt(jjp(1)**2 + jjp(2)**2)
             basis%bcoef(2,ibfn,il) = -jjp(1)/sqrt(jjp(1)**2 + jjp(2)**2)
             basis%bfn(1:basis%npts,ibfn,il) = &
                  basis%bcoef(1,ibfn,il) * basis%bessfn(1:basis%npts,isw,il) &
                  + basis%bcoef(2,ibfn,il) * basis%bessfn(1:basis%npts,isw+1,il)

          end do
       end do

    end if

  end subroutine atom_create_basis

  subroutine atom_destroy_basis(basis)

    !==========================================================================!
    ! This subroutine deallocates storage associated with a BASIS_TYPE object. !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   basis (input) : BASIS_TYPE object to be deallocated.                   !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(BASIS_TYPE) :: basis

    ! Local variables
    integer :: ierr

    ! deallocate internal arrays in basis type
    deallocate(basis%bfn,stat=ierr)
    call utils_dealloc_check('atom_destroy_basis','basis%bfn',ierr)
    deallocate(basis%bcoef,stat=ierr)
    call utils_dealloc_check('atom_destroy_basis','basis%bcoef',ierr)
    deallocate(basis%nbfn,stat=ierr)
    call utils_dealloc_check('atom_destroy_basis','basis%nbfn',ierr)
    deallocate(basis%bessnorm,stat=ierr)
    call utils_dealloc_check('atom_destroy_basis','basis%bessnorm',ierr)
    deallocate(basis%bessfn,stat=ierr)
    call utils_dealloc_check('atom_destroy_basis','basis%bessfn',ierr)
    deallocate(basis%qb,stat=ierr)
    call utils_dealloc_check('atom_destroy_basis','basis%qb',ierr)
    deallocate(basis%nsws,stat=ierr)
    call utils_dealloc_check('atom_destroy_basis','basis%nsws',ierr)
    deallocate(basis%rad,stat=ierr)
    call utils_dealloc_check('atom_destroy_basis','basis%rad',ierr)
    deallocate(basis%rmax_l,stat=ierr)
    call utils_dealloc_check('atom_destroy_basis','basis%rmax_l',ierr)

  end subroutine atom_destroy_basis

  subroutine atom_get_lmax(lmax,elem,nfunc_target,config,pseudo_sp,paw_sp)

    !==========================================================================!
    ! This subroutine calculates the configuration of an atom and investigates !
    ! the projectors to find the value of lmax that will be required for the   !
    ! basis.                                                                   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   lmax (output)  : The maximum angular momentum value needed in basis    !
    !   elem (input)   : ELEMENT type describing the atom.                     !
    !   nfunc_target (input) : minimum number of NGWF radial functions needed. !
    !   config (input) : string describing the configuration (eg "3s2.0 3p2.0")!
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2011.                                   !
    !==========================================================================!

    use ion, only: element
    use paw, only: paw_get_projector_info
    use pseudopotentials, only: pseudo_get_projector_info
    use rundat, only: pub_paw

    implicit none

    ! Arguments
    integer,intent(out) :: lmax
    type(ELEMENT),intent(in) :: elem
    character(80),intent(in) :: config
    integer,intent(inout) :: nfunc_target
    type(PSEUDO_SPECIES), intent(in), target :: pseudo_sp(:)
    type(PAW_SPECIES), intent(in), target :: paw_sp(:)

    ! Local Variables
    integer :: nshells,nproj_tot,lmax_proj,log_npts_max
    integer :: first_orb,last_orb,norbs,lmax_orb

    ! Get lmax of projectors
    if (.not.pub_paw) then
       call pseudo_get_projector_info(elem%pspecies_number,pseudo_sp, &
            nshells,nproj_tot,lmax_proj)
    else
       call paw_get_projector_info(elem%pspecies_number,paw_sp, &
            nshells,nproj_tot,lmax_proj,log_npts_max)
    end if

    ! Get lmax required for this configuration
    call atom_orb_range(first_orb,last_orb,norbs,lmax_orb, &
         real(elem%atomic_number,kind=DP),elem%ion_charge,config,nfunc_target)

    ! Basis lmax is max of orbital lmax and projector lmax
    lmax = max(lmax_orb,lmax_proj)

  end subroutine atom_get_lmax

  character function atom_l_char(orb_l)

    ! Arguments
    integer, intent(in) :: orb_l

    select case (orb_l)
       case (0); atom_l_char='s'
       case (1); atom_l_char='p'
       case (2); atom_l_char='d'
       case (3); atom_l_char='f'
       case (4); atom_l_char='g'
       case (5); atom_l_char='h'
       case (6); atom_l_char='i'
       case (7); atom_l_char='j'
       case (8); atom_l_char='k'
       case (9); atom_l_char='l'
    end select

  end function atom_l_char

  subroutine atom_get_term_occ(occ,orb,config,atomic_number)

    !==========================================================================!
    ! This subroutine returns the occupancy of a given orbital in the list     !
    ! given in gs_occ according to the atomic number and the configuration     !
    ! override string                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2011.                                   !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: occ
    real(kind=DP), intent(in) :: atomic_number
    integer, intent(in) :: orb
    character(80),intent(in) :: config

    ! Local Variables
    integer :: orb_l, orb_n
    integer :: ipos,jpos
    character :: l_char
    character(4) :: term
    character(80) :: term_str

    ! Create string for this term (eg "2s")
    orb_n = atom_orb_n(orb)
    orb_l = atom_orb_l(orb)
    l_char = atom_l_char(orb_l)
    write(term,'(i3,a1)') orb_n,l_char

    ! Find the term in the config string, if it exists
    ipos = index(config,trim(adjustl(term)))
    if (ipos>0) then
       jpos = index(config(ipos:),' ') + ipos - 1
       if (jpos==0) jpos = len(config)
       term_str = config(ipos:jpos)

       ! Check for ':' (denoting splitnorms)
       jpos = index(term_str,':')
       if (jpos>0) term_str = term_str(:jpos-1)

       ! Check for '|' (denoting polarisation)
       jpos = index(term_str,'|')
       if (jpos>0) term_str = term_str(:jpos-1)

       ! Check for 'X', which means remove this term from the set
       if ((index(term_str,'X')>0).or.(index(term_str,'x')>0)) then
          occ = -100.0_DP
          return
       end if

       ! Read the occupation number from the config string
       read(term_str(3:),*) occ
       if ((occ<-real(orb_l*4+2,DP)).or.(occ>real(orb_l*4+2,DP))) then
          call utils_abort('Error in atom_get_term_occ: Configuration string &
               &could not be parsed')
       end if
    else
       if (orb<=size(gs_occ,2)) then
          ! Use default value from gs_occ table
          occ = gs_occ(int(atomic_number),orb)
          if (occ<0.0_DP) occ=-100.0_DP
       else
          occ = 0.0_DP
       end if
    end if

  end subroutine atom_get_term_occ

  subroutine atom_get_splitnorm(splitnorm,orb,config,pol_in)

    !==========================================================================!
    ! This subroutine returns the occupancy of a given orbital in the list     !
    ! given in gs_occ according to the atomic number and the configuration     !
    ! override string                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2011.                                   !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: splitnorm(4)
    integer, intent(in) :: orb
    character(80),intent(in) :: config
    integer, intent(in), optional :: pol_in

    ! Local Variables
    integer :: orb_l, orb_n
    integer :: ipos,jpos,kpos
    integer :: nsn
    character :: l_char
    character(4) :: term
    character(80) :: term_str
    character(80) :: norm_str,norms_str
    integer :: pol, ipol

    ! Deal with optional flag to read pol splitnorms
    if (present(pol_in)) then
       pol = pol_in
    else
       pol = 0
    end if

    ! Create string for this term (eg "2s")
    orb_n = atom_orb_n(orb)
    orb_l = atom_orb_l(orb)
    l_char = atom_l_char(orb_l)
    write(term,'(i3,a1)') orb_n,l_char

    ! Default is not to split any orbitals
    splitnorm(:) = 0.0_DP

    ! Find this term in the config string, if it exists
    ipos = index(config,trim(adjustl(term)))

    ! No splitnorms present if term is not present
    if (ipos==0) return

    ! Find end of term string
    jpos = index(config(ipos:),' ') + ipos - 1
    if (jpos==0) jpos = len(config)
    term_str = config(ipos:jpos)

    kpos = 1
    do ipol=1,pol

       ! Find next pol division flag '|' in term string
       kpos = index(term_str,'|')
       term_str = term_str(kpos+1:len(term_str))

    end do

    ! Trim off any further pol divisions beyond current
    kpos = index(term_str,'|')
    if (kpos==0) kpos = len(term_str)
    term_str = term_str(1:kpos-1)

    ! Now find the section on splitnorms (anything after the first ':')
    ipos = index(term_str,':')
    if (ipos==0) return
    norms_str = term_str(ipos+1:)

    do nsn=1,4

       ! Check for next instance of ':' or '|'
       ipos = index(norms_str,':')
       jpos = index(norms_str,'|')
       ! If neither found, use whole string
       if ((ipos == 0).and.(jpos == 0)) then
          jpos = len(norms_str)
       else if (ipos>0) then
          jpos = ipos - 1
       else if (jpos>0) then
          jpos = jpos - 1
       end if
       ! Read this section of norms_str
       norm_str = norms_str(1:jpos)
       read(norm_str,*) splitnorm(nsn)

       ! Check for errors
       if (splitnorm(nsn)>1.0_DP) then
          call utils_abort('Error in atom_get_splitnorm: splitnorm string &
               &could not be parsed')
       end if
       if (any(splitnorm(nsn)>splitnorm(1:nsn))) then
          call utils_abort('Error in atom_get_splitnorm: split-valence &
               &norms must be descending')
       end if

       ! Exit if this is the last one, else find the rest of the string
       if (ipos==0) exit
       norms_str = norms_str(ipos+1:)

    end do


  end subroutine atom_get_splitnorm

  subroutine atom_get_pol(test_pol,orb,config)

    !==========================================================================!
    ! This subroutine returns an integer determining whether a given orbital   !
    ! is going to be perturbatively polarised after the atom has been solved.  !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in January 2012.                                !
    !==========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(out) :: test_pol
    integer, intent(in) :: orb
    character(80),intent(in) :: config

    ! Local Variables
    integer :: orb_l, orb_n
    integer :: ipos,jpos
    character :: l_char
    character(4) :: term
    character(80) :: term_str

    ! Create string for this term (eg "2s")
    orb_n = atom_orb_n(orb)
    orb_l = atom_orb_l(orb)
    l_char = atom_l_char(orb_l)
    write(term,'(i3,a1)') orb_n,l_char

    ! Default is not to split any orbitals
    test_pol = 0

    ! Find this term in the config string, if it exists
    ipos = index(config,trim(adjustl(term)))
    if (ipos>0) then
       jpos = index(config(ipos:),' ') + ipos - 1
       if (jpos==0) jpos = len(config)
       term_str = config(ipos:jpos)

       ! Now find sections on polarisation (anything after the '|')
       do
          ipos = index(term_str,'|P')
          ! No more pol sections if ipos==0
          if (ipos==0) return
          test_pol = test_pol + 1
          if (ipos>=len(term_str)-1) return
          term_str = term_str(ipos+1:len(term_str))
       end do

    end if

  end subroutine atom_get_pol

  integer function atom_orb_n(orb)

    implicit none

    ! Arguments
    integer,intent(in) :: orb

    ! Local Variables
    integer :: borb,orb_l,orb_n

    ! Test for orb being bigger than size of gs_occ table
    if (orb>size(gs_occ,2)) then
       borb = size(gs_occ,2)
       do orb_n=8,100
          do orb_l=0,min(orb_n-1,10)
             borb = borb + 1
             if (borb == orb) exit
          end do
          if (borb == orb) exit
       end do
    else
       ! Get values from table
       orb_l = modulo(gs_occ(0,orb),10)
       orb_n = int(gs_occ(0,orb)/10)
    end if
    atom_orb_n = orb_n

  end function atom_orb_n

  integer function atom_orb_l(orb)

    implicit none

    ! Arguments
    integer,intent(in) :: orb

    ! Local Variables
    integer :: borb,orb_l,orb_n

    ! Test for orb being bigger than size of gs_occ table
    if (orb>size(gs_occ,2)) then
       borb = size(gs_occ,2)
       do orb_n=8,100
          do orb_l=0,min(orb_n-1,10)
             borb = borb + 1
             if (borb == orb) exit
          end do
          if (borb == orb) exit
       end do
    else
       ! Get values from table
       orb_l = modulo(gs_occ(0,orb),10)
       orb_n = int(gs_occ(0,orb)/10)
    end if
    atom_orb_l = orb_l

  end function atom_orb_l

  subroutine atom_set_nfunc(nfunc,ngwf_set,elem,cond)

    !==========================================================================!
    ! This subroutine sets up the number of functions required, if this is     !
    ! being determined automatically. In the case of an automatically-set-up   !
    ! multiple-zeta+polarisation basis, it also sets the splitnorms and        !
    ! polarisation directives for each configuration term.                     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   elem (input)     : ELEMENT type describing the atom.                   !
    !   nfunc (inout)    : number of NGWF radial functions needed.             !
    !   ngwf_set (inout) : string describing the configuration                 !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2012.                                   !
    !==========================================================================!

    use ion, only: element
    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer,intent(inout) :: nfunc
    character(*), intent(inout) :: ngwf_set
    type(ELEMENT),intent(in) :: elem
    logical, optional, intent(in) :: cond

    ! Local Variables
    integer :: iorb, jorb
    integer :: first_orb,last_orb,norbs,lmax
    integer :: orb_l, orb_n
    integer :: ipos
    real(kind=DP) :: test_occ
    real(kind=DP) :: atomic_number, ion_charge
    character(len=80) :: solverstr
    character(len=80) :: config
    character(len=4) :: term
    character(len=40) :: termstr
    character(len=20) :: occstr
    character(len=20) :: splitstr
    character(len=20) :: polstr
    character :: l_char
    logical :: loc_cond
    logical :: qzet,tzet,dzet,szet
    logical :: pol,tpol,dpol,spol,npol
    logical :: higher_n

    loc_cond = .false.
    if (present(cond)) loc_cond = cond
    atomic_number = real(elem%atomic_number,kind=DP)
    ion_charge = real(elem%ion_charge,kind=DP)

    ! Find out what kind of basis we want
    qzet=.false.;tzet=.false.;dzet=.false.;szet=.false.;
    tpol=.false.;dpol=.false.;spol=.false.;npol=.false.;
    ngwf_set = trim(adjustl(ngwf_set))
    ipos = index(ngwf_set,'SOLVE')
    if (ipos>0) then
       szet = .true.
       npol = .true.
       ngwf_set = trim(adjustl(ngwf_set(ipos+5:)))
    else ! look for special flags
       ipos = index(ngwf_set,'SZ')
       if (ipos>0) then
          szet = .true.
          ngwf_set = trim(adjustl(ngwf_set(ipos+2:)))
       end if
       ipos = index(ngwf_set,'DZ')
       if (ipos>0) then
          dzet = .true.
          ngwf_set = trim(adjustl(ngwf_set(ipos+2:)))
       end if
       ipos = index(ngwf_set,'TZ')
       if (ipos>0) then
          tzet = .true.
          ngwf_set = trim(adjustl(ngwf_set(ipos+2:)))
       end if
       ipos = index(ngwf_set,'QZ')
       if (ipos>0) then
          qzet = .true.
          ngwf_set = trim(adjustl(ngwf_set(ipos+2:)))
       end if
       npol = .true.
       ipos = index(ngwf_set,'SP')
       if (index(ngwf_set,'P')==1) ipos = 1
       if (ipos>0) then
          spol = .true.; npol = .false.;
          ngwf_set = trim(adjustl(ngwf_set(ipos+2:)))
       end if
       ipos = index(ngwf_set,'DP')
       if (ipos>0) then
          dpol = .true.; npol = .false.;
          ngwf_set = trim(adjustl(ngwf_set(ipos+2:)))
       end if
       ipos = index(ngwf_set,'TP')
       if (ipos>0) then
          tpol = .true.; npol = .false.;
          ngwf_set = trim(adjustl(ngwf_set(ipos+2:)))
       end if
    end if

    ! Abort if none of these have been set
    if ((.not.szet).and.(.not.dzet).and.(.not.tzet).and.(.not.qzet)) then
       call utils_abort('Error in atom_set_nfunc: Could not parse basis &
            &definition string for zeta.')
    end if
    if ((.not.npol).and.(.not.spol).and.(.not.dpol).and.(.not.tpol)) then
       call utils_abort('Error in atom_set_nfunc: Could not parse basis &
            &definition string for pol.')
    end if

    solverstr = trim(adjustl(ngwf_set))

    ! Parse the string following "SOLVE"
    ! now allows partially-set config eg for core holes
    ipos = index(ngwf_set,'conf=')
    if (ipos>0) then
       config = ngwf_set(ipos:)
    else
       config = 'conf='
    end if

    call atom_orb_range(first_orb,last_orb,norbs,lmax, &
         atomic_number,ion_charge,config,nfunc)

    do iorb=first_orb,last_orb

       ! Find occupancy of this term
       call atom_get_term_occ(test_occ,iorb,config,atomic_number)

       ! Skip if <0 occupation (ie term is not meant to be included)
       if (test_occ<-99.0_DP) cycle

       ! Create string for this term (eg "2s")
       orb_n = atom_orb_n(iorb)
       orb_l = atom_orb_l(iorb)
       l_char = atom_l_char(orb_l)
       write(term,'(i3,a1)') orb_n,l_char
       write(occstr,'(f10.1)') test_occ

       ! For conduction NGWFS, if subshell is filled and there is a higher-n
       ! shell of the same angular momentum, then exclude this term from the
       ! final NGWFs by setting the splitnorm to -1
       if (loc_cond) then
          if (abs(test_occ-real(2*(2*orb_l+1),kind=DP))<epsilon(1.0_DP)) then
             higher_n = .false.
             do jorb=iorb+1,last_orb
                if (atom_orb_n(jorb)>orb_n) higher_n = .true.
             end do
             if (higher_n) then
                write(occstr,'(f10.1,a)') test_occ,':-1'
             end if
          end if
       end if

       polstr = ''
       splitstr = ''

       ! If term is occupied, set multiple-zeta splitting string
       if (test_occ>0.0_DP) then
          if (dzet) write(splitstr,'(1(a,f4.2))') ':',0.05_DP
          if (tzet) write(splitstr,'(2(a,f4.2))') ':',0.15_DP,':',0.06_DP
          if (qzet) write(splitstr,'(3(a,f4.2))') ':',0.15_DP,':',0.06_DP, &
               ':',0.02_DP

          ! Set polarisation string: true unless an l+1 term already exists
          pol = .true.
          do jorb=iorb+1,last_orb
             if (atom_orb_l(jorb)==orb_l+1) pol = .false.
          end do
          if (pol) then
             if (spol) write(polstr,'(a)') '|P'
             if (dpol) write(polstr,'(a)') '|P:0.15'
             if (tpol) write(polstr,'(a)') '|P:0.15:0.05'
          end if

       end if

       write(termstr,'(4a)') trim(adjustl(term)), &
            trim(adjustl(occstr)),trim(adjustl(splitstr)), &
            trim(adjustl(polstr))
       config = trim(adjustl(config))//' '//trim(adjustl(termstr))

    end do

    ! Determine orbital range & count again, this time with splitnorms and pol
    nfunc = -1

    call atom_orb_range(first_orb,last_orb,norbs,lmax, &
         atomic_number,ion_charge,config,nfunc)

    write(ngwf_set,'(4a)') 'SOLVE ',trim(adjustl(config)),' ', &
         trim(adjustl(solverstr))

  end subroutine atom_set_nfunc


  subroutine atom_set_nfunc_nl_cdft(nfunc,ngwf_set,elem)

    !==========================================================================!
    ! This subroutine sets up the number of functions required, if this is     !
    ! being determined automatically. In the case of an automatically-set-up   !
    ! multiple-zeta+polarisation basis, it also sets the splitnorms and        !
    ! polarisation directives for each configuration term.                     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   elem (input)     : ELEMENT type describing the atom.                   !
    !   nfunc (inout)    : number of NGWF radial functions needed.             !
    !   ngwf_set (inout) : string describing the configuration                 !
    !--------------------------------------------------------------------------!
    ! cDFT-adaptation (Jul 2013) of atom_set_nfunc by Nicholas Hine            !
    !==========================================================================!

    use ion, only: element
    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer,intent(inout) :: nfunc
    character(*), intent(in) :: ngwf_set
    type(ELEMENT),intent(inout) :: elem

    ! Local Variables
    integer :: iorb
    integer :: first_orb,last_orb,norbs,lmax
    integer :: orb_l, orb_n
    integer :: ipos
    integer :: kk, orb_l_old
    real(kind=DP) :: test_occ
    real(kind=DP) :: atomic_number, ion_charge
    character(len=80) :: config

    atomic_number = real(elem%atomic_number,kind=DP)
    ion_charge = real(elem%ion_charge,kind=DP)

    ! Parse the string following "SOLVE"
    ipos = index(ngwf_set,'conf=')
    if (ipos>0) then
       config = ngwf_set(ipos+5:)
    else
       config = ''
    end if

    call atom_orb_range(first_orb,last_orb,norbs,lmax, &
         atomic_number,ion_charge,config,nfunc)

    !gibo: for cdft_multi_proj, initialise counter for ang.mom of each function
    kk = 0
    orb_l_old = -1 ! dummy of 'old' ang.mom. number
    elem%num_ang_mom_shells = 0
    elem%orb_ang_mom(:) = -1

    do iorb=first_orb,last_orb
       ! Find occupancy of this term
       call atom_get_term_occ(test_occ,iorb,config,atomic_number)
       ! Skip if <0 occupation (ie term is not meant to be included)
       if (test_occ<-99.0_DP) cycle

       ! work out principal (_n) and angular momentum (_l) quantum numbers
       orb_n = atom_orb_n(iorb)
       orb_l = atom_orb_l(iorb)

       !gibo: for cdft_multi_proj, store ang.mom of this function
       !      (needed for hubbard_species_proj)
       if ( (orb_l .ne. orb_l_old) .AND. (test_occ >= 0.0_DP)) then
        elem%num_ang_mom_shells = elem%num_ang_mom_shells + 1
        kk = kk + 1
        elem%orb_ang_mom(kk) = orb_l
        orb_l_old = orb_l
       endif
    end do

  end subroutine atom_set_nfunc_nl_cdft


  subroutine atom_orb_range(first_orb,last_orb,norbs,lmax, &
       atomic_number,ion_charge,config,nfunc_target)

    !==========================================================================!
    ! This subroutine returns the range of terms in the gs_occ list required   !
    ! for the given configuration.                                             !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2011.                                   !
    !==========================================================================!

    use utils, only: utils_abort, utils_real_to_str

    implicit none

    ! Arguments
    integer,intent(out) :: first_orb,last_orb
    integer,intent(out) :: norbs,lmax
    real(kind=DP), intent(in) :: atomic_number,ion_charge
    integer, intent(inout) :: nfunc_target
    character(80),intent(in) :: config

    ! Local Variables
    integer :: iorb,norbs_tot
    integer :: nsplitnorm,nsplitnorm_pol,test_pol,ipol
    real(kind=DP) :: charge_left, test_occ
    real(kind=DP) :: test_splitnorm(4), test_splitnorm_pol(4)
    real(kind=DP) :: orb_charge
    logical :: set_nfunc_target
    logical :: sp_complete

    set_nfunc_target = .false.
    if (nfunc_target < 0) set_nfunc_target = .true.
    nsplitnorm_pol = 0

    ! Find the first term needed to get the right total charge
    charge_left = atomic_number
    norbs_tot = 0
    first_orb = 999
    do iorb=1,size(gs_occ,2)
       if (charge_left == ion_charge) then
          first_orb = iorb
          exit
       end if
       call atom_get_term_occ(orb_charge,iorb,config,atomic_number)
       if (orb_charge<0.0_DP) orb_charge = 0.0_DP
       charge_left = charge_left - orb_charge
    end do
    if (first_orb==999) then
       call utils_abort('Error in atom_orb_range: Failed to determine &
            &first orbital in NGWF setup. Atomic number: '//&
            trim(utils_real_to_str(atomic_number))//', ion charge: '//&
            trim(utils_real_to_str(ion_charge))//', config: "'//&
            trim(adjustl(config))//'".')
    end if

    ! Find the last term needed to reach the required number of orbitals
    charge_left = ion_charge
    do last_orb=first_orb,500

       ! Find occupancy of this term
       call atom_get_term_occ(test_occ,last_orb,config,atomic_number)

       ! Find splitnorms of this term and count the number of shells
       ! which will eventually be generated for this term
       call atom_get_splitnorm(test_splitnorm,last_orb,config)
       ! Count splitnorms, excluding any terms where splitnorm == -1
       nsplitnorm = count(test_splitnorm(:)>0.0_DP) - &
            count(test_splitnorm(:)==-1.0_DP) + 1
       ! Skip if -100 occupation (ie term is not meant to be included at all)
       if (test_occ<-99.0_DP) cycle

       ! Negative occupancy means count the charge, but do not include the
       ! function in the set
       if (test_occ<0.0_DP) then
          charge_left = charge_left + test_occ
       else
          charge_left = charge_left - test_occ
       end if

       ! Skip if <0 occupation (ie term is not meant to be included)
       if (test_occ<0.0_DP) cycle

       ! Add number of orbitals (including m terms) in this term to total
       norbs_tot = norbs_tot + (2*atom_orb_l(last_orb) + 1) * &
            nsplitnorm

       ! Find whether this orbital will be perturbatively polarised
       ! afterwards (flexibility for how many times, later)
       call atom_get_pol(test_pol,last_orb,config)
       do ipol=1,test_pol
          ! Find splitnorms of this polarisation orbital for this term
          call atom_get_splitnorm(test_splitnorm_pol,last_orb,config,ipol)
          ! Count pol splitnorms, excluding any terms where splitnorm == -1
          nsplitnorm_pol = count(test_splitnorm_pol(:)>0.0_DP) - &
               count(test_splitnorm_pol(:)==-1.0_DP) + 1
          ! Add number of orbitals which will later result from pertubative
          ! polarisation of this term to total
          norbs_tot = norbs_tot + (2*(atom_orb_l(last_orb) + ipol) + 1) * &
               test_pol * nsplitnorm_pol
       end do

       ! Check if we have enough orbitals already (if fixed target)
       if ((.not.set_nfunc_target).and.(norbs_tot >= nfunc_target)) then
          exit
       end if

       ! Check if this is an s-function for n>1, in which case we need
       ! to also add the p-function
       sp_complete=.true.
       if ((atom_orb_l(last_orb)==0).and.(atom_orb_n(last_orb)>1)) then
          sp_complete = .false.
       end if
       ! Hacky special case to avoid Palladium breaking with auto set-up
       ! (would be nice to fix in a more robust way)
       if (set_nfunc_target.and.(atomic_number==46)) then
          if ((atom_orb_l(last_orb)==2).and.(charge_left==0.0_DP)) then
             sp_complete = .false.
          end if
       end if
       ! Check if we have enough charge already (if determining target)
       if (set_nfunc_target.and.(charge_left==0.0_DP).and.sp_complete) then
          nfunc_target = norbs_tot
          exit
       end if

    end do
    if (norbs_tot<nfunc_target) then
       call utils_abort('Error in atom_orb_range: Failed to determine &
            &last orbital in NGWF setup. Check config.')
    end if
    if ((last_orb>100).and.(set_nfunc_target)) then
       call utils_abort('Error in atom_orb_range: Failed to determine &
            &last orbital in NGWF auto-setup. Check config.')
    end if

    ! Now count again, testing for any terms to throw away
    ! and also getting max angular momentum of any orbital as we go
    norbs = last_orb - first_orb + 1
    norbs_tot = norbs
    last_orb = first_orb
    lmax = -1
    do iorb=1,norbs_tot
       call atom_get_term_occ(test_occ,last_orb,config,atomic_number)
       call atom_get_pol(test_pol,last_orb,config)
       if (test_occ<0.0_DP) then
          norbs = norbs - 1
       else
          ! keep a running count on the maximum l value encountered
          lmax = max(lmax,atom_orb_l(last_orb) + test_pol)
          norbs = norbs + test_pol
       end if
       last_orb = last_orb + 1
    end do
    last_orb = last_orb - 1

  end subroutine atom_orb_range

  subroutine atom_create(atom,basis,elem,nfunc_target,config,cwidth,cscale, &
       pseudo_sp,paw_sp)

    !==========================================================================!
    ! This subroutine allocates storage for an ATOM_TYPE object and sets up    !
    ! sizes based on the element type passed in.                               !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (inout)   : ATOM_TYPE object to be initialised                    !
    !   basis (input)  : BASIS_TYPE object, already initialised.               !
    !   elem (input)   : ELEMENT type describing the atom.                     !
    !   config (input) : string describing the configuration (eg "3s2.0 3p2.0")!
    !   nfunc_target (input) : minimum number of NGWF radial functions needed. !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use ion, only: element
    use paw, only: paw_get_projector_info
    use pseudopotentials, only: pseudo_get_projector_info
    use rundat, only: pub_paw
    use utils, only: utils_alloc_check, utils_abort

    implicit none

    ! Arguments
    type(ATOM_TYPE),intent(inout) :: atom
    type(BASIS_TYPE),intent(in) :: basis
    type(ELEMENT),intent(in) :: elem
    character(80),intent(in) :: config
    integer,intent(inout) :: nfunc_target
    real(kind=DP),intent(in) :: cwidth(0:5), cscale
    type(PAW_SPECIES), intent(in), target :: paw_sp(:)
    type(PSEUDO_SPECIES), intent(in), target :: pseudo_sp(:)

    ! Local variables
    integer :: ierr
    integer :: npts
    integer :: iorb, jorb, korb, first_orb, last_orb, il, ipol
    real(kind=DP) :: test_occ

    ! Initialise variables based on inputs
    npts = basis%npts
    atom%ion_charge = elem%ion_charge
    atom%atomic_number = elem%atomic_number
    atom%val_charge = atom%ion_charge
    atom%pspecies_number = elem%pspecies_number
    allocate(atom%confpot_width(0:basis%lmax),stat=ierr)
    call utils_alloc_check('atom_create','atom%confpot_width',ierr)
    do il=0,basis%lmax
       if (il<=5) then
          atom%confpot_width(il) = cwidth(il)
       else
          atom%confpot_width(il) = cwidth(5)
       end if
    end do
    atom%confpot_scale = cscale

    ! Set up PAW or NCPSP information
    atom%paw_sp => paw_sp
    atom%pseudo_sp => pseudo_sp
    if (.not.pub_paw) then
       call pseudo_get_projector_info(atom%pspecies_number,atom%pseudo_sp, &
            atom%nshells,atom%npwtot,atom%proj_lmax)
    else
       call paw_get_projector_info(atom%pspecies_number,atom%paw_sp, &
            atom%nshells,atom%npwtot,atom%proj_lmax,atom%log_npts_max)
    end if
    if (atom%proj_lmax > basis%lmax) then
       call utils_abort('Error in atom_create: Angular momentum of &
              &a projector is greater than basis%lmax')
    end if

    ! lpl: Compensated charge stuff used in ddec_mod. These values
    !      are initialized to be inert unless altered.
    atom%rcomp = 0.0_DP
    atom%qcomp = 0.0_DP
    atom%include_rcomppot = .false.
    atom%update_rcomppot  = .false.
    atom%rcomp_extra_scf  = .false.
    allocate(atom%rcomppot(npts),stat=ierr)
    call utils_alloc_check('atom_create','atom%rcomppot',ierr)
    allocate(atom%locpspot_orig(npts),stat=ierr)
    call utils_alloc_check('atom_create','atom%locpspot_orig',ierr)

    allocate(atom%den(npts),stat=ierr)
    call utils_alloc_check('atom_create','atom%den',ierr)
    allocate(atom%comp_den(npts),stat=ierr)
    call utils_alloc_check('atom_create','atom%comp_den',ierr)
    allocate(atom%core_den(npts),stat=ierr)
    call utils_alloc_check('atom_create','atom%core_den',ierr)
    allocate(atom%locpspot(npts),stat=ierr)
    call utils_alloc_check('atom_create','atom%locpspot',ierr)
    allocate(atom%confpot(npts,0:5),stat=ierr)
    call utils_alloc_check('atom_create','atom%confpot',ierr)
    allocate(atom%effpot(npts),stat=ierr)
    call utils_alloc_check('atom_create','atom%effpot',ierr)
    allocate(atom%xcpot(npts),stat=ierr)
    call utils_alloc_check('atom_create','atom%xcpot',ierr)
    allocate(atom%work(npts),stat=ierr)
    call utils_alloc_check('atom_create','atom%work',ierr)
    allocate(atom%proj_ang_mom(atom%nshells),stat=ierr)
    call utils_alloc_check('atom_create','atom%proj_ang_mom',ierr)
    allocate(atom%proj_q(basis%nsws_max,atom%nshells),stat=ierr)
    call utils_alloc_check('atom_create','atom%proj_q',ierr)
    allocate(atom%dij(atom%nshells,atom%nshells),stat=ierr)
    call utils_alloc_check('atom_create','atom%dij',ierr)
    allocate(atom%dij0(atom%nshells,atom%nshells),stat=ierr)
    call utils_alloc_check('atom_create','atom%dij0',ierr)
    allocate(atom%dij_hfco(atom%nshells,atom%nshells),stat=ierr)
    call utils_alloc_check('atom_create','atom%dij_hfco',ierr)
    if (pub_paw) then
       allocate(atom%dij_xc(atom%nshells,atom%nshells),stat=ierr)
       call utils_alloc_check('atom_create','atom%dij_xc',ierr)
       allocate(atom%dij_hart(atom%nshells,atom%nshells),stat=ierr)
       call utils_alloc_check('atom_create','atom%dij_hart',ierr)
       allocate(atom%dijhat(atom%nshells,atom%nshells),stat=ierr)
       call utils_alloc_check('atom_create','atom%dijhat',ierr)
    end if

    call atom_orb_range(first_orb,last_orb,atom%norbs,atom%orb_lmax, &
         atom%atomic_number,atom%ion_charge,config,nfunc_target)

    ! Allocate arrays for wavefunction
    allocate(atom%orb_ang_mom(atom%norbs),stat=ierr)
    call utils_alloc_check('atom_create','atom%orb_ang_mom',ierr)
    allocate(atom%psi_cin(basis%nbfn_max,atom%norbs),stat=ierr)
    call utils_alloc_check('atom_create','atom%psi_cin',ierr)
    allocate(atom%occ(atom%norbs),stat=ierr)
    call utils_alloc_check('atom_create','atom%occ',ierr)
    allocate(atom%psi_r(npts,atom%norbs),stat=ierr)
    call utils_alloc_check('atom_create','atom%psi_r',ierr)
    allocate(atom%eigs(1:atom%norbs),stat=ierr)
    call utils_alloc_check('atom_create','atom%eigs',ierr)
    allocate(atom%splitnorm(4,1:atom%norbs),stat=ierr)
    call utils_alloc_check('atom_create','atom%splitnorm',ierr)
    allocate(atom%polarise(atom%norbs),stat=ierr)
    call utils_alloc_check('atom_create','atom%polarise',ierr)

    ! Get atom occupation numbers and angular momenta
    iorb = 1
    do jorb=first_orb,last_orb
       call atom_get_term_occ(test_occ,jorb,config,atom%atomic_number)
       if (test_occ < 0.0_DP) cycle
       atom%occ(iorb) = test_occ
       atom%orb_ang_mom(iorb) = atom_orb_l(jorb)
       if (atom%orb_ang_mom(iorb) > basis%lmax) then
          write(stdout,*) iorb,jorb,atom_orb_l(jorb),basis%lmax
          call utils_abort('Error in atom_create: Angular momentum of &
               &a state is greater than basis%lmax')
       end if
       call atom_get_splitnorm(atom%splitnorm(:,iorb),jorb,config)
       ! Check if we are going to do perturbative polarisation for this func
       call atom_get_pol(atom%polarise(iorb),jorb,config)
       iorb = iorb + 1
    end do

    ! Check again, adding any perturbative polarisation orbitals on to the end
    ! of the list
    korb = iorb
    iorb = 1
    do jorb=first_orb,last_orb
       call atom_get_term_occ(test_occ,jorb,config,atom%atomic_number)
       if (test_occ < 0.0_DP) cycle
       do ipol=1,atom%polarise(iorb)
          atom%occ(korb) = 0.0_DP
          atom%orb_ang_mom(korb) = atom_orb_l(jorb) + ipol
          atom%polarise(korb) = 0
          call atom_get_splitnorm(atom%splitnorm(:,korb),jorb,config,ipol)
          if (atom%orb_ang_mom(iorb) > basis%lmax) then
             write(stdout,*) iorb,jorb,atom_orb_l(jorb),basis%lmax
             call utils_abort('Error in atom_create: Angular momentum of &
                  &a pol state is greater than basis%lmax')
          end if
          korb = korb + 1
       end do
       iorb = iorb + 1
    end do

    ! Set total valence charge
    atom%val_charge = sum(atom%occ)

    atom%orb_lmax = maxval(atom%orb_ang_mom(:))

    atom%total_energy = 0.0_DP
    atom%hartree_energy = 0.0_DP
    atom%kinetic_energy = 0.0_DP
    atom%locpsp_energy = 0.0_DP
    atom%xc_energy = 0.0_DP
    atom%xc_dc_energy = 0.0_DP
    atom%old_total_energy = 0.0_DP

    ! Allocate PAW compensation density information
    allocate(atom%qijL(max(atom%nshells,1),max(atom%nshells,1),0:atom%proj_lmax),stat=ierr)
    call utils_alloc_check('atom_create','atom%qijL',ierr)
    allocate(atom%shape_l(basis%npts,0:atom%proj_lmax),stat=ierr)
    call utils_alloc_check('atom_create','atom%shape_l',ierr)
    allocate(atom%grad_shape_l(basis%npts,0:atom%proj_lmax),stat=ierr)
    call utils_alloc_check('atom_create','atom%grad_shape_l',ierr)
    allocate(atom%rhoij(atom%nshells,atom%nshells),stat=ierr)
    call utils_alloc_check('atom_create','atom%rhoij',ierr)
    allocate(atom%rhoij0(atom%nshells,atom%nshells),stat=ierr)
    call utils_alloc_check('atom_create','atom%rhoij0',ierr)
    allocate(atom%qfunc(basis%npts,atom%nshells,atom%nshells),stat=ierr)
    call utils_alloc_check('atom_create','atom%qfunc',ierr)

  end subroutine atom_create

  subroutine atom_destroy(atom)

    !==========================================================================!
    ! This subroutine deallocates memory for the ATOM_TYPE object.             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (input) : ATOM_TYPE object to be deallocated.                     !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use rundat, only: pub_paw
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(ATOM_TYPE),intent(inout) :: atom

    ! Local variables
    integer :: ierr

    ! deallocate all storage associated with the atom
    if (allocated(atom%guess_rhoij)) then
       deallocate(atom%guess_rhoij,stat=ierr)
       call utils_dealloc_check('atom_destroy','atom%guess_rhoij',ierr)
    end if
    if (allocated(atom%guess_occ)) then
       deallocate(atom%guess_occ,stat=ierr)
       call utils_dealloc_check('atom_destroy','atom%guess_occ',ierr)
    end if
    if (allocated(atom%guess_comp_den)) then
       deallocate(atom%guess_comp_den,stat=ierr)
       call utils_dealloc_check('atom_destroy','atom%guess_comp_den',ierr)
    end if
    if (allocated(atom%guess_den)) then
       deallocate(atom%guess_den,stat=ierr)
       call utils_dealloc_check('atom_destroy','atom%guess_den',ierr)
    end if
    deallocate(atom%qfunc,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%qfunc',ierr)
    deallocate(atom%rhoij0,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%rhoij0',ierr)
    deallocate(atom%rhoij,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%rhoij',ierr)
    deallocate(atom%grad_shape_l,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%grad_shape_l',ierr)
    deallocate(atom%shape_l,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%shape_l',ierr)
    deallocate(atom%qijL,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%qijL',ierr)
    deallocate(atom%splitnorm,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%splitnorm',ierr)
    deallocate(atom%eigs,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%eigs',ierr)
    deallocate(atom%psi_r,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%psi_r',ierr)
    deallocate(atom%occ,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%occ',ierr)
    deallocate(atom%psi_cin,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%psi_cin',ierr)
    deallocate(atom%polarise,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%polarise',ierr)
    deallocate(atom%orb_ang_mom,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%orb_ang_mom',ierr)
    if (pub_paw) then
       deallocate(atom%dijhat,stat=ierr)
       call utils_dealloc_check('atom_destroy','atom%dijhat',ierr)
       deallocate(atom%dij_hart,stat=ierr)
       call utils_dealloc_check('atom_destroy','atom%dij_hart',ierr)
       deallocate(atom%dij_xc,stat=ierr)
       call utils_dealloc_check('atom_destroy','atom%dij_xc',ierr)
    end if
    deallocate(atom%dij0,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%dij0',ierr)
    deallocate(atom%dij_hfco,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%dij_hfco',ierr)
    deallocate(atom%dij,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%dij',ierr)
    deallocate(atom%proj_q,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%proj_q',ierr)
    deallocate(atom%proj_ang_mom,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%proj_ang_mom',ierr)
    deallocate(atom%work,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%work',ierr)
    deallocate(atom%xcpot,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%xcpot',ierr)
    deallocate(atom%effpot,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%effpot',ierr)
    deallocate(atom%confpot,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%confpot',ierr)
    deallocate(atom%locpspot,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%locpspot',ierr)
    deallocate(atom%core_den,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%core_den',ierr)
    deallocate(atom%comp_den,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%comp_den',ierr)
    deallocate(atom%den,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%den',ierr)
    deallocate(atom%confpot_width,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%confpot_width',ierr)

    ! lpl: Compensated charge stuff used in ddec_mod
    deallocate(atom%rcomppot,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%rcomppot',ierr)
    deallocate(atom%locpspot_orig,stat=ierr)
    call utils_dealloc_check('atom_destroy','atom%locpspot_orig',ierr)

  end subroutine atom_destroy

  subroutine atom_initialise_ham(atom,basis)

    !==========================================================================!
    ! This subroutine initialises the arrays required to calculate the         !
    ! Kohn-Sham Hamiltonian for the atom. Interfaces with pseudopotential      !
    ! modules to get required information on local potentials, core densities, !
    ! projectors, etc.                                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (input) : ATOM_TYPE object describing the atom.                   !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use constants, only: PI, paw_en_exc_core
    use pseudopotentials, only: pseudo_get_locpot_rad, &
         pseudo_get_core_den_rad, pseudo_get_projectors_q, pseudo_get_aug_funcs
    use paw, only: paw_get_locpot_rad, paw_get_core_den_rad, &
         paw_get_projectors_q, paw_get_aug_funcs, paw_exc_core_atom
    use rundat, only: pub_paw
    use services, only: services_regular_integral

    implicit none

    ! Arguments
    type(ATOM_TYPE),intent(inout) :: atom
    type(BASIS_TYPE),intent(in) :: basis

    ! Local Variables
    integer :: ir
    integer :: ishell, isw, ibfn, il
    real(kind=DP) :: den_int
    real(kind=DP) :: r_onset, delta
    real(kind=DP) :: sph_harm_fac

    ! Get the local potential
    if (.not.pub_paw) then
       call pseudo_get_locpot_rad(atom%locpspot,basis%npts,basis%rad, &
            atom%pspecies_number,atom%pseudo_sp)
    else
       call paw_get_locpot_rad(atom%locpspot,basis%npts,basis%rad, &
            atom%pspecies_number,atom%paw_sp)
    end if

    ! lpl: Compensated charge stuff used in ddec_mod
    atom%rcomp_ipt = int(aint( atom%rcomp/basis%dr )) ! jd: int() silences warning
    atom%locpspot_orig(:) = atom%locpspot(:)

    ! Get the core density
    if (.not.pub_paw) then
       call pseudo_get_core_den_rad(atom%core_den,basis%npts,basis%rad, &
            atom%pspecies_number,atom%pseudo_sp)
    else
       call paw_get_core_den_rad(atom%core_den,basis%npts,basis%rad, &
            atom%pspecies_number,atom%paw_sp)
    end if

    ! Get the projectors (and shape functions and rhoij0, in PAW)
    if (.not.pub_paw) then
       call pseudo_get_projectors_q(atom%proj_q,atom%dij0, &
            atom%proj_ang_mom,atom%nshells,basis%nsws,maxval(basis%nsws), &
            basis%lmax,basis%qb,atom%pspecies_number,atom%pseudo_sp)
       call pseudo_get_aug_funcs(atom%qijL(:,:,0),atom%qfunc, &
            atom%nshells,basis%npts,basis%rad,atom%pspecies_number,atom%pseudo_sp)
       atom%rhoij0(:,:) = 0.0_DP
    else
       call paw_get_projectors_q(atom%proj_q,atom%dij0,atom%proj_ang_mom, &
            atom%nshells,basis%nsws,maxval(basis%nsws),basis%lmax,basis%qb, &
            atom%pspecies_number,atom%paw_sp)
       call paw_get_aug_funcs(atom%qijL,atom%shape_l,atom%grad_shape_l, &
            atom%rhoij0,atom%proj_lmax,atom%nshells,basis%npts,basis%rad, &
            atom%pspecies_number,atom%paw_sp)
       do il=0,atom%proj_lmax
          sph_harm_fac = sqrt(4.0_DP*PI)*sqrt(real(2*il+1,kind=DP)/(4.0_DP*PI))
          atom%shape_l(:,il) = atom%shape_l(:,il)*sph_harm_fac
          den_int = services_regular_integral(basis%npts,basis%dr, &
               atom%shape_l(:,il)*basis%rad(:)**(2+il))
          atom%shape_l(:,il) = atom%shape_l(:,il) / den_int
       end do
    end if

    ! Normalise the projectors by the integrals of the basis functions
    do ishell=1,atom%nshells
       il = atom%proj_ang_mom(ishell)
       do isw=1,basis%nsws(il)
          atom%proj_q(isw,ishell) = atom%proj_q(isw,ishell) / &
               basis%bessnorm(isw,il)
       end do
    end do

    ! Take linear combinations of proj_q's to match basis functions
    if (basis%nfixderiv>0) then
       do ishell=1,atom%nshells
          il = atom%proj_ang_mom(ishell)
          do ibfn=1,basis%nbfn(il)
             atom%proj_q(ibfn,ishell) = &
                  basis%bcoef(1,ibfn,il) * atom%proj_q(ibfn,ishell) &
                  + basis%bcoef(2,ibfn,il) * atom%proj_q(ibfn+1,ishell)
          end do
       end do
    end if

    ! Calculate the core xc energy in PAW
    if (pub_paw) then
       call paw_exc_core_atom(atom%paw_sphere_energies(paw_en_exc_core), &
            atom%pspecies_number,atom%paw_sp)
    end if

    ! Set up a confining potential (see V. Blum, R. Gehrke, F. Hanke,
    ! P. Havu, V. Havu, X. Ren, K. Reuter, and M. Scheffler, Comp. Phys.
    ! Comm. 180, 2175-2196 (2009) for more details).
    do il=0,min(atom%orb_lmax,5)
       r_onset = basis%rmax_l(il) - atom%confpot_width(il)
       delta = 0.00001_DP
       do ir=1,basis%npts
          if (basis%rad(ir)>r_onset) then
             atom%confpot(ir,il) = atom%confpot_scale * &
                  exp(-atom%confpot_width(il)/(basis%rad(ir)-r_onset)) * &
                  1.0_DP / (basis%rad(ir)-basis%rmax+delta)**2
          else
             atom%confpot(ir,il) = 0.0_DP
          end if
       end do
    end do

    ! Guess an initial density
    do ir=1,basis%npts
       atom%den(ir) = exp(-0.5_DP*basis%rad(ir)**2)
    end do

    ! Normalise the initial density
    atom%work(:) = atom%den(:)*basis%rad(:)**2
    den_int = services_regular_integral(basis%npts,basis%dr,atom%work)
    atom%den(:) = atom%val_charge * atom%den(:) / den_int

    ! Guess initial rhoij
    atom%rhoij(:,:) = atom%rhoij0(:,:)

    ! Calculate the compensation density based on initial guess rhoij0
    call atom_compden(atom,basis)

    ! Subtract the guessed compensation density off the guess valence density
    atom%den(:) = atom%den(:) - atom%comp_den(:)

    ! Calculate the initial guess local potential
    call atom_locpot(atom,basis)

  end subroutine atom_initialise_ham

  subroutine atom_update_ham(atom,basis,mix_alpha)

    !==========================================================================!
    ! This subroutine calculates the new density, mixes the new density with   !
    ! the old density and recalculates the new effective potential, to update  !
    ! the Kohn-Sham Hamiltonian for the atom.                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (input) : ATOM_TYPE object describing the atom.                   !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use services, only: services_regular_integral
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(ATOM_TYPE) :: atom
    type(BASIS_TYPE) :: basis
    real(kind=DP) :: mix_alpha

    ! Local Variables
    integer :: iorb
    integer :: ir
    real(kind=DP) :: den_int

    ! Calculate the new charge density in atom%work
    atom%work(:) = 0.0_DP
    do iorb=1,atom%norbs
       do ir=1,basis%npts
          atom%work(ir) = atom%work(ir) + atom%psi_r(ir,iorb)**2*atom%occ(iorb)
       end do
    end do

    ! Mix new ps-valence density with old ps-valence density
    atom%den(:) = atom%den(:)*(1.0_DP-mix_alpha) + atom%work(:)*mix_alpha

    ! Mix new rhoij with old rhoij
    atom%rhoij(:,:) = atom%rhoij0*(1.0_DP-mix_alpha) + atom%rhoij*mix_alpha

    call atom_compden(atom,basis)

    ! Check norm of new density, then rescale to make it exact
    den_int = services_regular_integral(basis%npts,basis%dr, &
         (atom%den(:)+atom%comp_den(:))*basis%rad(:)**2)
    if (abs(den_int-atom%val_charge)>0.00001_DP) then
       call utils_abort('Error in atom_update_ham: Incorrect valence charge &
            &calculated after wavefunction update')
    end if

    ! Update the local potential
    call atom_locpot(atom,basis)

  end subroutine atom_update_ham

  subroutine atom_locpot(atom,basis)

    !==========================================================================!
    ! This subroutine calculates the new effective potential, to update        !
    ! the Kohn-Sham Hamiltonian for the atom.                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (input) : ATOM_TYPE object describing the atom.                   !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in March 2011 out of bits of previous routines. !
    ! Modified for EMFT by Joseph Prentice, June 2018                          !
    !==========================================================================!

    use rundat, only: pub_paw, pub_nhat_in_xc, pub_turn_off_hartree, &
         pub_active_region, pub_emft
    use services, only: services_regular_integral
    use utils, only: utils_abort
    use xc, only: xc_embed_swap_functional

    implicit none

    ! Arguments
    type(ATOM_TYPE) :: atom
    type(BASIS_TYPE) :: basis

    ! Set the effective potential to the local pseudopotential
    atom%effpot(:) = atom%locpspot(:)

    ! Calculate the xc potential and add it to the effective potential
    ! jcap: if the current atom is in the active region, swap functionals
    if (pub_emft.and.(atom%ireg==pub_active_region)) then
       call xc_embed_swap_functional(.true.)
    end if
    if (pub_nhat_in_xc.or.(.not.pub_paw)) then
       ! Add in the compensation density if required
       atom%den(:) = atom%den(:) + atom%comp_den(:)
    end if
    call atom_xc(basis,atom%den,atom%core_den,atom%xcpot,atom%xc_energy, &
         atom%xc_dc_energy)
    atom%effpot(:) = atom%effpot(:) + atom%xcpot
    ! jcap: swap functional back
    if (pub_emft.and.(atom%ireg==pub_active_region)) then
       call xc_embed_swap_functional(.false.)
    end if

    ! Calculate the Hartree potential and add it to the effective potential
    if (.not.(pub_nhat_in_xc.or.(.not.pub_paw))) then
       ! Add in the compensation density if we did not already do so
       atom%den(:) = atom%den(:) + atom%comp_den(:)
    end if

    if (.not. pub_turn_off_hartree) then
       call atom_hartree(basis,atom%den,atom%work,atom%hartree_energy)
    else
       atom%work = 0.0_DP
       atom%hartree_energy = 0.0_DP
    endif

    call atom_hartree(basis,atom%den,atom%work,atom%hartree_energy)
    atom%effpot(:) = atom%effpot(:) + atom%work

    ! lpl: Compensated charge stuff used in ddec_mod
    !      Add comppot to locpspot if flagged to do so
    if(atom%include_rcomppot) then
       if(atom%update_rcomppot) then
          ! Get new qcomp due to total V(rcomp)
          atom%qcomp = (atom%effpot(atom%rcomp_ipt) &
               - atom%xcpot(atom%rcomp_ipt))*basis%rad(atom%rcomp_ipt)
       end if
       ! Add potential due to qcomp
       atom%effpot(:) = atom%effpot(:) - atom%qcomp*atom%rcomppot(:)
       atom%locpspot(:) = atom%locpspot(:) - atom%qcomp*atom%rcomppot(:)
    end if

    ! Subtract off the compensation density
    atom%den(:) = atom%den(:) - atom%comp_den(:)

    ! Find the energy of the density in the local potential
    ! In PAW, this is the pseudovalence and compensation density
    ! In USP/NCPP this is just the valence density
    if (.not.pub_paw) then
       atom%work = (atom%locpspot(:) + atom%confpot(:,0)) * &
            atom%den(:) * basis%rad(:)**2
    else if (pub_paw) then
       atom%work = (atom%locpspot(:) + atom%confpot(:,0)) * &
            (atom%den(:) + atom%comp_den(:)) * basis%rad(:)**2
    end if
    atom%locpsp_energy = services_regular_integral(basis%npts,basis%dr, &
         atom%work(:))

    ! lpl: Compensation charge stuff for ddec_mod
    !      Restore local pspot
    atom%locpspot(:) = atom%locpspot_orig(:)

  end subroutine atom_locpot

  subroutine atom_compden(atom,basis)

    !==========================================================================!
    ! This subroutine calculates the new compensation density, to update       !
    ! the Kohn-Sham Hamiltonian for the atom.                                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (input) : ATOM_TYPE object describing the atom.                   !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in May 2012 out of bits of previous routines.   !
    !==========================================================================!

    use rundat, only: pub_paw
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(ATOM_TYPE) :: atom
    type(BASIS_TYPE) :: basis

    ! Local variables
    integer :: ishell, jshell, il

    ! Calculate the compensation density
    atom%work(:) = 0.0_DP
    atom%comp_den_int = 0.0_DP
    if (pub_paw) then ! PAW version
       do ishell=1,atom%nshells
          do jshell=1,atom%nshells
             do il=0,0 !atom%proj_lmax
                atom%work(:) = atom%work(:) + atom%shape_l(:,il) * &
                     atom%qijL(ishell,jshell,il) * atom%rhoij(ishell,jshell)
                atom%comp_den_int = atom%comp_den_int + &
                    atom%qijL(ishell,jshell,il) * atom%rhoij(ishell,jshell)
             end do
          end do
       end do
    else ! USP version
       do ishell=1,atom%nshells
          do jshell=1,atom%nshells
             atom%work(:) = atom%work(:) + &
                  atom%qfunc(:,ishell,jshell) * atom%rhoij(ishell,jshell)
          end do
       end do
       atom%work(1) = atom%work(2) / basis%rad(2)**2
       atom%work(2:basis%npts) = atom%work(2:basis%npts) &
            / basis%rad(2:basis%npts)**2
    end if

    ! Set current compensation density
    atom%comp_den(:) = atom%work(:)

  end subroutine atom_compden

  subroutine atom_diagonalise(atom,basis)

    !==========================================================================!
    ! This subroutine creates and diagonalises the Hamiltonian for an isolated !
    ! atom in a Spherical Bessel function basis to find the Kohn-Sham orbitals !
    ! for that atom.                                                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (input) : ATOM_TYPE object describing the atom.                   !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use constants, only: paw_en_dijhat, paw_en_dij0, paw_en_ehart
    use linalg, only: linalg_dsygv_lt
    use rundat, only: pub_paw, pub_nhat_in_xc
    use services, only: services_regular_integral
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(ATOM_TYPE),intent(inout):: atom
    type(BASIS_TYPE),intent(in) :: basis

    ! Local Variables
    real(kind=DP), allocatable :: ham(:,:)
    real(kind=DP), allocatable :: overlap(:,:)
    real(kind=DP), allocatable :: eigs(:)
    real(kind=DP), allocatable :: work(:)
    integer :: il
    integer :: iorb
    integer :: iorb_l
    integer :: ibfn
    integer :: iproj,jproj
    integer :: ierr
    integer :: nbfn
    real(kind=DP), parameter :: tol = 1e-10_DP
    real(kind=DP) :: int_veff, int_vHxc

    ! BLAS subroutine
    external :: dgemv

    ! Allocate full basis-sized matrices for ham, overlap and eigs
    allocate(ham(basis%nbfn_max,basis%nbfn_max),stat=ierr)
    call utils_alloc_check('atom_diagonalise','ham',ierr)
    allocate(overlap(basis%nbfn_max,basis%nbfn_max),stat=ierr)
    call utils_alloc_check('atom_diagonalise','overlap',ierr)
    allocate(eigs(basis%nbfn_max),stat=ierr)
    call utils_alloc_check('atom_diagonalise','eigs',ierr)
    allocate(work(6*basis%nbfn_max),stat=ierr)
    call utils_alloc_check('atom_diagonalise','work',ierr)

    ! Set the nonlocal energies

    ! Calculate the PAW nonlocal energy terms if required
    if (pub_paw) then

       atom%paw_sphere_energies(paw_en_dij0) = 0.0_DP
       do iproj=1,atom%nshells
          do jproj=1,atom%nshells
             atom%paw_sphere_energies(paw_en_dij0) = &
                  atom%paw_sphere_energies(paw_en_dij0) + &
                  atom%dij0(iproj,jproj) * atom%rhoij(iproj,jproj)
          end do
       end do

       ! Find the dijhat terms \sum_L \int\tilde{v}_eff(r)\hat{Q}_ij^L (r) dr
       atom%dijhat = 0.0_DP
       atom%dij_hfco = 0.0_DP
       atom%paw_sphere_energies(paw_en_dijhat) = 0.0_DP
       do il=0,0!atom%proj_lmax
          ! Screened nonlocal energies
          if (pub_nhat_in_xc) then
             atom%work = atom%shape_l(:,il) * atom%effpot(:) * basis%rad(:)**2
          else
             atom%work = atom%shape_l(:,il) * (atom%effpot(:) - atom%xcpot(:)) &
                  * basis%rad(:)**2
          end if
          int_veff = services_regular_integral(basis%npts,basis%dr, &
               atom%work)
          ! for HF corrective term in nonlocal energies
          atom%work = atom%shape_l(:,il) * (atom%effpot(:) - atom%locpspot(:)) &
               * basis%rad(:)**2
          int_vHxc = services_regular_integral(basis%npts,basis%dr, &
               atom%work)
          do iproj=1,atom%nshells
             do jproj=1,atom%nshells
               atom%dijhat(iproj,jproj) = atom%dijhat(iproj,jproj) + &
                    int_veff * atom%qijL(iproj,jproj,il)
               atom%dij_hfco(iproj,jproj) = atom%dij_hfco(iproj,jproj) + &
                    int_vHxc * atom%qijL(iproj,jproj,il)
               atom%paw_sphere_energies(paw_en_dijhat) = &
                    atom%paw_sphere_energies(paw_en_dijhat) + &
                    atom%dijhat(iproj,jproj) * atom%rhoij(iproj,jproj)
            end do
         end do
       end do

       ! Add the dijhartree terms if in PAW
       call atom_paw_dij_hartree(atom)
       do iproj=1,atom%nshells
          do jproj=1,atom%nshells
             atom%paw_sphere_energies(paw_en_ehart) = &
                  atom%paw_sphere_energies(paw_en_ehart) + &
                  0.5_DP * atom%dij_hart(iproj,jproj) * atom%rhoij(iproj,jproj)
          end do
       end do

       ! Add the dijxc terms if in PAW
       call atom_paw_dij_xc(atom)

       ! Add up the total dij
       atom%dij = atom%dij0 + atom%dij_hart + atom%dij_xc + atom%dijhat

    else

       atom%dij(:,:) = atom%dij0(:,:)
       atom%dij_hfco(:,:) = 0.0_DP

       do iproj=1,atom%nshells
          do jproj=1,atom%nshells
             if (atom%proj_ang_mom(iproj)/=atom%proj_ang_mom(jproj)) cycle
             ! Screened nonlocal energies
             atom%work(:) = atom%qfunc(:,iproj,jproj)*atom%effpot(:)
             int_veff = services_regular_integral(basis%npts,basis%dr, &
                  atom%work)
             atom%dij(iproj,jproj) = atom%dij(iproj,jproj) + int_veff
             ! For Harris-Foulkes total energy
             atom%work(:) = atom%qfunc(:,iproj,jproj) * (atom%effpot(:) - &
                  atom%locpspot(:))
             int_vHxc = services_regular_integral(basis%npts,basis%dr, &
                  atom%work)
             atom%dij_hfco(iproj,jproj) = atom%dij_hfco(iproj,jproj) + int_vHxc
          end do
       end do

    end if

    atom%psi_cin(:,:) = 0.0_DP
    atom%psi_r(:,:) = 0.0_DP

    ! Loop over angular momentum channels
    do il=0,atom%orb_lmax

       if (.not.any(atom%orb_ang_mom(:)==il)) cycle

       ! Get Hamiltonian and overlap for this angular momentum
       call atom_ham_matrix(ham,overlap,il,atom,basis)

#if 0
    nbfn = basis%nbfn(il)
    print *,'overlap: ', il, nbfn
    do ibfn=1,nbfn
       do jbfn=1,nbfn
          write(stdout,'(f18.12)',advance='no') overlap(jbfn,ibfn)
       end do
       write(stdout,*)
    end do
    print *,'ham: ', il, nbfn
    do ibfn=1,nbfn
       do jbfn=1,nbfn
          if (jbfn>ibfn) then
             write(stdout,'(f18.12)',advance='no') ham(ibfn,jbfn)
          else
             write(stdout,'(f18.12)',advance='no') ham(jbfn,ibfn)
          end if
       end do
       write(stdout,*)
    end do
#endif

       ! Solve generalised eigenvalue problem
       eigs = 0.0_DP
       work = 0.0_DP
       nbfn = basis%nbfn(il)

       call linalg_dsygv_lt(ham,eigs,overlap,nbfn,ut=.true.)

       ! Copy the eigenvalues for this l into the relevant parts of the
       ! psi_cin array, then calculate psi_r in real space
       iorb_l = 0
       do iorb=1,atom%norbs
          if (atom%orb_ang_mom(iorb) == il) then
             iorb_l = iorb_l + 1
             atom%eigs(iorb) = eigs(iorb_l)
             atom%psi_cin(1:nbfn,iorb) = ham(1:nbfn,iorb_l)

             ! Use BLAS matrix-vector routine on matrix of evaluated Bessel
             ! functions on real space grid to transform psi to real space
             call dgemv('N',basis%npts,nbfn,1.0_DP,basis%bfn(1,1,il), &
                  basis%npts,atom%psi_cin(1,iorb),1,0.0_DP,atom%psi_r(1,iorb),1)

             ! If an eigenstate is -ve on more points than it is +ve
             ! multiply it by -1 for consistency between architectures
             if (count(atom%psi_r(1:(basis%npts),iorb)<0.0_DP) .gt. &
                  count(atom%psi_r(1:(basis%npts),iorb)>0.0_DP)) then
                atom%psi_r(:,iorb) = -1.0_DP * atom%psi_r(:,iorb)
             end if

          end if
       end do

    end do

    ! Deallocate full basis-sized matrices for ham, overlap and eigs
    deallocate(work,stat=ierr)
    call utils_dealloc_check('atom_diagonalise','work',ierr)
    deallocate(eigs,stat=ierr)
    call utils_dealloc_check('atom_diagonalise','eigs',ierr)
    deallocate(overlap,stat=ierr)
    call utils_dealloc_check('atom_diagonalise','overlap',ierr)
    deallocate(ham,stat=ierr)
    call utils_dealloc_check('atom_diagonalise','ham',ierr)

    ! Calculate the bandstructure energy
    atom%eigs(:) = atom%eigs(:)
    atom%band_energy = sum(atom%eigs(:)*atom%occ(:))

    ! Calculate the kinetic energy
    atom%kinetic_energy = 0.0_DP
    do iorb=1,atom%norbs
       il = atom%orb_ang_mom(iorb)
       nbfn = basis%nbfn(il)
       do ibfn=1,nbfn
          atom%kinetic_energy = atom%kinetic_energy + &
               atom%occ(iorb) * 0.5_DP * atom%psi_cin(ibfn,iorb)**2 &
               * basis%qb(ibfn,il)**2 * basis%bcoef(1,ibfn,il)**2
          if (basis%nfixderiv>0) then
             atom%kinetic_energy = atom%kinetic_energy + &
                  atom%occ(iorb) * 0.5_DP * atom%psi_cin(ibfn,iorb)**2 &
                  * basis%qb(ibfn+1,il)**2 * basis%bcoef(2,ibfn,il)**2
             if (ibfn>1) atom%kinetic_energy = atom%kinetic_energy + &
                  atom%occ(iorb) * 0.5_DP * atom%psi_cin(ibfn,iorb) &
                  * atom%psi_cin(ibfn-1,iorb) * basis%qb(ibfn,il)**2 &
                  * basis%bcoef(1,ibfn,il) * basis%bcoef(2,ibfn-1,il)
             if (ibfn<nbfn) atom%kinetic_energy = atom%kinetic_energy + &
                  atom%occ(iorb) * 0.5_DP * atom%psi_cin(ibfn,iorb) &
                  * atom%psi_cin(ibfn+1,iorb) * basis%qb(ibfn+1,il)**2  &
                  * basis%bcoef(2,ibfn,il) * basis%bcoef(1,ibfn+1,il)
          end if
       end do
    end do

    ! Recalculate the rhoij matrix
    call atom_rhoij(atom,basis)

    ! Calculate the nonlocal energy using the unscreened nonlocal energies
    atom%nonlocal_energy = 0.0_DP
    do iproj=1,atom%nshells
       do jproj=1,atom%nshells
          atom%nonlocal_energy = atom%nonlocal_energy + &
               atom%rhoij(iproj,jproj)*(atom%dij(iproj,jproj) - &
               atom%dij_hfco(iproj,jproj))
       end do
    end do

  end subroutine atom_diagonalise

  subroutine atom_ham_matrix(ham,overlap,il,atom,basis)

    !==========================================================================!
    ! This subroutine evaluates the matrix elements of the Hamiltonian and     !
    ! overlap matrix for the atom, for a given angular momentum channel il.    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   ham (output) : Hamiltonian matrix for the basis functions of this il.  !
    !   overlap (output) : overlap matrix for the basis functions of this il.  !
    !   il (input) : current angular momentum channel il.                      !
    !   atom (input) : ATOM_TYPE object describing the atom.                   !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    ! Moved to standalone routine in January 2012.                             !
    !==========================================================================!

    use services, only: services_regular_integral
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(ATOM_TYPE), intent(inout) :: atom
    type(BASIS_TYPE), intent(in) :: basis
    real(kind=DP), intent(inout) :: ham(basis%nbfn_max,basis%nbfn_max)
    real(kind=DP), intent(inout) :: overlap(basis%nbfn_max,basis%nbfn_max)
    integer, intent(in) :: il

    ! Local Variables
    integer :: nbfn,ibfn,jbfn
    integer :: iproj,jproj

    nbfn = basis%nbfn(il)
    ham(:,:) = 0.0_DP

    ! Kinetic Energy
    do ibfn=1,nbfn
       do jbfn=ibfn,nbfn
          if (basis%nfixderiv==0) then
             if (ibfn==jbfn) ham(ibfn,jbfn) = 0.5_DP*basis%qb(ibfn,il)**2
          else
             if (ibfn==jbfn) ham(ibfn,jbfn) = &
                  0.5_DP*basis%bcoef(1,ibfn,il)**2 * basis%qb(ibfn,il)**2 &
                  + 0.5_DP*basis%bcoef(2,ibfn,il)**2 * basis%qb(ibfn+1,il)**2
             if (ibfn==jbfn+1) ham(ibfn,jbfn) = &
                  0.5_DP*basis%bcoef(1,ibfn,il) * basis%bcoef(2,jbfn,il) * &
                  basis%qb(jbfn,il)**2
             if (ibfn+1==jbfn) ham(ibfn,jbfn) = &
                  0.5_DP*basis%bcoef(2,ibfn,il) * basis%bcoef(1,jbfn,il) * &
                  basis%qb(jbfn,il)**2
          end if
       end do
    end do

    ! Add the confining potential for this angular momentum
    atom%effpot(:) = atom%effpot(:) + atom%confpot(:,min(il,5))

    ! Local Potential
    do ibfn=1,nbfn
       atom%work(:) = atom%effpot(:) * basis%rad(:)**2 * basis%bfn(:,ibfn,il)
       do jbfn=ibfn,nbfn
          ham(ibfn,jbfn) = ham(ibfn,jbfn) + &
               services_regular_integral(basis%npts,basis%dr, &
               atom%work(:)*basis%bfn(:,jbfn,il))
       end do
    end do

    ! Subtract the confining potential for this angular momentum
    atom%effpot(:) = atom%effpot(:) - atom%confpot(:,min(il,5))

    ! Nonlocal Potential
    do ibfn=1,nbfn
       do jbfn=ibfn,nbfn
          do iproj=1,atom%nshells
             if (atom%proj_ang_mom(iproj)/=il) cycle
             do jproj=1,atom%nshells
                if (atom%proj_ang_mom(jproj)/=il) cycle
                ham(ibfn,jbfn) = ham(ibfn,jbfn) + atom%proj_q(ibfn,iproj) * &
                     atom%dij(iproj,jproj) * atom%proj_q(jbfn,jproj)
             end do
          end do
       end do
    end do

    ! Overlap matrix
    overlap(:,:) = 0.0_DP
    do ibfn=1,nbfn
       do jbfn=ibfn,nbfn
          if (ibfn==jbfn) overlap(ibfn,jbfn) = basis%bcoef(1,ibfn,il)**2 &
               + basis%bcoef(2,ibfn,il)**2
          if (basis%nfixderiv>0) then
             if (ibfn+1==jbfn) overlap(ibfn,jbfn) = basis%bcoef(2,ibfn,il) &
                  * basis%bcoef(1,jbfn,il)
             if (ibfn==jbfn+1) overlap(ibfn,jbfn) = basis%bcoef(1,ibfn,il) &
                  * basis%bcoef(2,jbfn,il)
          end if
       end do
    end do

    ! Augmentation of Overlap matrix
    do ibfn=1,nbfn
       do jbfn=1,nbfn
          do iproj=1,atom%nshells
             if (atom%proj_ang_mom(iproj)/=il) cycle
             do jproj=1,atom%nshells
                if (atom%proj_ang_mom(jproj)/=il) cycle
                overlap(ibfn,jbfn) = overlap(ibfn,jbfn) + atom%proj_q(ibfn,iproj) * &
                     atom%qijL(iproj,jproj,0) * atom%proj_q(jbfn,jproj)
             end do
          end do
       end do
    end do

    ! Symmetrise Hamiltonian matrix to ensure both LT & UT methods work
    do ibfn=1,nbfn
       do jbfn=1,ibfn-1
          ham(ibfn,jbfn) = ham(jbfn,ibfn)
          overlap(ibfn,jbfn) = overlap(jbfn,ibfn)
       end do
    end do

  end subroutine atom_ham_matrix

  subroutine atom_rhoij(atom,basis)

    !==========================================================================!
    ! This subroutine calculates the projector density kernel rhoij for an atom!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   basis (input) : BASIS_TYPE describing grid.                            !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in May 2012.                                    !
    !==========================================================================!

    implicit none

    ! Arguments
    type(ATOM_TYPE), intent(inout) :: atom
    type(BASIS_TYPE), intent(in) :: basis

    ! Local Variables
    integer :: iproj, jproj, nbfn, iorb
    real(kind=DP), external :: ddot

    ! Recalculate the rhoij matrix
    atom%rhoij0(:,:) = atom%rhoij(:,:)
    atom%rhoij(:,:) = 0.0_DP
    do iproj=1,atom%nshells
       do jproj=1,atom%nshells
          do iorb=1,atom%norbs
             nbfn = basis%nbfn(atom%orb_ang_mom(iorb))
             if ((atom%orb_ang_mom(iorb)==atom%proj_ang_mom(iproj)).and. &
                  (atom%orb_ang_mom(iorb)==atom%proj_ang_mom(jproj))) then
                atom%rhoij(iproj,jproj) = atom%rhoij(iproj,jproj) + &
                     ddot(nbfn,atom%psi_cin(1,iorb),1, &
                     atom%proj_q(1,iproj),1) * &
                     ddot(nbfn,atom%psi_cin(1,iorb),1, &
                     atom%proj_q(1,jproj),1) * &
                     atom%occ(iorb)
             end if
          end do
       end do
    end do

  end subroutine atom_rhoij

  subroutine atom_hartree(basis,den,hartree_pot,hartree_energy)

    !==========================================================================!
    ! This subroutine calculates the Hartree potential and Hartree energy for  !
    ! an atom, from the density on a regular real space grid.                  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   basis (input) : BASIS_TYPE describing grid.                            !
    !   den (input) : density on real space grid                               !
    !   hartree_pot (out) : Hartree potential of density on real space grid.   !
    !   hartree_energy (out) : Hartree energy of density.                      !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use services, only: services_regular_integral

    implicit none

    ! Arguments
    type(BASIS_TYPE) :: basis
    real(kind=DP) :: den(basis%npts)
    real(kind=DP) :: hartree_pot(basis%npts)
    real(kind=DP) :: hartree_energy

    ! Local Variables
    integer :: ipt
    real(kind=DP) :: vp
    real(kind=DP) :: vm

    ! Calculate Hartree potential
    vp = 0.0_DP
    vm = 0.0_DP
    hartree_pot(1) = 0.0_DP
    do ipt=2,basis%npts
       vp = vp + 2.0_DP*den(ipt)*basis%rad(ipt)**2*basis%dr
       vm = vm + 2.0_DP*den(ipt)*basis%rad(ipt)*basis%dr
       hartree_pot(ipt) = vp/basis%rad(ipt) - vm
    end do
    hartree_pot = hartree_pot + vm
    hartree_pot = 0.5_DP*hartree_pot

    ! Calculate Hartree energy
    hartree_energy = services_regular_integral(basis%npts,basis%dr, &
         hartree_pot(:)*den(:)*basis%rad(:)**2)
    hartree_energy = 0.5_DP*hartree_energy

  end subroutine atom_hartree

  subroutine atom_xc(basis,den,core_den,xc_pot,xc_energy,xc_dc_energy)

    !==========================================================================!
    ! This subroutine calculates the XC potential and XC energy for an atom,   !
    ! from the density on a regular real space grid.                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   basis (input)   : BASIS_TYPE describing grid.                          !
    !   den (input)     : density on real space grid                           !
    !   xc_pot (out)    : XC potential of density on real space grid.          !
    !   xc_energy (out) : XC energy of density.                                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use constants, only: PI
    use services, only: services_regular_integral, services_radial_derivative
    use xc, only: xc_radial
#if 0
    use xc, only: pub_xc_gradient_corrected
#endif
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(BASIS_TYPE),intent(in) :: basis
    real(kind=DP),intent(in) :: den(basis%npts)
    real(kind=DP),intent(in) :: core_den(basis%npts)
    real(kind=DP),intent(out) :: xc_pot(basis%npts)
    real(kind=DP),intent(out) :: xc_energy
    real(kind=DP),intent(out) :: xc_dc_energy

    ! Local Variables
    integer :: ierr
    real(kind=DP), allocatable :: den_pt(:)
    real(kind=DP), allocatable :: den_grad(:)
    real(kind=DP), allocatable :: dpot_dden(:)
    real(kind=DP), allocatable :: exc(:)

    allocate(den_pt(basis%npts),stat=ierr)
    call utils_alloc_check('atom_xc','den_pt',ierr)
    allocate(den_grad(basis%npts),stat=ierr)
    call utils_alloc_check('atom_xc','den_grad',ierr)
    allocate(dpot_dden(basis%npts),stat=ierr)
    call utils_alloc_check('atom_xc','dpot_dden',ierr)
    allocate(exc(basis%npts),stat=ierr)
    call utils_alloc_check('atom_xc','exc',ierr)

    ! Remove solid angle integration factor before passing density to xc_mod
    den_pt = (den+core_den) / (4.0_DP*PI)

#if 0
    if (pub_xc_gradient_corrected) then
       call services_radial_derivative(den_grad,den_pt,basis%npts,basis%rmax)
    else
       den_grad = 0.0_DP
    end if
#endif

    ! Calculate v_xc and e_xc
    call xc_radial(basis%npts,basis%npts,1,den_pt,den_grad,xc_pot,dpot_dden, &
         exc,rad=basis%rad,rmax=basis%rmax)

#if 0
    if (pub_xc_gradient_corrected) then
       do ipt=1,basis%npts
          if (abs(den_grad(ipt))>0.0_DP) then
             dpot_dden(ipt) = -dpot_dden(ipt) * den_grad(ipt) &
                  / abs(den_grad(ipt))
          else
             dpot_dden(ipt) = 0.0_DP
          end if
       end do

       call services_radial_derivative(den_grad,dpot_dden,basis%npts,basis%rmax)
       den_grad(2:) = den_grad(2:) + 2.0_DP*dpot_dden(2:)/basis%rad(2:)
       xc_pot(:) = xc_pot(:) + den_grad(:)

    end if
#endif

    ! Calculate XC energy
    xc_energy = services_regular_integral(basis%npts,basis%dr, &
         exc(:)*basis%rad(:)**2)*4.0_DP*PI
    xc_dc_energy = xc_energy - services_regular_integral(basis%npts,basis%dr, &
         xc_pot(:)*den(:)*basis%rad(:)**2)

    deallocate(exc,stat=ierr)
    call utils_dealloc_check('atom_xc','exc',ierr)
    deallocate(dpot_dden,stat=ierr)
    call utils_dealloc_check('atom_xc','dpot_dden',ierr)
    deallocate(den_grad,stat=ierr)
    call utils_dealloc_check('atom_xc','den_grad',ierr)
    deallocate(den_pt,stat=ierr)
    call utils_dealloc_check('atom_xc','den_pt',ierr)

  end subroutine atom_xc

  subroutine internal_derivative(grad,func,npts,xmax)

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP),intent(in) :: xmax
    real(kind=DP),intent(in) :: func(:)
    real(kind=DP),intent(out) :: grad(:)

    ! Local Variables
    integer :: ipt
    real(kind=DP) :: t1,t2,t3

    ! Gradient is zero at origin by symmetry
    grad(1) = 0.0_DP
    t1 = 0.0_DP; t2 = 0.0_DP; t3 = 0.0_DP

    ! Cubic interpolation
    do ipt=2,npts-2
       t1 = (6.0_DP*func(ipt+1) - 2.0_DP*func(ipt-1) - 3.0_DP*func(ipt) &
            - 1.0_DP*func(ipt+2)) / 6.0_DP
       t2 = (1.0_DP*func(ipt-1) + 1.0_DP*func(ipt+1) - 2.0_DP*func(ipt)) &
            / 2.0_DP
       t3 = (1.0_DP*func(ipt+2) - 1.0_DP*func(ipt-1) + 3.0_DP*func(ipt) &
            - 3.0_DP*func(ipt+1)) / 6.0_DP
       grad(ipt) = t1
    end do

    ! Last two points
    grad(npts-1) = t1 + 2.0_DP*t2 + 3.0_DP*t3
    grad(npts) = t1 + 4.0_DP*t2 + 12.0_DP*t3

    ! Normalise for dr
    grad(:) = grad(:)*real(npts-1,kind=DP)/xmax

  end subroutine internal_derivative


  subroutine atom_paw_dij_hartree(atom)

    !=================================================================!
    ! This subroutine calculates the XC terms in the nonlocal         !
    ! energies Dij for a PAW atom, and also the XC energies and       !
    ! double-counting terms.                                          !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 23/09/10.                           !
    !=================================================================!

    use constants, only: paw_en_ehart
    use paw, only: paw_dij_hartree_atom
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(ATOM_TYPE), intent(inout) :: atom

    ! Local Variables
    real(kind=DP), allocatable :: rhoij_at(:,:)
    real(kind=DP), allocatable :: dijhartree_at(:,:)
    integer :: ierr
    integer :: ipwtot,jpwtot
    integer :: ipw,jpw
    integer :: li,lj,di
    integer :: npw,npwtot

    ! Find array sizes
    npwtot = atom%npwtot
    npw = atom%nshells

    ! Allocate temporary arrays
    allocate(rhoij_at(npwtot,npwtot),stat=ierr)
    call utils_alloc_check('atom_paw_dij_hartree','rhoij_at',ierr)
    allocate(dijhartree_at(npwtot,npwtot),stat=ierr)
    call utils_alloc_check('atom_paw_dij_hartree','dijhartree_at',ierr)

    ! Expand rhoij (does not include m-degeneracy) to rhoij_at (does)
    ipwtot = 1
    rhoij_at(:,:) = 0.0_DP
    do ipw=1,npw
       li = atom%proj_ang_mom(ipw)
       jpwtot = 1
       do jpw=1,npw
          lj = atom%proj_ang_mom(jpw)
          if (li==lj) then
             do di=0,2*li
                rhoij_at(ipwtot+di,jpwtot+di) = &
                     rhoij_at(ipwtot+di,jpwtot+di) + &
                     atom%rhoij(ipw,jpw) / real(2*li+1,kind=DP)
             end do
          end if
          jpwtot = jpwtot + 2*lj + 1
       end do
       ipwtot = ipwtot + 2*li + 1
    end do

    ! Calculate exchange-correlation energy and potential for this atom
    call paw_dij_hartree_atom(atom%pspecies_number,npwtot,rhoij_at, &
         dijhartree_at,atom%paw_sp)

    ! Add this to the nonlocal energies
    ipwtot = 1
    atom%paw_sphere_energies(paw_en_ehart) = 0.0_DP
    atom%dij_hart = 0.0_DP
    do ipw=1,npw
       li = atom%proj_ang_mom(ipw)
       jpwtot = 1
       do jpw=1,npw
          lj = atom%proj_ang_mom(jpw)
          if (li==lj) then
             do di=0,2*li
                atom%dij_hart(ipw,jpw) = atom%dij_hart(ipw,jpw) + &
                     dijhartree_at(ipwtot+di,jpwtot+di) / real(2*li+1,kind=DP)
             end do
          end if
          jpwtot = jpwtot + 2*lj + 1
       end do
       ipwtot = ipwtot + 2*li + 1
    end do

    ! Deallocate temporary arrays
    deallocate(dijhartree_at,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_hartree','dijhartree_at',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_hartree','rhoij_at',ierr)

  end subroutine atom_paw_dij_hartree


  subroutine atom_paw_dij_xc(atom)

    !=================================================================!
    ! This subroutine calculates the XC terms in the nonlocal         !
    ! energies Dij for a PAW atom, and also the XC energies and       !
    ! double-counting terms.                                          !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 23/09/10.                           !
    ! Modified for EMFT by Joseph Prentice, June 2018                 !
    !=================================================================!

    use constants, only: max_spins, paw_en_exc, paw_en_exc_dc, &
         paw_en_etxc, paw_en_etxc_dc, paw_en_dijxc
    use paw, only: paw_dij_xc_atom
    use rundat, only: pub_emft, pub_active_region
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check
    use xc, only: xc_embed_swap_functional

    implicit none

    ! Arguments
    type(ATOM_TYPE), intent(inout) :: atom
    real(kind=DP) :: exc
    real(kind=DP) :: exc_dc
    real(kind=DP) :: etxc
    real(kind=DP) :: etxc_dc

    ! Local Variables
    real(kind=DP), allocatable :: vxc_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: density_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: rhoij_at(:,:,:)
    real(kind=DP), allocatable :: dijxc_at(:,:,:)
    real(kind=DP), allocatable :: dijtxc_at(:,:,:)
    real(kind=DP), allocatable :: pot_work(:,:,:,:)
    real(kind=DP), allocatable :: den_work(:,:,:,:)
    real(kind=DP), allocatable :: inter(:)
    integer :: ierr
    integer :: is
    integer :: ipwtot,jpwtot
    integer :: ipw,jpw
    integer :: li,lj,di
    integer :: npw,npwtot
    integer :: npts_max
    integer :: nLM_max
    integer :: nspins
    real(kind=DP) :: edijxc
    real(kind=DP) :: edijtxc
    real(kind=DP) :: total_nhat(max_spins)

    ! Find array sizes
    npts_max = atom%log_npts_max
    nLM_max = (atom%proj_lmax+1)**2
    nspins = 1
    npwtot = atom%npwtot
    npw = atom%nshells

    ! Allocate temporary arrays
    allocate(vxc_rad_LM(npts_max,nspins,nLM_max,2),stat=ierr)
    call utils_alloc_check('atom_paw_dij_xc','vxc_rad_LM',ierr)
    allocate(density_rad_LM(npts_max,nspins,nLM_max,3),stat=ierr)
    call utils_alloc_check('atom_paw_dij_xc','density_rad_LM',ierr)
    allocate(rhoij_at(npwtot,npwtot,nspins),stat=ierr)
    call utils_alloc_check('atom_paw_dij_xc','rhoij_at',ierr)
    allocate(dijxc_at(npwtot,npwtot,nspins),stat=ierr)
    call utils_alloc_check('atom_paw_dij_xc','dijxc_at',ierr)
    allocate(dijtxc_at(npwtot,npwtot,nspins),stat=ierr)
    call utils_alloc_check('atom_paw_dij_xc','dijtxc_at',ierr)
    allocate(pot_work(npts_max,nspins,nspins,6),stat=ierr)
    call utils_alloc_check('atom_paw_dij_xc','pot_work',ierr)
    allocate(den_work(npts_max,nspins,nspins,4),stat=ierr)
    call utils_alloc_check('atom_paw_dij_xc','den_work',ierr)
    allocate(inter(npts_max),stat=ierr)
    call utils_alloc_check('atom_paw_dij_xc','inter',ierr)

    ! Initialisations
    exc = 0.0_DP
    exc_dc = 0.0_DP
    etxc = 0.0_DP
    etxc_dc = 0.0_DP
    dijxc_at = 0.0_DP
    dijtxc_at = 0.0_DP

    ! Expand rhoij (does not include m-degeneracy) to rhoij_at (does)
    ipwtot = 1
    rhoij_at(:,:,:) = 0.0_DP
    do ipw=1,npw
       li = atom%proj_ang_mom(ipw)
       jpwtot = 1
       do jpw=1,npw
          lj = atom%proj_ang_mom(jpw)
          if (li==lj) then
             do di=0,2*li
                rhoij_at(ipwtot+di,jpwtot+di,1) = &
                     rhoij_at(ipwtot+di,jpwtot+di,1) + &
                     atom%rhoij(ipw,jpw) / real(2*li+1,kind=DP)
             end do
          end if
          jpwtot = jpwtot + 2*lj + 1
       end do
       ipwtot = ipwtot + 2*li + 1
    end do

    ! Calculate exchange-correlation energy and potential for this atom
    ! jcap: if the atom is in the active region, swap functionals
    if (pub_emft.and.(atom%ireg==pub_active_region)) then
       call xc_embed_swap_functional(.true.)
    end if
    call paw_dij_xc_atom(atom%pspecies_number,nspins,npts_max,nLM_max,npwtot, &
         rhoij_at,dijxc_at,dijtxc_at,exc,etxc,exc_dc,etxc_dc,total_nhat, &
         density_rad_LM,vxc_rad_LM,den_work,pot_work,inter,atom%paw_sp)
    ! jcap: swap back
    if (pub_emft.and.(atom%ireg==pub_active_region)) then
       call xc_embed_swap_functional(.false.)
    end if

    ! Add this to the nonlocal energy
    ipwtot = 1
    atom%dij_xc = 0.0_DP
    do ipw=1,npw
       li = atom%proj_ang_mom(ipw)
       jpwtot = 1
       do jpw=1,npw
          lj = atom%proj_ang_mom(jpw)
          if (li==lj) then
             do di=0,2*li
                atom%dij_xc(ipw,jpw) = atom%dij_xc(ipw,jpw) + &
                     (dijxc_at(ipwtot+di,jpwtot+di,1) - &
                     dijtxc_at(ipwtot+di,jpwtot+di,1)) / real(2*li+1,kind=DP)
             end do
          end if
          jpwtot = jpwtot + 2*lj + 1
       end do
       ipwtot = ipwtot + 2*li + 1
    end do

    atom%paw_sphere_energies(paw_en_exc) = exc
    atom%paw_sphere_energies(paw_en_etxc) = etxc
    atom%paw_sphere_energies(paw_en_exc_dc) = exc_dc
    atom%paw_sphere_energies(paw_en_etxc_dc) = etxc_dc

    edijxc = 0.0_DP
    edijtxc = 0.0_DP
    do is=1,nspins
       do ipwtot=1,npwtot
          do jpwtot=1,npwtot
             edijxc = edijxc + rhoij_at(ipwtot,jpwtot,is) * &
                  dijxc_at(ipwtot,jpwtot,is)
             edijtxc = edijtxc + rhoij_at(ipwtot,jpwtot,is) * &
                  dijtxc_at(ipwtot,jpwtot,is)
          end do
       end do
    end do
    atom%paw_sphere_energies(paw_en_dijxc) = edijxc - edijtxc

    ! Consistency Check
    if (abs(edijxc-exc_dc)>0.00001*abs(exc_dc)) then
       write(stdout,'(a,f20.12)') 'edijxc: ',edijxc
       write(stdout,'(a,f20.12)') 'exc_dc: ',exc_dc
       call utils_abort('Error in atom_paw_dij_xc: Consistency check failed: &
            &edijxc /= exc_dc')
    end if
    if (abs(edijtxc-etxc_dc)>0.00001_DP*abs(etxc_dc)) then
       write(stdout,'(a,f20.12)') 'edijtxc: ',edijtxc
       write(stdout,'(a,f20.12)') 'etxc_dc: ',etxc_dc
       call utils_abort('Error in atom_paw_dij_xc: Consistency check failed: &
            &edijtxc /= etxc_dc')
    end if

    ! Deallocate temporary arrays
    deallocate(inter,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_xc','inter',ierr)
    deallocate(den_work,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_xc','den_work',ierr)
    deallocate(pot_work,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_xc','pot_work',ierr)
    deallocate(dijtxc_at,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_xc','dijtxc_at',ierr)
    deallocate(dijxc_at,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_xc','dijxc_at',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_xc','rhoij_at',ierr)
    deallocate(density_rad_LM,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_xc','density_rad_LM',ierr)
    deallocate(vxc_rad_LM,stat=ierr)
    call utils_dealloc_check('atom_paw_dij_xc','vxc_rad_LM',ierr)

  end subroutine atom_paw_dij_xc

  subroutine atom_split_orbital(func,func_orig,nsplit,splitnorm,ang_mom,work, &
       basis)

    !==========================================================================!
    ! This subroutine splits a function to produce nsplitnorm different        !
    ! orbitals in the "split-valence" approach to multiple-zeta basis sets.    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   func (output) : New radial functions on grid defined by basis%rad      !
    !   func_orig (input) : Orig radial functions on grid defined by basis%rad !
    !   nsplit (input) : Number of functions to split original function into   !
    !   splitnorm (input) : Norm beyond matching radius for each split-off func!
    !   ang_mom (input) : Angular momentum of original (and new) functions     !
    !   work (output) : workspace array of size basis%npts                     !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in August 2011.                                 !
    ! Modified by Chris-Kriton Skylaris in December 2012 to change logic of    !
    ! multiple-zeta splitting.                                                 !
    !==========================================================================!

    use services, only: services_regular_integral, services_radial_derivative
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(BASIS_TYPE),intent(in) :: basis
    real(kind=DP), intent(out) :: func(:,:)
    real(kind=DP), intent(in) :: func_orig(1:basis%npts)
    real(kind=DP), intent(in) :: splitnorm(4)
    real(kind=DP), intent(inout) :: work(1:basis%npts)
    integer, intent(in) :: nsplit
    integer, intent(in) :: ang_mom

    ! Local Variables
    integer :: isplit
    integer :: ipt, iptm
    integer :: l
    real(kind=DP) :: rm, fm, fmp, fmc, fmpc, norm, norm_tot, al, bl
    real(kind=DP), parameter :: tol=1e-10_DP

    ! Copy original orbital into func(:,1)
    func(1:basis%npts,1) = func_orig(1:basis%npts)

    ! Now loop over the others that are to be split off
    do isplit=2,nsplit

       ! Workspace array holds original |psi(r)|^2 * r^2
       work = func(1:basis%npts,1)**2 * basis%rad(1:basis%npts)**2
       norm = 0.0_DP
       rm = 0.0_DP

       ! First calculate total norm
       norm_tot = 0.0_DP
       do ipt=1,basis%npts
          norm_tot = norm_tot + work(ipt)*basis%dr
       end do

       ! Find matching radius where remaining norm = splitnorm
       iptm = -999
       do ipt=1,basis%npts
          norm = norm + work(ipt)*basis%dr
          if (norm > (norm_tot - splitnorm(isplit-1))) then
             rm = basis%rad(ipt)
             iptm = ipt
             exit
          end if
       end do
       if (iptm == -999) then
          call utils_abort('Error in atom_split_orbital: Attempt to split &
               &function failed when finding matching radius')
       end if

       ! Find f(rm) and (df/dr)(rm) for original function
       call services_radial_derivative(work,func(:,1),basis%npts, &
            basis%rmax)
       fm = func(iptm,1)
       fmp = work(iptm)

       ! Solve for coefficients al,bl
       l = ang_mom
       al = fm / rm**l + l*fm/(2.0_DP*rm**l) - fmp/(2.0_DP*rm**(l-1))
       bl = l*fm/(2*rm**(l+2)) - fmp/(2.0_DP*rm**(l+1))
       !print '(i5,6f12.7)',isplit,splitnorm(isplit-1),rm,fm,fmp,al,bl

       ! Check solution has worked
       fmc = rm**l*(al-bl*rm**2)
       fmpc = rm**(l-1)*(l*al-(l+2)*bl*rm**2)
       if ((abs(fm-fmc)>tol).or.(abs(fmp-fmpc)>tol)) then
          call utils_abort('Error in atom_split_orbital: Attempt to split &
               &function failed when matching values at rm')
       end if

       ! Evaluate function from 1 to iptm-1
       do ipt=1,iptm-1
          func(ipt,isplit) = basis%rad(ipt)**l*(al-bl*basis%rad(ipt)**2)
       end do

       ! Copy original function from rm to rmax into new function
       func(iptm:basis%npts,isplit) = func(iptm:basis%npts,1)

       ! cks: Subtract new function from original one and put into new fucntion
       func(1:basis%npts,isplit) = &
            func(1:basis%npts,1) - func(1:basis%npts,isplit)


       ! Renormalise new function
       norm = services_regular_integral(basis%npts,basis%dr, &
            func(1:basis%npts,isplit)**2*basis%rad(1:basis%npts)**2)
       func(1:basis%npts,isplit) = func(1:basis%npts,isplit) / sqrt(norm)

    end do

  end subroutine atom_split_orbital


  subroutine atom_polarise_orbital(func,func_orig,eig_orig,ang_mom,atom,basis)

    !==========================================================================!
    ! This subroutine uses perturbation theory to determine the radial part    !
    ! of the polarisation orbitals of a given orbital. These have angular      !
    ! momentum ang_mom+1.                                                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   func (output) : New radial functions of angular momentum ang_mom+1     !
    !   func_orig (input) : Orig radial functions on grid defined by basis%rad !
    !   eig_orig (input) : Hamiltonian eigenvalue of original function         !
    !   ang_mom (input) : Angular momentum of original function                !
    !   atom (inout) : ATOM_TYPE object describing the atom                    !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in January 2012.                                !
    !==========================================================================!

    use services, only: services_regular_integral
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check
    use linalg, only: linalg_invert_sym_matrix

    implicit none

    ! Arguments
    type(BASIS_TYPE),intent(in) :: basis
    type(ATOM_TYPE),intent(inout) :: atom
    real(kind=DP), intent(out) :: func(1:basis%npts)
    real(kind=DP), intent(in) :: func_orig(1:basis%npts)
    real(kind=DP), intent(in) :: eig_orig
    integer, intent(in) :: ang_mom

    ! Local Variables
    integer :: ibfn,jbfn,nbfn,il
    integer :: ierr
    real(kind=DP) :: norm
    real(kind=DP), allocatable :: ham(:,:), ham_inv(:,:), overlap(:,:)

    ! BLAS subroutine
    external :: dgemv

    ! Allocate matrix workspaces
    allocate(ham(basis%nbfn_max,basis%nbfn_max),stat=ierr)
    call utils_alloc_check('atom_polarise_orbital','ham',ierr)
    allocate(ham_inv(basis%nbfn_max,basis%nbfn_max),stat=ierr)
    call utils_alloc_check('atom_polarise_orbital','ham_inv',ierr)
    allocate(overlap(basis%nbfn_max,basis%nbfn_max),stat=ierr)
    call utils_alloc_check('atom_polarise_orbital','overlap',ierr)

    ! Get Hamiltonian for L+1
    il = ang_mom + 1
    nbfn = basis%nbfn(il)
    call atom_ham_matrix(ham,overlap,il,atom,basis)

    ! Fill in lower triangle of hamiltonian and overlap
    do ibfn=1,nbfn
       do jbfn=min(ibfn+1,nbfn),nbfn
          ham(jbfn,ibfn) = ham(ibfn,jbfn)
          overlap(jbfn,ibfn) = overlap(ibfn,jbfn)
       end do
    end do

    ! Subtract eigenvalue times overlap matrix
    ham(:,:) = ham(:,:) - overlap(:,:) * eig_orig

    ! Invert Hamiltonian Matrix
    ham_inv = ham
    call linalg_invert_sym_matrix(ham_inv,nbfn)

    ! Calculate overlap of basis functions with original function times r
    ! Use overlap(:,1) as workspace
    do ibfn=1,nbfn
       atom%work(:) = basis%rad(:) * func_orig(:) * &
            basis%bfn(:,ibfn,il) * basis%rad(:)**2
       overlap(ibfn,1) = services_regular_integral(basis%npts,basis%dr, &
            atom%work(1:basis%npts))
    end do

    ! Apply inverse Hamiltonian to matrix of overlaps between basis functions
    ! and original function times r to get coeffs of polarisation function
    ! Use overlap(:,2) to hold resulting coefficients of bfns
    call dgemv('N',nbfn,nbfn,1.0_DP,ham_inv(1,1), &
         basis%nbfn_max,overlap(1,1),1,0.0_DP,overlap(1,2),1)

    ! Get new function on real space grid from bfn coefficients
    call dgemv('N',basis%npts,nbfn,1.0_DP,basis%bfn(1,1,il), &
         basis%npts,overlap(1,2),1,0.0_DP,func(1),1)

    ! Renormalise new function
    norm = services_regular_integral(basis%npts,basis%dr, &
         func(1:basis%npts)**2*basis%rad(1:basis%npts)**2)
    func(1:basis%npts) = func(1:basis%npts) / sqrt(norm)

    ! Deallocate matrix workspaces
    deallocate(overlap,stat=ierr)
    call utils_dealloc_check('atom_polarise_orbital','overlap',ierr)
    deallocate(ham_inv,stat=ierr)
    call utils_dealloc_check('atom_polarise_orbital','ham_inv',ierr)
    deallocate(ham,stat=ierr)
    call utils_dealloc_check('atom_polarise_orbital','ham',ierr)

  end subroutine atom_polarise_orbital

  subroutine atom_write_orbitals(atom,basis,filename)

    !==========================================================================!
    ! This subroutine writes the single-electron orbitals of an atom to a file.!
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (input) : ATOM_TYPE object describing the atom.                   !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !   filename (input) : file to write orbitals to.                          !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use utils, only: utils_open_unit_check, utils_close_unit_check, utils_unit

    implicit none

    ! Arguments
    type(ATOM_TYPE), intent(in) :: atom
    type(BASIS_TYPE),intent(in) :: basis
    character(len=80),intent(in) :: filename

    ! Local Variables
    integer :: ierr
    integer :: ipt, iorb
    integer :: iunit

    ! ndmh: open the file
    iunit = utils_unit()
    open(iunit,file=trim(filename),iostat=ierr)
    call utils_open_unit_check('atom_write_orbitals',trim(filename),ierr)

    ! ndmh: write the orbitals
    do ipt=1,basis%npts
       write(iunit,'(f22.15)',advance='no') basis%rad(ipt)
       do iorb=1,atom%norbs
          write(iunit,'(f22.15)',advance='no') atom%psi_r(ipt,iorb)
       end do
       write(iunit,*)
    end do

    ! ndmh: close the file
    close(iunit,iostat=ierr)
    call utils_close_unit_check('atom_write_orbitals',trim(filename),ierr)

  end subroutine atom_write_orbitals


  subroutine atom_write_fireball(atom,basis,filename)

    !==========================================================================!
    ! This subroutine writes a fireball file using the orbitals of an atom.    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (input) : ATOM_TYPE object describing the atom.                   !
    !   basis (input) : BASIS_TYPE object describing the basis.                !
    !   filename (input) : file to write orbitals to.                          !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                              !
    !==========================================================================!

    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_open_unit_check, utils_close_unit_check, utils_unit

    implicit none

    ! Arguments
    type(ATOM_TYPE), intent(in) :: atom
    type(BASIS_TYPE),intent(in) :: basis
    character(len=80),intent(in) :: filename

    ! Local Variables
    integer :: ierr
    integer :: ipt, iorb
    integer :: iunit
    character(len=1) :: tmp1,tmp2

    ! ndmh: open the file
    iunit = utils_unit()
    open(iunit,file=trim(filename),iostat=ierr)
    call utils_open_unit_check('atom_write_fireball',trim(filename),ierr)

    ! ndmh: write the angular momenta
    write(iunit,'(a)') 'ANGULAR_MOMENTA'
    do iorb=1,atom%norbs
       write(iunit,'(f22.15)',advance='no') atom%orb_ang_mom(iorb)
    end do
    write(iunit,'(a)') 'END_ANGULAR_MOMENTA'

    ! ndmh: write the radial grid positions
    write(iunit,'(a)') 'RADIAL_POSITIONS'
    do ipt=1,basis%npts
       write(iunit,'(f22.15)',advance='no') basis%rad(ipt)
    end do
    write(iunit,'(a)') 'END_RADIAL_POSITIONS'

    ! ndmh: write the orbitals
    do iorb=1,atom%norbs
       write(tmp1,'(i5)') iorb
       write(tmp2,'(i5)') atom%orb_ang_mom(iorb)
       write(iunit,'(a,a,a,a)') 'SHELL_',trim(adjustl(tmp1)), &
            '_ANGMOM_',trim(adjustl(tmp2))
       do ipt=1,basis%npts
          write(iunit,'(f22.15)',advance='no') atom%psi_r(ipt,iorb)
       end do
       write(iunit,'(a,a,a,a)') 'END_SHELL_',trim(adjustl(tmp1)), &
            '_ANGMOM_',trim(adjustl(tmp2))
    end do

    ! ndmh: close the file
    close(iunit,iostat=ierr)
    call utils_close_unit_check('atom_write_fireball',trim(filename),ierr)

  end subroutine atom_write_fireball

  subroutine atom_set_density_guess(atom,basis,charge,spin,report,num_report_lines)

    !==========================================================================!
    ! This subroutine guesses the radial electron density of an atom in a      !
    ! given charge and spin state.                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atom (in/out)  : ATOM_TYPE object describing the atom.                 !
    !   basis (input)  : BASIS_TYPE object describing the basis.               !
    !   charge (input) : Guess charge (currently unused).                      !
    !   spin (input)   : Guess spin (currently unused).                        !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in May 2012.                                    !
    !==========================================================================!

    use rundat, only: pub_num_spins
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(ATOM_TYPE), intent(inout) :: atom
    type(BASIS_TYPE),intent(in) :: basis
    real(kind=DP), intent(in) :: charge
    real(kind=DP), intent(in) :: spin
    character(len=128),intent(inout) :: report(400)
    integer, intent(inout) :: num_report_lines

    ! Local Variables
    integer :: ierr
    integer :: ir, iorb, is
    integer :: max_occ_iorb, min_unocc_iorb
    integer :: norbs_write
    integer :: counter
    real(kind=DP) :: maxocc
    real(kind=DP) :: changeleft
    real(kind=DP) :: change
    real(kind=DP) :: spin_fac
    real(kind=DP) :: max_eig, min_eig
    logical :: dbg=.false.
    character(len=80) :: fmt,tmp

    ! ndmh: allocate storage
    allocate(atom%guess_den(basis%npts,pub_num_spins),stat=ierr)
    call utils_alloc_check('atom_set_density_guess','atom%guess_den',ierr)
    allocate(atom%guess_comp_den(basis%npts,pub_num_spins),stat=ierr)
    call utils_alloc_check('atom_set_density_guess','atom%guess_comp_den',ierr)
    allocate(atom%guess_occ(atom%norbs,pub_num_spins),stat=ierr)
    call utils_alloc_check('atom_set_density_guess','atom%guess_occ',ierr)
    allocate(atom%guess_rhoij(atom%nshells,atom%nshells,pub_num_spins),stat=ierr)
    call utils_alloc_check('atom_set_density_guess','atom%guess_rhoij',ierr)

    if ((spin>0.0_DP).and.(pub_num_spins==1)) then
       call utils_abort(' ERROR in atom_set_density_guess: Nonzero initial &
            &spin requested, but calculation is not spin-polarised')
    end if

    spin_fac = 2.0_DP / real(pub_num_spins,DP)

    do is=1,pub_num_spins

       atom%guess_occ(:,is) = atom%occ(:) / real(pub_num_spins,DP)

       ! ndmh: adjust configuration occupancies according to guess charge / spin
       if (abs(charge)>0.0_DP) then
          if (is==1) changeleft = -charge / real(pub_num_spins,DP)
          if (is==2) changeleft = -charge / real(pub_num_spins,DP)

          if (dbg) print *
          if (dbg) print *,'CHARGE:'
          if (dbg) print *,'is,changeleft',is,changeleft
          if (dbg) print *,'occ:     ',atom%guess_occ(:,is)
          if (dbg) print *,'eigs:    ',atom%eigs(:)

          counter = 0
          do while (abs(changeleft)>tiny(1.0_DP))
             if (counter>9) dbg = .true.
             if (counter>10) then
                call utils_abort(' ERROR in atom_set_density_guess: Failed to &
                     &adjust total charge on species ',atom%pspecies_number)
             end if
             counter = counter + 1

             ! ndmh: find highest occupied orbital
             max_occ_iorb = 1
             do iorb=1,atom%norbs
                if ((atom%guess_occ(iorb,is)>tiny(1.0_DP)) .and. &
                     (atom%eigs(iorb)>atom%eigs(max_occ_iorb))) max_occ_iorb = iorb
             end do

             ! ndmh: find lowest unoccupied orbital
             min_unocc_iorb = atom%norbs
             do iorb=atom%norbs,1,-1
                if ((abs(atom%guess_occ(iorb,is))<tiny(1.0_DP)) .and. &
                    (atom%eigs(iorb)<atom%eigs(min_unocc_iorb))) min_unocc_iorb = iorb
             end do

             if (dbg) print *
             if (dbg) print *,'  OCC: ', max_occ_iorb,atom%eigs(max_occ_iorb), atom%guess_occ(max_occ_iorb,is)
             if (dbg) print *,'UNOCC: ',min_unocc_iorb,atom%eigs(min_unocc_iorb), atom%guess_occ(min_unocc_iorb,is)

             if (changeleft>0.0_DP) then
                maxocc = real(atom%orb_ang_mom(max_occ_iorb)*2+1,DP)*spin_fac
                if (abs(atom%guess_occ(max_occ_iorb,is)-maxocc)<tiny(1.0_DP)) then
                   ! highest occupied orbital is already full
                   iorb = min_unocc_iorb
                else
                   ! highest occupied orbital still has some space
                   iorb = max_occ_iorb
                end if
                change = changeleft
                maxocc = real(atom%orb_ang_mom(iorb)*2+1,DP)*spin_fac
                if ((atom%guess_occ(iorb,is)+change)>maxocc) then
                   change = maxocc - atom%guess_occ(max_occ_iorb,is)
                end if
             else
                iorb = max_occ_iorb
                change = changeleft
                if ((atom%guess_occ(iorb,is)+change)<tiny(1.0_DP)) then
                   change = -atom%guess_occ(max_occ_iorb,is)
                end if
             endif
             atom%guess_occ(iorb,is) = atom%guess_occ(iorb,is) + change
             changeleft = changeleft - change
             if (dbg) print *,'change = ',change
             if (dbg) print *,'changeleft = ',changeleft
          end do
          if (dbg) print *,'occ:     ',atom%guess_occ(:,is)
       end if

       ! ndmh: adjust configuration occupancies according to guess charge / spin
       if ((abs(spin)>0.0_DP)) then
          if (is==1) changeleft = spin*0.5_DP
          if (is==2) changeleft = -spin*0.5_DP
          if (dbg) print *
          if (dbg) print *,'SPIN:'
          if (dbg) print *,'is,changeleft',is,changeleft
          if (dbg) print *,'occ:     ',atom%guess_occ(:,is)
          if (dbg) print *,'eigs:    ',atom%eigs(:)

          counter = 0
          do while (abs(changeleft)>tiny(1.0_DP))
             if (counter>9) dbg = .true.
             if (counter>10) then
                call utils_abort(' ERROR in atom_set_density_guess: Failed to &
                     &adjust total spin on species ',atom%pspecies_number)
             end if
             counter = counter + 1

             ! ndmh: find highest occupied orbital
             ! ebl & kkbd: Find lowest non-fully occupied orbital
             max_occ_iorb = 1
             min_eig = MAXVAL(atom%eigs(:)) + 1.0_dp
             do iorb=1,atom%norbs
                if (atom%guess_occ(iorb,is)<real(atom%orb_ang_mom(iorb)*2+1,DP) .and. &
                     (atom%eigs(iorb) <= min_eig)) then
                      max_occ_iorb = iorb
                      min_eig = atom%eigs(iorb)
                end if
             end do

             ! ndmh: find lowest unoccupied orbital
             ! ebl & kkbd: Find highest partially occupied orbital
             min_unocc_iorb = atom%norbs
             max_eig = MINVAL(atom%eigs(:)) - 1.0_dp
             do iorb=atom%norbs,1,-1
                if ((abs(atom%guess_occ(iorb,is))>0.0_DP) .and. &
                    (atom%eigs(iorb)>max_eig)) then
                       min_unocc_iorb = iorb
                       max_eig = atom%eigs(iorb)
                end if
             end do

             if (dbg) print *,max_occ_iorb,atom%eigs(max_occ_iorb)
             if (dbg) print *,min_unocc_iorb,atom%eigs(min_unocc_iorb)
             if (dbg) print *,'occ:     ',atom%guess_occ(:,is)

             if (changeleft>0.0_DP) then
                iorb = max_occ_iorb
                change = changeleft
                maxocc = real(atom%orb_ang_mom(iorb)*2+1,DP)
                if ((atom%guess_occ(iorb,is)+change)>maxocc) then
                   change = maxocc - atom%guess_occ(iorb,is)
                end if
             else
                iorb = min_unocc_iorb
                change = changeleft
                if ((atom%guess_occ(iorb,is)+change)<tiny(1.0_DP)) then
                   change = -atom%guess_occ(iorb,is)
                end if
             endif
             atom%guess_occ(iorb,is) = atom%guess_occ(iorb,is) + change
             changeleft = changeleft - change
             if (dbg) print *,'change = ',change
             if (dbg) print *,'changeleft = ',changeleft
          end do
          if (dbg) print *,'occ:     ',atom%guess_occ(:,is)
       end if

    end do

    if ((abs(charge)>0.0_DP).or.(abs(spin)>0.0_DP)) then
       write(report(num_report_lines+1),'(a)') 'Adjusting orbital occupancies:'
       num_report_lines = num_report_lines + 1
       norbs_write = min(atom%norbs,15)
       write(tmp,'(i10)') norbs_write
       do is=1,pub_num_spins
          write(fmt,'(5a)') '(a,2i3,4x,',trim(adjustl(tmp)),'f5.2)'
          write(report(num_report_lines+1),fmt) 'Orbitals (num,spin,occ):', &
               atom%norbs,is,atom%guess_occ(1:norbs_write,is)
          num_report_lines = num_report_lines + 1
          write(fmt,'(5a)') '(a,2i3,4x,',trim(adjustl(tmp)),'i5)'
          write(report(num_report_lines+1),fmt) 'Orbitals   (num,spin,l):', &
               atom%norbs,is,atom%orb_ang_mom(1:norbs_write)
          num_report_lines = num_report_lines + 1
       end do
    end if

    if ((sum(atom%guess_occ(:,:))+charge-atom%ion_charge)>tiny(1.0_DP)) then
       call utils_abort(' ERROR: Atomsolver failed to initialise charge')
    end if

    if (abs(spin)>0.0_DP) then
       if ((sum(atom%guess_occ(:,1))-sum(atom%guess_occ(:,2))-spin)>tiny(1.0_DP)) then
       call utils_abort(' ERROR: Atomsolver failed to initialise spin')
       end if
    end if

    do is=1,pub_num_spins

       ! ndmh: recalculate density
       atom%work(:) = 0.0_DP
       do iorb=1,atom%norbs
          do ir=1,basis%npts
             atom%work(ir) = atom%work(ir) &
                  + atom%psi_r(ir,iorb)**2*atom%guess_occ(iorb,is)
          end do
       end do
       atom%den(1:basis%npts) = atom%work(1:basis%npts)

       ! ndmh: recalculate rhoij and compensation density
       atom%occ = atom%guess_occ(:,is)
       call atom_rhoij(atom,basis)
       atom%guess_rhoij(:,:,is) = atom%rhoij(:,:)
       call atom_compden(atom,basis)

       ! ndmh: store in (spin-polarised) guess density arrays
       atom%guess_den(1:basis%npts,is) = atom%den(1:basis%npts)
       atom%guess_comp_den(1:basis%npts,is) = atom%comp_den(1:basis%npts)

    end do

  end subroutine atom_set_density_guess

  integer function atom_gs_occ(atomic_number,nl_index,flag)

    !==========================================================================!
    ! This function simply accesses the atomic configuration array gs_occ(:,:) !
    ! Mostly written for use in ddec_mod                                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   atomic_number (input) : Integer describing the atomic number.          !
    !   nl_index (input) : Integer for the nl_index (2nd index) of gs_occ.     !
    !   flag (input) : 'Size' to return size of the 2nd index of gs_occ for a  !
    !                  given atomic_number (1st index). 'Occupancy' to return  !
    !                  the non-negative (-1 gets set to 0) value of the array  !
    !                  index gs_occ(atomic_number,nl_index). Leave blank to    !
    !                  simply access gs_occ(atomic_number,nl_index).           !
    !--------------------------------------------------------------------------!
    ! Written by Louis Lee in September 2013.                                  !
    !==========================================================================!

     use utils, only: utils_abort

     implicit none
     integer, intent(in) :: atomic_number, nl_index
     character(len=*), optional, intent(in) :: flag

     integer :: occ
     logical :: return_occ

     return_occ = .false.

     if ( present(flag) ) then
        select case (flag)
           case ('Size')
              atom_gs_occ = size(gs_occ,atomic_number)
              return
           case ('Occupancy')
              return_occ = .true.
           case default
              call utils_abort(' ERROR: Unknown flag in atom_gs_occ')
        end select
     end if

     if ( atomic_number > ubound(gs_occ,1) .or. &
          atomic_number < lbound(gs_occ,1) ) call utils_abort(' ERROR: &
          &atomic_number out of range in atom_gs_occ')

     if ( nl_index > ubound(gs_occ,2) .or. &
          nl_index < lbound(gs_occ,2) ) call utils_abort(' ERROR: &
          &nl_index out of range in atom_gs_occ')

     ! If 'return_occ' is true, set negative occupancies to zero
     occ = gs_occ(atomic_number,nl_index)
     if (return_occ) then
        if ( atomic_number == 0 ) occ = 0 ! no ion
        if ( occ < 0 ) occ = 0
     end if
     atom_gs_occ = occ

  end function atom_gs_occ

end module atom
