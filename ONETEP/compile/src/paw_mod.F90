! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!             Projector-Augmented Wave module                    !
!                                                                !
! This module contains routines associated with the use of the   !
! Projector-Augmented Wave method, including initialisation,     !
! calculation of the augmentation density, calculation of the    !
! nonlocal energy terms, and calculation of the forces.          !
!----------------------------------------------------------------!
! This module was created by Nicholas Hine in May 2010.          !
!================================================================!

module paw

  use constants, only: DP, PI, SQRT_PI, NORMAL, VERBOSE, stdout
  use paw_shape, only: PAW_SHAPE_INFO
  use rundat, only: pub_debug_on_root

  implicit none

  private

  ! Type describing radial grids
  type RADIAL_GRID

     ! ndmh: number of points in this radial grid
     integer :: npt

     ! ndmh: step sizes for radial grid
     real(kind=DP) :: rad_step, log_step

     ! ndmh: positions on radial grid
     real(kind=DP), allocatable :: r(:)

     ! ndmh: "spacings" of radial grid
     real(kind=DP), allocatable :: rab(:)

  end type RADIAL_GRID

  ! Type describing PAW species
  type PAW_SPECIES

     ! ndmh: dataset name
     character(len=64) :: ds_name

     ! ndmh: core wavefunctions filename
     character(len=64) :: core_wf_name

     ! ndmh: charge of core of atom
     real(kind=DP) :: ion_charge

     ! ndmh: atomic number of atom
     integer :: atomic_number

     ! ndmh: PAW cutoff radius
     real(kind=DP) :: rcut

     ! ndmh: PAW projector max radius
     real(kind=DP) :: proj_rcut

     ! ndmh: Number of projectors / partial waves (counting n, l only)
     integer :: npw

     ! ndmh: Number of projectors / partial waves (counting n, l, m)
     integer :: npw_tot

     ! ndmh: Angular momentum of each partial wave / projector pair
     integer, allocatable :: l_pw(:)
     integer, allocatable :: l_pw_tot(:),m_pw_tot(:)
     integer, allocatable :: ipw_tot(:)

     ! ndmh: Maximum angular momentum of any partial wave / projector pair
     integer :: lmax

     ! ndmh: Number of meshes defined
     integer :: ngrid

     ! ndmh: Which mesh the partial waves use
     integer :: phi_grid

     ! ndmh: Which mesh the projectors use
     integer :: proj_grid

     ! ndmh: Which mesh the core densities use
     integer :: core_den_grid

     ! ndmh: Which mesh the local potential uses
     integer :: vhntzc_grid

     ! ndmh: The grid to use for the shape function
     integer :: shape_grid

     ! ndmh: The format of the vhntzc potential
     integer :: vhntzc_format

     ! ndmh: Array to hold information about grids
     type(RADIAL_GRID), allocatable :: grid(:)

     ! ndmh: Type to hold information about shape function
     type(PAW_SHAPE_INFO) :: shape

     integer :: n_recip_pts
     real(kind=DP) :: g_max
     real(kind=DP) :: inv_g_spacing

     ! ndmh: AE and PS partial waves on radial grid
     real(kind=DP), allocatable :: phi_rad(:,:)
     real(kind=DP), allocatable :: tphi_rad(:,:)

     ! ndmh: projectors on radial grid
     real(kind=DP), allocatable :: tproj_rad(:,:)
     real(kind=DP), allocatable :: tproj_recip(:,:)

     ! ndmh: core densities on radial grids (real/recip)
     real(kind=DP), allocatable :: core_den_rad(:)
     real(kind=DP), allocatable :: tcore_den_rad(:)
     real(kind=DP), allocatable :: tcore_den_recip(:)

     ! ndmh: Frozen part of nonlocal energies Dij
     real(kind=DP), allocatable :: dij0(:,:)

     ! ndmh: Initial guess for projector density kernel
     real(kind=DP), allocatable :: rhoij0(:,:)

     ! ndmh: Hartree potential of core, real and recip radial grids
     real(kind=DP), allocatable :: vhntzc_rad(:)
     real(kind=DP), allocatable :: vhntzc_recip(:)

     ! ndmh: for augmentation charge:
     ! ndmh: \int (phi_i(r)\phi_j(r)-\tphi_i(r)\tphi_j(r))r^L dr
     real(kind=DP), allocatable :: aug_nLij(:,:,:)

     ! ndmh: for onsite Hartree energy
     real(kind=DP), allocatable :: e_ijkl(:,:,:,:)

     ! ndmh: whether AE core density is nonzero
     logical :: core_charge
     logical :: core_charge_calculated

     ! ndmh: whether PS core density is nonzero
     logical :: tcore_charge

     ! ndmh: XC energy of AE core density
     real(kind=DP) :: exc_core

     ! ---- FOR CORE LEVEL RECONSTRUCTION ----
     logical :: core_wvfns_exist

     ! ndmh: number core orbitals on radial grid (and number including m)
     integer :: n_core_wfs, n_core_wfs_tot

     ! ndmh: angular momentum values of each core wf
     integer, allocatable :: l_core_wf(:)
     integer, allocatable :: l_core_wf_tot(:)
     integer, allocatable :: m_core_wf_tot(:)
     integer, allocatable :: icore_wf_tot(:)

     ! ndmh: principle quantum number of each state (for labels)
     integer, allocatable :: n_core_wf(:)

     ! ndmh: Array to hold information about core wf grids
     integer :: ngrid_core
     integer :: core_wf_grid
     type(RADIAL_GRID), allocatable :: grid_core(:)
     real(kind=DP) :: rcut_core

     ! ndmh: charge of core of atom
     real(kind=DP) :: core_wf_charge

     ! ndmh: Core orbitals on radial grid
     real(kind=DP), allocatable :: core_wf_rad(:,:)
     real(kind=DP), allocatable :: core_wf_recip(:,:)
     real(kind=DP), allocatable :: core_wf_eig(:)
     real(kind=DP), allocatable :: core_wf_occ(:)

     ! lr408: \int (psi_c(r)\phi_j(r)-\psi_c(r)\tphi_j(r))r^L dr
     real(kind=DP), allocatable :: core_aug_nLij(:,:,:)

     ! ndmh: --- Internal state of the local openbc pseudo ---
     logical :: openbc_locps_initialised = .false.
     ! ndmh: Vloc(x) on radial grid, for all species
     real(kind=DP), allocatable :: vloc_lookup(:)
     ! ndmh: d/dx(Vloc(x)) on radial grid, for all species,
     !     only filled if pub_forces_needed
     real(kind=DP), allocatable :: vlocder_lookup(:)
     ! ndmh: -------------------------------------------------

  end type PAW_SPECIES

  public :: PAW_SPECIES

  ! Maximum number of partial waves on any atom (including m_i)
  integer :: max_paw_proj_tot
  integer :: max_core_wf_tot

  real(kind=DP), parameter :: sqrt_4pi = 2.0_DP*SQRT_PI
  real(kind=DP), parameter :: inv_sqrt_4pi = 1.0_DP/sqrt_4pi

  ! Public subroutines
  public :: paw_read_species
  public :: paw_all_species_exit
  public :: paw_tcore_hartree_on_grid
  public :: paw_tcore_density
  public :: paw_tcore_density_recip
  public :: paw_species_init_proj
  public :: paw_species_calc_proj_prec_mat
  public :: paw_projector_overlap
  public :: paw_position_operator
  public :: paw_grad_operator
  public :: paw_projector_denskern_init
  public :: paw_nonlocal_energies
  public :: paw_show_atomblocks
  public :: paw_tcore_hartree_calc_forces
  public :: paw_nlcc_calculate_forces
  public :: paw_nlcc_calculate_forces_recip
  public :: paw_sphere_density_on_grid
  public :: paw_atom_aug_den
  public :: paw_atom_aug_integrals
  public :: paw_atom_aug_force

  public :: paw_get_projector_info
  public :: paw_get_locpot_rad
  public :: paw_get_core_den_rad
  public :: paw_get_projectors_q
  public :: paw_get_aug_funcs
  public :: paw_exc_core_init
  public :: paw_exc_core_atom
  public :: paw_dij_hartree_atom
  public :: paw_dij_hartree
  public :: paw_dij_xc_atom
  public :: paw_dij_xc_tddft
  public :: paw_dij_so
  public :: paw_store_initial_proj_kern

  public :: paw_species_init_core_wvfns
  public :: paw_species_exit_core_wvfns
  public :: paw_core_pw_position_operator
  public :: paw_core_pw_overlap

  ! gcc32: for LR_PHONONS
  public :: paw_atom_grad_aug_integrals
  public :: paw_FO_nlcc_energy
contains

  subroutine paw_read_species(paw_sp,elements,species,par)

    !==================================================================!
    ! This subroutine reads the elements array and works out what PAW  !
    ! species are present, allocates storage for them in the paw_sp    !
    ! array and then read the PAW dataset file for each species into   !
    ! the paw_sp array                                                 !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  elements (inout) : list of atoms present                        !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 17/05/2010.                          !
    !==================================================================!

    use bibliography, only: bibliography_cite
    use comms, only: pub_on_root, pub_imroots_comm, comms_barrier
    use image_comms, only: pub_my_image
    use gaunt_coeff, only: gaunt_init
    use ion, only: ELEMENT, SPECIE
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_nlcc, pub_aug, pub_usp, pub_nhat_in_xc, &
                      pub_num_images
    use utils, only: utils_abort, utils_alloc_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout), pointer     :: paw_sp(:)
    type(PARAL_INFO), intent(inout)               :: par
    type(ELEMENT), intent(inout)                  :: elements(par%nat)
    type(SPECIE), intent(in)                      :: species(par%num_species)
    ! Local Variables
    logical :: found
    integer :: ierr, i, img
    integer :: iat, jat
    integer :: isp,nsp,msp
    integer :: lmax

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_read_species'

    call bibliography_cite('PAW')

    ! count unique species
    nsp = 0
    do iat=1,par%nat
       found = .false.
       do jat=1,iat-1
          if ((species(elements(iat)%species_number)%pseudo_name==&
               species(elements(jat)%species_number)%pseudo_name) .and. &
               (species(elements(iat)%species_number)%core_wf_name==&
               species(elements(jat)%species_number)%core_wf_name)) then
             found = .true.
             exit
          end if
       end do
       if (.not.found) then
          nsp = nsp + 1
       end if
    end do

    allocate(paw_sp(nsp),stat=ierr)
    call utils_alloc_check('paw_read_species','paw_sp',ierr)

    ! set dataset file names
    paw_sp(:)%ds_name = ''
    msp = 0
    do iat=1,par%nat
       found = .false.
       do isp=1,nsp
          if ((species(elements(iat)%species_number)%pseudo_name==&
               paw_sp(isp)%ds_name) .and. &
               (species(elements(iat)%species_number)%core_wf_name==&
               paw_sp(isp)%core_wf_name)) then
             found = .true.
             exit
          end if
       end do
       if (.not.found) then
          msp = msp + 1
          paw_sp(msp)%ds_name = &
             species(elements(iat)%species_number)%pseudo_name
          paw_sp(msp)%core_wf_name = &
             species(elements(iat)%species_number)%core_wf_name
       end if
    end do

    ! Consistency check
    if (msp/=nsp) then
       call utils_abort('Error: Inconsistent dataset names in paw_read_species')
    end if

    if (pub_on_root) then
       ! kkbd: only one image should do file manipulation at once
       do img=0, pub_num_images-1
          if ((pub_num_images == 1) .or.(img == pub_my_image)) then
             ! Read Datasets
             do isp=1,nsp
                call paw_read_dataset(paw_sp(isp))
             end do

             ! Read Core wavefunctions
             paw_sp(:)%core_wvfns_exist = .false.
             do isp=1,nsp
                if (adjustl(trim(paw_sp(isp)%core_wf_name))/='NONE') then

                   call paw_read_core_wfs(paw_sp(isp))

                   if (paw_sp(isp)%n_core_wfs > 0) then
                      paw_sp(isp)%core_wvfns_exist = .true.
                   else
                      paw_sp(isp)%n_core_wfs = 0
                      paw_sp(isp)%n_core_wfs_tot = 0
                      paw_sp(isp)%ngrid_core = 0
                      paw_sp(isp)%core_wf_grid = 0
                      paw_sp(isp)%rcut_core = 0.0_DP
                   end if
                else
                   paw_sp(isp)%n_core_wfs = 0
                   paw_sp(isp)%n_core_wfs_tot = 0
                   paw_sp(isp)%ngrid_core = 0
                   paw_sp(isp)%core_wf_grid = 0
                   paw_sp(isp)%rcut_core = 0.0_DP
                end if
             end do
          end if
          if (pub_num_images > 1) call comms_barrier(pub_imroots_comm)
       end do ! image mutex
    end if

    ! Broadcast data read from dataset from root proc to all other procs
    do isp=1,nsp
       call paw_bcast_dataset(paw_sp(isp))
    end do

    ! Copy relevant information back into the species array on all procs
    elements(:)%pspecies_number = -1
    do iat=1,par%nat

       ! Loop over species found and check if name matches
       do isp=1,nsp
          if ((paw_sp(isp)%ds_name == &
               species(elements(iat)%species_number)%pseudo_name) .and. &
               (paw_sp(isp)%core_wf_name == &
               species(elements(iat)%species_number)%core_wf_name)) then
             elements(iat)%pspecies_number = isp
             elements(iat)%atomic_number = paw_sp(isp)%atomic_number
             elements(iat)%ion_charge = paw_sp(isp)%ion_charge
             elements(iat)%npawpws = paw_sp(isp)%npw_tot
             elements(iat)%nprojectors = 0
             elements(iat)%max_core_radius = paw_sp(isp)%proj_rcut
             elements(iat)%ncorewfs = paw_sp(isp)%n_core_wfs_tot
             elements(iat)%max_core_wf_radius = paw_sp(isp)%rcut_core
             exit
          end if
       end do

       ! Check we have found a species number for every species
       if (elements(iat)%pspecies_number == -1) then
          call utils_abort('Error in paw_read_species: No species found to &
               &match'//species(elements(iat)%species_number)%pseudo_name)
       end if

    end do

    ! Count total number of PAW partial waves in system
    par%num_projectors = 0
    par%num_pawpws = 0
    par%num_corewfs = 0
    do iat=1,par%nat
       par%num_pawpws = par%num_pawpws + &
            paw_sp(elements(iat)%pspecies_number)%npw_tot

       par%num_corewfs = par%num_corewfs + &
            paw_sp(elements(iat)%pspecies_number)%n_core_wfs_tot
    end do

    ! Display species information in table
    if (pub_on_root) call internal_print_datasets

    ! Count highest angular momentum over all species
    lmax = 0
    do isp=1,par%num_pspecies
       lmax=max(paw_sp(isp)%lmax,lmax)
    end do

    ! Initialise Gaunt Coefficients
    call gaunt_init(lmax+1)

    ! Set module-level variables
    max_core_wf_tot = maxval(paw_sp(:)%n_core_wfs_tot) ! lr408
    max_paw_proj_tot = maxval(paw_sp(:)%npw_tot)
    if (all(paw_sp(:)%vhntzc_format==2)) then
       pub_nhat_in_xc = .false.
    else if (all(paw_sp(:)%vhntzc_format==1)) then
       pub_nhat_in_xc = .true.
    else
       call utils_abort('Error in paw_read_species: vhntzc_format must be &
            &same for all datasets')
    end if

    ! No NLCC unless we find a species with a nonzero PS core charge
    pub_nlcc = .false.

    ! Calculate all the pre-calculated information required to perform
    ! a PAW calculation, for each species
    do isp=1,nsp
       call paw_dataset_init(paw_sp(isp),paw_sp(isp)%core_wvfns_exist)
    end do

    ! Charge augmentation will always be active in PAW calculations
    pub_aug = .true.

    ! USP and PAW cannot be simultaneously active
    pub_usp = .false.

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_read_species'

contains

    subroutine internal_print_datasets

      !=========================================================-!
      ! This subroutine prints out the number of atoms, NGWFs    !
      ! and partial waves for each species.                      !
      !----------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 25/03/2007.          !
      ! Adapted for PAW by Nicholas Hine on 20/05/2010.          !
      !==========================================================!

      implicit none

      ! ndmh: Local Variables
      integer :: ipw

      character(len=64)    :: current_file


      write(stdout,'(a)') '<<<<<<<<<<<<<<<<<<<<<<<<<<< &
           &PAW Dataset information >>>>>>>>>>>>>>>>>>>>>>>>>>>'

      do isp=1,par%num_pspecies

         ! ndmh: get file name for this species
         current_file = paw_sp(isp)%ds_name

         ! ndmh: print basic information about this pseudopotential
         write(stdout,'(3a,f5.2,a,i3,a,f5.2,a)') 'File: ',trim(current_file), &
                 '; rc =',paw_sp(isp)%rcut,' bohr; shape type =', &
                 paw_sp(isp)%shape%shape_type, '; rshape =', &
                 paw_sp(isp)%shape%rshape,' bohr;'
         write(stdout,'(a,i3,a,f10.6)') '  Atomic number:', &
              paw_sp(isp)%atomic_number,';  ionic charge:', &
              paw_sp(isp)%ion_charge
         do ipw=1,paw_sp(isp)%npw
            write(stdout,'(2(a,i2))') '    Partial Wave',ipw,': l =', &
                 paw_sp(isp)%l_pw(ipw)
         end do
         if (paw_sp(isp)%tcore_charge) write(stdout,'(a)') &
              '  Core charge supplied for Nonlinear Core Corrections'

      end do

      write(stdout,'(a/)') '<<<<<<<<<<<<<<<<<<<<<<<<<<&
           &<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

    end subroutine internal_print_datasets

  end subroutine paw_read_species


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_read_dataset(species)

    !==================================================================!
    ! This subroutine reads one PAW dataset file into the PAW_SPECIES  !
    ! type 'species'.                                                  !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  species (inout)  : Dataset for the PAW species to load.         !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 17/05/2010.                          !
    !==================================================================!

    use services, only: services_open_pspot_file
    use utils, only: utils_abort, utils_unit, utils_alloc_check, &
         utils_close_unit_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: species

    ! Local Variables
    character(len=100) :: dummy,dummy2
    integer :: ierr
    integer :: iunit
    integer :: i

    integer :: ipw, ipw_tot
    integer :: mi,li

    ! data to read in
    character(len=32) :: pspfmt
    real(kind=DP) :: r2well
    real(kind=DP) :: atnum_real
    integer :: pspdat, pspcod, pspxc, lmax, lloc, mmax
    integer :: creatorID
    integer :: igrid, jgrid

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_read_dataset'

    ! ndmh: open the file
    iunit = utils_unit()
    call services_open_pspot_file(iunit,trim(species%ds_name),ierr)

    ! Read comment and header
    read(iunit,*) dummy
    read(iunit,*) atnum_real,species%ion_charge,pspdat
    species%atomic_number = int(atnum_real)
    read(iunit,*) pspcod,pspxc,lmax,lloc,mmax,r2well
    read(iunit,*) pspfmt,creatorID
    read(iunit,*) species%npw,species%npw_tot

    ! Allocate arrays to store angular momentum values of each partial wave
    ! numbered either by ni,li or by ni,li,mi. Also store the ni,li value
    ! corresponding to each ni,li,mi value (in ipwtot)
    allocate(species%l_pw(species%npw),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%l_pw',ierr)
    allocate(species%ipw_tot(species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%ipw_tot',ierr)
    allocate(species%l_pw_tot(species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%l_pw_tot',ierr)
    allocate(species%m_pw_tot(species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%m_pw_tot',ierr)

    read(iunit,*,iostat=ierr) (species%l_pw(i),i=1,species%npw)
    species%lmax = maxval(species%l_pw(:))
    species%lmax = max(species%lmax*2,species%lmax+1)

    ! Fill in values of ipw_tot, l_pw_tot, m_pw_tot
    ipw_tot = 1
    do ipw=1,species%npw
       li = species%l_pw(ipw)
       do mi=-li,li
          species%ipw_tot(ipw_tot) = ipw
          species%l_pw_tot(ipw_tot) = li
          species%m_pw_tot(ipw_tot) = mi
          ipw_tot = ipw_tot + 1
       end do
    end do

    ! Read number of grids
    read(iunit,*) species%ngrid
    allocate(species%grid(species%ngrid),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%grid',ierr)

    ! Read in and initialise the grids
    call paw_read_grids(species%ngrid,species%grid,iunit)

    ! Read augmentation region radius
    read(iunit,*) species%rcut

    ! Read compensation density shape type and radius
    read(iunit,*) species%shape%shape_type,species%shape%rshape
    if (species%shape%rshape<1e-8_DP) species%shape%rshape = species%rcut

    ! Read the AE partial Waves
    do ipw=1,species%npw

       read(iunit,*) dummy, dummy2
       if (index(dummy2,'PHI')==0) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Expected "PHI"')
       end if

       read(iunit,*) igrid
       if (ipw==1) then
          species%phi_grid = igrid
          allocate(species%phi_rad(species%grid(igrid)%npt,species%npw), &
               stat=ierr)
          call utils_alloc_check('paw_read_dataset','species%phi_rad',ierr)
          allocate(species%tphi_rad(species%grid(igrid)%npt,species%npw), &
               stat=ierr)
          call utils_alloc_check('paw_read_dataset','species%tphi_rad',ierr)
       else
          if (igrid/=species%phi_grid) then
             call utils_abort('Error in paw_read_dataset: Non-matching grids &
                  &for partial waves in '//trim(species%ds_name))
          end if
       end if

       read(iunit,*) (species%phi_rad(i,ipw),i=1,species%grid(igrid)%npt)

    end do

    ! Read the PS partial waves
    do ipw=1,species%npw

       read(iunit,*) dummy, dummy2
       if (index(dummy2,'TPHI')==0) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Expected "TPHI"')
       end if

       read(iunit,*) igrid
       if (igrid/=species%phi_grid) then
          call utils_abort('Error in paw_read_dataset: Non-matching grids &
               &for partial waves in '//trim(species%ds_name))
       end if

       read(iunit,*) (species%tphi_rad(i,ipw),i=1,species%grid(igrid)%npt)

    end do

    ! Read the projectors
    species%proj_rcut = -1.0_DP
    do ipw=1,species%npw

       read(iunit,*) dummy, dummy2
       if (index(dummy2,'TPROJ')==0) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Expected "TPROJ"')
       end if
       read(iunit,*) igrid

       if (ipw==1) then
          species%proj_grid = igrid
          allocate(species%tproj_rad(species%grid(igrid)%npt,species%npw), &
               stat=ierr)
          call utils_alloc_check('paw_read_dataset','species%tproj_rad',ierr)
       else
          if (igrid/=species%proj_grid) then
             call utils_abort('Error in paw_read_dataset: Non-matching grids &
                  &for projectors in '//trim(species%ds_name))
          end if
       end if

       read(iunit,*) (species%tproj_rad(i,ipw),i=1,species%grid(igrid)%npt)

    end do
    species%proj_rcut = species%rcut

    ! Find the AE Core density
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'CORE_DENSITY')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "CORE_DENSITY"')
    end if
    read(iunit,*) igrid
    species%core_den_grid = igrid

    ! Check that the Core density grid has more or equal points to the phi grid
    ! and that it has the same spacings
    jgrid = species%phi_grid
    if (species%grid(igrid)%npt<species%grid(jgrid)%npt) then
       call utils_abort('Error in paw_read_dataset: Core density grid has &
            &fewer points than partial wave grid')
    end if
    if (species%grid(igrid)%rad_step/=species%grid(jgrid)%rad_step) then
       call utils_abort('Error in paw_read_dataset: Core density grid has &
            &different spacings from partial wave grid')
    end if
    if (species%grid(igrid)%log_step/=species%grid(jgrid)%log_step) then
       call utils_abort('Error in paw_read_dataset: Core density grid has &
            &different spacings from partial wave grid')
    end if

    ! Allocate storage for core densities
    allocate(species%core_den_rad(species%grid(igrid)%npt),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%core_den_rad',ierr)
    allocate(species%tcore_den_rad(species%grid(igrid)%npt),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%tcore_den_rad',ierr)
    read(iunit,*) (species%core_den_rad(i),i=1,species%grid(igrid)%npt)

    ! Read the Pseudo Core density
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'CORE_DENSITY')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "PSEUDO_CORE_DENSITY" or "TCORE_DENSITY"')
    end if
    read(iunit,*) igrid
    if (igrid/=species%core_den_grid) then
       call utils_abort('Error in paw_read_dataset: Non-matching grids &
            &for core densities in '//trim(species%ds_name))
    end if
    read(iunit,*) (species%tcore_den_rad(i),i=1,species%grid(igrid)%npt)
    if (any(species%tcore_den_rad/=0.0_DP)) then
       species%tcore_charge = .true.
    else
       species%tcore_charge = .false.
    end if

    ! Read Dij0
    allocate(species%dij0(species%npw_tot,species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%dij0',ierr)
    species%dij0 = 0.0_DP
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'Dij0')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "Dij0"')
    end if
    do ipw=1,species%npw_tot
       read(iunit,*) (species%dij0(i,ipw),i=1,ipw)
    end do
    do ipw=1,species%npw_tot
       do i=ipw+1,species%npw_tot
          species%dij0(i,ipw) = species%dij0(ipw,i)
       end do
    end do

    ! Read rhoij0
    allocate(species%rhoij0(species%npw_tot,species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%rhoij0',ierr)
    species%rhoij0 = 0.0_DP
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'Rhoij0')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "Rhoij0"')
    end if
    do ipw=1,species%npw_tot
       read(iunit,*) (species%rhoij0(i,ipw),i=1,ipw)
    end do
    do ipw=1,species%npw_tot
       do i=ipw+1,species%npw_tot
          species%rhoij0(i,ipw) = species%rhoij0(ipw,i)
       end do
    end do

    ! Read the vhntzc potential
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'VHntZC')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "VHntZC"')
    end if
    read(iunit,'(a)') dummy
    if (index(pspfmt,'paw5')>0) then
       read(dummy,*) igrid, species%vhntzc_format
       if ((species%vhntzc_format<1).or.(species%vhntzc_format>2)) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "vhntzc format" to be 1 or 2')
       end if
    else
       species%vhntzc_format = 1
       read(dummy,*) igrid
    end if
    species%vhntzc_grid = igrid
    allocate(species%vhntzc_rad(species%grid(igrid)%npt),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%vhntzc_rad',ierr)
    read(iunit,*) (species%vhntzc_rad(i),i=1,species%grid(igrid)%npt)

    ! ndmh: close the file
    close(iunit,iostat=ierr)
    call utils_close_unit_check('paw_read_dataset',trim(species%ds_name),ierr)

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_read_dataset'

  end subroutine paw_read_dataset


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_read_core_wfs(species)

    !==================================================================!
    ! This subroutine reads one PAW dataset file into the PAW_SPECIES  !
    ! type 'species'.                                                  !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  species (inout)  : Dataset for the PAW species to load.         !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 17/05/2010.                          !
    !==================================================================!

    use services, only: services_open_pspot_file
    use utils, only: utils_abort, utils_unit, utils_alloc_check, &
         utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: species

    ! Local Variables
    character(len=100) :: dummy,dummy2
    integer :: ierr
    integer :: iunit
    integer :: i

    integer :: icore_wf, icore_wf_tot
    integer :: mi,li

    ! data to read in
    character(len=32) :: pspfmt
    integer :: method,nspinor,nsppol
    integer :: pspdat, pspcod, pspxc, lmax, lloc
    integer :: creatorID
    integer :: igrid
    real(kind=DP) :: atnum_real

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_read_core_wfs'

    ! ndmh: open the file
    iunit = utils_unit()
    call services_open_pspot_file(iunit,trim(species%core_wf_name),ierr)

    ! Read comment and header
    read(iunit,*) dummy
    read(iunit,*) method,nspinor,nsppol
    read(iunit,*) atnum_real,species%core_wf_charge,pspdat
    if (species%atomic_number /= int(atnum_real)) then
       call utils_abort('Error in paw_read_core_wfs: Atomic number of core &
            &wvfn file does not match PAW dataset')
    end if
    read(iunit,*) pspcod,pspxc,lmax
    read(iunit,*) pspfmt,creatorID
    read(iunit,*) species%n_core_wfs,species%n_core_wfs_tot

    allocate(species%l_core_wf(species%n_core_wfs),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%l_core_wf',ierr)
    allocate(species%l_core_wf_tot(species%n_core_wfs_tot),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%l_core_wf_tot',ierr)
    allocate(species%m_core_wf_tot(species%n_core_wfs_tot),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%m_core_wf_tot',ierr)
    allocate(species%icore_wf_tot(species%n_core_wfs_tot),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%icore_wf_tot',ierr)
    allocate(species%n_core_wf(species%n_core_wfs),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%n_core_wf',ierr)

    read(iunit,*,iostat=ierr) (species%l_core_wf(i),i=1,species%n_core_wfs)

    ! Fill in values of icore_wf_tot, l_core_wf_tot, m_core_wf_tot
    icore_wf_tot = 1
    do icore_wf=1,species%n_core_wfs
       li = species%l_core_wf(icore_wf)
       do mi=-li,li
          species%icore_wf_tot(icore_wf_tot) = icore_wf
          species%l_core_wf_tot(icore_wf_tot) = li
          species%m_core_wf_tot(icore_wf_tot) = mi
          icore_wf_tot = icore_wf_tot + 1
       end do
    end do

    if (icore_wf_tot-1/=species%n_core_wfs_tot) then
       call utils_abort('Error reading dataset '//trim(species%core_wf_name) &
            //': Core wavefunction total counts do not match')
    end if

    ! Read number of grids
    read(iunit,*) species%ngrid_core
    allocate(species%grid_core(species%ngrid_core),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%grid_core',ierr)

    ! Read and initialise the grids
    call paw_read_grids(species%ngrid_core,species%grid_core,iunit)

    ! Read augmentation region radius
    read(iunit,*) species%rcut_core

    ! Read the Core wavefunctions
    species%core_wf_grid = 0
    do icore_wf=1,species%n_core_wfs

       read(iunit,*) dummy, dummy2
       if (index(dummy2,'Core')==0) then
          call utils_abort('Error reading dataset '//trim(species%core_wf_name) &
               //': Expected "Core"')
       end if

       ! Read in grid and check it matches previous grids
       read(iunit,*) igrid
       if (icore_wf==1) then
          species%core_wf_grid = igrid
          allocate(species%core_wf_rad(species%grid_core(igrid)%npt, &
               species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%core_wf_rad',ierr)
          allocate(species%core_wf_eig(species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%core_wf_eig',ierr)
          allocate(species%core_wf_occ(species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%core_wf_occ',ierr)
       else
          if (igrid/=species%core_wf_grid) then
             call utils_abort('Error in paw_read_core_wfs: Non-matching grids &
                  &for core wavefunctions in '//trim(species%core_wf_name))
          end if
       end if

       read(iunit,*) species%n_core_wf(icore_wf), lloc, nsppol
       if (lloc/=species%l_core_wf(icore_wf)) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Core wavefunction angular momenta do not match')
       end if
       if (nsppol/=1) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Unexpected value of spin')
       end if

       read(iunit,*) species%core_wf_eig(icore_wf),species%core_wf_occ(icore_wf)

       read(iunit,*) (species%core_wf_rad(i,icore_wf),i=1, &
            species%grid_core(igrid)%npt)

    end do

    ! ndmh: close the file
    close(iunit,iostat=ierr)
    call utils_close_unit_check('paw_read_core_wfs', &
         trim(species%core_wf_name),ierr)

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_read_core_wfs'

  end subroutine paw_read_core_wfs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_read_grids(ngrid,grids,iunit)

    !==================================================================!
    ! This subroutine reads and initialises grids from PAW datasets    !
    ! and core wavefunction files.                                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  ngrid (in)     : Number of grids to be read                     !
    !  grids (inout)  : RADIAL_GRID type (already allocated) to create !
    !  iunit (in)     : Unit number to read from                       !
    !------------------------------------------------------------------!
    ! Originally ritten by Nicholas Hine on 17/05/2010 as part of      !
    ! paw_read_species.                                                !
    ! Moved to its own routine on 25/10/2011 so that it can also be    !
    ! used by paw_read_core_wfs                                        !
    !==================================================================!

    use utils, only: utils_abort, utils_alloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: ngrid
    type(RADIAL_GRID), intent(inout) :: grids(ngrid)
    integer, intent(in) :: iunit

    ! Local Variables
    integer :: grid_type
    integer :: igrid,jgrid
    integer :: ierr
    integer :: i
    character(len=100) :: dummy

    ! Read in and initialise the grids
    do igrid=1,ngrid

       read(iunit,'(a)') dummy
       read(dummy,*) jgrid, grid_type

       select case (grid_type)
       case (1)
          grids(igrid)%log_step = 0.0_DP
          read(dummy,*) jgrid, grid_type, grids(igrid)%npt, &
               grids(igrid)%rad_step
       case (2)
          read(dummy,*) jgrid, grid_type, grids(igrid)%npt, &
               grids(igrid)%rad_step, grids(igrid)%log_step
       case default
          call utils_abort('Error in paw_read_grids: Unsupported grid type')
       end select

       allocate(grids(igrid)%r(grids(igrid)%npt),stat=ierr)
       call utils_alloc_check('paw_read_grids','grids(igrid)%r',ierr)
       allocate(grids(igrid)%rab(grids(igrid)%npt),stat=ierr)
       call utils_alloc_check('paw_read_grids','grids(igrid)%rab', &
            ierr)

       select case (grid_type)
       case (1)
          ! Regular grid
          do i=1,grids(igrid)%npt
             grids(igrid)%r(i) = real(i-1,kind=DP)*grids(igrid)%rad_step
             grids(igrid)%rab(i) = grids(igrid)%rad_step
          end do
       case (2)
          ! Logarithmic grid
          do i=1,grids(igrid)%npt
             grids(igrid)%r(i) = (exp(grids(igrid)%log_step &
                  *real(i-1,kind=DP))-1.0_DP)*grids(igrid)%rad_step
             grids(igrid)%rab(i) = grids(igrid)%log_step* &
                  (grids(igrid)%r(i)+grids(igrid)%rad_step)
          end do
       case default
          call utils_abort('Error in paw_read_grids: Unsupported grid type')
       end select

    end do

  end subroutine paw_read_grids


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_bcast_dataset(species)

    !==================================================================!
    ! This subroutine broadcasts the contents of one PAW_SPECIES type  !
    ! from the root proc (which read it in) to all the other procs.    !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  species (inout) : Dataset for the PAW species to share between  !
    !                    procs.                                        !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 18/05/2010.                          !
    !==================================================================!

    use comms, only: comms_bcast, pub_root_proc_id, pub_on_root
    use utils, only: utils_abort, utils_alloc_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: species

    ! Local Variables
    integer :: ierr
    integer :: npt

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_bcast_dataset'

    call comms_bcast(pub_root_proc_id,species%atomic_number)
    call comms_bcast(pub_root_proc_id,species%ion_charge)
    call comms_bcast(pub_root_proc_id,species%npw)
    call comms_bcast(pub_root_proc_id,species%npw_tot)

    if (.not.pub_on_root) then
       allocate(species%l_pw(species%npw),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%l_pw',ierr)
       allocate(species%ipw_tot(species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%ipw_tot',ierr)
       allocate(species%l_pw_tot(species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%l_pw_tot',ierr)
       allocate(species%m_pw_tot(species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%m_pw_tot',ierr)
    end if
    call comms_bcast(pub_root_proc_id,species%l_pw,species%npw)
    call comms_bcast(pub_root_proc_id,species%ipw_tot,species%npw_tot)
    call comms_bcast(pub_root_proc_id,species%l_pw_tot,species%npw_tot)
    call comms_bcast(pub_root_proc_id,species%m_pw_tot,species%npw_tot)
    call comms_bcast(pub_root_proc_id,species%lmax)

    call comms_bcast(pub_root_proc_id,species%ngrid)
    if (.not.pub_on_root) then
       allocate(species%grid(species%ngrid),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%grid',ierr)
    end if
    call paw_bcast_grids(species%ngrid,species%grid)

    call comms_bcast(pub_root_proc_id,species%rcut)
    call comms_bcast(pub_root_proc_id,species%proj_rcut)
    call comms_bcast(pub_root_proc_id,species%shape%shape_type)
    call comms_bcast(pub_root_proc_id,species%shape%rshape)

    call comms_bcast(pub_root_proc_id,species%phi_grid)
    call comms_bcast(pub_root_proc_id,species%proj_grid)
    call comms_bcast(pub_root_proc_id,species%core_den_grid)
    call comms_bcast(pub_root_proc_id,species%vhntzc_grid)

    if (.not.pub_on_root) then
       allocate(species%phi_rad(species%grid(species%phi_grid)%npt, &
            species%npw),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%phi_rad',ierr)
       allocate(species%tphi_rad(species%grid(species%phi_grid)%npt, &
            species%npw),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%tphi_rad',ierr)
       allocate(species%tproj_rad(species%grid(species%proj_grid)%npt, &
            species%npw),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%tproj_rad',ierr)
       allocate(species%core_den_rad(species%grid(species%core_den_grid)%npt), &
            stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%core_den_rad',ierr)
       allocate(species%tcore_den_rad(species%grid(species%core_den_grid)%npt),&
            stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%tcore_den_rad',ierr)
    end if
    call comms_bcast(pub_root_proc_id,species%phi_rad)
    call comms_bcast(pub_root_proc_id,species%tphi_rad)
    call comms_bcast(pub_root_proc_id,species%tproj_rad)
    call comms_bcast(pub_root_proc_id,species%core_den_rad)
    call comms_bcast(pub_root_proc_id,species%tcore_den_rad)
    call comms_bcast(pub_root_proc_id,species%tcore_charge)

    if (.not.pub_on_root) then
       allocate(species%dij0(species%npw_tot,species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%dij0',ierr)
       allocate(species%rhoij0(species%npw_tot,species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%rhoij0',ierr)
    end if
    call comms_bcast(pub_root_proc_id,species%dij0)
    call comms_bcast(pub_root_proc_id,species%rhoij0)

    call comms_bcast(pub_root_proc_id,species%vhntzc_format)
    if (.not.pub_on_root) then
       allocate(species%vhntzc_rad(species%grid(species%vhntzc_grid)%npt), &
            stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%vhntzc_rad',ierr)
    end if
    call comms_bcast(pub_root_proc_id,species%vhntzc_rad)

    ! Broadcast core WF info
    call comms_bcast(pub_root_proc_id,species%core_wvfns_exist) ! lr408
    call comms_bcast(pub_root_proc_id,species%n_core_wfs)
    call comms_bcast(pub_root_proc_id,species%n_core_wfs_tot)
    call comms_bcast(pub_root_proc_id,species%ngrid_core)

    if (species%core_wvfns_exist) then
       if (.not.pub_on_root) then
          allocate(species%l_core_wf(species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%l_core_wf',ierr)
          allocate(species%l_core_wf_tot(species%n_core_wfs_tot),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%l_core_wf_tot', ierr)
          allocate(species%m_core_wf_tot(species%n_core_wfs_tot),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%m_core_wf_tot',ierr)
          allocate(species%icore_wf_tot(species%n_core_wfs_tot),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%icore_wf_tot',ierr)
          allocate(species%n_core_wf(species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%n_core_wf',ierr)
          allocate(species%grid_core(species%ngrid_core),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%grid_core',ierr)
       end if

       call comms_bcast(pub_root_proc_id,species%l_core_wf)
       call comms_bcast(pub_root_proc_id,species%l_core_wf_tot)
       call comms_bcast(pub_root_proc_id,species%m_core_wf_tot)
       call comms_bcast(pub_root_proc_id,species%icore_wf_tot)
       call comms_bcast(pub_root_proc_id,species%n_core_wf)
       call paw_bcast_grids(species%ngrid_core,species%grid_core)
       call comms_bcast(pub_root_proc_id,species%core_wf_grid)
       call comms_bcast(pub_root_proc_id,species%rcut_core)

       if (species%core_wf_grid>0) then
          npt = species%grid_core(species%core_wf_grid)%npt
       else
          npt = 0
       end if

       if (.not.pub_on_root) then
          allocate(species%core_wf_rad(npt,species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%core_wf_rad',ierr)
          allocate(species%core_wf_eig(species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%core_wf_eig', ierr)
          allocate(species%core_wf_occ(species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%core_wf_occ',ierr)
       end if
       call comms_bcast(pub_root_proc_id,species%core_wf_rad)
       call comms_bcast(pub_root_proc_id,species%core_wf_eig)
       call comms_bcast(pub_root_proc_id,species%core_wf_occ)
    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_bcast_dataset'

  end subroutine paw_bcast_dataset


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_bcast_grids(ngrid,grids)

    !==================================================================!
    ! This subroutine broadcasts grids from PAW datasets to all procs. !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  ngrid (in)     : Number of grids to be bcast                    !
    !  grids (inout)  : RADIAL_GRID type (already allocated) to bcast. !
    !------------------------------------------------------------------!
    ! Originally ritten by Nicholas Hine on 17/05/2010 as part of      !
    ! paw_bcast_dataset.                                               !
    ! Moved to its own routine on 25/10/2011 so that it can also be    !
    ! used for core wf grids .                                         !
    !==================================================================!

    use comms, only: comms_bcast, pub_on_root, pub_root_proc_id
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: ngrid
    type(RADIAL_GRID), intent(inout) :: grids(ngrid)

    ! Local Variables
    integer :: igrid
    integer :: ierr

    do igrid=1,ngrid
       call comms_bcast(pub_root_proc_id,grids(igrid)%npt)
       call comms_bcast(pub_root_proc_id,grids(igrid)%rad_step)
       call comms_bcast(pub_root_proc_id,grids(igrid)%log_step)
       if (.not.pub_on_root) then
          allocate(grids(igrid)%r(grids(igrid)%npt),stat=ierr)
          call utils_alloc_check('paw_read_dataset','grids(igrid)%r', &
               ierr)
          allocate(grids(igrid)%rab(grids(igrid)%npt),stat=ierr)
          call utils_alloc_check('paw_read_dataset','grids(igrid)%rab', &
               ierr)
       end if

       call comms_bcast(pub_root_proc_id,grids(igrid)%r)
       call comms_bcast(pub_root_proc_id,grids(igrid)%rab)
    end do

  end subroutine paw_bcast_grids


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dataset_init(species,core_wvfns_exist)

    !==================================================================!
    ! This subroutine allocates and initialises those arrays in the    !
    ! PAW_SPECIES type which are not loaded in directly but must be    !
    ! calculated, but which are nevertheless independent of any        !
    ! quantities depending on the system.                              !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  species (inout) : Dataset for the PAW species to calculate all  !
    !                    quantities required which are not stored in   !
    !                    the dataset file.                             !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/05/2010.                          !
    !==================================================================!

    use gaunt_coeff, only: realgaunt
    use paw_shape, only: paw_shape_init
    use services, only: services_radial_transform, services_locate_interp, &
         services_radial_integral, services_radial_integral_rmax, &
         services_radial_derivative
    use rundat, only: pub_nlcc
!$  use rundat, only: pub_threads_max
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_erf, utils_assert

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: species
    logical, intent(in) :: core_wvfns_exist ! lr408

    ! Local Variables
    integer :: ierr
    integer :: igrid
    integer :: iq,ir,ir_c
    integer :: ipw,jpw,kpw,lpw,ipwtot,jpwtot,kpwtot,lpwtot
    integer :: li,mi,lj,mj,lk,mk,ll,ml
    integer :: lup, mup
    integer :: npts
    integer :: nptsc, igridc ! lr408
    real(kind=DP) :: q, r
    real(kind=DP) :: rcmax, r2new
    real(kind=DP) :: rgij, rgkl
    real(kind=DP) :: int1, int2, int3, int4
    real(kind=DP) :: lfac
    real(kind=DP), allocatable :: work(:),work2(:),work3(:)
    real(kind=DP), allocatable :: inter1(:),inter2(:),inter3(:),inter4(:)
    real(kind=DP), allocatable :: phir_phjr(:),tphir_tphjr(:)
    real(kind=DP), allocatable :: phkr_phlr(:),tphkr_tphlr(:)
    real(kind=DP), allocatable :: rwork(:,:)

    ! ndmh: temporary arrays
    real(kind=DP), allocatable :: v_shape_L(:,:)
    real(kind=DP), allocatable :: int_v_L_g_L(:)
    real(kind=DP), allocatable :: Vhat_L_ij(:,:,:)
    real(kind=DP), allocatable :: V_L_ijkl(:,:,:,:,:)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_dataset_init'

    ! Start timer
    call timer_clock('paw_dataset_init',1)

    ! Prepare reciprocal space grid
    species%n_recip_pts = 3001
    species%g_max = 50.0_DP !13.145341380124_DP
    species%inv_g_spacing = real(species%n_recip_pts-1,kind=DP)/species%g_max

    ! Fix all grid sizes to be odd
    species%grid(:)%npt = species%grid(:)%npt + modulo(species%grid(:)%npt,2) -1

    ! Find interpolation points on phi_grid for integrations up to r_c
    igrid = species%phi_grid
    npts = species%grid(igrid)%npt
    ir_c = services_locate_interp(species%rcut,species%grid(igrid)%r,npts)
    call utils_assert(ir_c <= npts-3, 'Error in paw_dataset_init: Not enough &
         &grid points beyond PAW radius on partial wave grid for accurate &
         &integration of partial waves')

    !=========================================================================!
    ! PREPARE SHAPE FUNCTION IN REAL SPACE                                    !
    !=========================================================================!

    ! Initialise the shape functions
    igrid = species%phi_grid
    species%shape_grid = igrid
    call paw_shape_init(species%shape,species%grid(igrid)%r, &
         species%grid(igrid)%rab,species%grid(igrid)%npt,species%lmax)

    !=========================================================================!
    ! PREPARE PROJECTORS IN RECIPROCAL SPACE                                  !
    !=========================================================================!

    allocate(species%tproj_recip(species%n_recip_pts,species%npw),stat=ierr)
    call utils_alloc_check('paw_dataset_init','species%tproj_recip',ierr)

    igrid = species%proj_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    do ipw=1,species%npw
       call services_radial_transform(species%l_pw(ipw), 1, npts,  &
            species%grid(igrid)%r,species%grid(igrid)%rab, &
            species%n_recip_pts,species%g_max,species%tproj_rad(:,ipw), &
            species%tproj_recip(:,ipw))
    end do

    !=========================================================================!
    ! PREPARE HARTREE POTENTIAL OF CORE DENSITY n_Zc                          !
    !=========================================================================!

    allocate(species%vhntzc_recip(species%n_recip_pts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','species%vhntzc_recip',ierr)

    igrid = species%vhntzc_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    rcmax = 1.0_dp/sqrt(min(-100.0_dp/(4.0_dp*log(1.d-12)),0.7_dp))

    allocate(work(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work',ierr)
    allocate(work2(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work2',ierr)
    work(1) = 0.0_DP ! jd
    do ir=2,npts
       r = species%grid(igrid)%r(ir)
       work(ir) = species%vhntzc_rad(ir)*r + &
            species%ion_charge*utils_erf(r/rcmax)
    end do
    r2new = 0.25_dp*(rcmax)**2
    species%vhntzc_recip(1) = services_radial_integral(npts, &
         species%grid(igrid)%rab,work*species%grid(igrid)%r(1:npts)) &
         + species%ion_charge*r2new
    do iq=2,species%n_recip_pts
       q = real(iq-1,dp)/real(species%n_recip_pts-1,dp)*species%g_max
       do ir=1,npts
          work2(ir)=work(ir)*sin(species%grid(igrid)%r(ir)*q)
       enddo
       species%vhntzc_recip(iq) = &
            services_radial_integral(npts,species%grid(igrid)%rab,work2)/q &
            + species%ion_charge*(1.0_dp-exp(-r2new*q*q))/(q**2)
    end do
    deallocate(work2,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work2',ierr)
    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work',ierr)

    !=========================================================================!
    ! PREPARE TRANSFORM OF PSEUDO CORE DENSITY \tilde{n}_c                    !
    !=========================================================================!

    if (species%tcore_charge) then

       ! Set global flag if not yet already true
       pub_nlcc = .true.

       ! Storage for reciprocal space core density
       allocate(species%tcore_den_recip(species%n_recip_pts),stat=ierr)
       call utils_alloc_check('paw_dataset_init','species%tcore_den_recip',ierr)

       igrid = species%core_den_grid
       npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1

       ! Transform the real-space core charge to reciprocal space
       call services_radial_transform(0,2,npts,species%grid(igrid)%r, &
            species%grid(igrid)%rab,species%n_recip_pts,&
            species%g_max,species%tcore_den_rad,species%tcore_den_recip)

       species%tcore_den_recip = species%tcore_den_recip * 4.0_DP * PI

    else
       species%tcore_charge = .false.
    end if

    if (any(species%core_den_rad/=0.0_DP)) then
       species%core_charge = .true.
       species%core_charge_calculated = .false.
    else
       species%core_charge = .false.
       species%core_charge_calculated = .false.
    end if

    !=========================================================================!
    ! PREPARE n^L_{n_i l_i n_j l_j} ARRAY                                     !
    !=========================================================================!

    ! This needs to be stored in the species array, for calculating nhat later
    allocate(species%aug_nLij(species%npw,species%npw,-2:species%lmax), &
         stat=ierr)
    call utils_alloc_check('paw_dataset_init','species%aug_nLij',ierr)

    igrid = species%phi_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(work,inter1,rwork,ierr,lup,ipw,jpw) &
!$OMP SHARED(species,igrid,npts,ir_c,pub_threads_max)
    allocate(work(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work',ierr)
    allocate(rwork(npts,2),stat=ierr)
    call utils_alloc_check('paw_dataset_init','rwork',ierr)
    ! sum runs from -1 so that we get the (phi phj-tphi tphj)/r term for
    ! use in the grad operator
    do lup=-1,species%lmax
       ! check for divide by zero
       if ((lup<0).and.(species%grid(igrid)%r(1)==0.0_DP)) then
          rwork(1,1) = 0.0_DP
          rwork(2:npts,1) = species%grid(igrid)%r(2:npts)**lup
       else
          rwork(1:npts,1) = species%grid(igrid)%r(1:npts)**lup
       end if
!$OMP DO
       do ipw=1,species%npw
          do jpw=1,ipw
             work(1:npts) = (species%phi_rad(1:npts,ipw)*species%phi_rad(1:npts,jpw) - &
                  species%tphi_rad(1:npts,ipw)*species%tphi_rad(1:npts,jpw)) * &
                  rwork(1:npts,1)
             if (modulo(ir_c,2)==1) then
                species%aug_nLij(ipw,jpw,lup) = services_radial_integral( &
                     ir_c,species%grid(igrid)%rab(1:),work(1:))
             else
                species%aug_nLij(ipw,jpw,lup) = services_radial_integral( &
                     ir_c-1,species%grid(igrid)%rab(2:),work(2:))
             end if

             species%aug_nLij(jpw,ipw,lup) = species%aug_nLij(ipw,jpw,lup)
          end do
       end do
!$OMP END DO
    end do
    ! extra set of values, stored in L=-2, for use in grad operator
!$OMP DO
    do ipw=1,species%npw
       do jpw=1,species%npw
          call services_radial_derivative(rwork(:,1),species%phi_rad(:,jpw), &
               npts,real(npts,kind=DP))
          rwork(1:npts,1) = rwork(1:npts,1) / species%grid(igrid)%rab(1:npts)
          call services_radial_derivative(rwork(:,2),species%tphi_rad(:,jpw), &
               npts,real(npts,kind=DP))
          rwork(1:npts,2) = rwork(1:npts,2) / species%grid(igrid)%rab(1:npts)
          work(1:npts) = (species%phi_rad(1:npts,ipw)*rwork(1:npts,1) - &
               species%tphi_rad(1:npts,ipw)*rwork(1:npts,2))
          if (modulo(ir_c,2)==1) then
             species%aug_nLij(ipw,jpw,-2) = services_radial_integral( &
                  ir_c,species%grid(igrid)%rab(1:),work(1:))
          else
             species%aug_nLij(ipw,jpw,-2) = services_radial_integral( &
                  ir_c-1,species%grid(igrid)%rab(2:),work(2:))
          end if
       end do
    end do
!$OMP END DO

    deallocate(rwork,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','rwork',ierr)
    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work',ierr)
!$OMP END PARALLEL

    ! Allocate e_ijkl now, since subsequent temporary arrays will be
    ! deallocated at the end of the routine
    allocate(species%e_ijkl(species%npw_tot,species%npw_tot,species%npw_tot, &
         species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_dataset_init','species%e_ijkl',ierr)

    !=========================================================================!
    ! PREPARE v^L(r) FUNCTION                                                 !
    !=========================================================================!

    igrid = species%shape_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    allocate(v_shape_L(npts,0:species%lmax),stat=ierr)
    call utils_alloc_check('paw_dataset_init','v_shape_L',ierr)
    allocate(work(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work',ierr)
    allocate(work2(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work2',ierr)
    allocate(inter1(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','inter1',ierr)
    allocate(inter2(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','inter2',ierr)
    allocate(rwork(npts,4),stat=ierr)
    call utils_alloc_check('paw_dataset_init','rwork',ierr)
    do lup=0,species%lmax
       rwork(1:npts,1) = species%grid(igrid)%r(1:npts)**lup
       rwork(2:npts,2) = 1.0_DP/(species%grid(igrid)%r(2:npts)**(lup+1))
       rwork(1:npts,3) = species%grid(igrid)%r(1:npts)**(lup+2)
       rwork(2:npts,4) = species%grid(igrid)%r(2:npts)**(1-lup)
       if (species%grid(igrid)%r(1)==0.0_DP) then
          rwork(1,2) = 0.0_DP
          if (lup>1) then
             rwork(1,4) = 0.0_DP
          else
             rwork(1,4) = species%grid(igrid)%r(1)**(1-lup)
          end if
       else
          rwork(1,2) = 1.0_DP/(species%grid(igrid)%r(1)**(lup+1))
          rwork(1,4) = species%grid(igrid)%r(1)**(1-lup)
       end if
       ! work(:)  = g_L(r) * rp^L / r ^(L+1) * rp^2
       ! work2(:) = g_L(r) * r ^L / rp^(L+1) * rp^2
       work(:)  = rwork(:,3) * species%shape%shape_rad(:,lup)
       work2(:) = rwork(:,4) * species%shape%shape_rad(:,lup)
       int1 = services_radial_integral(npts,species%grid(igrid)%rab, &
            work,inter1)
       int2 = services_radial_integral(npts,species%grid(igrid)%rab, &
            work2,inter2)
       v_shape_L(:,lup) = inter1(:)* rwork(:,2) + (int2 - inter2(:)) * rwork(:,1)
       !do ir=1,npts
          !r = species%grid(igrid)%r(ir)
          !write(stdout,'(i5,7f20.12)') ir, r, species%shape%shape_rad(ir,lup), &
          !     inter1(ir),int2,inter2(ir),work2(ir),v_shape_L(ir,lup)
       !end do
       !write(stdout,*)
       v_shape_L(:,lup) = v_shape_L(:,lup) * 4.0_DP*PI / real(2*lup+1,kind=DP)
    end do
    deallocate(rwork,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','rwork',ierr)
    deallocate(inter2,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','inter2',ierr)
    deallocate(inter1,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','inter1',ierr)
    deallocate(work2,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work2',ierr)
    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work',ierr)

    !=========================================================================!
    ! PREPARE \int v^L(r) g_L(r) dr INTEGRALS                                 !
    !=========================================================================!

    igrid = species%shape_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    allocate(int_v_L_g_L(0:species%lmax),stat=ierr)
    call utils_alloc_check('paw_dataset_init','int_v_L_g_L',ierr)

    allocate(work(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work',ierr)
    do lup=0,species%lmax
       work(1:npts) = v_shape_L(1:npts,lup) * species%shape%shape_rad(1:npts,lup) &
            * species%grid(igrid)%r(1:npts)**2
       int_v_L_g_L(lup) = services_radial_integral(npts, &
            species%grid(igrid)%rab,work)
    end do
    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work',ierr)

    !=========================================================================!
    ! PREPARE \hat{V}^L_ij ARRAY                                              !
    !=========================================================================!

    igrid = species%phi_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    allocate(Vhat_L_ij(species%npw,species%npw,0:species%lmax),stat=ierr)
    call utils_alloc_check('paw_dataset_init','Vhat_L_ij',ierr)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(work,inter1,ierr,lup,ipw,jpw) &
!$OMP SHARED(species,igrid,npts,ir_c,v_shape_L,Vhat_L_ij,pub_threads_max)
    allocate(work(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work',ierr)
    allocate(inter1(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','inter1',ierr)
    do lup=0,species%lmax
!$OMP DO
       do ipw=1,species%npw
          do jpw=1,species%npw
             work(1:npts) = v_shape_L(1:npts,lup)*species%tphi_rad(1:npts,ipw)* &
                  species%tphi_rad(1:npts,jpw)
             Vhat_L_ij(ipw,jpw,lup) = &
                  services_radial_integral_rmax(npts,species%grid(igrid)%rab, &
                  species%grid(igrid)%r,species%rcut,work,inter1)
          end do
       end do
!$OMP END DO
    end do
    deallocate(inter1,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','inter1',ierr)
    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work',ierr)
!$OMP END PARALLEL

    !=========================================================================!
    ! PREPARE V^L_ijkl ARRAY                                                  !
    !=========================================================================!

    igrid = species%phi_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    allocate(V_L_ijkl(species%npw,species%npw,species%npw,species%npw, &
         0:species%lmax),stat=ierr)
    call utils_alloc_check('paw_dataset_init','V_L_ijkl',ierr)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(work,work2,work3,inter1,inter2,inter3,inter4,rwork,phir_phjr, &
!$OMP      tphir_tphjr,phkr_phlr,tphkr_tphlr,ierr,lup,ipw,jpw,lfac, &
!$OMP      int1,int2,int3,int4) &
!$OMP SHARED(species,igrid,npts,ir_c,V_L_ijkl,pub_threads_max)
    allocate(work(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work',ierr)
    allocate(work2(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work2',ierr)
    allocate(work3(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','work3',ierr)
    allocate(inter1(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','inter1',ierr)
    allocate(inter2(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','inter2',ierr)
    allocate(inter3(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','inter3',ierr)
    allocate(inter4(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','inter4',ierr)
    allocate(phir_phjr(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','phir_phjr',ierr)
    allocate(tphir_tphjr(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','tphir_tphjr',ierr)
    allocate(phkr_phlr(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','phkr_phlr',ierr)
    allocate(tphkr_tphlr(npts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','tphkr_tphlr',ierr)
    allocate(rwork(npts,2),stat=ierr)
    call utils_alloc_check('paw_dataset_init','rwork',ierr)

    ! jd: services_radial_integral_rmax() does not initialise all elements.
    !     Ensure there's no garbage in inter[1-4] before they are operated on.
    inter1(:) = 0.0_DP
    inter2(:) = 0.0_DP
    inter3(:) = 0.0_DP
    inter4(:) = 0.0_DP

    ! Loop over L and quadruple loop over partial waves
    do lup=0,species%lmax
       lfac = 4.0_DP * PI / real(2*lup+1,kind=DP)
       rwork(1:npts,1) = species%grid(igrid)%r(1:npts)**lup
       rwork(2:npts,2) = 1.0_DP / (species%grid(igrid)%r(2:npts)**(lup+1))
       if (species%grid(igrid)%r(1)==0.0_DP) then
          rwork(1,2) = 0.0_DP
       else
          rwork(1,2) = 1.0_DP / (species%grid(igrid)%r(1)**(lup+1))
       end if
!$OMP DO
       do ipw=1,species%npw
          do jpw=1,ipw
             phir_phjr(1:npts) = species%phi_rad(1:npts,ipw)*species%phi_rad(1:npts,jpw)
             tphir_tphjr(1:npts) = species%tphi_rad(1:npts,ipw)*species%tphi_rad(1:npts,jpw)

             ! ndmh_230316: new version of Hartree tensor calculation
             ! ndmh_230316: should be much faster

             ! Calculate Hartree potential of phi_i * phi_j

             ! Integrate (phi_i(r)/r).(phi_j(r)/r) * r^L+2
             ! from 0 to rc, storing intermediates in inter1
             work(:) = phir_phjr(:) * rwork(:,1)
             int1 = services_radial_integral_rmax(npts, &
                  species%grid(igrid)%rab,species%grid(igrid)%r, &
                  species%rcut,work,inter1)

             ! Integrate (tphi_i(r)/r).(tphi_j(r)/r) * r^L+2
             ! from 0 to rc, storing intermediates in inter3
             work(:) = tphir_tphjr(:) * rwork(:,1)
             int3 = services_radial_integral_rmax(npts, &
                  species%grid(igrid)%rab,species%grid(igrid)%r, &
                  species%rcut,work,inter3)

             ! Integrate (phi_i(r)/r).(phi_j(r)/r) * r^(1-L)
             ! from 0 to rc, storing intermediates in inter2
             work(:) = phir_phjr(:) * rwork(:,2)
             int2 = services_radial_integral_rmax(npts, &
                  species%grid(igrid)%rab,species%grid(igrid)%r, &
                  species%rcut,work,inter2)

             ! Integrate (tphi_i(r)/r).(tphi_j(r)/r) * r^(1-L)
             ! from 0 to rc, storing intermediates in inter4
             work(:) = tphir_tphjr(:) * rwork(:,2)
             int4 = services_radial_integral_rmax(npts, &
                  species%grid(igrid)%rab,species%grid(igrid)%r, &
                  species%rcut,work,inter4)

             do kpw=1,species%npw
                do lpw=1,kpw
                   phkr_phlr(1:npts) = species%phi_rad(1:npts,kpw) * &
                        species%phi_rad(1:npts,lpw)
                   tphkr_tphlr(1:npts) = species%tphi_rad(1:npts,kpw) * &
                        species%tphi_rad(1:npts,lpw)

                   ! inter1(r),inter3(r) are integrals from 0 to r
                   ! int2-inter2(r),int4-inter4(r) are integrals from r to rc
                   work(:) = inter1(:) * phkr_phlr(:) * rwork(:,2) - &
                             inter3(:) * tphkr_tphlr(:) * rwork(:,2) + &
                             (int2 - inter2(:)) * phkr_phlr(:) * rwork(:,1) - &
                             (int4 - inter4(:)) * tphkr_tphlr(:) * rwork(:,1)

                   ! integrate work3(r) from 0 to rc
                   V_L_ijkl(ipw,jpw,kpw,lpw,lup) = lfac * &
                        services_radial_integral_rmax(npts, &
                        species%grid(igrid)%rab,species%grid(igrid)%r, &
                        species%rcut,work,work2)

                end do  ! lpw
             end do  ! kpw

             ! Fill in remaining values from symmetric equivalents for this ij
             do kpw=1,species%npw
                do lpw=kpw+1,species%npw
                   V_L_ijkl(ipw,jpw,kpw,lpw,lup) = V_L_ijkl(ipw,jpw,lpw,kpw,lup)
                end do  ! lpw
             end do  ! kpw

          end do  ! jpw
       end do  ! ipw
!$OMP END DO
    end do  ! lup

    deallocate(rwork,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','rwork',ierr)
    deallocate(tphkr_tphlr,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','tphkr_tphlr',ierr)
    deallocate(phkr_phlr,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','phkr_phlr',ierr)
    deallocate(tphir_tphjr,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','tphir_tphjr',ierr)
    deallocate(phir_phjr,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','phir_phjr',ierr)
    deallocate(inter4,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','inter4',ierr)
    deallocate(inter3,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','inter3',ierr)
    deallocate(inter2,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','inter2',ierr)
    deallocate(inter1,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','inter1',ierr)
    deallocate(work3,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work3',ierr)
    deallocate(work2,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work2',ierr)
    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','work',ierr)
!$OMP END PARALLEL

    ! Fill in remaining values from symmetric equivalents
    do ipw=1,species%npw
       do jpw=ipw+1,species%npw
          V_L_ijkl(ipw,jpw,:,:,:) = V_L_ijkl(jpw,ipw,:,:,:)
       end do
    end do

    !=========================================================================!
    ! PREPARE e_ijkl TENSOR                                                   !
    !=========================================================================!

    species%e_ijkl(:,:,:,:) = 0.0_DP
!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(lup,mup,ipwtot,li,mi,ipw,jpwtot,lj,mj,jpw,kpwtot,lk,mk,kpw, &
!$OMP      lpwtot,ll,ml,lpw,rgij,rgkl) &
!$OMP SHARED(species,V_L_ijkl,Vhat_L_ij,int_v_L_g_L,pub_threads_max)
    do lup=0,species%lmax
       do mup=-lup,lup
!$OMP DO
          do ipwtot=1,species%npw_tot
             ipw = species%ipw_tot(ipwtot)
             li = species%l_pw_tot(ipwtot)
             mi = species%m_pw_tot(ipwtot)
             do jpwtot=1,species%npw_tot
                jpw = species%ipw_tot(jpwtot)
                lj = species%l_pw_tot(jpwtot)
                mj = species%m_pw_tot(jpwtot)
                rgij = realgaunt(lup,mup,li,mi,lj,mj)
                if (abs(rgij)<1e-16) cycle
                do kpwtot=1,species%npw_tot
                   kpw = species%ipw_tot(kpwtot)
                   lk = species%l_pw_tot(kpwtot)
                   mk = species%m_pw_tot(kpwtot)
                   do lpwtot=1,species%npw_tot
                      lpw = species%ipw_tot(lpwtot)
                      ll = species%l_pw_tot(lpwtot)
                      ml = species%m_pw_tot(lpwtot)
                      rgkl = realgaunt(lup,mup,lk,mk,ll,ml)
                      if (abs(rgkl)<1e-16) cycle
                      species%e_ijkl(ipwtot,jpwtot,kpwtot,lpwtot) = &
                           species%e_ijkl(ipwtot,jpwtot,kpwtot,lpwtot) &
                           + rgij*rgkl*(V_L_ijkl(ipw,jpw,kpw,lpw,lup) &
                           - species%aug_nLij(ipw,jpw,lup)*Vhat_L_ij(kpw,lpw,lup) &
                           - species%aug_nLij(kpw,lpw,lup)*Vhat_L_ij(ipw,jpw,lup) &
                           - species%aug_nLij(ipw,jpw,lup)*species%aug_nLij(kpw, &
                           lpw,lup)*int_v_L_g_L(lup))
                   end do
                end do
             end do
          end do
!$OMP END DO
       end do
    end do
!$OMP END PARALLEL

    do ipwtot=1,species%npw_tot
       do jpwtot=1,species%npw_tot
          do kpwtot=1,species%npw_tot
             do lpwtot=1,species%npw_tot
                !print '(4i4,f20.12)',ipwtot,jpwtot,kpwtot,lpwtot,species%e_ijkl(ipwtot,jpwtot,kpwtot,lpwtot)
             end do
          end do
       end do
    end do


    ! lr408: assuming same grid here so basically not changing that much
    if (core_wvfns_exist) then

       ! Fix all grid sizes to be odd
       species%grid_core(:)%npt = species%grid_core(:)%npt + modulo(species%grid_core(:)%npt,2) -1

       ! Find interpolation points on phi_grid for integrations up to r_c
       igrid = species%core_wf_grid
       npts = species%grid_core(igrid)%npt

       ir_c = services_locate_interp(species%rcut,species%grid_core(igrid)%r,npts)
       call utils_assert(ir_c <= npts-3, 'Error in paw_dataset_init: Not enough&
            &grid points beyond PAW radius on partial wave grid for accurate &
            & integration of partial waves')

       !=========================================================================!
       ! PREPARE CORE WAVEFUNCTIONS IN RECIPROCAL SPACE                          !
       !=========================================================================!

       allocate(species%core_wf_recip(species%n_recip_pts,species%n_core_wfs),stat=ierr)
       call utils_alloc_check('paw_dataset_init','species%core_wf_recip',ierr)
       igrid = species%core_wf_grid
       npts = species%grid_core(igrid)%npt + modulo(species%grid_core(igrid)%npt,2) - 1
       do ipw=1,species%n_core_wfs
          call services_radial_transform(species%l_core_wf(ipw), 1, npts,  &
               species%grid_core(igrid)%r,species%grid_core(igrid)%rab, &
               species%n_recip_pts,species%g_max,species%core_wf_rad(:,ipw), &
               species%core_wf_recip(:,ipw))
       end do
    end if

    if (core_wvfns_exist) then

       !=========================================================================!
       ! PREPARE n^L_{n_i l_i n_j l_j} ARRAY FOR CORE WAVEFUNCTIONS              !
       !=========================================================================!

       ! Check that grids have same spacing - fail if not
       if (species%grid(species%core_wf_grid)%rad_step /= &
            species%grid(species%phi_grid)%rad_step) then
          call utils_abort('Error in paw_read_dataset: Core wavefunction grid &
               &has different spacings from partial wave grid')
       end if
       if (species%grid(species%core_wf_grid)%log_step /= &
            species%grid(species%phi_grid)%log_step) then
          call utils_abort('Error in paw_read_dataset: Core wavefunction grid &
               &has different spacings from partial wave grid')
       end if

       ! This needs to be stored in the species array, for calculating nhat later

       ! unclear what to do with lmax - hopefully the same for core wvfns but may
       ! need to return to
       allocate(species%core_aug_nLij(species%n_core_wfs, &
            species%npw,0:species%lmax),stat=ierr)
       call utils_alloc_check('paw_dataset_init','species%core_aug_nLij',ierr)

       igrid = species%phi_grid
       igridc = species%core_wf_grid
       npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
       nptsc = species%grid_core(igridc)%npt + &
            modulo(species%grid_core(igridc)%npt,2) - 1
       npts = min(npts,nptsc)

       allocate(work(npts),stat=ierr)
       call utils_alloc_check('paw_dataset_init','work',ierr)
       allocate(inter1(npts),stat=ierr)
       call utils_alloc_check('paw_dataset_init','inter1',ierr)
       allocate(rwork(npts,1),stat=ierr)
       call utils_alloc_check('paw_dataset_init','rwork',ierr)
       do lup=0,species%lmax
          rwork(1:npts,1) = species%grid(igrid)%r(1:npts)**lup
          do ipw=1,species%n_core_wfs
             do jpw=1,species%npw
                work(1:npts) = (species%core_wf_rad(1:npts,ipw) * &
                     species%phi_rad(1:npts,jpw) - &
                     species%core_wf_rad(1:npts,ipw) * &
                     species%tphi_rad(1:npts,jpw)) * rwork(1:npts,1)
                if (modulo(ir_c,2)==1) then
                   species%core_aug_nLij(ipw,jpw,lup) = services_radial_integral( &
                        ir_c,species%grid(igrid)%rab(1:npts),work(1:npts))
                else
                   species%core_aug_nLij(ipw,jpw,lup) = services_radial_integral( &
                        ir_c-1,species%grid(igrid)%rab(2:npts),work(2:npts))
                end if

             end do
          end do
       end do
       deallocate(rwork,stat=ierr)
       call utils_dealloc_check('paw_dataset_init','rwork',ierr)
       deallocate(inter1,stat=ierr)
       call utils_dealloc_check('paw_dataset_init','inter1',ierr)
       deallocate(work,stat=ierr)
       call utils_dealloc_check('paw_dataset_init','work',ierr)

    end if

    deallocate(V_L_ijkl,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','V_L_ijkl',ierr)
    deallocate(Vhat_L_ij,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','Vhat_L_ij',ierr)
    deallocate(int_v_L_g_L,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','int_v_L_g_L',ierr)
    deallocate(v_shape_L,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','v_shape_L',ierr)

    ! Stop timer
    call timer_clock('paw_dataset_init',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_dataset_init'

  end subroutine paw_dataset_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_all_species_exit(paw_sp)

    !====================================================================!
    ! This subroutine deallocates the contents of the PAW species array, !
    ! deallocates the PAW projectors and cleans up storage for the Gaunt !
    ! Coefficients.                                                      !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  species (inout) : Dataset for the PAW species to share between    !
    !                    procs.                                          !
    !--------------------------------------------------------------------!
    ! Written by Nicholas Hine on 18/05/2010.                            !
    !====================================================================!

    use gaunt_coeff, only: gaunt_exit
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), pointer, intent(inout) :: paw_sp(:)
    ! Local Variables
    integer :: ierr
    integer :: isp
    integer :: num_pspecies

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_all_species_exit'

    call gaunt_exit

    ! Deallocate arrays in paw_sp type
    num_pspecies = size(paw_sp)
    do isp=num_pspecies,1,-1
       call paw_species_exit(paw_sp(isp))
    end do

    ! Deallocate paw_sp itself
    deallocate(paw_sp,stat=ierr)
    call utils_dealloc_check('paw_species_exit','paw_sp',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_all_species_exit'

  end subroutine paw_all_species_exit


  subroutine paw_species_exit(species)

    !====================================================================!
    ! This subroutine deallocates the contents of the PAW species array, !
    ! deallocates the PAW projectors and cleans up storage for the Gaunt !
    ! Coefficients.                                                      !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  species (inout) : Dataset for the PAW species to share between    !
    !                    procs.                                          !
    !--------------------------------------------------------------------!
    ! Written by Nicholas Hine on 18/05/2010.                            !
    !====================================================================!

    use paw_shape, only: paw_shape_exit
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: species

    ! Local Variables
    integer :: ierr
    integer :: igrid

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_species_exit'

    ! Deallocate reciprocal space arrays in species type
    deallocate(species%e_ijkl,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%e_ijkl',ierr)
    deallocate(species%aug_nLij,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%aug_nLij',ierr)
    if (allocated(species%tcore_den_recip)) then
       deallocate(species%tcore_den_recip,stat=ierr)
       call utils_dealloc_check('paw_species_exit', &
            'species%tcore_den_recip',ierr)
    end if
    deallocate(species%vhntzc_recip,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%vhntzc_recip',ierr)

    ! lr408
    if (species%core_wvfns_exist) then
       deallocate(species%core_wf_recip,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%core_wf_recip',ierr)
    end if

    deallocate(species%tproj_recip,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%tproj_recip',ierr)

    ! Deallocate the shape functions
    call paw_shape_exit(species%shape)

    ! Deallocate remaining arrays in paw_sp type
    deallocate(species%vhntzc_rad,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%vhntzc_rad',ierr)
    deallocate(species%rhoij0,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%rhoij0',ierr)
    deallocate(species%dij0,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%dij0',ierr)
    deallocate(species%tcore_den_rad,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%tcore_den_rad',ierr)
    deallocate(species%core_den_rad,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%core_den_rad',ierr)
    deallocate(species%tproj_rad,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%tproj_rad',ierr)
    deallocate(species%tphi_rad,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%tphi_rad',ierr)
    deallocate(species%phi_rad,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%phi_rad',ierr)
    do igrid=species%ngrid,1,-1
       deallocate(species%grid(igrid)%rab,stat=ierr)
       call utils_dealloc_check('paw_species_exit', &
            'grids(igrid)%rab',ierr)
       deallocate(species%grid(igrid)%r,stat=ierr)
       call utils_dealloc_check('paw_species_exit', &
            'grids(igrid)%r',ierr)
    end do
    deallocate(species%grid,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%grid',ierr)
    deallocate(species%m_pw_tot,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%m_pw_tot',ierr)
    deallocate(species%l_pw_tot,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%l_pw_tot',ierr)
    deallocate(species%ipw_tot,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%ipw_tot',ierr)
    deallocate(species%l_pw,stat=ierr)
    call utils_dealloc_check('paw_species_exit','species%l_pw',ierr)

    ! lr408
    if (species%core_wvfns_exist) then

       deallocate(species%core_wf_occ,stat=ierr)
       call utils_dealloc_check('paw_read_core_wfs','species%core_wf_occ',ierr)
       deallocate(species%core_wf_eig,stat=ierr)
       call utils_dealloc_check('paw_read_core_wfs','species%core_wf_eig',ierr)
       deallocate(species%core_wf_rad,stat=ierr)
       call utils_dealloc_check('paw_read_core_wfs','species%core_wf_rad',ierr)

       ! Deallocate Core Wavefunction grid(s)
       do igrid=species%ngrid_core,1,-1
          deallocate(species%grid_core(igrid)%rab,stat=ierr)
          call utils_dealloc_check('paw_species_exit', &
               'grids(igrid)%rab',ierr)
          deallocate(species%grid_core(igrid)%r,stat=ierr)
          call utils_dealloc_check('paw_species_exit', &
               'grids(igrid)%r',ierr)
       end do
       deallocate(species%grid_core,stat=ierr)
       call utils_dealloc_check('paw_read_core_wfs','species%grid_core',ierr)

       deallocate(species%n_core_wf,stat=ierr)
       call utils_dealloc_check('paw_read_core_wfs','species%n_core_wf',ierr)
       deallocate(species%icore_wf_tot,stat=ierr)
       call utils_dealloc_check('paw_read_core_wfs','species%icore_wf_tot',ierr)
       deallocate(species%m_core_wf_tot,stat=ierr)
       call utils_dealloc_check('paw_read_core_wfs','species%m_core_wf_tot',ierr)
       deallocate(species%l_core_wf_tot,stat=ierr)
       call utils_dealloc_check('paw_read_core_wfs','species%l_core_wf_tot',ierr)
       deallocate(species%l_core_wf,stat=ierr)
       call utils_dealloc_check('paw_read_core_wfs','species%l_core_wf',ierr)

       deallocate(species%core_aug_nLij,stat=ierr)
       call utils_dealloc_check('paw_dataset_init','species%core_aug_nLij',ierr)

    end if

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_species_exit'

  end subroutine paw_species_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! SUBROUTINES FOR INTERFACING WITH THE ATOM SOLVER !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_locpot_rad(locpot,npts,rad,isp,paw_sp)

    !=====================================================================!
    ! This subroutine fetches the local pseudopotential on a regular grid !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                         !
    !=====================================================================!

    use services, only: services_sbessj
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    integer, intent(in) :: npts
    real(kind=DP), intent(out) :: locpot(npts)
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp

    ! Local variables
    real(kind=DP), allocatable :: work(:)
    real(kind=DP) :: q,dq,Z
    integer :: npts_q
    integer :: ir,iq
    integer :: ierr

    npts_q = paw_sp(isp)%n_recip_pts
    allocate(work(npts_q),stat=ierr)
    call utils_alloc_check('paw_get_locpot_rad','work',ierr)

    dq = paw_sp(isp)%g_max / real(npts_q-1,kind=DP)
    Z = paw_sp(isp)%ion_charge

    do ir=1,npts
       do iq=1,npts_q
          q = real(iq-1,kind=DP)*dq
          work(iq) = (paw_sp(isp)%vhntzc_recip(iq)*q**2 - Z) * &
               services_sbessj(0,q*rad(ir))
       end do

       locpot(ir) = work(1)+work(npts_q)
       do iq = 2,npts_q-1,2
          locpot(ir) = locpot(ir) + 4.0_dp*work(iq)+2.0_dp*work(iq+1)
       enddo
       locpot(ir) = locpot(ir)*dq/3.0_dp
    end do

    locpot(:) = locpot(:)*(2.0_DP/PI)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_get_locpot_rad','work',ierr)

  end subroutine paw_get_locpot_rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_core_den_rad(core_den,npts,rad,isp,paw_sp)

    !============================================================!
    ! This subroutine fetches the core density on a regular grid !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    use services, only: services_sbessj
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    integer, intent(in) :: npts
    real(kind=DP), intent(out) :: core_den(npts)
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp

    ! Local variables
    real(kind=DP), allocatable :: work(:)
    real(kind=DP) :: q,dq
    integer :: npts_q
    integer :: ir,iq
    integer :: ierr

    if (.not.paw_sp(isp)%tcore_charge) then
       core_den(:) = 0.0_DP
       return
    end if

    npts_q = paw_sp(isp)%n_recip_pts
    allocate(work(npts_q),stat=ierr)
    call utils_alloc_check('paw_get_core_den_rad','work',ierr)

    dq = paw_sp(isp)%g_max / real(npts_q-1,kind=DP)

    do ir=1,npts
       do iq=1,npts_q
          q = real(iq-1,kind=DP)*dq
          work(iq) = paw_sp(isp)%tcore_den_recip(iq)*q**2 * &
                services_sbessj(0,q*rad(ir))
       end do

       core_den(ir) = work(1) + work(npts_q)
       do iq=2,npts_q-1,2
          core_den(ir) = core_den(ir) + 4.0_DP*work(iq) + 2.0_DP*work(iq+1)
       end do
       core_den(ir) = core_den(ir)*dq/3.0_DP
    end do

    core_den(:) = core_den(:)*(2.0_DP/PI)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_get_core_den_rad','work',ierr)

  end subroutine paw_get_core_den_rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_projector_info(isp,paw_sp,npw,npwtot,lmax,log_npts_max)

    !============================================================!
    ! This function fetches the number of shells of projectors.  !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    integer,intent(in) :: isp
    integer,intent(out) :: npw
    integer,intent(out) :: npwtot
    integer,intent(out) :: lmax
    integer,intent(out) :: log_npts_max

    ! Local Variables
    integer :: igrid

    npw = paw_sp(isp)%npw
    npwtot = paw_sp(isp)%npw_tot
    lmax = paw_sp(isp)%lmax

    log_npts_max = 0
    igrid = paw_sp(isp)%phi_grid
    log_npts_max = max(log_npts_max,paw_sp(isp)%grid(igrid)%npt)
    igrid = paw_sp(isp)%shape_grid
    log_npts_max = max(log_npts_max,paw_sp(isp)%grid(igrid)%npt)
    igrid = paw_sp(isp)%core_den_grid
    log_npts_max = max(log_npts_max,paw_sp(isp)%grid(igrid)%npt)

  end subroutine paw_get_projector_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_projectors_q(proj_q,dij0,ang_mom,npw,nsws, &
       nsws_max,lmax,qb,isp,paw_sp)

    !============================================================!
    ! This subroutine fetches the projectors at chosen q-points  !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    use services, only: services_1d_interpolation
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    integer, intent(in) :: npw
    integer, intent(in) :: lmax
    integer, intent(in) :: nsws(0:lmax)
    integer, intent(in) :: nsws_max
    integer, intent(in) :: isp
    integer, intent(out) :: ang_mom(npw)
    real(kind=DP), intent(in) :: qb(nsws_max,0:lmax)
    real(kind=DP), intent(out) :: proj_q(nsws_max,npw)
    real(kind=DP), intent(out) :: dij0(npw,npw)

    ! Local variables
    integer :: ipw,jpw,isw
    integer :: ipwtot,jpwtot
    integer :: li,lj,di
    real(kind=DP) :: qq

    ! Check matrices are right size
    if (npw/=paw_sp(isp)%npw) then
       call utils_abort('Error in paw_get_projectors_q: Wrong number of shells &
            &npw')
    end if

    ! Set the dij0 terms
    dij0(:,:) = 0.0_DP
    ipwtot = 1
    do ipw=1,npw
       li = paw_sp(isp)%l_pw(ipw)
       jpwtot = 1
       do jpw=1,npw
          lj = paw_sp(isp)%l_pw(jpw)
          if (li==lj) then
             do di=0,0
                dij0(ipw,jpw) = dij0(ipw,jpw) + &
                     paw_sp(isp)%dij0(ipwtot+di,jpwtot+di)
             end do
          end if
          jpwtot = jpwtot + 2*lj + 1
       end do
       ipwtot = ipwtot + 2*li + 1
    end do

    ! Interpolate the reciprocal-space projectors at the qb values
    do ipw=1,npw
       proj_q(:,ipw) = 0.0_DP
       ang_mom(ipw) = paw_sp(isp)%l_pw(ipw)
       do isw=1,nsws(ang_mom(ipw))
          qq = qb(isw,paw_sp(isp)%l_pw(ipw))
          proj_q(isw,ipw) = services_1d_interpolation( &
               paw_sp(isp)%tproj_recip(:,ipw), &
               paw_sp(isp)%n_recip_pts, &
               qq*paw_sp(isp)%inv_g_spacing, &
               paw_sp(isp)%l_pw(ipw))
       end do
    end do

  end subroutine paw_get_projectors_q


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_aug_funcs(qijL,shape_l,grad_shape_l,rhoij0,lmax, &
       npw,npts,rad,isp,paw_sp)

    !============================================================!
    ! This subroutine fetches the shape functions and the qijL   !
    ! terms for each angular momentum channel.                   !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    use gaunt_coeff, only: realgaunt
    use paw_shape, only: paw_shape_calculate
    use services, only: services_1d_interpolation
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    integer, intent(in) :: npts
    integer, intent(in) :: npw
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp
    integer, intent(inout) :: lmax
    real(kind=DP), intent(out) :: qijL(npw,npw,0:lmax)
    real(kind=DP), intent(out) :: shape_l(npts,0:lmax)
    real(kind=DP), intent(out) :: grad_shape_l(npts,0:lmax)
    real(kind=DP), intent(out) :: rhoij0(npw,npw)

    ! Local variables
    integer :: lup, mup
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,lj,mi,mj

    ! Ensure the guessed lmax value is the same PAW species' version
    if (lmax < paw_sp(isp)%lmax) then
       call utils_abort('Error paw_get_aug_funcs: value of lmax provided is &
            &too small')
    else if (lmax > paw_sp(isp)%lmax) then
       lmax = paw_sp(isp)%lmax
    end if

    ! Add up qijL factors to reduce npwtot x npwtot down to npw x npw
    qijL(:,:,:) = 0.0_DP
    do lup=0,lmax
       do mup=0,0
          do ipwtot=1,paw_sp(isp)%npw_tot
             li = paw_sp(isp)%l_pw_tot(ipwtot)
             mi = paw_sp(isp)%m_pw_tot(ipwtot)
             if (mi/=0) cycle
             ipw = paw_sp(isp)%ipw_tot(ipwtot)
             do jpwtot=1,paw_sp(isp)%npw_tot
                lj = paw_sp(isp)%l_pw_tot(jpwtot)
                mj = paw_sp(isp)%m_pw_tot(jpwtot)
                if (mj/=0) cycle
                jpw = paw_sp(isp)%ipw_tot(jpwtot)
                qijL(ipw,jpw,lup) = qijL(ipw,jpw,lup) + &
                     paw_sp(isp)%aug_nLij(ipw,jpw,lup) * &
                     realgaunt(lup,mup,li,mi,lj,mj) * sqrt_4pi
             end do
          end do
       end do
    end do

    ! Evaluate the shape function on the radial grid
    do lup=0,lmax
       shape_l(:,lup) = 0.0_DP
       grad_shape_l(:,lup) = 0.0_DP
       call paw_shape_calculate(shape_l(:,lup),grad_shape_l(:,lup), &
            rad,npts,lup,paw_sp(isp)%shape%shape_type, &
            paw_sp(isp)%shape%rshape, &
            paw_sp(isp)%shape%shape_alpha(:,lup), &
            paw_sp(isp)%shape%shape_q(:,lup), &
            paw_sp(isp)%shape%shape_lambda(lup), &
            paw_sp(isp)%shape%shape_sigma(lup))
       shape_l(:,lup) = shape_l(:,lup) / paw_sp(isp)%shape%norm(lup)
       grad_shape_l(:,lup) = grad_shape_l(:,lup) / paw_sp(isp)%shape%norm(lup)
    end do

    ! Add up rhoij for ipwtot,jpwtot contributing to each ipw,jpw to get
    ! initial guess projector density matrix
    rhoij0(:,:) = 0.0_DP
    do ipwtot=1,paw_sp(isp)%npw_tot
       ipw = paw_sp(isp)%ipw_tot(ipwtot)
       do jpwtot=1,paw_sp(isp)%npw_tot
          jpw = paw_sp(isp)%ipw_tot(jpwtot)

          rhoij0(ipw,jpw) = rhoij0(ipw,jpw) + paw_sp(isp)%rhoij0(ipwtot,jpwtot)
       end do
    end do

  end subroutine paw_get_aug_funcs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_tcore_hartree_on_grid(tcore_hartree_potential,&
        tcore_hartree_gzero_term,struct_fac,struct_fac_classical,paw_sp,grid,par)

    !===================================================================!
    ! This subroutine generates the Hartree potential of the pseudized  !
    ! core charge \tilde{n_Zc} on the simulation cell fine grid. This   !
    ! is the PAW equivalent of the local potential.                     !
    !-------------------------------------------------------------------!
    ! Arguments:                                                        !
    !  tcore_hartree_potential (out) : V_H [n_Zc] (r) in real space     !
    !  tcore_hartree_gzero_term (out)    : V_H [n_Zc] (G=0)             !
    !  struct_fac (in) : The structure factor for each species in       !
    !  reciprocal space.                                                !
    !  struct_fac_classical (in) : The structure factor for classical   !
    !  atoms in reciprocal space.                                       !
    !  grid (in) : The grid definition.                                 !
    !-------------------------------------------------------------------!
    ! Adapted on 20/05/2010.by Nicholas Hine from the routine           !
    ! pseudopotentials_local_on_fine, originally written by             !
    ! Chris-Kriton Skylaris on 21/2/2004 with subsequent modifications  !
    ! by Peter Haynes, Arash Mostofi and Nicholas Hine.                 !
    !===================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_bcast, pub_my_proc_id, pub_root_proc_id
    use classical_pot, only: classical_pot_recip
    use fourier, only: fourier_apply_cell_backward
    use parallel_strategy, only: PARAL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)  :: grid
    type(PARAL_INFO), intent(in)  :: par
    complex(kind=DP), intent(in) :: struct_fac(par%num_pspecies, &
         grid%ld3, grid%ld2, grid%max_slabs23)
    complex(kind=DP), intent(in) :: struct_fac_classical( &
         grid%ld3, grid%ld2, grid%max_slabs23)
    real(kind=DP), intent(out) :: tcore_hartree_potential(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    real(kind=DP), intent(out) :: tcore_hartree_gzero_term
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)

    ! Local Variables
    integer :: ierr              ! Error flag
    complex(kind=DP), allocatable :: tcore_hartree_recip(:,:,:)

    call timer_clock('paw_tcore_hartree_on_grid',1)

    ! ndmh: expand to 3D in reciprocal fine grid and sum together
    !       (along with the structure factor) the pseudo-core Hartree
    !       potential for each species to obtain the total pseudo-core Hartree
    !       potential in reciprocal representation on the fine grid
    allocate(tcore_hartree_recip(grid%ld3,grid%ld2, &
         grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_tcore_hartree_on_grid','tcore_hartree_recip', &
         ierr)

    call paw_tcore_hartree_rec(tcore_hartree_recip, &      ! output
         struct_fac,paw_sp,grid,par)                       ! input

    ! store g=0 term seperately for charged system correction
    tcore_hartree_gzero_term = 0.0_DP
    if (pub_my_proc_id==grid%proc_slab23(1)) then
       tcore_hartree_gzero_term = real(tcore_hartree_recip(1,1,1) &
            / (grid%n1 * grid%n2 * grid%n3), kind=DP)
       ! jd: Wrapped to work around warning on complex to real assignment
    end if
    call comms_bcast(pub_root_proc_id,tcore_hartree_gzero_term)

    ! ndmh: include external potential from "classical" atoms
    if (par%nat_classical > 0) then
       call classical_pot_recip(tcore_hartree_recip, &       ! output
            struct_fac_classical,grid)                       ! input
    endif

    ! FFT the local ionic potential from reciprocal to real space
    call fourier_apply_cell_backward(tcore_hartree_potential, &
         tcore_hartree_recip,grid)

    deallocate(tcore_hartree_recip,stat=ierr)
    call utils_dealloc_check('paw_tcore_hartree_on_grid', &
         'tcore_hartree_recip',ierr)

    call timer_clock('paw_tcore_hartree_on_grid',2)

  end subroutine paw_tcore_hartree_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_tcore_hartree_rec(fine_complex,struct_fac,paw_sp,grid,par)

    !=================================================================!
    ! This subroutine generates in reciprocal space the Hartree       !
    ! potential in the simulation cell due to the pseudized core      !
    ! density \tilde{n_Zc} of all ions.                               !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !  fine_complex (out) : the Hartree potential in reciprocal space !
    !  struct_fac (in) : The structure factor for each species in     !
    !  reciprocal space.                                              !
    !-----------------------------------------------------------------!
    ! Adapted on 20/05/2010.by Nicholas Hine from the routine         !
    ! pseudopotentials_local_on_fine, originally written by           !
    ! Chris-Kriton Skylaris in 2000 with subsequent modifications     !
    ! by Peter Haynes, Arash Mostofi and Nicholas Hine.               !
    !=================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
!$  use rundat, only: pub_threads_max
    use services, only: services_1d_interpolation
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in) :: par
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    type(GRID_INFO), intent(in) :: grid
    complex(kind=DP), intent(in) :: struct_fac(par%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    complex(kind=DP), intent(out) :: fine_complex(grid%ld3,&
         grid%ld2, grid%max_slabs23)

    ! Local variables
    integer :: species              ! Atomic species counter
    integer :: ipt,i3,i2,islab23    ! Reciprocal grid loop counters
    real(kind=DP) :: gvec(3)        ! G vector
    real(kind=DP) :: g_length       ! Length of this G vector
    real(kind=DP) :: v_loc_value    ! Local potential at this G
    !real(kind=DP) :: g_cut, alpha   ! For filtering locps
    real(kind=DP),parameter :: fourpi = 4.0_DP * PI ! Constant multiplier

    ! Loop over reciprocal space grid on this proc
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i2,i3,islab23,gvec,g_length,species,v_loc_value) &
!$OMP SHARED (pub_my_proc_id,grid,paw_sp,fine_complex,par, &
!$OMP      struct_fac,pub_threads_max)
    do ipt=1,grid%num_slabs23*grid%n2*grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1

       ! Get magnitude of this G-vector
       call cell_grid_recip_pt(gvec,islab23 + &
            grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)
       g_length = sqrt(sum(gvec(1:3)**2))

       ! ndmh: initialise (cache-efficiently)
       fine_complex(i3,i2,islab23) = cmplx(0.0_DP,0.0_DP,kind=DP)

       ! Loop over atomic species
       do species=1,par%num_pspecies

          ! Get potential at this G-vector
          v_loc_value = services_1d_interpolation( &
               paw_sp(species)%vhntzc_recip, &
               paw_sp(species)%n_recip_pts,&
               g_length*paw_sp(species)%inv_g_spacing,0)

          ! Add back Coulomb potential
          v_loc_value = v_loc_value - paw_sp(species)%ion_charge * &
                grid%coulomb_recip(i3,i2,islab23)

          ! ndmh: scale by 4 pi/weight
          v_loc_value = v_loc_value * fourpi / grid%weight

          fine_complex(i3,i2,islab23) = fine_complex(i3,i2,islab23) + &
               struct_fac(species,i3,i2,islab23) * v_loc_value

       end do    ! loop over species

    end do
!$OMP END PARALLEL DO

    ! G=0 element must be real
    if (pub_my_proc_id==grid%proc_slab23(1)) then
       call utils_assert(aimag(fine_complex(1,1,1)) == 0.0_DP, &
            'Error in paw_tcore_hartree_rec: potential not real')
    end if

    if (grid%num_slabs23 > 0) then
       ! Nyquist filter (fine grid is always going to be even)
       fine_complex(grid%n3/2+1,:,:) = (0.0_DP,0.0_DP)
       fine_complex(:,grid%n2/2+1,:) = (0.0_DP,0.0_DP)
    end if
    ! Nyquist filter for last slab
    if (pub_my_proc_id==grid%proc_slab23(grid%n1/2+1))&
         fine_complex(:,:,grid%num_slabs23) = (0.0_DP,0.0_DP)

  end subroutine paw_tcore_hartree_rec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_tcore_density_recip(tcore_density,struct_fac,paw_sp,grid,par)

    !==================================================================!
    ! This subroutine reconstructs the core density for the whole      !
    ! supercell on the fine real space grid using the information      !
    ! stored for each dataset in the array tcore_den_recip.            !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  tcore_density (out) : data-parallelised core density            !
    !  struct_fac  (in)    : data-parallelised structure factor        !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 28/05/10.           !
    !==================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_proc_id
    use fourier, only: fourier_apply_cell_backward
!$  use rundat, only: pub_threads_max
    use parallel_strategy, only: PARAL_INFO
    use services, only: services_1d_interpolation
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in) :: par
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out) :: tcore_density(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    complex(kind=DP), intent(in) :: struct_fac(par%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)

    ! Local variables
    integer :: ierr              ! Error flag
    complex(kind=DP), allocatable :: tcore_density_recip(:,:,:)
    real(kind=DP) :: gvec(3),g_length, tcore_den_value, factor
    integer :: species              ! Atomic species counter
    integer :: ipt,i3,i2,islab23    ! Reciprocal grid loop counters

    call timer_clock('paw_tcore_density_recip',1)

    ! ndmh: allocate storage for core density in reciprocal space
    allocate(tcore_density_recip(grid%ld3, &
         grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_tcore_density', &
         'tcore_density_recip',ierr)

    ! ndmh: loop over reciprocal space grid on this proc
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i2,i3,islab23,gvec,g_length,species,tcore_den_value) &
!$OMP SHARED (pub_my_proc_id,grid,paw_sp,tcore_density_recip,par, &
!$OMP      struct_fac,pub_threads_max)
    do ipt=1,grid%num_slabs23*grid%n2*grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1

       ! Initialise
       tcore_density_recip(i3,i2,islab23) = (0.0_DP,0.0_DP)

       ! Loop over atomic species
       do species=1,par%num_pspecies

          ! Check if we have a core charge for this species
          if (.not.paw_sp(species)%tcore_charge) cycle

          ! Get magnitude of this G-vector
          call cell_grid_recip_pt(gvec,islab23 + &
               grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)
          g_length = sqrt(sum(gvec(1:3)**2))

          ! Get core density at this G-vector
          tcore_den_value = services_1d_interpolation( &
               paw_sp(species)%tcore_den_recip, &
               paw_sp(species)%n_recip_pts,&
               g_length*paw_sp(species)%inv_g_spacing,0)

          tcore_density_recip(i3,i2,islab23) = &
               tcore_density_recip(i3,i2,islab23) + &
               struct_fac(species,i3,i2,islab23) * tcore_den_value

       end do    ! loop over species

    end do
!$OMP END PARALLEL DO

    ! FFT the core density from reciprocal to real space
    call fourier_apply_cell_backward(tcore_density,tcore_density_recip,grid)

    ! ndmh: scale with 1.0/weight
    factor = 1.0_DP / grid%weight
    tcore_density = factor * tcore_density

    ! ndmh: deallocate storage for core density in reciprocal space
    deallocate(tcore_density_recip,stat=ierr)
    call utils_dealloc_check('paw_tcore_density', &
         'tcore_density_recip',ierr)

    call timer_clock('paw_tcore_density_recip',2)

  end subroutine paw_tcore_density_recip




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !==================================================================!
  ! This subroutine initialises the box on which the core density is !
  ! calculated (this is used in paw_tcore_density and                !
  ! paw_nlcc_calculate_forces)                                       !
  !------------------------------------------------------------------!
  ! This subroutine was written by Nicholas Hine in July 2019        !
  !==================================================================!

  subroutine paw_make_core_den_box(core_den_box,grid,cell,paw_sp,par)

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: stdout,DP
    use geometry, only: geometry_magnitude
    use fft_box, only: FFTBOX_INFO, fftbox_init, fftbox_find_size, fftbox_exit
    use fourier, only: fourier_init_box, fourier_exit_box
    use parallel_strategy, only: PARAL_INFO
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_banner

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(GRID_INFO), intent(in) :: grid
    type(PARAL_INFO), intent(in)  :: par
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    type(FFTBOX_INFO), intent(inout) :: core_den_box

    ! Local Variables
    integer :: isp, igrid, npts
    integer :: box_n1,box_n2,box_n3
    real(kind=DP) :: max_rad
    logical :: extend(3)
    integer :: pref(3)
    real(kind=DP) :: halo

    ! Find radius of biggest core density over all species
    max_rad = 0.0_DP
    do isp=1,par%num_pspecies
       igrid = paw_sp(isp)%core_den_grid
       npts = paw_sp(isp)%grid(igrid)%npt
       max_rad = max(max_rad,paw_sp(isp)%grid(igrid)%r(npts))
    end do

    ! Find size of core density box
    extend = (/.false.,.false.,.false./)
    pref = (/0,0,0/)
    halo = 0.0_DP
    call fftbox_find_size(box_n1,box_n2,box_n3, &
       max_rad,1,extend,0.0_DP,pref, &
       cell%a1,cell%a2,cell%a3,cell%b1,cell%b2,cell%b3, &
       geometry_magnitude(grid%da1),geometry_magnitude(grid%da2), &
       geometry_magnitude(grid%da3),grid%n1,grid%n2,grid%n3)

    ! Report size of core den box
    if (pub_on_root) then
       write(stdout,'(a)') ''
       write(stdout,'(a)') utils_banner('=', 'Core Density')
       write(stdout,'(3(a,i4))') '               Core den box size: ',&
            box_n1,' x',box_n2,' x',box_n3
       write(stdout,'(a)') repeat('=',80)
    end if

    ! Initialise mdl%fftbox structure
    call fftbox_init(core_den_box,box_n1,box_n2,box_n3,cell)
    call fourier_init_box(core_den_box)

  end subroutine paw_make_core_den_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !==================================================================!
  ! This subroutine cleans up structures associated with the core    !
  ! density box.                                                     !
  !------------------------------------------------------------------!
  ! This subroutine was written by Nicholas Hine in July 2019        !
  !==================================================================!

  subroutine paw_exit_core_den_box(core_den_box)

    use comms, only: pub_on_root
    use constants, only: stdout,DP
    use fft_box, only: FFTBOX_INFO, fftbox_exit
    use fourier, only: fourier_exit_box

    implicit none

    ! Arguments
    type(FFTBOX_INFO), intent(inout) :: core_den_box

    call fourier_exit_box(core_den_box)
    call fftbox_exit(core_den_box)

  end subroutine paw_exit_core_den_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !==================================================================!
  ! This subroutine replaces the old paw_tcore_density subroutine. It!
  ! constructs the core density for the whole supercell on the fine  !
  ! real space grid using the information stored for each dataset in !
  ! the array tcore_den_recip. The approach is derived from          !
  ! augmentation_density_on_grid                                     !
  !------------------------------------------------------------------!
  ! Arguments:                                                       !
  !  atom_tcore (out) : core density on real-space augmentation box  !
  !------------------------------------------------------------------!
  ! This subroutine was written by Gabriel Constantinescu - 20/07/15 !
  ! Updates the old paw_tcore_density, written by Nicholas Hine      !
  ! Updated by Nicholas Hine on 20/06/19 to use a new box, whose     !
  ! size is set by the actual maximum core density radius, to        !
  ! perform all box FFTs for generation of core density.             !
  !==================================================================!

  subroutine paw_tcore_density(tcore_den,grid,cell,paw_sp,par, &
       direction1,weight1,direction2,weight2)

    use cell_grid, only: GRID_INFO, cell_grid_deposit_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: pub_my_proc_id, pub_on_root
    use constants, only: VERBOSE, stdout,DP
    use geometry, only: geometry_magnitude
    use fft_box, only: FFTBOX_INFO, fftbox_init, fftbox_find_size, fftbox_exit
    use fourier, only: fourier_init_box, fourier_exit_box
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_paw, pub_debug_on_root, &
         pub_threads_num_fftboxes
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_banner, utils_alloc_check, &
         utils_dealloc_check

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: tcore_den(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    type(PARAL_INFO), intent(in)  :: par
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    integer, intent(in), optional :: direction1(par%nat)
    real(kind=DP), intent(in), optional :: weight1(par%nat)
    integer, intent(in), optional :: direction2(par%nat)
    real(kind=DP), intent(in), optional :: weight2(par%nat)

    ! Local Variables
    integer :: iat, loc_iat
    integer :: ierr
    integer :: isp
    integer :: igrid, npts
    integer :: box_n1,box_n2,box_n3
    integer :: box_start1,box_start2,box_start3
    real(kind=DP), allocatable :: atom_tcore(:,:,:)
    real(kind=DP), allocatable :: buffer(:,:,:)
    real(kind=DP) :: max_rad
    type(FFTBOX_INFO) :: core_den_box
    complex(kind=DP), allocatable :: atom_tcore_recip(:,:,:)
    logical :: i_have_box
    logical :: extend(3)
    integer :: pref(3)
    real(kind=DP) :: halo

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &paw_tcore_density'

    ! Start Timer
    call timer_clock('paw_tcore_density',1)

    tcore_den(:,:,:) = 0.0_DP

    ! Set up structures for core density box
    call paw_make_core_den_box(core_den_box,grid,cell,paw_sp,par)
    box_n1 = core_den_box%total_pt1
    box_n2 = core_den_box%total_pt2
    box_n3 = core_den_box%total_pt3

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(atom_tcore,buffer,atom_tcore_recip,loc_iat,iat,isp, &
!$OMP      ierr,box_start1,box_start2,box_start3,i_have_box) &
!$OMP SHARED(pub_my_proc_id,grid,cell,par,core_den_box,paw_sp, &
!$OMP      box_n1,box_n2,box_n3,pub_threads_num_fftboxes,pub_paw, &
!$OMP      tcore_den,direction1,weight1,direction2,weight2)

    ! Allocate workspace
    allocate(atom_tcore(box_n1,box_n2,box_n3),stat=ierr)
    call utils_alloc_check('paw_tcore_density','atom_tcore',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('paw_tcore_density','buffer',ierr)
    allocate(atom_tcore_recip(box_n1,box_n2,box_n3),stat=ierr)
    call utils_alloc_check('paw_tcore_density','atom_tcore_recip',ierr)


!$OMP DO
    do loc_iat=1,par%max_atoms_on_proc

       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id)) then
          iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          isp = par%elements_on_proc(loc_iat)%pspecies_number

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               par%elements_on_proc(loc_iat)%centre, box_n1, box_n2, box_n3, &
               grid,cell)

          ! Call appropriate routine to generate the tcore density for
          ! this atom in the augmentation box
          if (paw_sp(isp)%tcore_charge) then
             if (.not.present(direction1)) then
                call paw_atom_tcore_den(par%elements_on_proc(loc_iat)%centre, &
                     grid,cell,paw_sp,core_den_box,box_start1, &
                     box_start2,box_start3,atom_tcore_recip, &
                     atom_tcore,isp,par%num_pspecies)
                i_have_box = .true.
             else
                if (present(direction2)) then
                   call paw_atom_tcore_den( &
                        par%elements_on_proc(loc_iat)%centre, &
                        grid,cell,paw_sp,core_den_box,box_start1, &
                        box_start2,box_start3,atom_tcore_recip, &
                        atom_tcore,isp,par%num_pspecies,&
                        direction1(par%orig_atom(iat)), &
                        direction2(par%orig_atom(iat)))

                    atom_tcore = atom_tcore * weight1(par%orig_atom(iat)) * &
                         weight2(par%orig_atom(iat))
                else
                   call paw_atom_tcore_den( &
                        par%elements_on_proc(loc_iat)%centre, &
                        grid,cell,paw_sp,core_den_box,box_start1, &
                        box_start2,box_start3,atom_tcore_recip, &
                        atom_tcore,isp,par%num_pspecies,&
                        direction1(par%orig_atom(iat)))

                   atom_tcore = atom_tcore * weight1(par%orig_atom(iat))
                end if
                i_have_box = .true.
             end if
          else
              atom_tcore(:,:,:) = 0.0_DP
              i_have_box = .false.
          end if

       else
          ! Nothing to deposit on this proc
          i_have_box = .false.
       end if

       ! Deposit this box to the simulation cell if present, or just wait for
       ! data from other procs if no box
!$OMP CRITICAL
       call cell_grid_deposit_box(tcore_den(:,:,:), &
            atom_tcore(:,:,:), buffer, grid, &
            box_n1, box_n2, box_n3, box_n1, box_n2, &
            box_start1, box_start2, box_start3, i_have_box, .false.)
!$OMP END CRITICAL
    end do
!$OMP END DO

    ! Deallocate temporary arrays and matrices
    deallocate(atom_tcore_recip,stat=ierr)
    call utils_dealloc_check('paw_tcore_density', &
         'atom_tcore_recip',ierr)
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('paw_tcore_density','buffer',ierr)
    deallocate(atom_tcore,stat=ierr)
    call utils_dealloc_check('paw_tcore_density','atom_tcore',ierr)

!$OMP END PARALLEL

    call paw_exit_core_den_box(core_den_box)

    ! Stop Timer
    call timer_clock('paw_tcore_density',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_tcore_density'

  end subroutine paw_tcore_density


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !==================================================================!
  ! This subroutine outputs the pseudo-core density for a given atom !
  ! of species isp, on the real-space augmentation box               !
  !------------------------------------------------------------------!
  ! Arguments:                                                       !
  !  tcore_den (out) : data-parallelised core density                !
  !------------------------------------------------------------------!
  ! This subroutine was written by Gabriel Constantinescu - 20/07/15 !
  !==================================================================!

  subroutine paw_atom_tcore_den(atom_centre,grid,cell,paw_sp,box, &
       box_start1,box_start2,box_start3, &
       atom_tcore_recip, atom_tcore,isp, num_pspecies, cart1, cart2)

    use constants, only: DP
    use cell_grid, only: GRID_INFO
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use gaunt_coeff, only: realgaunt
    use geometry, only: POINT
    use simulation_cell, only: CELL_INFO

    ! Arguments
    type(FFTBOX_INFO), intent(in) :: box
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(POINT), intent(in) :: atom_centre
    integer, intent(in) :: isp
    integer, intent(in) :: box_start1,box_start2,box_start3
    complex(kind=DP), intent(out) :: atom_tcore_recip(box%total_pt1, &
         box%total_pt2,box%total_pt3)
    real(kind=DP), intent(out) :: atom_tcore(box%total_pt1, &
         box%total_pt2,box%total_pt3)
    integer, intent(in) :: num_pspecies
    integer, intent(in), optional :: cart1
    integer, intent(in), optional :: cart2

    atom_tcore_recip = 0.0_DP

    ! Get tcore on augmentation box in reciprocal space

    if (.not.present(cart1)) then
       call paw_tcore_on_box_recip(atom_tcore_recip, box_start1,box_start2, &
            box_start3,grid,cell,box%total_pt1,box%total_pt2,box%total_pt3, &
            atom_centre,isp,paw_sp,num_pspecies)
    end if

    if (present(cart1).and.(.not.present(cart2))) then
       call paw_tcore_on_box_recip(atom_tcore_recip, box_start1,box_start2, &
            box_start3,grid,cell,box%total_pt1,box%total_pt2,box%total_pt3, &
            atom_centre,isp,paw_sp,num_pspecies,cart1)
    end if

    if (present(cart1).and.(present(cart2))) then
       call paw_tcore_on_box_recip(atom_tcore_recip, box_start1,box_start2, &
            box_start3,grid,cell,box%total_pt1,box%total_pt2,box%total_pt3, &
            atom_centre,isp,paw_sp,num_pspecies,cart1,cart2)
    end if




    ! Fourier transform to real space in augmentation box
    call fourier_apply_box('F','B',atom_tcore_recip(:,:,:),box=box)
    atom_tcore(:,:,:) = real(atom_tcore_recip(:,:,:),kind=DP) / grid%weight

  end subroutine paw_atom_tcore_den


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !==================================================================!
  ! This subroutine outputs the FFT of the pseudo-core density for a !
  ! given atom of species isp, on the augmentation box               !
  !------------------------------------------------------------------!
  ! Arguments:                                                       !
  !  tcore_den (out) : data-parallelised core density                !
  !------------------------------------------------------------------!
  ! This subroutine was written by Gabriel Constantinescu - 20/07/15 !
  !==================================================================!

  subroutine paw_tcore_on_box_recip(tcore_on_box_recip, cell_start1, &
      cell_start2,cell_start3,grid,cell,box_n1,box_n2,box_n3,atom_origin, &
      isp,paw_sp, num_pspecies, cart1, cart2)


    use basis, only: basis_box_origin_to_atom
    use cell_grid, only: GRID_INFO
    use constants, only: PI, cmplx_i
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_real_sph_harm, sw_bessel_accurate, sw_init
    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    complex(kind=DP),intent(out) :: tcore_on_box_recip(box_n1,box_n2, &
         box_n3)
    integer,intent(in) :: cell_start1
    integer,intent(in) :: cell_start2
    integer,intent(in) :: cell_start3
    type(CELL_INFO),intent(in) :: cell
    type(GRID_INFO),intent(in) :: grid
    type(POINT),intent(in) :: atom_origin
    integer, intent(in) :: num_pspecies
    type(PAW_SPECIES), intent(in) :: paw_sp(num_pspecies)
    integer, intent(in) :: isp
    integer, intent(in), optional :: cart1
    integer, intent(in), optional :: cart2

    ! Local Variables
    integer :: i1,i2,i3
    integer :: k1,k2,k3
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    real(kind=DP) :: gval
    real(kind=DP) :: slmval
    type(POINT) :: g_vector3,g_vector2,g_vector1
    type(POINT) :: pairvec3,pairvec2,pairvec1
    real(kind=DP) :: g_length
    complex(kind=DP) :: phase_fac

    call basis_box_origin_to_atom(box_to_atom,box_origin,atom_origin, &
         cell_start1,cell_start2,cell_start3,grid%da1,grid%da2,grid%da3, &
         cell)

    ! initialisations
    tcore_on_box_recip = 0.0_DP
    pairvec3=( real(grid%n3,DP)/real(box_n3,DP) ) * cell%b3
    pairvec2=( real(grid%n2,DP)/real(box_n2,DP) ) * cell%b2
    pairvec1=( real(grid%n1,DP)/real(box_n1,DP) ) * cell%b1
    phase_fac = 1.0_DP ! cmplx(0.0_DP,-1.0_DP,kind=DP)**0

    ! loop over reciprocal-space grid points of augmentation box
    do i3=1,box_n3
       k3 = i3 - 1
       if (k3>box_n3/2) k3 = k3 - box_n3
       g_vector3 = real(k3,kind=DP)*pairvec3

       do i2=1,box_n2
          k2 = i2 - 1
          if (k2>box_n2/2) k2 = k2 - box_n2
          g_vector2 = real(k2,kind=DP)*pairvec2 + g_vector3

          do i1=1,box_n1
             k1 = i1 - 1
             if (k1>box_n1/2) k1 = k1 - box_n1
             g_vector1 = real(k1,kind=DP)*pairvec1 + g_vector2

             g_length = geometry_magnitude(g_vector1)


            gval =  services_1d_interpolation( &
                paw_sp(isp)%tcore_den_recip, &
                paw_sp(isp)%n_recip_pts,&
                g_length*paw_sp(isp)%inv_g_spacing,0)

            ! multiply value of tcore_recip at each g-point by the
            ! appropriate phase factor (n.b. real part and imaginary part
            ! stored as separate consecutive elements in the x-direction
            ! of the array) and the appropriate real spherical harmonic
            !factor.

            slmval = sw_real_sph_harm(g_vector1%x,g_vector1%y, &
                     g_vector1%z, g_length, 0, 0)  * sqrt(4.0_DP*PI)

            tcore_on_box_recip(i1,i2,i3) = &
                 tcore_on_box_recip(i1,i2,i3) + &
                 gval * slmval * phase_fac * &
                 exp(-cmplx_i*(g_vector1.dot.box_to_atom))

            if (present(cart1)) then
               if (cart1 == 1) then
                  tcore_on_box_recip(i1,i2,i3) = tcore_on_box_recip(i1,i2,i3) *&
                       (-cmplx_i*g_vector1%x)
               else if (cart1 == 2) then
                  tcore_on_box_recip(i1,i2,i3) = tcore_on_box_recip(i1,i2,i3) *&
                       (-cmplx_i*g_vector1%y)
               else if (cart1 == 3) then
                  tcore_on_box_recip(i1,i2,i3) = tcore_on_box_recip(i1,i2,i3) *&
                       (-cmplx_i*g_vector1%z)
               end if
            end if

            if (present(cart2)) then
               if (cart2 == 1) then
                  tcore_on_box_recip(i1,i2,i3) = tcore_on_box_recip(i1,i2,i3) *&
                       (-cmplx_i*g_vector1%x)
               else if (cart2 == 2) then
                  tcore_on_box_recip(i1,i2,i3) = tcore_on_box_recip(i1,i2,i3) *&
                       (-cmplx_i*g_vector1%y)
               else if (cart2 == 3) then
                  tcore_on_box_recip(i1,i2,i3) = tcore_on_box_recip(i1,i2,i3) *&
                       (-cmplx_i*g_vector1%z)
               end if
            end if

          end do
       end do
    end do

  end subroutine paw_tcore_on_box_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_species_init_proj(proj_set,paw_sp,elements,par,is_cmplx)

    !=================================================================!
    ! This subroutine allocates and initialises all fftbox_proj_recip !
    ! for each paw_sp element. Each such fftbox_proj_recip is         !
    ! a full non-local projector in reciprocal space in the FFTbox.   !
    !-----------------------------------------------------------------!
    ! Adapted on 21/05/2010.by Nicholas Hine from the routine         !
    ! pseudopot_species_init_proj, originally written by              !
    ! Chris-Kriton Skylaris in 2004 with subsequent modifications     !
    ! by Peter Haynes, Arash Mostofi and Nicholas Hine.               !
    ! Modified to accept par as input by Robert Charlton, 29/08/2018. !
    !=================================================================!

    use comms, only: pub_on_root
    use geometry, only: OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use projectors, only: PROJECTOR_SET, projectors_allocate_set
    use rundat, only: pub_output_detail
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    type(PROJECTOR_SET), intent(inout) :: proj_set
    ! agrecocmplx
    logical, optional, intent(in) :: is_cmplx

    ! Local Variables
    integer :: iat, orig_iat
    integer :: isp
    integer :: shell
    integer :: proj_count
    ! agrecocmplx
    logical :: loc_cmplx

    ! agrecocmplx
    if (present(is_cmplx)) then
       loc_cmplx = is_cmplx
    else
       loc_cmplx = .false.
    end if

    if (pub_on_root  .and. pub_output_detail >= VERBOSE) &
         write(stdout,'(a)') '... PAW projector initialisation'

    ! ndmh: find number of unique projectors
    proj_count = 0
    do isp=1,par%num_pspecies
       proj_count = proj_count + paw_sp(isp)%npw_tot
    end do

    ! Set up the entries of proj_set
    proj_set%n_proj_species = par%num_pspecies
    ! agrecocmplx
    call projectors_allocate_set(proj_set, &
         maxval(paw_sp(:)%npw),maxval(paw_sp(:)%n_recip_pts),par=par, &
         is_cmplx=loc_cmplx)

    ! ndmh: set species_num_proj and species_first_proj values
    ! ndmh: also set gmax, n_rad_pts, n_shells, ang_mom and rad_proj_recip
    proj_count = 1
    do isp=1,proj_set%n_proj_species
       proj_set%species_num_proj(isp) = paw_sp(isp)%npw_tot
       proj_set%species_first_proj(isp) = proj_count
       proj_set%gmax(isp) = paw_sp(isp)%g_max
       proj_set%n_rad_pts(isp) = paw_sp(isp)%n_recip_pts
       proj_set%num_shells(isp) = paw_sp(isp)%npw
       proj_set%ang_mom(:,isp) = 0
       proj_set%rad_proj_recip(:,:,isp) = 0.0_DP
       do shell=1,paw_sp(isp)%npw
          proj_set%ang_mom(shell,isp) = paw_sp(isp)%l_pw(shell)
          proj_set%rad_proj_recip(1:paw_sp(isp)%n_recip_pts,shell,isp) = &
               paw_sp(isp)%tproj_recip(1:paw_sp(isp)%n_recip_pts,shell)
       end do
       proj_count = proj_count + proj_set%species_num_proj(isp)
    end do

    ! ndmh: copy projector centre and radius from elements array
    do iat=1,par%nat
       orig_iat = par%orig_atom(iat)
       proj_set%proj_centre(iat) = elements(orig_iat)%centre
       proj_set%proj_max_radius(iat) = elements(orig_iat)%max_core_radius
       proj_set%proj_species(iat) = elements(orig_iat)%pspecies_number
    end do

  end subroutine paw_species_init_proj


  subroutine paw_species_calc_proj_prec_mat(proj_set,cell,fftbox, &
       precond_func_recip,paw_sp)

    !==================================================================!
    ! This subroutine returns the matrix used in the PAW-correction to !
    ! the kinetic energy preconditioning:                              !
    !   O_prec = -O*(I+C*O)^-1                                         !
    ! with O - partial wave overlap difference                         !
    ! and C - <p^i|(1+T)^-1|p^j>                                       !
    !------------------------------------------------------------------!
    ! This subroutine was written by Gabriel Constantinescu in May'14. !
    !==================================================================!

    use constants, only: DP, stdout
    use fft_box, only: FFTBOX_INFO
    use linalg, only: linalg_mat_mul_serial, linalg_invert_serial
    use projectors, only: PROJECTOR_SET, projectors_proj_overlap
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    real(kind=DP), intent(in)  :: precond_func_recip(:,:,:)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    ! Local Variables
    integer :: i1,i2
    integer :: ierr
    integer :: isp
    integer :: max_val_species
    integer :: max_num_proj_tot ! max number of projectors for all the species
    real(kind=DP), allocatable :: proj_overlap(:,:,:) ! C^ij=<p^i|(1+T)^-1|p^j>
    real(kind=DP), allocatable :: O_mat(:,:,:) ! <phi_i|phi_j>-<tphi_i|tphi_j>
    real(kind=DP), allocatable :: prod_mat(:,:,:) ! C*O
    real(kind=DP), allocatable :: inv_mat(:,:,:) ! (I+CO)^-1 (I= identity)
    integer :: num_pspecies

    ! Start Timer
    call timer_clock('paw_species_calc_proj_prec_mat', 1)

    num_pspecies = size(paw_sp)

    ! gcc32: allocate proj_set%O_prec, i.e.  i.e. -O*(I+C*O)^-1
    if (.not.allocated(proj_set%O_prec)) then
       max_val_species = maxval(proj_set%species_num_proj(:))
       allocate(proj_set%O_prec(max_val_species, max_val_species, &
            num_pspecies),stat=ierr)
       call utils_alloc_check('paw_species_init_proj','proj_set%O_prec',ierr)
    end if

    ! Allocation
    max_num_proj_tot = maxval(proj_set%species_num_proj(:))
    allocate(proj_overlap(max_num_proj_tot, max_num_proj_tot, &
             num_pspecies),stat=ierr)
    call utils_alloc_check('paw_species_calc_proj_prec_mat','proj_overlap',ierr)
    allocate(O_mat(max_num_proj_tot, max_num_proj_tot, &
             num_pspecies),stat=ierr)
    call utils_alloc_check('paw_species_calc_proj_prec_mat','O_mat',ierr)
    allocate(prod_mat(max_num_proj_tot, max_num_proj_tot, &
             num_pspecies),stat=ierr)
    call utils_alloc_check('paw_species_calc_proj_prec_mat','prod_mat',ierr)
    allocate(inv_mat(max_num_proj_tot, max_num_proj_tot, &
             num_pspecies),stat=ierr)
    call utils_alloc_check('paw_species_calc_proj_prec_mat','inv_mat',ierr)

    ! Calculate O_ij(s)
    call paw_projector_overlap_basic(O_mat,paw_sp,num_pspecies)

    ! Calculate C^ij(s)
    call projectors_proj_overlap(proj_overlap,proj_set,precond_func_recip, &
         cell,fftbox)

    ! Perform C*O matrix-multiplication
    prod_mat = 0.0_DP
    do isp=1,num_pspecies
       call linalg_mat_mul_serial(prod_mat(:,:,isp),proj_overlap(:,:,isp), &
            O_mat(:,:,isp))
    end do

    ! gcc32: now, add identity matrix to P*O
    do isp=1,num_pspecies
       do i1 = 1,max_num_proj_tot
          prod_mat(i1,i1,isp) = prod_mat(i1,i1,isp) + 1.0_DP
       end do
    end do

    ! first, copy to prod_mat to inv_mat
    do isp=1,num_pspecies
       do i1 = 1,max_num_proj_tot
          do i2 = 1,max_num_proj_tot
             inv_mat(i1,i2,isp) = prod_mat(i1,i2,isp)
          end do
       end do
    end do

    ! gcc32: now, invert (I+C*O)
    do isp=1,num_pspecies
       call linalg_invert_serial(inv_mat(1:max_num_proj_tot, &
            1:max_num_proj_tot,isp),ierr)
       if(ierr/=0) then
          write(stdout,'(a)')'======== INVERSION FAILED ========'
       end if
    end do

    ! gcc32: now, finish by taking the product of O_mat and -inv_mat
    do isp=1,num_pspecies
       call linalg_mat_mul_serial(proj_set%O_prec(:,:,isp),O_mat(:,:,isp), &
            inv_mat(:,:,isp),-1.0_DP)
    end do

    ! Deallocate memory
    deallocate(proj_overlap,stat=ierr)
    call utils_dealloc_check('paw_species_calc_proj_prec_mat','proj_overlap', &
         ierr)
    deallocate(O_mat,stat=ierr)
    call utils_dealloc_check('paw_species_calc_proj_prec_mat','O_mat',ierr)
    deallocate(prod_mat,stat=ierr)
    call utils_dealloc_check('paw_species_calc_proj_prec_mat','prod_mat',ierr)
    deallocate(inv_mat,stat=ierr)
    call utils_dealloc_check('paw_species_calc_proj_prec_mat','inv_mat',ierr)

    ! Stop timer
    call timer_clock('paw_species_calc_proj_prec_mat', 2)

  end subroutine paw_species_calc_proj_prec_mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_species_init_core_wvfns(core_wvfns,elements,paw_sp,par,kpt, &
             is_cmplx)

    !=================================================================!
    ! This subroutine allocates and initialises all fftbox_proj_recip !
    ! for each paw_sp element. Each such fftbox_proj_recip is         !
    ! a full non-local projector in reciprocal space in the FFTbox.   !
    !-----------------------------------------------------------------!
    ! Adapted on 21/05/2010.by Nicholas Hine from the routine         !
    ! pseudopot_species_init_proj, originally written by              !
    ! Chris-Kriton Skylaris in 2004 with subsequent modifications     !
    ! by Peter Haynes, Arash Mostofi and Nicholas Hine.               !
    !=================================================================!

    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use projectors, only: PROJECTOR_SET, projectors_allocate_set
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in) :: par
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    type(PROJECTOR_SET), intent(inout) :: core_wvfns
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(POINT), intent(in), optional :: kpt
    logical, intent(in), optional :: is_cmplx

    ! Local Variables
    integer :: iat, orig_iat
    integer :: isp
    integer :: shell
    integer :: proj_count
    type(POINT) :: kpt_loc
    ! agrecocmplx
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &paw_init_core_wvfns'

    ! Check for optional arguments
    if (present(kpt)) then
       kpt_loc = kpt
    else
       kpt_loc%x = 0.0_DP ; kpt_loc%y = 0.0_DP ; kpt_loc%z = 0.0_DP
    end if

    if (present(is_cmplx)) then
       loc_cmplx = is_cmplx
    else
       loc_cmplx = .false.
    end if

    ! ndmh: find number of unique projectors
    proj_count = 0
    do isp=1,par%num_pspecies
       proj_count = proj_count + paw_sp(isp)%n_core_wfs_tot
    end do

    ! Set up the entries of paw_projectors
    core_wvfns%n_proj_species = par%num_pspecies
    ! agrecocmplx
    call projectors_allocate_set(core_wvfns, &
         maxval(paw_sp(:)%n_core_wfs),maxval(paw_sp(:)%n_recip_pts), par, &
         is_cmplx=loc_cmplx)

    ! ndmh: set species_num_proj and species_first_proj values
    proj_count = 1
    do isp=1,core_wvfns%n_proj_species
       core_wvfns%species_num_proj(isp) = paw_sp(isp)%n_core_wfs_tot
       core_wvfns%gmax(isp) = paw_sp(isp)%g_max
       core_wvfns%n_rad_pts(isp) = paw_sp(isp)%n_recip_pts
       core_wvfns%num_shells(isp) = paw_sp(isp)%n_core_wfs
       core_wvfns%ang_mom(:,isp) = 0
       core_wvfns%rad_proj_recip(:,:,isp) = 0.0_DP
       do shell=1,paw_sp(isp)%n_core_wfs
          core_wvfns%ang_mom(shell,isp) = paw_sp(isp)%l_core_wf(shell)
          core_wvfns%rad_proj_recip(1:paw_sp(isp)%n_recip_pts,shell,isp) = &
               paw_sp(isp)%core_wf_recip(1:paw_sp(isp)%n_recip_pts,shell)
       end do
       core_wvfns%species_first_proj(isp) = proj_count
       proj_count = proj_count + core_wvfns%species_num_proj(isp)
    end do

    ! ndmh: copy projector centre and radius from elements array
    do iat=1,par%nat
       orig_iat = par%orig_atom(iat)
       core_wvfns%proj_centre(iat) = elements(orig_iat)%centre
       core_wvfns%proj_max_radius(iat) = elements(orig_iat)%max_core_wf_radius
       core_wvfns%proj_species(iat) = elements(orig_iat)%pspecies_number
    end do

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &paw_init_core_wvfns'

  end subroutine paw_species_init_core_wvfns


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_species_exit_core_wvfns(core_wvfns,paw_sp)

    !====================================================================!
    ! This subroutine deallocates storage for the core wavefunctions     !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  core_wvfns (inout) : PROJECTOR_SET type describing core wvfns.    !
    !--------------------------------------------------------------------!
    ! Written by Nicholas Hine on 16/07/2012.                            !
    !====================================================================!

    use projectors, only: PROJECTOR_SET, projectors_deallocate_set

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(PROJECTOR_SET), intent(inout) :: core_wvfns

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &paw_species_exit_core_wvfns'

    if (any(paw_sp(:)%core_wvfns_exist)) then
       call projectors_deallocate_set(core_wvfns)
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &paw_species_exit_core_wvfns'

  end subroutine paw_species_exit_core_wvfns


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_projector_overlap(proj_overlap,paw_sp)

    !==================================================================!
    ! This subroutine returns the PAW-projector overlap matrix in      !
    ! SPAM3 format.                                                    !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  proj_overlap (inout) : The SPAM3 block-diagonal overlap matrix  !
    !  of the partial waves of each atom                               !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 06/06/10.           !
    ! Modified to get par from SPAM3 by Robert Charlton, 29/08/2018.   !
    !==================================================================!

    use comms, only: pub_my_proc_id
    use gaunt_coeff, only: realgaunt
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_put_block, sparse_get_par
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(SPAM3), intent(inout) :: proj_overlap

    ! Local Variables
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,mi,lj,mj
    integer :: ierr
    real(kind=DP),allocatable :: ovlp_block(:,:)
    ! agrecocmplx
    complex(kind=DP), allocatable :: ovlp_block_cmplx(:,:)
    logical :: iscmplx
    type(PARAL_INFO), pointer :: par

    ! agrecocmplx
    iscmplx = proj_overlap%iscmplx

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, proj_overlap)
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_projector_overlap: allocated parallel strategy is &
         &incompatible with paw_sp.')

    allocate(ovlp_block(max_paw_proj_tot,max_paw_proj_tot),stat=ierr)
    call utils_alloc_check('paw_projector_overlap','ovlp_block',ierr)

    ! agrecocmplx
    if (iscmplx) then
       allocate(ovlp_block_cmplx(max_paw_proj_tot,max_paw_proj_tot),stat=ierr)
       call utils_alloc_check('paw_projector_overlap','ovlp_block_cmplx',ierr)
    end if

    ! Loop over atoms on this proc
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       ! Double loop over partial waves i,j
       do ipwtot=1,paw_sp(isp)%npw_tot
          ipw = paw_sp(isp)%ipw_tot(ipwtot)
          li = paw_sp(isp)%l_pw_tot(ipwtot)
          mi = paw_sp(isp)%m_pw_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find O_ij = sqrt(4pi)*n_nilinjlj^0*G_limiljmj^00
             ovlp_block(ipwtot,jpwtot) = &
                  paw_sp(isp)%aug_nLij(ipw,jpw,0) * &
                  realgaunt(0,0,li,mi,lj,mj) * sqrt_4pi
          end do
       end do

       ! Put this atom's block into SPAM3
       ! agrecocmplx
       if (iscmplx) then
          ovlp_block_cmplx(:,:) = cmplx(ovlp_block(:,:),kind=DP)
          call sparse_put_block(ovlp_block_cmplx,proj_overlap,iat,iat)
       else
          call sparse_put_block(ovlp_block,proj_overlap,iat,iat)
       end if
    end do

    ! agrecocmplx
    if (iscmplx) then
       deallocate(ovlp_block_cmplx, stat=ierr)
       call utils_dealloc_check('paw_projector_overlap','ovlp_block_cmplx',ierr)
    end if
    nullify(par)

    deallocate(ovlp_block,stat=ierr)
    call utils_dealloc_check('paw_projector_overlap','ovlp_block',ierr)

  end subroutine paw_projector_overlap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_projector_overlap_basic(proj_overlap,paw_sp,num_pspecies)

    ! gcc32

    !====================================================================!
    ! This subroutine returns the partial-wave overlap-difference matrix !
    ! Oij in a species-dependent real matrix format                      !
    !--------------------------------------------------------------------!
    !====================================================================!

    use gaunt_coeff, only: realgaunt
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in)          :: num_pspecies
    real(kind=DP), intent(inout) :: proj_overlap(max_paw_proj_tot, &
                    max_paw_proj_tot,num_pspecies)
    type(PAW_SPECIES), intent(in) :: paw_sp(num_pspecies)

    ! Local Variables
    integer :: isp
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,mi,lj,mj


    proj_overlap = 0.0_DP

    ! Loop over species
    do isp=1,num_pspecies

       ! Double loop over partial waves i,j
       do ipwtot=1,paw_sp(isp)%npw_tot
          ipw = paw_sp(isp)%ipw_tot(ipwtot)
          li = paw_sp(isp)%l_pw_tot(ipwtot)
          mi = paw_sp(isp)%m_pw_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find O_ij = sqrt(4pi)*n_nilinjlj^0*G_limiljmj^00
             proj_overlap(ipwtot,jpwtot,isp) = &
                  paw_sp(isp)%aug_nLij(ipw,jpw,0) * &
                  realgaunt(0,0,li,mi,lj,mj) * sqrt_4pi
          end do
       end do
    end do

  end subroutine paw_projector_overlap_basic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_position_operator(pos_ij,paw_sp,axis)

    !==================================================================!
    ! This subroutine returns the PAW sphere part of the position      !
    ! operator, defined as <phi_i|r|phi_j> - <tphi_i|r|tphi_j> for     !
    ! each pair of partial waves i,j                                   !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  pos_ij(3) (inout) : The SPAM3 block-diagonal position matrix    !
    !  between the partial waves of each atom                          !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 06/10/11.           !
    ! Modified to get par from SPAM3 by Robert Charlton, 29/08/2018.   !
    !==================================================================!

    use comms, only: pub_my_proc_id
    use gaunt_coeff, only: realgaunt
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_put_block, sparse_get_par
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(SPAM3), intent(inout) :: pos_ij(3)
    integer, intent(in), optional :: axis

    ! Local Variables
    integer :: axmin,axmax
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,mi,lj,mj,lr,mr
    integer :: icart
    integer :: ierr
    real(kind=DP),allocatable :: pos_block(:,:,:)
    type(PARAL_INFO), pointer :: par

    ! Optional variable
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, pos_ij(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_position_operator: allocated parallel strategy is &
         &incompatible with paw_sp.')

    allocate(pos_block(max_paw_proj_tot,max_paw_proj_tot,3),stat=ierr)
    call utils_alloc_check('paw_position_operator','pos_block',ierr)

    ! Expansion of r in spherical harmonics has only L=1 components
    lr = 1

    ! Loop over atoms on this proc
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       ! Double loop over partial waves i,j
       do ipwtot=1,paw_sp(isp)%npw_tot
          ipw = paw_sp(isp)%ipw_tot(ipwtot)
          li = paw_sp(isp)%l_pw_tot(ipwtot)
          mi = paw_sp(isp)%m_pw_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find r_ij = sqrt(4pi)*n_nilinjlj^0*G_limiljmj^1M
             do icart=axmin,axmax
                if (icart==1) mr = 1
                if (icart==2) mr = -1
                if (icart==3) mr = 0
                pos_block(ipwtot,jpwtot,icart) = &
                     paw_sp(isp)%aug_nLij(ipw,jpw,lr) * &
                     realgaunt(lr,mr,li,mi,lj,mj) * sqrt_4pi / sqrt(3.0_DP)
             end do
          end do
       end do

       ! Put this atom's blocks into SPAM3 matrices
       do icart=axmin,axmax
          call sparse_put_block(pos_block(:,:,icart),pos_ij(icart),iat,iat)
       end do
    end do
    nullify(par)

    deallocate(pos_block,stat=ierr)
    call utils_dealloc_check('paw_position_operator','pos_block',ierr)

  end subroutine paw_position_operator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_grad_operator(grad_ij,paw_sp,axis)

    !==================================================================!
    ! This subroutine returns the PAW sphere part of the grad          !
    ! operator, defined as <phi_i|nabla|phi_j> - <tphi_i|nabla|tphi_j> !
    ! for each pair of partial waves i,j                               !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  grad_ij(3) (inout) : The SPAM3 block-diagonal grad operator     !
    !  matrix between the partial waves of each atom                   !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 02/12/11.           !
    ! Modified to get par from SPAM3 by Robert Charlton, 29/08/2018.   !
    !==================================================================!

    use comms, only: pub_my_proc_id
    use gaunt_coeff, only: realgaunt, gaunt_grad_integrals
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_put_block, sparse_get_par
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(SPAM3), intent(inout) :: grad_ij(3)
    integer, intent(in), optional :: axis

    ! Local Variables
    integer :: axmin,axmax
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,mi,lj,mj
    integer :: icart
    integer :: ierr
    real(kind=DP),allocatable :: grad_block(:,:,:)
    real(kind=DP) :: angint(3,3)
    type(PARAL_INFO), pointer :: par

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, grad_ij(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_position_operator: allocated parallel strategy is &
         &incompatible with paw_sp.')

    ! Optional variable
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    allocate(grad_block(max_paw_proj_tot,max_paw_proj_tot,3),stat=ierr)
    call utils_alloc_check('paw_grad_operator','grad_block',ierr)

    ! Loop over atoms on this proc
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       ! Double loop over partial waves i,j
       do ipwtot=1,paw_sp(isp)%npw_tot
          ipw = paw_sp(isp)%ipw_tot(ipwtot)
          li = paw_sp(isp)%l_pw_tot(ipwtot)
          mi = paw_sp(isp)%m_pw_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find grad_ij
             call gaunt_grad_integrals(angint,li,mi,lj,mj)
             do icart=axmin,axmax
                grad_block(ipwtot,jpwtot,icart) = &
                     paw_sp(isp)%aug_nLij(ipw,jpw,-2) * angint(icart,1) &
                     + paw_sp(isp)%aug_nLij(ipw,jpw,-1) * (-angint(icart,1) &
                     + angint(icart,2) + angint(icart,3))
             end do
          end do
       end do

       ! Put this atom's blocks into SPAM3 matrices
       do icart=axmin,axmax
          call sparse_put_block(grad_block(:,:,icart),grad_ij(icart),iat,iat)
       end do
    end do
    nullify(par)

    deallocate(grad_block,stat=ierr)
    call utils_dealloc_check('paw_grad_operator','grad_block',ierr)

  end subroutine paw_grad_operator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine paw_core_pw_overlap(core_pw_overlap,paw_sp)

    !==================================================================!
    ! This subroutine returns the overlap matrix between core orbitals !
    ! and PAW partial waves, in SPAM3 format.                          !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  core_pw_overlap (inout) : The SPAM3 block-diagonal overlap      !
    !  matrix between the AE partial waves and the AE core orbitals    !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 17/07/12.           !
    ! Modified to get par from SPAM3 by Robert Charlton, 29/08/2018.   !
    !==================================================================!

    use comms, only: pub_my_proc_id
    use gaunt_coeff, only: realgaunt
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_put_block, sparse_get_par
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(SPAM3), intent(inout) :: core_pw_overlap

    ! Local Variables
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: icwf,jpw
    integer :: icwftot,jpwtot
    integer :: li,mi,lj,mj
    integer :: ierr
    real(kind=DP),allocatable :: ovlp_block(:,:)
    type(PARAL_INFO), pointer :: par

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, core_pw_overlap)
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_position_operator: allocated parallel strategy is &
         &incompatible with paw_sp.')

    allocate(ovlp_block(max_core_wf_tot,max_paw_proj_tot),stat=ierr)
    call utils_alloc_check('paw_core_pw_overlap','ovlp_block',ierr)

    ! Loop over atoms on this proc
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       ! Double loop over partial waves i,j
       do icwftot=1,paw_sp(isp)%n_core_wfs_tot
          icwf = paw_sp(isp)%icore_wf_tot(icwftot)
          li = paw_sp(isp)%l_core_wf_tot(icwftot)
          mi = paw_sp(isp)%m_core_wf_tot(icwftot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find O_ij = sqrt(4pi)*n_nilinjlj^0*G_limiljmj^00
             ovlp_block(icwftot,jpwtot) = &
                  paw_sp(isp)%core_aug_nLij(icwf,jpw,0) * &
                  realgaunt(0,0,li,mi,lj,mj) * sqrt_4pi
          end do
       end do

       ! Put this atom's block into SPAM3 - check first if the block has a
       ! non-zero dimension
       if (paw_sp(isp)%n_core_wfs_tot==0 .or. paw_sp(isp)%npw_tot==0) cycle
       call sparse_put_block(ovlp_block,core_pw_overlap,iat,iat)
    end do
    nullify(par)

    deallocate(ovlp_block,stat=ierr)
    call utils_dealloc_check('paw_core_pw_overlap','ovlp_block',ierr)

  end subroutine paw_core_pw_overlap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_core_pw_position_operator(core_pw_pos,paw_sp)

    !==================================================================!
    ! This subroutine returns the PAW sphere part of the position      !
    ! operator between partial waves and core wavefunctions, defined   !
    ! as <phi_c|r|phi_j> - <phi_c|r|tphi_j> for each combination       !
    ! of partial wave j and core wavefunction c                        !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  pos_ij(3) (inout) : The SPAM3 block-diagonal position matrix    !
    !  between the partial waves and core wavefunctions of each atom   !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 06/10/11.           !
    ! Adapted for Core Wavefunctions by Laura Ratcliff.                !
    ! Modified to get par from SPAM3 by Robert Charlton, 29/08/2018.   !
    !==================================================================!

    use comms, only: pub_my_proc_id
    use gaunt_coeff, only: realgaunt
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_put_block, sparse_get_par
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    type(SPAM3), intent(inout) :: core_pw_pos(3)

    ! Local Variables
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: icwf,jpw
    integer :: icwftot,jpwtot
    integer :: li,mi,lj,mj,lr,mr
    integer :: icart
    integer :: ierr
    real(kind=DP),allocatable :: pos_block(:,:,:)
    type(PARAL_INFO), pointer :: par

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, core_pw_pos(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_position_operator: allocated parallel strategy is &
         &incompatible with paw_sp.')

    ! lr408: might need to change this - presumably max_paw_proj_tot will always be smaller
    ! than number of core wfs but really need to add some kind of check
    allocate(pos_block(max_core_wf_tot,max_paw_proj_tot,3),stat=ierr)
    call utils_alloc_check('paw_core_pw_position_operator','pos_block',ierr)

    ! Expansion of r in spherical harmonics has only L=1 components
    lr = 1

    ! Loop over atoms on this proc
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       ! but core wvfns don't exist on each atom so need to check that beforehand!
       if (.not. paw_sp(isp)%core_wvfns_exist) cycle

       ! Double loop over partial waves i,j
       do icwftot=1,paw_sp(isp)%n_core_wfs_tot
          icwf = paw_sp(isp)%icore_wf_tot(icwftot)
          li = paw_sp(isp)%l_core_wf_tot(icwftot)
          mi = paw_sp(isp)%m_core_wf_tot(icwftot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find r_ij = sqrt(4pi)*n_nilinjlj^0*G_limiljmj^1M
             do icart=1,3
                if (icart==1) mr = 1
                if (icart==2) mr = -1
                if (icart==3) mr = 0
                pos_block(icwftot,jpwtot,icart) = &
                     paw_sp(isp)%core_aug_nLij(icwf,jpw,lr) * &
                     realgaunt(lr,mr,li,mi,lj,mj) * sqrt_4pi / sqrt(3.0_DP)
             end do
          end do
       end do

       ! Put this atom's block into SPAM3 - check first if the block has a non-zero dimension
       if (paw_sp(isp)%n_core_wfs_tot==0 .or. paw_sp(isp)%npw_tot==0) cycle
       ! Put this atom's blocks into SPAM3 matrices
       do icart=1,3
          call sparse_put_block(pos_block(:,:,icart),core_pw_pos(icart),iat,iat)
       end do
    end do
    nullify(par)

    deallocate(pos_block,stat=ierr)
    call utils_dealloc_check('paw_core_pw_position_operator','pos_block',ierr)

  end subroutine paw_core_pw_position_operator



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_store_initial_proj_kern(radial_densities,rhoij0, &
       isp,paw_sp)

    !=================================================================!
    ! This subroutine stores the initial guess for the projector      !
    ! density kernel from the atomsolver in an array inside the       !
    ! RADIAL_DENSITY_TYPE structure for the particular species        !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !  rhoij0 (out) : The initial projector density kernel.           !
    !  radial_densities(inout) : RADIAL_DENSITY_TYPE structures       !
    !  isp : species number (not pspecies number) of this species     !
    !  paw_sp : type containing info on PAW species                   !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 17/02/16.                           !
    !=================================================================!

    use ion, only: RADIAL_DENSITY_TYPE
    use rundat, only: pub_num_spins

    implicit none

    ! Arguments
    type(RADIAL_DENSITY_TYPE), intent(inout) :: radial_densities(:)!size(paw_sp)
    real(kind=DP), intent(in) :: rhoij0(:,:,:)
    integer, intent(in) :: isp
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)

    ! Local Variables
    integer :: is
    integer :: ipsp
    integer :: ipw,jpw
    integer :: npw,npwtot
    integer :: ipwtot,jpwtot
    integer :: li,lj,dmi

    ! Initialisation
    ipsp = radial_densities(isp)%ipsp
    npw = paw_sp(ipsp)%npw
    npwtot = paw_sp(ipsp)%npw_tot
    radial_densities(isp)%rhoij0(:,:,:) = 0.0_DP

    ! Loop over spins
    do is=1,pub_num_spins

       ! Loop over partial wave pairs (not counting m-degeneracy)
       ipwtot = 1
       do ipw=1,npw
          li = paw_sp(ipsp)%l_pw(ipw)
          jpwtot = 1
          do jpw=1,npw
             lj = paw_sp(ipsp)%l_pw(jpw)
             if (li==lj) then
                ! expand to include m-degeneracy, store in radial_densities array
                do dmi=-li,li
                    radial_densities(isp)%rhoij0(ipwtot+li+dmi,jpwtot+lj+dmi,is) = &
                         rhoij0(ipw,jpw,is)/real(2*li+1,kind=DP)
                    !print '(9i4,f20.15)',isp,is,ipw,jpw,ipwtot,jpwtot,li,lj,dmi, &
                    !    radial_densities(isp)%rhoij0(ipwtot+li+dmi,jpwtot+lj+dmi,is)
                end do
             end if
             jpwtot = jpwtot + 2*lj + 1
          end do
          ipwtot = ipwtot + 2*li + 1
       end do

    end do

  end subroutine paw_store_initial_proj_kern

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_projector_denskern_init(rhoij0,paw_sp,radial_densities)

    !=================================================================!
    ! This subroutine calculates an initial guess for the projector   !
    ! density kernel, based on the atomic rhoij0 values found in the  !
    ! atomsolver calculation. This is used to calculate the nonlocal  !
    ! energies for the Hamiltonian used in the Palser-Manolopoulos    !
    ! canonical purification.                                         !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !  rhoij0 (out) : The initial projector density kernel.           !
    !  paw_sp (in)  : PAW species type                                !
    !  radial_densities (in) : RADIAL_DENSITY_TYPE storing initial    !
    !                          guess densities and projector denskern !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 04/06/10.                           !
    ! Re-written by Nicholas Hine, 17/02/16 to use atomsolver.        !
    ! Bugfix by Nicholas Hine, 20/02/17 to use dataset versions in    !
    ! some cases                                                      !
    ! Modified to get par from SPAM3 by Robert Charlton, 29/08/2018.  !
    !=================================================================!

    use comms, only: pub_my_proc_id
    use ion, only: RADIAL_DENSITY_TYPE
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins, pub_spin_fac, pub_devel_code, &
         pub_initial_dens_realspace
    use rundat, only: pub_num_spins, pub_spin_fac, pub_devel_code
    use sparse, only: SPAM3, sparse_put_block, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_devel_code, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: rhoij0(pub_num_spins)
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)
    type(RADIAL_DENSITY_TYPE), intent(in), allocatable :: &
         radial_densities(:)

    ! Local Variables
    integer :: loc_iat, iat
    integer :: isp, ipsp
    integer :: is
    integer :: ipwtot,jpwtot
    integer :: ierr
    real(kind=DP), allocatable :: rhoij0_sp(:,:,:,:)
    logical :: use_dataset_rhoij
    type(PARAL_INFO), pointer :: par

    ! Allow over-ride of use of rhoij0 from dataset with devel code
    use_dataset_rhoij = utils_devel_code(.false.,'PAW', &
         'USE_DATASET_RHOIJ0',pub_devel_code)

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, rhoij0(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_projector_denskern_init: allocated parallel strategy is &
         &incompatible with paw_sp.')
    call utils_assert(par%num_species == size(radial_densities), 'Error in &
         &paw_projector_denskern_init: allocated parallel strategy is &
         &incompatible with radial_densities.')

    ! Avoid breaking if real space densities not initialised
    ! (use dataset versions instead)
    if (.not.pub_initial_dens_realspace) then
       use_dataset_rhoij = .true.
    end if

    ! Initialisation
    allocate(rhoij0_sp(max_paw_proj_tot,max_paw_proj_tot, &
         pub_num_spins,par%num_species),stat=ierr)
    call utils_alloc_check('paw_projector_denskern_init','rhoij0_sp',ierr)
    rhoij0_sp(:,:,:,:) = 0.0_DP

    ! Loop over species, setting up array of rhoij0 terms
    do isp=1,par%num_species
       if (.not.use_dataset_rhoij) then
          ipsp = radial_densities(isp)%ipsp
       else
          ipsp = 0
          do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
             if (par%elements_on_proc(loc_iat)%species_number==isp) then
                ipsp = par%elements_on_proc(loc_iat)%pspecies_number
             end if
          end do
          if (ipsp==0) cycle
       end if
       do is=1,pub_num_spins
          do ipwtot=1,paw_sp(ipsp)%npw_tot
             do jpwtot=1,paw_sp(ipsp)%npw_tot
                if (use_dataset_rhoij) then
                   ! old version, used rhoij0 from dataset
                   rhoij0_sp(ipwtot,jpwtot,is,isp) = &
                        paw_sp(ipsp)%rhoij0(ipwtot,jpwtot) * 0.5_DP * pub_spin_fac
                else
                   ! new version, uses rhoij0 from atomsolver
                    rhoij0_sp(ipwtot,jpwtot,is,isp) = &
                         radial_densities(isp)%rhoij0(ipwtot,jpwtot,is)
                end if
             end do
          end do
       end do

    end do

    ! Loop over atoms
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
       isp = par%elements_on_proc(loc_iat)%species_number
       ! Put block of rhoij0 into SPAM3 matrix
       do is=1,pub_num_spins
          call sparse_put_block(rhoij0_sp(:,:,is,isp),rhoij0(is),iat,iat)
       end do
    end do

    ! Cleanup
    deallocate(rhoij0_sp,stat=ierr)
    call utils_dealloc_check('paw_projector_denskern_init','rhoij0_sp',ierr)
    nullify(par)

  end subroutine paw_projector_denskern_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_atom_aug_den(atom_nhat,atom_grad_nhat,total_nhat, &
       total_nhat_targ,rho_ij_block, &
       isp,atom_centre,grid,cell,paw_sp,box, &
       box_start1,box_start2,box_start3,atom_aug_func_real, &
       atom_aug_func_grad_real,atom_aug_func_recip,atom_aug_func_grad_recip, &
       atom_nhat_recip,atom_grad_nhat_recip,grad_cartgrad,input_cart)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use gaunt_coeff, only: realgaunt
    use geometry, only: POINT
    use paw_shape, only: paw_shape_gLSLM_recip, paw_shape_gLSLM_real, &
         paw_shape_grad_gLSLM_real, paw_shape_grad_gLSLM_recip, &
         paw_shape_grad_cartgrad_gLSLM_recip
    use rundat, only: pub_aug_funcs_recip, pub_num_spins
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort ! gcc32
    use xc, only: pub_xc_gradient_corrected

    implicit none

    ! Arguments
    type(FFTBOX_INFO) :: box
    real(kind=DP), intent(out) :: atom_nhat(box%total_pt1,box%total_pt2, &
         box%total_pt3,pub_num_spins)
    real(kind=DP), intent(out) :: atom_grad_nhat(box%total_pt1,box%total_pt2, &
         box%total_pt3,pub_num_spins,3)
    real(kind=DP), intent(out) :: total_nhat(max_spins)
    real(kind=DP), intent(out) :: total_nhat_targ(max_spins)
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)
    type(POINT), intent(in) :: atom_centre
    real(kind=DP), intent(in) :: rho_ij_block(:,:,:)
    integer, intent(in) :: isp
    integer, intent(in) :: box_start1,box_start2,box_start3
    real(kind=DP), intent(out), optional :: atom_aug_func_real(box%total_pt1, &
         box%total_pt2, box%total_pt3)
    real(kind=DP), intent(out), optional :: atom_aug_func_grad_real(box%total_pt1, &
         box%total_pt2, box%total_pt3,3)
    complex(kind=DP), intent(out), optional :: atom_aug_func_recip(box%total_pt1, &
         box%total_pt2, box%total_pt3)
    complex(kind=DP), intent(out), optional :: atom_aug_func_grad_recip(box%total_pt1, &
         box%total_pt2, box%total_pt3,3)
    complex(kind=DP), intent(out), optional :: atom_nhat_recip(box%total_pt1, &
         box%total_pt2, box%total_pt3,pub_num_spins)
    complex(kind=DP), intent(out), optional :: atom_grad_nhat_recip(box%total_pt1, &
         box%total_pt2, box%total_pt3,pub_num_spins,3)
    integer, intent(in), optional :: grad_cartgrad ! gcc32: optional tag
           ! used for calculating d/dr (d/dR_I)_cart of \hat{\rho_I},
           ! instead of the standard gradient ( d/dr \hat{\rho_I} )
    integer, intent(in), optional :: input_cart

    ! Local Variables
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: lup,mup
    integer :: li,mi,lj,mj
    integer :: is
    integer :: cart
    real(kind=DP) :: rgij
    real(kind=DP) :: qijLM
    real(kind=DP) :: qLM(max_spins)
    real(kind=DP) :: total_gLSLM

    if (pub_aug_funcs_recip) then
       atom_nhat_recip = 0.0_DP
       if (pub_xc_gradient_corrected) atom_grad_nhat_recip = 0.0_DP
    else
       atom_nhat = 0.0_DP
       if (pub_xc_gradient_corrected) atom_grad_nhat = 0.0_DP
    end if

    total_nhat = 0.0_DP
    total_nhat_targ = 0.0_DP

    do lup=0,paw_sp(isp)%lmax
       do mup=-lup,lup

          qLM(:) = 0.0_DP

          ! Loop over partial waves to find total target comp. charge
          ! and calculate q_LM^\sigma for this L,M,sigma channel
          do ipwtot=1,paw_sp(isp)%npw_tot
             ipw = paw_sp(isp)%ipw_tot(ipwtot)
             li = paw_sp(isp)%l_pw_tot(ipwtot)
             mi = paw_sp(isp)%m_pw_tot(ipwtot)
             do jpwtot=1,paw_sp(isp)%npw_tot
                jpw = paw_sp(isp)%ipw_tot(jpwtot)
                lj = paw_sp(isp)%l_pw_tot(jpwtot)
                mj = paw_sp(isp)%m_pw_tot(jpwtot)
                rgij = realgaunt(lup,mup,lj,mj,li,mi)
                qijLM = paw_sp(isp)%aug_nLij(ipw,jpw,lup)*rgij
                do is=1,pub_num_spins
                   qLM(is) = qLM(is) + &
                        rho_ij_block(ipwtot,jpwtot,is) * qijLM
                   if (lup==0) &
                        total_nhat_targ(is) = total_nhat_targ(is) + &
                        qijLM*rho_ij_block(ipwtot,jpwtot,is) * sqrt_4pi
                end do

             end do
          end do

          if (any(abs(qLM(:))>1e-20_DP)) then

             ! Get shape function g_L(r) multiplied by spherical harmonic
             ! S_LM(r) for this L,M pair
             if (pub_aug_funcs_recip) then
                ! Get augmentation function in reciprocal space
                call paw_shape_gLSLM_recip(atom_aug_func_recip, &
                     paw_sp(isp)%shape, &
                     lup,mup,box_start1,box_start2,box_start3,grid,cell, &
                     box%total_pt1, box%total_pt2, box%total_pt3,atom_centre)
                total_gLSLM = real(atom_aug_func_recip(1,1,1),kind=DP)

                if (.not.present(input_cart)) then

                   do is=1,pub_num_spins
                      atom_nhat_recip(:,:,:,is) = atom_nhat_recip(:,:,:,is) &
                           + qLM(is)*atom_aug_func_recip(:,:,:)
                      total_nhat(is) = total_nhat(is) + total_gLSLM*qLM(is)
                   end do

                end if

                ! Find gradient of \hat{n}(r) if required
                if (pub_xc_gradient_corrected) then
                  if(.not.present(grad_cartgrad)) then ! by gcc32
                     ! Get gradient of (shape function g_L(r) multiplied by
                     ! spherical harmonic S_LM(r) for this L,M pair)
                     call paw_shape_grad_gLSLM_recip(atom_aug_func_grad_recip, &
                          paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                          box_start3,grid,cell,box%total_pt1, &
                          box%total_pt2,box%total_pt3,atom_centre)
                     do cart=1,3
                        do is=1,pub_num_spins
                           atom_grad_nhat_recip(:,:,:,is,cart) = &
                                atom_grad_nhat_recip(:,:,:,is,cart) + &
                                qLM(is)*atom_aug_func_grad_recip(:,:,:,cart)
                        end do
                     end do
                   else
                     ! Get gradient of cartesian component grad_cartgrad of
                     ! the gradient of (shape function g_L(r) multiplied by
                     ! spherical harmonic S_LM(r) for this L,M pair)
                     call paw_shape_grad_cartgrad_gLSLM_recip( &
                          atom_aug_func_grad_recip, &
                          paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                          box_start3,grid,cell,box%total_pt1, &
                          box%total_pt2, box%total_pt3, &
                          atom_centre, grad_cartgrad)
                     do cart=1,3
                        do is=1,pub_num_spins
                           atom_grad_nhat_recip(:,:,:,is,cart) = &
                                atom_grad_nhat_recip(:,:,:,is,cart) + &
                                qLM(is)*atom_aug_func_grad_recip(:,:,:,cart)
                        end do
                     end do
                   endif ! if .not.present(grad_cartgrad)
                end if ! if (pub_xc_gradient_corrected)

                if ((.not.pub_xc_gradient_corrected).and.present(input_cart)) then
                   ! Get gradient of (shape function g_L(r) multiplied by
                   ! spherical harmonic S_LM(r) for this L,M pair)
                   call paw_shape_grad_gLSLM_recip(atom_aug_func_grad_recip, &
                        paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                        box_start3,grid,cell,box%total_pt1, &
                        box%total_pt2,box%total_pt3,atom_centre)
                   do is=1,pub_num_spins
                      atom_nhat_recip(:,:,:,is) = atom_nhat_recip(:,:,:,is) + &
                           qLM(is)*atom_aug_func_grad_recip(:,:,:,input_cart)
                   end do
                end if


             else

                ! Get augmentation function in real space
                call paw_shape_gLSLM_real(atom_aug_func_real,paw_sp(isp)%shape, &
                     lup,mup,box_start1,box_start2,box_start3,grid,cell, &
                     box%total_pt1,box%total_pt2,box%total_pt3, &
                     atom_centre)
                total_gLSLM = sum(atom_aug_func_real(:,:,:))*grid%weight

                do is=1,pub_num_spins
                   atom_nhat(:,:,:,is) = atom_nhat(:,:,:,is) &
                        + qLM(is)*atom_aug_func_real(:,:,:)
                   total_nhat(is) = total_nhat(is) + total_gLSLM*qLM(is)
                end do

                ! Find gradient of \hat{n}(r) if required
                if (pub_xc_gradient_corrected) then
                   if(.not.present(grad_cartgrad)) then ! by gcc32
                      ! Get gradient of (shape function g_L(r) multiplied by
                      ! spherical harmonic S_LM(r) for this L,M pair)
                      call paw_shape_grad_gLSLM_real(atom_aug_func_grad_real, &
                           paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                           box_start3,grid,cell,box%total_pt1, &
                           box%total_pt2,box%total_pt3,atom_centre)
                      do cart=1,3
                         do is=1,pub_num_spins
                            atom_grad_nhat(:,:,:,is,cart) = &
                                 atom_grad_nhat(:,:,:,is,cart) + &
                                 qLM(is)*atom_aug_func_grad_real(:,:,:,cart)
                         end do
                      end do
                  else
                     call utils_abort('Cannot compute gradient of first-order &
                          &augmentation charge in real-space')
                  endif
                end if

             end if
          end if
       end do
    end do

    if (pub_aug_funcs_recip) then

       ! Fourier transform to real space in augmentation box
       do is=1,pub_num_spins
          call fourier_apply_box('F','B',atom_nhat_recip(:,:,:,is),box=box)
          atom_nhat(:,:,:,is) = real(atom_nhat_recip(:,:,:,is),kind=DP) &
               / grid%weight
       end do

       ! Fourier transform gradient of \hat{n}(r) if required
       if (pub_xc_gradient_corrected) then
          do cart=1,3
             do is=1,pub_num_spins
                call fourier_apply_box('F','B', &
                     atom_grad_nhat_recip(:,:,:,is,cart),box=box)
                atom_grad_nhat(:,:,:,is,cart) = &
                     real(atom_grad_nhat_recip(:,:,:,is,cart),kind=DP)/grid%weight
             end do
          end do
       end if
    end if

  end subroutine paw_atom_aug_den


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_atom_aug_integrals(dijhat_at, &
       num_spins,isp,atom_centre,grid,cell,paw_sp,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3,locpot_box_real,atom_aug_func_real, &
       locpot_box_recip,atom_aug_func_recip)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use gaunt_coeff, only: realgaunt
    use geometry, only: POINT
    use paw_shape, only: paw_shape_gLSLM_recip, paw_shape_gLSLM_real
    use rundat, only: pub_aug_funcs_recip
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    integer, intent(in) :: num_spins
    type(CELL_INFO),intent(in) :: cell
    type(GRID_INFO), intent(in) :: grid
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)
    type(POINT), intent(in) :: atom_centre
    real(kind=DP), intent(inout) :: dijhat_at(max_paw_proj_tot, &
         max_paw_proj_tot,num_spins)
    integer, intent(in) :: isp
    integer, intent(in) :: box_start1,box_start2,box_start3
    real(kind=DP), intent(in), optional :: locpot_box_real(box_n1,box_n2, &
         box_n3,num_spins)
    real(kind=DP), intent(out), optional :: atom_aug_func_real(box_n1,box_n2, &
         box_n3)
    complex(kind=DP), intent(in), optional :: locpot_box_recip(box_n1,box_n2, &
         box_n3,num_spins)
    complex(kind=DP), intent(out), optional :: atom_aug_func_recip(box_n1, &
         box_n2,box_n3)

    ! Local Variables
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: lup,mup
    integer :: li,mi,lj,mj
    integer :: is
    real(kind=DP) :: rgij
    real(kind=DP) :: qijLM
    real(kind=DP) :: locpot_gLSLM_product(max_spins)

    ! Check arguments
    if (pub_aug_funcs_recip) then
       if ((.not.present(locpot_box_recip)) .or. &
            (.not.present(atom_aug_func_recip))) then
          call utils_abort('Error in paw_atom_aug_integrals: Incorrect &
               &arguments supplied')
       end if
    else
       if ((.not.present(locpot_box_real)) .or. &
            (.not.present(atom_aug_func_real))) then
          call utils_abort('Error in paw_atom_aug_integrals: Incorrect &
               &arguments supplied')
       end if
    end if

    ! Loop over angular momentum channels L,M
    do lup=0,paw_sp(isp)%lmax
       do mup=-lup,lup

          ! Get shape function g_L(r) multiplied by spherical harmonic
          ! S_LM(r) for this L,M pair and integrate with locpot

          if (pub_aug_funcs_recip) then

             ! Get augmentation function in reciprocal space
             call paw_shape_gLSLM_recip(atom_aug_func_recip,paw_sp(isp)%shape, &
                  lup,mup,box_start1,box_start2,box_start3,grid,cell, &
                  box_n1,box_n2,box_n3,atom_centre)

             ! Calculate integral \int veff(r) g_L(r) S_LM(r) dr
             ! as sum in reciprocal space: \sum_G veff(G) (g_L S_LM)(G)
             do is=1,num_spins
                locpot_gLSLM_product(is) = real(sum(locpot_box_recip(:,:,:,is) &
                     * atom_aug_func_recip(:,:,:)), kind=DP)
             end do

          else

             ! Get augmentation function in real space
             call paw_shape_gLSLM_real(atom_aug_func_real,paw_sp(isp)%shape, &
                  lup,mup,box_start1,box_start2,box_start3,grid,cell, &
                  box_n1,box_n2,box_n3,atom_centre)

             ! Calculate integral \int veff(r) g_L(r) S_LM(r) dr
             do is=1,num_spins
                locpot_gLSLM_product(is) = sum(locpot_box_real(:,:,:,is) &
                     * atom_aug_func_real(:,:,:)) * grid%weight
             end do

          end if

          ! Loop over spins and double loop over partial waves
          do is=1,num_spins
             do ipwtot=1,paw_sp(isp)%npw_tot
                ipw = paw_sp(isp)%ipw_tot(ipwtot)
                li = paw_sp(isp)%l_pw_tot(ipwtot)
                mi = paw_sp(isp)%m_pw_tot(ipwtot)
                do jpwtot=1,paw_sp(isp)%npw_tot
                   jpw = paw_sp(isp)%ipw_tot(jpwtot)
                   lj = paw_sp(isp)%l_pw_tot(jpwtot)
                   mj = paw_sp(isp)%m_pw_tot(jpwtot)
                   rgij = realgaunt(lup,mup,li,mi,lj,mj)
                   qijLM = paw_sp(isp)%aug_nLij(ipw,jpw,lup)*rgij

                   dijhat_at(ipwtot,jpwtot,is) = &
                        dijhat_at(ipwtot,jpwtot,is) + &
                        locpot_gLSLM_product(is) * qijLM

                end do  ! jpwtot
             end do  ! ipwtot
          end do  ! is

       end do  ! mup
    end do  ! lup

  end subroutine paw_atom_aug_integrals


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_atom_grad_aug_integrals(dijhat_at, &
       num_spins,isp,atom_centre,grid,cell,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3,cart,paw_sp,locpot_box_real, &
       atom_grad_aug_func_real,locpot_box_recip,atom_grad_aug_func_recip, &
       cart2)

    !==================================================================!
    ! This subroutine calculates the integral between a local          !
    ! potential box and the gradient (along direction 'cart') of the   !
    ! multipole moment Q_ij(\mathbf{r}) associated to an atom with     !
    ! species 'isp':                                                   !
    ! \sum_LM locpot_box(\mathbf{r}) [\nabla Q^LM_ij]_cart d\mathbf{r} !
    !------------------------------------------------------------------!
    !------------------------------------------------------------------!
    ! This subroutine was written by Gabriel Constantinescu in May     !
    ! 2015                                                             !
    !==================================================================!


    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use gaunt_coeff, only: realgaunt
    use geometry, only: POINT
    use paw_shape, only: paw_shape_grad_gLSLM_recip, paw_shape_grad_gLSLM_real
    use rundat, only: pub_aug_funcs_recip
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_abort

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    integer, intent(in) :: num_spins
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO),intent(in) :: cell
    type(POINT), intent(in) :: atom_centre
    real(kind=DP), intent(inout) :: dijhat_at(max_paw_proj_tot, &
         max_paw_proj_tot,num_spins)
    integer, intent(in) :: isp
    integer, intent(in) :: box_start1,box_start2,box_start3
    real(kind=DP), intent(in), optional :: locpot_box_real(box_n1,box_n2, &
         box_n3,num_spins)
    real(kind=DP), intent(out), optional :: atom_grad_aug_func_real(box_n1, &
         box_n2, box_n3,3)
    complex(kind=DP), intent(in), optional :: locpot_box_recip(box_n1,box_n2, &
         box_n3,num_spins)
    complex(kind=DP), intent(out), optional :: atom_grad_aug_func_recip(box_n1,&
         box_n2,box_n3,3)
    integer, intent(in) :: cart
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    integer, intent(in), optional :: cart2

    ! Local Variables
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: lup,mup
    integer :: li,mi,lj,mj
    integer :: is
    real(kind=DP) :: rgij
    real(kind=DP) :: qijLM
    real(kind=DP) :: locpot_grad_gLSLM_product(max_spins)


    ! Check arguments
    if (pub_aug_funcs_recip) then
       if ((.not.present(locpot_box_recip)) .or. &
            (.not.present(atom_grad_aug_func_recip))) then
          call utils_abort('Error in paw_atom_grad_aug_integrals: Incorrect &
               &arguments supplied')
       end if
    else
       if ((.not.present(locpot_box_real)) .or. &
            (.not.present(atom_grad_aug_func_real))) then
          call utils_abort('Error in paw_atom_grad_aug_integrals: Incorrect &
               &arguments supplied')
       end if
    end if

    ! Loop over angular momentum channels L,M
    do lup = 0,paw_sp(isp)%lmax
       do mup = -lup,lup

          ! Get shape function g_L(r) multiplied by spherical harmonic
          ! S_LM(r) for this L,M pair and integrate with locpot

          if (pub_aug_funcs_recip) then

             if (present(cart2)) then ! second-order perturbation of \hat{Q}_ab
                ! Get augmentation function in reciprocal space
                call paw_shape_grad_gLSLM_recip(atom_grad_aug_func_recip, &
                     paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                     box_start3,grid,cell,box_n1,box_n2,box_n3,atom_centre, &
                     cart2)
             else
                ! Get augmentation function in reciprocal space
                call paw_shape_grad_gLSLM_recip(atom_grad_aug_func_recip, &
                     paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                     box_start3,grid,cell,box_n1,box_n2,box_n3,atom_centre)
             end if

             ! Calculate integral \int veff(r) d/dr(g_L(r) S_LM(r)) dr
             ! as sum in reciprocal space: \sum_G veff(G) [d/dr(g_L S_LM)](G)
             do is = 1,num_spins
                locpot_grad_gLSLM_product(is) = &
                     real(sum(locpot_box_recip(:,:,:,is) &
                     * atom_grad_aug_func_recip(:,:,:,cart)), kind=DP)
             end do

          else

             ! Get augmentation function in real space
             call paw_shape_grad_gLSLM_real(atom_grad_aug_func_real, &
                  paw_sp(isp)%shape, lup,mup,box_start1,box_start2, &
                  box_start3,grid,cell,box_n1,box_n2,box_n3,atom_centre)

             ! Calculate integral \int veff(r) d/dr[g_L(r) S_LM(r)] dr
             do is = 1,num_spins
                locpot_grad_gLSLM_product(is) = sum(locpot_box_real(:,:,:,is) &
                     * atom_grad_aug_func_real(:,:,:,cart)) * grid%weight
             end do

          end if

          ! Loop over spins and double loop over partial waves
          do is = 1,num_spins
             do ipwtot=1,paw_sp(isp)%npw_tot
                ipw = paw_sp(isp)%ipw_tot(ipwtot)
                li = paw_sp(isp)%l_pw_tot(ipwtot)
                mi = paw_sp(isp)%m_pw_tot(ipwtot)
                do jpwtot = 1,paw_sp(isp)%npw_tot
                   jpw = paw_sp(isp)%ipw_tot(jpwtot)
                   lj = paw_sp(isp)%l_pw_tot(jpwtot)
                   mj = paw_sp(isp)%m_pw_tot(jpwtot)
                   rgij = realgaunt(lup,mup,li,mi,lj,mj)
                   qijLM = paw_sp(isp)%aug_nLij(ipw,jpw,lup) * rgij

                   dijhat_at(ipwtot,jpwtot,is) = &
                        dijhat_at(ipwtot,jpwtot,is) + &
                        locpot_grad_gLSLM_product(is) * qijLM

                end do  ! jpwtot
             end do  ! ipwtot
          end do  ! is

       end do  ! mup
    end do  ! lup

  end subroutine paw_atom_grad_aug_integrals


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_atom_aug_force(nhat_force, &
       rhoij_at,num_spins,isp,atom_centre,grid,cell,paw_sp,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3,locpot_box_real, &
       atom_grad_aug_func_real,locpot_box_recip,atom_grad_aug_func_recip)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use gaunt_coeff, only: realgaunt
    use geometry, only: POINT
    use paw_shape, only: paw_shape_grad_gLSLM_real, paw_shape_grad_gLSLM_recip
    use rundat, only: pub_aug_funcs_recip
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    real(kind=DP),intent(inout) :: nhat_force(3)
    integer,intent(in) :: box_n1,box_n2,box_n3
    integer,intent(in) :: num_spins
    type(CELL_INFO),intent(in) :: cell
    type(GRID_INFO),intent(in) :: grid
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)
    type(POINT),intent(in) :: atom_centre
    real(kind=DP),intent(in) :: rhoij_at(:,:,:)
    integer,intent(in) :: isp
    integer,intent(in) :: box_start1,box_start2,box_start3
    real(kind=DP),intent(out),optional :: atom_grad_aug_func_real(box_n1, &
         box_n2,box_n3,3)
    real(kind=DP),intent(in),optional :: locpot_box_real(box_n1,box_n2, &
         box_n3,num_spins)
    complex(kind=DP),intent(out),optional :: atom_grad_aug_func_recip(box_n1, &
         box_n2,box_n3,3)
    complex(kind=DP),intent(in),optional :: locpot_box_recip(box_n1,box_n2, &
         box_n3,num_spins)

    ! Local Variables
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: lup,mup
    integer :: li,mi,lj,mj
    integer :: is
    integer :: cart
    real(kind=DP) :: rgij
    real(kind=DP) :: qijLM
    real(kind=DP) :: locpot_grad_gLSLM_product(3,max_spins)

    ! Loop over angular momentum channels L,M
    do lup=0,paw_sp(isp)%lmax
       do mup=-lup,lup

          ! Construct gradient d/dR_I (g_L(r)*S_LM(\hat{r})) of shape
          ! function g_L(r) multiplied by spherical harmonic S_LM(r) for
          ! this L,M pair, then integrate it with part of locpot in aug box

          if (pub_aug_funcs_recip) then

             call paw_shape_grad_gLSLM_recip(atom_grad_aug_func_recip, &
                  paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                  box_start3,grid,cell,box_n1,box_n2,box_n3, &
                  atom_centre)

             ! Loop over spins and loop over cartesian directions
             do is=1,num_spins
                do cart=1,3
                   ! Calculate integral \int veff(r) d/dR(g_L(r) S_LM(r)) dr
                   locpot_grad_gLSLM_product(cart,is) = &
                        real(sum(locpot_box_recip(:,:,:,is) &
                        * atom_grad_aug_func_recip(:,:,:,cart)),kind=DP)
                end do
             end do

          else

             call paw_shape_grad_gLSLM_real(atom_grad_aug_func_real, &
                  paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                  box_start3,grid,cell,box_n1,box_n2,box_n3, &
                  atom_centre)

             ! Loop over spins and loop over cartesian directions
             do is=1,num_spins
                do cart=1,3
                   ! Calculate integral \int veff(r) d/dR(g_L(r) S_LM(r)) dr
                   locpot_grad_gLSLM_product(cart,is) = &
                        sum(locpot_box_real(:,:,:,is) &
                        * atom_grad_aug_func_real(:,:,:,cart)) * grid%weight
                end do
             end do

          end if

          ! Double loop over partial waves
          do ipwtot=1,paw_sp(isp)%npw_tot
             ipw = paw_sp(isp)%ipw_tot(ipwtot)
             li = paw_sp(isp)%l_pw_tot(ipwtot)
             mi = paw_sp(isp)%m_pw_tot(ipwtot)
             do jpwtot=1,paw_sp(isp)%npw_tot
                jpw = paw_sp(isp)%ipw_tot(jpwtot)
                lj = paw_sp(isp)%l_pw_tot(jpwtot)
                mj = paw_sp(isp)%m_pw_tot(jpwtot)
                rgij = realgaunt(lup,mup,li,mi,lj,mj)
                qijLM = paw_sp(isp)%aug_nLij(ipw,jpw,lup)*rgij

                do is=1,num_spins
                   nhat_force(:) = nhat_force(:) &
                        + (locpot_grad_gLSLM_product(:,is) * qijLM &
                        * rhoij_at(ipwtot,jpwtot,is))
                end do  ! is

             end do  ! jpwtot
          end do  ! ipwtot

       end do  ! mup
    end do  ! lup

  end subroutine paw_atom_aug_force


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine paw_nonlocal_energies(dij,rho_ij,paw_sp,par, &
       paw_sphere_energies,show_matrices,first_order,FO_rho_ij)

    !==================================================================!
    ! This subroutine creates the PAW nonlocal energies Dij, given by  !
    !  D_ij = \hat{D}_ij + D^1_ij - \tilde{D}^1_ij                     !
    !                                                                  !
    ! Can also provide first-order equivalent, if 'FO_rho_ij', 'atom'  !
    ! and 'cart' are present                                           !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  dij (out) : nonlocal matrix in Nproj x Nproj sized matrix       !
    !  rhoij (in) : projector-reduced density kernel (Nproj-sized)     !
    !  grid (in) : Grid definition for fine grid                       !
    !  paw_sphere_energies (inout) : array of energy contributions     !
    !    from sphere terms                                             !
    !  show_matrices (in,opt) : Whether it is a good time to display   !
    !    the intermediate matrices (eg during energy components report)!
    !                                                                  !
    !  FOR FIRST-ORDER:                                                !
    !  FO_rhoij (in) : first-order projector-reduced density kernel    !
    !  atom (in) : perturbed atom                                      !
    !  cart (un) : direction of perturbation                           !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.           !
    ! First-order capabilities were added by Gabriel Constantinescu in !
    ! May 2015                                                         !
    ! Modified to remove pub_par by Joseph Prentice, June 2018         !
    !==================================================================!

    use comms, only: pub_on_root
    use constants, only: paw_en_size, paw_en_dij0, &
         paw_en_ehart, paw_en_exc, paw_en_exc_dc, paw_en_etxc, paw_en_etxc_dc, &
         paw_en_dijxc, paw_en_exc_core
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_paw_output_detail, pub_xc_functional, pub_num_spins, &
         pub_debug_on_root
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_trace, &
         sparse_axpy, sparse_scale
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dij(pub_num_spins)
    type(SPAM3), intent(in) :: rho_ij(pub_num_spins)
    type(PARAL_INFO), intent(in) :: par
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    real(kind=DP), intent(inout), optional :: paw_sphere_energies(paw_en_size)
    logical, intent(in), optional :: show_matrices
    logical, intent(in), optional :: first_order
    type(SPAM3), intent(in), optional :: FO_rho_ij(pub_num_spins)

    ! Local Variables
    type(SPAM3), allocatable :: dij0(:)
    type(SPAM3), allocatable :: dij_hartree(:)
    type(SPAM3), allocatable :: dij_xc(:)
    real(kind=DP) :: exc,etxc,exc_dc,etxc_dc,edijxc,exc_core
    integer :: ierr
    integer :: is
    logical :: loc_show_matrices
    logical :: loc_first_order

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_nonlocal_energies'

    ! check whether we are doing first order response version
    loc_first_order = .false.
    if ( present(first_order).and.present(FO_rho_ij) ) then
       loc_first_order = first_order
    end if

    ! Handle optional arguments
    if (present(paw_sphere_energies)) paw_sphere_energies(:) = 0.0_DP
    loc_show_matrices = .false.
    if (present(show_matrices)) loc_show_matrices = show_matrices

    ! Allocate temporary matrices
    allocate(dij0(pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_nonlocal_energies','dij0',ierr)
    allocate(dij_hartree(pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_nonlocal_energies','dij_hartree',ierr)
    allocate(dij_xc(pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_nonlocal_energies','dij_xc',ierr)

    ! Allocate temporary matrices
    do is=1,pub_num_spins
       call sparse_create(dij0(is),rho_ij(is))
       call sparse_create(dij_hartree(is),rho_ij(is))
       call sparse_create(dij_xc(is),rho_ij(is))
    end do

    if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
       do is=1,pub_num_spins
          if (pub_on_root) write(stdout,'(a,i4)') 'rho_ij', is
          call paw_show_atomblocks(rho_ij(is),paw_sp)
       end do
    end if

    if (.not.loc_first_order) then

       ! Calculate the fixed kinetic and Hartree contributions to the nonlocal
       ! energies. Not relevant for first-order calculation!
       call paw_dij0(dij0,paw_sp)
       do is = 1,pub_num_spins
          if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
             if (pub_on_root) write(stdout,'(a,i4)') 'dij0', is
             call paw_show_atomblocks(dij0(is),paw_sp)
          end if
          call sparse_axpy(dij(is),dij0(is),1.0_DP)
          if (present(paw_sphere_energies)) then
             paw_sphere_energies(paw_en_dij0) = &
                  paw_sphere_energies(paw_en_dij0) + &
                  sparse_trace(rho_ij(is),dij0(is))
          end if
       end do

    end if ! if NOT a first-order calculation

    ! Calculate the remaining Hartree contribution to the nonlocal energies
    if (.not.loc_first_order) then

       call paw_dij_hartree(dij_hartree,rho_ij,paw_sp)
       do is = 1,pub_num_spins
          if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
             if (pub_on_root) write(stdout,'(a,i4)') 'dij_hartree', is
             call paw_show_atomblocks(dij_hartree(is),paw_sp)
          end if
          call sparse_axpy(dij(is),dij_hartree(is),1.0_DP)
          if (present(paw_sphere_energies)) then
             paw_sphere_energies(paw_en_ehart) = &
                  paw_sphere_energies(paw_en_ehart) &
                  + 0.5_DP*sparse_trace(rho_ij(is),dij_hartree(is))
          end if
       end do

    else ! if it is a first-order calculation

       call paw_dij_hartree(dij_hartree,FO_rho_ij,paw_sp)
       do is = 1,pub_num_spins
          call sparse_axpy(dij(is),dij_hartree(is),1.0_DP)
       end do

    end if

    ! Calculate the exchange-correlation part of the nonlocal energies
    if (pub_xc_functional/='NONE') then

       if (.not.loc_first_order) then

          call paw_exc_core_tot(exc_core,paw_sp,par)
          call paw_dij_xc(dij_xc,rho_ij,exc,exc_dc,etxc,etxc_dc,paw_sp)
          edijxc = 0.0_DP
          do is = 1,pub_num_spins
             if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
                if (pub_on_root) write(stdout,'(a,i4)') 'dijxc', is
                call paw_show_atomblocks(dij_xc(is),paw_sp)
             end if
             call sparse_axpy(dij(is),dij_xc(is),1.0_DP)
             edijxc = edijxc + sparse_trace(rho_ij(is),dij_xc(is))
          end do

       else ! if it is a first-order calculation

          call paw_dij_xc(dij_xc,rho_ij,exc,exc_dc,etxc,etxc_dc,paw_sp, &
               loc_first_order, FO_rho_ij)

          do is = 1,pub_num_spins
             call sparse_axpy(dij(is),dij_xc(is),1.0_DP)
          end do

       end if

    else  ! if NO xc contribution
       exc = 0.0_DP; exc_dc = 0.0_DP
       etxc = 0.0_DP; etxc_dc = 0.0_DP
       exc_core = 0.0_DP; edijxc = 0.0_DP
    end if
    if (present(paw_sphere_energies)) then
       paw_sphere_energies(paw_en_exc) = exc
       paw_sphere_energies(paw_en_exc_dc) = exc_dc
       paw_sphere_energies(paw_en_etxc) = etxc
       paw_sphere_energies(paw_en_etxc_dc) = etxc_dc
       paw_sphere_energies(paw_en_dijxc) = edijxc
       paw_sphere_energies(paw_en_exc_core) = exc_core
    end if

    ! Clean up temporary matrices
    do is=pub_num_spins,1,-1
       call sparse_destroy(dij_xc(is))
       call sparse_destroy(dij_hartree(is))
       call sparse_destroy(dij0(is))
    end do
    deallocate(dij_xc,stat=ierr)
    call utils_dealloc_check('paw_nonlocal_energies','dij_xc',ierr)
    deallocate(dij_hartree,stat=ierr)
    call utils_dealloc_check('paw_nonlocal_energies','dij_hartree',ierr)
    deallocate(dij0,stat=ierr)
    call utils_dealloc_check('paw_nonlocal_energies','dij0',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_nonlocal_energies'

  end subroutine paw_nonlocal_energies


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_show_atomblocks(mat,paw_sp)

    use comms, only: comms_reduce, pub_on_root, pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use services, only: services_flush
    use sparse, only: SPAM3, sparse_get_block, sparse_get_par
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert

    ! Arguments
    type(SPAM3), intent(in) :: mat
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    ! Local Variables
    integer :: ierr
    integer :: ipw,npw
    integer :: loc_iat,iat,orig_iat,global_iat,isp
    real(kind=DP), allocatable :: blk(:,:)
    character(20) :: fmt,tmp
    type(PARAL_INFO), pointer :: par

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, mat)
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_show_atomblocks: allocated parallel strategy is &
         &incompatible with paw_sp.')

    allocate(blk(max_paw_proj_tot,max_paw_proj_tot),stat=ierr)
    call utils_alloc_check('paw_show_atomblocks','blk',ierr)

    ! Loop over all species, printing block of first atom of each
    do isp=1,par%num_pspecies
       orig_iat = par%nat + 1
       iat = par%nat + 1
       npw = paw_sp(isp)%npw_tot

       ! Find first instance of this species on local proc
       do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
          iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          if (par%elements_on_proc(loc_iat)%pspecies_number==isp) then
             orig_iat = par%orig_atom(iat)
             exit
          end if
       end do

       ! Find first atom of this type over all procs
       global_iat = orig_iat
       call comms_reduce('MIN',global_iat)
       if (global_iat == par%nat + 1) then
          call utils_abort('ERROR in paw_show_atomblocks: Failed finding &
               &first instance of species')
       end if

       ! Collect to root proc
       blk = 0.0_DP
       if (global_iat == orig_iat) then
          call sparse_get_block(blk(:,:),mat,iat,iat)
       end if
       call comms_reduce('SUM',blk)

       ! Output block
       if (pub_on_root) then
          write(stdout,'(a,i3,a,i5)') 'Species ',isp,': atom ',orig_iat
          write(tmp,'(i5)') npw
          write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f14.10)'
          do ipw=1,npw
             write(stdout,fmt) blk(1:npw,ipw)
          end do
       end if
       call services_flush

    end do
    nullify(par)

    deallocate(blk,stat=ierr)
    call utils_dealloc_check('paw_show_atomblocks','blk',ierr)

  end subroutine paw_show_atomblocks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij0(dij0,paw_sp)

    !=================================================================!
    ! This subroutine calculates the constant terms in the nonlocal   !
    ! energies Dij, which are given in the PAW dataset as Dij0. These !
    ! result from the kinetic energy of the sphere parts of the all-  !
    ! electron wavefunctions, and the interaction of these with the   !
    ! core Hartree potential.                                         !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !   dij0(inout) : Constant contribution D_ij^0 to nonlocal energy !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    !=================================================================!

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_put_block, sparse_get_par
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dij0(pub_num_spins)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    ! Local Variables
    integer :: loc_iat, iat
    integer :: isp
    integer :: is
    integer :: ipwtot,jpwtot
    integer :: ierr
    real(kind=DP), allocatable :: dij0_sp(:,:,:)
    complex(kind=DP), allocatable :: dij0_sp_cmplx(:,:,:)
    type(PARAL_INFO), pointer :: par

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, dij0(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_show_atomblocks: allocated parallel strategy is &
         &incompatible with paw_sp.')

    allocate(dij0_sp(max_paw_proj_tot,max_paw_proj_tot, &
         par%num_pspecies),stat=ierr)
    call utils_alloc_check('paw_dij0','dij0_sp',ierr)

    ! Loop over species, setting up array of dij0 terms
    do isp=1,par%num_pspecies
       if (any(par%elements_on_proc(:)%pspecies_number==isp)) then
          do ipwtot=1,paw_sp(isp)%npw_tot
             do jpwtot=1,paw_sp(isp)%npw_tot
                dij0_sp(ipwtot,jpwtot,isp) = paw_sp(isp)%dij0(ipwtot,jpwtot)
             end do
          end do
       end if
    end do

    if (dij0(1)%iscmplx) then
       allocate(dij0_sp_cmplx(max_paw_proj_tot,max_paw_proj_tot, &
            par%num_pspecies),stat=ierr)
       call utils_alloc_check('paw_dij0','dij0_sp',ierr)
       dij0_sp_cmplx(:,:,:) = cmplx(dij0_sp(:,:,:),kind=DP)
    end if

    ! Loop over atoms
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number
       ! Put block of dij0 into SPAM3 matrix
       do is=1,pub_num_spins
          if (dij0(1)%iscmplx) then
             call sparse_put_block(dij0_sp_cmplx(:,:,isp),dij0(is),iat,iat)
          else
             call sparse_put_block(dij0_sp(:,:,isp),dij0(is),iat,iat)
          end if
       end do
    end do
    nullify(par)

    if (dij0(1)%iscmplx) then
       deallocate(dij0_sp_cmplx,stat=ierr)
       call utils_dealloc_check('paw_dij0','dij0_sp_cmplx',ierr)
    end if

    deallocate(dij0_sp,stat=ierr)
    call utils_dealloc_check('paw_dij0','dij0_sp',ierr)

  end subroutine paw_dij0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_hartree(dijhartree,rho_ij,paw_sp)

    !=================================================================!
    ! This subroutine calculates the Hartree terms in the nonlocal    !
    ! energies Dij, from the tensor e_ijkl calculated from the        !
    ! information in the PAW dataset, in paw_init_species.            !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !   dijhartree(inout) : Hartree contribution to nonlocal energy   !
    !                       term hat{D}_ij.                           !
    !   rho_ij(in)        : Projector density kernel rho^ij           !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    ! First-order capabilities added by Gabriel Constantinescu in May !
    ! 2015.                                                           !
    !=================================================================!

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dijhartree(pub_num_spins)
    type(SPAM3),intent(in) :: rho_ij(pub_num_spins)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    ! Local Variables
    integer :: loc_iat, iat, ierr
    integer :: isp
    integer :: is
    real(kind=DP), allocatable :: rhoij_at(:,:,:), dijh_at(:,:)
    type(PARAL_INFO), pointer :: par

    ! Start Timer
    call timer_clock('paw_dij_hartree',1)

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, dijhartree(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_dij_hartree: allocated parallel strategy is &
         &incompatible with paw_sp.')

    allocate(rhoij_at(max_paw_proj_tot,max_paw_proj_tot, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_hartree','rhoij_at',ierr)
    allocate(dijh_at(max_paw_proj_tot,max_paw_proj_tot),stat=ierr)
    call utils_alloc_check('paw_dij_hartree','dijh_at',ierr)

    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       ! Get rhoij_at for this atom (summed over spins if necessary
       rhoij_at(:,:,:) = 0.0_DP
       do is=1,pub_num_spins
          call sparse_get_block(rhoij_at(:,:,is),rho_ij(is),iat,iat)
       end do
       if (pub_num_spins==2) then
          rhoij_at(:,:,1) = rhoij_at(:,:,1) + rhoij_at(:,:,2)
       end if

       ! Calculate Hartree terms in nonlocal energies for this atom
       dijh_at(:,:) = 0.0_DP
       call paw_dij_hartree_atom(isp,max_paw_proj_tot,rhoij_at,dijh_at, &
            paw_sp)
       do is = 1,pub_num_spins
          call sparse_put_block(dijh_at,dijhartree(is),iat,iat)
       end do

    end do

    nullify(par)
    deallocate(dijh_at,stat=ierr)
    call utils_dealloc_check('paw_dij_hartree','dijh_at',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('paw_dij_hartree','rhoij_at',ierr)

    ! Stop Timer
    call timer_clock('paw_dij_hartree',2)

  end subroutine paw_dij_hartree


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_hartree_atom(isp,npwdim,rhoij,dijh,paw_sp)

    !=================================================================!
    ! This subroutine calculates the Hartree terms in the nonlocal    !
    ! energies Dij for a given atom, from the tensor e_ijkl and the   !
    ! projector density matrix rhoij.                                 !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !   dijh(inout)   : Compensation density contribution to nonlocal !
    !                   energy term \hat{D}_ij.                       !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    !=================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: isp
    integer, intent(in) :: npwdim
    real(kind=DP), intent(in) :: rhoij(npwdim,npwdim)
    real(kind=DP), intent(out) :: dijh(npwdim,npwdim)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    ! Local Variables
    integer :: ipwtot, jpwtot, kpwtot, lpwtot

    dijh(:,:) = 0.0_DP
    do ipwtot=1,paw_sp(isp)%npw_tot
       do jpwtot=1,paw_sp(isp)%npw_tot
          do kpwtot=1,paw_sp(isp)%npw_tot
             do lpwtot=1,paw_sp(isp)%npw_tot
                dijh(ipwtot,jpwtot) = dijh(ipwtot,jpwtot) + &
                     paw_sp(isp)%e_ijkl(ipwtot,jpwtot,kpwtot,lpwtot) * &
                     rhoij(kpwtot,lpwtot)
             end do
          end do
       end do
    end do

  end subroutine paw_dij_hartree_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_xc(dijxc,rho_ij,exc,exc_dc,etxc,etxc_dc,paw_sp,&
       first_order, FO_rho_ij)

    !=================================================================!
    ! This subroutine calculates the XC terms in the nonlocal         !
    ! energies Dij, and also the XC energies and double-counting      !
    ! terms.                                                          !
    !                                                                 !
    ! FO_rho_ij is the first-order projector density kernel, if a     !
    ! first order calculation is required (i.e. if first_order = T and!
    ! 'FO_rho_ij', is present).                                       !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    ! First-order capabilities added by Gabriel Constantinescu in May !
    ! 2015                                                            !
    !=================================================================!

    use comms, only: comms_reduce, pub_my_proc_id, pub_on_root
    use constants, only: max_spins
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_paw_output_detail, pub_aug_funcs_recip,&
         pub_num_spins
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dijxc(pub_num_spins)
    type(SPAM3),intent(in) :: rho_ij(pub_num_spins)
    real(kind=DP),intent(out) :: exc
    real(kind=DP),intent(out) :: exc_dc
    real(kind=DP),intent(out) :: etxc
    real(kind=DP),intent(out) :: etxc_dc
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    logical, intent(in), optional :: first_order
    type(SPAM3), intent(in), optional :: FO_rho_ij(pub_num_spins)

    ! Local Variables
    character(20) :: fmt,tmp
    real(kind=DP), allocatable :: vxc_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: density_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: rhoij_at(:,:,:)
    real(kind=DP), allocatable :: dijxc_at(:,:,:)
    real(kind=DP), allocatable :: dijtxc_at(:,:,:)
    real(kind=DP), allocatable :: pot_work(:,:,:,:)
    real(kind=DP), allocatable :: den_work(:,:,:,:)
    real(kind=DP), allocatable :: inter(:)
    integer :: loc_iat
    integer :: iat
    integer :: ierr
    integer :: is
    integer :: isp
    integer :: igrid
    integer :: ipwtot,jpwtot
    integer :: npts_max
    integer :: nLM_max
    real(kind=DP) :: exc_at
    real(kind=DP) :: exc_dc_at
    real(kind=DP) :: etxc_at
    real(kind=DP) :: etxc_dc_at
    real(kind=DP) :: edijxc_at
    real(kind=DP) :: edijtxc_at
    real(kind=DP) :: total_nhat(max_spins)
    real(kind=DP) :: total_nhat_at(max_spins)

    ! gcc32:
    logical :: loc_first_order
    real(kind=DP), allocatable :: FO_vxc_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: FO_density_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: FO_rhoij_at(:,:,:)
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_dij_xc'

    ! Start Timer
    call timer_clock('paw_dij_xc',1)

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, dijxc(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_dij_xc: allocated parallel strategy is &
         &incompatible with paw_sp.')

    ! check if we are doing first order response
    loc_first_order = .false.
    if (present(first_order).and.present(FO_rho_ij)) then
       loc_first_order = first_order
    end if

    ! Find array sizes
    npts_max = 0
    nLM_max = 0
    do isp=1,par%num_pspecies
       if (.not.any(par%elements_on_proc(:)%pspecies_number==isp)) cycle
       igrid = paw_sp(isp)%phi_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       igrid = paw_sp(isp)%shape_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       igrid = paw_sp(isp)%core_den_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       nLM_max = max(nLM_max,(paw_sp(isp)%lmax+1)**2)
    end do

    ! Allocate temporary arrays
    allocate(vxc_rad_LM(npts_max,pub_num_spins,nLM_max,2),stat=ierr)
    call utils_alloc_check('paw_dij_xc','vxc_rad_LM',ierr)
    allocate(density_rad_LM(npts_max,pub_num_spins,nLM_max,3),stat=ierr)
    call utils_alloc_check('paw_dij_xc','density_rad_LM',ierr)
    allocate(rhoij_at(max_paw_proj_tot,max_paw_proj_tot, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc','rhoij_at',ierr)
    allocate(dijxc_at(max_paw_proj_tot,max_paw_proj_tot, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc','dijxc_at',ierr)
    allocate(dijtxc_at(max_paw_proj_tot,max_paw_proj_tot, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc','dijtxc_at',ierr)
    allocate(pot_work(npts_max,pub_num_spins,pub_num_spins,6), &
         stat=ierr)
    call utils_alloc_check('paw_dij_xc','pot_work',ierr)
    allocate(den_work(npts_max,pub_num_spins,pub_num_spins,4), &
         stat=ierr)
    call utils_alloc_check('paw_dij_xc','den_work',ierr)
    allocate(inter(npts_max),stat=ierr)
    call utils_alloc_check('paw_dij_xc','inter',ierr)

    if (loc_first_order) then ! gcc32: first-order calculation
       allocate(FO_vxc_rad_LM(npts_max,pub_num_spins,nLM_max,2),stat=ierr)
       call utils_alloc_check('paw_dij_xc','FO_vxc_rad_LM',ierr)
       allocate(FO_density_rad_LM(npts_max,pub_num_spins,nLM_max,3),stat=ierr)
       call utils_alloc_check('paw_dij_xc','FO_density_rad_LM',ierr)
       allocate(FO_rhoij_at(max_paw_proj_tot,max_paw_proj_tot,pub_num_spins), &
            stat=ierr)
       call utils_alloc_check('paw_dij_xc','FO_rhoij_at',ierr)
    end if

    ! Initialisations
    exc = 0.0_DP
    exc_dc = 0.0_DP
    etxc = 0.0_DP
    etxc_dc = 0.0_DP
    total_nhat(:) = 0.0_DP
    total_nhat_at(:) = 0.0_DP

    ! Loop over atoms
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       ! Get block of projector density matrix for this atom
       do is = 1,pub_num_spins
          call sparse_get_block(rhoij_at(:,:,is),rho_ij(is),iat,iat)

          if (loc_first_order) then
             call sparse_get_block(FO_rhoij_at(:,:,is),FO_rho_ij(is),iat,iat)
          end if
       end do  ! is

       ! Calculate exchange-correlation energy and potential for this atom
       if (loc_first_order) then

          FO_vxc_rad_LM = 0.0_DP
          FO_density_rad_LM = 0.0_DP

          call paw_dij_xc_atom(isp,pub_num_spins,npts_max,nLM_max, &
               max_paw_proj_tot,rhoij_at,dijxc_at,dijtxc_at, &
               exc_at,etxc_at,exc_dc_at,etxc_dc_at,total_nhat_at, &
               density_rad_LM,vxc_rad_LM,den_work,pot_work,inter,paw_sp, &
              .false.,FO_rhoij_at,FO_density_rad_LM, &
               FO_vxc_rad_LM,.true.)
          ! the weight of the perturbation is implicit in FO_rhoij, since
          ! we do not do any physical displacements in this subroutine
       else

             call paw_dij_xc_atom(isp,pub_num_spins,npts_max,nLM_max, &
                  max_paw_proj_tot,rhoij_at,dijxc_at,dijtxc_at, &
                  exc_at,etxc_at,exc_dc_at,etxc_dc_at,total_nhat_at, &
                  density_rad_LM,vxc_rad_LM,den_work,pot_work,inter,paw_sp, &
                  .false.)
       end if

       ! Add up contributions to energy from this atom
       exc = exc + exc_at
       etxc = etxc + etxc_at
       exc_dc = exc_dc + exc_dc_at
       etxc_dc = etxc_dc + etxc_dc_at
       total_nhat(:) = total_nhat(:) + total_nhat_at(:)

       ! Create total dijxc and add it to matrix (then undo total)
       do is = 1,pub_num_spins
          dijxc_at(:,:,is) = dijxc_at(:,:,is) - dijtxc_at(:,:,is)
          call sparse_put_block(dijxc_at(:,:,is),dijxc(is),iat,iat)
          dijxc_at(:,:,is) = dijxc_at(:,:,is) + dijtxc_at(:,:,is)
       end do  ! is

       ! Add up AE dijxc energy \sum_ij \rho_ij D^xc_ij and
       ! and PS dijtxc energy \sum_ij \rho_ij \tilde{D}^xc_ij
       edijxc_at = 0.0_DP
       edijtxc_at = 0.0_DP
       do is=1,pub_num_spins
          do ipwtot=1,paw_sp(isp)%npw_tot
             do jpwtot=1,paw_sp(isp)%npw_tot
                edijxc_at = edijxc_at + rhoij_at(ipwtot,jpwtot,is) * &
                     dijxc_at(ipwtot,jpwtot,is)
                edijtxc_at = edijtxc_at + rhoij_at(ipwtot,jpwtot,is) * &
                     dijtxc_at(ipwtot,jpwtot,is)
             end do
          end do
       end do

       ! Consistency Check
       ! gcc32: also unnecessary for loc_first_order
       ! tjz21: removed loc_is_tddft_here
       if (.not. (loc_first_order) ) then
          if ((abs(edijxc_at-exc_dc_at)>0.0000001_DP*abs(exc_dc_at)) .and. &
               (exc_dc_at>0.0000001_DP)) then
             if (pub_on_root) then
                write(stdout,'(a,i6,a)') 'For atom ',iat,':'
                write(stdout,'(a,f20.12)') 'edijxc_at: ',edijxc_at
                write(stdout,'(a,f20.12)') 'exc_dc_at: ',exc_dc_at
             end if
             call utils_abort('Error in paw_dij_xc: Consistency check &
                  &failed: edijxc_at /= exc_dc_at')
          end if
          if ((abs(edijtxc_at-etxc_dc_at)>0.0000001_DP*abs(etxc_dc_at)).and. &
               (etxc_dc_at>0.0000001_DP)) then
             if (pub_on_root) then
                write(stdout,'(a,i6,a)') 'For atom ',iat,':'
                write(stdout,'(a,f20.12)') 'edijtxc_at: ',edijtxc_at
                write(stdout,'(a,f20.12)') 'etxc_dc_at: ',etxc_dc_at
             end if
             call utils_abort('Error in paw_dij_xc: Consistency check &
                  &failed: edijtxc_at /= etxc_dc_at')
          end if
       end if

    end do  ! loc_iat

    ! Sum energies across all procs
    call comms_reduce('SUM',exc)
    call comms_reduce('SUM',exc_dc)
    call comms_reduce('SUM',etxc)
    call comms_reduce('SUM',etxc_dc)

    ! Check and show total compensation charge on radial grid
    call comms_reduce('SUM',total_nhat(:))
    if (pub_on_root.and.(pub_paw_output_detail>=VERBOSE).and. &
         (.not.pub_aug_funcs_recip)) then
       write(tmp,'(i5)') pub_num_spins
       write(fmt,'(3a)') '(a,',trim(adjustl(tmp)),'f14.8)'
       write(stdout,fmt) 'Total Compensation Charge on Radial Grid : ',&
            total_nhat(1:pub_num_spins)
    end if

    ! Deallocate temporary arrays
    deallocate(inter,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','inter',ierr)
    deallocate(den_work,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','den_work',ierr)
    deallocate(pot_work,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','pot_work',ierr)
    deallocate(dijtxc_at,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','dijtxc_at',ierr)
    deallocate(dijxc_at,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','dijxc_at',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','rhoij_at',ierr)
    deallocate(density_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','density_rad_LM',ierr)
    deallocate(vxc_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','vxc_rad_LM',ierr)


    if (loc_first_order) then ! gcc32: first-order calculation
       deallocate(FO_vxc_rad_LM,stat=ierr)
       call utils_dealloc_check('paw_dij_xc','FO_vxc_rad_LM',ierr)
       deallocate(FO_density_rad_LM,stat=ierr)
       call utils_dealloc_check('paw_dij_xc','FO_density_rad_LM',ierr)
       deallocate(FO_rhoij_at,stat=ierr)
       call utils_dealloc_check('paw_dij_xc','FO_rhoij_at',ierr)
    end if
    nullify(par)

    ! Stop Timer
    call timer_clock('paw_dij_xc',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_dij_xc'

  end subroutine paw_dij_xc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine paw_dij_xc_triplet(dijxc,rho_ij,rho_ij_0,paw_sp)

    !=================================================================!
    ! This subroutine calculates the XC terms in the nonlocal         !
    ! energies Dij in the case that this is a TDDFT triplet           !
    ! calculation for spin-degenerate systems. Since the exchange-    !
    ! correlation potential for triplet states is given by [fxc(up,up)!
    ! -fxc(up,down)]*rho^{1}, the most convenient way of computing    !
    ! the PAW contribution is to call spin-polarised xc_routines even !
    ! for the spin-degenerate system and to use a finite difference   !
    ! approximation. This routine is adapted from paw_dij_xc.         !
    !-----------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff in April 2016.                       !
    !=================================================================!

    use comms, only: comms_reduce, pub_my_proc_id
    use constants, only: max_spins
    use rundat, only: pub_num_spins
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block,sparse_get_par,&
         sparse_create,sparse_destroy,sparse_copy,sparse_scale,sparse_axpy
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dijxc(pub_num_spins)
    type(SPAM3),intent(in) :: rho_ij(pub_num_spins) ! perturbed rho_ij
    type(SPAM3),intent(in) :: rho_ij_0(pub_num_spins) ! ground state rho_ij
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    ! Local Variables
    integer, parameter :: eff_num_spins=2 ! effectively we are treating this
        ! routine as if we have pub_num_spins=2
    type(SPAM3), allocatable :: eff_rho_ij(:)
    type(SPAM3), allocatable :: eff_dijxc(:)
    real(kind=DP), allocatable :: vxc_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: density_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: rhoij_at(:,:,:)
    real(kind=DP), allocatable :: dijxc_at(:,:,:)
    real(kind=DP), allocatable :: dijtxc_at(:,:,:)
    real(kind=DP), allocatable :: pot_work(:,:,:,:)
    real(kind=DP), allocatable :: den_work(:,:,:,:)
    real(kind=DP), allocatable :: inter(:)
    integer :: loc_iat
    integer :: iat
    integer :: ierr
    integer :: is
    integer :: isp
    integer :: igrid
    integer :: npts_max
    integer :: nLM_max
    real(kind=DP) :: exc_at
    real(kind=DP) :: exc_dc_at
    real(kind=DP) :: etxc_at
    real(kind=DP) :: etxc_dc_at
    real(kind=DP) :: total_nhat_at(eff_num_spins) ! should be max_spins. Is this right?
    type(PARAL_INFO), pointer :: par

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, dijxc(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_dij_xc_triplet: allocated parallel strategy is &
         &incompatible with paw_sp.')

    ! Find array sizes
    npts_max = 0
    nLM_max = 0
    do isp=1,par%num_pspecies
       if (.not.any(par%elements_on_proc(:)%pspecies_number==isp)) cycle
       igrid = paw_sp(isp)%phi_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       igrid = paw_sp(isp)%shape_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       igrid = paw_sp(isp)%core_den_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       nLM_max = max(nLM_max,(paw_sp(isp)%lmax+1)**2)
    end do

    ! Allocate temporary arrays
    allocate(eff_dijxc(eff_num_spins), stat=ierr)
    call utils_alloc_check('paw_dij_xc_triplet','eff_dijxc',ierr)
    allocate(eff_rho_ij(eff_num_spins), stat=ierr)
    call utils_alloc_check('paw_dij_xc_triplet','eff_rho_ij',ierr)
    allocate(vxc_rad_LM(npts_max,eff_num_spins,nLM_max,2),stat=ierr)
    call utils_alloc_check('paw_dij_xc_triplet','vxc_rad_LM',ierr)
    allocate(density_rad_LM(npts_max,eff_num_spins,nLM_max,3),stat=ierr)
    call utils_alloc_check('paw_dij_xc','density_rad_LM',ierr)
    allocate(rhoij_at(max_paw_proj_tot,max_paw_proj_tot, &
         eff_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc','rhoij_at',ierr)
    allocate(dijxc_at(max_paw_proj_tot,max_paw_proj_tot, &
         eff_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc','dijxc_at',ierr)
    allocate(dijtxc_at(max_paw_proj_tot,max_paw_proj_tot, &
         eff_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc','dijtxc_at',ierr)
    allocate(pot_work(npts_max,eff_num_spins,eff_num_spins,6), &
         stat=ierr)
    call utils_alloc_check('paw_dij_xc','pot_work',ierr)
    allocate(den_work(npts_max,eff_num_spins,eff_num_spins,4), &
         stat=ierr)
    call utils_alloc_check('paw_dij_xc','den_work',ierr)
    allocate(inter(npts_max),stat=ierr)
    call utils_alloc_check('paw_dij_xc','inter',ierr)

    ! create sparse matrix
    do is=1, eff_num_spins
       call sparse_create(eff_rho_ij(is),rho_ij(1))
       call sparse_create(eff_dijxc(is),dijxc(1))
    enddo

    ! fill eff_rho_ij
    ! eff_rho_ij(1)=rho_ij_0+0.5_DP*epsilon*rho_ij_1
    ! eff_rho_ij(2)=rho_ij_0
    call sparse_copy(eff_rho_ij(1),rho_ij(1))
    call sparse_scale(eff_rho_ij(1),0.5_DP)
    call sparse_copy(eff_rho_ij(2),rho_ij_0(1))

    ! Loop over atoms
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       ! Get block of projector density matrix for this atom
       do is = 1,eff_num_spins
          call sparse_get_block(rhoij_at(:,:,is),eff_rho_ij(is),iat,iat)
       end do  ! is

       call paw_dij_xc_atom(isp,eff_num_spins,npts_max,nLM_max, &
            max_paw_proj_tot,rhoij_at,dijxc_at,dijtxc_at, &
            exc_at,etxc_at,exc_dc_at,etxc_dc_at,total_nhat_at, &
            density_rad_LM,vxc_rad_LM,den_work,pot_work,inter,paw_sp, &
            .false.)

       ! Create total dijxc and add it to matrix (then undo total)
       do is = 1,eff_num_spins
          dijxc_at(:,:,is) = dijxc_at(:,:,is) - dijtxc_at(:,:,is)
          call sparse_put_block(dijxc_at(:,:,is),eff_dijxc(is),iat,iat)
          dijxc_at(:,:,is) = dijxc_at(:,:,is) + dijtxc_at(:,:,is)
       end do  ! is

    enddo ! iat

    ! create dijxc from eff_dijxc
    call sparse_copy(dijxc(1),eff_dijxc(1))
    call sparse_axpy(dijxc(1),eff_dijxc(2),-1.0_DP)

    ! deallocate data structures
    do is=1,eff_num_spins
       call sparse_destroy(eff_rho_ij(is))
    enddo
    nullify(par)
    deallocate(eff_rho_ij,stat=ierr)
    call utils_dealloc_check('paw_dij_xc_triplet','eff_rho_ij',ierr)
    deallocate(inter,stat=ierr)
    call utils_dealloc_check('paw_dij_xc_triplet','inter',ierr)
    deallocate(den_work,stat=ierr)
    call utils_dealloc_check('paw_dij_xc_triplet','den_work',ierr)
    deallocate(pot_work,stat=ierr)
    call utils_dealloc_check('paw_dij_xc_triplet','pot_work',ierr)
    deallocate(dijtxc_at,stat=ierr)
    call utils_dealloc_check('paw_dij_xc_triplet','dijtxc_at',ierr)
    deallocate(dijxc_at,stat=ierr)
    call utils_dealloc_check('paw_dij_xc_triplet','dijxc_at',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('paw_dij_xc_triplet','rhoij_at',ierr)
    deallocate(density_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_dij_xc_triplet','density_rad_LM',ierr)
    deallocate(vxc_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_dij_xc_triplet','vxc_rad_LM',ierr)

  end subroutine paw_dij_xc_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine paw_dij_xc_tddft(dijxc,rho_1_ij,rho_0_ij,paw_sp,&
       epsilon_value,triplet)

    !================================================================!
    ! This subroutine calculates the XC terms in the nonlocal        !
    ! energies Dij for the TDDFT case by using a finite difference   !
    ! approximation. This allows the reuse of the standard           !
    ! paw_dij_xc routine used in ground state dft                    !
    !----------------------------------------------------------------!
    ! Written by Tim Zuehlsdorff on 10/01/14.                        !
    !================================================================!

    use comms, only: comms_reduce
    use constants, only: stdout
    use rundat, only: pub_num_spins, pub_debug_on_root
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block, sparse_create,&
        sparse_destroy, sparse_copy, sparse_scale, sparse_axpy
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dijxc(pub_num_spins) ! correction term dij_xc
    type(SPAM3),intent(in) :: rho_1_ij(pub_num_spins) ! transition density rho_ij
    type(SPAM3),intent(in) :: rho_0_ij(pub_num_spins) ! ground state rho_ij
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)
    real(kind=DP), intent(in) :: epsilon_value ! finite difference parameter
    logical, intent(in), optional :: triplet   ! flag determining if triplet calculation

    ! local variables
    logical :: loc_triplet
    integer :: ierr, is
    type(SPAM3), allocatable, dimension(:) :: rho_ij_effective
    type(SPAM3), allocatable, dimension(:)  :: dij_xc_1, dij_xc_2
    real(kind=DP) :: exc ! dummy variables. The sphere energies calculated
    real(kind=DP) :: exc_dc ! through paw_dij_xc are not actually used in
    real(kind=DP) :: etxc   ! lr_tddft
    real(kind=DP) :: etxc_dc

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_dij_xc_tddft'

    ! Start Timer
    call timer_clock('paw_dij_xc_tddft',1)

    loc_triplet=.false.
    if(present(triplet)) loc_triplet=triplet

    ! allocate appropriate temporary storage space:
    allocate(rho_ij_effective(pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc_tddft','rho_ij_effective',ierr)
    allocate(dij_xc_1(pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc_tddft','dij_xc_1',ierr)
    allocate(dij_xc_2(pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc_tddft','dij_xc_2',ierr)

    do is=1, pub_num_spins
       call sparse_create(rho_ij_effective(is), rho_1_ij(is))
       call sparse_create(dij_xc_1(is),dijxc(is))
       call sparse_create(dij_xc_2(is),dijxc(is))
    enddo

    if(pub_num_spins==1) then
       ! create first effective perturbed density kernel inside the spheres through
       ! the finite difference approximation
       call sparse_copy(rho_ij_effective(1), rho_0_ij(1))
       call sparse_scale(rho_ij_effective(1), 2.0_DP)
       call sparse_axpy(rho_ij_effective(1), rho_1_ij(1), epsilon_value)

       ! calculate dij_xc through paw_dij_xc
       if(loc_triplet) then
          call paw_dij_xc_triplet(dij_xc_1,rho_ij_effective,rho_0_ij,paw_sp)
       else
          call paw_dij_xc(dij_xc_1,rho_ij_effective,exc,exc_dc,etxc,etxc_dc,&
              paw_sp)

       endif
       ! construct new effective rho as a second value for the central difference
       ! approximation
       call sparse_copy(rho_ij_effective(1), rho_0_ij(1))
       call sparse_scale(rho_ij_effective(1), 2.0_DP)
       call sparse_axpy(rho_ij_effective(1), rho_1_ij(1), -epsilon_value)

       ! calculate second part of dij
       if(loc_triplet) then
          call paw_dij_xc_triplet(dij_xc_2,rho_ij_effective,rho_0_ij,paw_sp)
       else
          call paw_dij_xc(dij_xc_2,rho_ij_effective,exc,exc_dc,etxc,etxc_dc,&
             paw_sp)
       endif

       ! find correct effective value through finite difference approximation
       call sparse_copy(dijxc(1),dij_xc_1(1))
       call sparse_axpy(dijxc(1),dij_xc_2(1), -1.0_DP)
       call sparse_scale(dijxc(1),1.0_DP/(2.0_DP*epsilon_value))
       ! scale for spin degeneracy:
       call sparse_scale(dijxc(1),2.0_DP)
    else  ! spin polarised case
       ! if this is a spin polarised case, just call normal paw dij_xc module but with
       ! rho_ij_effective treated in both spin channels
       call sparse_copy(rho_ij_effective(1),rho_0_ij(1))
       call sparse_copy(rho_ij_effective(2),rho_0_ij(2))

       do is=1, pub_num_spins
          call sparse_axpy(rho_ij_effective(is),rho_1_ij(is), 1.0_DP*epsilon_value)
       enddo

       call paw_dij_xc(dij_xc_1,rho_ij_effective,exc,exc_dc,etxc,etxc_dc,&
            paw_sp)

       call sparse_copy(rho_ij_effective(1),rho_0_ij(1))
       call sparse_copy(rho_ij_effective(2),rho_0_ij(2))
       do is=1, pub_num_spins
          call sparse_axpy(rho_ij_effective(is),rho_1_ij(is), -1.0_DP*epsilon_value)
       enddo

       call paw_dij_xc(dij_xc_2,rho_ij_effective,exc,exc_dc,etxc,etxc_dc,&
            paw_sp)

       do is=1, pub_num_spins
          call sparse_copy(dijxc(is),dij_xc_1(is))
          call sparse_axpy(dijxc(is),dij_xc_2(is),-1.0_DP)
          call sparse_scale(dijxc(is),1.0_DP/(2.0_DP*epsilon_value))
       enddo
    endif

    ! deallocate data storage
    do is=1, pub_num_spins
       call sparse_destroy(dij_xc_1(is))
       call sparse_destroy(dij_xc_2(is))
       call sparse_destroy(rho_ij_effective(is))
    enddo
    deallocate(dij_xc_1, stat=ierr)
    call utils_dealloc_check('paw_dij_xc_tddft','dij_xc_1',ierr)
    deallocate(dij_xc_2, stat=ierr)
    call utils_dealloc_check('paw_dij_xc_tddft','dij_xc_2',ierr)
    deallocate(rho_ij_effective, stat=ierr)
    call utils_dealloc_check('paw_dij_xc_tddft','rho_ij_effective',ierr)

    call timer_clock('paw_dij_xc_tddft',2)

  end subroutine paw_dij_xc_tddft


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_xc_atom(isp,nspins,npts_max,nLM_max,npw, &    ! in
       rhoij,dijxc,dijtxc,exc,etxc,exc_dc,etxc_dc,total_nhat_at, & ! out
       density_rad_LM,vxc_rad_LM,den_work,pot_work,inter,paw_sp,prt,& !workspace
       FO_rhoij,FO_density_rad_LM,FO_vxc_rad_LM,first_order) ! first-order param

    !=================================================================!
    ! This subroutine calculates the xc terms in the nonlocal         !
    ! energies Dij for a given atom, using the projector density      !
    ! matrix rhoij.                                                   !
    !                                                                 !
    ! For first-order terms, one also needs the first-order rhoij     !
    ! ('FO_rhoij'), as well as the direction ('cart') of perturbation !
    ! for a perturbed atom of species isp                             !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    ! First-order capabilities added by Gabriel Constantinescu in May !
    ! 2015                                                            !
    !=================================================================!

    use constants, only: max_spins, PI, stdout
    use gaunt_coeff, only: realgaunt
    use paw_xc, only: paw_xc_pot_rad_LM, paw_xc_exc_dc
    use rundat, only: pub_nhat_in_xc, pub_edft_grand_canonical
    use services, only: services_radial_integral_rmax
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: isp
    integer, intent(in) :: npw
    integer, intent(in) :: nspins
    integer, intent(in) :: npts_max, nLM_max
    real(kind=DP), intent(in) :: rhoij(npw,npw,nspins)
    real(kind=DP), intent(out) :: dijxc(npw,npw,nspins)
    real(kind=DP), intent(out) :: dijtxc(npw,npw,nspins)
    real(kind=DP), intent(out) :: exc, etxc
    real(kind=DP), intent(out) :: exc_dc, etxc_dc
    real(kind=DP), intent(out) :: density_rad_LM(npts_max,nspins,nLM_max,3)
    real(kind=DP), intent(out) :: vxc_rad_LM(npts_max,nspins,nLM_max,2)
    real(kind=DP), intent(out) :: den_work(npts_max,nspins,nspins,4)
    real(kind=DP), intent(out) :: pot_work(npts_max,nspins,nspins,6)
    real(kind=DP), intent(out) :: inter(npts_max)
    real(kind=DP), intent(out) :: total_nhat_at(max_spins)
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)
    logical, optional, intent(in) :: prt
    ! Optional Arguments for first order calculation
    real(kind=DP), intent(in), optional :: FO_rhoij(npw,npw,nspins)
    real(kind=DP), intent(out), optional :: FO_density_rad_LM(npts_max, &
         nspins,nLM_max,3)
    real(kind=DP), intent(out), optional :: FO_vxc_rad_LM(npts_max,nspins, &
         nLM_max,2)
    logical, intent(in), optional :: first_order


    ! Local Variables
    integer :: ipwtot, jpwtot
    integer :: igrid, is
    integer :: lup, mup
    integer :: ipw,jpw
    integer :: li,mi,lj,mj
    integer :: npts_shp
    integer :: npts_phi
    integer :: npts_den
    integer :: nLM
    integer :: iLM
    integer :: lmax
    logical :: loc_prt
    real(kind=DP) :: rgij
    real(kind=DP) :: total_nhat_at_targ
    logical :: loc_first_order

    ! gcc32:
    real(kind=DP), allocatable :: pert_density_mi_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: pert_density_pl_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: temp_vxc_rad_LM(:,:,:)
    real(kind=DP) :: delta
    integer :: iter,ierr

    ! check if we are doing a first order response
    loc_first_order = .false.
    if ( present(first_order) ) then
       loc_first_order = first_order
    end if

    if (loc_first_order) then
       allocate(pert_density_mi_rad_LM(npts_max,nspins,nLM_max,3),stat=ierr)
       call utils_alloc_check('paw_dij_xc_atom','pert_density_mi_rad_LM',ierr)
       allocate(pert_density_pl_rad_LM(npts_max,nspins,nLM_max,3),stat=ierr)
       call utils_alloc_check('paw_dij_xc_atom','pert_density_pl_rad_LM',ierr)
       allocate(temp_vxc_rad_LM(npts_max,nspins,nLM_max),stat=ierr)
       call utils_alloc_check('paw_dij_xc_atom','temp_vxc_rad_LM',ierr)

       delta = 0.001_DP
    end if

    loc_prt = .false.
    if (present(prt)) loc_prt = prt

    ! Find numbers of grid points and angular momenta for this atom
    igrid = paw_sp(isp)%phi_grid
    npts_phi = paw_sp(isp)%grid(igrid)%npt
    igrid = paw_sp(isp)%shape_grid
    npts_shp = paw_sp(isp)%grid(igrid)%npt
    igrid = paw_sp(isp)%core_den_grid
    npts_den = paw_sp(isp)%grid(igrid)%npt
    npts_den = max(npts_den,npts_shp,npts_phi)
    nLM = (paw_sp(isp)%lmax + 1)**2
    lmax = paw_sp(isp)%lmax

    ! Reset dijxc blocks for this atom
    dijxc(:,:,:) = 0.0_DP
    dijtxc(:,:,:) = 0.0_DP
    exc_dc = 0.0_DP
    etxc_dc = 0.0_DP

    ! Create (n^1 + n_c)(r) in density(1)
    ! Create (\tilde{n}^1 + \tilde{n}_c)(r) in density(2)
    ! Create \hat{n}  in density(3)
    call paw_ae_density_rad_LM(density_rad_LM(:,:,:,1),isp,rhoij, &
         npw,npts_max,nspins,nLM_max,paw_sp,den_work)
    call paw_ps_density_rad_LM(density_rad_LM(:,:,:,2),isp,rhoij, &
         npw,npts_max,nspins,nLM_max,paw_sp,den_work)
    call paw_aug_density_rad_LM(density_rad_LM(:,:,:,3),isp,rhoij, &
         npw,npts_max,nspins,nLM_max,paw_sp,den_work)

    ! Calculate total compensation density and compare to target
    do is=1,nspins
       den_work(1:npts_den,1,1,1) = density_rad_LM(1:npts_den,is,1,3) * &
            paw_sp(isp)%grid(igrid)%r(1:npts_den)**2 * inv_sqrt_4pi
       total_nhat_at(is) = services_radial_integral_rmax(npts_den, &
            paw_sp(isp)%grid(igrid)%rab,paw_sp(isp)%grid(igrid)%r, &
            paw_sp(isp)%rcut,den_work(:,1,1,1),inter) * 4.0_DP*PI
       total_nhat_at_targ = 0.0_DP
       do ipwtot=1,paw_sp(isp)%npw_tot
          ipw = paw_sp(isp)%ipw_tot(ipwtot)
          li = paw_sp(isp)%l_pw_tot(ipwtot)
          mi = paw_sp(isp)%m_pw_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)

             total_nhat_at_targ = total_nhat_at_targ + sqrt(4.0_DP*PI) * &
                  realgaunt(0,0,li,mi,lj,mj) * rhoij(ipwtot,jpwtot,is) * &
                  paw_sp(isp)%aug_nLij(ipw,jpw,0)
          end do  ! jpwtot
       end do  ! ipwtot

       ! Check total compensation density does not differ by more than 0.1%
       ! from its intended value
       if (.not.pub_edft_grand_canonical) then
          if (abs(total_nhat_at_targ-total_nhat_at(is)) > &
               0.001_DP*abs(total_nhat_at_targ)) then
               call utils_abort('Error in paw_dij_xc_atom: Total Compensation &
                    &charge on radial grid does not match target')
          end if
       end if
    end do  ! is

    ! Print AE density
    if (loc_prt) call print_rad(density_rad_LM(:,:,:,1), &
         paw_sp(isp)%grid(igrid)%r,npts_den,nspins,nLM,'n^1+n_c',55)

    ! Create d n^1 (r) / d_epsilon in FO_density(1)
    ! Create d \tilde{n}^1 (r) / d_epsilon in FO_density(2)
    ! Create d \hat{n} / d_epsilon in FO_density(3)
    if (loc_first_order) then
       call paw_ae_density_rad_LM(FO_density_rad_LM(:,:,:,1),isp,FO_rhoij, &
            npw,npts_max,nspins,nLM_max,paw_sp,den_work,.true.)
       call paw_ps_density_rad_LM(FO_density_rad_LM(:,:,:,2),isp,FO_rhoij, &
            npw,npts_max,nspins,nLM_max,paw_sp,den_work,.true.)
       call paw_aug_density_rad_LM(FO_density_rad_LM(:,:,:,3),isp,FO_rhoij, &
            npw,npts_max,nspins,nLM_max,paw_sp,den_work)

       ! calculate total perturbed density with + delta and - delta
       do iLM=1,nLM_max
          do is=1,nspins
             do iter=1,3
                pert_density_pl_rad_LM(:,is,iLM,iter) = density_rad_LM(:,is, &
                     iLM,iter) + delta * FO_density_rad_LM(:,is,iLM,iter)
                pert_density_mi_rad_LM(:,is,iLM,iter) = density_rad_LM(:,is, &
                     iLM,iter) - delta * FO_density_rad_LM(:,is,iLM,iter)
             end do
          end do
       end do

    end if ! if loc_first_order

    ! Calculate v_xc[density(1)]
    call paw_xc_pot_rad_LM(vxc_rad_LM(:,:,:,1),exc, &
         density_rad_LM(:,:,:,1), &
         paw_sp(isp)%grid(igrid)%r(:),paw_sp(isp)%grid(igrid)%rab(:), &
         paw_sp(isp)%rcut,npts_phi,npts_max,nspins,nLM_max,lmax, &
         den_work,pot_work,inter(:),loc_prt)


    if (loc_first_order) then

       ! Calculate v_xc[density(1)+delta*FO_density(1)]
       call paw_xc_pot_rad_LM(temp_vxc_rad_LM(:,:,:),exc, &
            pert_density_pl_rad_LM(:,:,:,1), &
            paw_sp(isp)%grid(igrid)%r(:),paw_sp(isp)%grid(igrid)%rab(:), &
            paw_sp(isp)%rcut,npts_phi,npts_max,nspins,nLM_max,lmax, &
            den_work,pot_work,inter(:),loc_prt)

       FO_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max,1) = &
            temp_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max) / (2.0_DP*delta)

       ! Calculate v_xc[density(1)-delta*FO_density(1)]
       call paw_xc_pot_rad_LM(temp_vxc_rad_LM(:,:,:),exc, &
            pert_density_mi_rad_LM(:,:,:,1), &
            paw_sp(isp)%grid(igrid)%r(:),paw_sp(isp)%grid(igrid)%rab(:), &
            paw_sp(isp)%rcut,npts_phi,npts_max,nspins,nLM_max,lmax, &
            den_work,pot_work,inter(:),loc_prt)

       FO_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max,1) = &
            FO_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max,1) - &
            temp_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max) / (2.0_DP*delta)
    end if


    ! Subtract off AE core density and calculate double-counting term
    ! not required if this is a tddft calculation or first-order
    ! tjz21: Removed lr_tddft from here
    if (.not. (loc_first_order)) then
       do is=1,nspins
          density_rad_LM(1:npts_den,is,1,1) =density_rad_LM(1:npts_den,is,1,1) &
               - paw_sp(isp)%core_den_rad(1:npts_den) * sqrt_4pi &
               / real(nspins,kind=DP)
       end do

       ! Calculate Double Counting correction to AE XC energy
       call paw_xc_exc_dc(exc_dc,vxc_rad_LM(:,:,:,1), &
            density_rad_LM(:,:,:,1),paw_sp(isp)%grid(igrid)%r(:), &
            paw_sp(isp)%grid(igrid)%rab(:),paw_sp(isp)%rcut, &
            npts_phi,npts_max,nspins,nLM,nLM_max,den_work(:,1,1,1),inter(:))
    endif

    ! Print AE potential
    if (loc_prt) call print_rad(vxc_rad_LM(:,:,:,1),paw_sp(isp)%grid(igrid)%r, &
         npts_den,nspins,nLM,'vxc(n^1+n_c)',56)

    ! Add compensation density nhat to valence density \tilde{n}^1
    ! (but only if we are using the nhat density in the xc terms)
    if (pub_nhat_in_xc) then
       density_rad_LM(1:npts_den,1:nspins,1:nLM,2) = &
            density_rad_LM(1:npts_den,1:nspins,1:nLM,2) + &
            density_rad_LM(1:npts_den,1:nspins,1:nLM,3)

       if (loc_first_order) then
          pert_density_pl_rad_LM(1:npts_den,1:nspins,1:nLM,2) = &
               pert_density_pl_rad_LM(1:npts_den,1:nspins,1:nLM,2) + &
               pert_density_pl_rad_LM(1:npts_den,1:nspins,1:nLM,3)

          pert_density_mi_rad_LM(1:npts_den,1:nspins,1:nLM,2) = &
               pert_density_mi_rad_LM(1:npts_den,1:nspins,1:nLM,2) + &
               pert_density_mi_rad_LM(1:npts_den,1:nspins,1:nLM,3)
       end if

    end if

    ! Print soft PS density
    if (loc_prt) call print_rad(density_rad_LM(:,:,:,2), &
         paw_sp(isp)%grid(igrid)%r,npts_den,nspins,nLM,'tn^1+nhat+tn_c',65)

    ! Calculate v_xc[density(2)]
    call paw_xc_pot_rad_LM(vxc_rad_LM(:,:,:,2),etxc, &
         density_rad_LM(:,:,:,2), &
         paw_sp(isp)%grid(igrid)%r(:),paw_sp(isp)%grid(igrid)%rab(:), &
         paw_sp(isp)%rcut,npts_phi,npts_max,nspins,nLM_max,lmax, &
         den_work,pot_work,inter(:))

    if (loc_first_order) then

       ! Calculate v_xc[density(2)+delta*FO_density(2)]
       call paw_xc_pot_rad_LM(temp_vxc_rad_LM(:,:,:),exc, &
            pert_density_pl_rad_LM(:,:,:,2), &
            paw_sp(isp)%grid(igrid)%r(:),paw_sp(isp)%grid(igrid)%rab(:), &
            paw_sp(isp)%rcut,npts_phi,npts_max,nspins,nLM_max,lmax, &
            den_work,pot_work,inter(:),loc_prt)

       FO_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max,2) = &
            temp_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max) / (2.0_DP*delta)

       ! Calculate v_xc[density(2)-delta*FO_density(2)]
       call paw_xc_pot_rad_LM(temp_vxc_rad_LM(:,:,:),exc, &
            pert_density_mi_rad_LM(:,:,:,2), &
            paw_sp(isp)%grid(igrid)%r(:),paw_sp(isp)%grid(igrid)%rab(:), &
            paw_sp(isp)%rcut,npts_phi,npts_max,nspins,nLM_max,lmax, &
            den_work,pot_work,inter(:),loc_prt)

       FO_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max,2) = &
            FO_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max,2) - &
            temp_vxc_rad_LM(1:npts_den,1:nspins,1:nLM_max) / (2.0_DP*delta)
    end if


    ! Subtract off PS core density and calculate double-counting term
    ! not required if this is a tddft calculation
    ! tjz21: removed tddft option here
    if (.not. (loc_first_order)) then
       do is=1,nspins
          density_rad_LM(1:npts_den,is,1,2) =density_rad_LM(1:npts_den,is,1,2) &
               - paw_sp(isp)%tcore_den_rad(1:npts_den) * sqrt_4pi &
               / real(nspins,kind=DP)
       end do

       ! Calculate Double Counting correction to PS XC energy
       call paw_xc_exc_dc(etxc_dc,vxc_rad_LM(:,:,:,2), &
            density_rad_LM(:,:,:,2),paw_sp(isp)%grid(igrid)%r(:), &
            paw_sp(isp)%grid(igrid)%rab(:),paw_sp(isp)%rcut, &
            npts_phi,npts_max,nspins,nLM,nLM_max,den_work(:,1,1,1),inter(:))
    endif

    ! Print PS potential
    if (loc_prt) call print_rad(vxc_rad_LM(:,:,:,2),paw_sp(isp)%grid(igrid)%r, &
         npts_den,nspins,nLM,'vxc(tn^1+nhat+tn_c)',66)

    ! Calculate matrix elements of dij_xc (first-order dij_xc if first_order)
    do is=1,nspins
       do lup=0,lmax
          do mup=-lup,lup
             iLM = lup**2 + lup + 1 + mup

             ! Double loop over projectors ipwtot,jpwtot
             do ipwtot=1,paw_sp(isp)%npw_tot
                ipw = paw_sp(isp)%ipw_tot(ipwtot)
                li = paw_sp(isp)%l_pw_tot(ipwtot)
                mi = paw_sp(isp)%m_pw_tot(ipwtot)
                do jpwtot=1,paw_sp(isp)%npw_tot
                   jpw = paw_sp(isp)%ipw_tot(jpwtot)
                   lj = paw_sp(isp)%l_pw_tot(jpwtot)
                   mj = paw_sp(isp)%m_pw_tot(jpwtot)

                   ! Find Gaunt coefficient for this combination of indices
                   ! G = G^LM_{l_i m_i l_j m_j}
                   rgij = realgaunt(lup,mup,li,mi,lj,mj)

                   ! Term 1: \sum_LM G \int_0^{r_c} dr
                   !           v_xc[n^1 + n_c](r) \phi_i(r) \phi_j(r)
                   ! if loc_first_order, use FO_v_xc instead of v_xc
                   if (loc_first_order) then
                      den_work(1:npts_phi,1,1,1) = &
                           paw_sp(isp)%phi_rad(1:npts_phi,ipw) * &
                           paw_sp(isp)%phi_rad(1:npts_phi,jpw) * &
                           FO_vxc_rad_LM(1:npts_phi,is,iLM,1)
                   else
                      den_work(1:npts_phi,1,1,1) = &
                           paw_sp(isp)%phi_rad(1:npts_phi,ipw) * &
                           paw_sp(isp)%phi_rad(1:npts_phi,jpw) * &
                           vxc_rad_LM(1:npts_phi,is,iLM,1)
                   end if

                   dijxc(ipwtot,jpwtot,is) = &
                        dijxc(ipwtot,jpwtot,is) + rgij * &
                        services_radial_integral_rmax(npts_phi, &
                        paw_sp(isp)%grid(igrid)%rab, &
                        paw_sp(isp)%grid(igrid)%r,paw_sp(isp)%rcut, &
                        den_work(1:npts_phi,1,1,1),inter)

                   ! Term 2: -\sum_LM G \int_0^{r_c} dr
                   !           v_xc[\tilde{n}^1 + \hat{n}* + \tilde{n}_c](r)
                   !          \tilde{\phi}_i(r) \tilde{\phi}_j(r)
                   ! * only if pub_nhat_in_xc is true
                   ! if loc_first_order, use FO_v_xc instead of v_xc
                   if (loc_first_order) then
                      den_work(1:npts_phi,1,1,2) = &
                           paw_sp(isp)%tphi_rad(1:npts_phi,ipw) * &
                           paw_sp(isp)%tphi_rad(1:npts_phi,jpw) * &
                           FO_vxc_rad_LM(1:npts_phi,is,iLM,2)
                   else
                      den_work(1:npts_phi,1,1,2) = &
                           paw_sp(isp)%tphi_rad(1:npts_phi,ipw) * &
                           paw_sp(isp)%tphi_rad(1:npts_phi,jpw) * &
                           vxc_rad_LM(1:npts_phi,is,iLM,2)
                   end if

                   dijtxc(ipwtot,jpwtot,is) = &
                        dijtxc(ipwtot,jpwtot,is) + rgij * &
                        services_radial_integral_rmax(npts_phi, &
                        paw_sp(isp)%grid(igrid)%rab, &
                        paw_sp(isp)%grid(igrid)%r,paw_sp(isp)%rcut, &
                        den_work(1:npts_phi,1,1,2),inter)

                   ! Last part not present if nhat is not included in XC
                   if (pub_nhat_in_xc) then

                      ! Term 3: -\sum_LM G \int_0^{r_c} dr
                      !           v_xc[\tilde{n}^1 + \hat{n} + \tilde{n}_c](r)
                      !           n_{n_i l_i n_j l_j}^L g_L(r) * r^2
                      ! if loc_first_order, use FO_v_xc instead of v_xc
                      if (loc_first_order) then
                         den_work(1:npts_phi,1,1,3) = &
                              paw_sp(isp)%aug_nLij(ipw,jpw,lup) * &
                              paw_sp(isp)%shape%shape_rad(1:npts_phi,lup) * &
                              FO_vxc_rad_LM(1:npts_phi,is,iLM,2) * &
                              paw_sp(isp)%grid(igrid)%r(1:npts_phi)**2.0_DP
                      else
                         den_work(1:npts_phi,1,1,3) = &
                              paw_sp(isp)%aug_nLij(ipw,jpw,lup) * &
                              paw_sp(isp)%shape%shape_rad(1:npts_phi,lup) * &
                              vxc_rad_LM(1:npts_phi,is,iLM,2) * &
                              paw_sp(isp)%grid(igrid)%r(1:npts_phi)**2.0_DP
                      end if
                      dijtxc(ipwtot,jpwtot,is) = &
                           dijtxc(ipwtot,jpwtot,is) + rgij * &
                           services_radial_integral_rmax(npts_phi, &
                           paw_sp(isp)%grid(igrid)%rab, &
                           paw_sp(isp)%grid(igrid)%r,paw_sp(isp)%rcut, &
                           den_work(1:npts_phi,1,1,3),inter)

                   end if

                end do  ! jpwtot
             end do  ! ipwtot
          end do  ! mup
       end do  ! lup

    end do  ! is

    return

    if (loc_first_order) then
       deallocate(pert_density_mi_rad_LM,stat=ierr)
       call utils_dealloc_check('paw_dij_xc_atom','pert_density_mi_rad_LM',ierr)
       deallocate(pert_density_pl_rad_LM,stat=ierr)
       call utils_dealloc_check('paw_dij_xc_atom','pert_density_pl_rad_LM',ierr)
       deallocate(temp_vxc_rad_LM,stat=ierr)
       call utils_dealloc_check('paw_dij_xc_atom','temp_vxc_rad_LM',ierr)
    end if

contains

    subroutine print_rad(rad_fun,rad,npt,ns,ncomp,name,unit)

       ! Arguments
       integer, intent(in) :: npt
       integer, intent(in) :: ns
       integer, intent(in) :: ncomp
       integer, intent(in) :: unit
       real(kind=DP), intent(in) :: rad_fun(npt,ns,ncomp)
       real(kind=DP), intent(in) :: rad(npt)
       character(*), intent(in) :: name

       ! Local Variables
       integer :: ipt
       integer :: is
       integer :: icomp

       do is=1,ns
          write(unit,*)
          write(unit,*) name
          do icomp=1,ncomp
             do ipt=1,npts_den
                write(unit,'(3i6,2f22.12)') is,icomp,ipt,rad(ipt),rad_fun(ipt,is,icomp)
             end do
             write(unit,*)
          end do
       end do

    end subroutine print_rad

  end subroutine paw_dij_xc_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_exc_core_tot(exc_core,paw_sp,par)

    !==================================================================!
    ! This subroutine adds up the core density XC energies for all the !
    ! atoms in the system                                              !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    ! Modified to remove pub_par by Robert Charlton, 29/08/2018.       !
    !==================================================================!

    use comms, only: comms_reduce, pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: exc_core
    type(PARAL_INFO), intent(in)  :: par
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)

    ! Local Variables
    integer :: isp
    integer :: loc_iat
    integer :: iat

    ! Add up exc core over all atoms in system
    exc_core = 0.0_DP
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number
       exc_core = exc_core + paw_sp(isp)%exc_core
    end do
    call comms_reduce('SUM',exc_core)

  end subroutine paw_exc_core_tot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_exc_core_init(paw_sp)

    !==================================================================!
    ! This subroutine adds up the core density XC energies for all the !
    ! atoms in the system                                              !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    !==================================================================!

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: paw_sp(:)

    ! Local Variables
    integer :: isp, num_pspecies

    ! Calculate exc core for each species if not already done
    num_pspecies = size(paw_sp)
    do isp=1,num_pspecies
       if (.not.paw_sp(isp)%core_charge_calculated) then
          call paw_exc_core_atom(paw_sp(isp)%exc_core,isp,paw_sp)
          paw_sp(isp)%core_charge_calculated = .true.
       end if
    end do

  end subroutine paw_exc_core_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_exc_core_atom(exc_core_atom,isp,paw_sp)

    !==================================================================!
    ! This subroutine calculates the core density XC energies for a    !
    ! Single species of atom.                                          !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    !==================================================================!

    use comms, only: comms_reduce
    use paw_xc, only: paw_xc_pot_rad_LM
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: exc_core_atom
    integer, intent(in) :: isp
    type(PAW_SPECIES), intent(inout) :: paw_sp(:) ! size(par%num_pspecies)

    ! Local Variables
    integer :: ierr
    integer :: is
    integer :: npts
    integer :: nspins
    integer :: igrid
    real(kind=DP), allocatable :: core_den(:,:,:)
    real(kind=DP), allocatable :: core_den_vxc(:,:,:)
    real(kind=DP), allocatable :: xc_pot_work(:,:,:,:)
    real(kind=DP), allocatable :: xc_den_work(:,:,:,:)
    real(kind=DP), allocatable :: inter(:)

    if (paw_sp(isp)%core_charge) then
       igrid = paw_sp(isp)%core_den_grid
       npts = paw_sp(isp)%grid(igrid)%npt
       nspins = 1
       allocate(core_den(npts,nspins,1),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','core_den',ierr)
       allocate(core_den_vxc(npts,nspins,1),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','core_den_vxc',ierr)
       allocate(xc_pot_work(npts,nspins,nspins,6),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','xc_pot_work',ierr)
       allocate(xc_den_work(npts,nspins,nspins,4),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','xc_den_work',ierr)
       allocate(inter(npts),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','inter',ierr)

       do is=1,nspins
          core_den(1:npts,is,1) = paw_sp(isp)%core_den_rad(1:npts) &
               * sqrt_4pi / real(nspins,kind=DP)
       end do

       ! Call radial XC function with just nLM_max = 1 (no moment expansion)
       call paw_xc_pot_rad_LM(core_den_vxc(:,:,:),exc_core_atom, &
            core_den(:,:,:),paw_sp(isp)%grid(igrid)%r, &
            paw_sp(isp)%grid(igrid)%rab,paw_sp(isp)%rcut,npts,npts, &
            nspins,1,0,xc_den_work,xc_pot_work,inter)

       deallocate(inter,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','inter',ierr)
       deallocate(xc_den_work,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','xc_den_work',ierr)
       deallocate(xc_pot_work,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','xc_pot_work',ierr)
       deallocate(core_den_vxc,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','core_den_vxc',ierr)
       deallocate(core_den,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','core_den',ierr)

    else

       paw_sp(isp)%exc_core = 0.0_DP
       paw_sp(isp)%core_charge_calculated = .true.

    end if

  end subroutine paw_exc_core_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_ae_density_rad_LM(ae_den,isp,rhoij_at,npw,npts_max,nspins, &
       nLM_max,paw_sp,den_work,first_order)

    !==================================================================!
    ! This subroutine calculates the all-electron sphere density in    !
    ! the PAW method, as a sum over the moments LM.                    !
    !  ae_den(:,is,iLM) = n^1 density + core density n_c               !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    ! Split off to a separate routine by Nicholas Hine on 09/05/11.    !
    !==================================================================!

    use gaunt_coeff, only: realgaunt
    use services, only: services_analytic_limit

    implicit none

    ! Arguments
    integer, intent(in) :: npts_max
    integer, intent(in) :: nspins
    integer, intent(in) :: npw
    integer, intent(in) :: nLM_max
    real(kind=DP), intent(out) :: ae_den(npts_max,nspins,nLM_max)
    real(kind=DP), intent(out) :: den_work(npts_max)
    real(kind=DP), intent(in) :: rhoij_at(npw,npw,nspins)
    integer, intent(in) :: isp
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)
    logical, intent(in), optional :: first_order

    ! Local Variables
    integer :: iLM
    integer :: npts_den
    integer :: npts_phi
    integer :: kpwtot,lpwtot
    integer :: kpw,lpw
    integer :: lk,mk,ll,ml
    integer :: lup,mup
    integer :: is
    integer :: igrid
    real(kind=DP) :: rgkl
    real(kind=DP),parameter :: gaunt_tol = 1.0e-14_DP
    real(kind=DP),parameter :: rho_tol = 1.0e-14_DP
    logical :: loc_first_order

    if (present(first_order)) then
       loc_first_order = first_order
    else
       loc_first_order = .false.
    end if

    ! Initialise, and find information about this species and its phi_grid
    ae_den(:,:,:) = 0.0_DP
    igrid = paw_sp(isp)%phi_grid
    npts_phi = paw_sp(isp)%grid(igrid)%npt
    igrid = paw_sp(isp)%core_den_grid
    npts_den = paw_sp(isp)%grid(igrid)%npt

    ! Double loop over partial waves on this atom
    do kpwtot=1,paw_sp(isp)%npw_tot
       kpw = paw_sp(isp)%ipw_tot(kpwtot)
       lk = paw_sp(isp)%l_pw_tot(kpwtot)
       mk = paw_sp(isp)%m_pw_tot(kpwtot)
       do lpwtot=1,paw_sp(isp)%npw_tot
          lpw = paw_sp(isp)%ipw_tot(lpwtot)
          ll = paw_sp(isp)%l_pw_tot(lpwtot)
          ml = paw_sp(isp)%m_pw_tot(lpwtot)

          if (all(abs(rhoij_at(kpwtot,lpwtot,:))<rho_tol)) cycle

          ! den_work(:) = tphi_i(r)*tphi_j(r) / r^2
          den_work(2:npts_phi) = paw_sp(isp)%phi_rad(2:npts_phi,kpw) &
               * paw_sp(isp)%phi_rad(2:npts_phi,lpw) &
               / paw_sp(isp)%grid(igrid)%r(2:npts_phi)**2
          ! Find analytic limit
          call services_analytic_limit(npts_phi,paw_sp(isp)%grid(igrid)%r(:), &
               den_work(:))

          ! Loop over angular momentum channels
          do lup=0,paw_sp(isp)%lmax
             do mup=-lup,lup
                rgkl = realgaunt(lup,mup,lk,mk,ll,ml)
                if (abs(rgkl)<gaunt_tol) cycle
                iLM = lup*lup + lup + 1 + mup
                do is=1,nspins
                   ae_den(1:npts_phi,is,iLM) = ae_den(1:npts_phi,is,iLM) + &
                        rhoij_at(kpwtot,lpwtot,is)*den_work(1:npts_phi)*rgkl
                end do
             end do
          end do
       end do
    end do

    if (.not.loc_first_order) then
    ! Add tcore density to tn1_{LM=00}(r)
       do is=1,nspins
          ae_den(1:npts_den,is,1) = ae_den(1:npts_den,is,1) + &
               paw_sp(isp)%core_den_rad(1:npts_den) &
               * sqrt_4pi / real(nspins,kind=DP)
       end do
    end if

  end subroutine paw_ae_density_rad_LM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_ps_density_rad_LM(ps_den,isp,rhoij_at,npw,npts_max,nspins, &
       nLM_max,paw_sp,den_work,first_order)

    !==================================================================!
    ! This subroutine calculates the soft pseudo sphere density in     !
    ! the PAW method, as a sum over the moments LM.                    !
    !  ps_den(:,is,iLM) = \tilde{n^1} density + psc density \tilde{n_c}!
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    ! Split off to a separate routine by Nicholas Hine on 09/05/11.    !
    !==================================================================!

    use gaunt_coeff, only: realgaunt
    use services, only: services_analytic_limit

    implicit none

    ! Arguments
    integer, intent(in) :: npts_max
    integer, intent(in) :: nspins
    integer, intent(in) :: npw
    integer, intent(in) :: nLM_max
    real(kind=DP), intent(out) :: ps_den(npts_max,nspins,nLM_max)
    real(kind=DP), intent(out) :: den_work(npts_max)
    real(kind=DP), intent(in) :: rhoij_at(npw,npw,nspins)
    integer, intent(in) :: isp
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)
    logical, intent(in),optional :: first_order

    ! Local Variables
    integer :: iLM
    integer :: npts_den
    integer :: npts_phi
    integer :: kpwtot,lpwtot
    integer :: kpw,lpw
    integer :: lk,mk,ll,ml
    integer :: lup,mup
    integer :: is
    integer :: igrid
    real(kind=DP) :: rgkl
    real(kind=DP),parameter :: gaunt_tol = 1.0e-14_DP
    real(kind=DP),parameter :: rho_tol = 1.0e-14_DP
    logical :: loc_first_order

    if (present(first_order)) then
       loc_first_order = first_order
    else
       loc_first_order = .false.
    end if

    ! Initialise, and find information about this species and its phi_grid
    ps_den(:,:,:) = 0.0_DP
    igrid = paw_sp(isp)%phi_grid
    npts_phi = paw_sp(isp)%grid(igrid)%npt
    igrid = paw_sp(isp)%core_den_grid
    npts_den = paw_sp(isp)%grid(igrid)%npt

    ! Double loop over partial waves on this atom
    do kpwtot=1,paw_sp(isp)%npw_tot
       kpw = paw_sp(isp)%ipw_tot(kpwtot)
       lk = paw_sp(isp)%l_pw_tot(kpwtot)
       mk = paw_sp(isp)%m_pw_tot(kpwtot)
       do lpwtot=1,paw_sp(isp)%npw_tot
          lpw = paw_sp(isp)%ipw_tot(lpwtot)
          ll = paw_sp(isp)%l_pw_tot(lpwtot)
          ml = paw_sp(isp)%m_pw_tot(lpwtot)

          if (all(abs(rhoij_at(kpwtot,lpwtot,:))<rho_tol)) cycle

          ! den_work(:) = tphi_i(r)*tphi_j(r) / r^2
          den_work(2:npts_phi) = paw_sp(isp)%tphi_rad(2:npts_phi,kpw) &
               * paw_sp(isp)%tphi_rad(2:npts_phi,lpw) &
               / paw_sp(isp)%grid(igrid)%r(2:npts_phi)**2
          ! Find analytic limit
          call services_analytic_limit(npts_phi,paw_sp(isp)%grid(igrid)%r(:), &
               den_work(:))

          ! Loop over angular momentum channels
          do lup=0,paw_sp(isp)%lmax
             do mup=-lup,lup
                rgkl = realgaunt(lup,mup,lk,mk,ll,ml)
                if (abs(rgkl)<gaunt_tol) cycle
                iLM = lup*lup + lup + 1 + mup
                do is=1,nspins
                   ps_den(1:npts_phi,is,iLM) = ps_den(1:npts_phi,is,iLM) + &
                        rhoij_at(kpwtot,lpwtot,is)*den_work(1:npts_phi)*rgkl
                end do
             end do
          end do
       end do
    end do

    if(.not.loc_first_order) then
       ! Add tcore density to tn1_{LM=00}(r)
       do is=1,nspins
          ps_den(1:npts_den,is,1) = ps_den(1:npts_den,is,1) + &
               paw_sp(isp)%tcore_den_rad(1:npts_den) &
               * sqrt_4pi / real(nspins,kind=DP)
       end do
    end if

  end subroutine paw_ps_density_rad_LM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine paw_aug_density_rad_LM(aug_den,isp,rhoij_at,npw,npts_max, &
       nspins,nLM_max,paw_sp,den_work)

    !==================================================================!
    ! This subroutine calculates the various sphere-only densities in  !
    ! the PAW as a sums over their moments LM.                         !
    !  den(:,is,iLM,1) = n^1 density + core density n_c                !
    !  den(:,is,iLM,2) = \tilde{n^1} density + psc density \tilde{n_c} !
    !  den(:,is,iLM,3) = \hat{n} compensation density                  !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    !==================================================================!

    use gaunt_coeff, only: realgaunt

    implicit none

    ! Arguments
    integer, intent(in) :: npts_max
    integer, intent(in) :: nspins
    integer, intent(in) :: npw
    integer, intent(in) :: nLM_max
    real(kind=DP), intent(out) :: aug_den(npts_max,nspins,nLM_max)
    real(kind=DP), intent(out) :: den_work(npts_max)
    real(kind=DP), intent(in) :: rhoij_at(npw,npw,nspins)
    integer, intent(in) :: isp
    type(PAW_SPECIES), intent(in) :: paw_sp(:) ! size(par%num_pspecies)

    ! Local Variables
    integer :: iLM
    integer :: npts_shp
    integer :: kpwtot,lpwtot
    integer :: kpw,lpw
    integer :: lk,mk,ll,ml
    integer :: lup,mup
    integer :: is
    integer :: igrid
    real(kind=DP) :: rgkl
    real(kind=DP),parameter :: gaunt_tol = 1.0e-14_DP
    real(kind=DP),parameter :: rho_tol = 1.0e-14_DP

    ! Initialise, and find information about this species and its phi_grid
    aug_den(:,:,:) = 0.0_DP
    igrid = paw_sp(isp)%shape_grid
    npts_shp = paw_sp(isp)%grid(igrid)%npt

    ! Double loop over partial waves on this atom
    do kpwtot=1,paw_sp(isp)%npw_tot
       kpw = paw_sp(isp)%ipw_tot(kpwtot)
       lk = paw_sp(isp)%l_pw_tot(kpwtot)
       mk = paw_sp(isp)%m_pw_tot(kpwtot)
       do lpwtot=1,paw_sp(isp)%npw_tot
          lpw = paw_sp(isp)%ipw_tot(lpwtot)
          ll = paw_sp(isp)%l_pw_tot(lpwtot)
          ml = paw_sp(isp)%m_pw_tot(lpwtot)

          if (all(abs(rhoij_at(kpwtot,lpwtot,:))<rho_tol)) cycle

          ! Loop over angular momentum channels
          do lup=0,paw_sp(isp)%lmax
             ! work(3) = g_L(r)*q_ij
             den_work(1:npts_shp) = paw_sp(isp)%aug_nLij(kpw,lpw,lup) * &
                  paw_sp(isp)%shape%shape_rad(1:npts_shp,lup)
             do mup=-lup,lup
                rgkl = realgaunt(lup,mup,lk,mk,ll,ml)
                if (abs(rgkl)<gaunt_tol) cycle
                iLM = lup*lup + lup + 1 + mup
                do is=1,nspins
                   aug_den(1:npts_shp,is,iLM) = aug_den(1:npts_shp,is,iLM) + &
                        rhoij_at(kpwtot,lpwtot,is)*den_work(1:npts_shp)*rgkl
                end do
             end do
          end do
       end do
    end do

  end subroutine paw_aug_density_rad_LM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_so(dijso,num_spin_dens,rho_ij,paw_sp,atom,cart)

    !=================================================================!
    ! This subroutine calculates the spin orbit terms in the nonlocal !
    ! energies Dij.                                                   !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !   dijso(inout) : spin-orbit contribution to nonlocal energy     !
    !                       term hat{D}_ij.                           !
    !   rho_ij(in)   : Projector density kernel rho^ij                !
    !   atom, cart   : if a first-order calculation if required, these!
    !        variables show the atom and direction of perturbation    !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 15/03/16.                           !
    ! Fixes by Gabriel Constantinescu and JM Escartin, Summer 2016.   !
    !=================================================================!

    use comms, only: pub_my_proc_id
    use constants, only: max_spins
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: num_spin_dens
    type(SPAM3),intent(inout) :: dijso(num_spin_dens)
    type(SPAM3),intent(in) :: rho_ij(pub_num_spins)
    type(PAW_SPECIES), intent(in) :: paw_sp(:)
    integer, intent(in), optional :: atom
    integer, intent(in), optional :: cart

    ! if 'atom' and 'cart' are present, rho_ij is actually the first-order
    ! projector density kernel

    ! Local Variables
    integer :: loc_iat, iat, ierr
    integer :: isp
    integer :: is
    integer :: igrid
    integer :: npts_den

    real(kind=DP), allocatable :: veff_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: density_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: rhoij_at(:,:,:)
    real(kind=DP), allocatable :: dijxc_at(:,:,:)
    real(kind=DP), allocatable :: dijtxc_at(:,:,:)
    real(kind=DP), allocatable :: pot_work(:,:,:,:)
    real(kind=DP), allocatable :: den_work(:,:,:,:)
    real(kind=DP), allocatable :: inter(:)
    real(kind=DP) :: exc_at
    real(kind=DP) :: exc_dc_at
    real(kind=DP) :: etxc_at
    real(kind=DP) :: etxc_dc_at
    real(kind=DP) :: total_nhat_at(max_spins)

    complex(kind=DP), allocatable :: dijso_at(:,:,:)
    integer :: npts_max
    integer :: nLM_max
    integer :: nspins
    logical :: do_step
    type(PARAL_INFO), pointer :: par

    ! Start Timer
    call timer_clock('paw_dij_so',1)

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, dijso(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_dij_so: allocated parallel strategy is &
         &incompatible with paw_sp.')

    ! Find array sizes
    npts_max = 0
    nLM_max = 0
    nspins = pub_num_spins
    do isp=1,par%num_pspecies
       if (.not.any(par%elements_on_proc(:)%pspecies_number==isp)) cycle
       igrid = paw_sp(isp)%phi_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       igrid = paw_sp(isp)%shape_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       igrid = paw_sp(isp)%core_den_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       nLM_max = max(nLM_max,(paw_sp(isp)%lmax+1)**2)
    end do

    ! Allocate temporary arrays
    allocate(veff_rad_LM(npts_max,nspins,nLM_max,2),stat=ierr)
    call utils_alloc_check('paw_dij_so','veff_rad_LM',ierr)
    allocate(density_rad_LM(npts_max,nspins,nLM_max,3),stat=ierr)
    call utils_alloc_check('paw_dij_so','density_rad_LM',ierr)
    allocate(rhoij_at(max_paw_proj_tot,max_paw_proj_tot, &
         nspins),stat=ierr)
    call utils_alloc_check('paw_dij_so','rhoij_at',ierr)
    allocate(dijxc_at(max_paw_proj_tot,max_paw_proj_tot, &
         nspins),stat=ierr)
    call utils_alloc_check('paw_dij_so','dijxc_at',ierr)
    allocate(dijtxc_at(max_paw_proj_tot,max_paw_proj_tot, &
         nspins),stat=ierr)
    call utils_alloc_check('paw_dij_so','dijtxc_at',ierr)
    allocate(pot_work(npts_max,nspins,nspins,6), &
         stat=ierr)
    call utils_alloc_check('paw_dij_so','pot_work',ierr)
    allocate(den_work(npts_max,nspins,nspins,4), &
         stat=ierr)
    call utils_alloc_check('paw_dij_so','den_work',ierr)
    allocate(inter(npts_max),stat=ierr)
    call utils_alloc_check('paw_dij_so','inter',ierr)

    allocate(dijso_at(max_paw_proj_tot,max_paw_proj_tot,num_spin_dens),stat=ierr)
    call utils_alloc_check('paw_dij_so','disjo_at',ierr)

    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       if ( present(atom).and.present(cart) ) then ! if first-order
          if (par%orig_atom(iat) == atom) then
             do_step = .true.
          else
             do_step = .false.
          end if
       else                                        ! if NOT first-order
          do_step = .true.
       end if

       if (do_step) then

          ! Get rhoij_at for this atom (summed over spins if necessary
          rhoij_at(:,:,:) = 0.0_DP
          do is=1,nspins
             call sparse_get_block(rhoij_at(:,:,is),rho_ij(is),iat,iat)
          end do
          if (nspins==2) then
             rhoij_at(:,:,1) = rhoij_at(:,:,1) + rhoij_at(:,:,2)
          end if

          call paw_dij_xc_atom(isp,nspins,npts_max,nLM_max, &
               max_paw_proj_tot,rhoij_at,dijxc_at,dijtxc_at, &
               exc_at,etxc_at,exc_dc_at,etxc_dc_at,total_nhat_at, &
               density_rad_LM,veff_rad_LM,den_work,pot_work,inter,paw_sp, &
               .false.)

          igrid = paw_sp(isp)%core_den_grid
          npts_den = paw_sp(isp)%grid(igrid)%npt
          do is = 1,nspins
             density_rad_LM(1:npts_den,is,1,1) = &
                  density_rad_LM(1:npts_den,is,1,1) &
                  + paw_sp(isp)%core_den_rad(1:npts_den) * sqrt_4pi &
                  / real(nspins,kind=DP)
          end do

          ! calculate hartree potential of n1 density plus core density
          call paw_vhartree_rad(veff_rad_LM(:,:,:,2),igrid,nspins,npts_max, &
               nLM_max,density_rad_LM(:,:,:,1),den_work,paw_sp(isp))


          ! Add up terms in total effective potential
          do is = 1,nspins

             ! Add Hartree potential of n1 density to xc potential of n1 density
             veff_rad_LM(1:npts_den,is,1,1) = veff_rad_LM(1:npts_den,is,1,1) + &
                  veff_rad_LM(1:npts_den,is,1,2)

             ! Add potential of nucleus with charge Z
             if (paw_sp(isp)%grid(igrid)%r(1)==0.0_DP) then
                veff_rad_LM(2:npts_den,is,1,1) = &
                     veff_rad_LM(2:npts_den,is,1,1) - &
                     paw_sp(isp)%atomic_number / &
                     paw_sp(isp)%grid(igrid)%r(2:npts_den) * sqrt_4pi &
                     / real(nspins,kind=DP)
             else
                veff_rad_LM(1:npts_den,is,1,1) = &
                     veff_rad_LM(1:npts_den,is,1,1) - &
                     paw_sp(isp)%atomic_number / &
                     paw_sp(isp)%grid(igrid)%r(1:npts_den) * sqrt_4pi &
                     / real(nspins,kind=DP)
             end if


          end do

          ! Calculate so terms in nonlocal energies for this atom
          call paw_dij_so_atom(dijso_at,npts_max,nspins,num_spin_dens, &
               max_paw_proj_tot,nLM_max,igrid,veff_rad_LM(:,:,:,1), &
               den_work,paw_sp(isp))

          do is = 1,num_spin_dens
             call sparse_put_block(dijso_at(:,:,is),dijso(is),iat,iat)
          end do

       end if ! if do_step

    end do


    ! Deallocate temporary arrays
    deallocate(dijso_at,stat=ierr)
    call utils_dealloc_check('paw_dij_so','dijso_at',ierr)

    deallocate(inter,stat=ierr)
    call utils_dealloc_check('paw_dij_so','inter',ierr)
    deallocate(den_work,stat=ierr)
    call utils_dealloc_check('paw_dij_so','den_work',ierr)
    deallocate(pot_work,stat=ierr)
    call utils_dealloc_check('paw_dij_so','pot_work',ierr)
    deallocate(dijtxc_at,stat=ierr)
    call utils_dealloc_check('paw_dij_so','dijtxc_at',ierr)
    deallocate(dijxc_at,stat=ierr)
    call utils_dealloc_check('paw_dij_so','dijxc_at',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('paw_dij_so','rhoij_at',ierr)
    deallocate(density_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_dij_so','density_rad_LM',ierr)
    deallocate(veff_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_dij_so','veff_rad_LM',ierr)
    nullify(par)

    ! Stop Timer
    call timer_clock('paw_dij_so',2)

  end subroutine paw_dij_so



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_vhartree_rad(vh_rad_LM,igrid,nspins,npts_max, &
       nLM_max,density_rad_LM,work,species)

    !=========================================================================!
    ! This subroutine calculates the hartree potential of a density on a PAW  !
    ! radial grid                                                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 18/03/2016.                                 !
    ! Some fixes by Gabriel Constantinescu and JM Escartin, Summer 2016.      !
    !=========================================================================!

    use services, only: services_radial_integral
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: npts_max
    integer, intent(in) :: nspins
    integer, intent(in) :: nLM_max
    real(kind=DP), intent(out) :: vh_rad_LM(npts_max,nspins,nLM_max)
    integer, intent(in) :: igrid
    real(kind=DP), intent(in) :: density_rad_LM(npts_max,nspins,nLM_max)
    real(kind=DP), intent(out) :: work(npts_max,nspins,nspins,4)
    type(PAW_SPECIES), intent(in) :: species

    ! Local Variables
    integer :: lup, mup
    integer :: is
    integer :: iLM
    integer :: ierr
    integer :: npts
    real(kind=DP) :: lfac
    real(kind=DP) :: int1,int2
    real(kind=DP), allocatable :: rwork(:,:)

    npts = species%grid(igrid)%npt
    allocate(rwork(npts,4),stat=ierr)
    call utils_alloc_check('paw_vhartree_rad','rwork',ierr)
    work(:,:,:,:) = 0.0_DP

    do is=1,nspins
       do lup=0,species%lmax
          do mup=-lup,lup
             iLM = lup*lup + lup + 1 + mup
             lfac = 4.0_DP*PI / real(2*lup+1,kind=DP)

             rwork(2:npts,1) = species%grid(igrid)%r(2:npts)**lup
             rwork(2:npts,2) = 1.0_DP/(species%grid(igrid)%r(2:npts)**(lup+1))
             rwork(1:npts,3) = species%grid(igrid)%r(1:npts)**(lup+2)
             rwork(2:npts,4) = species%grid(igrid)%r(2:npts)**(1-lup)

             ! handle r=0 separately depending on whether it is included in grid
             if (species%grid(igrid)%r(1)==0.0_DP) then
                rwork(1,2) = 0.0_DP
                rwork(1,1) = 0.0_DP
                if (lup.ge.1) then
                   rwork(1,4) = 0.0_DP
                else
                   rwork(1,4) = species%grid(igrid)%r(1)**(1-lup)
                end if
             else
                rwork(1,1) = species%grid(igrid)%r(1)**lup
                rwork(1,2) = 1.0_DP/(species%grid(igrid)%r(1)**(lup+1))
                rwork(1,4) = species%grid(igrid)%r(1)**(1-lup)
             end if

             ! work(:)  = n_LM(r) * rp^L / r ^(L+1) * rp^2
             ! work2(:) = n_LM(r) * r ^L / rp^(L+1) * rp^2
             work(1:npts,1,1,1) = rwork(1:npts,3) * density_rad_LM(1:npts,is,iLM)
             work(1:npts,1,1,2) = rwork(1:npts,4) * density_rad_LM(1:npts,is,iLM)
             int1 = services_radial_integral(npts,species%grid(igrid)%rab, &
                  work(1:npts,1,1,1),work(1:npts,1,1,3))
             int2 = services_radial_integral(npts,species%grid(igrid)%rab, &
                  work(1:npts,1,1,2),work(1:npts,1,1,4))
             vh_rad_LM(1:npts,is,iLM) = lfac * (work(1:npts,1,1,3) * rwork(1:npts,2) + &
                  (int2 - work(1:npts,1,1,4)) * rwork(1:npts,1))
          end do

       end do
    end do

    deallocate(rwork,stat=ierr)
    call utils_dealloc_check('paw_vhartree_rad','rwork',ierr)

  end subroutine paw_vhartree_rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_so_atom(dijso_at,npts,nspins,num_spin_dens,npwtot, &
       nLM_max,igrid,veff_rad_LM,work,species)

    !=========================================================================!
    ! This subroutine calculates the spin-orbit Dij terms                     !
    ! for an individual atom.                                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 18/03/2016.                                 !
    ! Multiple essential fixes and enhancements by Jose M Escartin and        !
    ! Gabriel Constantinescu, Summer 2016.                                    !
    !=========================================================================!

    use constants, only: FINE_STRUCTURE, cmplx_0, cmplx_i, SQRT_PI
    use services, only: services_radial_derivative, services_radial_integral_rmax
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    integer, intent(in) :: nspins
    integer, intent(in) :: num_spin_dens
    integer, intent(in) :: npts
    integer, intent(in) :: npwtot
    integer, intent(in) :: nLM_max
    integer, intent(in) :: igrid ! mesh the core densities use
    complex(kind=DP), intent(out) :: dijso_at(npwtot,npwtot,num_spin_dens)
    real(kind=DP), intent(in) :: veff_rad_LM(npts,nspins,nLM_max)
    real(kind=DP), intent(out) :: work(npts,nspins,nspins,4)
    type(PAW_SPECIES), intent(in) :: species

    ! Local Variables
    integer :: is
    integer :: iLM
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,lj,mi,mj
    integer :: mu,mup
    real(kind=DP) :: radint, fac
    real(kind=DP), parameter :: alpha2 = FINE_STRUCTURE**2
    complex(kind=DP) :: Slimi_Lop_Sljmj, c_mup_mi, c_mu_mj
    integer, parameter :: LZ=1,MINUSLZ=2,LPLUS=3,LMINUS=4
    integer :: lop
    integer :: lsgn
    integer :: npts_den, npts_phi
    real(kind=DP), allocatable :: veff_rad_00_average(:)
    integer :: ierr

    call utils_assert(num_spin_dens.le.4, 'Error in paw_dij_so_atom: too many &
         &components required.')

    allocate(veff_rad_00_average(npts), stat=ierr)
    call utils_alloc_check('paw_dij_so_atom', 'veff_rad_00_average', ierr)

    dijso_at(:,:,:) = cmplx_0

    npts_den = species%grid(igrid)%npt
    npts_phi =  species%grid(species%phi_grid)%npt

    ! jme: we average veff over spins.
    ! For polarised calculations, this should lead to Hermiticity of
    ! the full Hamiltonian (ud vs du blocks), at the cost of changing
    ! the uu and dd blocks.
    ! This is the same strategy that ABINIT uses.
    if (nspins == 1) then
       veff_rad_00_average(:) = veff_rad_LM(:,1,1)
    else
       veff_rad_00_average(:) = 0.5_DP * (veff_rad_LM(:,1,1) + &
            veff_rad_LM(:,2,1))
    end if

    ! jme: reinstate Y_{00} factor:
    veff_rad_00_average(:) = veff_rad_00_average(:) / (2.0_DP * SQRT_PI)

    ! Calculate dV_eff[n1] / dr
    call services_radial_derivative(work(:,1,1,1),veff_rad_00_average, &
         npts,real(npts,kind=DP))
    ! required scaling of derivative for logarithmic grid (implied)
    work(1:npts_den,1,1,1) = work(1:npts_den,1,1,1) / &
         species%grid(igrid)%rab(1:npts_den)

    ! Divide by r (avoiding divide by zero)
    if (species%grid(igrid)%r(1)==0.0_DP) then
       work(1,1,1,2) = 0.0_DP
       work(2:npts_den,1,1,2) = work(2:npts_den,1,1,1) / &
            species%grid(igrid)%r(2:npts_den)
    else
       work(1:npts_den,1,1,2) = work(1:npts_den,1,1,1) / &
            species%grid(igrid)%r(1:npts_den)
    end if

    ! ZORA correction
    work(1:npts_den,1,1,2) = work(1:npts_den,1,1,2) / &
         ( 1.0_DP - alpha2 * veff_rad_00_average(1:npts_den) )

    deallocate(veff_rad_00_average, stat=ierr)
    call utils_dealloc_check('paw_dij_so_atom', 'veff_rad_00_average', ierr)

    ! Loop over spins
    do is=1,num_spin_dens

       ! Choose operator for angular matrix element.
       select case (is)
       case (1)
          lop = LZ
       case (2)
          lop = MINUSLZ
       case (3)
          lop = LPLUS
       case (4)
          lop = LMINUS
       end select

       ! Double loop over projectors ipwtot,jpwtot
       do jpwtot=1,species%npw_tot
          jpw = species%ipw_tot(jpwtot)
          lj = species%l_pw_tot(jpwtot)
          mj = species%m_pw_tot(jpwtot)

          ! No spin-orbit couplings in ell=0 subspaces.
          if (lj==0) cycle

          do ipwtot=1,species%npw_tot
             ipw = species%ipw_tot(ipwtot)
             li = species%l_pw_tot(ipwtot)
             mi = species%m_pw_tot(ipwtot)

             ! Selection rule on nonzero coeffs since L op preserves l value
             if (lj/=li) cycle

             ! Angular matrix element
             Slimi_Lop_Sljmj = cmplx_0

             ! L_z or -L_z operators
             if ((lop==LZ).or.(lop==MINUSLZ)) then
                if (mi == -mj) then
                   if (lop==LZ) then
                      lsgn = 1
                   else
                      lsgn = -1
                   end if
                   Slimi_Lop_Sljmj = lsgn*cmplx_i*real(mj,kind=DP)
                end if

             ! L_+ or L_- operators
             else ! ((lop==LPLUS).or.(lop==LMINUS))

                ! Loop over mu = -m,+m (avoid repeating zero)
                do mu = -abs(mj), abs(mj), max(2*abs(mj),1)
                   if (lop==LPLUS) then
                      mup = mu+1
                   else ! (lop==LMINUS)
                      mup = mu-1
                   end if
                   if (abs(mup)/=abs(mi)) cycle
                   c_mu_mj  = c_mu_m(mu, mj)
                   c_mup_mi = c_mu_m(mup,mi)
                   fac = sqrt(real(lj*(lj+1)-mu*mup,kind=DP))
                   Slimi_Lop_Sljmj = Slimi_Lop_Sljmj + &
                           conjg(c_mup_mi) * c_mu_mj * fac
                end do
             end if

             ! cycle if angular element is zero
             if (abs(Slimi_Lop_Sljmj)<1e-10) cycle

             ! Radial integral
             work(1:npts_phi,1,1,3) = work(1:npts_phi,1,1,2) * &
                  species%phi_rad(1:npts_phi,ipw) * &
                  species%phi_rad(1:npts_phi,jpw)

             if (npts>npts_phi) then
                work((npts_phi+1):npts,1,1,3) = 0.0_DP
             end if

             radint = services_radial_integral_rmax(npts_den, &
                  species%grid(igrid)%rab,species%grid(igrid)%r, &
                  species%rcut,work(:,1,1,3),work(:,1,1,4))

             ! jme: (alpha^2/4) <r^1 Y_{00} dV0/dr> <L.sigma>
             ! (S = sigma / 2)
             dijso_at(ipwtot,jpwtot,is) = &
                  0.25_DP * alpha2 * radint * Slimi_Lop_Sljmj

          end do
       end do
    end do

contains

   pure complex(kind=DP) function c_mu_m(mu,m)

      use constants, only: cmplx_0, cmplx_1, cmplx_i, SQRT_TWO
      implicit none
      integer, intent(in) :: mu,m
      real(kind=DP), parameter :: invsq2 = 1.0_DP / SQRT_TWO

      c_mu_m = cmplx_0

      if (m>0) then
         if (mu==-m) c_mu_m = invsq2*cmplx_1
         if (mu==m)  c_mu_m = invsq2*cmplx_1*(-1.0_DP)**m
      else if (m==0) then
         if (mu==0)  c_mu_m = cmplx_1
      else if (m<0) then
         if (mu==m)  c_mu_m = invsq2*cmplx_i
         if (mu==-m) c_mu_m = -invsq2*cmplx_i*(-1.0_DP)**m
      end if

   end function c_mu_m

  end subroutine paw_dij_so_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_tcore_hartree_calc_forces(den_slabs12,grid,elements, &
       paw_sp,loc_forces,par)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the Hartree potential of the pseudized core charge.                !
    !  F_I = d/dR_I (\int (\tilde{n} + \hat{n})(r) v_H[\tilde{n}_{Zc}](r)) dr !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  den_slabs12  : input  : ground state data-parallelised charge density  !
    !  elements     : input  : list of elements and corresponding info        !
    !  loc_forces   : output : ionic forces due to local part                 !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 09/06/2010.                                 !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use cell_grid, only: cell_grid_recip_pt, GRID_INFO
    use comms, only: comms_reduce, pub_my_proc_id
    use constants, only: PI
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins
!$  use rundat, only: pub_threads_max
    use services, only: services_1d_interpolation
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in) :: par
    type(GRID_INFO), intent(in) :: grid
    real(kind=dp), intent(inout) :: den_slabs12(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    real(kind=dp), intent(out) :: loc_forces(1:3,par%nat)

    ! Local Variables
    integer :: ipt,i2,i3,islab23
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gvec(3)
    real(kind=DP) :: g_length,inv_gsq,gdotR
    real(kind=DP) :: factor,v_loc_value
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: iG_vden(:,:)
    complex(kind=DP), allocatable :: recip(:,:,:)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_tcore_hartree_calc_forces'

    ! Start timer
    call timer_clock('paw_tcore_hartree_calc_forces',1)

    ! Allocate
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_tcore_hartree_calc_forces','recip',ierr)

    ! Initialise
    loc_forces = 0.0_DP

    ! If spin polarised, put total density in up spin
    if (pub_num_spins == 2) den_slabs12(:,:,:,1) = &
         den_slabs12(:,:,:,1) + den_slabs12(:,:,:,2)

    ! Fourier transform the charge density to reciprocal space
    call fourier_apply_cell_forward(den_slabs12(:,:,:,1),recip,grid)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i2,i3,gvec,g_length,inv_gsq,species,v_loc_value, &
!$OMP      iG_vden,factor,islab23,atom,gdotR,eiGR,ierr) &
!$OMP SHARED (pub_my_proc_id,grid,paw_sp,elements,recip,par, &
!$OMP      pub_threads_max) &
!$OMP REDUCTION (+:loc_forces)

    allocate(iG_vden(1:3,par%num_pspecies),stat=ierr)
    call utils_alloc_check('paw_tcore_hartree_calc_forces','iG_vden',ierr)

    ! For components g with symmetry points at -g
    factor = 2.0_DP

!$OMP DO
    do ipt=1,grid%num_slabs23*grid%n2*grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1

       if ((islab23==1).and.(grid%proc_slab23(1) == pub_my_proc_id)) &
            factor = 1.0_DP

       ! g-vector, |g| and 1/|g|^2
       call cell_grid_recip_pt(gvec,islab23 + &
            grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)
       g_length = sqrt(sum(gvec(1:3)**2))
       inv_gsq = grid%coulomb_recip(i3,i2,islab23)

       ! ndmh: loop over species, calculating v(G) for each.
       ! ndmh: this bit is O(N) so as much as possible should be
       ! ndmh: pre-calculated here.
       do species=1,par%num_pspecies

          ! Interpolate value of local potential at current g
          v_loc_value = services_1d_interpolation( &
               paw_sp(species)%vhntzc_recip, &
               paw_sp(species)%n_recip_pts,g_length * &
               paw_sp(species)%inv_g_spacing,0)

          ! Add back the Coulomb potential; set g=0 term to zero
          if (g_length .gt. 0.0_DP) then
             ! pa: changed to allow fractional ionic charge
             v_loc_value = v_loc_value - &
                  paw_sp(species)%ion_charge*inv_gsq
          else
             v_loc_value = 0.0_DP
          endif

          ! ndmh: calculate iG.v(G).n*(G) for each species and
          ! ndmh: each Cartesian direction
          iG_vden(1:3,species) = factor * cmplx(0.0_DP,1.0_DP,kind=DP) &
               * gvec(1:3) * v_loc_value * conjg(recip(i3,i2,islab23))

       end do

       ! ndmh: loop over atoms: this bit is O(N^2) so should be made as
       ! ndmh: short and simple as possible
       do atom=1,par%nat

          ! e^{-ig.R}
          gdotR =  -(gvec(1)*elements(atom)%centre%x + &
               gvec(2)*elements(atom)%centre%y + &
               gvec(3)*elements(atom)%centre%z)
          eiGR = cmplx(COS(gdotR),SIN(gdotR),kind=DP)

          ! Sum force over g in each Cartesian direction i
          ! ==>  f = sum_{g} i.g.e^{-ig.R}.V_H[n_Zc](g).den^{*}(g)
          loc_forces(:,atom) = loc_forces(:,atom) + &
               real(iG_vden(:,elements(atom)%pspecies_number)*eiGR,kind=DP)

       end do     ! End loop over atoms

       ! For g1/=0 slabs
       factor=2.0_DP

    end do
!$OMP END DO

    deallocate(iG_vden,stat=ierr)
    call utils_dealloc_check('paw_tcore_hartree_calc_forces','iG_vden',ierr)

!$OMP END PARALLEL

    ! Sum the result over all procs
    call comms_reduce('SUM',loc_forces)

    ! Deallocate
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('paw_tcore_hartree_calc_forces','recip',ierr)

    ! Scale
    loc_forces = loc_forces * 4.0_dp * PI / (grid%n1*grid%n2*grid%n3)

    ! If spin polarised, restore up spin density
    if (pub_num_spins == 2) den_slabs12(:,:,:,1) = den_slabs12(:,:,:,1) - &
         den_slabs12(:,:,:,2)

    ! Stop timer
    call timer_clock('paw_tcore_hartree_calc_forces',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_tcore_hartree_calc_forces'

    return
  end subroutine paw_tcore_hartree_calc_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_nlcc_calculate_forces_recip(density,tcore_density,nhat_den_grad, &
       grid,cell,elements,paw_sp,nlcc_forces,par,active_density,&
       active_nhat_den_grad)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the nonlinear core correction core charge in PAW.                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  density : input : ground state charge density on grid                  !
    !  tcore_density   : input  : NLCC core charge on grid                    !
    !  elements        : input  : list of elements and corresponding info     !
    !  grid            : input  : GRID_INFO type describing the grid          !
    !  nlcc_forces     : output : ionic forces due to NLCC corrections        !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in June 2010.                                  !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use cell_grid, only: cell_grid_recip_pt, GRID_INFO
    use comms, only: comms_reduce, pub_my_proc_id
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_aug_den_dim, pub_nhat_in_xc, pub_num_spins, pub_emft
!$  use rundat, only: pub_threads_max
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    use xc, only: xc_energy_potential,xc_embed_swap_functional

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in) :: par
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=dp), intent(inout) :: density(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(in) :: tcore_density(grid%ld1, &
         grid%ld2,grid%max_slabs12)
    real(kind=DP), intent(in) :: nhat_den_grad(grid%ld1, &
         grid%ld2,grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    real(kind=dp), intent(out) :: nlcc_forces(1:3,par%nat)
    ! jcap: if doing emft, active subsystem density
    real(kind=dp), intent(in), optional :: active_density(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    ! jcap: if doing emft and pub_aug, active subsystem nhat_den_grad
    real(kind=dp), optional, intent(in) :: active_nhat_den_grad(grid%ld1, &
         grid%ld2,grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)

    ! Local Variables
    integer :: ipt,i2,i3,islab23,is
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gvec(3)
    real(kind=DP) :: g_length
    real(kind=DP) :: factor
    real(kind=DP) :: xc_energy
    real(kind=DP) :: xc_energy_emft
    real(kind=DP) :: coreden, GdotR
    real(kind=DP), allocatable :: total_density(:,:,:,:)
    real(kind=DP), allocatable :: total_active_density(:,:,:,:)
    real(kind=DP), allocatable :: xc_pot(:,:,:,:)
    real(kind=DP), allocatable :: xc_pot_emft(:,:,:,:)
    complex(kind=DP), allocatable :: iG_coreden_vxc(:,:)
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: recip(:,:,:)

    logical :: loc_emft

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_nlcc_calculate_forces_recip'

    ! Start timer
    call timer_clock('paw_nlcc_calculate_forces_recip',1)

    ! jcap: check to see if we are doing an emft calculation and if
    ! this matches pub_emft
    loc_emft=present(active_density)
    call utils_assert(loc_emft.eqv.pub_emft,'Arguments of &
         &pseudo_nlcc_calculate_forces inconsistent with pub_emft:',pub_emft)

    ! Allocate
    allocate(total_density(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces_recip','total_density', &
         ierr)
    allocate(xc_pot(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins), &
         stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces_recip','xc_pot',ierr)
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces_recip','recip',ierr)
    if (loc_emft) then
       allocate(total_active_density(grid%ld1,grid%ld2,grid%max_slabs12, &
            pub_num_spins),stat=ierr)
       call utils_alloc_check('paw_nlcc_calculate_forces_recip', &
            'total_active_density',ierr)
       allocate(xc_pot_emft(grid%ld1,grid%ld2,grid%max_slabs12, &
            pub_num_spins),stat=ierr)
       call utils_alloc_check('paw_nlcc_calculate_forces_recip','xc_pot_emft', &
            ierr)
    end if

    ! Calculate total density
    factor = 1.0_DP / pub_num_spins
    do is=1,pub_num_spins
       total_density(:,:,:,is) = density(:,:,:,is) + tcore_density*factor
       if (loc_emft) total_active_density(:,:,:,is) = &
            active_density(:,:,:,is) + tcore_density*factor
    end do

    if (.not.pub_nhat_in_xc) then
       total_density = total_density - nhat_den_grad(:,:,:,:,0)
       if (loc_emft) total_active_density = total_active_density - &
            active_nhat_den_grad(:,:,:,:,0)
    end if

    ! Calculate the exchange correlation potential
    call xc_energy_potential(total_density,xc_energy,xc_pot,grid,cell,&
         pub_aug_den_dim,nhat_den_grad)

    ! jcap: If we are using embedded mean field theory, we need to
    ! call this twice more, in order to get the xc energy and
    ! potential for the active region
    if (loc_emft) then

       ! jcap: xc part of the hamiltonian is the only part that
       ! changes between different levels of DFT theory. First,
       ! calculate xc for low level theory for active subregion

       call xc_energy_potential(total_active_density, xc_energy_emft, &
               xc_pot_emft, grid, cell, pub_aug_den_dim, active_nhat_den_grad)

       ! jcap: subtract the contribution of the active subregion at
       ! the low level of theory
       xc_energy=xc_energy-xc_energy_emft
       xc_pot=xc_pot-xc_pot_emft

       ! jcap: need to swap to high level xc functional
       call xc_embed_swap_functional(.true.)

       ! jcap: calculate xc for high level theory for active
       ! subregion
       call xc_energy_potential(total_active_density, xc_energy_emft, &
            xc_pot_emft, grid, cell, pub_aug_den_dim, active_nhat_den_grad)

       ! jcap: need to swap the xc functional back
       call xc_embed_swap_functional(.false.)

       ! jcap: add the contribution of the active subregion at the
       ! high level of theory
       xc_energy=xc_energy+xc_energy_emft
       xc_pot=xc_pot+xc_pot_emft
    end if

    ! Average the spin channels in (:,:,:,1)
    if (pub_num_spins == 2) then
       xc_pot(:,:,:,1) = factor*(xc_pot(:,:,:,1) &
            + xc_pot(:,:,:,2))
    end if

    ! Initialise
    nlcc_forces = 0.0_DP

    ! Fourier transform the xc potential to reciprocal space
    call fourier_apply_cell_forward(xc_pot(:,:,:,1),recip,grid)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(i2,i3,gvec,g_length,species,coreden, &
!$OMP      iG_coreden_vxc,factor,islab23,atom,gdotR,eiGR,ierr) &
!$OMP SHARED (pub_my_proc_id,grid,paw_sp,elements,recip,par, &
!$OMP      pub_threads_max) &
!$OMP REDUCTION (+:nlcc_forces)

    allocate(iG_coreden_vxc(1:3,par%num_pspecies),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces_recip','iG_coreden_vxc', &
         ierr)

    ! For components g with symmetry points at -g
    factor=2.0_DP

!$OMP DO
    do ipt=1,grid%num_slabs23*grid%n2*grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1

       if ((islab23==1).and.(grid%proc_slab23(1) == pub_my_proc_id)) &
            factor = 1.0_DP

       ! Get g-vector and |g|
       call cell_grid_recip_pt(gvec,islab23 + &
            grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)
       g_length = sqrt(sum(gvec(1:3)**2))

       ! Loop over atoms to find n(G) for each
       do species=1,par%num_pspecies

          ! Check if we actually have a core charge for this species
          if (.not.paw_sp(species)%tcore_charge) cycle

          ! Interpolate value of core density at current g
          coreden = services_1d_interpolation( &
               paw_sp(species)%tcore_den_recip, &
               paw_sp(species)%n_recip_pts, &
               g_length*paw_sp(species)%inv_g_spacing,0)

          ! ndmh: calculate iG.n(G).vxc*(G) for each species and
          ! ndmh: each Cartesian direction
          iG_coreden_vxc(1:3,species) = factor * cmplx(0.0_DP,1.0_DP,kind=DP) &
               * gvec(1:3) * coreden * conjg(recip(i3,i2,islab23))

       end do

       ! ndmh: loop over atoms: this bit is O(N^2) so should be made as
       ! ndmh: short and simple as possible
       do atom=1,par%nat

          species = elements(atom)%pspecies_number

          ! Check if we actually have a core charge for this atom
          if (.not.paw_sp(species)%tcore_charge) cycle

          ! e^{-ig.R}
          gdotR = -(gvec(1)*elements(atom)%centre%x + &
               gvec(2)*elements(atom)%centre%y + &
               gvec(3)*elements(atom)%centre%z)
          eiGR = cmplx(COS(gdotR),SIN(gdotR),kind=DP)

          ! Sum force over g in each Cartesian direction i
          ! ==>  f = sum_{g} i.g.e^{-ig.R}.den_core(g).Vxc^{*}(g)
          nlcc_forces(:,atom) = nlcc_forces(:,atom) + &
               real(iG_coreden_vxc(:,species)*eiGR,kind=DP)

       enddo     ! End loop over atoms

       ! For g1/=0 slabs
       factor=2.0_DP

    end do
!$OMP END DO

    deallocate(iG_coreden_vxc,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces_recip', &
         'iG_coreden_vxc',ierr)

!$OMP END PARALLEL

    ! Sum the result over all procs
    call comms_reduce('SUM',nlcc_forces)

    ! Deallocate
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces_recip','recip',ierr)
    deallocate(xc_pot,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces_recip','xc_pot',ierr)
    deallocate(total_density,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces_recip', &
         'total_density',ierr)

    ! Scale factor
    nlcc_forces = nlcc_forces / real(grid%n1*grid%n2*grid%n3,dp)

    ! Stop timer
    call timer_clock('paw_nlcc_calculate_forces_recip',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_nlcc_calculate_forces_recip'

    return
  end subroutine paw_nlcc_calculate_forces_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_nlcc_calculate_forces(density,tcore_density,nhat_den_grad, &
       grid,cell,elements,paw_sp,nlcc_forces,par,active_density,&
       active_nhat_den_grad)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the nonlinear core correction core charge in PAW.                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  density : input : ground state charge density on grid                  !
    !  tcore_density   : input  : NLCC core charge on grid                    !
    !  elements        : input  : list of elements and corresponding info     !
    !  grid            : input  : GRID_INFO type describing the grid          !
    !  nlcc_forces     : output : ionic forces due to NLCC corrections        !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in June 2010.                                  !
    ! Modified to remove pub_par by Joseph Prentice, May 2018                 !
    !=========================================================================!

    use cell_grid, only: cell_grid_recip_pt, GRID_INFO, cell_grid_extract_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: comms_barrier, comms_reduce, pub_my_proc_id
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use geometry, only: POINT, OPERATOR(.dot.), geometry_magnitude
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_aug_den_dim, pub_nhat_in_xc, pub_num_spins, pub_emft
!$  use rundat, only: pub_threads_max, pub_threads_num_fftboxes
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block
    use spherical_wave, only: sw_init
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, utils_assert
    use xc, only: xc_energy_potential,xc_embed_swap_functional

    use rundat, only: pub_paw, pub_usp, pub_aug_funcs_recip, &
         pub_num_spins, pub_debug_on_root


    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in) :: par
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=dp), intent(inout) :: density(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(in) :: tcore_density(grid%ld1, &
         grid%ld2,grid%max_slabs12)
    real(kind=DP), intent(in) :: nhat_den_grad(grid%ld1, &
         grid%ld2,grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    real(kind=dp), intent(out) :: nlcc_forces(1:3,par%nat)
    ! jcap: if doing emft, active subsystem density
    real(kind=dp), intent(in), optional :: active_density(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    ! jcap: if doing emft and pub_aug, active subsystem nhat_den_grad
    real(kind=dp), optional, intent(in) :: active_nhat_den_grad(grid%ld1, &
         grid%ld2,grid%max_slabs12,pub_num_spins,0:pub_aug_den_dim)

    ! Local Variables
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gvec(3)
    real(kind=DP) :: g_length
    real(kind=DP) :: factor
    real(kind=DP) :: xc_energy
    real(kind=DP) :: xc_energy_emft
    real(kind=DP) :: coreden, GdotR
    real(kind=DP), allocatable :: total_density(:,:,:,:)
    real(kind=DP), allocatable :: total_active_density(:,:,:,:)
    real(kind=DP), allocatable :: xc_pot(:,:,:,:)
    real(kind=DP), allocatable :: xc_pot_emft(:,:,:,:)
    complex(kind=DP), allocatable :: iG_coreden_vxc(:,:)
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: recip(:,:,:)
    logical :: loc_emft
    integer :: loc_iat, orig_iat, iat
    integer :: isp, is, cart
    integer :: box_n1
    integer :: box_n2
    integer :: box_n3
    integer :: box_start1
    integer :: box_start2
    integer :: box_start3
    logical :: i_need_box
    real(kind=DP), allocatable :: buffer(:,:,:)
    real(kind=DP), allocatable :: xc_pot_box(:,:,:)
    real(kind=DP), allocatable :: atom_grad_core_den(:,:,:,:)
    complex(kind=DP), allocatable :: xc_pot_box_recip(:,:,:)
    complex(kind=DP), allocatable :: atom_grad_core_den_recip(:,:,:)
    integer :: igrid, npts
    real(kind=DP), allocatable :: atom_tcore(:,:,:)
    real(kind=DP) :: max_rad
    type(FFTBOX_INFO) :: core_den_box
    logical :: extend(3)
    integer :: pref(3)
    real(kind=DP) :: halo

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_nlcc_calculate_forces'

    ! Start timer
    call timer_clock('paw_nlcc_calculate_forces',1)

    ! Set up structures for core density box
    call paw_make_core_den_box(core_den_box,grid,cell,paw_sp,par)
    box_n1 = core_den_box%total_pt1
    box_n2 = core_den_box%total_pt2
    box_n3 = core_den_box%total_pt3

    ! jcap: check to see if we are doing an emft calculation and if
    ! this matches pub_emft
    loc_emft=present(active_density)
    call utils_assert(loc_emft.eqv.pub_emft,'Arguments of &
         &pseudo_nlcc_calculate_forces inconsistent with pub_emft:',pub_emft)

    ! Allocate
    allocate(total_density(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces','total_density',ierr)
    allocate(xc_pot(grid%ld1,grid%ld2,grid%max_slabs12,pub_num_spins), &
         stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces','xc_pot',ierr)
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces','recip',ierr)
    if (loc_emft) then
       allocate(total_active_density(grid%ld1,grid%ld2,grid%max_slabs12, &
            pub_num_spins),stat=ierr)
       call utils_alloc_check('paw_nlcc_calculate_forces','total_active_density',&
            ierr)
       allocate(xc_pot_emft(grid%ld1,grid%ld2,grid%max_slabs12, &
            pub_num_spins),stat=ierr)
       call utils_alloc_check('paw_nlcc_calculate_forces','xc_pot_emft',ierr)
    end if

    ! Calculate total density
    factor = 1.0_DP / pub_num_spins
    do is=1,pub_num_spins
       total_density(:,:,:,is) = density(:,:,:,is) + tcore_density*factor
       if (loc_emft) total_active_density(:,:,:,is) = &
            active_density(:,:,:,is) + tcore_density*factor
    end do

    if (.not.pub_nhat_in_xc) then
       total_density = total_density - nhat_den_grad(:,:,:,:,0)
       if (loc_emft) total_active_density = total_active_density - &
            active_nhat_den_grad(:,:,:,:,0)
    end if

    ! Calculate the exchange correlation potential
    call xc_energy_potential(total_density,xc_energy,xc_pot,grid,cell,&
         pub_aug_den_dim,nhat_den_grad)

    ! jcap: If we are using embedded mean field theory, we need to
    ! call this twice more, in order to get the xc energy and
    ! potential for the active region
    if (loc_emft) then

       ! jcap: xc part of the hamiltonian is the only part that
       ! changes between different levels of DFT theory. First,
       ! calculate xc for low level theory for active subregion

       call xc_energy_potential(total_active_density, xc_energy_emft, &
               xc_pot_emft, grid, cell, pub_aug_den_dim, active_nhat_den_grad)

       ! jcap: subtract the contribution of the active subregion at
       ! the low level of theory
       xc_energy=xc_energy-xc_energy_emft
       xc_pot=xc_pot-xc_pot_emft

       ! jcap: need to swap to high level xc functional
       call xc_embed_swap_functional(.true.)

       ! jcap: calculate xc for high level theory for active
       ! subregion
       call xc_energy_potential(total_active_density, xc_energy_emft, &
            xc_pot_emft, grid, cell, pub_aug_den_dim, active_nhat_den_grad)

       ! jcap: need to swap the xc functional back
       call xc_embed_swap_functional(.false.)

       ! jcap: add the contribution of the active subregion at the
       ! high level of theory
       xc_energy=xc_energy+xc_energy_emft
       xc_pot=xc_pot+xc_pot_emft
    end if

    ! Average the spin channels in (:,:,:,1)
    if (pub_num_spins == 2) then
       xc_pot(:,:,:,1) = factor*(xc_pot(:,:,:,1) &
            + xc_pot(:,:,:,2))
    end if

    ! Initialise
    nlcc_forces = 0.0_DP

!$OMP PARALLEL NUM_THREADS(pub_threads_num_fftboxes) DEFAULT(none) &
!$OMP PRIVATE(xc_pot_box,buffer,xc_pot_box_recip, &
!$OMP      atom_grad_core_den_recip,loc_iat,iat,isp,is,ierr,box_start1, &
!$OMP      box_start2,box_start3,orig_iat,i_need_box) &
!$OMP SHARED (pub_my_proc_id,par,grid,cell,paw_sp,xc_pot,core_den_box, &
!$OMP      box_n1,box_n2,box_n3,nlcc_forces,pub_threads_num_fftboxes,stdout)

    ! Allocate temporary arrays
    allocate(xc_pot_box(box_n1,box_n2,box_n3),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces','xc_pot_box',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces','buffer',ierr)
    allocate(xc_pot_box_recip(box_n1,box_n2,box_n3),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces', &
         'xc_pot_box_recip',ierr)
    allocate(atom_grad_core_den_recip(box_n1,box_n2,box_n3),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces', &
         'atom_grad_core_den_recip',ierr)

    ! Loop over atoms
!$OMP DO
    do loc_iat=1,par%max_atoms_on_proc

       ! Only need to extract if there is an atom left on this proc
       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id)) then

          ! Find atom number in input file order and species number
          iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          isp = par%elements_on_proc(loc_iat)%pspecies_number
          orig_iat = par%orig_atom(iat)

          if (paw_sp(isp)%tcore_charge) then
             ! Find where box for this atom is located in simulation cell
             call cell_grid_box_start_wrt_atom( &
                  box_start1, box_start2, box_start3, &
                  par%elements_on_proc(loc_iat)%centre, &
                  box_n1, box_n2, box_n3, grid, cell)

             i_need_box = .true.
          else
             ! No need for xc_pot box as no core charge is present
             i_need_box = .false.
          end if
       else
          i_need_box = .false.
          orig_iat = -1  ! suppress warning
       end if ! loc_iat < loc_nat

       ! Extract box of data from xc potential over simulation
       ! cell for this atom
!$OMP CRITICAL
       call cell_grid_extract_box(xc_pot_box(:,:,:), &
            buffer, xc_pot(:,:,:,1), grid, &
            box_n1, box_n2, box_n3, box_n1, box_n2, &
            box_start1, box_start2, box_start3, i_need_box, .false.)
!$OMP END CRITICAL

       ! Only need to calculate force if there is an atom left on this proc
       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id).and.i_need_box) then

          ! Fourier transform the xc pot box
          xc_pot_box_recip(:,:,:) = xc_pot_box(:,:,:)
          call fourier_apply_box('F','B',xc_pot_box_recip(:,:,:), &
               box=core_den_box)

          ! Calculate the NLCC force on this atom in recip space
          do cart=1,3

             call paw_tcore_on_box_recip(atom_grad_core_den_recip, box_start1, &
                  box_start2,box_start3,grid,cell,box_n1,box_n2,box_n3, &
                  par%elements_on_proc(loc_iat)%centre, &
                  isp,paw_sp,par%num_pspecies,cart)

             nlcc_forces(cart,orig_iat) = nlcc_forces(cart,orig_iat) - &
                  real(sum(xc_pot_box_recip(:,:,:) &
                  * atom_grad_core_den_recip(:,:,:)),kind=DP)
          end do

       end if  ! loc_iat < loc_nat
    end do  ! loc_iat
!$OMP END DO

    ! Deallocate temporary arrays
    deallocate(atom_grad_core_den_recip,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces', &
         'atom_grad_core_den_recip',ierr)
    deallocate(xc_pot_box_recip,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces', &
         'xc_pot_box_recip',ierr)
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces','buffer',ierr)
    deallocate(xc_pot_box,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces','xc_pot_box',ierr)

!$OMP END PARALLEL

    ! Reduce result across procs
    call comms_barrier
    call comms_reduce('SUM',nlcc_forces,3*par%nat)

    ! Destroy structures for core density box
    call paw_exit_core_den_box(core_den_box)

    ! Stop Timer
    call timer_clock('paw_nlcc_calculate_forces',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_nlcc_calculate_forces'

  end subroutine paw_nlcc_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_sphere_density_on_grid(density_on_grid,grid,cell,aug_box,rhoij, &
       ae_coeff,ps_coeff,paw_sp)

    !=========================================================================!
    ! This subroutine calculates the PAW sphere terms of the density on a     !
    ! regular real space grid (not used in the main calculation, but useful   !
    ! for analysing results afterwards).                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  density_on_grid (inout) : on input, the soft PS density (or zero)      !
    !                            on output, the contents on input plus the    !
    !                            sphere ae density times ae_coeff plus the    !
    !                            sphere ps density times ps_coeff             !
    !  grid (input) : GRID_INFO type describing grid                          !
    !  rhoij (input) : projector density kernel                               !
    !  ae_coeff (input) : coefficient of ae density                           !
    !  ps_coeff (input) : coefficient of ps density                           !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 09/06/2011                                    !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_deposit_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: pub_my_proc_id
    use fft_box, only: FFTBOX_INFO
    use fourier, only: fourier_apply_box
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins
    use services, only: services_radial_transform, services_radial_integral_rmax
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_get_block, sparse_get_par
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_assert

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(GRID_INFO), intent(in) :: grid
    type(FFTBOX_INFO), intent(in) :: aug_box
    real(kind=DP), intent(inout) :: density_on_grid(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    type(SPAM3), intent(in) :: rhoij(pub_num_spins)
    real(kind=DP), intent(in) :: ae_coeff, ps_coeff
    type(PAW_SPECIES), intent(in) :: paw_sp(:)

    ! Local Variables
    real(kind=DP),allocatable :: density_box(:,:,:,:)
    real(kind=DP),allocatable :: buffer(:,:,:)
    real(kind=DP),allocatable :: rho_ij_block(:,:,:)
    real(kind=DP),allocatable :: ae_den_rad_LM(:,:,:)
    real(kind=DP),allocatable :: ps_den_rad_LM(:,:,:)
    real(kind=DP),allocatable :: den_work(:)
    integer :: iat, loc_iat
    integer :: ierr
    integer :: is
    integer :: isp
    integer :: npw, nspins
    integer :: box_n1,box_n2,box_n3
    integer :: box_start1,box_start2,box_start3
    integer :: nLM_max, npts_max, npts_recip, npts_phi
    integer :: igrid
    integer :: lup, mup, iLM
    logical :: i_have_box
    type(PARAL_INFO), pointer :: par

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_sphere_density_on_grid'

    ! Start timer
    call timer_clock('paw_sphere_density_on_grid',1)

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, rhoij(1))
    call utils_assert(par%num_pspecies == size(paw_sp), 'Error in &
         &paw_sphere_density_on_grid: allocated parallel strategy is &
         &incompatible with paw_sp.')

    ! Find sizes of box, shorthands, grid sizes tec
    box_n1 = aug_box%total_pt1
    box_n2 = aug_box%total_pt2
    box_n3 = aug_box%total_pt3
    npw = maxval(paw_sp(:)%npw_tot)
    nspins = pub_num_spins
    npts_max = 0
    nLM_max = 0
    npts_recip = 0
    do isp=1,par%num_pspecies
       if (.not.any(par%elements_on_proc(:)%pspecies_number==isp)) cycle
       npts_max = max(npts_max,paw_sp(isp)%grid(paw_sp(isp)%phi_grid)%npt)
       npts_max = max(npts_max,paw_sp(isp)%grid(paw_sp(isp)%core_den_grid)%npt)
       nLM_max = max(nLM_max,(paw_sp(isp)%lmax+1)**2)
       npts_recip = max(npts_recip,paw_sp(isp)%n_recip_pts)
    end do

    ! Allocate workspace
    allocate(rho_ij_block(npw,npw,nspins),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','rho_ij_block',ierr)
    allocate(density_box(box_n1,box_n2,box_n3,nspins),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','density_box',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','buffer',ierr)
    allocate(ae_den_rad_LM(npts_max,nspins,nLM_max),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','ae_den_rad_LM',ierr)
    allocate(ps_den_rad_LM(npts_max,nspins,nLM_max),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','ps_den_rad_LM',ierr)
    allocate(den_work(npts_max),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','den_work',ierr)

    do loc_iat=1,par%max_atoms_on_proc

       if (loc_iat<=par%num_atoms_on_proc(pub_my_proc_id)) then
          iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
          isp = par%elements_on_proc(loc_iat)%pspecies_number
          igrid = paw_sp(isp)%phi_grid
          npts_phi = paw_sp(isp)%grid(igrid)%npt

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               par%elements_on_proc(loc_iat)%centre, box_n1, box_n2, box_n3, &
               grid, cell)

          ! Get block of \rho_ij for this atom
          do is=1,nspins
             call sparse_get_block(rho_ij_block(:,:,is),rhoij(is),iat,iat)
          end do

          ! Generate the ae and ps densities for this atom in the augmentation
          ! box using rhoij
          call paw_ae_density_rad_LM(ae_den_rad_LM,isp,rho_ij_block,npw, &
               npts_max,nspins,nLM_max,paw_sp,den_work)
          call paw_ps_density_rad_LM(ps_den_rad_LM,isp,rho_ij_block,npw, &
               npts_max,nspins,nLM_max,paw_sp,den_work)

          density_box(:,:,:,:) = 0.0_DP
          do lup=0,paw_sp(isp)%lmax
             do mup=-lup,lup
                iLM = lup*lup + lup + 1 + mup

                do is=1,nspins

                   ! Create density to be created, from ae and ps coefficients
                   den_work(:) = ae_coeff * ae_den_rad_LM(:,is,iLM) + &
                        ps_coeff * ps_den_rad_LM(:,is,iLM)

                   call paw_den_to_box_real(density_box(:,:,:,is),cell,den_work, &
                        paw_sp(isp)%grid(igrid)%r,paw_sp(isp)%rcut,npts_phi,&
                        lup,mup,box_start1,box_start2,box_start3,grid,box_n1, &
                        box_n2,box_n3,par%elements_on_proc(loc_iat)%centre)
                end do  ! is
             end do  ! mup
          end do  ! lup

          i_have_box = .true.
       else
          ! Nothing to deposit on this proc
          i_have_box = .false.
       end if

       ! Deposit this box to the simulation cell if present, or just wait for
       ! data from other procs if no box
       do is=1,pub_num_spins
          call cell_grid_deposit_box(density_on_grid(:,:,:,is), &
               density_box(:,:,:,is), buffer, grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               box_start1, box_start2, box_start3, i_have_box, .false.)
       end do
    end do

    ! Deallocate
    deallocate(den_work,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','den_work',ierr)
    deallocate(ps_den_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','ps_den_rad_LM',ierr)
    deallocate(ae_den_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','ae_den_rad_LM',ierr)
    deallocate(density_box,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','density_box',ierr)
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','buffer',ierr)
    deallocate(rho_ij_block,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','rho_ij_block',ierr)
    nullify(par)

    ! Stop timer
    call timer_clock('paw_sphere_density_on_grid',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_sphere_density_on_grid'

    return
  end subroutine paw_sphere_density_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_den_to_box_real(box,cell,den_lm,r,rcut,npts,lup,mup, &
      cell_start1,cell_start2,cell_start3,grid,box_n1,box_n2,box_n3, &
      atom_origin)

    !==================================================================!
    ! This subroutine calculates a l,m component of the density in     !
    ! real space within a sphere of radius rcut, provided on a grid    !
    ! of points r(1:npts), where the density is given by               !
    !   n_LM(\vec r) = n_LM(r) S_LM(\vec r)                            !
    ! It does so on the simulation cell fine grid in a box centered on !
    ! the atom.                                                        !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  box_(out) :                                                     !
    !  den_lm (in) :                                                   !
    !  r (in) :                                                        !
    !  rcut (in) :                                                     !
    !  npts (in) :                                                     !
    !  lup (in) :                                                      !
    !  mup (in) :                                                      !
    !  cell_start1 (in) :                                              !
    !  cell_start2 (in) :                                              !
    !  cell_start3 (in) :                                              !
    !  grid (in) :                                                     !
    !  box_n1 (in) :                                                   !
    !  box_n2 (in) :                                                   !
    !  box_n3 (in) :                                                   !
    !  atom_origin (in) :                                              !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 09/05/11.           !
    !==================================================================!

    use basis, only: basis_box_origin_to_atom
    use cell_grid, only: GRID_INFO
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use services, only: services_locate_interp,services_linear_interpolation
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_real_sph_harm, sw_init

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    real(kind=DP),intent(inout) :: box(box_n1,box_n2,box_n3)
    integer,intent(in) :: npts
    real(kind=DP), intent(in) :: den_lm(npts)
    real(kind=DP), intent(in) :: r(npts), rcut
    integer,intent(in) :: lup
    integer,intent(in) :: mup
    integer,intent(in) :: cell_start1
    integer,intent(in) :: cell_start2
    integer,intent(in) :: cell_start3
    type(GRID_INFO),intent(in) :: grid
    type(POINT),intent(in) :: atom_origin

    ! Local Variables
    integer :: i1,i2,i3
    integer :: ipt
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    real(kind=DP) :: denval
    real(kind=DP) :: slmval
    type(POINT) :: r_cell, r_sphere
    real(kind=DP) :: rmag

    if (lup>4) call sw_init(lup+1,1)

    call basis_box_origin_to_atom(box_to_atom,box_origin,atom_origin, &
         cell_start1,cell_start2,cell_start3,grid%da1,grid%da2,grid%da3, &
         cell)

    do i3=1,box_n3
       do i2=1,box_n2
          do i1=1,box_n1

             r_cell = box_origin + (i1-1)*grid%da1 &
                                 + (i2-1)*grid%da2 &
                                 + (i3-1)*grid%da3

             r_sphere = r_cell - atom_origin
             rmag = geometry_magnitude(r_sphere)

             if (rmag>rcut) then
                box(i1,i2,i3) = 0.0_DP
                cycle
             end if

             ipt = services_locate_interp(rmag,r,npts)
             denval = services_linear_interpolation(rmag,den_lm(ipt), &
                  den_lm(ipt+1),r(ipt),r(ipt+1))

             ! multiply value of shapefunc at each g-point by the appropriate
             ! phase factor (n.b. real part and imaginary part stored as
             ! separate consecutive elements in the x-direction of the array),
             ! and the appropriate real spherical harmonic factor.
             slmval = sw_real_sph_harm(r_sphere%x, &
                  r_sphere%y,r_sphere%z,rmag,lup,mup)

             box(i1,i2,i3) = box(i1,i2,i3) + &
                  denval * slmval

          end do
       end do
    end do

  end subroutine paw_den_to_box_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine paw_FO_nlcc_energy(xc_pot,grid,cell,elements,paw_sp,par, &
       nlcc_energy,directions1,weights1,directions2,weights2)

    ! THIS SUBROUTINE COULD USE IMPROVEMENTS: ONLY WORK IN THE AUGMENTATION
    ! BOX OF THE SELECTED ATOM, MAYBE?

    use cell_grid, only: cell_grid_recip_pt, GRID_INFO
    use comms, only: comms_reduce, pub_my_proc_id
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins
!$  use rundat, only: pub_threads_max
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    real(kind=dp), intent(in) :: xc_pot(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    type(PARAL_INFO), intent(in) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(PAW_SPECIES), intent(in) :: paw_sp(par%num_pspecies)
    real(kind=dp), intent(out) :: nlcc_energy
    integer, intent(in) :: directions1(par%nat)
    real(kind=DP), intent(in) :: weights1(par%nat)
    integer, intent(in), optional :: directions2(par%nat)
    real(kind=DP), intent(in), optional :: weights2(par%nat)

    ! Local Variables
    integer :: ipt,i2,i3,islab23,is
    integer :: species, atom
    integer :: ierr
    real(kind=DP) :: gvec(3)
    real(kind=DP) :: g_length
    real(kind=DP) :: factor
    real(kind=DP) :: coreden, GdotR
    complex(kind=DP) :: iG_coreden_vxc
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: recip(:,:,:)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_FO_nlcc_energy'

    ! Start timer
    call timer_clock('paw_FO_nlcc_energy',1)

    ! Allocate
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_FO_nlcc_energy','recip',ierr)

!    factor = 1.0_DP / pub_num_spins
!    ! Average the spin channels in (:,:,:,1)
!    if (pub_num_spins == 2) then
!       xc_pot(:,:,:,1) = factor*(xc_pot(:,:,:,1) + xc_pot(:,:,:,2))
!    end if
!
!    gcc32: neglect for the moment

    ! Initialise
    nlcc_energy = 0.0_DP

    ! Fourier transform the xc potential to reciprocal space
    call fourier_apply_cell_forward(xc_pot(:,:,:,1),recip,grid)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(i2,i3,gvec,g_length,coreden, &
!$OMP      iG_coreden_vxc,factor,islab23,gdotR,eiGR,ierr,atom,species) &
!$OMP SHARED (pub_my_proc_id,grid,paw_sp,elements,recip,par, &
!$OMP      pub_threads_max,directions1,weights1,directions2,weights2) &
!$OMP REDUCTION (+:nlcc_energy)

    ! For components g with symmetry points at -g
    factor=2.0_DP

!$OMP DO
    do ipt=1,grid%num_slabs23*grid%n2*grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1

       if ((islab23==1).and.(grid%proc_slab23(1) == pub_my_proc_id)) &
            factor = 1.0_DP

       ! Get g-vector and |g|
       call cell_grid_recip_pt(gvec,islab23 + &
            grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)
       g_length = sqrt(sum(gvec(1:3)**2))

       do species = 1, par%num_pspecies

          ! Check if we actually have a core charge for this species
          if (.not.paw_sp(species)%tcore_charge) cycle

          ! Interpolate value of core density at current g
          coreden = services_1d_interpolation( &
               paw_sp(species)%tcore_den_recip, &
               paw_sp(species)%n_recip_pts, &
               g_length*paw_sp(species)%inv_g_spacing,0)

          ! calculate -iG_dir1.n(G).vxc*(G) if FO (i.e. direction2 NOT present)
          ! calculate -G_dir1.G_dir2.n(G).vxc*(G) if SO (direction2 present)

          do atom = 1, par%nat

             if (species == elements(atom)%pspecies_number) then

                if (.not.present(directions2)) then
                   iG_coreden_vxc = factor * cmplx(0.0_DP,-1.0_DP,kind=DP) * &
                        gvec(directions1(atom)) * coreden * &
                        conjg(recip(i3,i2,islab23)) * weights1(atom)
                else
                   iG_coreden_vxc = cmplx(-1.0_DP,0.0_DP,kind=DP) * factor * &
                        gvec(directions1(atom)) * gvec(directions2(atom)) * &
                        weights1(atom) * weights2(atom) * coreden * &
                        conjg(recip(i3,i2,islab23))
                end if

                ! e^{-ig.R}
                gdotR = -(gvec(1)*elements(atom)%centre%x + &
                     gvec(2)*elements(atom)%centre%y + &
                     gvec(3)*elements(atom)%centre%z)
                eiGR = cmplx(COS(gdotR),SIN(gdotR),kind=DP)

                ! Sum force over g in each Cartesian direction i
                ! ==>  f = sum_{g} i.g.e^{-ig.R}.den_core(g).Vxc^{*}(g)
                nlcc_energy = nlcc_energy + real(iG_coreden_vxc*eiGR,kind=DP)

             end if

          end do ! over atoms of that species

       end do ! over species

       ! For g1/=0 slabs
       factor = 2.0_DP

    end do
!$OMP END DO

!$OMP END PARALLEL

    call comms_reduce('SUM',nlcc_energy)

    ! Deallocate
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('paw_FO_nlcc_energy','recip',ierr)

    ! Scale factor
    nlcc_energy = nlcc_energy / real(grid%n1*grid%n2*grid%n3,dp)

    ! Stop timer
    call timer_clock('paw_FO_nlcc_energy',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_FO_nlcc_energy'

    return
  end subroutine paw_FO_nlcc_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module paw
