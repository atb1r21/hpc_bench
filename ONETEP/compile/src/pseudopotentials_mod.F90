! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Arash A. Mostofi, Nicholas D.M. Hine, and
!   Peter Haynes
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!   Subsequent additions and modifications by Laura Ratcliff, Jacek Dziedzic,
!   Gabriel Constantinescu, Andrea Greco, and Jose M Escartin.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module pseudopotentials

  use constants, only : DP, PI, stdout
  use geometry, only : POINT
  use projectors, only : PROJECTOR_SET
  use rundat, only: pub_debug_on_root

  implicit none

  private

  ! cks, 22/1/2004: type for containing pseudopotential information
  type PSEUDO_SPECIES

     ! cks: name of pseudopotential file
     character(len=64) :: pseudo_name

     ! cks: atomic number of atom
     integer :: atomic_number

     ! cks: charge of "pseudised" core of atom
     real(kind=DP) :: ion_charge   ! pa: changed from integer to real

     ! ndmh: whether to subtract the coulomb potential -Z/r
     logical :: subtract_coul

     ! cks: number of points in the radial grid
     integer :: n_rad_pts

     ! ndmh: spacing of points in the radial grid
     real(kind=DP) :: inv_g_spacing

     ! cks: number of angular momentum shells of projectors
     integer :: n_shells

     ! ndmh: number of projectors in total
     integer :: n_proj

     ! cks: angular momentum of each shell, this
     ! cks: array is allocated to be ang_mom(n_shells)
     integer, allocatable, dimension(:) :: ang_mom

     ! cks: maximum radial reciprocal k-vector up to which
     ! cks: the potential and projectors are defined
     real(kind=DP) :: ps_gmax

     ! cks: core radius for each projector shell,
     ! cks: core_radius(n_shells)
     real(kind=DP), allocatable, dimension(:) :: core_radius

     ! ndmh: Kleinman-Bylander energies for each pair of projector
     ! ndmh: shells: D0(n_shells,n_shells)
     real(kind=DP), allocatable, dimension(:,:) :: D0

     ! cks: the radial recip part of the local potential
     ! cks: rad_locpot_recip(n_rad_pts)
     real(kind=DP), allocatable, dimension(:) :: rad_locpot_recip

     ! cks: the values of each projector (shell) on the radial
     ! cks: grid. The last column of this array contains the
     ! cks: grid points, so the size of the array is
     ! cks: rad_proj_recip(n_rad_pts, n_shells +1 )
     real(kind=DP), allocatable, dimension(:,:) :: rad_proj_recip

     ! ndmh: whether a core charge is present for this species
     logical :: core_charge

     ! ndmh: if core charge is present, array to store core charge in
     ! ndmh: reciprocal space
     real(kind=DP), allocatable, dimension(:) :: core_charge_recip

     ! ndmh: the augmentation charges for each shell
     real(kind=DP), allocatable, dimension(:,:) :: aug_q

     ! ndmh: the L-independent function
     integer :: kkbeta, mesh
     real(kind=DP), allocatable, dimension(:) :: rlog
     real(kind=DP), allocatable, dimension(:) :: rab
     real(kind=DP), allocatable, dimension(:,:,:) :: qfunc
     real(kind=DP), allocatable, dimension(:) :: rinner

     ! ndmh: the L-dependent coefficients
     integer :: nqfcoef
     integer :: qf_lmax
     real(kind=DP), allocatable, dimension(:,:,:,:) :: qfcoef

     ! jd: --- Internal state of the local openbc pseudo ---
     logical :: openbc_locps_initialised = .false.
     ! jd: Vloc(x) on radial grid, for all species
     real(kind=DP), allocatable :: vloc_lookup(:)
     ! jd: d/dx(Vloc(x)) on radial grid, for all species,
     !     only filled if pub_forces_needed
     real(kind=DP), allocatable :: vlocder_lookup(:)
     ! jd: -------------------------------------------------

  end type PSEUDO_SPECIES

  public :: PSEUDO_SPECIES

  ! ndmh: subroutines relating to pseudopotential initialisation/exit
  public :: pseudopotentials_read_species
  public :: pseudo_species_init_proj
  public :: pseudopotentials_species_exit

  ! ndmh: routines for retrieving the pseudopotential for use by atom_mod
  public :: pseudo_get_locpot_rad
  public :: pseudo_get_core_den_rad
  public :: pseudo_get_projectors_q
  public :: pseudo_get_aug_funcs
  public :: pseudo_get_projector_info

  ! ndmh: subroutines relating to construction of local potential and forces
  public :: pseudo_make_structure_factor
  public :: pseudopotentials_local_on_grid
  public :: pseudopotentials_core_density
  public :: pseudopotentials_FO_nlcc_energy ! gcc32
  public :: pseudo_local_calculate_forces
  public :: pseudo_nlcc_calculate_forces

  ! ndmh: subroutines relating to nonlocal potentials
  public :: pseudo_get_dij
  public :: pseudopotentials_nonlocal_mat
  public :: pseudo_nl_calculate_forces
  public :: pseudo_aug_Q_matrix
  public :: pseudo_atom_aug_den
  public :: pseudo_atom_aug_integrals
  public :: pseudo_atom_aug_force
  public :: pseudo_nonlocal_commutator_mat
  public :: pseudo_nonlocal_com_mat_fd
  ! lr408 - for debugging purposes only
  public :: pseudo_nonlocal_com_mat_direct

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SUBROUTINES FOR PSEUDOPOTENTIAL INITIALISATION !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_read_species(p_species,elements,cell,species, &
          par,any_nl_proj)

    !=================================================================!
    ! This subroutine allocates and initialises the p_species array   !
    ! whose elements are PSEUDO_SPECIES types - one for each distinct !
    ! ionic species. It also sets the par%num_projectors value.  !
    ! The norm-conserving pseudopotential for each ionic species is   !
    ! read from a file in CASTEP .recpot format.                      !
    !-----------------------------------------------------------------!
    ! WARNING: The fftbox_proj_recip corresponding to each p_species  !
    !          are not initialised by this subroutine. They are       !
    !          initialised by a separate call to the                  !
    !          pseudo_species_init_proj subroutine.                   !
    !-----------------------------------------------------------------!
    ! Written  by Chris-Kriton Skylaris on 23/1/2004.                 !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                       !
    ! Modified to allow using this routine on a system-wide basis in  !
    ! embedding calculations too by Joseph Prentice, September 2018   !
    !=================================================================!

#ifdef ACCELRYS
    use comms, only: pub_my_proc_id
#endif
    use gaunt_coeff, only: gaunt_init
    use ion, only: ELEMENT, SPECIE
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_any_nl_proj, pub_nhat_in_xc, pub_usp
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in)   :: cell
    type(ELEMENT), intent(inout) :: elements(:)!par%nat)
    type(PSEUDO_SPECIES), intent(inout), pointer :: p_species(:)
    type(SPECIE), intent(in) :: species(:)
    type(PARAL_INFO), intent(inout) :: par
    logical, intent(inout)          :: any_nl_proj ! rc2013: any nl's in this rep?
    ! Local Variables
    integer :: ierr                ! pdh: Error flag
    integer :: species_counter
    integer :: sp
    integer :: atom
    integer :: shell
    integer :: mom
    integer :: p_number
    integer :: lmax
    integer :: nat
    logical :: we_have_it
    type(ELEMENT), allocatable :: element_buffer(:)

    ! jcap: get number of atoms from size of elements array, not par
    nat=size(elements)

    ! Allocate workspace
    allocate(element_buffer(nat),stat=ierr)
    call utils_alloc_check('pseudopotentials_read_species', &
         'element_buffer',ierr)

    species_counter = 1
    element_buffer(1)%species_number = elements(1)%species_number
    element_buffer(1)%atomic_number = elements(1)%atomic_number

    ! cks: loop over all atoms as they appear in the input file
    ! cks: and fill element_buffer with information about each
    ! cks: distinct species.
    do atom=2,nat

       we_have_it = .false.

       ! cks: loop over all species found so far and check if we have this
       !      species
       species_loop: do sp=1,species_counter
          if (species(element_buffer(sp)%species_number)%pseudo_name == &
                species(elements(atom)%species_number)%pseudo_name) then
             we_have_it = .true.
             exit species_loop
          end if
       end do species_loop

       ! cks: if we don't have it, let's get it!
       if (.not. we_have_it) then
          species_counter = species_counter + 1
          element_buffer(species_counter)%species_number = &
               elements(atom)%species_number
          element_buffer(species_counter)%atomic_number = &
               elements(atom)%atomic_number
       end if

    end do

    ! cks: allocate p_species type(pseudo_species) array
    if (.not.associated(p_species)) then
       allocate(p_species(species_counter),stat=ierr)
       call utils_alloc_check('pseudopotentials_read_species','p_species',ierr)
    else
       call utils_abort('Error in pseudopotentials_read_species: &
            &p_species already allocated')
    end if

    ! cks: do some p_species initialisation from the element_buffer
    do sp=1,species_counter
       p_species(sp)%pseudo_name = &
            species(element_buffer(sp)%species_number)%pseudo_name
       p_species(sp)%atomic_number = element_buffer(sp)%atomic_number
    end do

    deallocate(element_buffer,stat=ierr)
    call utils_dealloc_check('pseudopotentials_read_species', &
        'element_buffer',ierr)

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! cks: Now connect the elements array to the species array
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    elements(:)%pspecies_number = -1
    do atom=1,nat

       ! cks: loop over all species found so far and check if we have this
       !      species
       number_set_loop: do sp=1,species_counter
          if (p_species(sp)%pseudo_name == &
             species(elements(atom)%species_number)%pseudo_name) then
             elements(atom)%pspecies_number = sp
             exit number_set_loop
          end if
       end do number_set_loop

       ! cks: sanity check
       if (elements(atom)%pspecies_number == -1) then
          call utils_abort('Error in pseudopotentials_read_species: &
               &elements%pspecies_number = -1 for atom', atom)
       end if

    end do
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! cks: end connect the elements array to the species array
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    ! cks: allocate memory for the allocatable components of p_species array
    call pseudo_species_alloc(p_species)

    ! cks: complete the initialisation of the p_species array
    call pseudo_species_read_radial(p_species,cell)

    ! cks: find the total number of projectors
    ! jcap: these are assigned even when this is called for the whole
    ! subsystem - this will adjust regions(1)%par, and this routine
    ! should be called again AFTERWARDS for this region to get the
    ! correct values for that region
    par%num_projectors = 0
    par%num_pawpws = 0
    do atom=1,nat

#ifdef ACCELRYS
       ! ndmh: protection against dodgy compilers
       if(atom>size(elements))then
          call utils_abort('ERROR in pseudopotentials_read_species on proc i. &
               &Atom j does not fit element array (size k). Values of i, j, k &
               &follow.', pub_my_proc_id,atom,size(elements))
       endif
#endif

       p_number = elements(atom)%pspecies_number
       do shell=1,p_species(p_number)%n_shells

          mom = p_species(p_number)%ang_mom(shell)
          par%num_projectors = par%num_projectors + 2*mom+1

       end do
    end do

    ! ndmh: set flag to denote presence of nl projectors
    ! ndmh 04/01/10: only if core radius of any species is > 0
    ! rc2013: there may be projectors outside this region
    !pub_any_nl_proj = .false.
    any_nl_proj=.false.
    if (par%num_projectors > 0) then
       ! jcap: this will be incorrect for the whole system, but this
       ! is fine because we then check for nl projs in all regions
       ! later
       do sp=1,par%num_pspecies
          if (any(p_species(sp)%core_radius(:)>0.0_DP)) then
             any_nl_proj = .true.
          end if
       end do
    end if

    if (pub_usp) then

       ! Count highest angular momentum over all species
       lmax = 0
       do sp=1,size(p_species)
          lmax=max(p_species(sp)%qf_lmax,lmax)
       end do

       ! Initialise Gaunt Coefficients
       call gaunt_init(lmax+1)

    end if

    ! ndmh: set flag to include augmentation charge in XC calculation
    pub_nhat_in_xc = .true.

    ! cks: now set some redundant information in the elements
    ! cks: structure which we need for the time being but should
    ! cks: be removed in the future and be only available in elements
    do atom=1,nat

       elements(atom)%nprojectors     = 0
       elements(atom)%npawpws         = 0
       elements(atom)%ncorewfs        = 0
       elements(atom)%max_core_radius = 0.0_DP
       elements(atom)%max_core_wf_radius = 0.0_DP
       p_number = elements(atom)%pspecies_number

       elements(atom)%ion_charge = p_species(p_number)%ion_charge
       do shell=1,p_species(p_number)%n_shells

          mom = p_species(p_number)%ang_mom(shell)

          elements(atom)%nprojectors = elements(atom)%nprojectors + 2*mom+1

          elements(atom)%max_core_radius = &
               max(elements(atom)%max_core_radius, &
               p_species(p_number)%core_radius(shell))

       end do

    end do

  end subroutine pseudopotentials_read_species


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_species_alloc(p_species)

    !============================================================!
    ! This subroutine reads the pseudopotentials for each atomic !
    ! species and allocates the necessary memory for them as     !
    ! elements of the p_species array.                           !
    !------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/1/2004.             !
    ! Tidied up by Peter Haynes on 1/7/2004.                     !
    ! Modified to read multiple file formats (recpots and usps)  !
    ! by Nicholas Hine on 29/01/2009.                            !
    ! par removed as input by Joseph Prentice, June 2019         !
    !============================================================!

    use comms, only: comms_bcast, pub_on_root, pub_root_proc_id, &
                     pub_image_comm, pub_imroots_comm, comms_barrier
    use image_comms, only: pub_my_image
    use rundat, only: pub_nlcc, pub_xc_functional, pub_aug, pub_usp, &
                      pub_num_images
    use utils, only: utils_alloc_check, utils_abort

    implicit none

    ! Arguments
    type(PSEUDO_SPECIES), intent(inout) :: p_species(:)

    ! Local Variables
    integer, parameter  :: cfac=4
    integer             :: ierr
    integer             :: ps
    integer             :: tot_num_points
    integer             :: tot_num_projectors
    integer             :: num_projectors
    integer             :: n_shells
    integer             :: i
    integer             :: img
    !character(len=32)   :: string
    character(len=80)   :: line
    character(len=64)   :: current_file
    character(len=80)   :: current_job
    character(len=10)   :: psp_xc_functional
    logical             :: xc_mismatch
    logical             :: count_points
    real(kind=DP)       :: gmax
    real(kind=DP)       :: temp
    real(kind=DP)       :: fact1
    real(kind=DP)       :: ionic_charge
    integer             :: num_pspecies

    ! ndmh: determine if any species present contain core charges for NLCC
    pub_nlcc = .false.

    ! ndmh: Charge augmentation will not be active unless any species with
    ! ndmh: nonzero augmentation charges are found.
    pub_usp = .false.
    pub_aug = .false.

    ! rc2013: get number of elements from p_species
    num_pspecies = size(p_species)

    if (pub_on_root) then

       do ps=1,num_pspecies

          ! pdh: get file name for this species
          current_file = p_species(ps)%pseudo_name

          ! ndmh: determine what file format we are reading and scan it
          ! ndmh: to find sizes of various arrays

          ! kkbd: simple mutex on the file so that multiple
          !       images aren't manipulating a file simultaneously
          do img=0, pub_num_images-1
             if ((pub_num_images == 1).or.(img == pub_my_image)) then
                ! ndmh: Old-style CASTEP .recpot
                if (index(current_file,'recpot') > 0) then

                   call internal_scan_recpot

                ! ndmh: New-style CASTEP OTFG .usp (not necessarily ultrasoft!)
                else if (index(current_file,'usp') > 0) then

                   call internal_scan_usp

                else
                    call utils_abort("Error in pseudopotentials_read_species: &
                         &Unrecognised file format in  '"//trim(current_file)//"'&
                         &. NB: This can also happen if you forget to put quotes &
                         &around a filename that contains slashes.")
                end if
             end if
             if (pub_num_images > 1) call comms_barrier(pub_imroots_comm)
          end do ! image mutex

          ! ndmh: check if pseudopotential file contains a reference to a
          ! ndmh: different XC functional from the one we are using
          xc_mismatch = .false.
          select case (pub_xc_functional)
          case ('LDA','CAPZ')
             if (psp_xc_functional=='GGA') xc_mismatch = .true.
          case ('GGA','PBE','PW91','RPBE','WC')
             if (psp_xc_functional=='LDA') xc_mismatch = .true.
          end select

          ! ndmh: report if there is a mismatch
          if (xc_mismatch.and.pub_on_root) then
             write(stdout,'(3a)')'WARNING in pseudopotentials_read_species: &
                  &string "',trim(psp_xc_functional),'" found in pseudopotential'
             write(stdout,'(5a)')'file "',trim(current_file), &
                  '", yet xc_functional = "',trim(pub_xc_functional),'".'
             write(stdout,*)
          end if

          ! ndmh: check for NLCC core charge
          if (p_species(ps)%core_charge) pub_nlcc = .true.

       end do

    endif

    ! ndmh: let all procs know whether nonlinear core corrections are present
    call comms_bcast(pub_root_proc_id,pub_nlcc)

    ! cks: let all procs know what pub_root_proc_id just got
    do ps=1,num_pspecies
       call comms_bcast(pub_root_proc_id,p_species(ps)%n_shells)
       call comms_bcast(pub_root_proc_id,p_species(ps)%n_proj)
       call comms_bcast(pub_root_proc_id,p_species(ps)%n_rad_pts)
       call comms_bcast(pub_root_proc_id,p_species(ps)%ion_charge)
       call comms_bcast(pub_root_proc_id,p_species(ps)%ps_gmax)
       call comms_bcast(pub_root_proc_id,p_species(ps)%inv_g_spacing)
       call comms_bcast(pub_root_proc_id,p_species(ps)%core_charge)
       call comms_bcast(pub_root_proc_id,p_species(ps)%subtract_coul)
       call comms_bcast(pub_root_proc_id,p_species(ps)%kkbeta)
       call comms_bcast(pub_root_proc_id,p_species(ps)%mesh)
       call comms_bcast(pub_root_proc_id,p_species(ps)%nqfcoef)
       call comms_bcast(pub_root_proc_id,p_species(ps)%qf_lmax)
    end do


    ! cks: now allocate memory for all allocatable components
    do ps=1,num_pspecies

       allocate(p_species(ps)%ang_mom(p_species(ps)%n_shells),stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%ang_mom',ierr)

       allocate(p_species(ps)%core_radius(p_species(ps)%n_shells),stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%core_radius',ierr)

       allocate(p_species(ps)%D0(p_species(ps)%n_shells, &
            p_species(ps)%n_shells),stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%D0',ierr)

       allocate(p_species(ps)%rad_locpot_recip(p_species(ps)%n_rad_pts), &
            stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%rad_locpot_recip',ierr)

       allocate(p_species(ps)%rad_proj_recip(p_species(ps)%n_rad_pts, &
            p_species(ps)%n_shells+1),stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%rad_proj_recip',ierr)

       ! ndmh: allocate memory for core charge if present
       if (p_species(ps)%core_charge) then
          allocate(p_species(ps)%core_charge_recip(p_species(ps)%n_rad_pts), &
               stat=ierr)
          call utils_alloc_check('pseudo_species_alloc','p_species(ps)%core_charge_recip',ierr)
       end if

       ! ndmh: allocate memory for augmentation charges
       allocate(p_species(ps)%aug_q(p_species(ps)%n_shells, &
            p_species(ps)%n_shells), stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%aug_q',ierr)

       ! ndmh: allocate memory for L-independent function
       allocate(p_species(ps)%qfunc(p_species(ps)%kkbeta, &
            p_species(ps)%n_shells,p_species(ps)%n_shells), stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%qfunc',ierr)
       n_shells = p_species(ps)%n_shells
       allocate(p_species(ps)%qfcoef(p_species(ps)%nqfcoef, &
            0:p_species(ps)%qf_lmax,n_shells,n_shells),stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%qfcoef',ierr)
       allocate(p_species(ps)%rinner(0:p_species(ps)%qf_lmax),stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%rinner',ierr)
       allocate(p_species(ps)%rlog(p_species(ps)%mesh), stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%rlog',ierr)
       allocate(p_species(ps)%rab(p_species(ps)%mesh), stat=ierr)
       call utils_alloc_check('pseudo_species_alloc','p_species(ps)%rab',ierr)
       p_species(ps)%qfcoef(:,:,:,:) = 0.0_DP
       p_species(ps)%qfunc(:,:,:) = 0.0_DP
       p_species(ps)%rlog(:) = 0.0_DP
       p_species(ps)%rab(:) = 0.0_DP

    end do

    return

  contains

    subroutine internal_scan_recpot

      !============================================================!
      ! This subroutine reads the pseudopotential for a .recpot    !
      ! and stores the sizes required to allocate the necessary    !
      ! memory for this element in the p_species array.            !
      !------------------------------------------------------------!
      ! Assembled by Nicholas Hine on 29/01/09 based on code by    !
      ! Chris-Kriton Skylaris from 23/1/2004 and tidied up by      !
      ! Peter Haynes on 1/7/2004.                                  !
      !============================================================!

      use constants, only: ANGSTROM, HARTREE_IN_EVS
      use services, only: services_open_pspot_file
      use utils, only: utils_abort

      implicit none

      call services_open_pspot_file(2,current_file,ierr)

      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      psp_xc_functional = ''
      do while(index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line
         if (index(line,'LDA')>0) psp_xc_functional = 'LDA'
         if (index(line,'lda')>0) psp_xc_functional = 'LDA'
         if (index(line,'GGA')>0) psp_xc_functional = 'GGA'
         if (index(line,'gga')>0) psp_xc_functional = 'GGA'
         if (index(line,'PBE')>0) psp_xc_functional = 'GGA'
         if (index(line,'PW91')>0) psp_xc_functional = 'GGA'
      end do

      ! * Find out the number of projectors, and g-grid points
      !    --- stop at 1000 (or a flag of 4 characters) if the
      !        number of projectors is less than 32 (2 * s,p,d,f)

      count_points       = .true.
      tot_num_points     = 0
      tot_num_projectors = 0
      num_projectors     = 0

      current_job = 'counting projectors and points'
      line = repeat(' ',80)
      do while (len_trim(adjustl(line)) /= 4 .and. tot_num_projectors < 32)
         read(2,'(a80)',err=100,end=200) line

         if (len_trim(adjustl(line)) == 1) then
            read(line,*) i
            tot_num_projectors = tot_num_projectors + (2*i+1) ! 2l+1
            num_projectors = num_projectors + 1
            count_points = .false.
         end if

         if (count_points) then
            do i=1,len(line)
               if (line(i:i) == '.') tot_num_points = tot_num_points + 1
            end do
         end if
      end do

      ! cks: initialise
      p_species(ps)%n_shells  = num_projectors
      p_species(ps)%n_proj = tot_num_projectors
      p_species(ps)%n_rad_pts = tot_num_points - 1

      ! ndmh: skip to the end of the projectors
      current_job = 'seeking marker "1000"'
      do while ((index(line,'1000') == 0).and.(len_trim(adjustl(line)) /= 4))
         read(2,'(a80)',err=100,end=200) line
      end do

      ! ndmh: check if there are numbers here - if so there is non-linear
      ! ndmh: core correction data to read in
      read(2,'(a80)',iostat=ierr) line
      if (ierr /= 0) then
         p_species(ps)%core_charge = .false.
      else
         p_species(ps)%core_charge = .true.
      end if

      ! ** Rewind and find out the ionic charge
      rewind(2,iostat=ierr,err=400)
400   if (ierr /= 0) call utils_abort('Error in pseudo_species_alloc: rewinding&
           & file "'//trim(current_file)//'" failed with code ',ierr)

      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      do while (index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line
      end do

      current_job = 'skipping version string'
      read(2,*,err=100,end=200) line
      current_job = 'reading Gmax'
      read(2,*,err=100,end=200) gmax
      current_job = 'reading first two points of local potential'
      read(2,*,err=100,end=200) temp,temp

      ! * Calculate the ionic charge by looking at the g -> 0 part
      !   of the local potential
      fact1 = 4.0_DP*PI*HARTREE_IN_EVS/ANGSTROM

      ! pa: changed from
      ! ionic_charge = real(nint(-temp*(gmax/real(tot_num_points-1,DP))**2/fact1),kind=DP)
      ! pa: to allow fractional ionic charges
      ionic_charge = -temp*(gmax/real(tot_num_points-1,DP))**2/fact1

      !write(stdout,'(a,f9.6)') &
      !     'ionic_charge from psp before rounding: ', ionic_charge

      ! pa: round ionic charge to the nearest 1/cfac
      ionic_charge = (real(nint(cfac*ionic_charge), kind=DP))/real(cfac,DP)

      !write(stdout,'(a,f9.6)') &
      !            'ionic_charge from psp after rounding is: ', ionic_charge

      ! cks: convert gmax to atomic units
      p_species(ps)%ps_gmax = gmax/ANGSTROM

      ! ndmh: calculate spacing of recip space points
      p_species(ps)%inv_g_spacing = (p_species(ps)%n_rad_pts-1) / &
           p_species(ps)%ps_gmax

      p_species(ps)%ion_charge = ionic_charge

      ! ndmh: recpots do include the -Z/r, so flag to subtract it
      p_species(ps)%subtract_coul = .true.

      ! ndmh: recpots do not use kkbeta, mesh, qf_lmax or nqfcoef
      p_species(ps)%kkbeta = 0
      p_species(ps)%mesh = 0
      p_species(ps)%qf_lmax = 0
      p_species(ps)%nqfcoef = 0

      ! ** Close the file
      close(2,iostat=ierr,err=500)
500   if (ierr /= 0) call utils_abort('Error in pseudo_species_alloc: closing &
           & file "'//trim(current_file)//'" failed with code ',ierr)

      return

100   call utils_abort('Error in pseudo_species_alloc: reading file "'//&
           trim(current_file)//'" failed while '//trim(current_job))

200   call utils_abort('Error in pseudo_species_alloc: file "'//&
           trim(current_file)//'" ended unexpectedly while '//trim(current_job))

    end subroutine internal_scan_recpot

    subroutine internal_scan_usp

      !============================================================!
      ! This subroutine reads the pseudopotential for a .usp       !
      ! and stores the sizes required to allocate the necessary    !
      ! memory for this element in the p_species array.            !
      !------------------------------------------------------------!
      ! Written by Nicholas Hine on 29/01/09.                      !
      !============================================================!

      use constants, only: ANGSTROM
      use services, only: services_open_pspot_file
      use utils, only: utils_abort

      ! Local Variables
      integer, parameter :: max_proj=huge(1)
      integer :: version(2), usp_nlcc
      integer :: shell
      integer :: j,k
      integer :: lshell, lmax
      real(kind=DP) :: tmp1,tmp2

      call services_open_pspot_file(2,current_file,ierr)

      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      psp_xc_functional = ''
      do while(index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line
         if (index(line,'LDA')>0) psp_xc_functional = 'LDA'
         if (index(line,'lda')>0) psp_xc_functional = 'LDA'
         if (index(line,'GGA')>0) psp_xc_functional = 'GGA'
         if (index(line,'gga')>0) psp_xc_functional = 'GGA'
         if (index(line,'PBE')>0) psp_xc_functional = 'GGA'
         if (index(line,'PW91')>0) psp_xc_functional = 'GGA'
      end do

      ! ** Read the version
      read(2,*,err=100,end=200) version(1),version(2)

      ! ** Read the Ionic Charge
      read(2,*,err=100,end=200) ionic_charge

      p_species(ps)%ion_charge = ionic_charge

      ! ** Read the gmax, and nlcc flag
      read(2,*,err=100,end=200) gmax, usp_nlcc
      p_species(ps)%ps_gmax = gmax/ANGSTROM

      ! ** Determine whether to include nonlinear core corrections
      p_species(ps)%core_charge = (usp_nlcc==2)

      ! ** Find out the number of projectors, and g-grid points
      count_points = .true.
      tot_num_points = 0
      tot_num_projectors = 0
      num_projectors = 0
      line = ''
      lmax = -100
      do while((len_trim(adjustl(line))/=4).and.(num_projectors<max_proj))
         read(2,'(a)',err=100,end=200) line
         if(len_trim(adjustl(line))==1) then
            read(line,*,err=100,end=200) lshell
            if (lshell>lmax) lmax = lshell
            tot_num_projectors = tot_num_projectors + (2*lshell+1)
            num_projectors = num_projectors+1
            count_points = .false.
         end if
         if(count_points) then
            do i=1,len(line)
               if(line(i:i)=='.') tot_num_points = tot_num_points + 1
            end do
         end if
      end do

      p_species(ps)%n_rad_pts = tot_num_points
      p_species(ps)%n_proj = tot_num_projectors
      p_species(ps)%n_shells  = num_projectors

      ! Skip over the augmentation charges
      write(current_job,'(a80)') 'skipping augmentation charges'
      do shell=1,p_species(ps)%n_shells
         read(2,*,err=100,end=200) line
      end do

       ! ndmh: read in grid sizes for kkbeta and mesh
      write(current_job,'(a80)') 'reading grid sizes'
      if(p_species(ps)%core_charge) then
         read(2,*,err=100,end=200) p_species(ps)%kkbeta,p_species(ps)%mesh
      else
         read(2,*,err=100,end=200) p_species(ps)%kkbeta
         p_species(ps)%mesh = p_species(ps)%kkbeta
      end if

      ! Skip over the augmentation charges
      write(current_job,'(a80)') 'skipping grid points'
      read(2,*,err=100,end=200) (tmp1,i=1,p_species(ps)%mesh)
      read(2,*,err=100,end=200) (tmp2,i=1,p_species(ps)%mesh)

      ! ndmh: skip the L independent function
      write(current_job,'(a80)') 'skipping L independent function'
      do i=1,p_species(ps)%n_shells
         do j=1,i
            read(2,*,err=100,end=200) (tmp1,k=1,p_species(ps)%kkbeta)
         end do
      end do

      ! ndmh: read the number of qfunc coefficients
      read(2,*,err=100,end=200) p_species(ps)%nqfcoef

      ! ndmh: calculate maximum L value of qfcoefs
      p_species(ps)%qf_lmax = 2*lmax

      ! ndmh: calculate inverse of spacing of recip space points
      p_species(ps)%inv_g_spacing = (p_species(ps)%n_rad_pts-1) / &
           p_species(ps)%ps_gmax

      ! ndmh: usps do not include the -Z/r, so flag to not subtract it
      p_species(ps)%subtract_coul = .false.

      ! ** Close the file
      close(2,iostat=ierr,err=500)
500   if (ierr /= 0) call utils_abort('Error in pseudo_species_alloc: &
           &closing file "'//trim(current_file)//'" failed with code ',ierr)

      return

100   call utils_abort('Error in pseudo_species_alloc: reading file "'//&
           trim(current_file)//'" failed while '//trim(current_job))

200   call utils_abort('Error in pseudo_species_alloc: file "'//&
           trim(current_file)//'" ended unexpectedly while '//trim(current_job))

    end subroutine internal_scan_usp

  end subroutine pseudo_species_alloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_species_read_radial(p_species,cell)

    !=============================================================!
    ! This subroutine reads the radial parts (local and nonlocal) !
    ! of the pseudopotential for a particular ionic species.      !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/1/2004.              !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                   !
    ! Modified to read multiple file formats (recpots and usps)   !
    ! by Nicholas Hine on 29/01/2009.                             !
    ! par removed as input by Joseph Prentice, June 2019          !
    !=============================================================!

    use comms, only: comms_bcast, pub_on_root, pub_root_proc_id, &
                     pub_imroots_comm, comms_barrier
    use image_comms, only: pub_my_image
    use constants, only: ANGSTROM, HARTREE_IN_EVS
    use simulation_cell, only: CELL_INFO
    use rundat, only: pub_aug, pub_usp, pub_num_images
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_banner

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in)         :: cell
    type(PSEUDO_SPECIES), intent(inout) :: p_species(:)

    ! Local Variables
    integer :: i, img
    integer :: ps
    integer :: shell
    integer :: ierr
    real(kind=DP) :: delta
    character(len=64)    :: current_file
    character(len=80)    :: line
    character(len=3)     :: shell_string
    character(len=80)    :: current_job
    integer :: num_pspecies

    ! jcap: Get num_pspecies from p_species itself, not par
    num_pspecies = size(p_species)

    if (pub_on_root) then

       write(stdout,'(a)') utils_banner('<', 'Pseudopotential information')

       do ps=1,num_pspecies

          ! pdh: get file name for this species
          current_file = p_species(ps)%pseudo_name

          ! ndmh: determine what file format we are reading and scan it
          ! ndmh: to find sizes of various arrays

          do img=0, pub_num_images-1
             if ((pub_num_images == 1).or.(img == pub_my_image)) then
                ! ndmh: Old-style CASTEP .recpot
                if (index(current_file,'recpot') > 0) then

                   call internal_read_recpot

                ! ndmh: New-style CASTEP OTFG .usp (not necessarily ultrasoft!)
                else if (index(current_file,'usp') > 0) then

                   call internal_read_usp

                else
                   write(stdout,"(a)")"Error in pseudopotentials_read_species: &
                    & Unrecognised file format: ", current_file
                end if
             end if
             if (pub_num_images > 1) call comms_barrier(pub_imroots_comm)
          end do

          ! pdh: print basic information about this pseudopotential
          write(stdout,'(3a,i5,a,f6.1,a)') 'File: ',trim(current_file),' [', &
               p_species(ps)%n_rad_pts,' points up to Gmax=', &
               p_species(ps)%ps_gmax,' (1/bohr)]'
          ! pa: changed from i3 to f10.6
          write(stdout,'(a,i3,a,f10.6)') '  Atomic number:', &
               p_species(ps)%atomic_number,';  ionic charge:', &
               p_species(ps)%ion_charge
          do shell=1,p_species(ps)%n_shells
             write(stdout,'(2(a,i2),a,f5.2,a)') '    Shell',shell,': l =', &
                  p_species(ps)%ang_mom(shell),'; rc =', &
                  p_species(ps)%core_radius(shell),' bohr'
          end do
          if (p_species(ps)%n_shells < 1) write(stdout,'(a)') '    Local potential'
          if (p_species(ps)%core_charge) write(stdout,'(a)') &
               '  Core charge supplied for Nonlinear Core Corrections'

       end do

       write(stdout,'(a/)') utils_banner('<')

       ! Subtract Coloumb potential from any recpots
       call subtract_or_add_coul(p_species)

    end if


    ! cks: broadcast results from root to the rest of the processors
    do ps=1,num_pspecies
       call comms_bcast(pub_root_proc_id,p_species(ps)%core_radius)
       call comms_bcast(pub_root_proc_id,p_species(ps)%ang_mom)
       call comms_bcast(pub_root_proc_id,p_species(ps)%D0)
       call comms_bcast(pub_root_proc_id,p_species(ps)%rad_locpot_recip)
       call comms_bcast(pub_root_proc_id,p_species(ps)%rad_proj_recip)
       if (p_species(ps)%core_charge) then
          call comms_bcast(pub_root_proc_id,p_species(ps)%core_charge_recip)
       end if
       ! ndmh: Broadcast augmentation charge info
       call comms_bcast(pub_root_proc_id,p_species(ps)%aug_q)
       call comms_bcast(pub_root_proc_id,p_species(ps)%qfunc)
       call comms_bcast(pub_root_proc_id,p_species(ps)%rinner)
       if (p_species(ps)%nqfcoef > 0) then
          do i=1,p_species(ps)%n_shells
             call comms_bcast(pub_root_proc_id,p_species(ps)%qfcoef(:,:,:,i))
          end do
       end if
       call comms_bcast(pub_root_proc_id,p_species(ps)%rlog)
       call comms_bcast(pub_root_proc_id,p_species(ps)%rab)

    end do

    return

  contains

    subroutine internal_read_recpot

      !==============================================================!
      ! This subroutine reads the pseudopotential for a .recpot      !
      ! and stores the local and nonlocal components of the          !
      ! potential in the arrays allocated in pseudo_alloc_species !
      !--------------------------------------------------------------!
      ! Assembled by Nicholas Hine on 29/01/09 based on code by      !
      ! Chris-Kriton Skylaris from 23/1/2004 and tidied up by        !
      ! Peter Haynes on 1/7/2004.                                    !
      !==============================================================!

      use rundat, only: pub_smooth_projectors
      use services, only: services_open_pspot_file
      use utils, only: utils_abort, utils_open_unit_check, &
           utils_close_unit_check

      ! Local Variables
      ! cks: for scaling magic
      real(kind=DP) :: g, factor, denom

      call services_open_pspot_file(2,current_file,ierr)

      ! * Get to the top of the pseudopotential proper
      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      do while (index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line  ! To the end of the comments
      end do
      current_job = 'moving to end of header'
      do i=1,2
         read(2,'(a80)',err=100,end=200) line  ! To the end of the header
      end do

      ! * Read the local component of the pseudopotential
      current_job = 'reading local component'
      read(2,*,err=100,end=200) (p_species(ps)%rad_locpot_recip(i), &
           i=1,p_species(ps)%n_rad_pts)


      ! cks: scale radial locpot with suitable factors to make it
      !      consistent with atomic units
      ! cks, 28/3/2004: numerical treatment consistent with CASTEP
      factor = ((ANGSTROM**3)/HARTREE_IN_EVS) / (4.0_DP*PI)
      p_species(ps)%rad_locpot_recip(:) = &
           p_species(ps)%rad_locpot_recip(:) * factor

      ! cks, 28/3/2004: scaling near low g consistent with CASTEP
      g = p_species(ps)%ps_gmax / real(p_species(ps)%n_rad_pts-1,kind=DP)

      ! pa: removed real(...,DP) for fractional ion charges
      factor = 3.75_DP*p_species(ps)%ion_charge / (g*g)

      denom = p_species(ps)%rad_locpot_recip(3) - &
           4.0_DP * p_species(ps)%rad_locpot_recip(2) + &
           3.0_DP * p_species(ps)%rad_locpot_recip(1)

      if (abs(denom) > epsilon(1.0_DP)) factor = factor / denom

      p_species(ps)%rad_locpot_recip(:) = &
           p_species(ps)%rad_locpot_recip(:) * factor

      ! * Read the non-local component of the pseudopotential
      p_species(ps)%D0(:,:) = 0.0_DP
      do shell=1,p_species(ps)%n_shells
         write(shell_string,'(i3)') shell

         write(current_job,'(a80)') 'reading angular momentum of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) p_species(ps)%ang_mom(shell)

         write(current_job,'(a80)') 'reading KB denominator of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) p_species(ps)%D0(shell,shell)

         ! cks: convert to atomic units
         ! cks: scaling consistent with CASTEP
         p_species(ps)%D0(shell,shell) = p_species(ps)%D0(shell,shell) * &
              HARTREE_IN_EVS
         p_species(ps)%D0(shell,shell) = 1.0_DP / p_species(ps)%D0(shell,shell)

         write(current_job,'(a80)') 'reading projector of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) (p_species(ps)%rad_proj_recip(i,shell), &
              i=1,p_species(ps)%n_rad_pts)

         ! cks: convert projector units from A^{3/2}*eV to a0^{3/2}*Eh
         ! cks: scaling consistent with CASTEP
         p_species(ps)%rad_proj_recip(:,shell) = &
              p_species(ps)%rad_proj_recip(:,shell) * sqrt(ANGSTROM**3)

      end do

      ! cks: initialise the radial grid for the interpolation
      ! pdh: inline this call to avoid awkward copy of arguments
      delta = p_species(ps)%ps_gmax / real(p_species(ps)%n_rad_pts-1,DP)
      do i=1,p_species(ps)%n_rad_pts
         p_species(ps)%rad_proj_recip(i,p_species(ps)%n_shells+1) = &
              real(i-1,DP) * delta
      end do

      ! pdh: get core radii from projectors
      call internal_get_core_radii(p_species(ps))

      ! cks: Gauss filter projectors in reciprocal and real space
      if (pub_smooth_projectors >= 0.0_DP ) then
         call pseudo_gauss_filter_proj(p_species(ps),cell)
      endif

      ! ndmh: find end of projectors
      current_job = 'seeking marker "1000"'
      do while ((index(line,'1000') == 0).and.(len_trim(adjustl(line)) /= 4))
         read(2,'(a80)',err=100,end=200) line
      end do

      ! ndmh: get NLCC core charge if present
      if (p_species(ps)%core_charge) then
         read(2,*,err=100,end=200) (p_species(ps)%core_charge_recip(i), &
              i=1,p_species(ps)%n_rad_pts)
      end if

      ! ndmh: recpots do not use kkbeta or mesh
      p_species(ps)%rlog = 0.0_DP
      p_species(ps)%rab = 0.0_DP
      p_species(ps)%qfunc = 0.0_DP
      p_species(ps)%aug_q = 0.0_DP

      close(2,iostat=ierr,err=400)
400   call utils_close_unit_check('internal_read_recpot &
           &(pseudo_species_read_radial)',trim(current_file),ierr)

      return

100   call utils_abort('Error in internal_read_recpot (pseudo_species_read_&
           &radial): reading file "'//trim(current_file)//'" failed while '//&
           trim(current_job))

200   call utils_abort('Error in internal_read_recpot (pseudo_species_read_&
           &radial): file "'//trim(current_file)//&
           '" ended unexpectedly while '//trim(current_job))

    end subroutine internal_read_recpot

    subroutine internal_read_usp

      !==============================================================!
      ! This subroutine reads the pseudopotential for a .usp         !
      ! and stores the local and nonlocal components of the          !
      ! potential in the arrays allocated in pseudo_alloc_species !
      !--------------------------------------------------------------!
      ! Written by Nicholas Hine on 29/01/09                         !
      !==============================================================!

      use rundat, only: pub_smooth_projectors
      use services, only: services_radial_transform, services_open_pspot_file
      use utils, only: utils_abort, utils_open_unit_check, &
           utils_close_unit_check

      ! Local Variables
      integer,parameter :: max_proj=huge(1)
      integer :: j, k, l
      integer :: lshell
      real(kind=DP) :: factor
      real(kind=DP), allocatable :: core_charge_real(:)

      call services_open_pspot_file(2,current_file,ierr)

      ! ndmh: get to the top of the pseudopotential proper
      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      do while (index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line  ! To the end of the comments
      end do
      current_job = 'moving to end of header'
      do i=1,3
         read(2,'(a80)',err=100,end=200) line  ! To the end of the header
      end do

      ! ndmh: read the local component of the pseudopotential
      read(2,*,err=100,end=200) (p_species(ps)%rad_locpot_recip(i), &
           i=1,p_species(ps)%n_rad_pts)

      ! ndmh: Scale radial locpot with suitable factors to convert it
      ! ndmh: to atomic units
      factor = ((ANGSTROM**3)/HARTREE_IN_EVS) / (4.0_DP*PI)
      p_species(ps)%rad_locpot_recip(:) = &
           p_species(ps)%rad_locpot_recip(:) * factor

      ! ndmh: Read the non-local component of the pseudopotential
      p_species(ps)%D0(:,:) = 0.0_DP
      do shell=1,p_species(ps)%n_shells

         write(shell_string,'(i3)') shell

         write(current_job,'(a80)') 'reading angular momentum of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) p_species(ps)%ang_mom(shell)

         ! Determine lower limit of shells of this l
         l = p_species(ps)%ang_mom(shell)
         do lshell=1,p_species(ps)%n_shells
            if (p_species(ps)%ang_mom(lshell)==l) exit
         end do

         ! Read KB energies
         write(current_job,'(a80)') 'reading KB energy of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) (p_species(ps)%D0(i,shell), &
              i=lshell,shell)

         ! Fill in other half of D0 matrix
         do i=lshell,shell
            p_species(ps)%D0(shell,i) = p_species(ps)%D0(i,shell)
         end do

         write(current_job,'(a80)') 'reading projector of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) (p_species(ps)%rad_proj_recip(i,shell), &
              i=1,p_species(ps)%n_rad_pts)

         ! cks: convert projector units from A^{3/2}*eV to a0^{3/2}*Eh
         ! cks: scaling consistent with CASTEP
         p_species(ps)%rad_proj_recip(:,shell) = &
              p_species(ps)%rad_proj_recip(:,shell) * sqrt(ANGSTROM**3)

      end do

      ! Convert units of D0
      p_species(ps)%D0(:,:) = p_species(ps)%D0(:,:) / HARTREE_IN_EVS

      ! cks: initialise the radial grid for the interpolation
      ! pdh: inline this call to avoid awkward copy of arguments
      delta = p_species(ps)%ps_gmax / real(p_species(ps)%n_rad_pts-1,DP)
      do i=1,p_species(ps)%n_rad_pts
         p_species(ps)%rad_proj_recip(i,p_species(ps)%n_shells+1) = &
              real(i-1,DP) * delta
      end do

      ! ndmh: check we are at the right point of the file
      !    -- There is no '1000' flag if there are 2 projectors per
      !       all channels up to l=3 (the assumed maximum)
      write(current_job,'(a80)') 'reading "1000" '
      if(p_species(ps)%n_shells < max_proj) then
         read(2,'(a)',err=100,end=200) line
         if(trim(adjustl(line))/='1000') then
            call utils_abort('Error in internal_read_usp '//&
                '(pseudo_species_read_radial): file "'//trim(current_file)//&
                '" ended unexpectedly while '//trim(current_job))
         end if
      end if

      ! ndmh: Read in augmentation charges
      write(current_job,'(a80)') 'reading augmentation charges'
      p_species(ps)%aug_q(:,:) = 0.0_DP
      do shell=1,p_species(ps)%n_shells
         l = p_species(ps)%ang_mom(shell)
         ! Determine lower limit of shells of this l
         do lshell=1,p_species(ps)%n_shells
            if (p_species(ps)%ang_mom(lshell)==l) exit
         end do
         ! Read in the q-values
         read(2,*,err=100,end=200) (p_species(ps)%aug_q(i,shell), &
              i=lshell,shell)
         ! Fill in other half of matrix
         do i=lshell,shell
            p_species(ps)%aug_q(shell,i) = p_species(ps)%aug_q(i,shell)
         end do
      end do
      if (any(abs(p_species(ps)%aug_q(:,:))>tiny(1.0_DP))) then
         !call utils_abort('Error in internal_read_usp (pseudo_species_read_&
         !     &radial): Nonzero augmentation charge found. Ultrasoft &
         !     &Pseudopotentials with augmentation charges are not supported.')
         pub_usp = .true.
         pub_aug = .true.
      end if

      ! ndmh: read in grid sizes for kkbeta and mesh
      write(current_job,'(a80)') 'reading grid sizes'
      if(p_species(ps)%core_charge) then
         read(2,*,err=100,end=200) p_species(ps)%kkbeta,p_species(ps)%mesh
      else
         read(2,*,err=100,end=200) p_species(ps)%kkbeta
         p_species(ps)%mesh = p_species(ps)%kkbeta
      end if

      ! * Allocate temporary storage for core_charge_real if required
      if (p_species(ps)%core_charge) then
         allocate(core_charge_real(p_species(ps)%mesh),stat=ierr)
         call utils_alloc_check('internal_read_usp '//&
              '(pseudo_species_read_radial)','core_charge_real',ierr)
      end if

      ! ndmh: read in the grid points
      write(current_job,'(a80)') 'reading grid points'
      read(2,*,err=100,end=200) (p_species(ps)%rlog(i),i=1,p_species(ps)%mesh)
      read(2,*,err=100,end=200) (p_species(ps)%rab(i),i=1,p_species(ps)%mesh)

      ! ndmh: read the L independent function
      ! ndmh: this is zero for ncpp's
      write(current_job,'(a80)') 'reading L independent function'
      do i=1,p_species(ps)%n_shells
         do j=1,i
            read(2,*,err=100,end=200) (p_species(ps)%qfunc(k,i,j),k=1, &
                 p_species(ps)%kkbeta)
         end do
      end do
      do i=1,p_species(ps)%n_shells
         do j=i+1,p_species(ps)%n_shells
            p_species(ps)%qfunc(:,i,j) = p_species(ps)%qfunc(:,j,i)
         end do
      end do

      ! ndmh: the L dependent coeffs
      ! ndmh: these are zero for ncpp's
      p_species(ps)%qfcoef(:,:,:,:) = 0.0_DP
      write(current_job,'(a80)') 'reading L dependent coeffs'
      p_species(ps)%qf_lmax = 2*maxval(p_species(ps)%ang_mom( &
           1:p_species(ps)%n_shells))
      read(2,*,err=100,end=200) p_species(ps)%nqfcoef
      do l=0,p_species(ps)%qf_lmax
         read(2,*,err=100,end=200) p_species(ps)%rinner(l)
         do i=1,p_species(ps)%n_shells
            do j=1,i
               read(2,*,err=100,end=200) (p_species(ps)%qfcoef(k,l,i,j), &
                    k=1,p_species(ps)%nqfcoef)
            end do
         end do
      end do
      do l=0,p_species(ps)%qf_lmax
         do i=1,p_species(ps)%n_shells
            do j=i+1,p_species(ps)%n_shells
               p_species(ps)%qfcoef(:,l,i,j) = p_species(ps)%qfcoef(:,l,j,i)
            end do
         end do
      end do


      ! ndmh: finally read in the core charge on the log grid
      if (p_species(ps)%core_charge) then
         read(2,*,err=100,end=200) &
              (core_charge_real(i),i=1,p_species(ps)%mesh)

         ! ndmh: transform the real-space core charge to reciprocal space
         call services_radial_transform(0,0,p_species(ps)%mesh,p_species(ps)%rlog, &
              p_species(ps)%rab,p_species(ps)%n_rad_pts,&
              p_species(ps)%ps_gmax,core_charge_real,&
              p_species(ps)%core_charge_recip)

      end if

      ! pdh: get core radii from projectors
      call internal_get_core_radii(p_species(ps))

      ! cks: Gauss filter projectors in reciprocal and real space
      if (pub_smooth_projectors >= 0.0_DP ) then
         call pseudo_gauss_filter_proj(p_species(ps),cell)
      endif

      ! ndmh: deallocate temporary storage
      if (p_species(ps)%core_charge) then
         deallocate(core_charge_real,stat=ierr)
         call utils_dealloc_check('internal_read_usp '//&
              '(pseudo_species_read_radial)','core_charge_real',ierr)
      end if

      ! ndmh: close the file
      close(2,iostat=ierr,err=400)
400   call utils_close_unit_check('internal_read_usp &
           &(pseudo_species_read_radial)',trim(current_file),ierr)

      return

100   call utils_abort('Error in internal_read_usp '//&
           '(pseudo_species_read_radial): reading file "'//trim(current_file)//&
           '" failed while '//trim(current_job))

200   call utils_abort('Error in internal_read_usp '//&
           '(pseudo_species_read_radial): file "'//trim(current_file)//&
           '" ended unexpectedly while '//trim(current_job))

    end subroutine internal_read_usp

    subroutine internal_get_core_radii(p_species)

      !==================================================================!
      ! This subroutine calculates the core radii of the pseudopotential !
      ! species by examining the projectors in real space.               !
      !------------------------------------------------------------------!
      ! Written by Peter D. Haynes in early 2004.                        !
      ! Tidied up by Peter D. Haynes on 1/7/2004.                        !
      !==================================================================!

      use services, only: services_sbessj
      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argument
      type(PSEUDO_SPECIES), intent(inout) :: p_species

      ! Local variables
      integer :: ierr                           ! Error flag
      integer :: i_shell                        ! Shell counter
      integer :: l                              ! Angular momentum of shell
      integer :: g_rad_pt,r_rad_pt              ! Grid point counters
      integer :: n_r_rad_pts                    ! Num points on real radial grid
      real(kind=DP), parameter :: rmax = 5.0_DP ! Maximum core radius
      real(kind=DP), parameter :: halfpi = 0.5_DP * PI
      real(kind=DP) :: gnorm,rnorm              ! Norms in recip/real space
      real(kind=DP) :: proj,proj_max            ! Projector (max value) in rspc
      real(kind=DP) :: val,fac                  ! Variables for integration
      real(kind=DP) :: delr,r,delg,g            ! Grid variables
      real(kind=DP), allocatable :: integrand(:)! gspc integrand for FT

      ! Set parameters for radial grid
      delg = p_species%ps_gmax / (p_species%n_rad_pts - 1)
      delr = pi / p_species%ps_gmax
      n_r_rad_pts = nint(rmax / delr) + 1
      fac = (7.0_DP/17280.0_DP) * delg / sqrt(halfpi)

      ! Allocate workspace for integration
      allocate(integrand(p_species%n_rad_pts),stat=ierr)
      call utils_alloc_check('internal_get_core_radii &
           &(pseudo_species_read_radial)','integrand',ierr)

      ! Loop over shells: angular momentum is l
      do i_shell=1,p_species%n_shells
         l = p_species%ang_mom(i_shell)

         ! Calculate norm from reciprocal space representation
         gnorm = 0.0_DP
         do g_rad_pt=1,p_species%n_rad_pts
            g = (g_rad_pt-1) * delg
            val = p_species%rad_proj_recip(g_rad_pt,i_shell) * g
            gnorm = gnorm + val*val
         end do
         gnorm = gnorm * delg

         if (gnorm < epsilon(1.0_DP)) then
            p_species%core_radius(i_shell) = 0.0_DP
            cycle
         end if

         ! Calculate projector on real-space radial grid, and cumulative
         ! norm of projector in real space
         rnorm = 0.0_DP
         proj_max = 0.0_DP

         do r_rad_pt=1,n_r_rad_pts
            r = (r_rad_pt-1) * delr
            p_species%core_radius(i_shell) = 1.8_DP
            integrand(1) = 0.0_DP
            do g_rad_pt=2,p_species%n_rad_pts
               g = (g_rad_pt-1) * delg
               integrand(g_rad_pt) = services_sbessj(l,g*r) * g*g * &
                    p_species%rad_proj_recip(g_rad_pt,i_shell)
            end do

            ! Due to the oscillatory nature of the integrand,
            ! an eight point Newton-Cotes integration method is used.
            ! See  Abramowitz and Stegun Eq. 25.4.17

            proj = 0.0_DP
            do g_rad_pt=1,p_species%n_rad_pts-7,7
               proj = proj + &
                    751.0_DP*(integrand(g_rad_pt)+integrand(g_rad_pt+7)) + &
                    3577.0_DP*(integrand(g_rad_pt+1)+integrand(g_rad_pt+6)) + &
                    1323.0_DP*(integrand(g_rad_pt+2)+integrand(g_rad_pt+5)) + &
                    2989.0_DP*(integrand(g_rad_pt+3)+integrand(g_rad_pt+4))
            end do
            proj = proj * fac
            val = r * proj
            rnorm = rnorm + val*val * delr

            ! Set core radius if magnitude of projector is less then 1% of
            ! maximum value and cumulative norm is at least 99.9% of
            ! reciprocal space norm

            proj_max = max(proj_max,abs(proj))
            if (r_rad_pt > 1) then
               if (abs(proj/proj_max) < 1.0e-2_DP .and. &
                    rnorm/gnorm > 0.999_DP) then
                  p_species%core_radius(i_shell) = r
                  exit
               end if
            end if

         end do

      end do

      ! Free workspace
      deallocate(integrand,stat=ierr)
      call utils_dealloc_check('internal_get_core_radii &
           &(pseudo_species_read_radial)','integrand',ierr)

    end subroutine internal_get_core_radii

  end subroutine pseudo_species_read_radial


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_gauss_filter_proj(current_species,cell)

    !=========================================================+=========!
    ! This subroutine applies a Gaussian smoothing filter to a          !
    ! non-local projector while it is still in radial form.             !
    !-------------------------------------------------------------------!
    !   *Description of the algorithm*                                  !
    ! First the smoothing is applied in reciprocal space.               !
    ! The part of the projector from 0.667g_grid to g_max is multiplied !
    ! by a Gaussian which begins with the value 1.0 at 0.667*g_grid     !
    ! and decays with 0.1 of an exponent whose half-width which is      !
    ! equal to pub_smooth_projectors*g_grid.                            !
    ! Then the projector is Fourier transformed to radial real space    !
    ! and again filtered with a Gaussian in the region r_core -> infty. !
    ! The halfwidth of the Gaussian is now equal to                     !
    ! pub_smooth_projectors*r_core.                                     !
    ! Finally the projector is put back to reciprocal space.            !
    !-------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 20/01/2005 using parts of     !
    ! Peter Haynes' subroutine "internal_get_core_radii".               !
    !===================================================================!

    use rundat, only: pub_smooth_projectors
    use services, only: services_sbessj
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Argument
    type(CELL_INFO), intent(in) :: cell
    type(PSEUDO_SPECIES), intent(inout) :: current_species

    ! Local variables
    integer :: ierr                           ! Error flag
    integer :: i_shell                        ! Shell counter
    integer :: l                              ! Angular momentum of shell
    integer :: g_rad_pt,r_rad_pt              ! Grid point counters
    integer :: n_r_rad_pts                    ! Num points on real radial grid
    integer :: max_npts         ! maximum number of real and recip gridpoints
    real(kind=DP), parameter :: rmax = 10.0_DP ! Maximum real space distance
    real(kind=DP), parameter :: halfpi = 0.5_DP * PI
    real(kind=DP) :: proj                     ! Projector in rspc
    real(kind=DP) :: fac                      ! Variables for integration
    real(kind=DP) :: delr, delg               ! Grid variables
    real(kind=DP), allocatable :: integrand(:)! gspc integrand for FT
    real(kind=DP), allocatable :: rspace_projector (:)
    real(kind=DP) :: rval      ! grid distance in real space
    real(kind=DP) :: gval      ! grid distance in reciprocal space
    real(kind=DP) :: g_grid    ! 1-D maximum g-vector of psinc-grid
    real(kind=DP) :: fac_real  ! factor to multiply real space integral
    real(kind=DP) :: alpha     ! exponent of Gaussian
    real(kind=DP) :: core_rad  ! core radius of projector

    ! Set parameters for radial grid
    delg = current_species%ps_gmax / (current_species%n_rad_pts - 1)
    delr = PI / current_species%ps_gmax

    ! cks: Overkill to retain all precission in real-> recip FT
    delr =0.1_DP*delr

    n_r_rad_pts = nint(rmax / delr) + 1

    max_npts =max(n_r_rad_pts, current_species%n_rad_pts)

    g_grid =3.0_DP*PI/(cell%d1 + cell%d2 +cell%d3)

    fac = (7.0_DP/17280.0_DP) * delg  / halfpi
    fac_real =(7.0_DP/17280.0_DP) * delr

    if (pub_smooth_projectors > tiny(1.0_DP) ) then
       alpha = 0.1_DP*log(2.0_DP)/( (pub_smooth_projectors*g_grid)**2 )
    else
       alpha = 0.0_DP
    end if

    ! Allocate workspace for integration
    allocate(integrand(max_npts),stat=ierr)
    call utils_alloc_check('pseudo_gauss_filter_proj','integrand',ierr)

    ! cks: Allocate rspace_projector
    allocate(rspace_projector(n_r_rad_pts),stat=ierr)
    call utils_alloc_check('pseudo_gauss_filter_proj', &
         'rspace_projector',ierr)

    ! Loop over shells: angular momentum is l
    do i_shell=1, current_species%n_shells
       l = current_species%ang_mom(i_shell)


       ! cks: ++++++++ Gauss filter projector in reciprocal space +++++
       do g_rad_pt =1, current_species%n_rad_pts

          gval =current_species%rad_proj_recip(g_rad_pt, current_species%n_shells+1)

          if ( gval .gt. (0.667_DP*g_grid ) ) then
             if (pub_smooth_projectors > tiny(1.0_DP) ) then
                current_species%rad_proj_recip(g_rad_pt, i_shell) = &
                     exp( -alpha*(( gval -0.667_DP*g_grid )**2) ) &
                     * current_species%rad_proj_recip(g_rad_pt, i_shell)
             else
                current_species%rad_proj_recip(g_rad_pt, i_shell) =0.0_DP
             endif
          endif

       enddo
       ! cks: +++++ END Gauss filter projector in reciprocal space +++++




       ! cks: ====== RECIP to REAL transform ===========================

       ! cks: Calculate projector on real-space radial grid
       do r_rad_pt =1, n_r_rad_pts

          rval = (r_rad_pt-1) * delr
          integrand(1) = 0.0_DP
          do g_rad_pt=2, current_species%n_rad_pts

             gval =current_species%rad_proj_recip(g_rad_pt, current_species%n_shells+1)

             integrand(g_rad_pt) = services_sbessj(l,gval*rval) * gval*gval * &
                  current_species%rad_proj_recip(g_rad_pt,i_shell)

          end do

          ! Due to the oscillatory nature of the integrand,
          ! an eight point Newton-Cotes integration method is used.
          ! See  Abramowitz and Stegun Eq. 25.4.17
          proj = 0.0_DP
          do g_rad_pt=1, current_species%n_rad_pts-7,7
             proj = proj + &
                  751.0_DP*(integrand(g_rad_pt)+integrand(g_rad_pt+7)) + &
                  3577.0_DP*(integrand(g_rad_pt+1)+integrand(g_rad_pt+6)) + &
                  1323.0_DP*(integrand(g_rad_pt+2)+integrand(g_rad_pt+5)) + &
                  2989.0_DP*(integrand(g_rad_pt+3)+integrand(g_rad_pt+4))
          end do
          proj = proj * fac

          rspace_projector(r_rad_pt) =proj

       end do

       ! cks: == END RECIP to REAL transform ===========================



       ! cks: ++++++++ Gauss filter projector in real space +++++
       core_rad =current_species%core_radius(i_shell)

       if (pub_smooth_projectors > tiny (1.0_DP) ) then
          alpha =log(2.0_DP)/( (pub_smooth_projectors*core_rad)**2 )
       endif


       do r_rad_pt =1, n_r_rad_pts

          rval = (r_rad_pt-1) * delr

          if ( rval .gt. core_rad ) then
             if (pub_smooth_projectors > tiny (1.0_DP) ) then
                rspace_projector(r_rad_pt) = &
                     exp( -alpha*(( rval -core_rad )**2) ) &
                     * rspace_projector(r_rad_pt)
             else
                rspace_projector(r_rad_pt) =0.0_DP
             endif
          endif

       enddo
       ! cks: ++++ END Gauss filter projector in real space +++++



       ! cks: ====== REAL to RECIP transform ===========================

       ! cks: put projector back to reciprocal space
       do g_rad_pt =1, current_species%n_rad_pts

          gval =current_species%rad_proj_recip(g_rad_pt, current_species%n_shells+1)
          integrand(1) = 0.0_DP
          do r_rad_pt= 2, n_r_rad_pts
             rval = (r_rad_pt-1) * delr
             integrand(r_rad_pt) = services_sbessj(l,gval*rval) *rval *rval * &
                  rspace_projector(r_rad_pt)
          end do


          ! Due to the oscillatory nature of the integrand,
          ! an eight point Newton-Cotes integration method is used.
          ! See  Abramowitz and Stegun Eq. 25.4.17
          proj = 0.0_DP
          do r_rad_pt= 1, n_r_rad_pts -7,7
             proj = proj + &
                  751.0_DP*(integrand(r_rad_pt)+integrand(r_rad_pt+7)) + &
                  3577.0_DP*(integrand(r_rad_pt+1)+integrand(r_rad_pt+6)) + &
                  1323.0_DP*(integrand(r_rad_pt+2)+integrand(r_rad_pt+5)) + &
                  2989.0_DP*(integrand(r_rad_pt+3)+integrand(r_rad_pt+4))
          end do

          current_species%rad_proj_recip(g_rad_pt,i_shell) =proj *fac_real
       end do
       ! cks: == END REAL to RECIP transform ===========================

    end do

    ! Free workspace
    deallocate(rspace_projector, stat=ierr)
    call utils_dealloc_check('pseudo_gauss_filter_proj', &
         'rspace_projector',ierr)

    ! Free workspace
    deallocate(integrand,stat=ierr)
    call utils_dealloc_check('pseudo_gauss_filter_proj','integrand',ierr)

  end subroutine pseudo_gauss_filter_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_species_init_proj(proj_set,p_species,elements,par,is_cmplx)

    !=================================================================!
    ! This subroutine allocates and initialises all fftbox_proj_recip !
    ! for each p_species element. Each such fftbox_proj_recip is      !
    ! a full non-local projector in reciprocal space in the FFTbox.   !
    ! Optionally, this subroutine can also calculate instead Del_G of !
    ! the projectors along a given cartesian direction given by cart  !
    !-----------------------------------------------------------------!
    ! Original version by Chris-Kriton Skylaris, January 2004.        !
    ! Tidied up by Peter D. Haynes, July 2004.                        !
    ! Re-write by Nicholas Hine, creating PROJECTOR_SET type in       !
    ! July 2009.                                                      !
    ! Del_G code added by Laura Ratcliff in March 2011 as a separate  !
    ! routine (pseudo_species_calc_grad_proj)                         !
    ! Merged back into this routine to avoid code duplication by      !
    ! Nicholas Hine in April 2011.                                    !
    !=================================================================!

    use comms, only: pub_on_root
    use constants, only: VERBOSE
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
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(PSEUDO_SPECIES), intent(inout) :: p_species(par%num_pspecies)
    ! agrecocmplx
    logical, optional, intent(in) :: is_cmplx

    ! Local Variables
    integer :: iat, orig_iat
    integer :: ps
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
         write(stdout,'(a)') &
         '... Nonlocal pseudopotential projector initialisation'

    ! ndmh: find number of unique projectors
    proj_count = 0
    do ps=1,par%num_pspecies
       do shell=1,p_species(ps)%n_shells
          proj_count = proj_count + 2*p_species(ps)%ang_mom(shell) + 1
       end do
    end do

    ! Set up the entries of proj_set
    proj_set%n_proj_species = par%num_pspecies

    if (.not. allocated(proj_set%fftbox_proj_recip)) then
       ! agrecocmplx
       call projectors_allocate_set(proj_set, &
            maxval(p_species(:)%n_shells),maxval(p_species(:)%n_rad_pts), &
            is_cmplx=loc_cmplx, par=par)
    else
       call utils_abort('Error in pseudo_species_init_proj: &
            &proj_set already allocated')
    end if

    ! ndmh: set species_num_proj and species_first_proj values
    proj_count = 1
    do ps=1,proj_set%n_proj_species
       proj_set%species_num_proj(ps) = 0
       proj_set%species_first_proj(ps) = proj_count
       proj_set%gmax(ps) = p_species(ps)%ps_gmax
       proj_set%n_rad_pts(ps) = p_species(ps)%n_rad_pts
       proj_set%num_shells(ps) = p_species(ps)%n_shells
       proj_set%ang_mom(:,ps) = 0
       proj_set%rad_proj_recip(:,:,ps) = 0.0_DP
       do shell=1,p_species(ps)%n_shells
          proj_set%ang_mom(shell,ps) = p_species(ps)%ang_mom(shell)
          proj_set%species_num_proj(ps) = &
               proj_set%species_num_proj(ps) + &
               2*p_species(ps)%ang_mom(shell) + 1
          proj_set%rad_proj_recip(1:p_species(ps)%n_rad_pts,shell,ps) = &
               p_species(ps)%rad_proj_recip(1:p_species(ps)%n_rad_pts,shell)
       end do
       proj_count = proj_count + proj_set%species_num_proj(ps)
    end do

    ! ndmh: copy projector centre and radius from elements array
    do iat=1,par%nat
       orig_iat = par%orig_atom(iat)
       proj_set%proj_centre(iat) = elements(orig_iat)%centre
       proj_set%proj_max_radius(iat) = elements(orig_iat)%max_core_radius
       proj_set%proj_species(iat) = elements(orig_iat)%pspecies_number
    end do

    ! ndmh: Initialise projectors in reciprocal-space fftbox representation
    !call projectors_init_fftbox_recip(proj_set,kpt,delta,cart,swap_rc)

  end subroutine pseudo_species_init_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_species_exit(p_species)

    !==========================================================!
    ! This subroutine deallocates all the memory allocated for !
    ! the p_species array.                                     !
    !----------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/1/2004.           !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                !
    ! Modified to also deallocate nlps_projectors by Nicholas  !
    ! Hine on 02/11/2009.                                      !
    !==========================================================!

    use gaunt_coeff, only: gaunt_exit
    use rundat, only: pub_usp
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PSEUDO_SPECIES), intent(inout), pointer :: p_species(:)

    ! Local Variables
    integer :: ierr
    integer :: ps
    !character(len=32) :: string

    if (pub_usp) then
       call gaunt_exit
    end if

    if (associated(p_species)) then

       do ps=size(p_species),1,-1

          if (allocated(p_species(ps)%aug_q)) then
             deallocate(p_species(ps)%aug_q,stat=ierr)
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  'p_species(ps)%aug_q',ierr)
          end if

          if (allocated(p_species(ps)%qfunc)) then
             deallocate(p_species(ps)%qfunc,stat=ierr)
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  'p_species(ps)%qfunc',ierr)
          end if

          if (allocated(p_species(ps)%qfcoef)) then
             deallocate(p_species(ps)%qfcoef,stat=ierr)
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  'p_species(ps)%qfcoef',ierr)
          end if

          if (allocated(p_species(ps)%rinner)) then
             deallocate(p_species(ps)%rinner,stat=ierr)
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  'p_species(ps)%rinner',ierr)
          end if

          if (allocated(p_species(ps)%rlog)) then
             deallocate(p_species(ps)%rlog,stat=ierr)
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  'p_species(ps)%rlog',ierr)
          end if

          if (allocated(p_species(ps)%rab)) then
             deallocate(p_species(ps)%rab,stat=ierr)
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  'p_species(ps)%rab',ierr)
          end if

          if (p_species(ps)%core_charge) then
             deallocate(p_species(ps)%core_charge_recip,stat=ierr)
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  'p_species(ps)%core_charge_recip',ierr)
          end if

          deallocate(p_species(ps)%rad_proj_recip,stat=ierr)
          call utils_dealloc_check('pseudopotentials_species_exit', &
               'p_species(ps)%rad_proj_recip',ierr)

          deallocate(p_species(ps)%rad_locpot_recip,stat=ierr)
          call utils_dealloc_check('pseudopotentials_species_exit', &
               'p_species(ps)%rad_locpot_recip',ierr)

          deallocate(p_species(ps)%D0,stat=ierr)
          call utils_dealloc_check('pseudopotentials_species_exit', &
               'p_species(ps)%D0',ierr)

          deallocate(p_species(ps)%core_radius,stat=ierr)
          call utils_dealloc_check('pseudopotentials_species_exit', &
               'p_species(ps)%core_radius',ierr)

          deallocate(p_species(ps)%ang_mom,stat=ierr)
          call utils_dealloc_check('pseudopotentials_species_exit', &
               'p_species(ps)%ang_mom',ierr)

       end do

       deallocate(p_species,stat=ierr)
       call utils_dealloc_check('pseudopotentials_species_exit', &
            'p_species',ierr)

    end if

  end subroutine pseudopotentials_species_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! SUBROUTINES FOR INTERFACING WITH THE ATOM SOLVER !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_locpot_rad(locpot,npts,rad,isp,p_species)

    !=====================================================================!
    ! This subroutine fetches the local pseudopotential on a regular grid !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                         !
    ! References to par removed by Robert Charlton, 14/12/2016.           !
    !=====================================================================!

    use constants, only: stdout
    use services, only: services_sbessj
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(out) :: locpot(npts)
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp
    type(PSEUDO_SPECIES), intent(in) :: p_species(:) ! size(par%num_pspecies)

    ! Local Variables
    real(kind=DP), allocatable :: work(:)
    real(kind=DP) :: q,dq,Z
    integer :: npts_q
    integer :: ir,iq
    integer :: ierr

    npts_q = p_species(isp)%n_rad_pts
    allocate(work(npts_q),stat=ierr)
    call utils_alloc_check('pseudo_get_locpot_rad','work',ierr)

    dq = p_species(isp)%ps_gmax / real(npts_q-1,kind=DP)
    Z = p_species(isp)%ion_charge

    do ir=1,npts
       do iq=1,npts_q
          q = real(iq-1,kind=DP)*dq
          work(iq) = (p_species(isp)%rad_locpot_recip(iq)*q**2 - Z) * &
               services_sbessj(0,q*rad(ir))
       end do

       locpot(ir) = work(1)+work(npts_q)
       do iq = 2,npts_q-1,2
           locpot(ir) = locpot(ir) + 4.0_DP*work(iq)+2.0_DP*work(iq+1)
       enddo
       locpot(ir) = locpot(ir)*dq/3.0_DP
    end do

    locpot(:) = locpot(:)*(2.0_DP/PI)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('pseudo_get_locpot_rad','work',ierr)

  end subroutine pseudo_get_locpot_rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_core_den_rad(core_den,npts,rad,isp,p_species)

    !============================================================!
    ! This subroutine fetches the core density on a regular grid !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    use services, only: services_sbessj
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(out) :: core_den(npts)
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp
    type(PSEUDO_SPECIES), intent(in) :: p_species(:) ! size par%num_pspecies)

    ! Local variables
    real(kind=DP), allocatable :: work(:)
    real(kind=DP) :: q,dq
    integer :: npts_q
    integer :: ir,iq
    integer :: ierr

    if (.not.p_species(isp)%core_charge) then
       core_den(:) = 0.0_DP
       return
    end if

    npts_q = p_species(isp)%n_rad_pts
    allocate(work(npts_q),stat=ierr)
    call utils_alloc_check('pseudo_get_core_den_rad','work',ierr)

    dq = p_species(isp)%ps_gmax / real(npts_q-1,kind=DP)

    do ir=1,npts
       do iq=1,npts_q
          q = real(iq-1,kind=DP)*dq
          work(iq) = p_species(isp)%core_charge_recip(iq)*q**2 * &
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
    call utils_dealloc_check('pseudo_get_core_den_rad','work',ierr)

  end subroutine pseudo_get_core_den_rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_projector_info(isp,p_species,nshells,nproj_tot,lmax)

    !============================================================!
    ! This function fetches the number of shells of projectors.  !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    ! Arguments
    integer,intent(in) :: isp
    type(PSEUDO_SPECIES), intent(in) :: p_species(:) ! size(par%num_pspecies)
    integer,intent(out) :: nshells
    integer,intent(out) :: nproj_tot
    integer,intent(out) :: lmax

    ! Local Variables
    integer :: ishell

    nshells = p_species(isp)%n_shells
    nproj_tot = 0
    do ishell=1,nshells
       nproj_tot = nproj_tot + 2*p_species(isp)%ang_mom(ishell) + 1
    end do
    if (nshells>0) then
       lmax = maxval(p_species(isp)%ang_mom)
    else
       lmax = 0
    end if

  end subroutine pseudo_get_projector_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_aug_funcs(qijl,qfunc,nshells,npts,rad,isp,p_species)

    !============================================================!
    ! This subroutine fetches the q function and the qijL        !
    ! terms for each angular momentum channel.                   !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in February 2011.                 !
    !============================================================!

    use services, only: services_locate_interp, services_linear_interpolation

    ! Arguments
    integer, intent(in) :: npts
    integer, intent(in) :: nshells
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp
    type(PSEUDO_SPECIES), intent(in) :: p_species(:)! size(par%num_pspecies)
    real(kind=DP), intent(out) :: qijL(nshells,nshells)
    real(kind=DP), intent(out) :: qfunc(npts,nshells,nshells)

    ! Local variables
    integer :: ipt,jpt
    integer :: ishell,jshell

    ! Set the augmentation charges
    qijL(:,:) = p_species(isp)%aug_q(:,:)

    if (p_species(isp)%kkbeta==0) then
       qfunc(:,:,:) = 0.0_DP
       return
    end if

    ! Set the L-independent function
    do ishell=1,nshells
       do jshell=1,nshells
          do ipt=1,npts
             if (rad(ipt)>p_species(isp)%rlog(p_species(isp)%kkbeta-1)) then
                qfunc(ipt,ishell,jshell) = 0.0_DP
                cycle
             end if
             jpt = services_locate_interp(rad(ipt), &
                  p_species(isp)%rlog,p_species(isp)%mesh)
             qfunc(ipt,ishell,jshell) = services_linear_interpolation(rad(ipt), &
                  p_species(isp)%qfunc(jpt,ishell,jshell), &
                  p_species(isp)%qfunc(jpt+1,ishell,jshell), &
                  p_species(isp)%rlog(jpt),p_species(isp)%rlog(jpt+1))
          end do
       end do
    end do

  end subroutine pseudo_get_aug_funcs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_projectors_q(proj_q,dij,ang_mom,nshells,nsws, &
       nsws_max,lmax,qb,isp,p_species)

    !============================================================!
    ! This subroutine fetches the projectors at chosen q-points  !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    use services, only: services_1d_interpolation

    implicit none

    ! Arguments
    integer, intent(in) :: nshells
    integer, intent(in) :: lmax
    integer, intent(in) :: nsws(0:lmax)
    integer, intent(in) :: nsws_max
    integer, intent(in) :: isp
    type(PSEUDO_SPECIES), intent(in) :: p_species(:) ! size(par%num_pspecies)
    integer, intent(out) :: ang_mom(nshells)
    real(kind=DP), intent(in) :: qb(nsws_max,0:lmax)
    real(kind=DP), intent(out) :: proj_q(nsws_max,nshells)
    real(kind=DP), intent(out) :: dij(nshells,nshells)

    ! Local variables
    integer :: ishell, isw
    real(kind=DP) :: qq

    ! Set the KB denominators
    dij(:,:) = p_species(isp)%D0(:,:)

    ! Interpolate the reciprocal-space projectors at the qb values
    do ishell=1,nshells
       proj_q(:,ishell) = 0.0_DP
       ang_mom(ishell) = p_species(isp)%ang_mom(ishell)
       do isw=1,nsws(ang_mom(ishell))
          qq = qb(isw,p_species(isp)%ang_mom(ishell))
          proj_q(isw,ishell) = services_1d_interpolation( &
               p_species(isp)%rad_proj_recip(:,ishell), &
               p_species(isp)%n_rad_pts, &
               qq*p_species(isp)%inv_g_spacing, &
               p_species(isp)%ang_mom(ishell))
       end do
    end do

  end subroutine pseudo_get_projectors_q



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! SUBROUTINES FOR THE LOCAL PSEUDOPOTENTIAL !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_make_structure_factor(struct_fac,    &  ! output
       elements, grid, num_pspecies,                     &  ! input
       directions1, weights1, directions2, weights2)        ! input (optional)

    !==============================================================!
    ! This subroutine generates the structure factor each distinct !
    ! ionic pseudopotential species.                               !
    !--------------------------------------------------------------!
    ! Written by Arash A. Mostofi 19/2/2004.                       !
    ! Modified by Chris-Kriton Skylaris on 22/2/2004 so that       !
    ! struct_fac is parallelised in memory.                        !
    ! Re-written by Peter D. Haynes on 30/6/2004 to use data       !
    ! parallelisation according to the new Fourier routines.       !
    ! Modified by Nicholas Hine on 19/09/2008 to optimise cache    !
    ! performance by re-ordering array                             !
    ! Added first-order and second-order capabilities by Gabriel   !
    ! Constantinescu in Nov 2016                                   !
    ! Modified to remove pub_par by Robert Charlton, 22/08/2018.   !
    !==============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_barrier, pub_my_proc_id
    use constants, only: cmplx_0
    use ion, only: ELEMENT
!$  use rundat, only: pub_threads_max
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)   :: grid
    type(ELEMENT), intent(in)     :: elements(:) ! size(nat)
    integer, intent(in)           :: num_pspecies
    complex(kind=DP), intent(out) :: struct_fac(num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    integer, intent(in), optional :: directions1(:)    ! size(nat)
    real(kind=DP), intent(in), optional :: weights1(:) ! size(nat)
    integer, intent(in), optional :: directions2(:)    ! size(nat)
    real(kind=DP), intent(in), optional :: weights2(:) ! size(nat)

    ! Local variables
    integer :: ipt,i3,i2,islab23
    integer :: species, atom, atom_sp, nat
    integer :: ierr
    real(kind=DP) :: gpt(3)
    real(kind=DP),allocatable :: species_pos(:,:,:)
    complex(kind=DP),allocatable :: phase_components(:,:,:,:)
    integer, allocatable :: species_nat(:)
    real(kind=DP) :: g1r1, g2r2, g3r3
    complex(kind=DP) :: accum

    ! Start timer
    call timer_clock('pseudo_make_structure_factor',1)

    ! rc2013: get size of elements and pspecies arrays
    nat = size(elements)

    ! jd: Take care of padding, the loop ignores it
    struct_fac = cmplx_0

    ! ndmh: allocate storage for number of atoms in each species
    allocate(species_nat(num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_make_structure_factor','species_nat',ierr)

    ! ndmh: count number of atoms in each species
    species_nat = 0
    do atom=1,nat
       species = elements(atom)%pspecies_number
       species_nat(species) = species_nat(species) + 1
    end do

    ! ndmh: allocate temporary storage for species atom positions
    allocate(species_pos(3,maxval(species_nat),num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_make_structure_factor','species_pos',ierr)

    ! jd: allocate temporary for components of the phase factor
    allocate(phase_components(3,maxval(species_nat),num_pspecies, &
         max(grid%num_slabs23,grid%n2,grid%n3)),stat=ierr)
    call utils_alloc_check('pseudo_make_structure_factor','phase_components', &
         ierr)

    ! ndmh: copy atom positions to temporary array for cache-efficiency
    species_nat = 0
    do atom=1,nat
       species = elements(atom)%pspecies_number
       species_nat(species) = species_nat(species) + 1
       species_pos(1,species_nat(species),species) = elements(atom)%centre%x
       species_pos(2,species_nat(species),species) = elements(atom)%centre%y
       species_pos(3,species_nat(species),species) = elements(atom)%centre%z
    end do

    ! jd: Determine exp(-i g1*r1) for all atoms
    do islab23=1,grid%num_slabs23                ! along b1
       do species=1,num_pspecies
          do atom=1,species_nat(species)
             call cell_grid_recip_pt(gpt(1:3),islab23 + &
                  grid%first_slab23(pub_my_proc_id) - 1,1,1,grid)
             g1r1 = sum(gpt(1:3)*species_pos(1:3,atom,species))
             phase_components(1,atom,species,islab23) = &
                  cmplx(cos(g1r1),-sin(g1r1),kind=DP)
          end do
       end do
    end do

    ! jd: Determine exp(-i g2*r2) for all atoms
    do i2=1,grid%n2                              ! along b2
       do species=1,num_pspecies
          do atom=1,species_nat(species)
             call cell_grid_recip_pt(gpt(1:3),1,i2,1,grid)
             g2r2 = sum(gpt(1:3)*species_pos(1:3,atom,species))
             phase_components(2,atom,species,i2) = &
                  cmplx(cos(g2r2),-sin(g2r2),kind=DP)
          end do
       end do
    end do

    ! jd: Determine exp(-i g3*r3) for all atoms
    do i3=1,grid%n3                             ! along b3
       do species=1,num_pspecies
          do atom=1,species_nat(species)
             call cell_grid_recip_pt(gpt(1:3),1,1,i3,grid)
             g3r3 = sum(gpt(1:3)*species_pos(1:3,atom,species))
             phase_components(3,atom,species,i3) = &
                  cmplx(cos(g3r3),-sin(g3r3),kind=DP)
          end do
       end do
    end do

    ! Loop over reciprocal space on this proc
    ! jd: Use exp(-i g.r) = exp(-i g1*r1) * exp(-i g2*r2) * exp(-i g3*r3)
    !     to avoid calculating trigs every time.

!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i2,i3,islab23,species,gpt,atom,accum,atom_sp) &
!$OMP SHARED(species_nat,phase_components,struct_fac,grid,nat,num_pspecies, &
!$OMP      pub_threads_max,weights1,directions1,weights2,directions2, &
!$OMP      pub_my_proc_id, elements)
    do ipt=1,grid%num_slabs23*grid%n2*grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1

       call cell_grid_recip_pt(gpt(1:3),islab23 + &
            grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)

       do species=1,num_pspecies
          accum = cmplx_0
          atom_sp = 0
          do atom = 1, nat
             if (species == elements(atom)%pspecies_number) then
                atom_sp = atom_sp + 1
                if (present(directions1)) then
                   if (present(directions2)) then ! second_order
                      accum = accum + &
                           phase_components(1,atom_sp,species,islab23) * &
                           phase_components(2,atom_sp,species,i2) * &
                           phase_components(3,atom_sp,species,i3) * &
                           cmplx(-gpt(directions1(atom)) * &
                           weights1(atom) * gpt(directions2(atom)) * &
                           weights2(atom), 0.0_DP, kind=DP)
                   else      ! first-order
                      accum = accum + &
                           phase_components(1,atom_sp,species,islab23) * &
                           phase_components(2,atom_sp,species,i2) * &
                           phase_components(3,atom_sp,species,i3) * &
                           cmplx(0.0, -gpt(directions1(atom)) * &
                           weights1(atom), kind=DP)
                   end if
                else ! zero-order
                   accum = accum + &
                        phase_components(1,atom_sp,species,islab23) * &
                        phase_components(2,atom_sp,species,i2) * &
                        phase_components(3,atom_sp,species,i3)
                end if
             end if
          end do
          struct_fac(species,i3,i2,islab23) = accum
       end do
    end do
!$OMP END PARALLEL DO

    ! ndmh: this routine has no comms, so is not necessarily synchronised.
    ! ndmh: synchronise here before returning, otherwise timers report
    ! ndmh: inaccurate timings for next function called, due to load imbalance
    ! ndmh: in this routine
    call comms_barrier

    ! ndmh: deallocate temporary storage
    deallocate(species_pos,stat=ierr)
    call utils_dealloc_check('pseudo_make_structure_factor','species_pos',ierr)
    deallocate(species_nat,stat=ierr)
    call utils_dealloc_check('pseudo_make_structure_factor','species_nat',ierr)
    deallocate(phase_components,stat=ierr)
    call utils_dealloc_check('pseudo_make_structure_factor','phase_components',ierr)

    ! Stop timer
    call timer_clock('pseudo_make_structure_factor',2)

  end subroutine pseudo_make_structure_factor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_local_on_grid(v_local_fine,locps_gzero_term, &
        struct_fac,struct_fac_classical,grid,cell,p_species,par,nat_classical)

    !==================================================================!
    ! This subroutine generates the local part of the total ionic      !
    ! pseudopotential on the fine grid simulation cell.                !
    !------------------------------------------------------------------!
    !                                                                  !
    !------------------------------------------------------------------!
    ! This subroutine was originally written by Chris-Kriton Skylaris  !
    ! in 2000.                                                         !
    ! Modified by A. A. Mostofi in 19/02/2004 so that it uses          !
    ! complex-complex FFTs.                                            !
    ! Modified by Chris-Kriton Skylaris on 22/2/2004 so that           !
    ! struct_fac is memory-parallelised                                !
    ! Modified for parallelisation with new fourier routines by        !
    ! Peter D. Haynes on 30/6/2004.                                    !
    ! Modified to include smeared-ion contribution if required, which  !
    ! necessitated passing elements, by Jacek Dziedzic on 13/05/2010.  !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use classical_pot, only: classical_pot_recip
    use comms, only: comms_bcast, pub_my_proc_id, pub_root_proc_id
    use fourier, only: fourier_apply_cell_backward
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(GRID_INFO), intent(in) :: grid
    type(PARAL_INFO), intent(in) :: par
    complex(kind=DP), intent(in) :: struct_fac(par%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    complex(kind=DP), intent(in) :: struct_fac_classical( &
         grid%ld3,grid%ld2,grid%max_slabs23)
    real(kind=DP), intent(out) :: v_local_fine(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    real(kind=DP), intent(out)    :: locps_gzero_term
    type(PSEUDO_SPECIES), intent(in) :: p_species(par%num_pspecies)
    integer, intent(in)         :: nat_classical

    ! Local variables
    integer :: ierr              ! Error flag
    complex(kind=DP), allocatable :: v_local_recip(:,:,:)

    call timer_clock('pseudopotentials_local_on_grid',1)

    ! jd: Usual PBC method

    ! cks: expand to 3D in reciprocal fine grid and sum together
    !      (along with the structure factor) the local pseudopotentials
    !      for each species to obtain the total local ionic potential in
    !      reciprocal representation on the fine grid
    allocate(v_local_recip(grid%ld3,grid%ld2, &
         grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pseudopotentials_local_on_grid', &
         'v_local_recip',ierr)

    call pseudopotentials_sum_local_rec(v_local_recip, & ! output
         struct_fac,grid,cell,p_species,par)    ! input

    ! store g=0 term seperately for charged system correction
    locps_gzero_term = 0.0_DP
    if (pub_my_proc_id==grid%proc_slab23(1)) then
       locps_gzero_term = &
            real(v_local_recip(1,1,1) / (grid%n1 * grid%n2 * grid%n3),kind=DP)
            ! jd: Wrapped to work around warning on complex to real assignment
    end if
    call comms_bcast(pub_root_proc_id,locps_gzero_term)

    ! cks: include external potential from "classical" atoms
    if (nat_classical > 0) then
       call classical_pot_recip(v_local_recip, & ! output
            struct_fac_classical,grid)    ! input
    endif

    ! FFT the local ionic potential from reciprocal to real space
    call fourier_apply_cell_backward(v_local_fine,v_local_recip,grid)

    deallocate(v_local_recip,stat=ierr)
    call utils_dealloc_check('pseudopotentials_local_on_grid', &
         'v_local_recip',ierr)


    call timer_clock('pseudopotentials_local_on_grid',2)

  end subroutine pseudopotentials_local_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_sum_local_rec(fine_complex,struct_fac,grid,cell, &
       p_species, par)

    !===============================================================!
    ! This subroutine generates in reciprocal space the total local !
    ! potential in the simulation cell due to the ionic             !
    ! pseudopotentials of all ions.                                 !
    !---------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris in 2000.          !
    ! Improvements by Chris-Kriton Skylaris on 30/5/2001.           !
    ! Modified by Chris-Kriton Skylaris on 24/1/2004 to work with   !
    ! the pseudo_species type.                                      !
    ! Modified by Arash A. Mostofi in Aprl 2003 to be compatible    !
    ! with complex-to-complex FFTs (rather than complex-to-real).   !
    ! Modified by Arash A. Mostofi on 19/2/2004 to work with        !
    ! structure factor as input.                                    !
    ! Modified by Chris-Kriton Skylaris on 22/2/2004 to work with   !
    ! data-parallel structure factor.                               !
    ! Modified for parallelisation with new fourier routines by     !
    ! Peter D. Haynes on 30/6/2004.                                 !
    ! Modified by Nicholas Hine in December 2009 to allow less than !
    ! one slab23 per proc.                                          !
    ! OpenMP parallelised by Nicholas Hine in June 2013             !
    !===============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use pbc_corrections, only: pbc_corr_vloc, pbc_corr_initialise, &
         pbc_corr_is_initialised
    use rundat, only: pub_mt_cutoff, pub_smooth_loc_pspot
!$  use rundat, only: pub_threads_max
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in)   :: cell
    type(GRID_INFO), intent(in)   :: grid
    type(PARAL_INFO), intent(in)  :: par
    complex(kind=DP), intent(in)  :: struct_fac(par%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    complex(kind=DP), intent(out) :: fine_complex(grid%ld3,&
         grid%ld2, grid%max_slabs23)
    type(PSEUDO_SPECIES), intent(in) :: p_species(par%num_pspecies)

    ! Local variables
    integer :: species              ! Atomic species counter
    integer :: ipt,i3,i2,islab23    ! Reciprocal grid loop counters
    real(kind=DP) :: gvec(3)        ! G-vector
    real(kind=DP) :: g_length       ! Length of this G vector
    real(kind=DP) :: v_loc_value    ! Local potential at this G
    real(kind=DP) :: fourpi         ! Constant multiplier
    real(kind=DP) :: g_cut, alpha   ! For filtering locps

    fourpi = 4.0_DP * PI

    if (pub_smooth_loc_pspot > 0.0_DP) then
       g_cut = 3.0_DP*PI/(cell%d1+cell%d2+cell%d3)*(4.0_DP/7.0_DP)
       alpha = 0.1_DP*log(2.0_DP)/((pub_smooth_loc_pspot*g_cut)**2)
    else
       g_cut = huge(1.0_DP)
       alpha = 0.0_DP
    end if

    ! ndmh: Initialise pbc corrections if not already done so
    if (pub_mt_cutoff /= 0.0_DP) then
       if (.not. pbc_corr_is_initialised) call pbc_corr_initialise(grid,cell)
    end if

    ! Loop over reciprocal space grid on this proc
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i2,i3,islab23,species,gvec,g_length,v_loc_value) &
!$OMP SHARED(grid,fine_complex,p_species,struct_fac,par, &
!$OMP      g_cut,alpha,pub_mt_cutoff,fourpi,pub_my_proc_id,pub_threads_max)
    do ipt=1,grid%num_slabs23*grid%n2*grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1

       ! Get G-vector
       call cell_grid_recip_pt(gvec,islab23 + &
            grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)

       ! Get magnitude of this G-vector
       g_length = sqrt(sum(gvec(1:3)**2))

       ! ndmh: initialise (cache-efficiently)
       fine_complex(i3,i2,islab23) = cmplx(0.0_DP,0.0_DP,kind=DP)

       ! Loop over atomic species
       do species=1,par%num_pspecies

          ! Get pseudopotential at this G-vector
          v_loc_value = services_1d_interpolation( &
               p_species(species)%rad_locpot_recip, &
               p_species(species)%n_rad_pts,&
               g_length*p_species(species)%inv_g_spacing,0)

          ! Filter high G components if required
          if (g_length > g_cut) then
             v_loc_value = v_loc_value * exp(-alpha*(g_length-g_cut)**2)
          end if

          ! Add back Coulomb potential
          v_loc_value = v_loc_value - p_species(species)%ion_charge * &
                grid%coulomb_recip(i3,i2,islab23)

          ! ndmh: scale by 4 pi/fine_weight
          ! jd: split the scaling into two steps
          v_loc_value = v_loc_value * fourpi

          ! jd: Apply the Martyna-Tuckerman correction, if requested
          if (pub_mt_cutoff /= 0.0_DP) then
             call pbc_corr_vloc(i3,i2,islab23, &
                  p_species(species)%ion_charge,v_loc_value)
          end if

          ! jd: remaining part of scaling
          v_loc_value = v_loc_value / grid%weight

          fine_complex(i3,i2,islab23) = fine_complex(i3,i2,islab23) + &
               struct_fac(species,i3,i2,islab23) * v_loc_value

       end do    ! loop over species

    end do
!$OMP END PARALLEL DO


    ! G=0 element must be real
    if (pub_my_proc_id==grid%proc_slab23(1)) then
       if (aimag(fine_complex(1,1,1)) /= 0.0_DP) then
          call utils_abort('Error in pseudopotentials_sum_local_rec: &
               &potential not real')
       end if
    end if

    if (grid%num_slabs23 > 0) then
       ! Nyquist filter (fine grid is always going to be even)
       fine_complex(grid%n3/2+1,:,:) = (0.0_DP,0.0_DP)
       fine_complex(:,grid%n2/2+1,:) = (0.0_DP,0.0_DP)
    end if
    ! Nyquist filter for last slab
    if (pub_my_proc_id==grid%proc_slab23(grid%n1/2+1)) &
         fine_complex(:,:,grid%num_slabs23) = (0.0_DP,0.0_DP)

  end subroutine pseudopotentials_sum_local_rec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_core_density(core_density_fine,struct_fac, &
       p_species,grid,par)

    !==================================================================!
    ! This subroutine reconstructs the core density for the whole      !
    ! supercell on the fine real space grid using the information      !
    ! stored for each pseudopotential species in the array             !
    ! core_charge_recip.                                               !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !                                                                  !
    ! 1) core_density_fine  : output  : data-parallelised core density !
    ! 2) struct_fac : input : data-parallelised structure factor       !
    ! 3) grid : input : description of whole-cell grid                 !
    !                                                                  !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 29/01/09.           !
    ! OpenMP parallelised by Nicholas Hine in June 2013                !
    !==================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_proc_id
    use fourier, only: fourier_apply_cell_backward
    use parallel_strategy, only: PARAL_INFO
!$  use rundat, only: pub_threads_max
    use services, only: services_1d_interpolation
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(PARAL_INFO), intent(in) :: par
    real(kind=DP), intent(out) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    complex(kind=DP), intent(in) :: struct_fac(par%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    type(PSEUDO_SPECIES), intent(in) :: p_species(par%num_pspecies)

    ! Local variables
    integer :: ierr              ! Error flag
    complex(kind=DP), allocatable :: core_density_recip(:,:,:)
    real(kind=DP) :: gvec(3), g_length, core_den_value, factor
    integer :: species              ! Atomic species counter
    integer :: ipt,i3,i2,islab23    ! Reciprocal grid loop counters

    call timer_clock('pseudopotentials_core_density',1)

    ! ndmh: allocate storage for core density in reciprocal space
    allocate(core_density_recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pseudopotentials_core_density', &
         'core_density_recip',ierr)

    ! ndmh: loop over reciprocal space grid on this proc
!$OMP PARALLEL DO NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i2,i3,islab23,species,gvec,g_length,core_den_value) &
!$OMP SHARED(grid,core_density_recip,p_species,struct_fac,par, &
!$OMP      pub_my_proc_id,pub_threads_max)
    do ipt=1,grid%num_slabs23*grid%n2*grid%n3

       i3 = modulo(ipt-1,grid%n3) + 1
       i2 = modulo((ipt-i3)/grid%n3,grid%n2) + 1
       islab23 = (ipt-(i2-1)*grid%n3-i3) / (grid%n3*grid%n2) + 1

       ! Initialise
       core_density_recip(i3,i2,islab23) = (0.0_DP,0.0_DP)

       ! Loop over atomic species
       do species=1,par%num_pspecies

          ! Check if we have a core charge for this species
          if (.not.p_species(species)%core_charge) cycle

          ! Get length of this G-vector
          call cell_grid_recip_pt(gvec,islab23 + &
               grid%first_slab23(pub_my_proc_id) - 1,i2,i3,grid)
          g_length = sqrt(sum(gvec(1:3)**2))

          ! Get core density at this G-vector
          core_den_value = services_1d_interpolation( &
               p_species(species)%core_charge_recip, &
               p_species(species)%n_rad_pts,&
               g_length*p_species(species)%inv_g_spacing,0)

          core_density_recip(i3,i2,islab23) = &
               core_density_recip(i3,i2,islab23) + &
               struct_fac(species,i3,i2,islab23) * core_den_value

       end do    ! loop over species

    end do
!$OMP END PARALLEL DO

    ! FFT the core density from reciprocal to real space
    call fourier_apply_cell_backward(core_density_fine,core_density_recip,grid)

    ! ndmh: scale with 1.0/weight
    factor = 1.0_DP / grid%weight
    core_density_fine = factor * core_density_fine

    ! ndmh: deallocate storage for core density in reciprocal space
    deallocate(core_density_recip,stat=ierr)
    call utils_dealloc_check('pseudopotentials_core_density', &
         'core_density_recip',ierr)

    call timer_clock('pseudopotentials_core_density',2)

  end subroutine pseudopotentials_core_density


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pseudopotentials_FO_nlcc_energy(xc_pot,grid,cell,elements, &
       pseudo_sp,nlcc_energy,directions1,weights1,directions2,weights2)

    use cell_grid, only: cell_grid_recip_pt, GRID_INFO
    use comms, only: comms_reduce, pub_my_proc_id
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: par=>pub_par
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
    real(kind=DP), intent(in) :: xc_pot(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_num_spins)
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(PSEUDO_SPECIES), intent(in) :: pseudo_sp(par%num_pspecies)
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
         'DEBUG: Entering pseudopotentials_FO_nlcc_energy'

    ! Start timer
    call timer_clock('pseudopotentials_FO_nlcc_energy',1)

    ! Allocate
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pseudopotentials_FO_nlcc_energy','recip',ierr)

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

!!!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!!!$OMP PRIVATE(i2,i3,gvec,g_length,coreden, &
!!!$OMP      iG_coreden_vxc,factor,islab23,gdotR,eiGR,ierr,species) &
!!!$OMP SHARED (pub_my_proc_id,grid,pseudo_sp,elements,recip,par, &
!!!$OMP      pub_threads_max,directions1,weights1,directions2,weights2) &
!!!$OMP REDUCTION (+:nlcc_energy)

    ! For components g with symmetry points at -g
    factor=2.0_DP

!!!$OMP DO
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
          if (.not.pseudo_sp(species)%core_charge) cycle

          ! Interpolate value of core density at current g
          coreden = services_1d_interpolation( &
               pseudo_sp(species)%core_charge_recip, &
               pseudo_sp(species)%n_rad_pts, &
               g_length*pseudo_sp(species)%inv_g_spacing,0)

          do atom = 1, par%nat
             if (species == elements(atom)%pspecies_number) then

                ! calculate -iG_cart1.n(G).vxc*(G) if FO (dir2 NOT present)
                ! calculate -G_cart1.G_cart2.n(G).vxc*(G) if SO (dir2 present)

                if (.not.present(directions2)) then
                   iG_coreden_vxc = factor * cmplx(0.0_DP,-1.0_DP,kind=DP) * &
                        gvec(directions1(atom)) * coreden * &
                        conjg(recip(i3,i2,islab23)) * weights1(atom)
                else
                   iG_coreden_vxc = cmplx(-1.0_DP,0.0_DP,kind=DP) * factor * &
                        gvec(directions1(atom)) * gvec(directions2(atom)) * &
                        coreden * conjg(recip(i3,i2,islab23)) * &
                        weights1(atom) * weights2(atom)
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
          end do ! over atoms
       end do ! over species

       ! For g1/=0 slabs
       factor=2.0_DP

    end do
!!!$OMP END DO

!!!$OMP END PARALLEL

    call comms_reduce('SUM',nlcc_energy)

    ! Deallocate
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('pseudopotentials_FO_nlcc_energy','recip',ierr)

    ! Scale factor
    nlcc_energy = nlcc_energy / real(grid%n1*grid%n2*grid%n3,dp)

    ! Stop timer
    call timer_clock('pseudopotentials_FO_nlcc_energy',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudopotentials_FO_nlcc_energy'

    return
  end subroutine pseudopotentials_FO_nlcc_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pseudo_local_calculate_forces(den_slabs12,grid,cell,elements, &
       p_species,locps_forces,par)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the local part of the ionic pseudopotential.                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) den_slabs12  : in/out : ground state data-parallelised charge density!
    !                            in/out as spin normalisation is done in situ !
    ! 2) grid         : input  : the grid on which the density is represented !
    ! 3) cell         : input  : the simulation cell                          !
    ! 4) elements     : input  : list of elements and corresponding info      !
    ! 5) p_species    : input  : pseudo information of all ionic species      !
    ! 6) locps_forces : output : ionic forces due to local part of pseudopot  !
    !-------------------------------------------------------------------------!
    ! Written by Arash A. Mostofi, v0.0, 28th June 2004                       !
    ! Modified by Nicholas Hine on 04/09/2009 to re-order the loops for cache !
    ! efficiency and to pre-calculate v(G) for each species for each G vector.!
    ! Modified by Nicholas Hine in December 2009 to allow for less than one   !
    ! 23-slab per processor.                                                  !
    ! OpenMP parallelised by Nicholas Hine in June 2013                       !
    ! Fixed to correctly support Martyna-Tuckerman correction, jd, 09.2013.   !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_reduce, pub_my_proc_id
    use constants, only: pi
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use pbc_corrections, only: pbc_corr_vloc
    use rundat, only : pub_smooth_loc_pspot, pub_mt_cutoff, &
         pub_num_spins
!$  use rundat, only: pub_threads_max
    use services, only: services_1d_interpolation
    use simulation_cell, only: CELL_INFO
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in) :: cell
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: den_slabs12(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)
    type(PARAL_INFO), intent(in) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(PSEUDO_SPECIES), intent(in) :: p_species(par%num_pspecies)
    real(kind=DP), intent(out) :: locps_forces(1:3,par%nat)

    ! Local Variables
    integer :: ipt,i2,i3,islab23
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gvec(3)
    real(kind=DP) :: g_length,gdotR
    real(kind=DP) :: factor,v_loc_value
    real(kind=DP) :: g_cut
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: iG_vden(:,:)
    complex(kind=DP), allocatable :: recip(:,:,:)
    logical :: filter_locps
    real(kind=DP), parameter :: fourpi = 4.0_DP*PI

    ! -------------------------------------------------------------------------


    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_local_calculate_forces'

    ! Start timer
    call timer_clock('pseudo_local_calculate_forces',1)

    ! Allocate
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pseudo_local_calculate_forces','recip',ierr)

    ! Initialise
    locps_forces = 0.0_DP

    if (pub_smooth_loc_pspot > 0.0_DP) then
       filter_locps = .true.
       g_cut = 3.0_DP*PI/(cell%d1+cell%d2+cell%d3)*(4.0_DP/7.0_DP)
    else
       filter_locps = .false.
       g_cut = huge(1.0_DP)
    end if

    ! If spin polarised, put total density in up spin
    if (pub_num_spins == 2) den_slabs12(:,:,:,1) = den_slabs12(:,:,:,1) + &
         den_slabs12(:,:,:,2)

    ! Fourier transform the charge density to reciprocal space
    call fourier_apply_cell_forward(den_slabs12(:,:,:,1),recip,grid)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i2,i3,islab23,gvec,g_length,species,v_loc_value, &
!$OMP      iG_vden,factor,atom,gdotR,eiGR,ierr) &
!$OMP SHARED (pub_my_proc_id,grid,p_species,elements,recip, &
!$OMP      filter_locps,g_cut,pub_smooth_loc_pspot,pub_threads_max,par, &
!$OMP      pub_mt_cutoff) &
!$OMP REDUCTION (+:locps_forces)

    allocate(iG_vden(1:3,par%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_local_calculate_forces','iG_vden',ierr)

    ! For components g with symmetry points at -g
    factor=2.0_DP

    ! Loop over all g points
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

       ! ndmh: loop over species, calculating v(G) for each.
       ! ndmh: this bit is O(N) so as much as possible should be
       ! ndmh: pre-calculated here.
       do species=1,par%num_pspecies

          ! Interpolate value of local potential at current g
          v_loc_value = services_1d_interpolation( &
               p_species(species)%rad_locpot_recip, &
               p_species(species)%n_rad_pts,g_length * &
               p_species(species)%inv_g_spacing,0)

          ! Filter high G components if required
          if (filter_locps.and.(g_length > g_cut)) then
             v_loc_value = v_loc_value * &
                  exp(-pub_smooth_loc_pspot*(g_length/g_cut-1.0_DP)**2)
          end if

          ! Add back the Coulomb potential; set g=0 term to zero
          if (g_length .gt. 0.0_DP) then
             ! pa: changed to allow fractional ionic charge
             v_loc_value = v_loc_value - &
                  p_species(species)%ion_charge * &
                  grid%coulomb_recip(i3,i2,islab23)
          else
             v_loc_value = 0.0_DP
          endif

          ! jd: Apply the Martyna-Tuckerman correction, if requested
          if (pub_mt_cutoff /= 0.0_DP) then
             v_loc_value = v_loc_value * fourpi
             call pbc_corr_vloc(i3,i2,islab23, &
                  p_species(species)%ion_charge,v_loc_value)
             v_loc_value = v_loc_value / fourpi
          end if

          ! ndmh: calculate iG.v(G).n*(G) for each species and
          ! ndmh: each Cartesian direction
          iG_vden(1:3,species) = factor * cmplx(0.0_DP,1.0_DP,kind=DP) &
               * gvec * v_loc_value * conjg(recip(i3,i2,islab23))

       end do

       ! ndmh: loop over atoms: this bit is O(N^2) so should be made as
       ! ndmh: short and simple as possible
       do atom=1,par%nat

          ! e^{-ig.R}
          gdotR = -(gvec(1)*elements(atom)%centre%x + &
               gvec(2)*elements(atom)%centre%y + &
               gvec(3)*elements(atom)%centre%z)
          eiGR = cmplx(COS(gdotR),SIN(gdotR),kind=DP)

          ! Sum force over g in each Cartesian direction i
          ! ==>  f = sum_{g} i.g.e^{-ig.R}.Vlocps(g).den^{*}(g)
          locps_forces(:,atom) = locps_forces(:,atom) + &
               real(iG_vden(:,elements(atom)%pspecies_number)*eiGR,kind=DP)

       enddo     ! End loop over atoms

       ! For g1/=0 slabs
       factor=2.0_DP

    end do
!$OMP END DO

    deallocate(iG_vden,stat=ierr)
    call utils_dealloc_check('pseudo_local_calculate_forces','iG_vden',ierr)
!$OMP END PARALLEL

    ! Sum the result over all procs
    call comms_reduce('SUM',locps_forces)

    ! Deallocate
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('pseudo_local_calculate_forces','recip',ierr)

    ! Scale
    locps_forces = locps_forces * 4.0_DP * PI / (grid%n1*grid%n2*grid%n3)

    ! If spin polarised, restore up spin density
    if (pub_num_spins == 2) den_slabs12(:,:,:,1) = den_slabs12(:,:,:,1) - &
         den_slabs12(:,:,:,2)

    ! Stop timer
    call timer_clock('pseudo_local_calculate_forces',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_local_calculate_forces'

    return
  end subroutine pseudo_local_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_nlcc_calculate_forces(density_fine,core_density_fine, &
       grid,cell,elements,p_species,nlcc_forces,par,active_density_fine)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the nonlinear core correction core charge.                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) density_fine : input  : ground state charge density                  !
    ! 2) core_density_fine  : input  : NLCC core charge                       !
    ! 3) elements     : input  : list of elements and corresponding info      !
    ! 4) nlcc_forces : output : ionic forces due to NLCC corrections          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) p_species    : pseudopotential information of all ionic species      !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 6th February 2009                             !
    ! based on pseudo_local_calculate_forces by Arash Mostofi                 !
    ! Modified by Nicholas Hine on 04/09/2009 to re-order the loops for cache !
    ! efficiency and to pre-calculate n(G) for each species for each G vector.!
    ! Modified by Nicholas Hine in December 2009 to allow for less than one   !
    ! 23-slab per processor.                                                  !
    ! OpenMP parallelised by Nicholas Hine in June 2013                       !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_reduce, pub_my_proc_id
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: PARAL_INFO
    use rundat, only: pub_num_spins, pub_spin_fac, pub_emft
!$  use rundat, only: pub_threads_max
    use simulation_cell, only: CELL_INFO
    use services, only: services_1d_interpolation
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert
    use xc, only : xc_energy_potential,xc_emft_calculate

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(PARAL_INFO), intent(in) :: par
    real(kind=DP), intent(inout) :: density_fine(grid%ld1,grid%ld2,&
         grid%max_slabs12,pub_num_spins)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(PSEUDO_SPECIES), intent(in) :: p_species(par%num_pspecies)
    real(kind=DP), intent(out) :: nlcc_forces(1:3,par%nat)
    ! jcap: if doing emft, active subsystem density
    real(kind=dp), intent(in), optional :: active_density_fine(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_num_spins)

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
    real(kind=DP), allocatable :: total_density_fine(:,:,:,:)
    real(kind=DP), allocatable :: total_active_density_fine(:,:,:,:)
    real(kind=DP), allocatable :: xc_pot_fine(:,:,:,:)
    real(kind=DP), allocatable :: xc_pot_fine_emft(:,:,:,:)
    complex(kind=DP), allocatable :: iG_coreden_vxc(:,:)
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: recip(:,:,:)

    logical :: loc_emft

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nlcc_calculate_forces'

    ! Start timer
    call timer_clock('pseudo_nlcc_calculate_forces',1)

    ! jcap: check to see if we are doing an emft calculation and if
    ! this matches pub_emft
    loc_emft=present(active_density_fine)
    call utils_assert(loc_emft.eqv.pub_emft,'Arguments of &
         &pseudo_nlcc_calculate_forces inconsistent with pub_emft:',pub_emft)

    ! Allocate
    allocate(total_density_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('pseudo_nlcc_calculate_forces','total_density_fine',&
         ierr)
    allocate(xc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_num_spins),stat=ierr)
    call utils_alloc_check('pseudo_nlcc_calculate_forces','xc_pot_fine',ierr)
    if (loc_emft) then
       allocate(total_active_density_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
            pub_num_spins),stat=ierr)
       call utils_alloc_check('pseudo_nlcc_calculate_forces','total_active_density_fine',&
            ierr)
       allocate(xc_pot_fine_emft(grid%ld1,grid%ld2,grid%max_slabs12, &
            pub_num_spins),stat=ierr)
       call utils_alloc_check('pseudo_nlcc_calculate_forces','xc_pot_fine_emft',ierr)
    end if
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pseudo_nlcc_calculate_forces','recip',ierr)

    ! Calculate total density
    factor = 0.5_DP * pub_spin_fac
    do is=1,pub_num_spins
       total_density_fine(:,:,:,is) = density_fine(:,:,:,is) &
            + core_density_fine*factor
    end do

    ! Calculate the exchange correlation potential
    call xc_energy_potential(total_density_fine,xc_energy,xc_pot_fine,grid,cell,0)

    ! jcap: If we are using embedded mean field theory, we need to
    ! call this twice more, in order to get the xc energy and
    ! potential for the active region
    if (loc_emft) then

       do is=1,pub_num_spins
          total_active_density_fine(:,:,:,is) = &
               active_density_fine(:,:,:,is) + core_density_fine*factor
       end do

       ! rc2013: calculate the EMFT correction
       call xc_emft_calculate(total_active_density_fine, xc_energy_emft, &
            xc_pot_fine_emft, grid, cell)

       ! jcap: add the contribution of the active subregion at the
       ! high level of theory
       xc_pot_fine = xc_pot_fine + xc_pot_fine_emft
    end if

    ! Average the spin channels in (:,:,:,1)
    if (pub_num_spins == 2) then
       xc_pot_fine(:,:,:,1) = factor*(xc_pot_fine(:,:,:,1) &
            + xc_pot_fine(:,:,:,2))
    end if

    ! Initialise
    nlcc_forces = 0.0_DP

    ! Fourier transform the xc potential to reciprocal space
    call fourier_apply_cell_forward(xc_pot_fine(:,:,:,1),recip,grid)

!$OMP PARALLEL NUM_THREADS(pub_threads_max) DEFAULT(none) &
!$OMP PRIVATE(ipt,i2,i3,gvec,g_length,species,coreden, &
!$OMP      iG_coreden_vxc,factor,islab23,atom,gdotR,eiGR,ierr) &
!$OMP SHARED (pub_my_proc_id,grid,p_species,elements,recip,par, &
!$OMP      pub_threads_max) &
!$OMP REDUCTION (+:nlcc_forces)

    allocate(iG_coreden_vxc(1:3,par%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_nlcc_calculate_forces','iG_coreden_vxc',ierr)

    ! For components g with symmetry points at -g
    factor=2.0_DP

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

       ! Loop over atoms to find n(G) for each
       do species=1,par%num_pspecies

          ! Check if we actually have a core charge for this species
          if (.not.p_species(species)%core_charge) cycle

          ! Interpolate value of core density at current g
          coreden = services_1d_interpolation( &
               p_species(species)%core_charge_recip, &
               p_species(species)%n_rad_pts, &
               g_length*p_species(species)%inv_g_spacing,0)

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
          if (.not.p_species(species)%core_charge) cycle

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
    call utils_dealloc_check('pseudo_nlcc_calculate_forces','iG_coreden_vxc',ierr)

!$OMP END PARALLEL

    ! Sum the result over all procs
    call comms_reduce('SUM',nlcc_forces)

    ! Deallocate
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('pseudo_nlcc_calculate_forces','recip',ierr)
    deallocate(xc_pot_fine,stat=ierr)
    call utils_dealloc_check('pseudo_nlcc_calculate_forces','xc_pot_fine',ierr)
    deallocate(total_density_fine,stat=ierr)
    call utils_dealloc_check('pseudo_nlcc_calculate_forces', &
         'total_density_fine',ierr)
    if (loc_emft) then
       deallocate(total_active_density_fine,stat=ierr)
       call utils_dealloc_check('pseudo_nlcc_calculate_forces', &
            'total_active_density_fine',ierr)
       deallocate(xc_pot_fine_emft,stat=ierr)
       call utils_dealloc_check('pseudo_nlcc_calculate_forces','xc_pot_fine_emft',ierr)
    end if

    ! Scale factor
    nlcc_forces = nlcc_forces / real(grid%n1*grid%n2*grid%n3,DP)

    ! Stop timer
    call timer_clock('pseudo_nlcc_calculate_forces',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nlcc_calculate_forces'

    return
  end subroutine pseudo_nlcc_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! SUBROUTINES FOR THE NON-LOCAL PSEUDOPOTENTIAL !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_dij(dij,p_species)

    !=======================================================================!
    ! This subroutine finds the Kleinman-Bylander denominators for the      !
    ! nonlocal pseudopotential projectors on this proc and puts 1/kb as the !
    ! diagonal elements of the matrix kb_denominators.                      !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine on 15 June 2009.                             !
    ! Edited for general parallel strategies by Robert Charlton, 15/12/2016.!
    !=======================================================================!

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_put_block, sparse_get_par
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dij
    type(PSEUDO_SPECIES), intent(in) :: p_species(:)! size(par%num_pspecies)

    ! Local Variables
    integer :: iat, loc_iat
    integer :: isp
    integer :: ishell, jshell
    integer :: iproj, jproj
    integer :: li, lj
    integer :: mi, mj
    integer :: ierr
    integer :: max_atom_proj
    real(kind=DP), allocatable :: dij0_sp(:,:,:)
    complex(kind=DP), allocatable :: dij0_sp_cmplx(:,:,:)
    logical :: iscmplx
    ! rc2013: parallel strategy
    type(PARAL_INFO), pointer :: par

    ! ndmh: initialisations
    iscmplx = dij%iscmplx
    max_atom_proj = maxval(p_species(:)%n_proj)

    ! rc2013: get parallel strategy from matrix structure, assume row_par=col_par
    ! pseudo_get_dij should only be used with diagonal matrix blocks, so this
    ! should always be true
    call sparse_get_par(par,dij,'R')
    call utils_assert(par%num_pspecies == size(p_species), 'Error in &
         &pseudo_get_dij: allocated parallel strategy is &
         &incompatible with p_species.')

    allocate(dij0_sp(max_atom_proj,max_atom_proj, &
         par%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_get_dij','dij0_sp',ierr)
    dij0_sp(:,:,:) = 0.0_DP

    ! Loop over species, setting up array of dij0 terms
    do isp=1,par%num_pspecies
       if (any(par%elements_on_proc(:)%pspecies_number==isp)) then
          iproj = 0
          do ishell=1,p_species(isp)%n_shells
             li = p_species(isp)%ang_mom(ishell)
             do mi=-li,li
                iproj = iproj + 1
                jproj = 0
                do jshell=1,p_species(isp)%n_shells
                   lj = p_species(isp)%ang_mom(jshell)
                   do mj=-lj,lj
                      jproj = jproj + 1
                      if (mj==mi) dij0_sp(iproj,jproj,isp) = &
                           p_species(isp)%D0(ishell,jshell)
                   end do
                end do
             end do
          end do
       end if
    end do

    if (iscmplx) then
       allocate(dij0_sp_cmplx(max_atom_proj,max_atom_proj, &
            par%num_pspecies),stat=ierr)
       call utils_alloc_check('pseudo_get_dij','dij0_sp_cmplx',ierr)
       dij0_sp_cmplx(:,:,:) = cmplx(dij0_sp(:,:,:),kind=DP)
    end if

    ! Loop over atoms
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = loc_iat + par%first_atom_on_proc(pub_my_proc_id) - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number
       ! Put block of dij0 into SPAM3 matrix if there is one
       if (p_species(isp)%n_shells>0) then
          if (iscmplx) then
             call sparse_put_block(dij0_sp_cmplx(:,:,isp),dij,iat,iat)
          else
             call sparse_put_block(dij0_sp(:,:,isp),dij,iat,iat)
          end if
       end if
    end do

    if (iscmplx) then
       deallocate(dij0_sp_cmplx,stat=ierr)
       call utils_dealloc_check('pseudo_get_dij','dij0_sp_cmplx',ierr)
    end if

    deallocate(dij0_sp,stat=ierr)
    call utils_dealloc_check('pseudo_get_dij','dij0_sp',ierr)

  end subroutine pseudo_get_dij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_nonlocal_mat(nonlocpot, sp_overlap, dij)

    !==================================================================!
    ! This subroutine calculates the sparse nonlocal potential matrix  !
    ! from the ngwf-projector overlap matrix.                          !
    !------------------------------------------------------------------!
    ! Original version written by Arash A. Mostofi on 12/01/2004.      !
    ! Modified by Chris-Kriton Skylaris on 25/1/2004 to work with      !
    ! the pseudo_species type array.                                   !
    ! Modified by Chris-Kriton Skylaris on 30/1/2004 to use the        !
    ! sparsity pattern of nonlocpot for determining which elements     !
    ! to calculate and store.                                          !
    ! Modified by Peter Haynes to use parallel SPAM 2, July 2006       !
    ! Modified by Nick Hine to use sparse_symmetric_expand to even out !
    ! load between processes                                           !
    ! Rewritten using SPAM3 matrices by Nicholas Hine, June 2009       !
    ! Modified by Nicholas Hine on 03/11/2009 to no longer need        !
    ! elements array, since pseudo_get_dij no longer does              !
    ! Modified by Laura Ratcliff to allow for calculation of a         !
    ! non-local matrix between 2 different sets of NGWFs, October 2010 !
    ! Modified for embedding by Robert Charlton, 23/09/2017.           !
    ! Modified to accept dij as input so that multiple lists of        !
    ! nonlocal projectors can be used, 11/10/2017.                     !
    !==================================================================!

    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
        sparse_embed_destroy,sparse_embed_product, &
        sparse_embed_transpose
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: nonlocpot
    type(SPAM3_EMBED), intent(in)    :: sp_overlap
    type(SPAM3_EMBED), intent(in)    :: dij

    ! Local Variables
    type(SPAM3_EMBED) :: sp_overlap_dij, ps_overlap

    call timer_clock('pseudopotentials_nonlocal_mat',1)

    ! Create the matrix structures (set as real or complex depending on
    ! whether sp_overlap is real or complex)
    call sparse_embed_create(sp_overlap_dij,sp_overlap)!,iscmplx=sp_overlap%iscmplx)

    ! ndmh: Create ps_overlap (modifying structure code if required)
    ! rc2013: Create ps_overlap from the transpose of ps_overlap
    call sparse_embed_create(ps_overlap, sp_overlap, trans=.true.)

    ! Calculate the matrix <NGWF_a|Proj_i> * D_ij
    call sparse_embed_product(sp_overlap_dij,sp_overlap,dij)

    ! Create the transpose of the NGWF-projector overlap matrix
    call sparse_embed_transpose(ps_overlap,sp_overlap)

    ! Calculate the matrix \sum_i (<NGWF_a|Proj_i> D_ij <Proj_j|NGWF_b> )
    call sparse_embed_product(nonlocpot,sp_overlap_dij,ps_overlap)

    ! Clean up temporary matrices
    call sparse_embed_destroy(ps_overlap)
    call sparse_embed_destroy(sp_overlap_dij)

    call timer_clock('pseudopotentials_nonlocal_mat',2)

  end subroutine pseudopotentials_nonlocal_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pseudo_nonlocal_commutator_mat(nonlocpot_com, &
       proj_basis, proj_set, ngwf_basis, ngwfs_on_grid, sp_overlap, &
       cell, fftbox, dij, delta_in)

    !==================================================================!
    ! This subroutine calculates the commutator between the nonlocal   !
    ! potential and the position operator for the 3 Cartesian          !
    ! directions.                                                      !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff February 2011.                         !
    ! Modified by Nicholas Hine in April 2011 to not use NGWF_REP and  !
    ! to simplify and cleanup existing code, and for changes to        !
    ! projector initialisation routines.                               !
    ! Modified by Laura Ratcliff Dec 2011 to fix bug involving complex !
    ! matrices and phase factors.                                      !
    ! Modified for embedding by Joseph Prentice, July 2018             !
    !==================================================================!

    use datatypes, only: FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: operator(*)
    use ion, only: ELEMENT
    use projectors, only: PROJECTOR_SET, projectors_commutator
    use simulation_cell, only: CELL_INFO
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, &
         sparse_embed_transpose, sparse_embed_transpose_structure
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: nonlocpot_com(3)
    type(PROJECTOR_SET), intent(inout) :: proj_set(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid(:)
    type(SPAM3_EMBED), intent(in) :: sp_overlap
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(SPAM3_EMBED), intent(in) :: dij
    real(kind=DP), intent(in) :: delta_in ! finite difference shift

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nonlocal_commutator_mat'

    call timer_clock('pseudo_nonlocal_commutator_mat',1)

    ! Matrix of Kleinman-Bylander denominators precalculated and
    ! passed in as dij

    call projectors_commutator(nonlocpot_com, proj_basis, &
         ngwf_basis, ngwfs_on_grid, sp_overlap, proj_set, &
         cell, fftbox, delta_in, dij)

    call timer_clock('pseudo_nonlocal_commutator_mat',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nonlocal_commutator_mat'

  end subroutine pseudo_nonlocal_commutator_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pseudo_nonlocal_com_mat_fd(nonlocpot_com, &
       proj_basis, proj_set, ngwf_basis, ngwfs_on_grid, sp_overlap, &
       cell, fftbox, dij, nonlocpot, delta)

    !==================================================================!
    ! This subroutine calculates the commutator between the nonlocal   !
    ! potential and the position operator for the 3 Cartesian          !
    ! directions.                                                      !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff February 2011.                         !
    ! Modified by Nicholas Hine in April 2011 to not use NGWF_REP      !
    ! Modified by Laura Ratcliff Dec 2011 to fix bug involving complex !
    ! matrices and phase factors and changed from centred difference   !
    ! to forward difference.                                           !
    ! Modified for embedding by Joseph Prentice, July 2018             !
    !==================================================================!

    use datatypes, only: FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*)
    use ion, only: ELEMENT
    use projectors, only: PROJECTOR_SET, projectors_func_ovlp_box
    use simulation_cell, only: CELL_INFO
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_transpose, &
         sparse_embed_copy, sparse_embed_transpose_structure, &
         sparse_embed_axpy, sparse_embed_scale
    use sparse, only: sparse_axpy, sparse_scale, sparse_copy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: nonlocpot_com(3)
    type(SPAM3_EMBED), intent(in) :: nonlocpot
    type(PROJECTOR_SET), intent(inout) :: proj_set(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid(:)
    type(SPAM3_EMBED), intent(in) :: sp_overlap
    real(kind=DP), intent(in) :: delta ! finite difference shift
    type(SPAM3_EMBED), intent(in) :: dij
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    type(SPAM3_EMBED) :: sp_overlap_dij, ps_overlap, dij_ps_overlap
    type(SPAM3_EMBED) :: sDGp_overlap, DGps_overlap, sp_overlap_cmplx
    type(SPAM3_EMBED) :: nonlocpot_com2
    type(SPAM3_EMBED) :: ps_overlap_rc, sDGp_overlap_rc
    type(POINT) :: q_vec, q_cart(3)
    real(kind=DP) :: inv_delta
    integer :: cart, isub, jsub
    logical :: ovlp_iscmplx

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nonlocal_com_mat_fd'

    call timer_clock('pseudo_nonlocal_com_mat_fd',1)

    ! NB not allowing for sp_overlap or ps_overlap to be at diff k-points
    inv_delta = 1.0_dp / delta

    q_vec%x = delta
    q_vec%y = delta
    q_vec%z = delta

    do cart=1,3
      q_cart(cart)%x = 0.0_dp
      q_cart(cart)%y = 0.0_dp
      q_cart(cart)%z = 0.0_dp
    end do

    q_cart(1)%x = q_vec%x
    q_cart(2)%y = q_vec%y
    q_cart(3)%z = q_vec%z

    ovlp_iscmplx = sp_overlap%iscmplx

    call sparse_embed_create(nonlocpot_com2,nonlocpot_com(1),iscmplx=.true.)

    call sparse_embed_create(sp_overlap_dij,sp_overlap,iscmplx=.true.)
    call sparse_embed_create(sDGp_overlap,sp_overlap,iscmplx=.true.)

    ! ndmh: Create ps_overlap
    call sparse_embed_transpose_structure(ps_overlap%structure,sp_overlap)
    ps_overlap_rc%structure = ps_overlap%structure
    dij_ps_overlap%structure = ps_overlap%structure
    DGps_overlap%structure = ps_overlap%structure

    call sparse_embed_create(dij_ps_overlap,iscmplx=ovlp_iscmplx)
    call sparse_embed_create(DGps_overlap,iscmplx=.true.)
    call sparse_embed_create(ps_overlap,iscmplx=ovlp_iscmplx)

    call sparse_embed_transpose(ps_overlap,sp_overlap)

    ! The matrix of Kleinman-Bylander denominators is now
    ! precalculated and passed in as dij

    call sparse_embed_product(sp_overlap_dij,sp_overlap,dij)
    call sparse_embed_product(dij_ps_overlap,dij,ps_overlap)

    call sparse_embed_destroy(ps_overlap)

    ! Loop over Cartesian directions
    do cart=1,3

       q_vec = q_cart(cart)

       ! Calculate <phi|Del_G(proj(G))> overlap matrix for ngwf_basis
       ! jcap: loop over regions
       do jsub=1,size(proj_basis(:))
          do isub=1,size(ngwf_basis(:))
             call projectors_func_ovlp_box(sDGp_overlap%m(isub,jsub), &
                  ngwfs_on_grid(isub), ngwf_basis(isub), &
                  proj_basis(jsub), proj_set(jsub), &
                  fftbox, cell, q_vec)
          end do
       end do

       ! Transpose <phi|Del_G(proj(G))> to get <Del_G(proj(G))|phi>
       call sparse_embed_transpose(DGps_overlap,sDGp_overlap)

       ! Calculate the matrix \sum_ij <phi_a|Del_G(p_i(G))>D_ij<p_j|phi_b>
       call sparse_embed_product(nonlocpot_com2,sDGp_overlap,dij_ps_overlap, &
            allow_mix_types=(.not.ovlp_iscmplx))

       ! Calculate the matrix \sum_ij <phi_a|p_i>D_ij<Del_G(p_j(G))|phi_b>
       call sparse_embed_product(nonlocpot_com(cart),sp_overlap_dij,DGps_overlap, &
            allow_mix_types=(.not.ovlp_iscmplx))

       ! Sum the two
       call sparse_embed_axpy(nonlocpot_com(cart),nonlocpot_com2,(1.0_DP,0.0_DP))

       ! Subtract 2 * nonlocpot
       call sparse_embed_axpy(nonlocpot_com(cart),nonlocpot,(-2.0_DP,0.0_DP))

       ! Multiply by -i/|q|
       call sparse_embed_scale(nonlocpot_com(cart),cmplx(0.0,-inv_delta,kind=DP))

    end do

    ! Clean up temporary matrices
    call sparse_embed_destroy(DGps_overlap)
    call sparse_embed_destroy(dij_ps_overlap)
    call sparse_embed_destroy(sDGp_overlap)
    call sparse_embed_destroy(sp_overlap_dij)
    call sparse_embed_destroy(nonlocpot_com2)

    call timer_clock('pseudo_nonlocal_com_mat_fd',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nonlocal_com_mat_fd'

  end subroutine pseudo_nonlocal_com_mat_fd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pseudo_nonlocal_com_mat_direct(nonlocpot_com, proj_set, &
       proj_basis, ngwf_basis, ngwfs_on_grid, sp_overlap, dij, cell, fftbox)

    !==================================================================!
    ! This subroutine calculates the commutator between the nonlocal   !
    ! potential and the position operator for the 3 Cartesian          !
    ! directions.                                                      !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff February 2011.                         !
    ! Modified for embedding by Joseph Prentice, July 2018             !
    ! (on the off-chance this will be used again)                      !
    !==================================================================!

    use datatypes, only: FUNCTIONS
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*)
    use ion, only: ELEMENT
    use projectors, only: PROJECTOR_SET, projectors_func_pos_ovlp_box
    use simulation_cell, only: CELL_INFO
    use sparse_embed, only: SPAM3_EMBED, sparse_embed_create, sparse_embed_destroy, &
         sparse_embed_product, sparse_embed_transpose, sparse_embed_axpy, &
         sparse_embed_scale, sparse_embed_copy, sparse_embed_transpose_structure
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_destroy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3_EMBED), intent(inout) :: nonlocpot_com(3)
    type(PROJECTOR_SET), intent(inout) :: proj_set(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid(:)
    type(SPAM3_EMBED), intent(in) :: sp_overlap
    type(SPAM3_EMBED), intent(in) :: dij
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox

    ! Local Variables
    type(SPAM3_EMBED) :: sp_overlap_KB, ps_overlap, KB_ps_overlap
    type(SPAM3_EMBED) :: sDGp_overlap(3), DGps_overlap
    type(SPAM3_EMBED) :: nonlocpot_com2
    type(SPAM3) :: temp_spam(3)

    integer :: cart, isub, jsub

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nonlocal_com_mat_direct'

    call timer_clock('pseudo_nonlocal_com_mat_direct',1)

    ! Assuming all matrices are real

    call sparse_embed_create(nonlocpot_com2,nonlocpot_com(1),iscmplx=.false.)
    call sparse_embed_create(sp_overlap_KB,sp_overlap)
    do cart=1,3
       call sparse_embed_create(sDGp_overlap(cart),sp_overlap)
    end do

    ! ndmh: Create ps_overlap (modifying structure code if required)
    call sparse_embed_transpose_structure(ps_overlap%structure,sp_overlap)
    KB_ps_overlap%structure = ps_overlap%structure
    DGps_overlap%structure = ps_overlap%structure

    call sparse_embed_create(KB_ps_overlap, iscmplx=.false.)
    call sparse_embed_create(DGps_overlap, iscmplx=.false.)
    call sparse_embed_create(ps_overlap, iscmplx=.false.)

    ! Create the transpose of the NGWF-projector overlap matrix
    call sparse_embed_transpose(ps_overlap,sp_overlap)

    ! The matrix of Kleinman-Bylander denominators is precalculated
    ! and passed in as dij

    call sparse_embed_product(sp_overlap_KB,sp_overlap,dij)
    call sparse_embed_product(KB_ps_overlap,dij,ps_overlap)

    call sparse_embed_destroy(ps_overlap)

    ! Calculate <phi|r|proj> overlap matrix
    ! jcap: loop over regions
    do isub=1,size(ngwf_basis(:))
       do jsub=1,size(proj_basis(:))

          do cart=1,3
             call sparse_create(temp_spam(cart),sDGp_overlap(cart)%m(isub,jsub))
             call sparse_copy(temp_spam(cart),sDGp_overlap(cart)%m(isub,jsub))
          end do

          call projectors_func_pos_ovlp_box(temp_spam, &
               ngwfs_on_grid(isub),ngwf_basis(isub),&
               proj_basis(jsub),proj_set(jsub),fftbox,cell)

          do cart=1,3
             call sparse_copy(sDGp_overlap(cart)%m(isub,jsub),temp_spam(cart))
             call sparse_destroy(temp_spam(cart))
          end do

       end do
    end do

    ! Loop over Cartesian directions
    do cart=1,3

       ! Transpose <phi|r|proj> to get <proj|r|phi> for ngwf_basis
       call sparse_embed_transpose(DGps_overlap,sDGp_overlap(cart))

       ! Calculate the matrix \sum_i (<NGWF_a|r|Proj_i><Proj_i|NGWF_b>/D_i)
       call sparse_embed_product(nonlocpot_com2,sDGp_overlap(cart),KB_ps_overlap)

       ! Calculate the matrix \sum_i (<NGWF_a|Proj_i><Proj_i|r|NGWF_b>/D_i)
       call sparse_embed_product(nonlocpot_com(cart),sp_overlap_KB,DGps_overlap)

       ! Sum the two
       call sparse_embed_axpy(nonlocpot_com(cart),nonlocpot_com2,-1.0_dp)

    end do

    ! Clean up temporary matrices
    call sparse_embed_destroy(DGps_overlap)
    call sparse_embed_destroy(KB_ps_overlap)
    do cart=3,1,-1
       call sparse_embed_destroy(sDGp_overlap(cart))
    end do
    call sparse_embed_destroy(sp_overlap_KB)
    call sparse_embed_destroy(nonlocpot_com2)

    call timer_clock('pseudo_nonlocal_com_mat_direct',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nonlocal_com_mat_direct'

  end subroutine pseudo_nonlocal_com_mat_direct


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_aug_Q_matrix(aug_Q,p_species)

    !==================================================================!
    ! This subroutine returns the Q matrix for all the USP atoms in    !
    ! the simulation cell.                                             !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  proj_overlap (inout) : The SPAM3 block-diagonal overlap matrix  !
    !  of the partial waves of each atom                               !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 11/02/11.           !
    ! Modified to remove pub_par by Joseph Prentice, June 2018         !
    !==================================================================!

    use comms, only: pub_my_proc_id
    use parallel_strategy, only: PARAL_INFO
    use sparse, only: SPAM3, sparse_put_block, sparse_get_par
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: aug_Q
    type(PSEUDO_SPECIES), intent(in) :: p_species(:)

    ! Local Variables
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: ishell,jshell
    integer :: iproj_at,jproj_at
    integer :: li,lj,mi,mj
    integer :: ierr
    integer :: max_proj
    real(kind=DP),allocatable :: aug_Q_block(:,:)
    ! agrecocmplx
    complex(kind=DP), allocatable :: aug_Q_block_cmplx(:,:)
    logical :: iscmplx
    type(PARAL_INFO), pointer :: par

    iscmplx = aug_Q%iscmplx

    ! rc2013: get the parallel strategy for this diagonal matrix
    call sparse_get_par(par, aug_Q)
    call utils_assert(par%num_pspecies == size(p_species), 'Error in &
         &pseudo_aug_Q_matrix: allocated parallel strategy is &
         &incompatible with p_species.')

    max_proj = maxval(p_species(:)%n_proj)
    allocate(aug_Q_block(max_proj,max_proj),stat=ierr)
    call utils_alloc_check('pseudo_aug_Q_matrix','aug_Q_block',ierr)

    if (iscmplx) then
       allocate(aug_Q_block_cmplx(max_proj,max_proj),stat=ierr)
       call utils_alloc_check('pseudo_aug_Q_matrix','aug_Q_block_cmplx',ierr)
    end if

    ! Loop over atoms on this proc
    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1
       isp = par%elements_on_proc(loc_iat)%pspecies_number

       ! Double loop over partial waves i,j
       iproj_at = 0
       do ishell=1,p_species(isp)%n_shells
          li = p_species(isp)%ang_mom(ishell)
          do mi=-li,li
             iproj_at = iproj_at + 1
             jproj_at = 0
             do jshell=1,p_species(isp)%n_shells
                lj = p_species(isp)%ang_mom(jshell)
                do mj=-lj,lj
                   jproj_at = jproj_at + 1

                   ! Find Q_ij for this projector
                   aug_Q_block(iproj_at,jproj_at) = &
                        p_species(isp)%aug_q(ishell,jshell)
                enddo
             end do
          end do
       end do

       ! Put this atom's block into SPAM3 matrix
       ! agrecocmplx
       if (iscmplx) then
          aug_Q_block_cmplx(:,:) = cmplx(aug_Q_block(:,:),kind=DP)
          call sparse_put_block(aug_Q_block_cmplx,aug_Q,iat,iat)
       else
          call sparse_put_block(aug_Q_block,aug_Q,iat,iat)
       end if
    end do

    ! agrecocmplx
    if (iscmplx) then
       deallocate(aug_Q_block_cmplx, stat=ierr)
       call utils_dealloc_check('pseudo_aug_Q_matrix','aug_Q_block_cmplx',ierr)
    end if

    deallocate(aug_Q_block,stat=ierr)
    call utils_dealloc_check('pseudo_aug_Q_matrix','aug_Q_block',ierr)

  end subroutine pseudo_aug_Q_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_atom_aug_den(atom_nhat,atom_grad_nhat,atom_aug_func, &
       atom_aug_func_grad,total_nhat,total_nhat_targ,rho_ij_block, &
       isp,atom_centre,grid,cell,p_species,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use fourier, only: fourier_apply_box
    use geometry, only: POINT
    use rundat, only: pub_num_spins, pub_cmplx_ngwfs
    use simulation_cell, only: CELL_INFO
    use xc, only: pub_xc_gradient_corrected
    ! agrecocmplx
    use utils, only: utils_assert

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    real(kind=DP), intent(out) :: atom_aug_func(box_n1,box_n2,box_n3)
    real(kind=DP), intent(out) :: atom_aug_func_grad(box_n1,box_n2,box_n3,3)
    real(kind=DP), intent(out) :: atom_nhat(box_n1,box_n2,box_n3, &
         pub_num_spins)
    real(kind=DP), intent(out) :: atom_grad_nhat(box_n1,box_n2,box_n3, &
         pub_num_spins,3)
    real(kind=DP), intent(out) :: total_nhat(max_spins)
    real(kind=DP), intent(out) :: total_nhat_targ(max_spins)
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(PSEUDO_SPECIES), intent(in) :: p_species(:) !size(par%num_pspecies)
    type(POINT), intent(in) :: atom_centre
    ! agrecocmplx: need to be complex with complex NGWFs
    real(kind=DP), intent(in) :: rho_ij_block(:,:,:)
    integer, intent(in) :: isp
    !integer, intent(in) :: iat
    integer, intent(in) :: box_start1,box_start2,box_start3

    ! Local Variables
    integer :: ishell,jshell
    integer :: iproj,jproj
    integer :: li,mi,lj,mj
    integer :: is
    integer :: cart
    ! need to be complex with complex NGWFs
    real(kind=DP) :: Qij(max_spins)
    real(kind=DP) :: total_aug_func

    ! agrecocmplx
    call utils_assert(pub_cmplx_ngwfs.eqv..false., &
         'Subroutine pseudo_atom_aug_den not ready yet for &
         & complex NGWFs.')

    ! Initialisation
    atom_nhat = 0.0_DP
    atom_grad_nhat = 0.0_DP
    total_nhat = 0.0_DP
    total_nhat_targ = 0.0_DP

    ! Loop over projectors i,j
    iproj = 0
    do ishell=1,p_species(isp)%n_shells
       li = p_species(isp)%ang_mom(ishell)
       do mi=-li,li
          iproj = iproj + 1
          jproj = 0
          do jshell=1,p_species(isp)%n_shells
             lj = p_species(isp)%ang_mom(jshell)
             do mj=-lj,lj
                jproj = jproj + 1

                Qij(:) = 0.0_DP

                ! Find contribution to total target augmentation charge
                ! for this i,j pair
                do is=1,pub_num_spins
                   Qij(is) = Qij(is) + rho_ij_block(iproj,jproj,is)
                   !write(*,*) 'rho: ',isp,iat,iproj,jproj,rho_ij_block(iproj,jproj,is)
                end do

                !if (any(abs(Qij(:))>1e-20_DP)) then

                   ! Get shape function g_L(r) multiplied by spherical harmonic
                   ! S_LM(r) for this L,M pair
                   atom_aug_func = 0.0_DP
                   call pseudo_get_aug_func(atom_aug_func, &
                        isp,p_species,ishell,jshell,mi,mj, &
                        box_start1,box_start2,box_start3, &
                        grid,cell,box_n1,box_n2,box_n3,atom_centre)

                   total_aug_func = sum(atom_aug_func(:,:,:))*grid%weight

                   do is=1,pub_num_spins
                      atom_nhat(:,:,:,is) = atom_nhat(:,:,:,is) &
                           + Qij(is)*atom_aug_func(:,:,:)
                      total_nhat(is) = total_nhat(is) + total_aug_func*Qij(is)
                   end do

                   ! Find gradient of \hat{n}(r) if required
                   if (pub_xc_gradient_corrected) then
                      ! Get gradient of (shape function g_L(r) multiplied by
                      ! spherical harmonic S_LM(r) for this L,M pair)
                      atom_aug_func_grad = 0.0_DP
                      call pseudo_get_aug_func_grad(atom_aug_func_grad, &
                           isp,p_species,ishell,jshell,mi,mj, &
                           box_start1,box_start2,box_start3, &
                           grid,cell,box_n1,box_n2,box_n3,atom_centre)
                      do cart=1,3
                         do is=1,pub_num_spins
                            atom_grad_nhat(:,:,:,is,cart) = &
                                 atom_grad_nhat(:,:,:,is,cart) + &
                                 Qij(is)*atom_aug_func_grad(:,:,:,cart)
                         end do
                      end do
                   end if
                !end if

             end do
          end do
       end do
    end do

  end subroutine pseudo_atom_aug_den


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_atom_aug_integrals(locpot_box,atom_aug_func, &
       dij_at,num_spins,isp,atom_centre,grid,cell,p_species,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use geometry, only: POINT
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    integer, intent(in) :: num_spins
    real(kind=DP), intent(out) :: atom_aug_func(box_n1,box_n2,box_n3)
    real(kind=DP), intent(in) :: locpot_box(box_n1,box_n2,box_n3,num_spins)
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(POINT), intent(in) :: atom_centre
    real(kind=DP), intent(inout) :: dij_at(:,:,:)
    integer, intent(in) :: isp
    type(PSEUDO_SPECIES), intent(in) :: p_species(:) ! size(par%num_pspecies)
    integer, intent(in) :: box_start1,box_start2,box_start3

    ! Local Variables
    integer :: iproj, jproj
    integer :: ishell, jshell
    integer :: li, lj
    integer :: mi, mj
    integer :: is
    real(kind=DP) :: locpot_Qij_product(max_spins)

    ! Loop over projectors i,j
    iproj = 0
    do ishell=1,p_species(isp)%n_shells
       li = p_species(isp)%ang_mom(ishell)
       do mi=-li,li
          iproj = iproj + 1
          jproj = 0
          do jshell=1,p_species(isp)%n_shells
             lj = p_species(isp)%ang_mom(jshell)
             do mj=-lj,lj
                jproj = jproj + 1

                ! Get the augmentation function Qij(r) corresponding
                ! to this pair of shells i,j and this mi,mj
                atom_aug_func = 0.0_DP
                call pseudo_get_aug_func(atom_aug_func, &
                     isp,p_species,ishell,jshell,mi,mj, &
                     box_start1,box_start2,box_start3, &
                     grid,cell,box_n1,box_n2,box_n3,atom_centre)

                ! Calculate integral \int veff(r) Qij(r) dr
                do is=1,num_spins
                   locpot_Qij_product(is) = sum(locpot_box(:,:,:,is) &
                        * atom_aug_func(:,:,:)) * grid%weight

                   dij_at(iproj,jproj,is) = dij_at(iproj,jproj,is) + &
                        locpot_Qij_product(is)

                end do

             end do
          end do
       end do
    end do

  end subroutine pseudo_atom_aug_integrals


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_atom_aug_force(nhat_force,locpot_box,atom_grad_aug_func, &
       rhoij_at,num_spins,isp,atom_centre,grid,cell,p_species,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use geometry, only: POINT
    use simulation_cell, only: CELL_INFO

    implicit none

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    integer, intent(in) :: num_spins
    real(kind=DP), intent(out) :: atom_grad_aug_func(box_n1,box_n2,box_n3,3)
    real(kind=DP),intent(inout) :: nhat_force(3)
    real(kind=DP), intent(in) :: locpot_box(box_n1,box_n2,box_n3,num_spins)
    type(GRID_INFO), intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(POINT), intent(in) :: atom_centre
    real(kind=DP), intent(in) :: rhoij_at(:,:,:)
    integer, intent(in) :: isp
    type(PSEUDO_SPECIES), intent(in) :: p_species(:)! size(par%num_pspecies)
    integer, intent(in) :: box_start1,box_start2,box_start3

    ! Local Variables
    integer :: ishell,jshell
    integer :: iproj,jproj
    integer :: li,mi,lj,mj
    integer :: is
    integer :: cart
    real(kind=DP) :: locpot_Qij_grad_product(3,max_spins)

    ! Loop over projectors i,j
    iproj = 0
    do ishell=1,p_species(isp)%n_shells
       li = p_species(isp)%ang_mom(ishell)
       do mi=-li,li
          iproj = iproj + 1
          jproj = 0
          do jshell=1,p_species(isp)%n_shells
             lj = p_species(isp)%ang_mom(jshell)
             do mj=-lj,lj
                jproj = jproj + 1

                ! Construct gradient d/dR_I (Qij(r)) of augmentation
                ! function for this i,j pair
                atom_grad_aug_func = 0.0_DP
                call pseudo_get_aug_func_grad(atom_grad_aug_func, &
                     isp,p_species,ishell,jshell,mi,mj,box_start1,box_start2, &
                     box_start3,grid,cell,box_n1,box_n2,box_n3, &
                     atom_centre)

                ! Loop over spins and loop over cartesian directions
                do is=1,num_spins
                   do cart=1,3
                      ! Calculate integral \int veff(r) d/dR(Qij(r)) dr
                      locpot_Qij_grad_product(cart,is) = &
                           sum(locpot_box(:,:,:,is) &
                           * atom_grad_aug_func(:,:,:,cart)) * grid%weight

                   end do
                end do

                do is=1,num_spins
                   nhat_force(:) = nhat_force(:) &
                        + (locpot_Qij_grad_product(:,is) &
                        * rhoij_at(iproj,jproj,is))
                end do  ! is

             end do
          end do
       end do
    end do

  end subroutine pseudo_atom_aug_force


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_aug_func(aug_func,isp,p_species,ishell,jshell,mi,mj, &
       box_start1,box_start2,box_start3,grid,cell,box_n1,box_n2,box_n3, &
       atom_origin)

    !==================================================================!
    ! This subroutine calculates the function                          !
    !    Qnm(r) = \sum_LM c_LM^nm Y_LM(\hat{r}) Q_nm^rad(r)            !
    ! according to Eq 26-30 of Laasonen et al, PRB 47 10142 (1993)     !
    ! which is required to form the compensation charge \hat{n}.       !
    ! It does so on the simulation cell fine grid in a box centered on !
    ! the atom.                                                        !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  box_start1 (in) :                                               !
    !  box_start2 (in) :                                               !
    !  box_start3 (in) :                                               !
    !  grid (in) :                                                     !
    !  box_n1 (in) :                                                   !
    !  box_n2 (in) :                                                   !
    !  box_n3 (in) :                                                   !
    !  atom_origin (in) :                                              !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 15/02/11.           !
    ! NOT FINISHED YET - NON-FUNCTIONAL!                               !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: PI
    use gaunt_coeff, only: realgaunt
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use services, only: services_locate_interp, services_linear_interpolation
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_real_sph_harm

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    real(kind=DP),intent(out) :: aug_func(box_n1,box_n2,box_n3)
    integer, intent(in) :: isp
    type(PSEUDO_SPECIES), intent(in) :: p_species(:) ! size(par%num_pspecies)
    integer,intent(in) :: ishell
    integer,intent(in) :: jshell
    integer,intent(in) :: mi,mj
    integer,intent(in) :: box_start1
    integer,intent(in) :: box_start2
    integer,intent(in) :: box_start3
    type(GRID_INFO),intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(POINT),intent(in) :: atom_origin

    ! Local Variables
    integer :: li, lj
    integer :: i1,i2,i3
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    type(POINT) :: r_cell
    type(POINT) :: r_sphere
    real(kind=DP) :: rmag
    real(kind=DP) :: slmval
    integer :: ipt,jcoef, L, M
    real(kind=DP) :: QLnm(p_species(isp)%kkbeta,0:p_species(isp)%qf_lmax)
    real(kind=DP) :: rr, QLnm_rad

    ! Find vector to origin of box
    box_origin = (box_start1-1)*grid%da1 + &
                 (box_start2-1)*grid%da2 + &
                 (box_start3-1)*grid%da3

    ! ndmh: vector from origin of box to centre of atom
    box_to_atom = atom_origin - box_origin

    li = p_species(isp)%ang_mom(ishell)
    lj = p_species(isp)%ang_mom(jshell)

    ! ndmh: check components of this vector in terms of lattice vectors
    ! ndmh: are all positive. If not, the box has been looped back into the
    ! ndmh: cell, so shift the origin back accordingly
    if ((box_to_atom.DOT.cell%b1) < 0.0_DP) then
       box_origin = box_origin - cell%a1
    end if
    if ((box_to_atom.DOT.cell%b2) < 0.0_DP) then
       box_origin = box_origin - cell%a2
    end if
    if ((box_to_atom.DOT.cell%b3) < 0.0_DP) then
       box_origin = box_origin - cell%a3
    end if

    ! Generate QLnm(ipt,L)
    QLnm(:,:) = 0.0_DP
    do L=0,p_species(isp)%qf_lmax
       do ipt=1,p_species(isp)%kkbeta
          if(p_species(isp)%rlog(ipt) < p_species(isp)%rinner(L)) then
             rr = p_species(isp)%rlog(ipt)**2.0_DP
             QLnm(ipt,L) = p_species(isp)%qfcoef(1,L,ishell,jshell)
             do jcoef=2,p_species(isp)%nqfcoef
                QLnm(ipt,L) = QLnm(ipt,L) + &
                     p_species(isp)%qfcoef(jcoef,L,ishell,jshell)*rr**real(jcoef-1,dp)
             end do
             QLnm(ipt,L) = QLnm(ipt,L)*p_species(isp)%rlog(ipt)**real(L+2,dp)
          else
             QLnm(ipt,L) = p_species(isp)%qfunc(ipt,ishell,jshell)
          end if
       end do
       !write(62+isp,*) '#',isp,L,ishell,jshell,p_species(isp)%qfcoef(:,L,ishell,jshell)
       !do ipt=1,p_species(isp)%kkbeta
       !   write(62+isp,*) L,ishell,jshell,p_species(isp)%rlog(ipt),QLnm(ipt,L)
       !end do
   end do

    do i3=1,box_n3
       do i2=1,box_n2
          do i1=1,box_n1

             r_cell = box_origin + (i1-1)*grid%da1 &
                                 + (i2-1)*grid%da2 &
                                 + (i3-1)*grid%da3

             r_sphere = r_cell - atom_origin
             rmag = geometry_magnitude(r_sphere)

             if (rmag>p_species(isp)%core_radius(ishell)) then
                aug_func(i1,i2,i3) = 0.0_DP
                cycle
             end if

             do L=0,p_species(isp)%qf_lmax

                ipt = services_locate_interp(rmag,p_species(isp)%rlog,p_species(isp)%kkbeta)
                QLnm_rad = services_linear_interpolation(rmag,QLnm(ipt,L), &
                     QLnm(ipt+1,L),p_species(isp)%rlog(ipt),p_species(isp)%rlog(ipt+1))

                do M=-L,L

                   slmval = sw_real_sph_harm(r_sphere%x, &
                        r_sphere%y,r_sphere%z,rmag,L,M)

                   aug_func(i1,i2,i3) = aug_func(i1,i2,i3) + &
                        slmval*realgaunt(L,M,li,mi,lj,mj)*QLnm_rad*sqrt(4.0_DP*PI)
                end do
             end do

          end do
       end do
    end do

  end subroutine pseudo_get_aug_func


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_aug_func_grad(aug_func_grad,isp,p_species, &
       ishell,jshell,mi,mj,box_start1,box_start2,box_start3,grid,cell, &
       box_n1,box_n2,box_n3,atom_origin)

    !==================================================================!
    ! This subroutine calculates the function g_L(r)*S_LM(\hat{r})     !
    ! which is required to form the compensation charge \hat{n}.       !
    ! It does so on the simulation cell fine grid in a box centered on !
    ! the atom.                                                        !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  box_start1 (in) :                                               !
    !  box_start2 (in) :                                               !
    !  box_start3 (in) :                                               !
    !  grid (in) :                                                     !
    !  box_n1 (in) :                                                   !
    !  box_n2 (in) :                                                   !
    !  box_n3 (in) :                                                   !
    !  atom_origin (in) :                                              !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 15/02/11.           !
    ! NOT FINISHED YET - NON-FUNCTIONAL!                               !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use simulation_cell, only: CELL_INFO
    use spherical_wave, only: sw_real_sph_harm, sw_grad_real_sph_harm
    use utils, only: utils_use_var, utils_abort

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    real(kind=DP),intent(out) :: aug_func_grad(box_n1,box_n2,box_n3,3)
    integer, intent(in) :: isp
    integer,intent(in) :: ishell
    integer,intent(in) :: jshell
    integer,intent(in) :: mi,mj
    integer,intent(in) :: box_start1
    integer,intent(in) :: box_start2
    integer,intent(in) :: box_start3
    type(GRID_INFO),intent(in) :: grid
    type(CELL_INFO), intent(in) :: cell
    type(POINT),intent(in) :: atom_origin
    type(PSEUDO_SPECIES), intent(in) :: p_species(:) ! size(par%num_pspecies)

    ! Local Variables
    integer :: li, lj
    integer :: lup, mup
    integer :: i1,i2,i3
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    type(POINT) :: r_cell
    type(POINT) :: r_sphere
    type(POINT) :: r_sphere_unit
    real(kind=DP) :: rmag
    real(kind=DP) :: slmval
    real(kind=DP) :: grad_SLM_R(3)

    call utils_use_var(mi)
    call utils_use_var(mj)
    lup=0 ! this prevents 'uninitialized variable' compilation warnings
    mup=0 ! this prevents 'uninitialized variable' compilation warnings
    call utils_abort('Error in pseudo_get_aug_func_grad: &
         &routine not yet implemented')

    ! Find vector to origin of box
    box_origin = (box_start1-1)*grid%da1 + &
                 (box_start2-1)*grid%da2 + &
                 (box_start3-1)*grid%da3

    ! ndmh: vector from origin of box to centre of atom
    box_to_atom = atom_origin - box_origin

    li = p_species(isp)%ang_mom(ishell)
    lj = p_species(isp)%ang_mom(jshell)

    ! ndmh: check components of this vector in terms of lattice vectors
    ! ndmh: are all positive. If not, the box has been looped back into the
    ! ndmh: cell, so shift the origin back accordingly
    if ((box_to_atom.DOT.cell%b1) < 0.0_DP) then
       box_origin = box_origin - cell%a1
    end if
    if ((box_to_atom.DOT.cell%b2) < 0.0_DP) then
       box_origin = box_origin - cell%a2
    end if
    if ((box_to_atom.DOT.cell%b3) < 0.0_DP) then
       box_origin = box_origin - cell%a3
    end if

    do i3=1,box_n3
       do i2=1,box_n2
          do i1=1,box_n1

             r_cell = box_origin + (i1-1)*grid%da1 &
                                 + (i2-1)*grid%da2 &
                                 + (i3-1)*grid%da3

             r_sphere = r_cell - atom_origin
             rmag = geometry_magnitude(r_sphere)

             if (rmag>1d-10) then
                r_sphere_unit = (1.0_DP/rmag)*r_sphere
             else
                r_sphere_unit = 0.0_DP*r_sphere
             end if

             if (rmag>p_species(isp)%core_radius(ishell)) then
                aug_func_grad(i1,i2,i3,:) = 0.0_DP
                cycle
             end if

             ! Not functional yet!

             ! Get spherical harmonic at this \hat{r}
             slmval = sw_real_sph_harm(r_sphere%x,r_sphere%y,r_sphere%z, &
                  rmag,lup,mup)

             ! Get gradients of spherical harmonic at this \hat{r}
             call sw_grad_real_sph_harm(grad_SLM_R(:),r_sphere%x,r_sphere%y, &
                  r_sphere%z,rmag,lup,mup)

             ! S_LM(\hat{r}) * d/dR(g(|r-R|) part
             aug_func_grad(i1,i2,i3,1) = 1.0_DP * slmval * r_sphere_unit%x
             aug_func_grad(i1,i2,i3,2) = 1.0_DP * slmval * r_sphere_unit%y
             aug_func_grad(i1,i2,i3,3) = 1.0_DP * slmval * r_sphere_unit%z

             ! g(|r-R|) * d S_LM(\hat{r}) / dR(g(|r-R|) part
             aug_func_grad(i1,i2,i3,:) = aug_func_grad(i1,i2,i3,:) + &
                  grad_SLM_R(:) * 1.0_DP

          end do
       end do
    end do

  end subroutine pseudo_get_aug_func_grad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_nl_calculate_forces(nlps_forces, &
       sp_overlap,p_species,ngwfs_on_grid,ngwf_basis,proj_basis,proj_set, &
       cell,fftbox,pur_denskern,nsub,nat,first_atom_on_proc,num_atoms_on_proc,&
       orig_atom,global_atom_number,kb_denom, kpt) ! agrecokpt

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the nonlocal part of the ionic pseudopotential.                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1)  nlps_forces     : output : nonlocal forces                          !
    ! 2)  sp_overlap      : input  : ngwf-projector overlap matrix            !
    ! 3)  ngwfs_on_grid   : input  : NGWF data in ppd format                  !
    ! 4)  ngwf_basis      : input  : Function basis description for NGWFs     !
    ! 5)  proj_basis      : input  : Function basis description for projectors!
    ! 6)  proj_set        : input  : Nonlocal projector definitions           !
    ! 7)  pur_denskern    : input  : purified density kernel SPAM3            !
    !-------------------------------------------------------------------------!
    ! Written by Arash A. Mostofi, V0.0, 28th June 2004                       !
    !-------------------------------------------------------------------------!
    ! Modified by Arash A. Mostofi, 15th September 2004                       !
    !   - Distribution of the derivative with respect to atomic positions of  !
    !     the non-local potential matrix (nlfmtx_proc) amongst procs          !
    ! Modified by Peter Haynes for parallel SPAM3, July 2006                  !
    !   - Now all procs calculate contributions to the forces on all atoms    !
    !     and nlfmtx_proc is eliminated                                       !
    ! Modified by Nicholas Hine to improve efficiency of inner loop,          !
    !     November 2008.                                                      !
    ! Modified by Nicholas Hine in July 2009 to work with SPAM3.              !
    ! Modified by Andrea Greco in May 2016 to allow specification of k-point  !
    ! and compatibility with complex NGWFs.                                   !
    ! Modified to deal with embedding by Joseph Prentice, May 2018            !
    !=========================================================================!

    use datatypes, only: FUNCTIONS
    use comms, only: comms_barrier, comms_reduce, pub_my_proc_id, &
         pub_total_num_procs
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: PARAL_INFO
    use projectors, only: projectors_func_grad_ovlp_box
    use rundat, only: pub_num_spins
    use simulation_cell, only: CELL_INFO
    use sparse, only: sparse_get_element
    use sparse_embed, only: SPAM3_EMBED,sparse_embed_axpy,sparse_embed_create, &
         sparse_embed_destroy, sparse_embed_product, sparse_embed_transpose, &
         sparse_embed_scale, sparse_embed_transpose_structure, sparse_embed_copy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: nlps_forces(:,:)!size(1:3,par%nat)
    type(SPAM3_EMBED), intent(in) :: sp_overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis(:)
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid(:)
    type(FUNC_BASIS), intent(in) :: proj_basis(:)
    type(PROJECTOR_SET), intent(inout) :: proj_set(:)
    type(SPAM3_EMBED), intent(in) :: pur_denskern(pub_num_spins)
    type(PSEUDO_SPECIES), intent(in) :: p_species(:)!par%num_pspecies)
    type(CELL_INFO), intent(in) :: cell
    type(FFTBOX_INFO), intent(in) :: fftbox
    integer, intent(in) :: nsub,nat
    integer, intent(in) :: first_atom_on_proc(nsub,0:pub_total_num_procs-1)
    integer, intent(in) :: num_atoms_on_proc(nsub,0:pub_total_num_procs-1)
    integer, intent(in) :: orig_atom(nsub,nat)
    integer, intent(in) :: global_atom_number(nsub,nat)
    type(SPAM3_EMBED), intent(inout) :: kb_denom
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    type(SPAM3_EMBED) :: siGp_overlap,iGps_overlap
    type(SPAM3_EMBED) :: sp_overlap_dij,kb_ps_overlap
    type(SPAM3_EMBED) :: nl_force_mat(3)
    type(SPAM3_EMBED) :: kq,rkq
    integer :: cart, is, isub, jsub
    integer :: iat, orig_iat
    integer :: atom_proj, global_proj
    real(kind=DP) :: proj_force
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecocmplx: complex version for products
    type(SPAM3_EMBED) :: nl_force_mat_cmplx
    ! agrecokpt
    type(POINT) :: loc_kpt

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nl_calculate_forces'

    ! Start timer
    call timer_clock('pseudo_nl_calculate_forces',1)

    ! Initialise
    nlps_forces  = 0.0_DP

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid(1)%iscmplx

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! agrecokpt: specification of kpt allowed only in complex case
    if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.&
        loc_kpt%z/=0.0_DP) then
       if (.not.loc_cmplx) then
          call utils_abort('Error in pseudo_nl_calculate_forces: &
               &non-Gamma k-point requires complex projectors.')
       end if
    end if

    ! Create result matrices to hold nonlocal forces, projector-by-projector
    do cart=1,3
       nl_force_mat(cart)%structure='E'
       call sparse_embed_create(nl_force_mat(cart))
    end do

    ! agrecocmplx: in complex case, need temporary complex matrix to perform
    ! products between complex matrices
    if (loc_cmplx) then
       call sparse_embed_create(nl_force_mat_cmplx,nl_force_mat(1), &
               iscmplx=loc_cmplx)
    end if

    ! Create matrices to hold <NGWF_a|Proj_i> / D_i and <Proj_i|NGWF_a> / D_i
    ! agrecocmplx
    call sparse_embed_create(kb_ps_overlap,sp_overlap,trans=.true.,&
         iscmplx=loc_cmplx)
    call sparse_embed_create(sp_overlap_dij,sp_overlap)

    ! Get 1/KB denominators in diagonal matrix
    ! jcap: these are now precalculated and passed in

    ! Calculate the matrices <NGWF_a|Proj_i> / D_i and <Proj_i|NGWF_a> / D_i
    call sparse_embed_product(sp_overlap_dij,sp_overlap,kb_denom)
    call sparse_embed_transpose(kb_ps_overlap,sp_overlap_dij)

    ! Destroy kb_denom matrix to free up memory
    call sparse_embed_destroy(kb_denom)

    ! Create matrices to hold <phi|iG*proj> overlap matrix and transpose
    call sparse_embed_create(siGp_overlap,sp_overlap)
    call sparse_embed_create(iGps_overlap,kb_ps_overlap)

    ! Create temporary matrices kq and rkq
    call sparse_embed_create(kq,pur_denskern(1),sp_overlap_dij)
    ! agrecocmplx
    call sparse_embed_create(rkq,nl_force_mat(1),iscmplx=loc_cmplx)

    ! Loop over Cartesian directions
    do cart=1,3

       ! Calculate <phi|iG*proj> overlap matrix
       ! agrecokpt: at specified k-point
       ! jcap: loop over regions
       do isub=1,nsub
          do jsub=1,nsub
             call projectors_func_grad_ovlp_box(siGp_overlap%m(isub,jsub),&
                  ngwfs_on_grid(isub),ngwf_basis(isub),proj_basis(jsub),&
                  proj_set(jsub),fftbox,cell,cart,kpt=loc_kpt)
          end do
       end do

       ! Transpose it to get <iG*proj|phi> overlap matrix
       call sparse_embed_transpose(iGps_overlap,siGp_overlap)

       ! agrecocmplx
       ! complex case
       if (loc_cmplx) then
          ! agrecocmplx: reset nl_force_mat_cmplx for the
          ! new cartesian component being computed
          call sparse_embed_scale(nl_force_mat_cmplx, 0.0_DP)

          do is=1,pub_num_spins

             ! Calculate <proj_j|phi_b> K^ab <phi_a|iG.proj_i>
             call sparse_embed_product(kq,pur_denskern(is),siGp_overlap)
             call sparse_embed_product(rkq,kb_ps_overlap,kq)
             call sparse_embed_axpy(nl_force_mat_cmplx,rkq,1.0_DP)

             ! Calculate <iG.proj_j|phi_b> K^ab <phi_a|proj_i>
             call sparse_embed_product(kq,pur_denskern(is),sp_overlap_dij)
             call sparse_embed_product(rkq,iGps_overlap,kq)
             call sparse_embed_axpy(nl_force_mat_cmplx,rkq,1.0_DP)

          end do

          ! agrecocmplx: convert to real matrix since forces are real
          call sparse_embed_copy(nl_force_mat(cart),nl_force_mat_cmplx)

       ! real case
       else

          do is=1,pub_num_spins

             ! Calculate <proj_j|phi_b> K^ab <phi_a|iG.proj_i>
             call sparse_embed_product(kq,pur_denskern(is),siGp_overlap)
             call sparse_embed_product(rkq,kb_ps_overlap,kq)
             call sparse_embed_axpy(nl_force_mat(cart),rkq,1.0_DP)

             ! Calculate <iG.proj_j|phi_b> K^ab <phi_a|proj_i>
             call sparse_embed_product(kq,pur_denskern(is),sp_overlap_dij)
             call sparse_embed_product(rkq,iGps_overlap,kq)
             call sparse_embed_axpy(nl_force_mat(cart),rkq,1.0_DP)

          end do

       end if
    end do

    ! Destroy temporary matrices
    call sparse_embed_destroy(rkq)
    call sparse_embed_destroy(kq)
    call sparse_embed_destroy(iGps_overlap)
    call sparse_embed_destroy(siGp_overlap)
    call sparse_embed_destroy(kb_ps_overlap)
    call sparse_embed_destroy(sp_overlap_dij)
    ! agrecocmplx
    if (loc_cmplx) then
       call sparse_embed_destroy(nl_force_mat_cmplx)
    end if

    ! jcap: loop over regions
    do isub=1,nsub
       ! Loop over atoms
       do iat=first_atom_on_proc(isub,pub_my_proc_id), &
            first_atom_on_proc(isub,pub_my_proc_id)+num_atoms_on_proc(isub,pub_my_proc_id)-1

          ! Find atom number in input file order
          orig_iat = global_atom_number(isub,orig_atom(isub,iat))

          ! Loop over projectors on this atom
          do atom_proj=1,proj_basis(isub)%num_on_atom(iat)
             global_proj = proj_basis(isub)%first_on_atom(iat) + atom_proj - 1

             ! Loop over Cartesian co-ordinates
             do cart=1,3

                ! Find contribution of this projector to force on this atom
                ! from diagonal elements of nl_force_mat for this coordinate
                call sparse_get_element(proj_force,nl_force_mat(cart)%m(isub,isub), &
                     global_proj, global_proj)
                nlps_forces(cart,orig_iat) = nlps_forces(cart,orig_iat) + proj_force

             end do  ! cart

          end do  ! atom_proj

       end do   ! loc_iat

    end do ! isub
    ! Reduce result across procs
    call comms_barrier
    call comms_reduce('SUM',nlps_forces,3*nat)

    ! Destroy temporary matrices
    do cart=1,3
       call sparse_embed_destroy(nl_force_mat(cart))
    end do

    ! Stop timer
    call timer_clock('pseudo_nl_calculate_forces',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nl_calculate_forces'

  end subroutine pseudo_nl_calculate_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine subtract_or_add_coul(p_species,factor)

    !=================================================================!
    ! This subroutine is used to subtract the Coulomb (Z/r) potential !
    ! from the local part of each pseudopotential species before      !
    ! interpolation.                                                  !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 22/8/2001.                  !
    ! Modified on 24/1/2004 by Chris-Kriton Skylaris to work with the !
    ! PSEUDO_SPECIES type.                                            !
    ! Generalized on 27/5/2010 by Jacek Dziedzic to allow for adding  !
    ! 1/r back. Default behaviour is unchanged.                       !
    ! par removed as input by Joseph Prentice, June 2019              !
    !=================================================================!

    implicit none

    ! jd: Arguments
    real(kind=DP), optional, intent(in) :: factor
    type(PSEUDO_SPECIES), intent(inout) :: p_species(:)

    ! Local Variables
    integer :: row, species
    integer :: nsp1
    ! pa: changed ion_charge from int to real
    real(kind=DP) :: g_value, ion_charge
    real(kind=DP) :: fac
    integer :: num_pspecies

    ! ----------------------------------------------------------------------

    ! jd: Default behaviour is to subtract the Coulombic potential but alllow
    !     the caller to change it by passing a factor in an optional argument
    if (present(factor)) then
       fac = factor
    else
       fac = -1.0_DP
    end if

    ! jcap: Get num_pspecies from size of p_species, not par
    num_pspecies=size(p_species)

    ! cks: loop over different atomic species
    do species=1,num_pspecies

       ! ndmh: no need to subtract Coulomb potential from usps
       if (.not.p_species(species)%subtract_coul) cycle

       ion_charge = p_species(species)%ion_charge
       nsp1       = p_species(species)%n_shells + 1

       ! cks: loop over radial points of current species
       do row=2,p_species(species)%n_rad_pts

          g_value = p_species(species)%rad_proj_recip(row,nsp1)

          ! cks: subtract the attractive coulomb potential
          ! cks, 28/3/2004: modification for numerical consistency with CASTEP
          ! pa: changed to allow fractional ionic charge
          ! jd: changed to allow adding the potential back
          !     fac = -1.0 -> subtraction (default), fac = 1.0 -> addition
          p_species(species)%rad_locpot_recip(row) = &
               p_species(species)%rad_locpot_recip(row) - &
               fac * ion_charge / (g_value*g_value)

       end do

    end do

  end subroutine subtract_or_add_coul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pseudopotentials

