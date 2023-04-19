! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*- !
!================================================================!
!                                                                !
!             Hubbard DFT+U build module                         !
!                                                                !
!----------------------------------------------------------------!
! Written by David O'Regan in April 2009 and maintained          !
! subsequently by David O'Regan and Nicholas Hine.               !
! Re-worked to use the new HUBBARD_MODEL type by Nicholas Hine   !
! in 2011.                                                       !
! Addition of Constrained DFT code by Gilberto Teobaldi, 2011.   !
! Modifications by David Turban to include tensor correction 4,  !
! ie calculation of projector duals over the set of NGWFs        !
! in each CDFT-region                                            !
!================================================================!
! This module contains all information needed for DFT+U and cDFT !
! calculations.                                                  !
!                                                                !
! It defines the HUBBARD_MODEL type which contains all data      !
! relating to an instance of the Hubbard model, defining the     !
! correlated sites, the projectors used to define the correlated !
! manifold, the overlaps between NGWFs and Hubbard projectors    !
! of both hydrogenic and self-consistent NGWF projector type,    !
! the tensorial correction matrix needed to compute the          !
! occupation matrix when non-orthogonal projectors are used, and !
! the Hubbard Hamiltonian for the system in the projector        !
! representation.                                                !
!                                                                !
! Computation of the DFT+U energy and Hamiltonian is             !
! carried out in this module, as is mixing of the projectors     !
! to obtain projector-self-consistency.                          !
!                                                                !
! For details, see the following papers:                         !
!                                                                !
!  Linear-scaling DFT+U with full local orbital optimization     !
!  D. D. O'Regan, N. D. M. Hine, M. C. Payne and A. A. Mostofi,  !
!  Phys. Rev. B 85, 085107 (2012).                               !
!                                                                !
!  Optimised Projections for the Ab Initio Simulation of Large   !
!  and Strongly Correlated Systems,                              !
!  D. D. O'Regan, (Springer, Berlin, Heidelberg, 2012) 1st Ed.,  !
!  Springer Theses XVI, 225p., ISBN 978-3-642-23237-4.           !
!                                                                !
!  Subspace representations in ab initio methods for strongly    !
!  correlated systems,                                           !
!  D. D. O'Regan, M. C. Payne and A. A. Mostofi,                 !
!  Phys. Rev. B 83, 245124 (2011)                                !
!                                                                !
!  Projector self-consistent DFT+U using non-orthogonal          !
!  generalized Wannier functions                                 !
!  D. D. O'Regan, N. D. M. Hine, M. C. Payne and A. A. Mostofi,  !
!  Phys. Rev. B 82, 081102(R) (2010).                            !
!                                                                !
!================================================================!
! Minor adjustments to make this compatible with new embedding   !
! structures by Joseph Prentice, September 2018                  !
!================================================================!

module hubbard_build

  use datatypes, only: FUNCTIONS
  use constants, only: DP, stdout, VERBOSE
  use rundat, only: pub_debug_on_root
  use sparse, only: SPAM3
  use projectors, only: PROJECTOR_SET
  ! agreco: add dependence on SPAM3_ARRAY
  !use sparse_array, only: SPAM3_ARRAY

  implicit none

  private

  ! ndmh: "Container" type to store all information relating to an instance
  ! ndmh: of a Hubbard model projected from DFT+U orbitals
  type HUBBARD_MODEL

     ! Original atom number of each Hubbard atom
     integer, allocatable       :: orig(:)

     ! The hubbard_species index of each Hubbard atom
     integer, allocatable       :: species_number(:)

     ! A list of the NGWFs with greatest overlap with localised
     ! hydrogenic orbitals on this atom. Ordered according to
     ! the selection criterion.
     integer, allocatable       :: localisedngwfs(:,:) !(2l+1,nat_hub)

     ! ddor: The total occupancy for each site and spin
     real(kind=DP), allocatable :: population(:,:)

     ! ddor: The DFT+U energy for each site and spin
     real(kind=DP), allocatable :: energy(:,:)

     ! ddor: The total projection of the current set of NGWF projectors onto
     !       hydrogenic orbitals for each site
     real(kind=DP), allocatable :: current_projection(:)

     ! gibo: the cDFT gradient for each site, spin, and cdft-mode
     real(kind=DP), allocatable :: cdft_gradient(:,:,:)

     !gibo: complementary arrays for atom-resolved (not species-) cDFT U-pot.
     real(kind=DP), allocatable :: cdft_u_charge_up(:)
     real(kind=DP), allocatable :: cdft_u_charge_down(:)
     real(kind=DP), allocatable :: cdft_u_spin(:)

     ! ndmh: projector set type
     type(PROJECTOR_SET)        :: projectors

     ! ddor: arrays to store NGWFs, projectors and energies for projector
     !       self-consistency
     !real(kind=DP), allocatable :: consistency_ngwfs(:)
     ! agrecocmplx: use new FUNCTIONS type
     type(FUNCTIONS) :: consistency_ngwfs
     !real(kind=DP), allocatable :: consistency_projs(:)
     ! agrecocmplx: use new FUNCTIONS type
     type(FUNCTIONS) :: consistency_projs
     real(kind=DP), allocatable :: consistency_energies(:)
     real(kind=DP), allocatable :: consistency_projections(:)

     type(SPAM3) :: o_matrix    !ddor: Hubbard projector overlap matrix
     type(SPAM3) :: u_matrix    !ddor: DFT+U prefactor (U-J)/2
     type(SPAM3) :: j_matrix    !ddor: Hund's exchange J, not J/2
     type(SPAM3) :: up_matrix   !ddor: Spin up potential with alpha/spin-splitting
     type(SPAM3) :: down_matrix !ddor: Spin down potential with alpha/spin-splitting
     ! agreco: should these depend on k-points as well? use SPAM3_ARRAY instead
     !type(SPAM3_ARRAY) :: occupancy_matrix
     type(SPAM3), allocatable :: occupancy_matrix(:)
     !type(SPAM3_ARRAY) :: projector_ham
     type(SPAM3), allocatable :: projector_ham(:)

     ! ddor: The final NGWF overlap matrix of a Hubbard projector
     !       iteration is needed for the next
     !       in order to construct the covariant metric on the correlated sites
     type(SPAM3) :: consistency_overlap

     ! ddor: The current iteration number of  the Hubbard projector
     !       self-consistency cycle.
     integer :: consistency_iteration

     ! ddor: True if we wish to carry out DFT+U with hydrogenic projectors
     !       on the first HUBBARDSCF iteration, false for plain DFT.
     logical :: dftu_on_first_hubbardscf

     ! ddor: Make a copy of pub_maxit_ngwf_cg for restart of HUBBARDSCF
     integer :: store_maxit_ngwf_cg

     ! gom: the DFT+nu gradient for each site, spin, and cdft-mode
     real(kind=DP), allocatable :: nu_gradient(:,:,:)

     ! gom: the square of the total occupancy
     real(kind=DP), allocatable :: population_2(:,:)

     ! gom : the DFT+nu parameters
     type(SPAM3) :: u_matrix_up
     type(SPAM3) :: u_matrix_down
     real(kind=DP), allocatable :: nu_u1_up(:)
     real(kind=DP), allocatable :: nu_u1_down(:)
     real(kind=DP), allocatable :: nu_u2_up(:)
     real(kind=DP), allocatable :: nu_u2_down(:)
     real(kind=DP), allocatable :: nu_target_n_up(:)
     real(kind=DP), allocatable :: nu_target_n_down(:)
     real(kind=DP), allocatable :: nu_nu_up(:)
     real(kind=DP), allocatable :: nu_nu_down(:)
     real(kind=DP), allocatable :: group_pop(:,:)
     real(kind=DP), allocatable :: group_pop_2(:,:)
     integer, allocatable       :: nu_group(:)
     integer, allocatable       :: nu_projectors(:)
     integer                    :: num_groups

  end type HUBBARD_MODEL

  public :: HUBBARD_MODEL

  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  public :: hubbard_model_init
  public :: hubbard_model_exit
  public :: hubbard_build_matrices
  public :: hubbard_build_matrices_exit
  public :: cdft_build_matrices
  public :: hubbard_build_consist
  public :: hubbard_build_consist_exit
  public :: hubbard_spin_splitting_zero
  public :: hubbard_ham_matrix
  public :: hubbard_projector_ham
  public :: hubbard_energy_total
  public :: hubbard_energy_info
  public :: cdft_energy_total
  public :: cdft_energy_info
  public :: hubbard_projection_mtx
  public :: hubbard_projector_update
  public :: hubbard_test_convergence
  public :: hubbard_species_proj
  public :: hubbard_species_exit_proj
  public :: hubbard_calculate_forces
  public :: dft_nu_energy_total
  public :: dft_nu_energy_info
  ! ebl
  public :: hubbard_diago_occupancy
  public :: hubbard_atoms_occupancy
  public :: hubbard_cmp_rot

  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_model_init(hub,elements)

    !==========================================================!
    ! This subroutine allocates memory for the HUBBARD_MODEL   !
    ! type.                                                    !
    !----------------------------------------------------------!
    ! Written by Nicholas Hine in November 2011, based on code !
    ! by David O'Regan written in April 2009                   !
    !==========================================================!

    use bibliography, only: bibliography_cite
    use hubbard_init, only: h_species
    use ion, only: ELEMENT
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_cdft, &
         pub_cdft_atom_charge, pub_cdft_atom_spin, &
         pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor, &
         pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor, &
         pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
         pub_cdft_modes, pub_cdft_group_charge_up_only, &
         pub_cdft_group_charge_down_only, pub_cdft_multi_proj, &
         pub_num_spins, pub_dft_nu, pub_hubbard_unify_sites
    use utils, only: utils_alloc_check, utils_assert

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(ELEMENT), intent(in) :: elements(par%nat)

    ! Local Variables
    integer :: ierr
    integer :: dummy_h, iat_orig
    integer :: channels, sp
    integer :: group

    call bibliography_cite('DFTU')

    allocate(hub%orig(par%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%orig',ierr)
    allocate(hub%species_number(par%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%species_number',ierr)

    ! Set up orig and species_number arrays
    dummy_h = 0
    do sp=1,par%num_hub_species
       do iat_orig=1,par%nat ! The total number of atoms in the cell
          if (elements(iat_orig)%species_id == h_species(sp)%hub_species) then
             dummy_h = dummy_h + 1
             ! ddor: Original atom number of Hubbard atom
             hub%orig(dummy_h) = iat_orig
             hub%species_number(dummy_h) = sp
          end if
       end do
    end do

    call utils_assert(dummy_h .eq. par%nat_hub, &
         'Error in hubbard_model_init: final dummy_h does not match &
         &par%nat_hub. The two values were',dummy_h,par%nat_hub)

    !gibo: special counting of angular momentum channels for multi-proj cDFT
    if (pub_cdft_multi_proj) then

       channels = maxval(h_species(:)%cdft_num_proj)

       !gibo-23.07.13
       ! gibo: dummy transfer of NGWFs ang. momentum from elements to h_species
       do sp=1,par%num_hub_species
          do iat_orig=1,par%nat ! The total number of atoms in the cell
             if (elements(iat_orig)%species_id == h_species(sp)%hub_species) then
                h_species(sp)%num_ang_mom_shells = elements(iat_orig)%num_ang_mom_shells
                h_species(sp)%orb_ang_mom(:) = elements(iat_orig)%orb_ang_mom(:)
             end if
          end do
       end do

    else
       !gibo: standard DFT+U, single-channel cDFT treatment...
       channels = maxval(2*h_species(:)%hub_ang_mom + 1)

    endif

    allocate(hub%population(pub_num_spins,par%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%population',ierr)
    allocate(hub%energy(pub_num_spins,par%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%energy',ierr)
    allocate(hub%localisedngwfs(channels,par%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%localisedngwfs',ierr)
    allocate(hub%current_projection(par%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%current_projection',ierr)
    allocate(hub%population_2(pub_num_spins,par%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%population_2',ierr)

    ! ddor: Zero populations here, needed for DFT+U+J
    hub%population = 0.0_DP

    !gibo: for cDFT, allocate cDFT-gradient components
    CDFT: if (pub_cdft) then

       ! alocate cdft_gradient with one component for each cDFT-mode
       allocate(hub%cdft_gradient(pub_num_spins,par%nat_hub, &
            pub_cdft_modes), stat=ierr)
       call utils_alloc_check('hubbard_model_init','hub%cdft_gradient',ierr)

       ! charge-constrained modes
       if (pub_cdft_atom_charge .OR. pub_cdft_group_charge_donor .OR. &
            pub_cdft_group_charge_acceptor .OR. pub_cdft_group_charge_diff) then

          allocate(hub%cdft_u_charge_up(par%nat_hub),stat=ierr)
          call utils_alloc_check('hubbard_model_init',&
               'hub%cdft_u_charge_up',ierr)
          allocate(hub%cdft_u_charge_down(par%nat_hub),stat=ierr)
          call utils_alloc_check('hubbard_model_init',&
               'hub%cdft_u_charge_down',ierr)

          !initialise atom-specific U-potentials from species-specific input
          !Why? To be able to optimise different cDFT U-potentials for cDFT-atoms
          ! of the same species (thence targeted population) in different
          ! chemical environments...
          dummy_h = 0
          do sp=1,par%num_hub_species
             do iat_orig=1,par%nat ! The total number of atoms in the cell
                if (elements(iat_orig)%species_id == h_species(sp)%hub_species) then
                   dummy_h = dummy_h + 1

                   !gibo: modified with charge_UP/DOWN_only
                   if (pub_cdft_group_charge_up_only) then
                      hub%cdft_u_charge_up(dummy_h)   =  h_species(sp)%cdft_u_charge_up
                      hub%cdft_u_charge_down(dummy_h) =  0.0_DP ! ignore DOWN_potential/electrons

                   elseif (pub_cdft_group_charge_down_only) then
                      hub%cdft_u_charge_up(dummy_h)   =  0.0_DP ! ignore UP_potential/electrons
                      hub%cdft_u_charge_down(dummy_h) =  h_species(sp)%cdft_u_charge_down

                   else
                      hub%cdft_u_charge_up(dummy_h)   =  h_species(sp)%cdft_u_charge_up
                      hub%cdft_u_charge_down(dummy_h) =  h_species(sp)%cdft_u_charge_down

                   endif

                end if
             end do
          end do
       endif

       ! spin-constrained modes
       if (pub_cdft_atom_spin .OR. pub_cdft_group_spin_acceptor .OR. &
            pub_cdft_group_spin_donor .OR. pub_cdft_group_spin_diff) then

          allocate(hub%cdft_u_spin(par%nat_hub),stat=ierr)
          call utils_alloc_check('hubbard_model_init',&
               'hub%cdft_u_spin',ierr)

          !initialise atom-specific U-potentials from species-specific input
          !Why? To be able to optimise different cDFT U-potentials for cDFT-atoms
          ! of the same species (thence targeted population) in different
          ! chemical environments...
          dummy_h = 0
          do sp=1,par%num_hub_species
             do iat_orig=1,par%nat ! The total number of atoms in the cell
                if (elements(iat_orig)%species_id == h_species(sp)%hub_species) then
                   dummy_h = dummy_h + 1
                   hub%cdft_u_spin(dummy_h) =  h_species(sp)%cdft_u_spin
                end if
             end do
          end do
       endif
    endif CDFT

    !gom: for DFT+nu, allocate DFT_nu`-gradient components
    DFT_nu: if (pub_dft_nu) then

       !For single sites allocate DFT+nu_gradient with
       !two components (U1,U2) for each DFT+nu atom and spin
       if (.not.pub_hubbard_unify_sites) then

           allocate(hub%nu_gradient(pub_num_spins,par%nat_hub, 2), stat=ierr)
           call utils_alloc_check('hubbard_model_init','hub%nu_gradient',ierr)

           allocate(hub%nu_u1_up(par%nat_hub),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_u1_up',ierr)
           allocate(hub%nu_u1_down(par%nat_hub),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_u1_down',ierr)
           allocate(hub%nu_u2_up(par%nat_hub),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_u2_up',ierr)
           allocate(hub%nu_u2_down(par%nat_hub),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_u2_down',ierr)
           allocate(hub%nu_target_n_up(par%nat_hub),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_target_n_up',ierr)
           allocate(hub%nu_target_n_down(par%nat_hub),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_target_n_down',ierr)
           allocate(hub%nu_nu_up(par%nat_hub),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_nu_up',ierr)
           allocate(hub%nu_nu_down(par%nat_hub),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_nu_down',ierr)
           allocate(hub%nu_group(par%nat_hub),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_group',ierr)

           ! gom : assign values to the hubbard atom-resolved parameters
           dummy_h = 0
           do sp=1,par%num_hub_species
              do iat_orig=1,par%nat ! The total number of atoms in the cell
                 if (elements(iat_orig)%species_id == h_species(sp)%hub_species) then
                    dummy_h                       = dummy_h + 1
                    hub%nu_u1_up(dummy_h)         =  h_species(sp)%nu_u1_up
                    hub%nu_u1_down(dummy_h)       =  h_species(sp)%nu_u1_down
                    hub%nu_u2_up(dummy_h)         =  h_species(sp)%nu_u2_up
                    hub%nu_u2_down(dummy_h)       =  h_species(sp)%nu_u2_down
                    hub%nu_target_n_up(dummy_h)   =  h_species(sp)%nu_target_n_up
                    hub%nu_target_n_down(dummy_h) =  h_species(sp)%nu_target_n_down
                    hub%nu_nu_up(dummy_h)         =  h_species(sp)%nu_nu_up
                    hub%nu_nu_down(dummy_h)       =  h_species(sp)%nu_nu_down
                    hub%nu_group(dummy_h)         =  h_species(sp)%nu_group
                 end if
              end do
           end do
       else if (pub_hubbard_unify_sites) then

           !gom: DFT+nu_gradient with two components (U1,U2) for each DFT+nu
           !group and spin

           !gom: Number of groups is the maximum number appearing in groups column
           hub%num_groups=maxval(h_species(:)%nu_group)

           !gom: allocate memory for group_pop, group_pop_2
           allocate(hub%group_pop(pub_num_spins,hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init','hub%group_pop',ierr)
           allocate(hub%group_pop_2(pub_num_spins,hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init','hub%group_pop_2',ierr)

           allocate(hub%nu_projectors(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_projectors',ierr)

           allocate(hub%nu_gradient(pub_num_spins,hub%num_groups, 2), stat=ierr)
           call utils_alloc_check('hubbard_model_init','hub%nu_gradient',ierr)

           allocate(hub%nu_u1_up(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_u1_up',ierr)
           allocate(hub%nu_u1_down(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_u1_down',ierr)
           allocate(hub%nu_u2_up(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_u2_up',ierr)
           allocate(hub%nu_u2_down(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_u2_down',ierr)
           allocate(hub%nu_target_n_up(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_target_n_up',ierr)
           allocate(hub%nu_target_n_down(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_target_n_down',ierr)
           allocate(hub%nu_nu_up(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_nu_up',ierr)
           allocate(hub%nu_nu_down(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_nu_down',ierr)
           allocate(hub%nu_group(hub%num_groups),stat=ierr)
           call utils_alloc_check('hubbard_model_init',&
                'hub%nu_group',ierr)

          !gom : assign values to the hubbard group-resolved parameters
           hub%nu_projectors(:)=0
           do sp=1,par%num_hub_species
               do iat_orig=1,par%nat ! The total number of atoms in the cell
                   if (elements(iat_orig)%species_id == &
                                h_species(sp)%hub_species)then
                       do group=1,hub%num_groups ! The total number of groups in the cell
                         if (group == h_species(sp)%nu_group) then
                             hub%nu_u1_up(group)       = h_species(sp)%nu_u1_up
                             hub%nu_u1_down(group)     = h_species(sp)%nu_u1_down
                             hub%nu_u2_up(group)       = h_species(sp)%nu_u2_up
                             hub%nu_u2_down(group)     = h_species(sp)%nu_u2_down
                             hub%nu_target_n_up(group) = h_species(sp)%nu_target_n_up
                             hub%nu_target_n_down(group)=h_species(sp)%nu_target_n_down
                             hub%nu_nu_up(group)        =h_species(sp)%nu_nu_up
                             hub%nu_nu_down(group)      =h_species(sp)%nu_nu_down
                             hub%nu_group(group)        =h_species(sp)%nu_group

                             if (pub_cdft_multi_proj) then
                                 hub%nu_projectors(group)=hub%nu_projectors(group) &
                                                      +h_species(sp)%cdft_num_proj
                             else
                                 hub%nu_projectors(group)=hub%nu_projectors(group)&
                                                      +h_species(sp)%hub_ang_mom*2+ 1
                             endif

                         end if
                       end do
                   endif
               end do
           end do
       endif

    endif DFT_nu

  end subroutine hubbard_model_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_model_exit(hub)

    !==========================================================!
    ! This subroutine deallocates memory for the HUBBARD_MODEL !
    ! type.                                                    !
    !----------------------------------------------------------!
    ! Written by Nicholas Hine in November 2011, based on code !
    ! by David O'Regan written in April 2009                   !
    !==========================================================!

    use rundat, only: pub_cdft,&
         pub_cdft_atom_charge, pub_cdft_atom_spin, &
         pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor, &
         pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor, &
         pub_cdft_group_charge_diff, pub_cdft_group_spin_diff,&
         pub_dft_nu, pub_hubbard_unify_sites
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub

    ! Local Variables
    integer :: ierr

    ! Deallocate allocatable arrays in HUBBARD_MODEL

    !gibo: deallocate cDFT-specific arrays
    if (pub_cdft) then
       deallocate(hub%cdft_gradient,stat=ierr)
       call utils_dealloc_check('hubbard_model_exit','hub%cdft_gradient',ierr)
       !and atom specific potentials/target arrays
       if (pub_cdft_atom_charge .OR. pub_cdft_group_charge_acceptor .OR. &
            pub_cdft_group_charge_donor .OR. pub_cdft_group_charge_diff) then
          deallocate(hub%cdft_u_charge_up,stat=ierr)
          call utils_dealloc_check('hubbard_model_exit',&
               'hub%cdft_u_charge_up',ierr)
          deallocate(hub%cdft_u_charge_down,stat=ierr)
          call utils_dealloc_check('hubbard_model_exit',&
               'hub%cdft_u_charge_down',ierr)
       endif

       if (pub_cdft_atom_spin .OR. pub_cdft_group_spin_acceptor .OR. &
            pub_cdft_group_spin_donor .OR. pub_cdft_group_spin_diff) then
          deallocate(hub%cdft_u_spin,stat=ierr)
          call utils_dealloc_check('hubbard_model_exit',&
               'hub%cdft_u_spin',ierr)
       endif
    endif

    if (pub_dft_nu) then

        deallocate(hub%nu_group,stat=ierr)
        call utils_dealloc_check('hubbard_model_exit',&
             'hub%nu_group',ierr)
        deallocate(hub%nu_nu_down,stat=ierr)
        call utils_dealloc_check('hubbard_model_exit',&
             'hub%nu_nu_down',ierr)
        deallocate(hub%nu_nu_up,stat=ierr)
        call utils_dealloc_check('hubbard_model_exit',&
             'hub%nu_nu_up',ierr)
        deallocate(hub%nu_target_n_down,stat=ierr)
        call utils_dealloc_check('hubbard_model_exit',&
             'hub%nu_target_n_down',ierr)
        deallocate(hub%nu_target_n_up,stat=ierr)
        call utils_dealloc_check('hubbard_model_exit',&
             'hub%nu_target_n_up',ierr)
        deallocate(hub%nu_u2_down,stat=ierr)
        call utils_dealloc_check('hubbard_model_exit',&
             'hub%nu_u2_down',ierr)
        deallocate(hub%nu_u2_up,stat=ierr)
        call utils_dealloc_check('hubbard_model_exit',&
             'hub%nu_u2_up',ierr)
        deallocate(hub%nu_u1_down,stat=ierr)
        call utils_dealloc_check('hubbard_model_exit',&
             'hub%nu_u1_down',ierr)
        deallocate(hub%nu_u1_up,stat=ierr)
        call utils_dealloc_check('hubbard_model_exit',&
             'hub%nu_u1_up',ierr)
        deallocate(hub%nu_gradient, stat=ierr)
        call utils_dealloc_check('hubbard_model_init',&
                'hub%nu_gradient',ierr)

        if(pub_hubbard_unify_sites) then
            deallocate(hub%nu_projectors,stat=ierr)
            call utils_dealloc_check('hubbard_model_exit',&
                'hub%nu_projectors',ierr)
            deallocate(hub%group_pop_2,stat=ierr)
            call utils_dealloc_check('hubbard_model_exit',&
                'hub%group_pop_2',ierr)
            deallocate(hub%group_pop,stat=ierr)
            call utils_dealloc_check('hubbard_model_exit',&
                'hub%group_pop',ierr)
        endif


    endif

    deallocate(hub%population_2,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%population_2',ierr)
    deallocate(hub%current_projection,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%current_projection',ierr)
    deallocate(hub%localisedngwfs,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%localisedngwfs',ierr)
    deallocate(hub%energy,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%energy',ierr)
    deallocate(hub%population,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%population',ierr)
    deallocate(hub%species_number,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%species_number',ierr)
    deallocate(hub%orig,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%orig',ierr)

  end subroutine hubbard_model_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_build_matrices(hub,hub_proj_basis,is_cmplx)

    !==========================================================!
    ! This subroutine creates DFT+U SPAM3 matrices             !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in August 2009 !
    ! originally as hubbard_build_spam.                        !
    ! Modified by Andrea Greco on 14/01/2016 to allow use of   !
    ! complex matrices for complex NGWFs.                      !
    ! This was later simplified by Jose M Escartin (Sep 2017). !
    !==========================================================!

    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    ! agreco: added dependence on number of k-points
    use rundat, only: pub_cdft, pub_cdft_hubbard, &
         pub_task, pub_hubbard_restart, pub_hubbard_atomsolve, &
         pub_hub_on_the_fly, pub_num_spins, pub_num_kpoints, &
         pub_hubbard_unify_sites, pub_dft_nu, pub_hub_tensor_corr
    use sparse, only: sparse_create, sparse_destroy, sparse_put_element
    use utils, only: utils_alloc_check
    ! agreco: dependence on sparse_array routines
    !use sparse_array, only: sparse_array_create

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    logical, intent(in), optional :: is_cmplx

    ! Local Variables
    integer :: sp, hat_on_proc, hat, theatom
    integer :: hub_proj
    integer :: is, ierr
    real(kind=DP) :: half_u_minus_j, potential_up, potential_down
    real(kind=DP) :: half_u2_up, half_u2_down
    type(SPAM3) :: dummy1, dummy2, dummy3, dummy4
    ! agrecocmplx: local variable to set use of real/complex NGWFs
    logical :: loc_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_build_matrices'

    ! agrecocmplx
    if (present(is_cmplx)) then
       loc_cmplx = is_cmplx
    else
       loc_cmplx = .false.
    end if

    ! ndmh: Create projector overlap matrix
    ! gom: Combine all projectors into one Hubbard site
    if (pub_hub_tensor_corr .eq. 4) then
       hub%o_matrix%structure = 'I'
       ! agrecocmplx: not sure if we need complex Hubbard
       ! projectors when using complex NGWFs? need to check
       call sparse_create(hub%o_matrix,iscmplx=loc_cmplx)

    else if (pub_hub_tensor_corr .eq. 5) then

       ! dhpt: use dummy matrices to create WDV structure
       dummy1%structure = 'W'
       ! dhpt: ONLY PRELIMINARY, 'K' is only dense if
       ! dhpt: system < kernel cutoff
       dummy2%structure = 'K'
       dummy3%structure = 'V'
       ! agrecocmplx: need to be complex when complex
       ! NGWFs since they involve overlaps for NGWFs
       call sparse_create(dummy1,iscmplx=loc_cmplx)
       call sparse_create(dummy2,iscmplx=loc_cmplx)
       call sparse_create(dummy3,iscmplx=loc_cmplx)
       call sparse_create(dummy4, dummy2, dummy3)
       call sparse_create(hub%o_matrix, dummy1, dummy4)
       call sparse_destroy(dummy1)
       call sparse_destroy(dummy2)
       call sparse_destroy(dummy3)
       call sparse_destroy(dummy4)
    else
       hub%o_matrix%structure = 'G'
       ! agrecocmplx: need to check if this must be complex
       ! when using complex NGWFs
       call sparse_create(hub%o_matrix, iscmplx=loc_cmplx)
    endif

    ! ndmh: Create U and alpha matrices
    ! gom: check if DFT+nu is being used
    if(.not. pub_dft_nu) then
       ! gom: if not proceed with single u matrix
       if (pub_hubbard_unify_sites) then
          hub%u_matrix%structure ='I'
       else
          hub%u_matrix%structure = 'G'
       endif
       ! agrecocmplx: all these matrices should stay real anyway
       call sparse_create(hub%u_matrix)
       call sparse_create(hub%j_matrix,hub%u_matrix)
       call sparse_create(hub%up_matrix,hub%u_matrix)
       call sparse_create(hub%down_matrix,hub%u_matrix)
    else
       ! gom: if it is then create 2 u matrices
       if (pub_hubbard_unify_sites) then
          hub%u_matrix_up%structure   ='I'
          hub%u_matrix_down%structure ='I'
       else
          hub%u_matrix_up%structure   ='G'
          hub%u_matrix_down%structure ='G'
       endif
       ! agrecocmplx: all these matrices should stay real anyway
       call sparse_create(hub%u_matrix_up)
       call sparse_create(hub%u_matrix_down)
       call sparse_create(hub%up_matrix,hub%u_matrix_up)
       call sparse_create(hub%down_matrix,hub%u_matrix_down)
    endif

    ! agreco: include dependence on number of k-points?
    !allocate(hub%occupancy_matrix(pub_num_spins,pub_num_kpoints),stat=ierr)
    allocate(hub%occupancy_matrix(pub_num_spins),stat=ierr)
    call utils_alloc_check('hubbard_build_matrices','hub%occupancy_matrix',ierr)
    !allocate(hub%projector_ham(pub_num_spins,pub_num_kpoints),stat=ierr)
    allocate(hub%projector_ham(pub_num_spins),stat=ierr)
    call utils_alloc_check('hubbard_build_matrices','hub%projector_ham',ierr)
    ! agreco: spin and k-point dependent version using sparse_array routines
    !if (pub_hubbard_unify_sites) then
    !   call sparse_array_create(hub%occupancy_matrix, n_spins=pub_num_spins, &
    !          n_kpoints=pub_num_kpoints, structure='I')
    !   call sparse_array_create(hub%projector_ham, n_spins=pub_num_spins, &
    !          n_kpoints=pub_num_kpoints, structure='I')
    !else
    !   call sparse_array_create(hub%occupancy_matrix, n_spins=pub_num_spins, &
    !          n_kpoints=pub_num_kpoints, structure='G')
    !   call sparse_array_create(hub%projector_ham, n_spins=pub_num_spins, &
    !          n_kpoints=pub_num_kpoints, structure='G')
    !end if
    do is=1,pub_num_spins
       ! gom: Combine all projectors into one Hubbard site
       if (pub_hubbard_unify_sites) then
          hub%occupancy_matrix(is)%structure  = 'I'
          hub%projector_ham(is)%structure     = 'I'
       else
          hub%occupancy_matrix(is)%structure  = 'G'
          hub%projector_ham(is)%structure     = 'G'
       endif
       ! agreco: the occupancy matrix should stay real anyway
       call sparse_create(hub%occupancy_matrix(is))
       ! agreco: this should stay real anyway since product of real matrices
       call sparse_create(hub%projector_ham(is))
    end do

    ! gibo: Hubbard potentials only for pure Hubbard-only .OR. cdft_hubbard simulations
    if ((.not.pub_cdft) .or. (pub_cdft_hubbard)) then

       ! ddor: Loop over Hubbard atoms on my proc
       do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
          sp = hub%species_number(hat)
          theatom = par%distr_atom(hub%orig(hat))

          if (.not. pub_dft_nu) then  !Standard CDFT
             half_u_minus_j = (h_species(sp)%hub_u - h_species(sp)%hub_j) * 0.5_DP
             potential_up = h_species(sp)%hub_alpha - &
                  ( h_species(sp)%hub_spin_splitting * 0.5_DP )
             potential_down = h_species(sp)%hub_alpha + &
                  ( h_species(sp)%hub_spin_splitting * 0.5_DP )

             ! ddor: Loop over Hubbard projectors on my proc
             do hub_proj = hub_proj_basis%first_on_atom(theatom), &
                  hub_proj_basis%first_on_atom(theatom) + &
                  hub_proj_basis%num_on_atom(theatom) - 1

                call sparse_put_element(half_u_minus_j,&
                     hub%u_matrix,hub_proj,hub_proj)
                call sparse_put_element(h_species(sp)%hub_j,&
                     hub%j_matrix,hub_proj,hub_proj)
                call sparse_put_element(potential_up,&
                     hub%up_matrix,hub_proj,hub_proj)
                call sparse_put_element(potential_down,&
                     hub%down_matrix,hub_proj,hub_proj)
             end do

             ! gom : else populate DFT+nu matrices
          elseif (pub_dft_nu) then

             half_u2_up      = h_species(sp)%nu_u2_up * (-1.0_DP)
             half_u2_down    = h_species(sp)%nu_u2_down * (-1.0_DP)

             ! gom: V_up = U1/2
             potential_up    = h_species(sp)%nu_u1_up * 0.5_DP
             potential_down  = h_species(sp)%nu_u1_down * 0.5_DP

             ! gom: Loop over Hubbard projectors on my proc
             do hub_proj = hub_proj_basis%first_on_atom(theatom), &
                  hub_proj_basis%first_on_atom(theatom) + &
                  hub_proj_basis%num_on_atom(theatom) - 1

                call sparse_put_element(half_u2_up,&
                     hub%u_matrix_up,hub_proj,hub_proj)
                call sparse_put_element(half_u2_down,&
                     hub%u_matrix_down,hub_proj,hub_proj)

                call sparse_put_element(potential_up,&
                     hub%up_matrix,hub_proj,hub_proj)
                call sparse_put_element(potential_down,&
                     hub%down_matrix,hub_proj,hub_proj)

             end do
          endif

       end do

    end if

    ! gibo: for constrained-DFT create and add cDFT terms
    ! agreco: these matrices should be real anyway
    if (pub_cdft) call cdft_build_matrices(hub,hub_proj_basis)

    if ((pub_task == 'HUBBARDSCF') .or. pub_hubbard_restart &
         & .or. pub_hubbard_atomsolve .or. pub_hub_on_the_fly) then
       hub%consistency_overlap%structure = 'S'
       ! agreco: needs to be complex if complex NGWFs are used
       call sparse_create(hub%consistency_overlap,iscmplx=loc_cmplx)
    endif

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_build_matrices'

  end subroutine hubbard_build_matrices



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_build_matrices_exit(hub)

    !==========================================================!
    ! This subroutine destroys DFT+U SPAM3 matrices            !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in August 2009 !
    ! originally as hubbard_build_spam_exit.                   !
    !==========================================================!

    ! agreco: added dependence on number of k-points
    use rundat, only: pub_task, pub_hubbard_restart, pub_hubbard_atomsolve, &
         pub_hub_on_the_fly, pub_dft_nu, pub_num_spins, pub_num_kpoints
    use sparse, only: sparse_destroy
    use utils, only: utils_dealloc_check
    ! agreco: dependence on sparse_array routines
    !use sparse_array, only: sparse_array_destroy

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub

    ! Local Variables
    integer :: ierr
    integer :: is

    if (pub_task == 'HUBBARDSCF' .or. pub_hubbard_restart &
         & .or. pub_hubbard_atomsolve .or. pub_hub_on_the_fly) then
       call sparse_destroy(hub%consistency_overlap)
       call hubbard_build_consist_exit(hub)
    endif

    call sparse_destroy(hub%o_matrix)
    call sparse_destroy(hub%up_matrix)
    call sparse_destroy(hub%down_matrix)

    if (.not. pub_dft_nu) then
       call sparse_destroy(hub%j_matrix)
       call sparse_destroy(hub%u_matrix)
    else
       call sparse_destroy(hub%u_matrix_up)
       call sparse_destroy(hub%u_matrix_down)
    endif

    ! agreco: spin and k-point dependent version using
    ! routines in sparse_array module
    !call sparse_array_destroy(hub%occupancy_matrix)
    !call sparse_destroy(hub%projector_ham)
    do is=pub_num_spins,1,-1
       call sparse_destroy(hub%occupancy_matrix(is))
       call sparse_destroy(hub%projector_ham(is))
    end do
    deallocate(hub%projector_ham,stat=ierr)
    call utils_dealloc_check('hubbard_build_matrices_exit', &
         'hub%projector_ham',ierr)
    deallocate(hub%occupancy_matrix,stat=ierr)
    call utils_dealloc_check('hubbard_build_matrices_exit', &
         'hub%occupancy_matrix',ierr)

  end subroutine hubbard_build_matrices_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cdft_build_matrices(hub,hub_proj_basis)

    !==========================================================!
    ! This subroutine creates CDFT SPAM3 matrices              !
    !----------------------------------------------------------!
    ! Adapted to CDFT from hubbard_build_spam in November 2011 !
    ! by Gilberto Teobaldi.                                    !
    !==========================================================!

    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_cdft_atom_charge, pub_cdft_atom_spin, &
         pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor, &
         pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor, &
         pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
         pub_cdft_modes
    use sparse, only: SPAM3, sparse_create, sparse_destroy,&
         sparse_put_element, sparse_axpy, sparse_get_element

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    type(SPAM3) :: cdft_potential_up
    type(SPAM3) :: cdft_potential_down
    integer :: sp, hat_on_proc, hat, theatom
    integer :: hub_proj
    real(kind=DP) :: u_charge_up, u_charge_down
    real(kind=DP) :: u_spin_up, u_spin_down

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &cdft_build_matrices'

    ! gibo: create local spam matrixes for CDFT UP/DOWN potential
    call sparse_create(cdft_potential_up,   hub%up_matrix)
    call sparse_create(cdft_potential_down, hub%down_matrix)

    !====select CDFT-mode and iterate...
    ! gibo: atom_charge .OR. one group_charge mode only
    CHARGE_CDFT: if (pub_cdft_atom_charge .OR. pub_cdft_group_charge_diff .OR. &
         (pub_cdft_group_charge_acceptor.OR.pub_cdft_group_charge_donor) ) then
       !====select CDFT-mode and iterate...

       ! ddor: Loop over Hubbard atoms on my proc
       cdft_atoms_1 : do hat_on_proc = 1, &
            par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
          sp = hub%species_number(hat)
          theatom = par%distr_atom(hub%orig(hat))

          u_charge_up   = hub%cdft_u_charge_up(hat)
          u_charge_down = hub%cdft_u_charge_down(hat)

          ! ddor: Loop over Hubbard projectors on my proc
          hubprojs_cdft_1: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
               hub_proj_basis%first_on_atom(theatom) + &
               hub_proj_basis%num_on_atom(theatom) - 1

             call sparse_put_element(u_charge_up,&
                  cdft_potential_up,hub_proj,hub_proj)
             call sparse_put_element(u_charge_down,&
                  cdft_potential_down,hub_proj,hub_proj)

          end do hubprojs_cdft_1

       end do cdft_atoms_1

       !====select CDFT-mode and iterate...
       ! gibo: atom_charge .OR. one group_charge mode only
    endif CHARGE_CDFT
    !====select CDFT-mode and iterate...

    !====select CDFT-mode and iterate...
    ! gibo: atom_spin .OR. one group_spin mode only
    SPIN_CDFT: if (pub_cdft_atom_spin .OR. pub_cdft_group_spin_diff .OR. &
         (pub_cdft_group_spin_acceptor.OR.pub_cdft_group_spin_donor) ) then
       !====select CDFT-mode and iterate...

       ! depending on how many cdft-modes are active, either populate or augment
       ! the cdft_potential matrices with spin-terms.
       ! Spin-only single-cDFT mode case
       CDFT_SINGLE_MULTIPLE: if (pub_cdft_modes == 1) then

          ! ddor: Loop over Hubbard atoms on my proc
          cdft_atoms_2 : do hat_on_proc = 1, &
               par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             sp = hub%species_number(hat)
             theatom = par%distr_atom(hub%orig(hat))

             ! MIND THE SIGNS: Ecdft_spin = Us*[Nup -Ndown - DN_target]
             ! thus V_up = Us, and V_down = -Us
             u_spin_up   =  hub%cdft_u_spin(hat)
             u_spin_down = -hub%cdft_u_spin(hat)

             ! ddor: Loop over Hubbard projectors on my proc
             hubprojs_cdft_2: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
                  hub_proj_basis%first_on_atom(theatom) + &
                  hub_proj_basis%num_on_atom(theatom) - 1

                call sparse_put_element(u_spin_up,&
                     cdft_potential_up,hub_proj,hub_proj)
                call sparse_put_element(u_spin_down,&
                     cdft_potential_down,hub_proj,hub_proj)

             end do hubprojs_cdft_2

          end do cdft_atoms_2

          ! multiple-cDFt modes case
       elseif (pub_cdft_modes == 2) then

          ! ddor: Loop over Hubbard atoms on my proc
          cdft_atoms_3 : do hat_on_proc = 1, &
               par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             sp = hub%species_number(hat)
             theatom = par%distr_atom(hub%orig(hat))

             ! ddor: Loop over Hubbard projectors on my proc
             hubprojs_cdft_3: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
                  hub_proj_basis%first_on_atom(theatom) + &
                  hub_proj_basis%num_on_atom(theatom) - 1

                ! step_1: retrieve (charge-)potential terms
                call sparse_get_element(u_spin_up,&
                     cdft_potential_up,hub_proj,hub_proj)
                call sparse_get_element(u_spin_down,&
                     cdft_potential_down,hub_proj,hub_proj)

                ! step_2: AUGMENT (with spin-contributions) (charge-)potential terms
                ! MIND THE SIGNS: Ecdft_spin = Us*[Nup -Ndown - DN_target]
                ! thus V_up = Us, and V_down = -Us
                u_spin_up   = u_spin_up   + hub%cdft_u_spin(hat)
                u_spin_down = u_spin_down - hub%cdft_u_spin(hat)

                ! step_3: update potential matrices with augmented terms
                call sparse_put_element(u_spin_up,&
                     cdft_potential_up,hub_proj,hub_proj)
                call sparse_put_element(u_spin_down,&
                     cdft_potential_down,hub_proj,hub_proj)

             end do hubprojs_cdft_3

          end do cdft_atoms_3

       endif CDFT_SINGLE_MULTIPLE

       !====select CDFT-mode and iterate...
       ! gibo: atom_spin .OR. one group_spin mode only
    endif SPIN_CDFT
    !====select CDFT-mode and iterate...

    ! gibo: add cdft_potential_UP/DOWN to hub_UP/DOWN_matrix
    call sparse_axpy(hub%up_matrix,   cdft_potential_up,   1.0_DP)
    call sparse_axpy(hub%down_matrix, cdft_potential_down, 1.0_DP)

    ! gibo: destroy local spam matrixes
    call sparse_destroy(cdft_potential_up)
    call sparse_destroy(cdft_potential_down)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving cdft_build_matrices'

  end subroutine cdft_build_matrices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_build_consist(hub, hub_proj_basis, ngwf_basis, iscmplx)

    !==========================================================!
    ! This subroutine allocates arrays necessary for carrying  !
    ! out self-consistency over Hubbard projectors             !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in Oct 2009    !
    !==========================================================!

    use datatypes, only: data_functions_alloc, data_set_to_zero
    use bibliography, only: bibliography_cite
    use function_basis, only: FUNC_BASIS
    use utils, only: utils_alloc_check
    use rundat, only: pub_hub_conv_win

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    logical, intent(in) :: iscmplx

    ! Local Variables
    integer :: ierr

    call bibliography_cite('PROJSC')

    ! ddor: initialise matrices used for DFT+U
    !       with self-consistent projectors
    allocate(hub%consistency_energies(pub_hub_conv_win),stat=ierr)
    call utils_alloc_check('hubbard_build_consist', &
         'hub%consistency_energies',ierr)
    hub%consistency_energies = 0.0_DP
    allocate(hub%consistency_projections(pub_hub_conv_win),stat=ierr)
    call utils_alloc_check('hubbard_build_consist', &
         'hub%consistency_projections',ierr)
    hub%consistency_projections = 0.0_DP
    call data_functions_alloc(hub%consistency_ngwfs, &
         ngwf_basis%size_on_grid, iscmplx)
    call data_functions_alloc(hub%consistency_projs, &
         hub_proj_basis%size_on_grid, iscmplx)

    call data_set_to_zero(hub%consistency_projs)
    call data_set_to_zero(hub%consistency_ngwfs)
    hub%consistency_iteration = 1
    hub%dftu_on_first_hubbardscf = .true.

  end subroutine hubbard_build_consist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_build_consist_exit(hub)

    !==========================================================!
    ! This subroutine deallocates arrays for carrying          !
    ! out self-consistency over Hubbard projectors             !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in Oct 2009    !
    !==========================================================!

    use datatypes, only: data_functions_dealloc
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub

    ! Local Variables
    integer :: ierr

    ! ddor: Deallocate arrays used for DFT+U with self-consistent projectors
    call data_functions_dealloc(hub%consistency_projs)
    call data_functions_dealloc(hub%consistency_ngwfs)
    deallocate(hub%consistency_projections,stat=ierr)
    call utils_dealloc_check('hubbard_build_consist_exit', &
         'hub%consistency_projections',ierr)
    deallocate(hub%consistency_energies,stat=ierr)
    call utils_dealloc_check('hubbard_build_consist_exit', &
         'hub%consistency_energies',ierr)

  end subroutine hubbard_build_consist_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_spin_splitting_zero(hub,hub_proj_basis)

    !==========================================================!
    ! This subroutine sets any spin-splitting back to zero.    !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in August 2009 !
    !==========================================================!

    use comms, only: pub_my_proc_id, pub_on_root
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    use sparse, only: SPAM3, sparse_create, &
         sparse_scale, sparse_put_element

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: sp, hat_on_proc, hat, theatom
    integer :: hub_proj
    real(kind=DP) :: potential_up, potential_down

    if (SUM(ABS(h_species(:)%hub_spin_splitting)) .gt. 1.0e-10_DP) then

       if (pub_on_root) write(stdout,*) &
            'WARNING : Removing DFT+U spin-splitting potential.'

       call sparse_scale(hub%up_matrix,0.0_DP)
       call sparse_scale(hub%down_matrix,0.0_DP)

       h_species(:)%hub_spin_splitting = 0.0_DP

       ! ddor: Loop over Hubbard atoms on my proc
       do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
          sp = hub%species_number(hat)
          theatom = par%distr_atom(hub%orig(hat))

          potential_up = h_species(sp)%hub_alpha
          potential_down = h_species(sp)%hub_alpha

          ! ddor: Loop over Hubbard projectors on my proc
          do hub_proj = hub_proj_basis%first_on_atom(theatom), &
               hub_proj_basis%first_on_atom(theatom) + &
               hub_proj_basis%num_on_atom(theatom) - 1

             call sparse_put_element(potential_up,&
                  hub%up_matrix,hub_proj,hub_proj)
             call sparse_put_element(potential_down,&
                  hub%down_matrix,hub_proj,hub_proj)

          end do

       end do

    end if

  end subroutine hubbard_spin_splitting_zero


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_ham_matrix(hub,hubbard_ham,hub_overlap,hub_overlap_t)

    !==========================================================!
    ! This subroutine calculates the Hubbard DFT+U Hamiltonian !
    ! in SPAM3 form, including the components due to the alpha !
    ! and spin-splitting parameters, used for determining the  !
    ! value of U which is  consistent with the screening in    !
    ! the system and breaking magnetic symmetry, respectively. !
    !----------------------------------------------------------!
    ! Written by David O'Regan in August 2009                  !
    ! Modified by Nicholas Hine in April 2011 to separate      !
    ! creation of projector Ham from creation of Ham           !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use rundat, only: pub_num_spins
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_product, sparse_copy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(in) :: hub
    type(SPAM3), intent(inout) :: hubbard_ham(pub_num_spins)
    type(SPAM3), intent(in) :: hub_overlap
    type(SPAM3), intent(in) :: hub_overlap_t

    ! Local Variables
    type(SPAM3) :: v_hamiltonian_buffer
    integer :: is

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_ham_matrix'

    call timer_clock('hubbard_ham_matrix',1)

    ! dhpt: create buffer
    call sparse_create(v_hamiltonian_buffer,hub_overlap,hub%projector_ham(1))

    do is = 1, pub_num_spins

       ! ddor: V [ U/2 (1-2N) + a]  blanked to 'VG'
       call sparse_product(v_hamiltonian_buffer, hub_overlap, &
            hub%projector_ham(is), allow_mix_types=hub_overlap%iscmplx)

       ! ddor: Hamiltonian
       call sparse_product(hubbard_ham(is),&
            v_hamiltonian_buffer,hub_overlap_t)

    end do

    call sparse_destroy(v_hamiltonian_buffer)

    call timer_clock('hubbard_ham_matrix',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_ham_matrix'

  end subroutine hubbard_ham_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_projector_ham(hub,kernel,hub_overlap, &
       hub_overlap_t,hub_proj_basis)

    !==========================================================!
    ! Calculates the DFT+U Hamiltonian in the Hubbard projector!
    ! basis for a given spin channel and for the quadratic     !
    ! functional.                                              !
    !----------------------------------------------------------!
    ! Written by David O'Regan in September 2009.              !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_hub_functional, pub_cdft, pub_cdft_hubbard, &
         pub_hubbard_unify_sites, pub_dft_nu, pub_num_spins, &
         pub_hubbard_j_minority_term, pub_imag_thr
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_product, sparse_scale, sparse_axpy, sparse_copy, &
         sparse_put_element, sparse_take_real_part
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(SPAM3), intent(in) :: kernel(pub_num_spins)
    type(SPAM3), intent(in) :: hub_overlap
    type(SPAM3), intent(in) :: hub_overlap_t
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    type(SPAM3) :: site_buffer
    type(SPAM3), allocatable, dimension(:) :: occupancy_buffer
    type(SPAM3) :: occupancy_buffer2, occupancy_buffer3
    integer :: is, js
    integer :: ierr
    integer :: hat_on_proc, hat, sp, theatom
    integer :: channels,i,ihub
    ! agrecocmplx
    logical :: loc_cmplx
    type(SPAM3) :: occupancy_matrix_cmplx


    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_projector_ham'

    call timer_clock('hubbard_projector_ham',1)

    allocate(occupancy_buffer(pub_num_spins),stat=ierr)
    call utils_alloc_check('hubbard_projector_ham','occupancy_buffer',ierr)

    ! agrecocmplx
    loc_cmplx = kernel(1)%iscmplx

    ! gibo: calculate Hubbard terms only for pure Hubbard-only
    !       .OR. cdft_hubbard simulations
    ! *** WARNING: this if-construct is closed on line 1142
    HUBBARD_CDFT: if ((.not.pub_cdft) .OR. (pub_cdft_hubbard) ) then

       call sparse_create(site_buffer,hub_overlap_t,kernel(1))     ! WK
       !gom: Combine all projectors into one Hubbard site
       do is=1,pub_num_spins
          if (pub_hubbard_unify_sites) then
             occupancy_buffer(is)%structure = 'I'
          else
             occupancy_buffer(is)%structure = 'G'
          endif
          call sparse_create(occupancy_buffer(is))                 !WKV=G
       enddo

       ! agrecocmplx
       if (loc_cmplx) then
          call sparse_create(occupancy_matrix_cmplx, occupancy_buffer(1), &
               iscmplx=loc_cmplx)
       end if

       SPIN_loop: do is=1,pub_num_spins

          ! Calculate Hubbard term G
          ! WK
          call sparse_product(site_buffer,hub_overlap_t,kernel(is))
          ! WKV blanked to 'G'
          ! agrecocmplx
          if (loc_cmplx) then
             call sparse_product(occupancy_matrix_cmplx, site_buffer, hub_overlap)
             ! convert complex to real
             ! agrecocmplx: use safe routine to convert complex to real
             !call sparse_copy(occupancy_buffer(is), occupancy_matrix_cmplx)
             call sparse_take_real_part(occupancy_buffer(is), occupancy_matrix_cmplx, &
                  pub_imag_thr)
          ! standard real case
          else
             call sparse_product(occupancy_buffer(is),site_buffer,hub_overlap)
          end if

          ! ddor: The Hamiltonian for the quadratic energy penalty functional
          if (.not. pub_dft_nu .and. pub_hub_functional == 1 ) then

             ! (1 - 2 WKV) blanked to 'G'
             call sparse_scale(occupancy_buffer(is),-2.0_DP,1.0_DP)

             ! ddor: The Hamiltonian for the quartic energy penalty functional
          elseif (.not.pub_dft_nu .and. pub_hub_functional == 2) then

             ! Create temporary matrices with 'G' structure
             call sparse_create(occupancy_buffer2,occupancy_buffer(is))    !G
             call sparse_create(occupancy_buffer3,occupancy_buffer(is))    !G

             ! (WKV)^2
             call sparse_product(occupancy_buffer2,occupancy_buffer(is), &
                  occupancy_buffer2)
             ! -(WKV) (1 - WKV)
             call sparse_axpy(occupancy_buffer2,occupancy_buffer(is),-1.0_DP)
             ! (1 - 2 WKV)
             call sparse_copy(occupancy_buffer3,occupancy_buffer(is))
             call sparse_scale(occupancy_buffer3,-2.0_DP,1.0_DP)
             ! -(WKV) (1 - WKV) (1 - 2 WKV)
             call sparse_product(occupancy_buffer(is),occupancy_buffer2, &
                  occupancy_buffer3)
             ! 2 (WKV) (1 - WKV) (1 - 2 WKV)
             call sparse_scale(occupancy_buffer(is),-2.0_DP)

             ! Destroy temporary matrices
             call sparse_destroy(occupancy_buffer3)
             call sparse_destroy(occupancy_buffer2)

          elseif (.not. pub_dft_nu) then
             call utils_abort('Invalid value for pub_hub_functional.')
          end if

          if (.not.pub_dft_nu) then !Standard DFT+U
             ! ddor: DFT+U term in Hamiltonian, proportional to (U-J)/2
             call sparse_product(hub%projector_ham(is),hub%u_matrix,occupancy_buffer(is))
          elseif (pub_dft_nu) then !spin-dependent DFT+nu
             if (is==1) then
                ! Add diagonal U_2 ptential for DFT_nu
                call sparse_product(hub%projector_ham(is),hub%u_matrix_up, &
                     occupancy_buffer(is))
                ! Add diagonal U_1 ptential for DFT_nu
                call sparse_axpy(hub%projector_ham(is),hub%up_matrix,1.0_DP)
             elseif(is==2) then
                ! Add diagonal U_2 ptential for DFT_nu
                call sparse_product(hub%projector_ham(is),hub%u_matrix_down, &
                     occupancy_buffer(is))
                ! Add diagonal U_1 ptential for DFT_nu
                call sparse_axpy(hub%projector_ham(is),hub%down_matrix,1.0_DP)
             endif
          else
             call utils_abort('Subroutine hubbard_projector_ham: problems with&
                  &DFT_nu Hamiltonian!')
          endif

          ! Add contribution due to hub_alpha or spin-splitting if necessary
          ! gibo:*******  WARNING:
          ! gibo:*******  for (pub_cdft_hubbard), add also the cDFT terms
          !               [in hub_up/down_matrix]
          ! gibo:*******  (see also sbrtne hubbard_build_spam &
          !               cdft_build_matrices)
          ! gibo:*******  WARNING:
          if ( ( ( SUM(ABS(h_species(:)%hub_spin_splitting)) + &
               &SUM(ABS(h_species(:)%hub_alpha)) ) .gt. 0.0_DP ) &
               & .OR. (pub_cdft_hubbard)) then
             if (is == 2) then
                call sparse_axpy(hub%projector_ham(is),hub%down_matrix,1.0_DP)
             else
                call sparse_axpy(hub%projector_ham(is),hub%up_matrix,1.0_DP)
             end if
          end if

       end do SPIN_loop


       ! ddor: Carry out a DFT+U+J calculation
       !     : see B. Himmetoglu et al., Phys. Rev. B 84, 115108 (2011)
       if ((SUM(ABS(h_species(:)%hub_j)) .gt. 0.0_DP) .and. &
            (pub_hub_functional == 1) .and. (.not.pub_dft_nu)) then

          ! A temporary buffer
          call sparse_create(occupancy_buffer2,occupancy_buffer(1))    !G
          if (pub_hubbard_j_minority_term) then
             call sparse_create(occupancy_buffer3,occupancy_buffer(1))    !G
          endif

          ! ddor: Use a separate spin loop to ensure the occuapncies are gathered
          SPIN_loop_J: do is=1,pub_num_spins

             ! ddor: Make sure that we have the occupancy for the opposite spin
             js = MOD(is,2)+1

             ! ddor: Restore occupancy, 0.5 - 0.5(1 - 2 WKV) = WKV blanked to 'G'
             ! ddor: Make sure that it's the occupancy for the opposite spin
             call sparse_scale(occupancy_buffer(js),-0.5_DP,0.5_DP)

             ! ddor: Add extra potential corresponding to -J n^{minority}
             !     : giving total V^{sigma} = J (n^{-sigma} - 1^{minority})
             if (pub_hubbard_j_minority_term) then

                ! ddor: Check that we have already computed populations
                if (SUM(hub%population(:,:)) == 0.0_DP) then
                   call utils_abort('Error in hubbard_projector_ham: &
                        &missing populations for hubbard_j_minority_term')
                end if

                ! ddor: Loop over Hubbard atoms on my proc
                hubatoms_J: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

                   hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
                   theatom = par%distr_atom(hub%orig(hat))
                   channels = hub_proj_basis%num_on_atom(theatom)
                   sp = hub%species_number(hat)

                   if ( channels .ne. (2 * h_species(sp)%hub_ang_mom + 1) ) then
                      call utils_abort('Error in hubbard_projector_ham: &
                           &hub_proj_basis%num_on_atom and hub_ang_mom mismatch')
                   end if

                   ! ddor: If this spin is minority for this site
                   if ( hub%population(is,hat) .le. hub%population(js,hat) ) then

                      ! ddor: enter a Kronecker delta for this spin and atom
                      do i=1,channels
                         ihub = hub_proj_basis%first_on_atom(theatom) + i - 1
                         call sparse_put_element(1.0_DP,occupancy_buffer3, &
                              ihub,ihub)
                      end do

                   endif

                end do hubatoms_J

                ! ddor: Subtract minority delta function from majority occupancy
                !     : Note that occupancy_buffer no longer returns occupancies.
                call sparse_axpy(occupancy_buffer(js),occupancy_buffer3,-1.0_DP)

             endif

             ! ddor: Hund's exchange J term in Hamiltonian, J n^{-sigma}
             call sparse_product(occupancy_buffer2,hub%j_matrix,&
                  occupancy_buffer(js))

             ! ddor: Add to the Hamiltonian
             call sparse_axpy(hub%projector_ham(is),occupancy_buffer2,1.0_DP)

          end do SPIN_loop_J

          ! ddor: Destroy temporary buffer
          if (pub_hubbard_j_minority_term) then
             call sparse_destroy(occupancy_buffer3)
          endif
          call sparse_destroy(occupancy_buffer2)

       endif


       do is=pub_num_spins,1,-1
          call sparse_destroy(occupancy_buffer(is))
       enddo
       call sparse_destroy(site_buffer)

       ! agrecocmplx
       if (loc_cmplx) then
          call sparse_destroy(occupancy_matrix_cmplx)
       end if

       ! gibo: for CDFT-only simulations, add cDFT contribution to 'Hubbard' hamiltonian
       ! [in hub_up/down_matrix following call to hubbard_build_spam & cdft_build_matrices]
    elseif ((pub_cdft) .AND. (.not. pub_cdft_hubbard)) then

       do is=1,pub_num_spins
          if (is == 2) then
             call sparse_copy(hub%projector_ham(is),hub%down_matrix)
          else
             call sparse_copy(hub%projector_ham(is),hub%up_matrix)
          end if
       end do

    end if HUBBARD_CDFT

    deallocate(occupancy_buffer,stat=ierr)
    call utils_dealloc_check('hubbard_projector_ham','occupancy_buffer',ierr)

    call timer_clock('hubbard_projector_ham',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_projector_ham'

  end subroutine hubbard_projector_ham


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_energy_total(hub,hubbard_e, &
       pur_denskern, hub_proj_basis, hub_overlap, hub_overlap_t, &
       hubbard_e_alpha, hubbard_e_spin)  !gibo: optional DFT+alpha/spin energy

    !==========================================================!
    ! Calculates the Hubbard DFT+U energy contribution to the  !
    ! total energy, and writes out the occupancy matrix for    !
    ! each correlated site.                                    !
    !----------------------------------------------------------!
    ! Written by David O'Regan in August 2009                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_my_proc_id, comms_reduce
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_hub_functional, pub_output_detail, &
         pub_cdft, pub_cdft_hubbard, pub_hubbard_unify_sites, &
         pub_dft_nu, pub_num_spins, pub_hubbard_j_minority_term, &
         pub_imag_thr
    use sparse, only: SPAM3, sparse_axpy, sparse_copy, &
         sparse_product, sparse_create, sparse_destroy, &
         sparse_get_element, sparse_take_real_part
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(out) :: hubbard_e  ! Hubbard DFT+U energy
    !gibo: optional DFT+alpha (\sum a*Tr[n]) energy
    real(kind=DP), optional, intent(out) :: hubbard_e_alpha
    !gibo: optional DFT+spin-splitting (\sum s*(Tr[n_UP]-Tr[n_DN]) energy
    real(kind=DP), optional, intent(out) :: hubbard_e_spin
    type(SPAM3), intent(in) :: pur_denskern(pub_num_spins)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(SPAM3), intent(in) :: hub_overlap, hub_overlap_t

    ! Local Variables
    type(SPAM3) :: site_buffer
    type(SPAM3) :: energy_buffer
    type(SPAM3) :: energy_buffer2
    integer :: channels
    integer :: hat_on_proc, hat, sp, theatom
    integer :: is, js
    integer :: i,ihub
    real(kind=DP) :: occupancyelement, energyelement
    real(kind=DP) :: u,j,a,s
    !gibo: local DFT+alpha (\sum a*Tr[n]) energy
    real(kind=DP) :: local_hubbard_e_alpha
    !gibo: local DFT+spin-splitting (\sum s*(Tr[n_UP]-Tr[n_DN]) energy
    real(kind=DP) :: local_hubbard_e_spin
    ! agrecocmplx
    logical :: loc_cmplx
    type(SPAM3) :: occupancy_matrix_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering hubbard_energy_total'

    call timer_clock('hubbard_energy_total',1)

    ! gibo: Initialise energy and occupancies here not to mess with
    !       cdft_energy_total for CDFT
    hubbard_e = 0.0_DP
    hub%energy = 0.0_DP
    hub%population = 0.0_DP
    hub%population_2 = 0.0_DP

    ! agrecocmplx
    loc_cmplx = pur_denskern(1)%iscmplx

    ! for U-(J-)calculate runs, initialise DFT+alpha and DFT+sigma energy
    if (present(hubbard_e_alpha)) local_hubbard_e_alpha = 0.0_DP
    if (present(hubbard_e_spin) )  local_hubbard_e_spin = 0.0_DP
    ! gom: Calculate the U1, U2 energies for DFT+nu
    HUBBARD_CDFT_out: if (pub_dft_nu) then

       call dft_nu_energy_total(hub,hubbard_e, &
            pur_denskern, hub_proj_basis, hub_overlap, hub_overlap_t)

       ! gibo: calculate Hubbard terms only for pure Hubbard-only
       !       .OR. cdft_hubbard simulations
    elseif ((.not.pub_cdft) .OR. (pub_cdft_hubbard) ) then

       ! ddor: Sanity check for spin-symmetry breaking
       if ( (sum(abs(h_species(:)%hub_spin_splitting)) .ne. 0.0_DP) .and. &
            & (pub_num_spins .eq. 1) ) then
          call utils_abort(&
               'Error in hubbard_energy_total: &
               &An attempt has been made to break spin symmetry &
               &while performing an unpolarised calculation.')
       end if

       call sparse_create(site_buffer,hub_overlap_t,pur_denskern(1))
       ! ddor: Blank occupancy matrix and products thereof to
       !       block-diagonal matrices in order to prevent
       !       Hubbard sites with null overlap interfering

       ! agrecocmplx
       if (loc_cmplx) then
          call sparse_create(occupancy_matrix_cmplx, hub%occupancy_matrix(1), &
               iscmplx=loc_cmplx)
       end if

       ! gom: Combine all projectors into one Hubbard site
       if (pub_hubbard_unify_sites) then
          energy_buffer%structure   = 'I'
          energy_buffer2%structure  = 'I'
       else
          energy_buffer%structure   = 'G'
          energy_buffer2%structure  = 'G'
       endif

       call sparse_create(energy_buffer)
       call sparse_create(energy_buffer2)

       spin: do is=1,pub_num_spins

          ! WK
          call sparse_product(site_buffer,hub_overlap_t,pur_denskern(is))
          ! WKV blanked to 'G'
          ! agrecocmplx
          if (loc_cmplx) then
             call sparse_product(occupancy_matrix_cmplx, site_buffer, hub_overlap)
             ! convert complex to real
             ! agrecocmplx: use safe routine to convert complex to real
             !call sparse_copy(hub%occupancy_matrix(is), occupancy_matrix_cmplx)
             call sparse_take_real_part(hub%occupancy_matrix(is), occupancy_matrix_cmplx, &
                  pub_imag_thr)
          ! default real case
          else
             call sparse_product(hub%occupancy_matrix(is),site_buffer,hub_overlap)
          end if

          ! WKV.WKV blanked to 'G'
          call sparse_product(energy_buffer2,hub%occupancy_matrix(is), &
               hub%occupancy_matrix(is))

          call sparse_copy(energy_buffer,hub%occupancy_matrix(is))

          ! ddor: N - N^2
          call sparse_axpy(energy_buffer,energy_buffer2,-1.0_DP)

          ! ddor: Use quartic penalty functional instead of quadratic
          ! ddor: ( N - N^2 )^2
          if (pub_hub_functional == 2) then

             call sparse_copy(energy_buffer2,energy_buffer)
             call sparse_product(energy_buffer,energy_buffer2,energy_buffer2)

          end if

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             theatom = par%distr_atom(hub%orig(hat))
             channels = hub_proj_basis%num_on_atom(theatom)
             sp = hub%species_number(hat)

             if ( channels .ne. (2 * h_species(sp)%hub_ang_mom + 1) ) then
                call utils_abort('Error in hubbard_energy_total: &
                     &hub_proj_basis%num_on_atom and hub_ang_mom mismatch')
             end if

             ! ddor: Local copies of Hubbard parameters
             ! ddor: When J is non-zero, subtract it from U
             u = h_species(sp)%hub_u - h_species(sp)%hub_j
             a = h_species(sp)%hub_alpha
             s = h_species(sp)%hub_spin_splitting

             do i=1,channels
                ihub = hub_proj_basis%first_on_atom(theatom) + i - 1
                call sparse_get_element(energyelement,energy_buffer, &
                     ihub,ihub)
                hub%energy(is,hat) = hub%energy(is,hat) &
                     + energyelement
                call sparse_get_element(occupancyelement, &
                     hub%occupancy_matrix(is),ihub,ihub)
                hub%population(is,hat) = hub%population(is,hat) &
                     + occupancyelement
             end do

             ! ddor: For spin-polarised systems
             !       Calculate energy
             !       U/2 Tr[N - N^2] + alpha Tr[N] +
             !       spin_spitting/2 Tr[N_down - N_up]
             !       That is (U/2 + alpha + spin_splitting/2) Tr[N_down] +
             !            (U/2 + alpha - spin_splitting/2) Tr[N_up] - U/2 Tr[N^2]
             if (pub_num_spins == 2) then
                if (is == 1) then
                   hub%energy(is,hat) = &
                        ( ((-0.5_DP*s)+a)* hub%population(is,hat) ) + &
                        ( (u*0.5_DP)* hub%energy(is,hat) )

                   if (present(hubbard_e_alpha)) then
                      local_hubbard_e_alpha = local_hubbard_e_alpha &
                           + a*hub%population(is,hat)
                      local_hubbard_e_spin  = local_hubbard_e_spin  &
                           - 0.5_DP*s*hub%population(is,hat)
                   endif

                elseif (is == 2) then
                   hub%energy(is,hat) = &
                        ( ((0.5_DP*s)+a)* hub%population(is,hat) ) + &
                        ( (u*0.5_DP)* hub%energy(is,hat) )

                   if (present(hubbard_e_alpha)) then
                      local_hubbard_e_alpha = local_hubbard_e_alpha &
                           + a*hub%population(is,hat)
                      local_hubbard_e_spin  = local_hubbard_e_spin  &
                           + 0.5_DP*s*hub%population(is,hat)
                   endif
                end if
                ! ddor: For un-spin-polarised systems
                !       Calculate energy U Tr[N - N^2] + 2*alpha Tr[N]
                !       That is (U + 2*alpha) Tr[N] - U Tr[N^2]
                !       Since we are dealing with a one-electron density matrix
             elseif  (pub_num_spins == 1) then
                hub%energy(is,hat) = &
                     ( (2.0_DP*a)* hub%population(is,hat) ) + &
                     (u* hub%energy(is,hat) )

                if (present(hubbard_e_alpha)) local_hubbard_e_alpha = &
                     local_hubbard_e_alpha + 2.0_DP*a*hub%population(is,hat)

             end if

             hubbard_e = hubbard_e + hub%energy(is,hat)

          end do hubatoms

       end do spin


       ! ddor: Carry out a DFT+U+J calculation
       !     : see B. Himmetoglu et al., Phys. Rev. B 84, 115108 (2011)
       if ((SUM(ABS(h_species(:)%hub_j)) .gt. 0.0_DP) .and. &
            (pub_hub_functional == 1) .and. (.not.pub_dft_nu)) then

          ! ddor: Sanity check for Hund's exchange
          if (pub_num_spins == 1) &
               call utils_abort(&
               'Error in hubbard_energy_total: &
               &An attempt has been made to run DFT+U+J &
               &while performing an unpolarised calculation.')

          ! ddor: Take product of up and down occupancy matrices
          !       be careful of your metric!
          call sparse_product(energy_buffer,hub%occupancy_matrix(1), &
               hub%occupancy_matrix(2))

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms_J: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             theatom = par%distr_atom(hub%orig(hat))
             channels = hub_proj_basis%num_on_atom(theatom)
             sp = hub%species_number(hat)

             if ( channels .ne. (2 * h_species(sp)%hub_ang_mom + 1) ) &
                  call utils_abort('Error in hubbard_energy_total: &
                  &hub_proj_basis%num_on_atom and hub_ang_mom mismatch')

             ! ddor: Local copies of Hund's exchange parameter
             j = h_species(sp)%hub_j

             ! ddor: Reset J energy counter for this atom
             energyelement = 0.0_DP

             ! ddor: Trace spin off-diagonal occupancy for this site
             do i=1,channels
                ihub = hub_proj_basis%first_on_atom(theatom) + i - 1
                call sparse_get_element(occupancyelement,energy_buffer, &
                     ihub,ihub)
                energyelement = energyelement + occupancyelement
             end do

             energyelement = energyelement * j * 0.5_DP

             spin_J: do is=1,pub_num_spins

                ! ddor: Find the occupancy for the opposite spin
                js = MOD(is,2)+1

                ! ddor: Add extra energy corresponding to -J n^{minority}
                if (pub_hubbard_j_minority_term .and. &
                     ( hub%population(is,hat) .le. &
                     hub%population(js,hat) )) then

                   ! ddor: Having made sure we have the minority channel,
                   ! ddor: add term that further penalises its occupancy
                   energyelement = energyelement - ( j * hub%population(is,hat) )

                endif

                ! ddor: Add this atoms's J energy contribution to
                !     : the atomic and global energy counts
                hub%energy(is,hat) = hub%energy(is,hat) + energyelement
                hubbard_e = hubbard_e + energyelement

             end do spin_J

          end do hubatoms_J

       endif ! ddor: DFT+U+J

       ! gibo: calculate DFT+U occupancy matrices and energies for pure
       ! gibo: Hubbard-only simulation
       HUBBARD_CDFT_in: if (.not.pub_cdft) then

          call comms_reduce('SUM',hub%energy)
          call comms_reduce('SUM',hub%population)

          ! ddor: Write out the DFT+U occupancy matrices and energies
          call comms_reduce('SUM', hubbard_e )
          if (pub_output_detail >= VERBOSE) &
               call hubbard_energy_info(hub,hub_proj_basis)

          !gibo: for U-(J-)-calculate runs, reduce also hubbard_e_alpha/spin
          if (present(hubbard_e_alpha)) then
             call comms_reduce('SUM',local_hubbard_e_alpha)
             hubbard_e_alpha = local_hubbard_e_alpha
          endif
          if (present(hubbard_e_spin))  then
             call comms_reduce('SUM',local_hubbard_e_spin)
             hubbard_e_spin = local_hubbard_e_spin
          endif

       end if HUBBARD_CDFT_in

       ! agrecocmplx
       if (loc_cmplx) then
          call sparse_destroy(occupancy_matrix_cmplx)
       end if

       call sparse_destroy(energy_buffer2)
       call sparse_destroy(energy_buffer)
       call sparse_destroy(site_buffer)

    end if HUBBARD_CDFT_out

    ! gibo: for cDFT runs, calculate cDFT energy energy and add it to hubbard_e
    CDFT_ENERGY: if (pub_cdft) then

       call cdft_energy_total(hub,hubbard_e, &
            pur_denskern, hub_proj_basis, hub_overlap, hub_overlap_t)

    end if CDFT_ENERGY

    call timer_clock('hubbard_energy_total',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_energy_total'

  end subroutine hubbard_energy_total


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_energy_info(hub,hub_proj_basis)

    !==========================================================!
    ! Writes out occupancies and energies of Hubbard sites.    !
    !----------------------------------------------------------!
    ! Written by David O'Regan in November 2009                !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: comms_bcast, pub_on_root, &
         pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_hubbard_compute_u_or_j, pub_num_spins
    use sparse, only: SPAM3, sparse_get_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(in) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: species, this_species_atom ! ddor: Loop over species
    integer :: jj,hat,project,spin,theatom, proc
    integer :: ierr
    real(kind=DP) :: local_population, local_moment, local_energy
    !gibo: local variables for species-averaged population*magnetic moment
    real(kind=DP) :: av_species_population, av_species_moment
    real(kind=DP),allocatable :: occupancy_block(:,:,:)
    character(10) :: fmt, tmp

    allocate(occupancy_block(hub_proj_basis%max_on_atom, &
         hub_proj_basis%max_on_atom,pub_num_spins),stat=ierr)
    call utils_alloc_check('hubbard_energy_info','occupancy_block',ierr)

    do species=1,par%num_hub_species

       this_species_atom = 0

       if (pub_hubbard_compute_u_or_j) then
          !gibo: initialise species-averages to zero
          av_species_population = 0.0_DP
          av_species_moment     = 0.0_DP
       endif

       atomloop: do hat=1,par%nat_hub

          if (hub%species_number(hat) == species) then

             this_species_atom = this_species_atom + 1
             if ((this_species_atom .gt. 2).AND.&
                  (.not.pub_hubbard_compute_u_or_j)) exit atomloop

             theatom = par%distr_atom(hub%orig(hat))
             project = hub_proj_basis%num_on_atom(theatom)
             proc = par%proc_of_atom(hub%orig(hat))

             if (proc==pub_my_proc_id) then
                do spin=1,pub_num_spins
                   call sparse_get_block(occupancy_block(:,:,spin), &
                        hub%occupancy_matrix(spin),theatom,theatom)
                end do
             end if
             call comms_bcast(proc,occupancy_block)

             if (.not.pub_on_root) cycle

             write(stdout,'(/a)') repeat('#',80)

             ! CW
             if(proc==pub_my_proc_id) then
               call hubbard_diago_occupancy(hubbard_atoms_occupancy(hub_proj_basis,hub,hat,1),project,hat)
             endif
             ! end CW

             ! ddor: TELL ME ABOUT THIS ATOM
             write(stdout,'(a,i6,a,a)') 'DFT+U information on atom ', &
                  this_species_atom,' of Hubbard species ',&
                  h_species(species)%hub_species

             do spin=1,pub_num_spins

                write(stdout,'(a)') repeat('#',80)

                if (pub_num_spins == 2) then
                   ! ddor: WRITE OUT INFORMATION ON OCCUPANCY MATRIX
                   write(stdout,'(a,2(i6,a))') 'Occupancy matrix of &
                        &Hubbard site ',hat,' and spin ',spin,' is '
                else

                end if

                ! ddor: write header line
                write(stdout,'(a)',advance='no') '  m_l ='
                do jj=-(project-1)/2,(project-1)/2
                   write(stdout,'(i5,7x)',advance='no') jj
                end do
                write(stdout,'(1x)')

                ! gom: occupancy output
                write(tmp,'(i6)') project
                write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f12.8)'
                do jj=1,project
                   write(stdout,fmt) occupancy_block(jj,1:project,spin)
                end do

                if (MAXVAL(occupancy_block(:,:,spin)) .gt. 1.0_DP) then
                   write(stdout,'(a,2(i6,a))') 'WARNING: OCCUPANCY MATRIX &
                        &of Hubbard site ',&
                        hat,' and spin ',spin,' exceeeds 1.'
                end if
                if (hub%energy(spin,hat) .lt. 0.0_DP) then
                   write(stdout,'(a,2(i6,a))') 'WARNING: DFT+U ENERGY &
                        &of Hubbard site ',&
                        hat,' and spin ',spin,' is negative.'
                end if

             end do

             write(stdout,'(a)') repeat('#',80)

             ! ddor: WRITE OUT INFORMATION ON OCCUPANCY OF SITE
             if (pub_num_spins==2) then
                local_population = sum(hub%population(:,hat))
                local_energy = hub%energy(1,hat) + &
                     hub%energy(2,hat)
             else
                local_population = 2.0_DP*hub%population(1,hat)
                local_energy = hub%energy(1,hat)
             end if

             if (pub_hubbard_compute_u_or_j) then
                !gibo: augment species-averaged population
                av_species_population = av_species_population + local_population
             endif

             write(stdout,'(a,i6,a,f12.8,a)') 'Total occupancy &
                  &of Hubbard site ', hat,' is       ', &
                  &local_population, ' e'

             if (pub_num_spins==2) then
                local_moment =  &
                     &hub%population(1,hat) - hub%population(2,hat)
                write(stdout,'(a,i6,a,f12.8,a)') 'Local magnetic moment &
                     &of Hubbard site ', hat,' is ', local_moment,' mu_B'
             end if

             if (pub_hubbard_compute_u_or_j) then
                !gibo: augment species-averaged population
                av_species_moment = av_species_moment + local_moment
             endif

             ! ddor: WRITE OUT INFORMATION ON DFT+U ENERGY
             write(stdout,'(a,i6,a,f12.8,a)') 'DFT+U energy of Hubbard &
                  &site ',hat,' is          ', local_energy,' Ha'

             write(stdout,'(a/)') repeat('#',80)

          end if

       end do atomloop

       if (pub_hubbard_compute_u_or_j) then
          if (pub_on_root) then
             !gibo: print-species resolved average&cumulative
             !      populations and magnetic moments

             write(stdout,'(a)')'#alpha' // repeat('#',74)

             write(stdout,'(a,a,a,i4)') 'Number of atoms of Hubbard species    ', &
                  h_species(species)%hub_species,'     ', this_species_atom
             write(stdout,'(a,a,a,f12.8,a)') 'Average population of Hubbard species ', &
                  h_species(species)%hub_species,'      ', &
                  av_species_population/this_species_atom,' e'
             write(stdout,'(a,a,a,f12.8,a)') 'Cumulative population of Hubbard species ', &
                  h_species(species)%hub_species,'    ', av_species_population,' e'
             write(stdout,'(a,a,f12.8,a)') 'Average magnetic moment of Hubbard species  ', &
                  h_species(species)%hub_species, &
                  av_species_moment/this_species_atom,' mu_B'
             write(stdout,'(a,a,f12.8,a)') 'Cumulative magnetic mom. of Hubbard species ', &
                  h_species(species)%hub_species, av_species_moment,' mu_B'

             write(stdout,'(a)') repeat('#',80)
          endif
       endif

    end do

    deallocate(occupancy_block,stat=ierr)
    call utils_dealloc_check('hubbard_energy_info','occupancy_block',ierr)

  end subroutine hubbard_energy_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cdft_energy_total(hub,hubbard_e, &
       pur_denskern, hub_proj_basis, hub_overlap, hub_overlap_t)

    !==========================================================!
    ! Calculates the cDFT energy contribution to the Hubbard   !
    ! (thence total) energy,                                   !
    ! and writes out the occupancy matrix for                  !
    ! each correlated site.                                    !
    !----------------------------------------------------------!
    ! Adapted from hubbard_energy_total by gibo in Nov. 2011   !
    !==========================================================!

    use comms, only: pub_on_root, pub_my_proc_id, comms_reduce
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_output_detail, pub_cdft_hubbard, &
         pub_cdft_atom_charge, pub_cdft_atom_spin, &
         pub_cdft_group_charge_acceptor, pub_cdft_group_charge_donor, &
         pub_cdft_group_spin_acceptor, pub_cdft_group_spin_donor, &
         pub_cdft_group_charge_acceptor_u, pub_cdft_group_charge_donor_u, &
         pub_cdft_group_spin_acceptor_u, pub_cdft_group_spin_donor_u, &
         pub_cdft_group_charge_acceptor_target, pub_cdft_group_charge_donor_target, &
         pub_cdft_group_spin_acceptor_target, pub_cdft_group_spin_donor_target, &
         pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
         pub_cdft_group_charge_diff_u, pub_cdft_group_spin_diff_u, &
         pub_cdft_group_charge_diff_target, pub_cdft_group_spin_diff_target, &
         pub_cdft_group_charge_up_only, &
         pub_cdft_group_charge_down_only, pub_cdft_multi_proj, pub_num_spins, &
         pub_imag_thr
    use sparse, only: SPAM3, sparse_axpy, sparse_copy, &
         sparse_product, sparse_create, sparse_destroy, &
         sparse_get_element, sparse_take_real_part
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(inout) :: hubbard_e  ! Hubbard DFT+U energy
    type(SPAM3), intent(in) :: pur_denskern(pub_num_spins)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(SPAM3), intent(in) :: hub_overlap, hub_overlap_t

    ! Local Variables
    type(SPAM3) :: site_buffer

    integer :: channels
    integer :: hat_on_proc, hat, sp, theatom
    integer :: is
    integer :: i,ihub
    integer :: mode_counter  !gibo: counter for cDFT-mode
    real(kind=DP) :: occupancyelement
    real(kind=DP) :: cdft_e      ! total cDFT energy
    real(kind=DP) :: tmp_cdft_e  ! group contribution(s) to cDFT energy
    real(kind=DP) :: u_charge_up, u_charge_down
    real(kind=DP) :: u_spin
    real(kind=DP) :: n_t_up, n_t_down, delta_n_t
    real(kind=DP) :: cdft_penalty
    real(kind=DP) :: pop_acceptor, pop_donor
    real(kind=DP) :: spinmom_acceptor, spinmom_donor
    ! agrecocmplx
    logical :: loc_cmplx
    type(SPAM3) :: occupancy_matrix_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &cdft_energy_total'

    call timer_clock('cdft_energy_total',1)

    call sparse_create(site_buffer,hub_overlap_t,pur_denskern(1))
    ! ddor: Blank occupancy matrix and products thereof to
    !       block-diagonal matrices in order to prevent
    !       Hubbard sites with null overlap interfering

    ! initialise cDFT-energy to zero
    cdft_e = 0.0_DP

    ! initialise cdft_gradient to zero
    hub%cdft_gradient = 0.0_DP

    ! initialise cDFT-mode counter to zero
    mode_counter = 0

    ! initialise cDFT-group energy contribution to zero
    tmp_cdft_e = 0.0_DP

    ! agrecocmplx
    loc_cmplx = pur_denskern(1)%iscmplx

    ! NEW-30.11.11. For cDFT+U runs [pub_cdft_hubbard=.T.],
    ! we have already calculated the populations in hubbard_energy_total
    ! so, no need to recalculate them here...

    CDFT_HUBBARD: if (.not.pub_cdft_hubbard) then

       ! Calculate population first,
       ! then work out energy depending on given cdft-mode.
       ! Slightly less efficient (we repeat the spin_ and hubatoms_ loops) for
       ! cdft_atom_charge = T, yet save a lot of repeated coding

       ! agrecocmplx
       if (loc_cmplx) then
          call sparse_create(occupancy_matrix_cmplx, hub%occupancy_matrix(1), &
               iscmplx=loc_cmplx)
       end if

       spin: do is = 1, pub_num_spins

          ! WK
          call sparse_product(site_buffer,hub_overlap_t,pur_denskern(is))
          ! WKV blanked to 'G'
          ! agrecocmplx
          if (loc_cmplx) then
             call sparse_product(occupancy_matrix_cmplx,site_buffer,hub_overlap)
             ! convert complex to real
             ! agrecocmplx: use safe routine to convert complex to real
             !call sparse_copy(hub%occupancy_matrix(is), occupancy_matrix_cmplx)
             call sparse_take_real_part(hub%occupancy_matrix(is), occupancy_matrix_cmplx, &
                  pub_imag_thr)
          ! standard default case
          else
             call sparse_product(hub%occupancy_matrix(is),site_buffer,hub_overlap)
          end if

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms : do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             theatom = par%distr_atom(hub%orig(hat))
             channels = hub_proj_basis%num_on_atom(theatom)
             sp = hub%species_number(hat)

             if  ((.not.(pub_cdft_multi_proj)) .AND. (channels .ne.  &
                  (2 * h_species(sp)%hub_ang_mom + 1)))  then
                call utils_abort('Error in cdft_energy_total: &
                     &hub_proj_basis%num_on_atom and hub_ang_mom mismatch')
             end if

             do i=1,channels
                ihub = hub_proj_basis%first_on_atom(theatom) + i - 1
                call sparse_get_element(occupancyelement,hub%occupancy_matrix(is), &
                     ihub,ihub)
                hub%population(is,hat) = hub%population(is,hat) &
                     + occupancyelement
             end do

          end do hubatoms

       end do spin

       ! agrecocmplx
       if (loc_cmplx) then
          call sparse_destroy(occupancy_matrix_cmplx)
       end if

    endif CDFT_HUBBARD

    !******************************************************************************
    ! Once the cDFT populations are calculated,                                   *
    ! deal with them according to the selected cDFT-mode                          *
    !******************************************************************************

    !====select CDFT-ATOM-mode and iterate...
    ! gibo: ATOM-CHARGE constrained run
    CDFT_ATOM: if (pub_cdft_atom_charge) then
       !====select CDFT-ATOM-mode and iterate...

       ! increase mode_counter (we need it for the cdft_gradient)
       mode_counter = mode_counter + 1

       ! ddor: Loop over Hubbard atoms on my proc
       hubatoms_1 : do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
          sp = hub%species_number(hat)

          ! Cycle loop if present species is not cDFT-active.
          ! For mixed DFT+U and cDFT+U make sure hubbard_energy contributions
          ! in hub%energy(is,hat) is not lost
          if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) then
             do is = 1, pub_num_spins
                cdft_e = cdft_e + hub%energy(is,hat)
             enddo
             CYCLE hubatoms_1
          endif

          !gibo: local copies of CONSTRAINED_DFT CHARGE parameters
          u_charge_up    = hub%cdft_u_charge_up(hat)
          u_charge_down  = hub%cdft_u_charge_down(hat)

          n_t_up      = h_species(sp)%cdft_target_up
          n_t_down    = h_species(sp)%cdft_target_down

          spin_1: do is = 1, pub_num_spins

             if (is == 1) then

                ! U_q (Tr[n_up] - n_t_up)
                cdft_penalty = hub%population(is,hat) - n_t_up
                hub%energy(is,hat) =  hub%energy(is,hat) + u_charge_up * cdft_penalty

                cdft_e = cdft_e + hub%energy(is,hat)

             elseif (is == 2) then

                ! U_q (Tr[n_down] - n_t_down)
                cdft_penalty = hub%population(is,hat) - n_t_down
                hub%energy(is,hat) =  hub%energy(is,hat) + u_charge_down * cdft_penalty

                cdft_e = cdft_e + hub%energy(is,hat)

             end if

             ! store the gradient
             hub%cdft_gradient(is,hat,mode_counter) = cdft_penalty

          end do spin_1

       end do hubatoms_1


       !====select CDFT-ATOM-mode and iterate...
       ! gibo: ATOM-SPIN constrained run
    elseif (pub_cdft_atom_spin) then
       !====select CDFT-ATOM-mode and iterate...

       ! increase mode_counter (we need it for the cdft_gradient)
       mode_counter = mode_counter + 1

       !gibo: add spin penalization  U_s * (Tr[Nup] - Tr[Ndown] - Delta_n_t)
       hubatoms_cdft_2 : do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

          sp = hub%species_number(hat)

          ! Cycle loop if present species is not cDFT-active.
          ! For mixed DFT+U and cDFT+U make sure hubbard_energy contributions
          ! in hub%energy(is,hat) is not lost
          if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) then
             do is = 1, pub_num_spins
                cdft_e = cdft_e + hub%energy(is,hat)
             enddo
             CYCLE hubatoms_cdft_2
          endif

          ! gibo: local copies of CONSTRAINED_DFT SPIN parameters
          u_spin      = hub%cdft_u_spin(hat)
          delta_n_t   = h_species(sp)%cdft_target_spin

          ! hub%population(is,hat) has been calculated in SPIN_2: do-end do loop above
          cdft_penalty = hub%population(1,hat) - hub%population(2,hat) - delta_n_t

          ! distribute equally the spin-penalization between the UP- and DOWN-spins...
          do is = 1, pub_num_spins

             hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP * u_spin * cdft_penalty

             cdft_e = cdft_e + hub%energy(is,hat)

             ! store the gradient
             ! [equally distributed between the UP- and DOWN- spin channels]
             hub%cdft_gradient(is,hat, mode_counter) = 0.5_DP * cdft_penalty

          end do

       end do hubatoms_cdft_2

       !====select CDFT-ATOM-mode and iterate...
    endif CDFT_ATOM
    !====select CDFT-ATOM-mode and iterate...


    !####select CDFT-GROUP-mode and iterate...
    ! GROUP-CHARGE-ACCEPTOR mode
    GROUP_CHARGE_ACCEPTOR: if (pub_cdft_group_charge_acceptor) then

       mode_counter = 1 ! charge-modes as first cDFT-gradient component

       ! gibo: initialise total population of acceptor to zero
       pop_acceptor = 0._DP

       ! group_charge_acceptor: UP_electrons only
       SPLIT_UP_DOWN_acceptor: if (pub_cdft_group_charge_up_only) then

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms_3 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_3

             ! gibo: no need to loop over spin for pub_cdft_group_charge_up_only
             !for group-charge-constrained cDFT runs, keep track of total populations
             if (h_species(sp)%cdft_charge_acceptor) &
                  pop_acceptor = pop_acceptor + hub%population(1,hat)  ! spin-UP only
          end do hubatoms_3

          ! collect and sum acceptor/donor populations
          call comms_reduce('SUM', pop_acceptor)

          ! calculate group-constrained energy
          cdft_penalty = pop_acceptor - pub_cdft_group_charge_acceptor_target
          tmp_cdft_e = pub_cdft_group_charge_acceptor_u * cdft_penalty

          !augment total cDFT_energy
          cdft_e = cdft_e + tmp_cdft_e

          ! split equally e_cdft among cDFT charge-acceptor atoms...
          tmp_cdft_e = tmp_cdft_e/real(par%nat_hub_charge_acceptor,kind=DP)

          ! split equally gradient (cdft_penalty) among cDFT charge-acceptor atoms...
          cdft_penalty = cdft_penalty/real(par%nat_hub_charge_acceptor,kind=DP)

          hubatoms_3_1: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_3_1

             ! single out (charge-acceptor) atoms
             if (h_species(sp)%cdft_charge_acceptor) then
                ! no need to split between spin UP and down
                !...and spin channels (0.5 factor)
                hub%energy(1,hat) =  hub%energy(1,hat) + tmp_cdft_e     ! UP-channel only
                hub%cdft_gradient(1,hat,mode_counter) =  cdft_penalty   ! UP-channel only
             endif

          end do hubatoms_3_1

          ! group_charge_acceptor: DOWN_electrons only
       elseif (pub_cdft_group_charge_down_only) then

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms_4 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_4

             ! gibo: no need to loop over spin for pub_cdft_group_charge_down_only
             !for group-charge-constrained cDFT runs, keep track of total populations
             if (h_species(sp)%cdft_charge_acceptor) &
                  pop_acceptor = pop_acceptor + hub%population(2,hat)  ! spin-DOWN only
          end do hubatoms_4

          ! collect and sum acceptor/donor populations
          call comms_reduce('SUM', pop_acceptor)

          ! calculate group-constrained energy
          cdft_penalty = pop_acceptor - pub_cdft_group_charge_acceptor_target
          tmp_cdft_e = pub_cdft_group_charge_acceptor_u * cdft_penalty

          !augment total cDFT_energy
          cdft_e = cdft_e + tmp_cdft_e

          ! split equally e_cdft among cDFT charge-acceptor atoms...
          tmp_cdft_e = tmp_cdft_e/real(par%nat_hub_charge_acceptor,kind=DP)

          ! split equally gradient (cdft_penalty) among cDFT charge-acceptor atoms...
          cdft_penalty = cdft_penalty/real(par%nat_hub_charge_acceptor,kind=DP)

          hubatoms_4_1: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_4_1

             ! single out (charge-acceptor) atoms
             if (h_species(sp)%cdft_charge_acceptor) then
                ! no need to split between spin UP and down
                !...and spin channels (0.5 factor)
                hub%energy(2,hat) =  hub%energy(2,hat) + tmp_cdft_e     !DOWN-channel only
                hub%cdft_gradient(2,hat,mode_counter) =  cdft_penalty   !DOWN-channel only
             endif

          end do hubatoms_4_1

          ! group_charge_acceptor: UP+DOWN electrons
       else

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms_5 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_5

             ! gibo: loop over spin
             spin_5: do is = 1, pub_num_spins

                !for group-charge-constrained cDFT runs, keep track of total populations
                if (h_species(sp)%cdft_charge_acceptor) &
                     pop_acceptor = pop_acceptor + hub%population(is,hat)

             end do spin_5

          end do hubatoms_5

          ! collect and sum acceptor/donor populations
          call comms_reduce('SUM', pop_acceptor)

          ! calculate group-constrained energy
          cdft_penalty = pop_acceptor - pub_cdft_group_charge_acceptor_target
          tmp_cdft_e = pub_cdft_group_charge_acceptor_u * cdft_penalty

          !augment total cDFT_energy
          cdft_e = cdft_e + tmp_cdft_e

          ! split equally e_cdft among cDFT charge-acceptor atoms...
          tmp_cdft_e = tmp_cdft_e/real(par%nat_hub_charge_acceptor,kind=DP)

          ! split equally gradient (cdft_penalty) among cDFT charge-acceptor atoms...
          cdft_penalty = cdft_penalty/real(par%nat_hub_charge_acceptor,kind=DP)

          hubatoms_5_1: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_5_1

             ! single out (charge-acceptor) atoms
             if (h_species(sp)%cdft_charge_acceptor) then

                !...and spin channels (0.5 factor)
                do is = 1, pub_num_spins
                   hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*tmp_cdft_e
                   hub%cdft_gradient(is,hat,mode_counter) =  0.5_DP*cdft_penalty
                end do

             endif

          end do hubatoms_5_1

       endif SPLIT_UP_DOWN_acceptor

    endif GROUP_CHARGE_ACCEPTOR
    !####select CDFT-GROUP-mode and iterate...


    !####select CDFT-GROUP-mode and iterate...
    ! GROUP-CHARGE-DONOR mode
    GROUP_CHARGE_DONOR: if (pub_cdft_group_charge_donor) then
       !####select CDFT-GROUP-mode and iterate...

       mode_counter = 1 ! charge-modes as first cDFT-gradient component

       ! gibo: initialise total population of acceptor to zero
       pop_donor    = 0._DP

       ! group_charge_donor: UP_electrons only
       SPLIT_UP_DOWN_donor: if (pub_cdft_group_charge_up_only) then

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms_6 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_6

             ! no need to loop over spin
             !for group-charge-constrained cDFT runs, keep track of total populations
             if (h_species(sp)%cdft_charge_donor)    &
                                !pop_donor = pop_donor + hub%population(is,hat)
                  pop_donor = pop_donor + hub%population(1,hat)  ! UP-electrons only
             !pop_donor = pop_donor + hub%population(2,hat)  ! DOWN-electrons only

          end do hubatoms_6

          ! collect and sum acceptor/donor populations
          call comms_reduce('SUM', pop_donor)

          ! calculate group-constrained energy
          cdft_penalty = pop_donor - pub_cdft_group_charge_donor_target
          tmp_cdft_e = pub_cdft_group_charge_donor_u * cdft_penalty

          !augment total cDFT_energy
          cdft_e = cdft_e + tmp_cdft_e

          ! split equally e_cdft among cDFT charge-donor atoms
          tmp_cdft_e = tmp_cdft_e/real(par%nat_hub_charge_donor,kind=DP)

          ! split equally gradient (cdft_penalty) among cDFT charge-donor atoms
          cdft_penalty = cdft_penalty/real(par%nat_hub_charge_donor,kind=DP)

          hubatoms_6_1: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_6_1

             ! single out (charge-donor) atoms
             if (h_species(sp)%cdft_charge_donor) then

                ! no need to split between spin-channels
                !...and spin channels (0.5 factor)
                hub%energy(1,hat) =  hub%energy(1,hat) + tmp_cdft_e     ! UP-electrons only
                hub%cdft_gradient(1,hat,mode_counter) =  cdft_penalty   ! UP-electrons only

             endif

          end do hubatoms_6_1

          ! group_charge_donor: DOWN_electrons only
       elseif (pub_cdft_group_charge_down_only) then

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms_7 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_7

             ! no need to loop over spin

             !for group-charge-constrained cDFT runs, keep track of total populations
             if (h_species(sp)%cdft_charge_donor)    &
                  pop_donor = pop_donor + hub%population(2,hat)  ! DOWN-electrons only

          end do hubatoms_7

          ! collect and sum acceptor/donor populations
          call comms_reduce('SUM', pop_donor)

          ! calculate group-constrained energy
          cdft_penalty = pop_donor - pub_cdft_group_charge_donor_target
          tmp_cdft_e = pub_cdft_group_charge_donor_u * cdft_penalty

          !augment total cDFT_energy
          cdft_e = cdft_e + tmp_cdft_e

          ! split equally e_cdft among cDFT charge-donor atoms
          tmp_cdft_e = tmp_cdft_e/real(par%nat_hub_charge_donor,kind=DP)

          ! split equally gradient (cdft_penalty) among cDFT charge-donor atoms
          cdft_penalty = cdft_penalty/real(par%nat_hub_charge_donor,kind=DP)

          hubatoms_7_1: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_7_1

             ! single out (charge-donor) atoms
             if (h_species(sp)%cdft_charge_donor) then

                ! no need to split between spin-channels
                !...and spin channels (0.5 factor)
                hub%energy(2,hat) =  hub%energy(2,hat) + tmp_cdft_e     ! DOWN-electrons only
                hub%cdft_gradient(2,hat,mode_counter) =  cdft_penalty   ! DOWN-electrons only

             endif

          end do hubatoms_7_1

          ! group_charge_donor: UP+DOWN electrons
       else

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms_8 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_8

             ! gibo: loop over spin
             spin_8: do is = 1, pub_num_spins

                !for group-charge-constrained cDFT runs, keep track of total populations
                if (h_species(sp)%cdft_charge_donor)    &
                     pop_donor = pop_donor + hub%population(is,hat)

             end do spin_8

          end do hubatoms_8

          ! collect and sum acceptor/donor populations
          call comms_reduce('SUM', pop_donor)

          ! calculate group-constrained energy
          cdft_penalty = pop_donor - pub_cdft_group_charge_donor_target
          tmp_cdft_e = pub_cdft_group_charge_donor_u * cdft_penalty

          !augment total cDFT_energy
          cdft_e = cdft_e + tmp_cdft_e

          ! split equally e_cdft among cDFT charge-donor atoms
          tmp_cdft_e = tmp_cdft_e/real(par%nat_hub_charge_donor,kind=DP)

          ! split equally gradient (cdft_penalty) among cDFT charge-donor atoms
          cdft_penalty = cdft_penalty/real(par%nat_hub_charge_donor,kind=DP)

          hubatoms_8_1: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_8_1

             ! single out (charge-donor) atoms
             if (h_species(sp)%cdft_charge_donor) then

                !...and spin channels (0.5 factor)
                do is = 1, pub_num_spins
                   hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*tmp_cdft_e
                   hub%cdft_gradient(is,hat,mode_counter) =  0.5_DP*cdft_penalty
                end do

             endif

          end do hubatoms_8_1

       endif SPLIT_UP_DOWN_donor

       !####select CDFT-GROUP-mode and iterate...
    endif GROUP_CHARGE_DONOR
    !####select CDFT-GROUP-mode and iterate...


    !####select CDFT-GROUP-mode and iterate...
    GROUP_SPIN_ACCEPTOR: if (pub_cdft_group_spin_acceptor) then
       !####select CDFT-GROUP-mode and iterate...

       !for charge- and spin-constrained runs, spin is the second-cDFT mode
       if (pub_cdft_group_charge_acceptor.OR.pub_cdft_group_charge_donor) then
          mode_counter = 2
       else
          mode_counter = 1
       endif

       ! gibo: initialise total population of acceptor and donor (to zero)
       spinmom_acceptor = 0._DP

       ! work out spin-pipulation difference (local spin_moment) for acceptor and donor
       hubatoms_cdft_9 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

          sp = hub%species_number(hat)

          ! Cycle loop if present species is not cDFT-active.
          if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_cdft_9

          if (h_species(sp)%cdft_spin_acceptor)     &
               spinmom_acceptor = spinmom_acceptor + &
               hub%population(1,hat) - hub%population(2,hat)

       end do hubatoms_cdft_9

       ! collect and sum acceptor/donor spin_moments
       call comms_reduce('SUM', spinmom_acceptor)

       ! calculate group-constrained energy
       cdft_penalty = spinmom_acceptor - pub_cdft_group_spin_acceptor_target
       tmp_cdft_e = pub_cdft_group_spin_acceptor_u * cdft_penalty

       !augment total cDFT_energy
       cdft_e = cdft_e + tmp_cdft_e

       ! split equally e_cdft among cDFT spin-acceptor atoms
       tmp_cdft_e = tmp_cdft_e/real(par%nat_hub_spin_acceptor,kind=DP)

       ! split equally gradient (cdft_penalty) among cDFT spin-acceptor atoms
       cdft_penalty = cdft_penalty/real(par%nat_hub_spin_acceptor,kind=DP)

       hubatoms_9_1: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
          sp = hub%species_number(hat)

          ! Cycle loop if present species is not cDFT-active.
          if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_9_1

          ! single out (spin-acceptor) atoms
          if (h_species(sp)%cdft_spin_acceptor)  then

             !...and spin channels (0.5 factor)
             do is = 1, pub_num_spins
                hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*tmp_cdft_e
                hub%cdft_gradient(is,hat,mode_counter) =  0.5_DP*cdft_penalty
             end do

          endif

       end do hubatoms_9_1

       !####select CDFT-GROUP-mode and iterate...
    endif GROUP_SPIN_ACCEPTOR
    !####select CDFT-GROUP-mode and iterate...


    !####select CDFT-GROUP-mode and iterate...
    GROUP_SPIN_DONOR: if (pub_cdft_group_spin_donor) then
       !####select CDFT-GROUP-mode and iterate...

       !for charge- and spin-constrained runs, spin is the second-cDFT mode
       if (pub_cdft_group_charge_acceptor.OR.pub_cdft_group_charge_donor) then
          mode_counter = 2
       else
          mode_counter = 1
       endif

       ! gibo: initialise total population of acceptor and donor (to zero)
       spinmom_donor    = 0._DP

       ! work out spin-pipulation difference (local spin_moment) for acceptor and donor
       hubatoms_cdft_10 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

          sp = hub%species_number(hat)

          ! Cycle loop if present species is not cDFT-active.
          if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_cdft_10

          if (h_species(sp)%cdft_spin_donor)   &
               spinmom_donor = spinmom_donor + &
               hub%population(1,hat) - hub%population(2,hat)

       end do hubatoms_cdft_10

       ! collect and sum acceptor/donor spin_moments
       call comms_reduce('SUM', spinmom_donor)

       ! calculate group-constrained energy
       cdft_penalty = spinmom_donor - pub_cdft_group_spin_donor_target
       tmp_cdft_e = pub_cdft_group_spin_donor_u * cdft_penalty

       !augment total cDFT_energy
       cdft_e = cdft_e + tmp_cdft_e

       ! split equally e_cdft among cDFT spin-donor atoms
       tmp_cdft_e = tmp_cdft_e/real(par%nat_hub_spin_donor,kind=DP)

       ! split equally gradient (cdft_penalty) among cDFT spin-donor atoms
       cdft_penalty = cdft_penalty/real(par%nat_hub_spin_donor,kind=DP)

       hubatoms_10_1: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
          sp = hub%species_number(hat)

          ! Cycle loop if present species is not cDFT-active.
          if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_10_1

          ! single out (spin-donor) atoms
          if (h_species(sp)%cdft_spin_donor)  then

             !...and spin channels (0.5 factor)
             do is = 1, pub_num_spins
                hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*tmp_cdft_e
                hub%cdft_gradient(is,hat,mode_counter) =  0.5_DP*cdft_penalty
             end do

          endif

       end do hubatoms_10_1

       !####select CDFT-GROUP-mode and iterate...
    endif GROUP_SPIN_DONOR
    !####select CDFT-GROUP-mode and iterate...


    !####select CDFT-GROUP-mode and iterate...
    ! gibo: GROUP-CHARGE_DIFFERENCE constrained run
    GROUP_CHARGE_DIFF: if (pub_cdft_group_charge_diff) then
       !####select CDFT-GROUP-mode and iterate...

       mode_counter = 1 ! charge-modes as first cDFT-gradient component

       ! gibo: initialise total population of acceptor and donor (to zero)
       pop_acceptor = 0._DP
       pop_donor    = 0._DP

       ! group_charge_diff: UP_electrons only
       SPLIT_UP_DOWN_diff: if (pub_cdft_group_charge_up_only) then

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms_11 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_11

             ! No need to loop over spins
             !for group-charge-constrained cDFT runs, keep track of total populations
             if (h_species(sp)%cdft_charge_acceptor) &
                  pop_acceptor = pop_acceptor + hub%population(1,hat) ! UP-electrons ONLY

             if (h_species(sp)%cdft_charge_donor)    &
                  pop_donor = pop_donor + hub%population(1,hat)      ! UP-electrons only

          end do hubatoms_11

          ! collect and sum acceptor/donor populations
          call comms_reduce('SUM', pop_acceptor)
          call comms_reduce('SUM', pop_donor)

          ! calculate group-constrained energy
          !gibo: penalty inverted to be consistent with -ve (+ve) potentials at
          !       acceptor (donor) sites [pub_cdft_group_diff_u has been set +ve in
          !       rundat_blocks_mod.F90]
          cdft_penalty = pub_cdft_group_charge_diff_target - pop_acceptor + pop_donor
          tmp_cdft_e = pub_cdft_group_charge_diff_u * cdft_penalty

          !augment total cDFT_energy
          cdft_e = cdft_e + tmp_cdft_e

          ! split equally e_cdft among CDFT_atoms...
          tmp_cdft_e = tmp_cdft_e / & ! mind the division-operation!!!
               real(par%nat_hub_charge_donor+par%nat_hub_charge_acceptor,kind=DP)

          ! split equally gradient (cdft_penalty) among CDFT_atoms...
          cdft_penalty = cdft_penalty / & ! mind the division-operation!!!
               real(par%nat_hub_charge_donor+par%nat_hub_charge_acceptor,kind=DP)

          hubatoms_11_1: do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_11_1

             ! no need to loop over spins
             !...and spin channels (0.5 factor)
             hub%energy(1,hat) =  hub%energy(1,hat) + tmp_cdft_e     ! UP-electrons only
             hub%cdft_gradient(1,hat, mode_counter) =  cdft_penalty  ! UP-electrons only

             ! gibo: for acceptors, invert sign of gradient to be consistent
             !      with (+ve) sign of pub_cdft_group_charge_diff_u
             if (h_species(sp)%cdft_charge_acceptor) &
                  hub%cdft_gradient(1,hat,mode_counter) = &    ! UP-electron only
                  -hub%cdft_gradient(1,hat, mode_counter)

          enddo hubatoms_11_1

          ! group_charge_diff: DOWN_electrons only
       elseif (pub_cdft_group_charge_down_only) then

          hubatoms_12 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_12

             ! No need to loop over spins
             ! gibo: loop over spin
             !for group-charge-constrained cDFT runs, keep track of total populations
             if (h_species(sp)%cdft_charge_acceptor) &
                  pop_acceptor = pop_acceptor + hub%population(2,hat) ! DOWN-electrons ONLY

             if (h_species(sp)%cdft_charge_donor)    &
                  pop_donor = pop_donor + hub%population(2,hat)      ! DOWN-electrons only

          end do hubatoms_12

          ! collect and sum acceptor/donor populations
          call comms_reduce('SUM', pop_acceptor)
          call comms_reduce('SUM', pop_donor)

          ! calculate group-constrained energy
          !gibo: penalty inverted to be consistent with -ve (+ve) potentials at
          !       acceptor (donor) sites [pub_cdft_group_diff_u has been set +ve in
          !       rundat_blocks_mod.F90]
          cdft_penalty = pub_cdft_group_charge_diff_target - pop_acceptor + pop_donor
          tmp_cdft_e = pub_cdft_group_charge_diff_u * cdft_penalty

          !augment total cDFT_energy
          cdft_e = cdft_e + tmp_cdft_e

          ! split equally e_cdft among CDFT_atoms...
          tmp_cdft_e = tmp_cdft_e / & ! mind the division-operation!!!
               real(par%nat_hub_charge_donor+par%nat_hub_charge_acceptor,kind=DP)

          ! split equally gradient (cdft_penalty) among CDFT_atoms...
          cdft_penalty = cdft_penalty / & ! mind the division-operation!!!
               real(par%nat_hub_charge_donor+par%nat_hub_charge_acceptor,kind=DP)

          hubatoms_12_1: do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_12_1

             ! no need to loop over spins
             !...and spin channels (0.5 factor)
             hub%energy(2,hat) =  hub%energy(2,hat) + tmp_cdft_e     ! DOWN-electrons only
             hub%cdft_gradient(2,hat, mode_counter) =  cdft_penalty  ! DOWN-electrons only

             ! gibo: for acceptors, invert sign of gradient to be consistent
             !      with (+ve) sign of pub_cdft_group_charge_diff_u
             if (h_species(sp)%cdft_charge_acceptor) &
                  hub%cdft_gradient(2,hat,mode_counter) = &    ! DOWN-electron only
                  -hub%cdft_gradient(2,hat, mode_counter)

          end do hubatoms_12_1

          ! group_charge_diff: UP+DOWN electrons
       else

          ! ddor: Loop over Hubbard atoms on my proc
          hubatoms_13 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_13

             ! gibo: loop over spin
             spin_13: do is = 1, pub_num_spins

                !for group-charge-constrained cDFT runs, keep track of total populations
                if (h_species(sp)%cdft_charge_acceptor) &
                     pop_acceptor = pop_acceptor + hub%population(is,hat)

                if (h_species(sp)%cdft_charge_donor)    &
                     pop_donor = pop_donor + hub%population(is,hat)

             end do spin_13

          end do hubatoms_13

          ! collect and sum acceptor/donor populations
          call comms_reduce('SUM', pop_acceptor)
          call comms_reduce('SUM', pop_donor)

          ! calculate group-constrained energy
          !gibo: penalty inverted to be consistent with -ve (+ve) potentials at
          !       acceptor (donor) sites [pub_cdft_group_diff_u has been set +ve in
          !       rundat_blocks_mod.F90]
          cdft_penalty = pub_cdft_group_charge_diff_target - pop_acceptor + pop_donor
          tmp_cdft_e = pub_cdft_group_charge_diff_u * cdft_penalty

          !augment total cDFT_energy
          cdft_e = cdft_e + tmp_cdft_e

          ! split equally e_cdft among CDFT_atoms...
          tmp_cdft_e = tmp_cdft_e / & ! mind the division-operation!!!
               real(par%nat_hub_charge_donor+par%nat_hub_charge_acceptor,kind=DP)

          ! split equally gradient (cdft_penalty) among CDFT_atoms...
          cdft_penalty = cdft_penalty / & ! mind the division-operation!!!
               real(par%nat_hub_charge_donor+par%nat_hub_charge_acceptor,kind=DP)

          hubatoms_13_1: do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

             hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
             sp = hub%species_number(hat)

             ! Cycle loop if present species is not cDFT-active.
             if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_13_1

             !...and spin channels (0.5 factor)
             do is = 1, pub_num_spins
                hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*tmp_cdft_e
                hub%cdft_gradient(is,hat, mode_counter) =  0.5_DP*cdft_penalty

                ! gibo: for acceptors, invert sign of gradient to be consistent
                !      with (+ve) sign of pub_cdft_group_charge_diff_u
                if (h_species(sp)%cdft_charge_acceptor) &
                     hub%cdft_gradient(is,hat,mode_counter) = &
                     -hub%cdft_gradient(is,hat, mode_counter)
             end do

          end do hubatoms_13_1

       endif SPLIT_UP_DOWN_diff

       !====select CDFT-mode and iterate...
    endif GROUP_CHARGE_DIFF
    !====select CDFT-mode and iterate...


    !====select CDFT-mode and iterate...
    ! gibo: GROUP-SPIN_DIFFERENCE constrained run
    GROUP_SPIN_DIFF: if (pub_cdft_group_spin_diff) then
       !====select CDFT-mode and iterate...

       !for charge- and spin-constrained runs, spin is the second-cDFT mode
       if (pub_cdft_group_charge_diff) then
          mode_counter = 2
       else
          mode_counter = 1
       endif

       ! gibo: initialise total population of acceptor and donor (to zero)
       spinmom_acceptor = 0._DP
       spinmom_donor    = 0._DP

       ! work out spin-pipulation difference (local spin_moment) for acceptor and donor
       hubatoms_cdft_14 : do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

          sp = hub%species_number(hat)

          ! Cycle loop if present species is not cDFT-active.
          if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_cdft_14

          if (h_species(sp)%cdft_spin_acceptor)     &
               spinmom_acceptor = spinmom_acceptor + &
               hub%population(1,hat) - hub%population(2,hat)

          if (h_species(sp)%cdft_spin_donor)   &
               spinmom_donor = spinmom_donor + &
               hub%population(1,hat) - hub%population(2,hat)

       end do hubatoms_cdft_14

       ! collect and sum acceptor/donor spin_moments
       call comms_reduce('SUM', spinmom_acceptor)
       call comms_reduce('SUM', spinmom_donor)

       ! calculate group-constrained energy
       !gibo: penalty inverted to be consistent with -ve (+ve) potentials at
       !       acceptor (donor) sites [pub_cdft_group_diff_u has been set +ve in
       !       rundat_blocks_mod.F90]
       cdft_penalty = pub_cdft_group_spin_diff_target - spinmom_acceptor + spinmom_donor
       tmp_cdft_e = pub_cdft_group_spin_diff_u * cdft_penalty

       !augment total cDFT_energy
       cdft_e = cdft_e + tmp_cdft_e

       ! split equally e_cdft among CDFT_atoms...
       !cdft_e = cdft_e/real(par%nat_hub,kind=DP)
       tmp_cdft_e = tmp_cdft_e / &  ! MIND THE DIVISION!!!
            real(par%nat_hub_spin_acceptor+par%nat_hub_spin_donor,kind=DP)

       ! split equally gradient (cdft_penalty) among CDFT_atoms...
       cdft_penalty = cdft_penalty / &  ! MIND THE DIVISION!!!
            real(par%nat_hub_spin_acceptor+par%nat_hub_spin_donor,kind=DP)

       hubatoms_14_1: do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
          sp = hub%species_number(hat)

          ! Cycle loop if present species is not cDFT-active.
          if (pub_cdft_hubbard .AND. (.not.h_species(sp)%cdft_active)) CYCLE hubatoms_14_1

          !...and spin channels (0.5 factor)
          do is = 1, pub_num_spins
             hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*tmp_cdft_e
             hub%cdft_gradient(is,hat,mode_counter) =  0.5_DP*cdft_penalty

             !gibo: regardless of the sign of pub_cdft_group_charge_diff_target,
             !      invert sign of gradient for acceptors [alpha, is=1 channel]
             !      to be consistent with (+ve) sign of  pub_cdft_group_diff_u
             ![MIND no need to invert [beta, is=2 channel] as that is done
             ! automatically in sbrtn internal_cdft_update_u [cdft_cg_mod.F90]

             ! invert only UP-gradient for acceptor-atoms
             if ((h_species(sp)%cdft_spin_acceptor).AND.(is==1)) &
                                !hub%cdft_gradient(is,hat) = -hub%cdft_gradient(is,hat)
                  hub%cdft_gradient(is,hat,mode_counter) = -hub%cdft_gradient(is,hat,mode_counter)

             !! invert only DOWN-gradient for donor-atoms
             !if ((h_species(sp)%cdft_spin_donor).AND.(is==2)) &
             !   hub%cdft_gradient(is,hat,mode_counter) = -hub%cdft_gradient(is,hat,mode_counter)
          end do

       end do hubatoms_14_1

       !====select CDFT-mode and iterate...
    endif GROUP_SPIN_DIFF
    !====select CDFT-mode and iterate...

    ! gibo: for ATOM-CHARGE/SPIN-constrained simulation,
    ! 1. augment hubbard_e with cdft_e; 2. reduce across procs
    ! update hubbard_e and write out the cDFT(+U) occupancy matrices and energies
    ATOMS: if (pub_cdft_atom_charge .OR. pub_cdft_atom_spin) then

       ! cDFT+U: copy cDFT-augmented Hubbard energy into hubbard_e
       ! [hub%energy(:,:)-> cdft_e already contains hubbard_contributions]
       if (pub_cdft_hubbard) then
          hubbard_e = cdft_e

          ! cDFT-only: copy cDFT-augmented Hubbard energy into hubbard_e
          ! [hubbard_e and hub%energy(:,:) = 0 for cdft_hubbard=.FALSE.
       else
          hubbard_e = hubbard_e + cdft_e

       endif

       ! gibo: MIND that at this point hubbard_e contains also cdft_e (see above)
       call comms_reduce('SUM', hubbard_e )


       ! if running GROUP_CHARGE/SPIN_* simulation behave differently:
       ! [cdft_e is de facto already reduced across all procs by use of
       ! "call comms_reduce('SUM',hub%population)"->tmp_cdft_e-> cdft_e ]
       ! 1. reduce hubbard_e across procs; 2. augment hubbard_e with cdft_e
    else

       ! collect and sum hubbard_e
       call comms_reduce('SUM', hubbard_e )

       hubbard_e = hubbard_e + cdft_e

    end if ATOMS


    ! Write out the cDFT(+U) occupancy matrices and energies=== START
    call comms_reduce('SUM',hub%energy)!
    call comms_reduce('SUM',hub%population)

    !if (pub_on_root .and. pub_output_detail == VERBOSE) &
    if (pub_output_detail >= VERBOSE) &
         call cdft_energy_info(hub,hub_proj_basis)
    ! Write out the cDFT(+U) occupancy matrices and energies=== END


    ! gibo: for cdft_group_charge/spin[_difference] print out
    ! info on total population/spin-moment [difference]
    TOT_POP: if ((.not.pub_cdft_atom_charge) .AND. (.not.pub_cdft_atom_spin)) then
       ROOT: if (pub_on_root) then

          write(stdout,'(/a)') repeat('#',80)

          if (pub_cdft_group_charge_up_only) then
             write(stdout,'(a)') 'MIND: UP-populations only '
          elseif (pub_cdft_group_charge_down_only) then
             write(stdout,'(a)') 'MIND: DOWN-populations only '
          else
             write(stdout,'(a)') 'MIND: UP+DOWN-populations '
          endif

          if (pub_cdft_group_charge_acceptor) then
             write(stdout,'(2(a,f13.8))') 'Total ACCEPTOR population:            ',&
                  &pop_acceptor,&
                  &' e    / Target: ', pub_cdft_group_charge_acceptor_target
          endif

          if (pub_cdft_group_charge_donor) then
             write(stdout,'(2(a,f13.8))') 'Total DONOR population:               ',&
                  &pop_donor,&
                  &' e    / Target: ', pub_cdft_group_charge_donor_target
          endif

          if (pub_cdft_group_charge_acceptor .AND. &
               pub_cdft_group_charge_donor) then
             write(stdout,'(2(a,f13.8))') 'ACCEPTOR-DONOR population difference: ',&
                  &pop_acceptor - pop_donor,&
                  &' e    / Target: ', pub_cdft_group_charge_acceptor_target &
                  & -pub_cdft_group_charge_donor_target
          endif

          if (pub_cdft_group_spin_acceptor) then
             write(stdout,'(2(a,f13.8))') 'Total ACCEPTOR magnetic-moment:       ',&
                  &spinmom_acceptor,&
                  &' mu_B / Target: ', pub_cdft_group_spin_acceptor_target
          endif

          if (pub_cdft_group_spin_donor) then
             write(stdout,'(2(a,f13.8))') 'Total DONOR magnetic-moment:          ',&
                  &spinmom_donor,&
                  &' mu_B / Target: ', pub_cdft_group_spin_donor_target
          endif

          if (pub_cdft_group_spin_acceptor .AND. &
               pub_cdft_group_spin_donor) then
             write(stdout,'(2(a,f13.8))') 'ACCEPTOR-DONOR magnetic-moment diff.: ',&
                  &spinmom_acceptor - spinmom_donor,&
                  &' e    / Target: ', pub_cdft_group_spin_acceptor_target &
                  & -pub_cdft_group_spin_donor_target
          endif

          if (pub_cdft_group_charge_diff) then
             write(stdout,'(a,f13.8,a)')    'Total ACCEPTOR population:            ',&
                  &pop_acceptor, ' e'
             write(stdout,'(a,f13.8,a)')    'Total DONOR population:               ',&
                  &pop_donor,    ' e'
             write(stdout,'(2(a,f13.8))') 'ACCEPTOR-DONOR population diff.:      ',&
                  &pop_acceptor - pop_donor,&
                  &' e    / Target: ', pub_cdft_group_charge_diff_target
          endif

          if (pub_cdft_group_spin_diff) then
             write(stdout,'(a,f13.8,a)')    'Total ACCEPTOR magnetic-moment:       ',&
                  &spinmom_acceptor,' mu_B'
             write(stdout,'(a,f13.8,a)')    'Total DONOR magnetic-moment:          ',&
                  &spinmom_donor,   ' mu_B'
             write(stdout,'(2(a,f13.8))') 'ACCEPTOR-DONOR magnetic-moment diff.: ',&
                  &spinmom_acceptor - spinmom_donor,&
                  &' mu_B / Target: ', pub_cdft_group_spin_diff_target
          endif


          write(stdout,'(a/)') repeat('#',80)

       end if ROOT
    end if TOT_POP

    ! gibo: destroy local SPAM3 matrixes
    call sparse_destroy(site_buffer)

    call timer_clock('cdft_energy_total',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &cdft_energy_total'

  end subroutine cdft_energy_total


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cdft_energy_info(hub,hub_proj_basis)

    !==========================================================!
    ! Writes out occupancies and energies of Hubbard sites.    !
    !----------------------------------------------------------!
    ! Written by David O'Regan in November 2009                !
    ! Adapted to cDFT by gibo in November 2011                 !
    !==========================================================!

    use comms, only: comms_bcast, pub_on_root, pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_cdft_print_all_occ, pub_num_spins
    use sparse, only: SPAM3, sparse_get_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(in) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: species, this_species_atom ! ddor: Loop over species
    integer :: jj,hat,project,spin,theatom, proc
    integer :: ierr
    real(kind=DP) :: local_population, local_moment, local_energy
    real(kind=DP) :: local_population_up, local_population_down
    real(kind=DP),allocatable :: occupancy_block(:,:,:)
    character(10) :: fmt, tmp

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &cdft_energy_info'

    allocate(occupancy_block(hub_proj_basis%max_on_atom, &
         hub_proj_basis%max_on_atom,pub_num_spins),stat=ierr)
    call utils_alloc_check('hubbard_energy_info','occupancy_block',ierr)

    do species=1,par%num_hub_species

       this_species_atom = 0

       atomloop: do hat=1,par%nat_hub

          if (hub%species_number(hat) == species) then

             this_species_atom = this_species_atom + 1

             ! Unless cdft_print_all_occ=T,
             ! print occupancies for only the first 2-atoms of each cDFT-species
             if ( (this_species_atom .gt. 2) .AND. &
                  (.not. pub_cdft_print_all_occ) ) exit atomloop

             theatom = par%distr_atom(hub%orig(hat))
             project = hub_proj_basis%num_on_atom(theatom)
             proc = par%proc_of_atom(hub%orig(hat))

             if (proc==pub_my_proc_id) then
                do spin=1,pub_num_spins
                   call sparse_get_block(occupancy_block(:,:,spin), &
                        hub%occupancy_matrix(spin),theatom,theatom)
                end do
             end if
             call comms_bcast(proc,occupancy_block)

             if (.not.pub_on_root) cycle

             write(stdout,'(/a)') repeat('#',80)


             ! ddor: TELL ME ABOUT THIS ATOM
             write(stdout,'(a,i6,a,a)') 'cDFT information on &
                  &atom ',this_species_atom,' of constrained species ',&
                  h_species(species)%hub_species

             if (pub_num_spins == 2) then

                do spin=1, pub_num_spins

                   write(stdout,'(a)') repeat('#',80)

                   ! ddor: WRITE OUT INFORMATION ON OCCUPANCY MATRIX
                   write(stdout,'(a,2(i6,a))') 'Occupancy matrix of &
                        &cDFT site ', hat,' and spin ',spin,' is '
                   write(tmp,'(i6)') project
                   write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f12.8)'
                   do jj=1,project
                      write(stdout,fmt) occupancy_block(jj,1:project,spin)
                   end do

                   if (MAXVAL(occupancy_block(:,:,spin)) .gt. 1.0_DP) then
                      write(stdout,'(a,2(i6,a))') 'WARNING: OCCUPANCY MATRIX &
                           &of cDFT site ',&
                           hat,' and spin ',spin,' exceeeds 1.'
                   end if

                end do

                write(stdout,'(a)') repeat('#',80)

                ! ddor: WRITE OUT INFORMATION ON OCCUPANCY OF SITE
                local_population_up   =  hub%population(1,hat)
                local_population_down =  hub%population(2,hat)

                local_population = hub%population(1,hat) &
                     + hub%population(2,hat)
                local_moment = hub%population(1,hat) &
                     - hub%population(2,hat)
                local_energy = hub%energy(1,hat) &
                     + hub%energy(2,hat)
                write(stdout,'(a,i6,a,f12.8,a)') 'Total occupancy &
                     &of cDFT site       ',  hat,' is ', &
                     &local_population, ' e'
                write(stdout,'(a,i6,a,f12.8,a)') 'UP occupancy &
                     &of cDFT site          ',  hat,' is ', &
                     &local_population_up, ' e'
                write(stdout,'(a,i6,a,f12.8,a)') 'DOWN occupancy &
                     &of cDFT site        ',  hat,' is ', &
                     &local_population_down, ' e'
                write(stdout,'(a,i6,a,f12.8,a)') 'Local magnetic moment &
                     &of cDFT site ',  hat,' is ', local_moment,' mu_B'
                ! ddor: WRITE OUT INFORMATION ON DFT+U ENERGY
                write(stdout,'(a,i6,a,f12.8,a)') 'cDFT energy of cDFT &
                     &site           ', hat,' is ', local_energy,' Ha'

                write(stdout,'(a/)') repeat('#',80)

             end if

          end if

       end do atomloop

    end do

    deallocate(occupancy_block,stat=ierr)
    call utils_dealloc_check('hubbard_energy_info','occupancy_block',ierr)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &cdft_energy_info'

  end subroutine cdft_energy_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dft_nu_energy_total(hub,hubbard_e, &
       pur_denskern, hub_proj_basis, hub_overlap, hub_overlap_t)

    !==========================================================!
    ! Writes out occupancies and energies of Hubbard sites.    !
    !----------------------------------------------------------!
    ! This module is a modified version of hubbard_energy_info.!
    ! Written by Glenn Moynihan in June/July 2014.             !
    !                                                          !
    !==========================================================!

    use comms, only: pub_my_proc_id, comms_barrier, comms_reduce
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_output_detail, &
         pub_cdft_multi_proj, &
         pub_hubbard_unify_sites, pub_num_spins, pub_imag_thr
    use sparse, only: SPAM3, sparse_axpy, sparse_copy, &
         sparse_product, sparse_create, sparse_destroy, &
         sparse_get_element, sparse_trace, sparse_show_matrix, &
         sparse_scale, sparse_take_real_part
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
        dense_put_element, dense_normal_eigensolve, dense_write

    implicit none

    !Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(inout) :: hubbard_e  ! Hubbard DFT+U energy
    type(SPAM3), intent(in) :: pur_denskern(pub_num_spins)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(SPAM3), intent(in) :: hub_overlap, hub_overlap_t

    ! Local Variables
    type(SPAM3) :: site_buffer
    type(SPAM3) :: energy_buffer
    type(SPAM3) :: energy_buffer2
    integer :: channels
    integer :: hat_on_proc, hat, sp, theatom
    integer :: is
    integer :: i,ihub
    real(kind=DP) :: occupancyelement
    integer :: natoms

    real(kind=DP) :: u1_up, u1_down, u2_up, u2_down
    real(kind=DP) :: nu_up, nu_down, target_n_up, target_n_down
    real(kind=DP) :: nu_up_2, nu_down_2

    integer      :: jj,kk,num_projs_in_group, ierr, group
    integer      :: fe,le
    type(DEM)    :: dense_occ_matrix
    type(DEM)    :: dense_occ_eigenvecs
    real(kind=DP), allocatable :: dense_occ_eigenvals(:)
    ! agrecocmplx
    logical      :: loc_cmplx
    type(SPAM3)  :: occupancy_matrix_cmplx

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering dft_nu_energy_total'

    call timer_clock('dft_nu_energy_total',1)

    call sparse_create(site_buffer,hub_overlap_t,pur_denskern(1))
    ! ddor: Blank occupancy matrix and products thereof to
    !       block-diagonal matrices in order to prevent
    !       Hubbard sites with null overlap interfering

    !gom: Initialise gradient to 0.0_DP
    hub%nu_gradient(:,:,:) = 0.0_DP

    ! agrecocmplx
    loc_cmplx = pur_denskern(1)%iscmplx

    if (loc_cmplx) then
       call sparse_create(occupancy_matrix_cmplx, hub%occupancy_matrix(1), &
            iscmplx=loc_cmplx)
    end if

    !gom: If in single sites loop over atomic sites
    if (.not.pub_hubbard_unify_sites) then

       energy_buffer%structure   = 'G'
       energy_buffer2%structure  = 'G'

       call sparse_create(energy_buffer)
       call sparse_create(energy_buffer2)

       spin_du: do is=1,pub_num_spins

           !gom create occupancy matrix
           ! WK
           call sparse_product(site_buffer,hub_overlap_t,&
                  pur_denskern(is))
           ! WKV blanked to 'G'
           ! agrecocmplx: occupancy_matrix is real anyway, so use
           ! occupancy_matrix_cmplx when using complex NGWFs, and then
           ! convert to real matrix
           if (loc_cmplx) then
              call sparse_product(occupancy_matrix_cmplx, site_buffer, &
                   hub_overlap)
              ! convert complex to real
              ! agrecocmplx: use safe conversion from complex to real
              !call sparse_copy(hub%occupancy_matrix(is), occupancy_matrix_cmplx)
              call sparse_take_real_part(hub%occupancy_matrix(is), occupancy_matrix_cmplx, &
                   pub_imag_thr)
           ! standard real case
           else
              call sparse_product(hub%occupancy_matrix(is),site_buffer,&
                   hub_overlap)
           end if

           !gom: for non-spin-polarised, double the occupancy matrix
           if (pub_num_spins==1) then
               call sparse_scale(hub%occupancy_matrix(is),2.0_DP)
           endif

           ! WKV.WKV blanked to 'G'
           call sparse_product(energy_buffer2,&
                  hub%occupancy_matrix(is),hub%occupancy_matrix(is))
           call sparse_copy(energy_buffer,hub%occupancy_matrix(is))


           ! ddor: Loop over Hubbard atoms on my proc
           natoms = par%num_hub_atoms_on_proc(pub_my_proc_id)

           hubatoms: do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

              hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
              theatom = par%distr_atom(hub%orig(hat))
              channels = hub_proj_basis%num_on_atom(theatom)
              sp = hub%species_number(hat)


              if  ((.not.(pub_cdft_multi_proj)) .AND. (channels .ne.&
                   (2 * h_species(sp)%hub_ang_mom + 1)))  then
                 call utils_abort('Error in hubbard_energy_total: &
                      &hub_proj_basis%num_on_atom and hub_ang_mom &
                      &mismatch')
              end if

              u1_up         = hub%nu_u1_up(hat)
              u1_down       = hub%nu_u1_down(hat)
              u2_up         = hub%nu_u2_up(hat)
              u2_down       = hub%nu_u2_down(hat)
              nu_up         = hub%nu_nu_up(hat)
              nu_down       = hub%nu_nu_down(hat)
              nu_up_2       = nu_up*nu_up
              nu_down_2     = nu_down*nu_down
              target_n_up   = hub%nu_target_n_up(hat)
              target_n_down = hub%nu_target_n_down(hat)


              do i=1,channels
                 ihub = hub_proj_basis%first_on_atom(theatom) + i - 1

                 ! Tr[n]
                 call sparse_get_element(occupancyelement, &
                      energy_buffer,ihub,ihub)
                 hub%population(is,hat) = hub%population(is,hat) &
                      + occupancyelement

                 ! Tr[n^2]
                 call sparse_get_element(occupancyelement, &
                      energy_buffer2,ihub,ihub)
                 hub%population_2(is,hat) = hub%population_2(is,hat) &
                      + occupancyelement
              end do

              ! add U1, U2 terms for each spin
              if (is == 1) then

                  hub%nu_gradient(is,hat,1) = &
                    0.5_DP * (hub%population(is,hat)-&
                              (target_n_up + nu_up))

                  hub%nu_gradient(is,hat,2) = &
                    0.5_DP * ((target_n_up +nu_up_2)-&
                              hub%population_2(is,hat))

                  hub%energy(is,hat) = &
                    u1_up * hub%nu_gradient(is,hat,1) + &
                    u2_up * hub%nu_gradient(is,hat,2)

              elseif (is == 2) then

                  hub%nu_gradient(is,hat,1) = &
                    0.5_DP * (hub%population(is,hat)-&
                              (target_n_down + nu_down))

                  hub%nu_gradient(is,hat,2) = &
                    0.5_DP * ((target_n_down+nu_down_2)-&
                              hub%population_2(is,hat))

                  hub%energy(is,hat) = &
                    u1_down * hub%nu_gradient(is,hat,1) + &
                    u2_down * hub%nu_gradient(is,hat,2)

             endif

             hubbard_e = hubbard_e + hub%energy(is,hat)

           end do hubatoms

       end do spin_du

       call sparse_destroy(energy_buffer2)
       call sparse_destroy(energy_buffer)
       call sparse_destroy(site_buffer)

   !gom: If unified sites is on, perform calculation over the groups
   else if(pub_hubbard_unify_sites) then

       !gom: Initialise 'group' population
       hub%group_pop(:,:)   = 0.0_DP
       hub%group_pop_2(:,:) = 0.0_DP

       spin_us: do is=1,pub_num_spins

           !gom create occupancy matrix
           ! WK
           call sparse_product(site_buffer,hub_overlap_t,&
                  pur_denskern(is))
           ! WKV blanked to 'G'
           ! agrecocmplx: convert complex matrix to real if using complex
           ! NGWFs
           if (loc_cmplx) then
              call sparse_product(occupancy_matrix_cmplx, site_buffer, &
                   hub_overlap)
              ! convert complex to real
              ! agrecocmplx: safe conversion of complex to real
              !call sparse_copy(hub%occupancy_matrix(is), occupancy_matrix_cmplx)
              call sparse_take_real_part(hub%occupancy_matrix(is), occupancy_matrix_cmplx, &
                   pub_imag_thr)
           ! standard real case
           else
              call sparse_product(hub%occupancy_matrix(is),site_buffer,&
                   hub_overlap)
           end if

           !gom: in unpolarised case, double occupancy matrix
           if (pub_num_spins==1) then
               call sparse_scale(hub%occupancy_matrix(is),2.0_DP)
           endif

           !gom start with the first element of occupancy matrix
           fe = 1

           !gom: Can parallelise over the groups?
           grouploop: do group=1,hub%num_groups

              !gom: create dense group matrices
              num_projs_in_group = hub%nu_projectors(group)
              allocate(dense_occ_eigenvals(num_projs_in_group),&
                   stat=ierr)
              call utils_alloc_check('dft_nu_energy_info',&
                   'dense_occ_eigenvals',ierr)

              call dense_create(dense_occ_matrix,&
                   num_projs_in_group,num_projs_in_group)
              call dense_create(dense_occ_eigenvecs,&
                   num_projs_in_group,num_projs_in_group)

              le = fe + num_projs_in_group - 1

              !gom: make the dense occupancy matrix
              do jj=fe,le
                  do kk=fe,le
                      call sparse_get_element(occupancyelement,&
                           hub%occupancy_matrix(is),jj,kk)
                      call dense_put_element(occupancyelement,&
                           dense_occ_matrix,jj-fe+1,kk-fe+1)
                  end do
              end do

              fe = fe + num_projs_in_group

              !gom: solve for the e'vectors and e'values of the projectors
              call dense_normal_eigensolve(num_projs_in_group,&
                   dense_occ_eigenvals,dense_occ_matrix, &
                   dense_occ_eigenvecs)

              !gom: Here my "Hubbard population" must be the population
              !of the group, i.e. the trace of the total occupancy matrix
              do kk=1,num_projs_in_group
                 hub%group_pop(is,group)= &
                      hub%group_pop(is,group)+dense_occ_eigenvals(kk)
                 hub%group_pop_2(is,group)=&
                      hub%group_pop_2(is,group)+&
                      dense_occ_eigenvals(kk)*dense_occ_eigenvals(kk)
              end do

              u1_up         = hub%nu_u1_up(group)
              u1_down       = hub%nu_u1_down(group)
              u2_up         = hub%nu_u2_up(group)
              u2_down       = hub%nu_u2_down(group)
              nu_up         = hub%nu_nu_up(group)
              nu_down       = hub%nu_nu_down(group)
              nu_up_2       = nu_up*nu_up
              nu_down_2     = nu_down*nu_down
              target_n_up   = hub%nu_target_n_up(group)
              target_n_down = hub%nu_target_n_down(group)

              if (is == 1) then

                  hub%nu_gradient(is,group,1) = &
                    0.5_DP * (hub%group_pop(is,group)-&
                              (target_n_up+nu_up))

                  hub%nu_gradient(is,group,2) = &
                    0.5_DP * ((target_n_up+nu_up_2)-&
                              hub%group_pop_2(is,group))

                  hub%energy(is,group) =&
                    u1_up * hub%nu_gradient(is,group,1) + &
                    u2_up * hub%nu_gradient(is,group,2)

              elseif (is==2) then

                  hub%nu_gradient(is,group,1) = &
                    0.5_DP * (hub%group_pop(is,group)-&
                              (target_n_down+nu_down))

                  hub%nu_gradient(is,group,2) = &
                    0.5_DP * ((target_n_down+nu_down_2)-&
                              hub%group_pop_2(is,group))

                  hub%energy(is,group) = &
                    u1_down * hub%nu_gradient(is,group,1) + &
                    u2_down * hub%nu_gradient(is,group,2)

             endif

             hubbard_e = hubbard_e + hub%energy(is,group)

             call dense_destroy(dense_occ_eigenvecs)
             call dense_destroy(dense_occ_matrix)
             deallocate(dense_occ_eigenvals,stat=ierr)
             call utils_dealloc_check('dft_nu_energy_info',&
                  'dense_occ_eigenvals',ierr)

          end do grouploop

        end do spin_us

    endif


    ! gibo_to_done: make sure this is compatible with cDFT should we aim
    ! for cDFT[cDFT_DU] runs...
    call comms_reduce('SUM',hub%energy)
    call comms_reduce('SUM',hub%population)
    call comms_reduce('SUM',hub%population_2)
    call comms_reduce('SUM',hub%group_pop)
    call comms_reduce('SUM',hub%group_pop_2)

    ! agrecocmplx
    if (loc_cmplx) then
       call sparse_destroy(occupancy_matrix_cmplx)
    end if

    call sparse_destroy(site_buffer)

    ! gom:  print DFT+nu info
    call comms_reduce('SUM', hubbard_e )

    if (pub_output_detail >= VERBOSE) &
        call dft_nu_energy_info(hub,hub_proj_basis)

    call timer_clock('dft_nu_energy_total',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &dft_nu_energy_total'

  end subroutine dft_nu_energy_total


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dft_nu_energy_info(hub,hub_proj_basis)

    !==========================================================!
    ! Writes out occupancies and energies of Hubbard sites.    !
    !----------------------------------------------------------!
    ! Version of hubbard_energy info edited by Glenn Moynihan  !
    ! in July 2014 to use in the DFT+nu mode.                  !
    !==========================================================!

    use comms, only: comms_bcast, pub_on_root, &
         pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_num_spins, pub_hubbard_unify_sites
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_get_block, sparse_show_matrix, &
         sparse_create, sparse_destroy, &
         sparse_get_element
    use utils, only: utils_alloc_check, utils_dealloc_check
    use timer, only: timer_clock
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
        dense_put_element, dense_normal_eigensolve, dense_write
    use constants, only: HARTREE_IN_EVS


    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(in) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: species, this_species_atom ! ddor: Loop over species
    integer :: jj,hat,project,spin,theatom, proc
    integer :: ierr, sp
    real(kind=DP) :: local_population, local_moment, local_energy
    real(kind=DP) :: u1, u2, N, nu           !gom : local variables
    real(kind=DP),allocatable :: occupancy_block(:,:,:)
    character(10) :: fmt, tmp

    integer      :: fe,le     ! first/last element to get from occupancy matrix
    integer      :: kk, num_projs_in_group, group
    real(kind=DP):: occupancyelement
    type(DEM)    :: dense_occ_matrix
    type(DEM)    :: dense_occ_eigenvecs
    real(kind=DP), allocatable :: dense_occ_eigenvals(:)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
       &dft_nu_energy_info'

    call timer_clock('dft_nu_energy_info',1)

    allocate(occupancy_block(hub_proj_basis%max_on_atom, &
             hub_proj_basis%max_on_atom,pub_num_spins),stat=ierr)
    call utils_alloc_check('dft_nu_energy_info','occupancy_block',ierr)

    !gom: Print Occupancies for Single Sites by looping over each
    !spin of each atom of each species
    if (.not. pub_hubbard_unify_sites) then

        speciesloop: do species=1,par%num_hub_species

           this_species_atom = 0

           atomloop: do hat=1,par%nat_hub
              if (hub%species_number(hat) == species) then
                 this_species_atom = this_species_atom + 1

                 if (this_species_atom .gt. 2) exit atomloop

                 theatom = par%distr_atom(hub%orig(hat))
                 project = hub_proj_basis%num_on_atom(theatom)
                 proc = par%proc_of_atom(hub%orig(hat))

                 if (proc==pub_my_proc_id) then
                    do spin=1,pub_num_spins
                       call sparse_get_block(occupancy_block(:,:,spin), &
                            hub%occupancy_matrix(spin),theatom,theatom)
                    end do
                 end if
                 call comms_bcast(proc,occupancy_block)

                 if (.not.pub_on_root) cycle

                 write(stdout,'(/a)') repeat('#',80)

                 ! ddor: TELL ME ABOUT THIS ATOM
                 write(stdout,'(a,i6,a,a)') ' DFT+nu information on atom ', &
                      this_species_atom,' of Hubbard species ',&
                      h_species(species)%hub_species

                 do spin=1,pub_num_spins

                    write(stdout,'(a)') repeat('#',80)

                    if (pub_num_spins == 2) then
                       ! ddor: WRITE OUT INFORMATION ON OCCUPANCY MATRIX
                       write(stdout,'(a,2(i6,a))') ' Occupancy matrix of &
                            &Hubbard site ',hat,' and spin ',spin,' is '
                    endif

                    ! ddor: write header line
                    write(stdout,'(a)',advance='no') '  m_l ='
                    do jj=-(project-1)/2,(project-1)/2
                       write(stdout,'(i5,7x)',advance='no') jj
                    end do
                    write(stdout,'(1x)')

                    write(tmp,'(i6)') project
                    write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f12.8)'

                    do jj=1,project
                       write(stdout,fmt) occupancy_block(jj,1:project,spin)
                    end do

                   ! gom : print Hubbard constraint info
                   write(stdout,'(a)') repeat('-',80)
                   write(stdout,'(a,2(i1,a))') ' Constraint information of &
                        &Hubbard site ',hat,' and spin ',spin,' is '

                    ! U1 U2 nu
                    if (spin==1) then
                       u1=hub%nu_u1_up(hat)
                       u2=hub%nu_u2_up(hat)
                       nu=hub%nu_nu_up(hat)
                       N =hub%nu_target_n_up(hat)
                    else
                       u1=hub%nu_u1_down(hat)
                       u2=hub%nu_u2_down(hat)
                       nu=hub%nu_nu_down(hat)
                       N =hub%nu_target_n_down(hat)
                   endif

                    write(stdout,'(a,2(i1,a),f12.8)') ' H(',hat,',',spin,')     &
                        & =    ',hub%energy(spin,hat) * HARTREE_IN_EVS
                    write(stdout,'(a,2(i1,a),f12.8)') ' U1(',hat,',',spin,')    &
                        & =    ',u1 * HARTREE_IN_EVS
                    write(stdout,'(a,2(i1,a),f12.8)') ' U2(',hat,',',spin,')    &
                        & =    ',u2 * HARTREE_IN_EVS
                    write(stdout,'(a,2(i1,a),f12.8)') ' N(',hat,',',spin,')     &
                        & =    ',N
                    write(stdout,'(a,2(i1,a),f12.8)') ' nu(',hat,',',spin,')    &
                        & =    ',nu
                    write(stdout,'(a,2(i1,a),f12.8)') ' Tr(',hat,',',spin,')    &
                        & =    ',hub%population(spin,hat)
                    write(stdout,'(a,2(i1,a),f12.8)') ' Tr2(',hat,',',spin,')   &
                        & =    ',hub%population_2(spin,hat)
                    write(stdout,'(a,2(i1,a),f12.8)') ' dE/dU1(',hat,',',spin,')&
                        & =    ',hub%nu_gradient(spin,hat,1)
                    write(stdout,'(a,2(i1,a),f12.8)') ' dE/dU2(',hat,',',spin,')&
                        & =    ',hub%nu_gradient(spin,hat,2)

                    if (MAXVAL(occupancy_block(:,:,spin)) .gt. 1.0_DP) then
                       write(stdout,'(a,2(i6,a))') 'WARNING: OCCUPANCY MATRIX &
                            &of Hubbard site ',&
                            hat,' and spin ',spin,' exceeds 1.'
                    end if

                 end do
                write(stdout,'(a)') repeat('#',80)

                 ! ddor: WRITE OUT INFORMATION ON OCCUPANCY OF SITE
                 local_population = sum(hub%population(:,hat))
                 local_energy = sum(hub%energy(:,hat))

                 write(stdout,'(a,i6,a,f12.8,a)') ' Total occupancy &
                      &of Hubbard site ', hat,' is       ', &
                      &local_population, ' e'

                 if (pub_num_spins==2) then
                    local_moment =  &
                         &hub%population(1,hat) - hub%population(2,hat)
                    write(stdout,'(a,i6,a,f12.8,a)') ' Local magnetic moment &
                         &of Hubbard site ', hat,' is ', local_moment,' mu_B'
                 end if

                 ! ddor: WRITE OUT INFORMATION ON DFT+U ENERGY
                 write(stdout,'(a,i6,a,f12.8,a)') ' DFT+U energy of Hubbard &
                      &site ',hat,' is          ', local_energy,' Ha'

                 write(stdout,'(a/)') repeat('#',80)
              end if
           end do atomloop
        end do speciesloop

    !gom: Print Occupancies for Unified Sites
    else if(pub_hubbard_unify_sites) then

        !gom: start at the first element
        fe = 1

        grouploop: do group=1,hub%num_groups

            !gom: create dense group matrices: occupancy matrix,
            num_projs_in_group = hub%nu_projectors(group)

            allocate(dense_occ_eigenvals(num_projs_in_group),stat=ierr)
            call utils_alloc_check('dft_nu_energy_info',&
                    'dense_occ_eigenvals',ierr)
            call dense_create(dense_occ_matrix,&
                     num_projs_in_group,num_projs_in_group)
            call dense_create(dense_occ_eigenvecs,&
                    num_projs_in_group,num_projs_in_group)

            le = fe + num_projs_in_group - 1

            write(stdout,'(/a)') repeat('#',80)

            spinloop: do spin=1,pub_num_spins

                write(stdout,'(a,i2,a,i1)') ' DFT+nu information of Group ', &
                     group,' Spin ',spin

                write(stdout,'(a)') repeat('#',80)

                write(stdout,'(a)') ' Species in Group:'
                !gom: print species in group
                do sp=1,par%num_hub_species
                    if (group == h_species(sp)%nu_group) then
                        write(stdout,'(1x,a)',advance='no') h_species(sp)%hub_species
                    end if
                end do
               write(stdout,'(a)')' '

               write(stdout,'(a)') repeat('-',80)

                !gom: Copy the projectors in the group we are interested
                !in from occupancy_matrix to dense_occ_matrix
                do jj=fe,le
                    do kk=fe,le
                        call sparse_get_element(occupancyelement,&
                                hub%occupancy_matrix(spin),jj,kk)
                        call dense_put_element(occupancyelement,&
                                dense_occ_matrix,jj-fe+1,kk-fe+1)
                    end do
                end do

                if(spin==1) then
                    call dense_write(dense_occ_matrix,"occupancy_matrix_1")
                else
                    call dense_write(dense_occ_matrix,"occupancy_matrix_2")
                endif

                !gom: solve for the eigenvectors and eigenvalues of the
                !group in that spin channel
                call dense_normal_eigensolve(num_projs_in_group,&
                        dense_occ_eigenvals,dense_occ_matrix,&
                dense_occ_eigenvecs)


                !gom: print eigenvalues as list
                do kk=1,num_projs_in_group
                    write(stdout,'(a,3(i2,a),f12.8)')' Occupancy &
                         &Eigenvalue(',group,',',spin,',',kk,') = ',&
                         dense_occ_eigenvals(kk)
                end do
               write(stdout,'(a)') repeat('-',80)

                if (spin==1) then
                   u1=hub%nu_u1_up(group)
                   u2=hub%nu_u2_up(group)
                   nu=hub%nu_nu_up(group)
                   N =hub%nu_target_n_up(group)
                else
                   u1=hub%nu_u1_down(group)
                   u2=hub%nu_u2_down(group)
                   nu=hub%nu_nu_down(group)
                   N =hub%nu_target_n_down(group)
                endif

                write(stdout,'(a,2(i1,a),f12.8)') ' H(',group,',',spin,')     &
                        & =    ',hub%energy(spin,group) * HARTREE_IN_EVS
                write(stdout,'(a,2(i1,a),f12.8)') ' U1(',group,',',spin,')    &
                        & =    ',u1 * HARTREE_IN_EVS
                write(stdout,'(a,2(i1,a),f12.8)') ' U2(',group,',',spin,')    &
                        & =    ',u2 * HARTREE_IN_EVS
                write(stdout,'(a,2(i1,a),f12.8)') ' N(',group,',',spin,')     &
                        & =    ',N
                write(stdout,'(a,2(i1,a),f12.8)') ' nu(',group,',',spin,')    &
                        & =    ',nu
                write(stdout,'(a,2(i1,a),f12.8)') ' Tr(',group,',',spin,')    &
                        & =    ',hub%group_pop(spin,group)
                write(stdout,'(a,2(i1,a),f12.8)') ' Tr2(',group,',',spin,')   &
                        & =    ',hub%group_pop_2(spin,group)
                write(stdout,'(a,2(i1,a),f12.8)') ' dE/dU1(',group,',',spin,')&
                        & =    ',hub%nu_gradient(spin,group,1)
                write(stdout,'(a,2(i1,a),f12.8)') ' dE/dU2(',group,',',spin,')&
                        & =    ',hub%nu_gradient(spin,group,2)

                if (MAXVAL(dense_occ_eigenvals(:)) .gt. 1.0_DP &
                    .and. pub_num_spins==2) then
                   write(stdout,'(a,2(i1,a))') 'WARNING: Eigenvalue of Hubbard group ',&
                        group,', spin ',spin,' exceeds 1.'
                end if

                write(stdout,'(/a)') repeat('#',80)

                !gom: print eigenvectors to file
                if(spin==1) then
                    call dense_write(dense_occ_eigenvecs,"eigenvectors_1")
                else if(spin==2) then
                    call dense_write(dense_occ_eigenvecs,"eigenvectors_2")
                end if

            end do spinloop

            fe = fe + num_projs_in_group

            local_population = sum(hub%group_pop(:,group))
            local_energy = sum(hub%energy(:,group))

            write(stdout,'(a,i6,a,f12.8,a)') 'Total occupancy &
                 &of Hubbard group ', group,' is       ', &
                 &local_population, ' e'

            if (pub_num_spins==2) then
               local_moment =  &
                    &hub%population(1,group) - hub%population(2,group)
               write(stdout,'(a,i6,a,f12.8,a)') 'Local magnetic moment &
                    &of Hubbard group ', group,' is ', local_moment,' mu_B'
            end if

            ! ddor: WRITE OUT INFORMATION ON DFT+U ENERGY
            write(stdout,'(a,i6,a,f12.8,a)') 'DFT+U energy of Hubbard &
                 &group ',group,' is          ', local_energy,' Ha'

            write(stdout,'(a/)') repeat('#',80)

            !gom: deallocate memory
            deallocate(dense_occ_eigenvals,stat=ierr)
            call utils_dealloc_check('hubbard_energy_info',&
                    'dense_occ_eigenvals',ierr)
            call dense_destroy(dense_occ_eigenvecs)
            call dense_destroy(dense_occ_matrix)

        end do grouploop


    endif

    deallocate(occupancy_block,stat=ierr)
    call utils_dealloc_check('hubbard_energy_info','occupancy_block',ierr)

    call timer_clock('dft_nu_energy_info',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &dft_nu_energy_info'

  end  subroutine dft_nu_energy_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_ngwf_select(hub,ngwf_basis,hub_proj_basis, &
       hub_consist_matrix,hub_consist_tmatrix)

    !==========================================================!
    ! This subroutine determines which NGWFs have the greatest !
    ! projection onto a set of localised hydrogenic orbitals   !
    ! which are chosen by the user.                            !
    !----------------------------------------------------------!
    ! Written by David O'Regan in April 2009.                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_my_proc_id, comms_barrier, &
         comms_reduce, comms_bcast, pub_on_root
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use ion, only: ELEMENT
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_output_detail, pub_cdft_multi_proj ,pub_dmft_order_proj, &
         pub_dmft_switch_off_proj_order
    use sparse, only: SPAM3, sparse_get_element
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    ! agrecocmplx: add dependence on COEF
    use datatypes, only: COEF

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(SPAM3), intent(in) :: hub_consist_matrix
    type(SPAM3), intent(in) :: hub_consist_tmatrix

    ! Local Variables
    integer, allocatable, dimension(:) :: ll, mm
    integer :: ierr
    integer :: ijk
    integer :: hat_on_proc, hat, theatom
    integer :: olaps, projs
    integer :: ngwf_we_are_on, hub_proj_we_are_on
    integer :: ngwf_count, hub_proj_count
    integer :: kk, jj
    !real(kind=DP) :: overlap_element, transpose_overlap_element
    ! agrecocmplx: use COEF type instead
    type(COEF) :: overlap_element, transpose_overlap_element
    real(kind=DP), allocatable :: ngwfsquareovlp(:,:,:)
    real(kind=DP), allocatable :: ngwfsquareoverlapsum(:)
    ! agrecocmplx
    logical :: iscmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_ngwf_select'

    call timer_clock('hubbard_ngwf_select',1)

    ! Allocate workspace
    allocate(ll(1),stat=ierr)
    call utils_alloc_check('hubbard_ngwf_select','ll',ierr)
    allocate(mm(1),stat=ierr)
    call utils_alloc_check('hubbard_ngwf_select','mm',ierr)
    allocate(ngwfsquareovlp(hub_proj_basis%max_on_atom,ngwf_basis%max_on_atom, &
         par%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_ngwf_select','ngwfsquareovlp',ierr)
    allocate(ngwfsquareoverlapsum(ngwf_basis%max_on_atom),stat=ierr)
    call utils_alloc_check('hubbard_ngwf_select','ngwfsquareoverlapsum',ierr)

    ngwfsquareovlp(:,:,:) = 0.0_DP

    ! agrecocmplx
    iscmplx = hub_consist_matrix%iscmplx
    overlap_element%iscmplx = iscmplx
    transpose_overlap_element%iscmplx = iscmplx

    ! dhpt: Initialise localised NGWFs with zero
    do hat = 1, par%nat_hub
       hub%localisedngwfs(:,hat) = 0
    end do

    ! ddor: Loop over Hubbard atoms on my proc
    do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

       hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

       ! Atom number of Hubbard atom
       theatom = par%distr_atom(hub%orig(hat))
       ! ddor: Number of NGWFs on this Hubbard atom
       olaps = ngwf_basis%num_on_atom(theatom)
       ! ddor: Number of Hubbard projectors on this Hubbard atom
       projs = hub_proj_basis%num_on_atom(theatom)

       !==============================================================
       ! LOCATE THE MOST LOCALISED NGWFS
       !==============================================================

       ! ddor: Initially set output to zero
       hub%localisedngwfs(:,hat) = 0
       ngwfsquareoverlapsum = 0.0_DP

       ! ddor: First do a loop over local row ngwfs
       !       in the hydrogenic-ngwf overlap matrix
       !       to see which NGWFs \psi on the Hubbard atom
       !       maximise the
       !       sum_{localised atomic projectors \phi} | < \phi | \psi > |**2

       ! ddor: loop over the NGWFs on this Hubbard atom
       do ngwf_count = 1, olaps

          ngwf_we_are_on = ngwf_count + ngwf_basis%first_on_atom(theatom) - 1

          do hub_proj_count = 1, projs

             hub_proj_we_are_on = hub_proj_count + &
                  hub_proj_basis%first_on_atom(theatom) - 1

             ! ddor: Add up the square of the overlap between
             !       the current NGWF and
             !       the contravariant projectors.
             !       As well as correcting for nonorthogonality,
             !       this fixes the normalisation of the NGWF projectors.
             ! agrecocmplx
             if (iscmplx) then
                call sparse_get_element(overlap_element%z,hub_consist_matrix,&
                     ngwf_we_are_on, hub_proj_we_are_on)
                call sparse_get_element(transpose_overlap_element%z,&
                     hub_consist_tmatrix,hub_proj_we_are_on, ngwf_we_are_on)
             else
                call sparse_get_element(overlap_element%d,hub_consist_matrix,&
                     ngwf_we_are_on, hub_proj_we_are_on)
                call sparse_get_element(transpose_overlap_element%d,&
                     hub_consist_tmatrix,hub_proj_we_are_on, ngwf_we_are_on)
             end if

             ! ddor: The elements on the trace of the
             !       \sum_i <NGWF_|\phi_i><\phi_i|NGWF^> matrix
             ! agrecocmplx
             if (iscmplx) then
                ngwfsquareovlp(hub_proj_count,ngwf_count,hat) = &
                     real(overlap_element%z * transpose_overlap_element%z, kind=DP)
             else
                ngwfsquareovlp(hub_proj_count,ngwf_count,hat) = &
                     overlap_element%d * transpose_overlap_element%d
             end if

             ngwfsquareoverlapsum(ngwf_count) = &
                  ngwfsquareoverlapsum(ngwf_count) + &
                  ngwfsquareovlp(hub_proj_count,ngwf_count,hat)

          end do

       end do

       ! ddor: Find out which NGWFs produce the greatest sum of
       !       square overlaps with the atomic projectors
       projectortest: do kk = 1, projs

          ! ddor: find the location of the maximum overlap**2
          ll = MAXLOC(ngwfsquareoverlapsum(1:olaps) - (/( pub_dmft_order_proj*real(ijk,kind=DP),ijk=1,olaps )/) )
          ! ddor: take note of the NGWF number of the projector
          hub%localisedngwfs(kk,hat) = ll(1) + &
               ngwf_basis%first_on_atom(theatom) - 1

          ! ddor: destroy this overlap
          ngwfsquareoverlapsum(ll(1)) = 0.0_DP

       end do projectortest

       ! ddor: Assuming that projectors are always grouped together,
       !       re-sort them to preserve magnetic quantum number
       mm(1) = MINVAL(hub%localisedngwfs(1:projs,hat))

       if(.not.pub_dmft_switch_off_proj_order)then
          projectorsort: do jj = 1, projs

             hub%localisedngwfs(jj,hat) = mm(1) + jj - 1

          end do projectorsort
       end if

    end do ! ddor: End loop over Hubbard atoms on my proc

    ! ddor: Reduce, Allow other procs to catch up
    call comms_reduce('SUM',ngwfsquareovlp)
    ! dhpt: Synchronise localisedngwfs as well
    call comms_reduce('SUM',hub%localisedngwfs)

    call comms_barrier
    ! ddor: Write out the | < \phi | \psi > |**2 matrix
    if (pub_on_root .and. pub_output_detail >= VERBOSE) &
         &call internal_ngwf_character_info

    ! Deallocate temporary storage
    deallocate(ngwfsquareoverlapsum,stat=ierr)
    call utils_dealloc_check('hubbard_ngwf_select','ngwfsquareoverlapsum',ierr)
    deallocate(ngwfsquareovlp,stat=ierr)
    call utils_dealloc_check('hubbard_ngwf_select','ngwfsquareovlp',ierr)
    deallocate(mm,stat=ierr)
    call utils_dealloc_check('hubbard_ngwf_select','mm',ierr)
    deallocate(ll,stat=ierr)
    call utils_dealloc_check('hubbard_ngwf_select','ll',ierr)

    call timer_clock('hubbard_ngwf_select',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_ngwf_select'

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_ngwf_character_info

      implicit none

      integer :: species, this_species_atom ! ddor: Loop over species
      integer :: project, jj
      integer :: kk, curr_proj
      character(10) :: fmt, tmp

      do species=1,par%num_hub_species

         this_species_atom = 0

         atomloop : do hat=1,par%nat_hub

            if (hub%species_number(hat) == species) then

               this_species_atom = this_species_atom + 1
               theatom = par%distr_atom(hub%orig(hat))
               if (this_species_atom .gt. 2) exit atomloop

               write(stdout,'(/a)') repeat('#',80)

               ! ddor: TELL ME ABOUT THIS ATOM
               write(stdout,'(a,i6,a,a)') 'DFT+U information on &
                    &atom ',this_species_atom,' of Hubbard species ',&
                    h_species(species)%hub_species

               ! ddor: Number of NGWFs on this Hubbard atom
               olaps = ngwf_basis%num_on_atom(theatom)

               ! ddor: Number of Hubbard projectors on this Hubbard atom
               project = hub_proj_basis%num_on_atom(theatom)

               write(stdout,'(a)') repeat('#',80)

               ! ddor: WRITE OUT INFORMATION ON | < \phi | \psi > |**2 MATRIX
               write(stdout,'(a)') 'Using the on-atom contravariant metric, &
                    &the NGWF hydrogenic character'
               write(stdout,'(a,i6,a)') '| < \phi | \psi > |**2 &
                    &matrix of Hubbard site ',hat,' is '


               ! ndmh: write header line
               if (.not.pub_cdft_multi_proj) then
                  write(stdout,'(a)',advance='no') '  m_l ='
                  do jj=-(project-1)/2,(project-1)/2
                     write(stdout,'(i5,7x)',advance='no') jj
                  end do
                  write(stdout,'(a)') '  SUM'

                  ! gibo: keep the m_l(l) labelling also for cdft_multi_proj
               elseif (pub_cdft_multi_proj) then
                  write(stdout,'(a)',advance='no') 'm_l(l)='

                  do kk = 1,  h_species(species)%num_ang_mom_shells
                     curr_proj =  2*h_species(species)%orb_ang_mom(kk) + 1
                     do jj=-(curr_proj-1)/2,(curr_proj-1)/2
                        write(stdout,'(i5,6x)',advance='no') jj
                     end do
                  enddo
                  write(stdout,'(a)') '  SUM'

               endif

               ! ndmh: write format string
               write(tmp,'(i6)') project+1
               write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f12.8)'

               ! ndmh: write matrix
               do jj=1,olaps
                  write(stdout,fmt) ngwfsquareovlp(1:project,jj,hat), &
                       SUM(ngwfsquareovlp(1:project,jj,hat))
               end do

               write(stdout,'(a,f12.8)') 'Total hydrogenic projection of &
                    &NGWFs on this atom           ', &
                    &SUM(ngwfsquareovlp(1:project,1:olaps,hat))

               write(stdout,'(a/)') repeat('#',80)

            end if

         end do atomloop

      end do

    end subroutine internal_ngwf_character_info

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine hubbard_ngwf_select


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_tensorial_correction(hub,ngwf_basis,proj_basis, &
       consistency,inv_overlap,loc_o_matrix)

    !==========================================================!
    ! This subroutine calculates for a given set of NGWFs the  !
    ! appropriate tensorial correction matrix, which is stored !
    ! in hubbard_o_matrix or hubbard_consist_o_matrix.         !
    !----------------------------------------------------------!
    ! Written by David O'Regan in April 2009.                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    ! Modified by David Turban in March 2014 to include        !
    ! pub_hub_tensor_corr = 4, ie projector duals              !
    ! over set of NGWFs in each CDFT region                    !
    !==========================================================!

    use bibliography, only: bibliography_cite
    use comms, only: pub_my_proc_id, pub_on_root, comms_barrier, &
         comms_reduce
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species ! dhpt: To check whether atoms belong to same hubbard species
    use linalg, only: linalg_invert_sym_matrix, linalg_invert_sym_cmatrix
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_hub_tensor_corr, pub_maxit_hotelling, &
         pub_max_resid_hotelling ! dhpt: need for matrix inversion
    use sparse, only: SPAM3, sparse_get_element, sparse_put_element, &
         sparse_scale, sparse_hotelling_init, sparse_hotelling_invert, &
         sparse_create, sparse_copy, sparse_destroy,sparse_num_rows
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    use dense, only: DEM, dense_create, dense_destroy, dense_convert,&
         dense_invert
    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    logical, intent(in) :: consistency
    type(SPAM3), intent(in) :: inv_overlap
    type(SPAM3), intent(inout) :: loc_o_matrix

    ! Local Variables
    integer :: hat_on_proc, hat, theatom
    integer :: hat_global, hat_global_min, hat_global_max, theatom2
    integer :: sp, sp2
    integer :: iii, jjj ! Used for tensorial correction
    integer :: num_on_atom_row, num_on_atom_col
    integer :: ierr
    real(kind=DP), allocatable :: o_matrix(:,:)
    ! agrecocmplx
    complex(kind=DP), allocatable :: o_matrix_cmplx(:,:)
    logical :: iscmplx
    type(SPAM3) :: o_matrix_tmp ! buffer loc_o_matrix for hotelling
    type (DEM)  :: o_matrix_dns


    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_tensorial_correction'

    call timer_clock('hubbard_tensorial_correction',1)

    call bibliography_cite('SUBSPACE')

    ! agrecocmplx
    iscmplx = hub%consistency_overlap%iscmplx

    if (consistency) then
       num_on_atom_row = ngwf_basis%max_on_atom
    else
       num_on_atom_row = proj_basis%max_on_atom
    end if

    allocate(o_matrix(num_on_atom_row,num_on_atom_row),stat=ierr)
    call utils_alloc_check('hubbard_tensorial_correction','o_matrix',ierr)

    ! agrecocmplx
    if (iscmplx) then
       allocate(o_matrix_cmplx(num_on_atom_row,num_on_atom_row),stat=ierr)
       call utils_alloc_check('hubbard_tensorial_correction','o_matrix_cmplx',ierr)
    end if

    ! ddor: Loop over Hubbard atoms on my proc
    do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

       hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)

       theatom = par%distr_atom(hub%orig(hat)) !Atom number of Hubbard atom

       ! dhpt: Second loop over all Hubbard atoms if pub_hub_tensor_corr==4
       if ((pub_hub_tensor_corr == 4 .or. pub_hub_tensor_corr == 5) .and. (.not. consistency)) then
          hat_global_min = 1
          hat_global_max = par%nat_hub ! total number of Hubbard atoms, length of hub%orig?
       else
          hat_global_min = hat
          hat_global_max = hat
       end if

       do hat_global = hat_global_min, hat_global_max

          if ((pub_hub_tensor_corr == 4 .or. pub_hub_tensor_corr == 5)&
              .and. (.not. consistency)) then

             theatom2 = par%distr_atom(hub%orig(hat_global)) !Atom number of Hubbard atom 2
          else
             theatom2 = theatom
          end if

          ! dhpt: skip rest of second loop if atom and atom2 belong to different CDFT species
          sp  = hub%species_number(hat)
          sp2 = hub%species_number(hat_global)

          if (.not.((h_species(sp)%cdft_charge_donor  .eqv.h_species(sp2)%cdft_charge_donor)&
               .and.(h_species(sp)%cdft_charge_acceptor.eqv.h_species(sp2)%cdft_charge_acceptor)&
               .and.(h_species(sp)%cdft_spin_donor     .eqv.h_species(sp2)%cdft_spin_donor)&
               .and.(h_species(sp)%cdft_spin_acceptor  .eqv.h_species(sp2)%cdft_spin_acceptor))) then

             ! dhpt: do not cycle if global duals are selected
             if (.not.(pub_hub_tensor_corr == 5)) then
                cycle
             endif
          end if

          !==============================================================
          ! CALCULATE THE CONTRAVARIANT METRIC ON EACH HUBBARD SITE
          !==============================================================

          ! ddor: If using nonorthogonal Hubbard projectors,
          !       calculate the necessary tensorial correction matrix
          !       tensorial: if (tensor) then

          o_matrix(:,:) = 0.0_DP

          ! agrecocmplx
          if (iscmplx) then
             o_matrix_cmplx(:,:) = 0.0_DP
          end if

          ! dhpt: Matrices stored row-wise, hence take row index from theatom2
          ! dhpt: which runs over all hubbard atoms
          if (consistency) then
             num_on_atom_row = ngwf_basis%num_on_atom(theatom2)
             num_on_atom_col = ngwf_basis%num_on_atom(theatom)
          else
             num_on_atom_row = proj_basis%num_on_atom(theatom2)
             num_on_atom_col = proj_basis%num_on_atom(theatom)
          end if

          ! ddor: Loop over column Hubbard projectors
          tensor_col_ngwf_loop1: do jjj = 1, num_on_atom_col

             ! ddor: Loop over row Hubbard projectors
             tensor_row_ngwf_loop1: do iii = 1,  num_on_atom_row

                ! ddor: The tensorial correction for all the NGWFs on an atom
                if (consistency) then

                   ! ddor: ONLY RETRIEVE ELEMENTS FROM hub%consistency_overlap WHEN
                   !       FORMING PROJECTOR DUALS IN THE CORRELATED SUBSPACE
                   if ((pub_hub_tensor_corr == 1) .or. (pub_hub_tensor_corr == 4) &
                       .or. (pub_hub_tensor_corr == 5)) then

                      ! agrecocmplx
                      if (iscmplx) then
                         call sparse_get_element(&
                              o_matrix_cmplx(iii,jjj), hub%consistency_overlap,&
                              iii + ngwf_basis%first_on_atom(theatom) - 1,&
                              jjj + ngwf_basis%first_on_atom(theatom) - 1)
                      else
                         call sparse_get_element(&
                              o_matrix(iii,jjj), hub%consistency_overlap,&
                              iii + ngwf_basis%first_on_atom(theatom) - 1,&
                              jjj + ngwf_basis%first_on_atom(theatom) - 1)
                      end if

                      ! ddor: WHEN FORMING PROJECTOR DUALS IN THE SIMULATION CELL,
                      !       USE FULL OVERLAP INVERSE FOR THE TENSOR CORRECTION.
                   else if (pub_hub_tensor_corr == 2) then

                      ! agrecocmplx
                      if (iscmplx) then
                         call sparse_get_element(&
                              o_matrix_cmplx(iii,jjj), inv_overlap,&
                              iii + ngwf_basis%first_on_atom(theatom) - 1,&
                              jjj + ngwf_basis%first_on_atom(theatom) - 1)
                      else
                         call sparse_get_element(&
                              o_matrix(iii,jjj), inv_overlap,&
                              iii + ngwf_basis%first_on_atom(theatom) - 1,&
                              jjj + ngwf_basis%first_on_atom(theatom) - 1)
                      end if

                      ! ddor: WHEN NOT PERFORMING TENSORIAL CORRECTION,
                      !       JUST USE A KRONECKER DELTA FUNCTION,
                      !       WITH A SCALING FACTOR WHICH ASSUMES ORTHOGONALITY
                      !       (TOTAL NUMBER OF CORRELATED ELECTRONS NOT PRESERVED).
                   else if (pub_hub_tensor_corr == 3) then

                      if (iii .eq. jjj) then

                         ! agrecocmplx
                         if (iscmplx) then
                            call sparse_get_element(&
                               o_matrix_cmplx(iii,jjj), hub%consistency_overlap,&
                               iii + ngwf_basis%first_on_atom(theatom) - 1,&
                               jjj + ngwf_basis%first_on_atom(theatom) - 1)

                            o_matrix_cmplx(iii,jjj) = o_matrix_cmplx(iii,jjj)**(-1.0_DP)
                         else
                            call sparse_get_element(&
                                 o_matrix(iii,jjj), hub%consistency_overlap,&
                                 iii + ngwf_basis%first_on_atom(theatom) - 1,&
                                 jjj + ngwf_basis%first_on_atom(theatom) - 1)

                            o_matrix(iii,jjj) = o_matrix(iii,jjj)**(-1.0_DP)
                         end if

                      else
                         ! agrecocmplx
                         if (iscmplx) then
                            o_matrix_cmplx(iii,jjj) = 0.0_DP
                         else
                            o_matrix(iii,jjj) = 0.0_DP
                         end if
                      end if

                   else if (pub_hub_tensor_corr .gt. 5) then
                      call utils_abort('Invalid value for pub_hub_tensor_corr.')
                   end if

                else

                   ! ddor: ONLY RETRIEVE ELEMENTS FROM hub%consistency_overlap WHEN
                   !       FORMING PROJECTOR DUALS IN THE CORRELATED SUBSPACE
                   if (pub_hub_tensor_corr == 1) then

                      ! agrecocmplx
                      if (iscmplx) then
                         call sparse_get_element(&
                              o_matrix_cmplx(iii,jjj), hub%consistency_overlap,&
                              hub%localisedngwfs(iii,hat),&
                              hub%localisedngwfs(jjj,hat))
                      else
                         call sparse_get_element(&
                              o_matrix(iii,jjj), hub%consistency_overlap,&
                              hub%localisedngwfs(iii,hat),&
                              hub%localisedngwfs(jjj,hat))
                      end if

                      ! ddor: WHEN FORMING PROJECTOR DUALS IN THE SIMULATION CELL,
                      !       USE FULL OVERLAP INVERSE FOR THE TENSOR CORRECTION.
                   else if (pub_hub_tensor_corr == 2) then

                      ! agrecocmplx
                      if (iscmplx) then
                         call sparse_get_element(&
                              o_matrix_cmplx(iii,jjj), inv_overlap,&
                              hub%localisedngwfs(iii,hat),&
                              hub%localisedngwfs(jjj,hat))
                      else
                         call sparse_get_element(&
                              o_matrix(iii,jjj), inv_overlap,&
                              hub%localisedngwfs(iii,hat),&
                              hub%localisedngwfs(jjj,hat))
                      end if

                      ! ddor: WHEN NOT PERFORMING TENSORIAL CORRECTION,
                      !       JUST USE A KRONECKER DELTA FUNCTION.
                   else if (pub_hub_tensor_corr == 3) then

                      if (iii .eq. jjj) then

                         ! agrecocmplx
                         if (iscmplx) then
                            call sparse_get_element(&
                                 o_matrix_cmplx(iii,jjj), hub%consistency_overlap,&
                                 hub%localisedngwfs(iii,hat),&
                                 hub%localisedngwfs(jjj,hat))

                            o_matrix_cmplx(iii,jjj) = o_matrix_cmplx(iii,jjj)**(-1.0_DP)
                         else
                            call sparse_get_element(&
                                 o_matrix(iii,jjj), hub%consistency_overlap,&
                                 hub%localisedngwfs(iii,hat),&
                                 hub%localisedngwfs(jjj,hat))

                            o_matrix(iii,jjj) = o_matrix(iii,jjj)**(-1.0_DP)
                         end if

                      else
                         ! agrecocmplx
                         if (iscmplx) then
                            o_matrix_cmplx(iii,jjj) = 0.0_DP
                         else
                            o_matrix(iii,jjj) = 0.0_DP
                         end if

                      end if

                   else if (pub_hub_tensor_corr == 4 .or. pub_hub_tensor_corr == 5) then

                      ! agrecocmplx
                      if (iscmplx) then
                         call sparse_get_element(&
                              o_matrix_cmplx(iii,jjj), hub%consistency_overlap,&
                              hub%localisedngwfs(iii,hat_global),&
                              hub%localisedngwfs(jjj,hat))
                      else
                         call sparse_get_element(&
                              o_matrix(iii,jjj), hub%consistency_overlap,&
                              hub%localisedngwfs(iii,hat_global),&
                              hub%localisedngwfs(jjj,hat))
                      end if

                   else if (pub_hub_tensor_corr .gt. 5) then
                      call utils_abort('Invalid value for pub_hub_tensor_corr.')
                   end if

                end if

             end do tensor_row_ngwf_loop1

          end do tensor_col_ngwf_loop1

          !write(stdout,*) 'Covariant metric on correlated space on site ',&
          !     hat,' is '
          !do mm=1,num_on_atom_col
          !   write(*,*) hub%o_matrix(mm,1:num_on_atom_row)
          !end do

          ! ddor: ONLY INVERT PROJECTOR OVERLAP MATRIX WHEN
          !       FORMING PROJECTOR DUALS IN THE CORRELATED SUBSPACE
          if (pub_hub_tensor_corr == 1) then
             ! ddor: If we want the tensorial correction for all
             !       the NGWFs on an atom
             if (consistency) then
                ! ddor: Perform matrix inverse to find the
                !       contravariant metric on the correlated subspace
                ! agrecocmplx
                if (iscmplx) then
                   call linalg_invert_sym_cmatrix(o_matrix_cmplx(:,:),&
                        ngwf_basis%num_on_atom(theatom) )
                else
                   call linalg_invert_sym_matrix(o_matrix(:,:),&
                        ngwf_basis%num_on_atom(theatom) )
                end if

             else

                ! ddor: Perform matrix inverse to find the
                !       contravariant metric on the correlated subspace
                ! agrecocmplx
                if (iscmplx) then
                   call linalg_invert_sym_cmatrix(o_matrix_cmplx(:,:),&
                        proj_basis%num_on_atom(theatom) )
                else
                   call linalg_invert_sym_matrix(o_matrix(:,:),&
                        proj_basis%num_on_atom(theatom) )
                end if

             end if
          end if

          if ((pub_hub_tensor_corr == 4 .or. pub_hub_tensor_corr ==5)&
              .and. consistency) then
                ! agrecocmplx
                if (iscmplx) then
                   call linalg_invert_sym_cmatrix(o_matrix_cmplx(:,:),&
                        ngwf_basis%num_on_atom(theatom) )
                else
                   call linalg_invert_sym_matrix(o_matrix(:,:),&
                        ngwf_basis%num_on_atom(theatom) )
                end if

          endif


          ! ddor: Loop over column Hubbard projectors
          tensor_col_ngwf_loop2: do jjj = 1,  num_on_atom_col

             ! ddor: Loop over row Hubbard projectors
             tensor_row_ngwf_loop2: do iii = 1, num_on_atom_row

                ! ddor: If we want the tensorial correction for all
                !       the NGWFs on an atom
                if (consistency) then

                   ! ddor: Get the corresponding element of old S
                   !       and put it in hub_consist_o_matrix
                   ! agrecocmplx
                   if (iscmplx) then
                      call sparse_put_element(&
                           o_matrix_cmplx(iii,jjj), loc_o_matrix,&
                           iii + ngwf_basis%first_on_atom(theatom2) - 1,&
                           jjj + ngwf_basis%first_on_atom(theatom) - 1)
                   else
                      call sparse_put_element(&
                           o_matrix(iii,jjj), loc_o_matrix,&
                           iii + ngwf_basis%first_on_atom(theatom2) - 1,&
                           jjj + ngwf_basis%first_on_atom(theatom) - 1)
                   end if

                else

                   ! ddor: Get the corresponding element of old S
                   !       and put it in hub%o_matrix
                   ! agrecocmplx
                   if (iscmplx) then
                      call sparse_put_element(&
                           o_matrix_cmplx(iii,jjj), loc_o_matrix,&
                           iii + proj_basis%first_on_atom(theatom2) - 1,&
                           jjj + proj_basis%first_on_atom(theatom) - 1)
                   else
                      call sparse_put_element(&
                           o_matrix(iii,jjj), loc_o_matrix,&
                           iii + proj_basis%first_on_atom(theatom2) - 1,&
                           jjj + proj_basis%first_on_atom(theatom) - 1)
                   end if

                end if

             end do tensor_row_ngwf_loop2

          end do tensor_col_ngwf_loop2

          !write(stdout,*) 'Contravariant metric on correlated space on site ',&
          !     hat,' is '
          !do mm=1,num_on_atom_col
          !   write(*,*) o_matrix(mm,1:num_on_atom_row)
          !end do

       end do ! dhpt: End of loop over all Hubbard atoms (relevant if hub_tensor_corr==4)
    end do ! ddor: End loop over Hubbard atoms on my proc

    ! dhpt: invert loc_o_matrix for tensorial correction 4
    if ((pub_hub_tensor_corr == 4 .or. pub_hub_tensor_corr == 5)&
        .and. (.not. consistency)) then

       call sparse_create(o_matrix_tmp, loc_o_matrix)
       call sparse_copy(o_matrix_tmp, loc_o_matrix)

       if (pub_maxit_hotelling > 0) then
           call sparse_hotelling_init(loc_o_matrix,o_matrix_tmp)

           ! dhpt: Invert hub_consist_o_matrix (using sparse_hotelling_invert)
           if (pub_on_root) then
              write(stdout,'(a)')'======= Calculation of HUBBARD PROJECTOR S^-1 using &
                   &Hotelling algorithm ======== '
              write(stdout,'(a)')'In hubbard_tensorial_correction.'
              write(stdout,'(a)')'         Iteration  |  Frobenius norm and Maximum &
                   &abs value    (of I-S*S_n^-1)   '
           end if

           call sparse_hotelling_invert(loc_o_matrix, & !inout
                o_matrix_tmp,show_output=.true., &
                max_resid_converged=pub_max_resid_hotelling, &
                num_iter=pub_maxit_hotelling)

           ! dhpt: Should check whether matrix inversion is successful!

           !call sparse_create(mim,hub_consist_o_matrix,overlap_tmp)
           !call sparse_product(mim,hub_consist_o_matrix,overlap_tmp)
           !call sparse_scale(mim,-1.0_DP,1.0_DP)
           !call sparse_create(mim2,mim,hub_consist_o_matrix)
           !call sparse_product(mim2,mim,hub_consist_o_matrix)
           !frob_norm = sparse_rms_element(mim2) * &
           !     sqrt(sparse_num_element(mim2))
           !call sparse_destroy(mim)
           !call sparse_destroy(mim2)
           !if (pub_on_root) then
           !   write(stdout,'(a,e16.8)') 'Frobenius norm of (I-S*S^-1)*S^-1: ',&
           !        &frob_norm
           !endif

            if (pub_on_root) then
               write(stdout,'(a)')'===================================&
                    &============================================='
            end if

            call sparse_destroy(o_matrix_tmp)

        else if (pub_maxit_hotelling == 0) then
            num_on_atom_row = sparse_num_rows(loc_o_matrix)

            ! agrecocmplx
            call dense_create(o_matrix_dns,num_on_atom_row,num_on_atom_row, &
                 iscmplx=iscmplx)
            call dense_convert(o_matrix_dns, loc_o_matrix)
            call dense_invert(o_matrix_dns)
            call dense_convert(loc_o_matrix,o_matrix_dns)
            call dense_destroy(o_matrix_dns)

        else
            call utils_abort('Error in hamiltonian_dens_indep_matrices: &
                 &negative maxit_hotelling ',pub_maxit_hotelling)
        endif

    endif

    ! agrecocmplx
    if (iscmplx) then
       deallocate(o_matrix_cmplx,stat=ierr)
       call utils_dealloc_check('hubbard_tensorial_correction','o_matrix_cmplx',ierr)
    end if

    deallocate(o_matrix,stat=ierr)
    call utils_dealloc_check('hubbard_tensorial_correction','o_matrix',ierr)

    call timer_clock('hubbard_tensorial_correction',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_tensorial_correction'

  end subroutine hubbard_tensorial_correction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_ngwf_on_grid(hub,ngwf_basis,hub_proj_basis,cell)

    !==========================================================!
    ! This subroutine copies the NGWFs which have been chosen  !
    ! as Hubbard projectors into an array which is indexed in  !
    ! the same order as how the <NGWF|projector> matrix        !
    ! elements are sorted.                                     !
    !----------------------------------------------------------!
    ! Written by David O'Regan in September 2009.              !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use datatypes, only: data_set_to_zero
    use comms, only: pub_my_proc_id
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: par=>pub_par
    use simulation_cell, only: CELL_INFO
    use sparse, only : SPAM3, sparse_create, sparse_destroy
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(CELL_INFO), intent(in) :: cell

    ! Local Variables
    integer :: ngwf_count
    integer :: ngwf_offset
    integer :: hubbard_offset
    integer :: loc_iat
    integer :: iat
    integer :: hat_on_proc, hat
    integer :: ngwf, ngwf_skip
    integer :: ngwf_fix_up, ngwf_fix_down
    integer :: hub_proj
    integer :: index, projector_index, ngwf_index

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_ngwf_on_grid'

    call timer_clock('hubbard_ngwf_on_grid',1)

    call data_set_to_zero(hub%consistency_projs)

    ngwf_count = 1
    ngwf_offset = 1
    hubbard_offset = 1

    do loc_iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       iat = par%first_atom_on_proc(pub_my_proc_id) + loc_iat - 1

       do hat_on_proc=1,par%num_hub_atoms_on_proc(pub_my_proc_id)

          hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
          if (iat .eq. par%distr_atom(hub%orig(hat)) ) then

             do hub_proj=1,hub_proj_basis%num_on_atom(iat)

                do ngwf=ngwf_basis%first_on_atom(iat),&
                     ngwf_basis%first_on_atom(iat)+ngwf_basis%num_on_atom(iat)-1

                   if (hub%localisedngwfs(hub_proj,hat) .eq. ngwf) then

                      do ngwf_fix_up=1,(ngwf - ngwf_basis%first_on_atom(iat))

                         ngwf_offset = ngwf_offset + &
                              cell%n_pts*ngwf_basis%spheres(&
                              &ngwf_count)%n_ppds_sphere

                         ngwf_count = ngwf_count + 1

                      end do

                      do index=1,cell%n_pts*&
                           &ngwf_basis%spheres(ngwf_count)%n_ppds_sphere

                         projector_index = index + hubbard_offset - 1
                         ngwf_index = index + ngwf_offset - 1

                         if (hub%consistency_projs%iscmplx) then
                            hub%consistency_projs%z(projector_index) = &
                                 hub%consistency_ngwfs%z(ngwf_index)
                         else
                            hub%consistency_projs%d(projector_index) = &
                                 hub%consistency_ngwfs%d(ngwf_index)
                         end if

                      end do

                      hubbard_offset = hubbard_offset + &
                           &cell%n_pts*&
                           &ngwf_basis%spheres(ngwf_count)%n_ppds_sphere

                      do ngwf_fix_down =1, &
                           (ngwf - ngwf_basis%first_on_atom(iat))

                         ngwf_count = ngwf_count - 1

                         ngwf_offset = ngwf_offset - &
                              cell%n_pts*ngwf_basis%spheres(&
                              &ngwf_count)%n_ppds_sphere

                      end do

                   end if

                end do

             end do

          end if

       end do

       do ngwf_skip=1,ngwf_basis%num_on_atom(iat)

          if (ngwf_count .lt. ngwf_basis%num) then
             ngwf_offset = ngwf_offset + &
                  cell%n_pts*ngwf_basis%spheres(ngwf_count)%n_ppds_sphere
          elseif (ngwf_count .gt. ngwf_basis%num) then
             call utils_abort('Something wrong in hubbard_ngwf_on_grid.')
          end if

          ngwf_count = ngwf_count + 1

       end do

    end do

    call timer_clock('hubbard_ngwf_on_grid',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_ngwf_on_grid'

  end subroutine hubbard_ngwf_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_projection_mtx(tensor,consistency, &
       hub_overlap, hub_overlap_t, &
       hub_o_matrix, hub_consist_matrix, hub_consist_tmatrix)

    !==========================================================!
    ! This matrix takes the  <NGWF|projector> matrix and the   !
    ! tensorial correction matrix if (tensor == .true.) and    !
    ! from then calculates the SPAM3 Hubbard projection.       !
    ! WARNING: in the tensor=true case, the sparsity pattern   !
    ! of hub_overlap_t is changed as a side-effect (at least   !
    ! in the extended metric approach)                         !
    !----------------------------------------------------------!
    ! Written by David O'Regan in April 2009.                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    ! Modified by David Turban in March 2014 to ensure         !
    ! consistency with CDFT-region duals (tensor correction 4) !
    !==========================================================!

    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_transpose, sparse_product
    use timer, only: timer_clock

    implicit none

    ! Arguments
    ! ddor: Apply tensorial correction for non-orthogonal Hubbard projectors
    logical, intent(in) :: tensor
    ! ddor: Apply tensorial correction for non-orthogonal Hubbard projectors
    !       for the sake of determining which NGWFs to use as projectors
    logical, intent(in) :: consistency
    type(SPAM3), intent(in) :: hub_overlap
    type(SPAM3), intent(inout) :: hub_overlap_t
    type(SPAM3), intent(in) :: hub_o_matrix
    type(SPAM3), optional, intent(in) :: hub_consist_matrix
    type(SPAM3), optional, intent(inout) :: hub_consist_tmatrix

    ! Local Variables
    type(SPAM3) :: temp_overlap_t

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_projection_mtx'

    call timer_clock('hubbard_projection_mtx',1)

    ! ddor: generate hub_consist_tmatrix instead of hub_overlap_t
    if (consistency) then

       ! ddor: Carry out tensorial correction for nonorthogonal projectors
       if (tensor) then
          ! ddor: Raise projector index in <Projector_NGWF> matrix
          !       Here the NGWF index is raised, not the projector, so the
          !       correction matrix appears on the right.
          call sparse_create(temp_overlap_t,hub_consist_tmatrix)
          call sparse_transpose(temp_overlap_t,hub_consist_matrix)
          ! dhpt: make sure that sparsity pattern of hub_consist_tmatrix is
          ! correct
          call sparse_destroy(hub_consist_tmatrix)
          call sparse_create(hub_consist_tmatrix,temp_overlap_t,hub_o_matrix)
          call sparse_product(hub_consist_tmatrix,temp_overlap_t,hub_o_matrix)
          call sparse_destroy(temp_overlap_t)
       else
          call sparse_transpose(hub_consist_tmatrix,hub_consist_matrix)
       end if

    else

       ! ddor: Carry out tensorial correction for nonorthogonal projectors
       if (tensor) then
          ! ddor: Raise projector index in <Projector_NGWF> matrix
          call sparse_create(temp_overlap_t,hub_overlap_t)
          call sparse_transpose(temp_overlap_t,hub_overlap)
          ! dhpt: make sure that sparsity pattern of hub_overlap_t is
          ! correct
          call sparse_destroy(hub_overlap_t)
          call sparse_create(hub_overlap_t,hub_o_matrix,temp_overlap_t)
          call sparse_product(hub_overlap_t,hub_o_matrix,temp_overlap_t)
          call sparse_destroy(temp_overlap_t)
       else
          call sparse_transpose(hub_overlap_t,hub_overlap)
       end if

    end if

    call timer_clock('hubbard_projection_mtx',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_projection_mtx'

  end subroutine hubbard_projection_mtx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_projector_character(hub,hub_proj_basis, &
       hub_consist_matrix)

    !==========================================================!
    ! This subroutine recomputes the hydrogenic character of   !
    ! NGWF DFT+U projectors using the updated                  !
    ! tensorial correction.                                    !
    !----------------------------------------------------------!
    ! Written by David O'Regan in October 2009.                !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    ! agrecocmplx: add dependence on COEF
    use datatypes, only: COEF
    use comms, only: pub_my_proc_id, comms_barrier, comms_reduce, pub_on_root
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_output_detail, pub_cdft_multi_proj, &
         pub_hub_tensor_corr
    use sparse, only: SPAM3, sparse_get_element, sparse_product
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(SPAM3), intent(in) :: hub_consist_matrix

    ! Local Variables
    integer :: hat_on_proc, hat, theatom
    integer :: projs
    integer :: ierr
    integer :: hub_proj_row_we_are_on, hub_proj_col_we_are_on
    integer :: hub_proj_count
    integer :: hub_hydro_we_are_on
    integer :: hub_tensorial_index_row, hub_tensorial_index_col
    integer :: hub_hydro_count, hub_tensorial_count
    !real(kind=DP) :: overlap_element
    !real(kind=DP) :: transpose_overlap_element, tensorial_element
    ! agrecocmplx: change these to COEF types
    type(COEF) :: overlap_element
    type(COEF) :: transpose_overlap_element, tensorial_element
    real(kind=DP), allocatable :: ngwfsquareovlp(:,:,:)
    ! agrecocmplx
    logical :: iscmplx

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_projector_character'

    call timer_clock('hubbard_projector_character',1)

    ! agrecocmplx
    iscmplx = hub_consist_matrix%iscmplx

    overlap_element%iscmplx = iscmplx
    transpose_overlap_element%iscmplx = iscmplx
    tensorial_element%iscmplx = iscmplx

    allocate(ngwfsquareovlp(hub_proj_basis%max_on_atom, &
         hub_proj_basis%max_on_atom,par%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_ngwf_select','ngwfsquareovlp',ierr)
    ngwfsquareovlp = 0.0_DP

    ! gibo: to avoid summing up contribution from different HUBBARDSCF cycles
    ! [otherwise output printing problems in internal_ngwf_character_info2]
    hub%current_projection(:) = 0.0_DP

    ! ddor: Loop over Hubbard atoms on my proc
    do hat_on_proc = 1, par%num_hub_atoms_on_proc(pub_my_proc_id)

       hat = par%hub_atoms_on_proc(hat_on_proc,pub_my_proc_id)
       theatom = par%distr_atom(hub%orig(hat)) !Atom number of Hubbard atom
       ! ddor: Number of Hubbard projectors on this Hubbard atom
       projs = hub_proj_basis%num_on_atom(theatom)

       ! ddor: loop over the row projectors on this Hubbard atom
       do hub_proj_count = 1, projs

          ! ddor: The rows of <NGWF_|phi> and columns of <NGWF^|NGWF^>
          !       are indexed differently
          hub_proj_row_we_are_on = &
               hub%localisedngwfs(hub_proj_count,hat)
          hub_proj_col_we_are_on = hub_proj_count + &
               ( hub_proj_basis%first_on_atom(theatom) - 1 )

          ! ddor: loop over the hydrogenic projectors on this Hubbard atom
          do hub_hydro_count = 1, projs

             hub_hydro_we_are_on = hub_hydro_count + &
                  hub_proj_basis%first_on_atom(theatom) - 1

             ! ddor: Add up the square of the overlap between
             !       the current NGWF and
             !       the contravariant projectors.
             !       As well as correcting for nonorthogonality,
             !       this fixes the normalisation of the NGWF projectors.
             ! agrecocmplx
             if (iscmplx) then
                call sparse_get_element(overlap_element%z,hub_consist_matrix,&
                     hub_proj_row_we_are_on, hub_hydro_we_are_on)
             else
                call sparse_get_element(overlap_element%d,hub_consist_matrix,&
                     hub_proj_row_we_are_on, hub_hydro_we_are_on)
             end if

             ! ddor: Loop over rows of the contravariant matrix on the site
             do hub_tensorial_count = 1, projs

                hub_tensorial_index_col = &
                     hub%localisedngwfs(hub_tensorial_count,hat)

                hub_tensorial_index_row = &
                     hub_tensorial_count + &
                     hub_proj_basis%first_on_atom(theatom) - 1

                ! ddor: Hydrogenic rows, NGWF projector columns,
                !       but use the transpose
                ! agrecocmplx
                if (iscmplx) then
                   call sparse_get_element(transpose_overlap_element%z,&
                        hub_consist_matrix,&
                        hub_tensorial_index_col, hub_hydro_we_are_on)

                   call sparse_get_element(tensorial_element%z,&
                        hub%o_matrix,&
                        hub_tensorial_index_row, hub_proj_col_we_are_on)
                else
                   call sparse_get_element(transpose_overlap_element%d,&
                        hub_consist_matrix,&
                        hub_tensorial_index_col, hub_hydro_we_are_on)

                   call sparse_get_element(tensorial_element%d,&
                        hub%o_matrix,&
                        hub_tensorial_index_row, hub_proj_col_we_are_on)
                end if

                ! ddor: The elements (alpha, i)
                !                 { <NGWF_alpha|\phi_i>    }
                !       sum_beta  { <\phi_i|NGWF_beta>     }
                !                 { <NGWF^beta|NGWF^alpha> }
                ! agrecocmplx
                if (iscmplx) then
                   ngwfsquareovlp(hub_hydro_count,&
                        hub_proj_count,hat) = &
                        ngwfsquareovlp(hub_hydro_count,&
                        hub_proj_count,hat) + &
                        real(overlap_element%z * transpose_overlap_element%z * &
                        tensorial_element%z, kind=DP)
                else
                   ngwfsquareovlp(hub_hydro_count,&
                        hub_proj_count,hat) = &
                        ngwfsquareovlp(hub_hydro_count,&
                        hub_proj_count,hat) + &
                        overlap_element%d * transpose_overlap_element%d * &
                        tensorial_element%d
                end if

             end do

          end do

       end do

       hub%current_projection(hat) = &
            SUM(ngwfsquareovlp(1:projs,1:projs,hat))

    end do ! ddor: End loop over Hubbard atoms on my proc

    ! ddor: Reduce, Allow other procs to catch up
    call comms_reduce('SUM',ngwfsquareovlp)
    call comms_reduce('SUM',hub%current_projection)

    call comms_barrier

    ! ddor: Write out the | < \phi | \psi > |**2 matrix
    if (pub_on_root .and. pub_output_detail >= VERBOSE) &
         &call internal_ngwf_character_info2

    deallocate(ngwfsquareovlp,stat=ierr)
    call utils_dealloc_check('hubbard_ngwf_select','ngwfsquareovlp',ierr)

    call timer_clock('hubbard_projector_character',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_projector_character'

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_ngwf_character_info2

      implicit none

      ! Local Variables
      integer :: project, jj
      integer :: kk, curr_proj, species
      character(10) :: fmt, tmp

      do hat=1,par%nat_hub

         ! ddor: Number of Hubbard projectors on this Hubbard atom
         theatom = par%distr_atom(hub%orig(hat))
         project = hub_proj_basis%num_on_atom(theatom)

         !gibo: added for cdft_multi_proj
         if (pub_cdft_multi_proj) species = hub%species_number(hat)

         write(stdout,'(/a)') repeat('#',80)

         ! ddor: Write a warning to ignore this in certain circumstances
         if (pub_hub_tensor_corr == 4 .or. pub_hub_tensor_corr == 5) &
              write(stdout,'(a)') 'The following character matrices&
                   & are unsuitable for hubbard_tensor_corr : 4 or 5.'

         ! ddor: WRITE OUT INFORMATION ON | < \phi | \psi > |**2 MATRIX
         write(stdout,'(a)') 'Using the subspace &
              &contravariant metric, &
              &the NGWF hydrogenic character'
         write(stdout,'(a,i6,a)') '| < \phi | \psi > |**2 &
              &matrix of Hubbard site ',hat,' is '

         ! ndmh: write header line
         if (.not.pub_cdft_multi_proj) then
            write(stdout,'(a)',advance='no') '  m_l ='
            do jj=-(project-1)/2,(project-1)/2
               write(stdout,'(i5,7x)',advance='no') jj
            end do
            write(stdout,'(a)') '  SUM'

         elseif (pub_cdft_multi_proj) then
            write(stdout,'(a)',advance='no') 'm_l(l)='

            do kk = 1,  h_species(species)%num_ang_mom_shells
               curr_proj =  2*h_species(species)%orb_ang_mom(kk) + 1
               do jj=-(curr_proj-1)/2,(curr_proj-1)/2
                  write(stdout,'(i5,6x)',advance='no') jj
               end do
            enddo

            write(stdout,'(a)') '  SUM'

         endif

         ! ndmh: write format string
         write(tmp,'(i6)') project+1
         write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f12.8)'

         ! ndmh: write matrix
         do jj=1,project
            write(stdout,fmt) ngwfsquareovlp(1:project,jj,hat), &
                 SUM(ngwfsquareovlp(1:project,jj,hat))
         end do

         write(stdout,'(a,f12.8)') 'Total hydrogenic projection of &
              &NGWFs on this correlated site', hub%current_projection(hat)

         write(stdout,'(a/)') repeat('#',80)

      end do

    end subroutine internal_ngwf_character_info2

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine hubbard_projector_character


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_projector_update( &
       ngwfs_on_grid, ngwf_basis, nl_projectors, projector_basis, &
       hub_proj_basis, hub, inv_overlap, hub_overlap, hub_overlap_t, &
       sp_overlap, hub_proj_paw_overlap, hub_ngwf_paw_overlap, mdl, kpt)

    !==================================================================!
    ! This subroutine generates (on the first HUBBARDSCF iteration) or !
    ! updates (on subsequent HUBBARDSCF iterations) the consistency    !
    ! NGWFs, chooses which of these should become the hubbard          !
    ! projectors, and generates the required overlap matrices          !
    !------------------------------------------------------------------!
    ! Written by David O'Regan in July 2009.                           !
    ! Modified by Nicholas Hine in November 2011 to use the            !
    ! HUBBARD_MODEL type                                               !
    ! Modified by D. D. O'Regan for on-the-fly HUBBARDSCF in July 2013.!
    ! Modified by Andrea Greco on 21/06/2015 to use complex NGWFs.     !
    !==================================================================!

    use augmentation, only: augmentation_overlap
    ! agrecocmplx
    use datatypes, only: FUNCTIONS, data_functions_mix, data_functions_copy
    use function_basis, only: FUNC_BASIS
    use function_ops, only: function_ops_brappd_ketppd
    use geometry, only: POINT ! agrecokpt
    use ion, only: element
    use model_type, only: MODEL
    use ngwfs, only: ngwfs_initialise_ngwfs
    use projectors, only: PROJECTOR_SET, projectors_func_ovlp_box
    use restart, only : restart_ngwfs_tightbox_output, &
         restart_ngwfs_tightbox_input
    use rundat, only: pub_hub_proj_mixing, &
         pub_write_tightbox_ngwfs, pub_maxit_ngwf_cg, &
         pub_hubbard_atomsolve, pub_hub_on_the_fly, pub_aug, &
         pub_task, pub_hub_proj_read_only, pub_kpoint_method
    use sparse, only: SPAM3, sparse_scale, sparse_axpy,&
         sparse_copy, sparse_create, sparse_destroy
    use timer, only: timer_clock
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: ngwf_basis
    type(FUNC_BASIS), intent(in)    :: projector_basis
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis
    !real(kind=DP), intent(in)       :: ngwfs_on_grid(&
    !     &ngwf_basis%size_on_grid)
    ! agrecocmplx: use new FUNCTIONS type
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    type(PROJECTOR_SET), intent(inout) :: nl_projectors
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(SPAM3), intent(in)         :: inv_overlap
    type(SPAM3), intent(inout)      :: hub_overlap
    type(SPAM3), intent(inout)      :: hub_overlap_t
    ! ddor: Needed for PAW+U
    type(SPAM3), intent(in)         :: sp_overlap
    type(SPAM3), intent(inout)      :: hub_proj_paw_overlap
    type(SPAM3), intent(inout)      :: hub_ngwf_paw_overlap
    type(MODEL), optional, intent(inout) :: mdl
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    ! ddor: Contravariant <projector|projector> matrix for tensorial correction
    !       In this case it is for all the NGWFs on a Hubbard atom, so that
    !       the most localised ones can be chosen in tensorially consistent way.
    type(SPAM3) :: hub_consist_o_matrix
    ! ddor: <NGWF|projector> Matrix used to determine which NGWFs
    !        to use as projectors
    type(SPAM3) :: hub_consist_matrix
    ! ddor: Corrected <projector|NGWF> Matrix used to determine
    !       which NGWFs to use as projectors
    type(SPAM3) :: hub_consist_tmatrix
    ! agrecokpt
    type(POINT) :: loc_kpt


    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_projector_update'

    call timer_clock('hubbard_projector_update',1)

    ! agrecokpt
    if (present(kpt)) then
       loc_kpt = kpt
    else
       loc_kpt%x = 0.0_DP
       loc_kpt%y = 0.0_DP
       loc_kpt%z = 0.0_DP
    end if

    ! agrecokpt: currently only KP method is being implemented
    call utils_assert(pub_kpoint_method == 'KP', &
         'Subroutine hubbard_projector_update currently supports&
         & only KP method for BZ sampling')

    ! agrecokpt: currently only ensuring compatibility with complex
    ! quantities, plus including KP dependence in projectors when
    ! using proj_overlap_box routine; need to check KP dependence in
    ! Hubbard formalism
    if (loc_kpt%x/=0.0_DP.or.loc_kpt%y/=0.0_DP.or.loc_kpt%z/=0.0_DP) then
       if (pub_debug_on_root) write(stdout,'(a)') &
           'WARNING: update of Hubbard projectors in &
           &hubbard_projector_update not yet checked for &
           &non-Gamma k-point'
    end if

    ! ddor: On first projector optimisation iteration do a plain
    !       ONETEP calculation, unless we have projectors stored in
    !       tightbox representation, in which case we go straight to
    !       DFT+U with NGWF projectors.
    if (( hub%consistency_iteration > 1) .or. &
         & (pub_hub_proj_mixing .lt. 0.0_DP) .or. &
         & pub_hubbard_atomsolve .or. pub_hub_on_the_fly ) then

       ! ddor: If necessary, mix in some of previous projectors
       !       to aid convergence
       if ((pub_hubbard_atomsolve) .and. (hub%consistency_iteration==1) .and. &
            & (.not. pub_hub_on_the_fly)) then

          ! ndmh: initialise to get original PAO projectors back
          ngwf_basis%name = 'hub_ngwfs'
          ! rc2013: only valid for 1 region
          call ngwfs_initialise_ngwfs(hub%consistency_ngwfs,ngwf_basis,mdl,1)
          ngwf_basis%name = 'ngwfs'

       else if ( ( pub_hub_proj_mixing .eq. 0.0_DP )  &
            & .or. ( hub%consistency_iteration == 2 ) &
            & .or. ((pub_hub_proj_mixing .eq. -1.0_DP) .and. &
            & ( hub%consistency_iteration > 1 )) &
            & .or. pub_hub_on_the_fly ) then

          ! ddor: Nothing to mix
          ! ndmh_pointerfun
          call data_functions_copy(hub%consistency_ngwfs,ngwfs_on_grid)
          !hub%consistency_ngwfs = ngwfs_on_grid

          ! ddor: Restore normal value of pub_maxit_ngwf_cg
          if ((pub_hub_proj_mixing .lt. 0.0_DP).and.(.not. pub_hub_proj_read_only)) then
             pub_maxit_ngwf_cg = hub%store_maxit_ngwf_cg
          end if

       else if (( pub_hub_proj_mixing .lt. 0.0_DP ) .and. &
            ( hub%consistency_iteration .eq. 1)) then

          ! ddor: If mixing is negative, take this to
          !       mean read Hubbard NGWF projectors from
          !       file in tightbox representation
          call restart_ngwfs_tightbox_input(&
               hub%consistency_ngwfs,ngwf_basis,mdl%cell,mdl%fftbox, &
               mdl%elements,'tightbox_hub_projs',mdl%regions(1))

          ! ddor: Do no NGWF iterations on first HUBBARDSCF iteration
          !       when restarting.
          hub%store_maxit_ngwf_cg = pub_maxit_ngwf_cg
          !gibo: allow for cDFT-optimisations with cDFT-projectors from file
          !pub_maxit_ngwf_cg = 0
          if ((pub_task.NE.'SINGLEPOINT')          .AND.  &
               (pub_task.NE.'GEOMETRYOPTIMIZATION') .AND.  &
               (pub_task.NE.'MOLECULARDYNAMICS')    .AND.  &
               (pub_task.NE.'PHONON')               .AND.  &
               (pub_task.NE.'COND')                 .AND.  &
               (pub_task.NE.'TRANSITIONSTATESEARCH') .AND. &
               (.not. pub_hub_proj_read_only)) pub_maxit_ngwf_cg = 0

       else

          ! ddor: Mix the projectors for iterations > 2
          call data_functions_mix(hub%consistency_ngwfs, ngwfs_on_grid, &
               ABS(pub_hub_proj_mixing))

       end if

       ! ddor: Overlap of Hubbard NGWF projectors
       call function_ops_brappd_ketppd(hub%consistency_overlap, & !output
            hub%consistency_ngwfs, ngwf_basis, &               !input
            hub%consistency_ngwfs, ngwf_basis, mdl%cell)

       ! ddor: Augment the <NGWF|NGWF> matrix for PAW+U or USP calculations
       if (pub_aug) then

          ! ddor: Calculate the overlap matrix of the Hubbard
          ! ddor: and PAW projectors
          ! agrecokpt: add kpoint dependence
          ! (NEED TO CHECK THIS!!!!)
          call projectors_func_ovlp_box(hub_ngwf_paw_overlap, &
               hub%consistency_ngwfs,ngwf_basis,projector_basis, &
               nl_projectors,mdl%fftbox,mdl%cell)!,kshift=loc_kpt)

          ! ndmh: Calculate the augmentation of the overlap matrix due to the
          ! ndmh: augmentation region part of the overlap operator
          call augmentation_overlap(hub%consistency_overlap,mdl%pseudo_sp, &
               mdl%paw_sp,hub_ngwf_paw_overlap)

       end if

       call sparse_create(hub_consist_matrix,hub_overlap)
       call sparse_create(hub_consist_tmatrix,hub_overlap_t)

       ! ddor: Renew hubbard_consistency_matrix
       !       (the overlap of NGWFs and atomic orbitals)
       !       needed to select which NGWFs have most localised character
       !       on this iteration.
       ! agrecokpt: add kpoint dependence
       ! (NEED TO CHECK THIS!!!!)
       call projectors_func_ovlp_box(hub_consist_matrix, &
            hub%consistency_ngwfs,ngwf_basis,hub_proj_basis,hub%projectors, &
            mdl%fftbox,mdl%cell)!,kshift=loc_kpt

       hub_consist_o_matrix%structure = 'D'
       ! agrecocmplx: need to be complex if using complex NGWFs
       call sparse_create(hub_consist_o_matrix, iscmplx=hub_consist_matrix%iscmplx)

       ! ddor: Find the tensorial correction matrix for all candidate
       !       NGWF projectors on Hubbard atoms
       call hubbard_tensorial_correction(hub,ngwf_basis,hub_proj_basis,.true.,&
            inv_overlap, hub_consist_o_matrix)

       ! ddor: With hub_consist_matrix and the tensorial correction matrix,
       !       calculate the Hubbard projection operator
       !       for all candidate projectors
       ! ndmh: this call does not use hub_overlap or hub_overlap_t
       call hubbard_projection_mtx(.true.,.true.,hub_overlap,hub_overlap_t, &
            hub_consist_o_matrix, hub_consist_matrix, hub_consist_tmatrix)

       call sparse_destroy(hub_consist_o_matrix)

       ! ddor: Select the most localised NGWFs to use as Hubbard
       !       projectors on each correlated site.
       !       Using the hub_overlap from previous Hubbard iteration
       !       (e.g. hubbard_consistency_overlap) calculate
       !       the new tensorial correction matrix, the contravariant
       !       covariant metric on each correlated site.
       call hubbard_ngwf_select(hub, ngwf_basis, hub_proj_basis,  &
            hub_consist_matrix, hub_consist_tmatrix)

       ! ddor: Move the Hubbard NGWF projectors on the grid
       !       to an array containing the minimal set.
       call hubbard_ngwf_on_grid(hub,ngwf_basis, hub_proj_basis,mdl%cell)

       ! ddor: Compute new NGWF-projector overlap matrix
       call function_ops_brappd_ketppd(hub_overlap, &              !input-output
            ngwfs_on_grid, ngwf_basis, hub%consistency_projs, & !input
            hub_proj_basis, mdl%cell)                           !input

       ! ddor: Augment the <NGWF|proj> matrix for PAW+U or USP calculations
       if (pub_aug) then

          ! ddor: Calculate the overlap matrix of the Hubbard
          ! ddor: and PAW projectors
          ! agrecokpt: not sure if kpoint dependence needs to be added here
          ! (CHECK THIS!!!)
          call projectors_func_ovlp_box(hub_proj_paw_overlap, &
               hub%consistency_projs,hub_proj_basis,&
               projector_basis,nl_projectors,mdl%fftbox,mdl%cell)!,kshift=loc_kpt

          ! ndmh: Calculate the augmentation of the overlap matrix due to the
          ! ndmh: augmentation region part of the overlap operator
          call augmentation_overlap(hub_overlap,mdl%pseudo_sp,mdl%paw_sp, &
               sp_overlap,hub_proj_paw_overlap)

       end if

       ! ddor: write new Hubbard NGWF projectors to file in tightbox rep.
       !       write out all NGWFs, so we won't need projector tightboxes,
       !       or write out hub%consistency_overlap.
       ! ebl:  do not write to file in case of pub_hub_proj_read_only (for DMFT)
       if ( (.not. pub_hub_proj_read_only) .and. (hub%consistency_iteration &
            .gt. 0) .and. pub_write_tightbox_ngwfs) then
          call restart_ngwfs_tightbox_output(&
               hub%consistency_ngwfs,ngwf_basis,mdl%cell,mdl%fftbox, &
               mdl%elements, 'tightbox_hub_projs',mdl%regions(1))
       end if

       ! ddor: Find the tensorial correction matrix for
       !       the chosen NGWF Hubbard projectors
       ! gom: When combining all projectors into one site, and in case where
       !      the number of projectors equals the number of NGWFs, the subspace
       !      metric O is just the overlap matrix S of the current NGWFs

       call hubbard_tensorial_correction(hub,ngwf_basis,hub_proj_basis,&
            .false.,inv_overlap,hub%o_matrix)

       ! ddor: Recalculate the hydrogenic character of the NGWF projectors
       !       using the appropriate tensorial correction
       call hubbard_projector_character(hub,hub_proj_basis, &
            hub_consist_matrix)

       call sparse_destroy(hub_consist_tmatrix)
       call sparse_destroy(hub_consist_matrix)

       ! ddor: With hub_overlap and the tensorial correction matrix,
       !       calculate the Hubbard projection operator for each site
       call hubbard_projection_mtx(.true.,.false.,hub_overlap,hub_overlap_t, &
            hub%o_matrix)

    end if

    call timer_clock('hubbard_projector_update',2)

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_projector_update'

  end subroutine hubbard_projector_update


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function hubbard_test_convergence(hub,total_energy)

    !==========================================================!
    ! When carrying out a DFT+U calculation with               !
    ! self-consistent NGWF projectors, this function tests     !
    ! the convergence by seeing if the energy difference       !
    ! between ground-state energies with is less than          !
    ! hub_energy_tol over hub_conv_win projector optimisation  !
    ! steps.                                                   !
    !                                                          !
    ! The self-consistent scheme also halts if the sum of the  !
    ! projections onto hydrogenic orbitals of the NGWF Hubbard !
    ! projectors remains constant over hub_conv_win projector  !
    ! optimisation steps.                                      !
    !----------------------------------------------------------!
    ! Written by David O'Regan in April 2009.                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_on_root
    use parallel_strategy, only: par=>pub_par
    use rundat, only: pub_hub_energy_tol,pub_hub_conv_win

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(in) :: total_energy

    ! Local Variables
    integer ::  hub_energies_flag,  hub_projections_flag
    integer ::  hub_iterations_counter

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_test_convergence'

    ! ddor: Zero if the energy convergence criterion has been met
    hub_energies_flag = 0
    ! ddor: Zero if the hydrogenic projection convergence criterion has been met
    hub_projections_flag = 0

    ! ddor: Check the differences of historical energies and projections
    do hub_iterations_counter = 1, (pub_hub_conv_win - 1)

       ! ddor: Energies
       if (ABS( hub%consistency_energies(hub_iterations_counter) - &
            hub%consistency_energies(hub_iterations_counter + 1 ) )&
            & .gt. pub_hub_energy_tol ) then
          ! ddor: The difference of energies exceeds the threshold
          hub_energies_flag = 1
       end if
       ! ddor: Move energies down in the array
       hub%consistency_energies(hub_iterations_counter) = &
            hub%consistency_energies(hub_iterations_counter + 1 )

       ! ddor: Projections
       if (ABS( hub%consistency_projections(hub_iterations_counter) - &
            hub%consistency_projections(hub_iterations_counter + 1 ) )&
            & .gt. 1.0e-8_DP ) then
          ! ddor: The difference of projections exceeds the threshold
          hub_projections_flag = 1
       end if
       ! ddor: Move projections down in the array
       hub%consistency_projections(hub_iterations_counter) = &
            hub%consistency_projections(hub_iterations_counter + 1 )

    end do

    ! ddor: fill the top of the array with the current energy
    hub%consistency_energies(pub_hub_conv_win) = total_energy

    ! ddor: fill the top of the array with the current
    !       sum of hydrogenic projections of Hubbard NGWFs of all sites
    if (hub%consistency_iteration .gt. 1) then
       hub%consistency_projections(pub_hub_conv_win) = &
            SUM(hub%current_projection(1:par%nat_hub))
    end if

    ! ddor: Check the difference between current and previous energy
    if (ABS( hub%consistency_energies(pub_hub_conv_win - 1 ) - &
         hub%consistency_energies(pub_hub_conv_win) ) .gt. &
         &pub_hub_energy_tol ) then
       hub_energies_flag = 1
    end if

    ! ddor: Check the difference between current and previous projection
    if (ABS( hub%consistency_projections(pub_hub_conv_win - 1 ) - &
         hub%consistency_projections(pub_hub_conv_win ) ) .gt. &
         & 1.0e-8_DP ) then
       hub_projections_flag = 1
    end if

    ! ddor: If we still think that the energy has converged, and we've
    !       done sufficiently many iterations, then pass the convergence test
    if ( ( hub_projections_flag .eq. 0) &
         & .and. (hub%consistency_iteration .gt. pub_hub_conv_win+1)) then
       hubbard_test_convergence = .true.
       ! ddor: Write out a congratulatory message
       if (pub_on_root) then
          write(stdout,'(/a)') repeat('#',80)
          write(stdout,'(a,i6)')'DFT+U projectors optimised on iteration ', &
               hub%consistency_iteration
          if (hub_energies_flag .eq. 0) then
             write(stdout,'(a,f12.8,a,i6,a)') 'Energy tolerance of ', &
                  pub_hub_energy_tol,' maintained over ',pub_hub_conv_win, &
                  ' iterations'
          end if
          if (hub_projections_flag .eq. 0) then
             write(stdout,'(a,f12.8,a,i6,a)') 'Projector character &
                  &tolerance of ', 1.0e-8_DP,' maintained over ',&
                  pub_hub_conv_win, ' iterations'
          end if
          write(stdout,'(a/)') repeat('#',80)
       end if
    else
       ! ddor: Otherwise, keep going.
       hubbard_test_convergence = .false.
    end if

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_test_convergence'

  end function hubbard_test_convergence


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_species_proj(hub,elements,fftbox,is_cmplx)

    !=================================================================!
    ! This subroutine allocates all the memory for                    !
    ! h_proj data on hydrogenic projectors used by DFT+U and          !
    ! calculates this data for each Hubbard species.                  !
    ! This subroutine allocates and initialises hub_fftbox_proj_recip !
    ! for each h_species element. Each such hub_fftbox_proj_recip is  !
    ! a full non-local projector in reciprocal space in the FFTbox.   !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/1/2004.                  !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                       !
    ! Modified for DFT+U module by David O'Regan in July 2008.        !
    ! Modified by Nicholas Hine in November 2011 to use the           !
    ! HUBBARD_MODEL type                                              !
    !=================================================================!

    use comms, only: pub_on_root
    use constants, only: DP
    use fft_box, only: FFTBOX_INFO
    use geometry, only: POINT, OPERATOR(.dot.)
    use hubbard_init, only: h_species
    use ion, only: element
    use parallel_strategy, only: par=>pub_par
    use projectors, only: projectors_allocate_set
    use rundat, only: pub_output_detail, &
         pub_cdft_multi_proj !gibo: 28.12.12
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(FFTBOX_INFO), intent(in) :: fftbox
    ! agrecocmplx
    logical, optional, intent(in) :: is_cmplx

    ! Local Variables
    integer :: jj, channels  !gibo: for multi-proj cDFT
    integer :: n1, n2, n3, ld1, ld2
    integer :: ps
    integer :: proj_count, hat
    integer :: iat, orig_iat
    real(kind=DP) :: spacing,deltar ! Grid variables
    type(POINT) :: kpt_loc
    ! agrecocmplx
    logical :: loc_cmplx

    ! ndmh: removed hubbard_proj type, so these variables need to be local
    ! ddor: number of points in the real radial grid
    integer :: n_r_rad_pts
    ! ddor: number of points in the reciprocal radial grid
    integer :: n_g_rad_pts, max_n_g_rad_pts
    ! ddor: maximum radial reciprocal k-vector up to which
    !      the hydrogenic projectors are defined
    real(kind=DP) :: gmax
    ! ddor: hydrogenic projector radius in real space
    real(kind=DP) :: rmax

    if (pub_on_root  .and. pub_output_detail >= VERBOSE) &
         write(stdout,'(a)') &
         '... Hubbard projector initialisation'

    if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_species_proj'

    call timer_clock('hubbard_species_proj',1)

    ! agrecocmplx
    if (present(is_cmplx)) then
       loc_cmplx = is_cmplx
    else
       loc_cmplx = .false.
    end if

    ! Local copies of FFT box dimensions
    n1 = fftbox%total_pt1
    n2 = fftbox%total_pt2
    n3 = fftbox%total_pt3
    ld1 = fftbox%total_ld1
    ld2 = fftbox%total_ld2

    kpt_loc%x = 0.0_DP ; kpt_loc%y = 0.0_DP ; kpt_loc%z = 0.0_DP

    ! ndmh: find max value of recip grid points
    max_n_g_rad_pts = -1
    spacing = min(fftbox%d1 ,fftbox%d2, fftbox%d3)
    do ps=1,par%num_hub_species
       do hat=1,par%nat_hub
          orig_iat = hub%orig(hat)
          if (hub%species_number(hat) == ps) exit
       end do
       rmax = elements(orig_iat)%radius
       n_r_rad_pts = int(100.0_DP*rmax/spacing)
       deltar = rmax / real(n_r_rad_pts - 1,kind=DP)
       n_g_rad_pts = nint(rmax/deltar) + 1
       max_n_g_rad_pts = max(max_n_g_rad_pts,n_g_rad_pts)
    end do

    ! ndmh: count unique projectors
    hub%projectors%n_proj_species = par%num_hub_species
    proj_count = 0
    do ps=1,par%num_hub_species
       !gibo: projectors of different ang.mom for multi-proj cDFT
       if (pub_cdft_multi_proj) then
          proj_count = proj_count + h_species(ps)%cdft_num_proj
       else
          !gibo: standard DFT+U, single ang.mom channel cDFT counting
          proj_count = proj_count + h_species(ps)%hub_ang_mom*2 + 1
       endif
    end do

    ! ndmh: allocate hub%projectors type
    !call projectors_allocate_set(hub%projectors,1,max_n_g_rad_pts)
    !gibo: accounting for projectors of multiple ang.mom on cDFT-sites
    if (pub_cdft_multi_proj) then
       ! set 'channels' to the the maximum total number of cDFT-proj=NGWFs
       channels = maxval(h_species(:)%cdft_num_proj)
    else
       channels = 1
    endif
    ! agrecocmplx
    call projectors_allocate_set(hub%projectors,channels,max_n_g_rad_pts, &
         par, is_cmplx=loc_cmplx)

    hub%projectors%normalise = .true.

    !gibo: just to be sure not to mess with radial part of cDFT-projectors
    if (pub_cdft_multi_proj) hub%projectors%rad_proj_recip(:,:,:) = 0.0_DP

    ! ndmh: set up entries in hub%projectors type
    proj_count = 1
    do ps=1,par%num_hub_species
       !gibo: allow for multi-proj cDFT
       if (pub_cdft_multi_proj) then
          hub%projectors%species_num_proj(ps) = &
               h_species(ps)%cdft_num_proj
       else
          hub%projectors%species_num_proj(ps) = h_species(ps)%hub_ang_mom*2 + 1
       endif

       ! ndmh: find Gmax for this species
       do hat=1,par%nat_hub
          orig_iat = hub%orig(hat)
          if (hub%species_number(hat) == ps) exit
       end do
       rmax = elements(orig_iat)%radius
       n_r_rad_pts = int(100.0_DP*rmax/spacing)
       deltar = rmax / real(n_r_rad_pts - 1,kind=DP)
       n_g_rad_pts = nint(rmax/deltar) + 1

       ! ndmh: calculate projector in reciprocal space for this species
       !gibo: modified for multi-proj cDFT
       if (pub_cdft_multi_proj) then
          call hubbard_internal_radial_multi_proj(h_species(ps))

       else
          call hubbard_internal_radial_proj(h_species(ps))

       endif

       hub%projectors%n_rad_pts(ps) = n_g_rad_pts
       hub%projectors%gmax(ps) = gmax

       !-----
       !gibo: modified for multi-proj cDFT
       !hub%projectors%num_shells(ps) = 1
       !hub%projectors%ang_mom(1,ps) = h_species(ps)%hub_ang_mom
       if (pub_cdft_multi_proj) then
          hub%projectors%num_shells(ps) =  h_species(ps)%num_ang_mom_shells
          do jj = 1, hub%projectors%num_shells(ps)
             hub%projectors%ang_mom(jj,ps) = h_species(ps)%orb_ang_mom(jj)
          enddo

       else
          hub%projectors%num_shells(ps) = 1
          hub%projectors%ang_mom(1,ps) = h_species(ps)%hub_ang_mom

       endif

       !       !hub%projectors%num_shells(ps) = 1
       !       !hub%projectors%ang_mom(1,ps) = h_species(ps)%hub_ang_mom
       !       !gibo: number of different angular-momentum shell
       !       hub%projectors%num_shells(ps) =  h_species(ps)%hub_ang_mom + 1
       !       do jj = 1, hub%projectors%num_shells(ps)
       !        hub%projectors%ang_mom(jj,ps) = h_species(ps)%hub_ang_mom - (jj-1)
       !       enddo
       !-----

       hub%projectors%species_first_proj(ps) = proj_count
       proj_count = proj_count + hub%projectors%species_num_proj(ps)
    end do

    ! ndmh: copy species number, atom centre and NGWF radius for Hubbard atoms
    hub%projectors%proj_species(:) = -1
    hub%projectors%proj_centre(:)%X = 0.0_DP
    hub%projectors%proj_centre(:)%Y = 0.0_DP
    hub%projectors%proj_centre(:)%Z = 0.0_DP
    hub%projectors%proj_max_radius(:) = -999.999_DP
    do hat=1,par%nat_hub
       orig_iat = hub%orig(hat)
       iat = par%distr_atom(orig_iat)
       hub%projectors%proj_species(iat) = hub%species_number(hat)
       hub%projectors%proj_max_radius(iat) = elements(orig_iat)%radius
       hub%projectors%proj_centre(iat) = elements(orig_iat)%centre
    end do

    call timer_clock('hubbard_species_proj',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving hubbard_species_proj'

  contains


    subroutine hubbard_internal_radial_proj(hd_species)

      !==================================================================!
      ! This module builds the radial part of hydrogenic wavefunctions   !
      ! in reciprocal space.                                             !
      !------------------------------------------------------------------!
      ! Written by David D. O'Regan in April 2009.                       !
      !==================================================================!

      use comms, only: pub_on_root
      use constants, only: PI
      use hubbard_init, only: HUBBARD_SPECIES
      use utils, only: utils_abort

      implicit none

      ! Arguments
      type(HUBBARD_SPECIES), intent(in) :: hd_species

      ! Local variables
      integer :: g_rad_pt                       ! Grid point counters
      real(kind=DP) :: gnorm                    ! Norm in recip space
      real(kind=DP) :: delg                     ! Grid variables
      real(kind=DP) :: gradius, denom, gath
      real(kind=DP) :: proj_max

      ! ddor: Set parameters for radial grid
      delg = PI / rmax
      proj_max = 0.0_DP
      gnorm = 0.0_DP
      gmax = 0.0_DP

      gath = ( hd_species%hub_charge )**( 0.5_DP*( &
           &real(2*hd_species%hub_ang_mom + 5,kind=DP) ) )

      ! ddor: Construct radial component in real space
      do g_rad_pt=1,n_g_rad_pts

         gradius =  real(g_rad_pt-1,DP) * delg

         denom = ( (hd_species%hub_charge)**2.0_DP + &
              (gradius * real(hd_species%hub_ang_mom + 1))**2.0_DP )**&
              (-1.0_DP*real(hd_species%hub_ang_mom) - 2.0_DP)

         select case (hd_species%hub_ang_mom)

         case(0) ! 1s

            hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
                 SQRT(32.0_DP / PI) * gath * denom

         case(1) ! 2p

            hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
                 (128.0_DP / SQRT(3.0_DP*PI) ) * gradius * gath * denom

         case(2) ! 3d

            hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
                 1728.0_DP * SQRT( 3.0_DP / (5.0_DP*PI) ) * &
                 &(gradius**2.0_DP) * gath * denom

         case(3) ! 4f

            hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
                 131072.0_DP * SQRT( 2.0_DP / (35.0_DP*PI) ) * &
                 &(gradius**3.0_DP) * gath * denom

         case default
            call utils_abort('Radial or angular quantum number too high in &
                 &hubbard_internal_get_radial_proj. Only up to r=4,l=3 &
                 &atomic radial functions are currently supported.')
         end select

         gnorm = gnorm + delg*(gradius * &
              hub%projectors%rad_proj_recip(g_rad_pt,1,ps))**2.0_DP

         ! ddor: Set core radius if magnitude of
         !       projector is less then 0.01% of
         !       maximum value and cumulative norm is at least 99.99%.
         proj_max = max(proj_max,abs(hub%projectors%rad_proj_recip(g_rad_pt,1,ps)))
         if (g_rad_pt > 1) then
            if ( (abs(hub%projectors%rad_proj_recip(g_rad_pt,1,ps)/proj_max) < &
                 &1.0e-4_DP) .and. (gnorm > 0.9999_DP) )then
               gmax = gradius ! The maximum necessary G
               n_g_rad_pts = g_rad_pt
               exit
            end if
         end if

      end do

      ! ddor: Sanity check on Fourier-Bessel transform
      if  ( (abs(hub%projectors%rad_proj_recip(n_g_rad_pts,1,ps)/proj_max) > &
           &1.0e-4_DP) .or. (gnorm < 0.9999_DP) ) then
         gmax = real(n_g_rad_pts-1,DP) * delg
         if (pub_on_root) write(stdout,*) 'Fourier-Bessel normalisation',&
              ' not reached in', &
              ' hubbard_internal_get_radial_proj. Normalisation',&
              ' attained was.',gnorm*100,'%',' gmax set to',gmax
      end if

    end subroutine hubbard_internal_radial_proj


    !gibo: new subroutine for multi-proj cDFT
    subroutine hubbard_internal_radial_multi_proj(hd_species)

      !==================================================================!
      ! This module builds the radial part of hydrogenic wavefunctions   !
      ! in reciprocal space.                                             !
      !------------------------------------------------------------------!
      ! Written by David D. O'Regan in April 2009.                       !
      ! Modified for multiple-projector cDFT by G. Teobaldi in Oct. 2012 !
      !==================================================================!
      ! This subroutine follows the same logic as                        !
      ! hubbard_internal_radial_proj: it that finds the smallest Gmax    !
      ! sufficient to recover normalisation of the radial part of the    !
      ! hydrogenic wavefunctions (HW). The only difference is that,      !
      ! since Gmax depend on the angular momentum of the given HW        !
      ! [the radial extent of the HW depends on their angular momentum:  !
      ! s>p>d>f], care is needed  to select (and return) the largest Gmax!
      ! sufficient to allow normalisation of all the HWs                 !
      ! Note also that, to save troubles, the HWs have been initialised  !
      ! to zero on line 4368 (i.e. they are zero beyong Gmax...)         !
      !------------------------------------------------------------------!

      use comms, only: pub_on_root
      use constants, only: PI
      use hubbard_init, only: HUBBARD_SPECIES
      use utils, only: utils_abort

      implicit none

      ! Arguments
      type(HUBBARD_SPECIES), intent(in) :: hd_species

      ! Local variables
      integer :: tmp_ang_mom                 ! temporary variable for angular moment channels
      integer :: n_g_rad_pts_initial         ! initial value of n_g_rad_pts
      integer :: g_rad_pt                    ! Grid point counters
      real(kind=DP) :: gnorm                 ! Norm in recip space
      real(kind=DP) :: delg                  ! Grid variables
      real(kind=DP) :: gradius, denom, gath
      real(kind=DP) :: proj_max
      real(kind=DP) :: gradius_tmp           ! temporary variable for gradius (depends on ang.mom)
      real(kind=DP) :: gmax_tmp              ! temporary variable for Gmax (depends on ang.mom)


      ! keep track of initial n_g_rad_pts value
      n_g_rad_pts_initial = n_g_rad_pts

      ! dummy initialisation of n_g_rad_pts
      !(we target the largest value to have all the HWs normalised)
      n_g_rad_pts = 0

      ! dummy initialisation of giradius
      !(we target the largest value to have all the HWs normalised)
      gradius = 0.0_DP
      gmax    = 0.0_DP

      !==============================================================================
      !gibo: backwards loop to mantain the s,p,d,f order for tmp_ang_mom (below)
      MULTI_PROJ: do jj = hd_species%num_ang_mom_shells, 1, -1
         !==============================================================================

         ! ddor: Set parameters for radial grid
         delg = PI / rmax
         proj_max = 0.0_DP
         gnorm = 0.0_DP
         gmax_tmp = 0.0_DP

         !gibo: istantanous angular momentum channel
         tmp_ang_mom = hd_species%orb_ang_mom(jj)

         !gath = ( hd_species%hub_charge )**( 0.5_DP*( &
         !     &real(2*hd_species%hub_ang_mom + 5,kind=DP) ) )
         gath = ( hd_species%hub_charge )**( 0.5_DP*( &
              &real(2*tmp_ang_mom + 5,kind=DP) ) )

         ! ddor: Construct radial component in real space
         !do g_rad_pt=1,n_g_rad_pts
         !gibo: mind dependency on agular momentum for n_g_rad_pts below
         RADIAL: do g_rad_pt=1,n_g_rad_pts_initial

            !gradius =  real(g_rad_pt-1,DP) * delg
            gradius_tmp =  real(g_rad_pt-1,DP) * delg

            !denom = ( (hd_species%hub_charge)**2.0_DP + &
            !     (gradius * real(hd_species%hub_ang_mom + 1))**2.0_DP )**&
            !     (-1.0_DP*real(hd_species%hub_ang_mom) - 2.0_DP)
            denom = ( (hd_species%hub_charge)**2.0_DP + &
                 (gradius_tmp * real(tmp_ang_mom + 1))**2.0_DP )**&
                 (-1.0_DP*real(tmp_ang_mom) - 2.0_DP)

            !select case (hd_species%hub_ang_mom)
            select case (tmp_ang_mom)

            case(0) ! 1s

               !hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
               !     SQRT(32.0_DP / PI) * gath * denom
               hub%projectors%rad_proj_recip(g_rad_pt,jj,ps) = &
                    SQRT(32.0_DP / PI) * gath * denom

            case(1) ! 2p

               !hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
               !     (128.0_DP / SQRT(3.0_DP*PI) ) * gradius * gath * denom
               hub%projectors%rad_proj_recip(g_rad_pt,jj,ps) = &
                    (128.0_DP / SQRT(3.0_DP*PI) ) * gradius_tmp * gath * denom

            case(2) ! 3d

               !hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
               !     1728.0_DP * SQRT( 3.0_DP / (5.0_DP*PI) ) * &
               !     &(gradius**2.0_DP) * gath * denom
               hub%projectors%rad_proj_recip(g_rad_pt,jj,ps) = &
                    1728.0_DP * SQRT( 3.0_DP / (5.0_DP*PI) ) * &
                    &(gradius_tmp**2.0_DP) * gath * denom

            case(3) ! 4f

               !hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
               !     131072.0_DP * SQRT( 2.0_DP / (35.0_DP*PI) ) * &
               !     &(gradius**3.0_DP) * gath * denom
               hub%projectors%rad_proj_recip(g_rad_pt,jj,ps) = &
                    131072.0_DP * SQRT( 2.0_DP / (35.0_DP*PI) ) * &
                    &(gradius_tmp**3.0_DP) * gath * denom

            case default
               call utils_abort('Radial or angular quantum number too high in &
                    &hubbard_internal_radial_multi_proj. Only up to r=4,l=3 &
                    &atomic radial functions are currently supported.')
            end select

            !gnorm = gnorm + delg*(gradius * &
            !     hub%projectors%rad_proj_recip(g_rad_pt,1,ps))**2.0_DP
            gnorm = gnorm + delg*(gradius_tmp * &
                 hub%projectors%rad_proj_recip(g_rad_pt,jj,ps))**2.0_DP

            ! ddor: Set core radius if magnitude of
            !       projector is less then 0.01% of
            !       maximum value and cumulative norm is at least 99.99%.
            !proj_max = max(proj_max,abs(hub%projectors%rad_proj_recip(g_rad_pt,1,ps)))
            proj_max = max(proj_max,abs(hub%projectors%rad_proj_recip(g_rad_pt,jj,ps)))
            if (g_rad_pt > 1) then
               !if ( (abs(hub%projectors%rad_proj_recip(g_rad_pt,1,ps)/proj_max) < &
               if ( (abs(hub%projectors%rad_proj_recip(g_rad_pt,jj,ps)/proj_max) < &
                    &1.0e-4_DP) .and. (gnorm > 0.9999_DP) )then

                  !gmax = gradius ! The maximum necessary G
                  gmax_tmp = gradius_tmp ! The maximum necessary G
                  if (gradius_tmp > gmax) gmax = gradius_tmp !gibo: keep only the largest

                  !n_g_rad_pts = g_rad_pt
                  ! make sure to return the largest n_g_rad_pts value
                  if (g_rad_pt > n_g_rad_pts) n_g_rad_pts = g_rad_pt

                  exit RADIAL

               end if
            end if

         end do RADIAL

         ! dhpt: problem arises if n_g_rad_pts is still 0 here
         ! set it back to initial value, ask gibo what is supposed to happen
         if (n_g_rad_pts == 0) n_g_rad_pts = n_g_rad_pts_initial

         ! ddor: Sanity check on Fourier-Bessel transform
         !if  ( (abs(hub%projectors%rad_proj_recip(g_rad_pt,1,ps)/proj_max) > &
         if  ( (abs(hub%projectors%rad_proj_recip(n_g_rad_pts,jj,ps)/proj_max) > &
              &1.0e-4_DP) .or. (gnorm < 0.9999_DP) ) then
            gmax = real(n_g_rad_pts-1,DP) * delg
            if (pub_on_root) write(stdout,*) 'Fourier-Bessel normalisation',&
                 ' not reached in', &
                 ' hubbard_internal_get_radial_multi_proj for ang.mom channel',jj, &
                 ' of angular quantum number', tmp_ang_mom,'. Normalisation',&
                 ' attained was.',gnorm*100,'%',' gmax set to',gmax
         end if

         !==============================================================================
         ! gibo: loop over angular_momentum channels
      enddo MULTI_PROJ
      !==============================================================================

    end subroutine hubbard_internal_radial_multi_proj


  end subroutine hubbard_species_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_species_exit_proj(hub)

    !==========================================================!
    ! This subroutine deallocates memory allocated for         !
    ! data on projecors.                                       !
    !----------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/1/2004.           !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                !
    ! Modified for DFT+U module by David O'Regan in July 2008  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use projectors, only: projectors_deallocate_set

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub

    call projectors_deallocate_set(hub%projectors)

  end subroutine hubbard_species_exit_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_calculate_forces(hub_forces, &
       ngwfs_on_grid,ngwf_basis,hub_proj_basis,hub,cell,fftbox,pur_denskern, &
       hub_overlap,hub_overlap_t,kpt) ! agrecokpt

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the Hubbard DFT+U correction.                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !    hub_forces      : output : DFT+U forces                              !
    !    ngwfs_on_grid   : input  : all of the tight boxes                    !
    !    ngwf_basis      : input  : all of the ngwf spheres                   !
    !    hub_proj_basis  : input  : all of the projector spheres              !
    !    pur_denskern    : input  : purified density kernel SPAM3             !
    !-------------------------------------------------------------------------!
    ! Written by David D. O'Regan in September 2009 based on                  !
    ! pseudo_nl_calculate_forces of the ODG.                                  !
    ! Modified by Nicholas Hine in November 2011 to use the                   !
    ! HUBBARD_MODEL type                                                      !
    !=========================================================================!

    use datatypes, only: FUNCTIONS
    use comms, only: comms_barrier, comms_reduce, pub_my_proc_id
    use fft_box, only: FFTBOX_INFO
    use function_basis, only: FUNC_BASIS
    ! agrecokpt
    use geometry, only: POINT
    use integrals, only : integrals_grad
    use ion, only: ELEMENT
    use parallel_strategy, only: par=>pub_par
    use projectors, only: projectors_func_grad_ovlp_box
    use rundat, only: pub_task, pub_hubbard_restart, pub_hubbard_atomsolve, &
         pub_hub_on_the_fly, pub_cdft, pub_cdft_group_charge_up_only, &
         pub_cdft_group_charge_down_only, pub_cdft_hubbard, &
         pub_num_spins, pub_hub_tensor_forces, pub_imag_thr
    use simulation_cell, only: CELL_INFO
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_get_element, sparse_product, sparse_transpose, &
         sparse_scale, sparse_copy, sparse_take_real_part
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: hub_forces(1:3,par%nat)
    type(FUNC_BASIS) :: ngwf_basis
    type(FUNCTIONS), intent(in) :: ngwfs_on_grid
    type(FUNC_BASIS) :: hub_proj_basis
    type(HUBBARD_MODEL) :: hub
    type(FFTBOX_INFO), intent(in) :: fftbox
    type(CELL_INFO), intent(in) :: cell
    type(SPAM3), intent(in) :: pur_denskern(pub_num_spins)
    type(SPAM3), intent(in) :: hub_overlap
    type(SPAM3), intent(in) :: hub_overlap_t
    ! agrecokpt
    type(POINT), optional, intent(in) :: kpt

    ! Local Variables
    type(SPAM3) :: siGp_overlap(3)
    type(SPAM3) :: hub_ps_overlap
    type(SPAM3) :: piGp_overlap(3)    !ddor: for tensor force
    !type(SPAM3) :: piGp_overlap_local !ddor: for tensor force
    type(SPAM3) :: sppiGp_overlap(3), sppiGp_overlap_buffer
    type(SPAM3) :: hub_force_mat(3)
    type(SPAM3) :: hk, hkiG
    integer :: is
    integer :: cart
    integer :: iat, orig_iat
    integer :: atom_proj, global_proj
    real(kind=DP) :: hub_force
    ! agrecocmplx
    logical :: loc_cmplx
    ! agrecocmplx: temp matrices
    type(SPAM3) :: hub_proj_cmplx
    type(SPAM3) :: hub_force_mat_real
    ! agrecokpt
    type(POINT) :: loc_kpt

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Entering hubbard_calculate_forces'

    ! Start timer
    call timer_clock('hubbard_calculate_forces',1)

    ! agrecocmplx
    loc_cmplx = ngwfs_on_grid%iscmplx

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
          call utils_abort('Error in hubbard_calculate_forces: &
               &non-Gamma k-point requires complex projectors.')
       end if
    end if

    ! Initialise
    hub_forces  = 0.0_DP

    ! ddor: Create matrix to hold DFT+U forces
    ! ddor: This may truncate, but only the diagonal is of interest
    hub_force_mat(:)%structure = 'G'
    do cart=1,3
       call sparse_create(hub_force_mat(cart))
    end do

    ! ddor: G = U/2 [ 1 - 2 (Occupancy matrix in projector representation) ]
    !       Call it hubbard_term

    ! Create matrices to hold G^_ <Proj^i|NGWF_a>
    call sparse_create(hub_ps_overlap,hub%projector_ham(1),hub_overlap_t)

    ! Create matrices to hold <phi|iG*proj> overlap matrices in each direction
    do cart=1,3
       call sparse_create(siGp_overlap(cart),hub_overlap)
       if (pub_hub_tensor_forces) then
          ! ddor: <proj_|iG.proj_> has sparsity O__ or  O^^
          call sparse_create(piGp_overlap(cart),hub%o_matrix)
          ! ddor: <phi|proj_> O^^ <proj_|iG.proj_> has sparsity if
          !       <phi|proj_> O^^ O__ has sparsity of
          !       <phi|proj_> O^^
          call sparse_create(sppiGp_overlap(cart),hub_overlap,&
               hub%o_matrix)
       endif
    end do

    ! Calculate <phi|iG*proj_> overlap matrix
    if (((pub_task == 'HUBBARDSCF') .and. (hub%consistency_iteration > 1)) &
         & .or. pub_hubbard_restart .or. pub_hubbard_atomsolve .or. &
         & pub_hub_on_the_fly) then
       call integrals_grad(siGp_overlap,ngwfs_on_grid,ngwf_basis,&
            hub%consistency_projs,hub_proj_basis,cell,fftbox)
       ! ddor: compute gradient of covariant subspace metric if needed
       if (pub_hub_tensor_forces) then
          call integrals_grad(piGp_overlap,hub%consistency_projs,hub_proj_basis,&
               hub%consistency_projs,hub_proj_basis,cell,fftbox)
          ! ddor: compute temporary transpose of hub_overlap_t
          call sparse_create(sppiGp_overlap_buffer,sppiGp_overlap(1))
          call sparse_transpose(sppiGp_overlap_buffer,hub_overlap_t)
          ! ddor: compute NGWF overlap with gradient matrix
          do cart=1,3
             ! ddor: First remove the local (atom-diagonal) part
             !call sparse_copy(piGp_overlap_local,piGp_overlap(cart))
             !call sparse_axpy(piGp_overlap(cart),piGp_overlap_local,-1.0_DP)
             ! ddor: now overlap on the left with NGWFs
             call sparse_product(sppiGp_overlap(cart),&
                  sppiGp_overlap_buffer,piGp_overlap(cart))
          enddo
          call sparse_destroy(sppiGp_overlap_buffer)
       endif
    else
       do cart=1,3
          ! agrecokpt: need to specify k-point here? or not necessary?
          call projectors_func_grad_ovlp_box(siGp_overlap(cart), &
               ngwfs_on_grid,ngwf_basis,hub_proj_basis,hub%projectors, &
               fftbox,cell,cart) !kpt=loc_kpt
       end do
    end if

    ! agrecocmplx
    if (loc_cmplx) then
       ! create temp cmplx matrix
       call sparse_create(hub_proj_cmplx, hub%projector_ham(1), &
            iscmplx=loc_cmplx)
       ! create real version for temp copy back
       call sparse_create(hub_force_mat_real, hub_force_mat(1))
    end if

    spinloop: do is= 1, pub_num_spins

       !gibo: for cdft_group_charge_up_only .OR. cdft_group_charge_down_only,
       !       make sure \nabla Tr[n] of the ignored spin is not added to forces
       !gibo: for cdft_group_charge_up_only, skip \nable Tr[n_down]
       if (pub_cdft .AND. (.not.pub_cdft_hubbard) .AND. &
            pub_cdft_group_charge_up_only .AND. (is==2)) CYCLE spinloop

       !gibo: for cdft_group_charge_down_only, skip \nable Tr[n_up]
       if (pub_cdft .AND. (.not.pub_cdft_hubbard) .AND. &
            pub_cdft_group_charge_down_only .AND. (is==1)) CYCLE spinloop

       ! Calculate the matrix H^_ <Proj^i|NGWF_a>
       if (loc_cmplx) then
          ! agrecocmplx: copy to temp cmplx matrix before taking product
          call sparse_copy(hub_proj_cmplx, hub%projector_ham(is))
          call sparse_product(hub_ps_overlap,hub_proj_cmplx,hub_overlap_t) !W
       else
          call sparse_product(hub_ps_overlap,hub%projector_ham(is),hub_overlap_t) !W
       end if

       ! Create temporary matrices hk and hkiG
       call sparse_create(hk,hub_ps_overlap,pur_denskern(is))
       call sparse_create(hkiG,hub_force_mat(1),iscmplx=loc_cmplx)

       ! Compute hk : H^_ <proj^j|phi_b> K^ab
       call sparse_product(hk,hub_ps_overlap,pur_denskern(is))

       ! Loop over Cartesian co-ordinates
       ! Force mat has projector tensor format F^m_m'
       do cart=1,3
          ! Calculate H^_ <proj^j|phi_b> K^ab <phi_a|iG.proj_i>
          call sparse_product(hkiG,hk,siGp_overlap(cart))
          if (loc_cmplx) then
             ! jme: take real part and accumulate it
             ! agrecocmplx: safe conversion of complex to real
             !call sparse_copy(hub_force_mat_real, hkiG)
             call sparse_take_real_part(hub_force_mat_real, hkiG, pub_imag_thr)
             call sparse_axpy(hub_force_mat(cart),hub_force_mat_real,2.0_DP)
          else
             call sparse_axpy(hub_force_mat(cart),hkiG,2.0_DP)
                ! ddor: one minus due to the definition of force
                ! ddor: one minus due to |iG.proj_i> overlapping
                ! ddor: with rho rather than H, for ease of coding
                ! ddor: Factor of two for complex conjugate
          end if
       end do

       ! ddor: compute gradient of covariant subspace metric if needed
       if ((((pub_task == 'HUBBARDSCF') .and. (hub%consistency_iteration > 1)) &
            & .or. pub_hubbard_restart .or. pub_hubbard_atomsolve .or. &
            & pub_hub_on_the_fly) .and. (pub_hub_tensor_forces)) then

          ! Loop over Cartesian co-ordinates
          ! Force mat has projector tensor format F^m_m'
          do cart=1,3
             ! ddor: Calculate - H^_ <proj^j|phi_b> K^ab <phi_a|proj^k> <proj_k|iG.proj_i>
             ! ddor: Here we use d/dR [O^^] = - O^^ d/dR [O__] O^^
             call sparse_product(hkiG,hk,sppiGp_overlap(cart))
             if (loc_cmplx) then
                ! agrecocmplx: save conversion of complex to real
                !call sparse_copy(hub_force_mat_real, hkiG)
                call sparse_take_real_part(hub_force_mat_real, hkiG, pub_imag_thr)
                call sparse_axpy(hub_force_mat(cart),hub_force_mat_real,-2.0_DP)
             else
                call sparse_axpy(hub_force_mat(cart),hkiG,-2.0_DP)
                ! ddor: one minus due to the definition of force
                ! ddor: one minus due to the negative exponent of the overlap
                ! ddor: one minus due to |proj^k> <proj_k|iG.proj_i> overlapping
                ! ddor: with rho rather than H, for ease of coding
                ! ddor: Factor of two for complex conjugate
             end if
          end do

       endif

       ! Destroy temporary matrices
       call sparse_destroy(hkiG)
       call sparse_destroy(hk)

    end do spinloop

    ! Destroy temporary matrices
    if (loc_cmplx) then
       call sparse_destroy(hub_force_mat_real)
       call sparse_destroy(hub_proj_cmplx)
    end if
    !if (pub_hub_tensor_forces) then
    !   call sparse_destroy(piGp_overlap_local)
    !endif
    do cart=3,1,-1
       if (pub_hub_tensor_forces) then
          call sparse_destroy(sppiGp_overlap(cart))
          call sparse_destroy(piGp_overlap(cart))
       endif
       call sparse_destroy(siGp_overlap(cart))
    end do
    call sparse_destroy(hub_ps_overlap)

    ! ddor: Force depends only on centre of Hubbard projector, that
    !       is the position of the Hubbard atom. There is no contribution
    !       to the force on atoms overlapping with Hubbard projectors.

    ! Loop over atoms
    do iat=par%first_atom_on_proc(pub_my_proc_id), &
         par%first_atom_on_proc(pub_my_proc_id + 1) - 1

       ! Find atom number in input file order
       orig_iat = par%orig_atom(iat)

       ! Loop over Hubbard projectors on this atom
       do atom_proj=1,hub_proj_basis%num_on_atom(iat)
          global_proj = hub_proj_basis%first_on_atom(iat) + atom_proj - 1

          ! Loop over Cartesian co-ordinates
          do cart=1,3

             ! Find contribution of this projector to force on this atom
             ! from diagonal elements of hub_force_mat for this coordinate
             call sparse_get_element(hub_force,hub_force_mat(cart), &
                  global_proj, global_proj)
             hub_forces(cart,orig_iat) = hub_forces(cart,orig_iat) + hub_force

          end do  ! cart

       end do  ! atom_proj

    end do   ! loc_iat

    ! Reduce result across procs
    call comms_barrier
    call comms_reduce('SUM',hub_forces)

    ! Destroy temporary matrices
    do cart=1,3
       call sparse_destroy(hub_force_mat(cart))
    end do

    ! Stop timer
    call timer_clock('hubbard_calculate_forces',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving hubbard_calculate_forces'

  end subroutine hubbard_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! ebl: the following subroutines are required by the DMFT-related
  ! alterations to this module

  ! ebl: Calculate the rotated occupancy matrix
  subroutine hubbard_diago_occupancy(occupancy, chan, hat, rot_vec_angle, &
      tmp_dens)

    use linalg, only: linalg_dsyev_lt, linalg_mat_mul_serial

    implicit none

    ! Arguments
    real(kind=DP),           intent(in   ) :: occupancy(:,:)
    integer,                 intent(in   ) :: chan
    integer,                 intent(in   ) :: hat
    real(kind=DP), optional, intent(in   ) :: rot_vec_angle(3,3)
    real(kind=DP), optional, intent(inout) :: tmp_dens(chan,chan)

    ! Local variables
    real(kind=DP)               :: trace
    real(kind=DP)               :: occupancy_rot(chan,chan)
    real(kind=DP)               :: tmp(chan,chan)
    real(kind=DP)               :: diagdens(chan)
    real(kind=DP)               :: rotation(chan,chan)
    logical                     :: read_in_file
    integer                     :: i

    read_in_file=.false.

    if (maxval(abs(occupancy))<1.e-8_DP) then
       return
    endif

    if (present(rot_vec_angle)) then
       read_in_file=maxval(abs(rot_vec_angle))>1.e-4_DP
    endif

    ! ebl: find the rotation that diagonalises the occupancy matrix
    if (.not.read_in_file) then
       rotation = (occupancy+transpose(occupancy))/2.0_DP
       call linalg_dsyev_lt(rotation, diagdens, chan)
    else
       rotation=hubbard_cmp_rot(rot_vec_angle,(chan-1)/2,chan,hat)
    endif

    ! ebl: rotate the occupancy matrix
    call linalg_mat_mul_serial(tmp, rotation, occupancy, opA="T")
    call linalg_mat_mul_serial(occupancy_rot, tmp, rotation)

    ! ebl: optionally return the rotated occupancy matrix
    if(present(tmp_dens)) tmp_dens=occupancy_rot

    ! ebl: calculate trace and diagonal elements
    trace=0.0_DP
    do i=1,chan
       trace=trace+occupancy(i,i)
       diagdens(i)=occupancy_rot(i,i)
    enddo

    ! ebl: printing of the results is currently suppressed, may return to this
    !      in the future
    ! write(*,*) '================================================='
    ! call write_array(occupancy_rot,' DIAGONAL DENSITIES ')
    ! write(*,*) 'TOTAL (trace rotated matrix)  = ', sum(diagdens)
    ! write(*,*) 'CHECK (trace original matrix) = ', occupancy
    ! write(*,*) '================================================='

    ! occupancy_rot=occupancy_rot/sum(diagdens)
    ! entropy=0.d0
    ! do i=1,chan
    !    n=2.0*occupancy_rot(i,i)
    !    d=n**2/4.0
    !    l(1)= 1.0 - n  + d
    !    l(2)= n/2.0 -d
    !    l(3)= n/2.0 -d
    !    l(4)= d
    !    write(*,'(a,i3,2f15.3)') 'ORBITAL / SUM P / ENT : ', i, sum(l), -sum(l*log(l))
    !    entropy = entropy + sum(l*log(l))
    ! enddo

    ! write(*,*) '================================================='
    ! write(*,*) 'ENTROPY 4 MANY BODY STATES (U=0)  = ', -entropy
    ! write(*,*) '================================================='

    ! entropy = 1.d0 - sum(diag(matmul(occupancy_rot,occupancy_rot)))
    ! write(*,*) '================================================='
    ! write(*,*) 'NORMALIZATION                   = ', sum(diag(occupancy_rot))
    ! write(*,*) 'LINEAR ENTROPY                  = ', entropy
    ! write(*,*) 'WITH SPIN FACTOR 2              = ', 2.d0*entropy
    ! write(*,*) '================================================='

    ! occupancy_rot=occupancy_rot*2.d0
    ! occupancy_rot=occupancy_rot/sum(diag(occupancy_rot))
    ! entropy = 1.d0 - sum(diag(matmul(occupancy_rot,occupancy_rot)))
    ! write(*,*) '================================================='
    ! write(*,*) 'NORMALIZATION                   = ', sum(diag(occupancy_rot))
    ! write(*,*) 'LINEAR ENTROPY FULL RHO (WITH S)= ', entropy
    ! write(*,*) '================================================='

  end subroutine

  function hubbard_atoms_occupancy(hub_proj_basis, hub, hat, spin)

    use function_basis,    only: FUNC_BASIS
    use sparse,            only: sparse_get_block
    use parallel_strategy, only: par=>pub_par
    use utils, only: utils_assert, utils_int_to_str

    implicit none

    ! Arguments
    type(FUNC_BASIS),    intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(in) :: hub
    integer,             intent(in) :: hat
    integer,             intent(in) :: spin

    ! Local variables
    integer                    :: n_on_at
    real(kind=DP)              :: h_atoms_occupancy_tmp( &
         hub_proj_basis%max_on_atom, hub_proj_basis%max_on_atom)
    real(kind=DP), allocatable :: hubbard_atoms_occupancy(:, :)

    n_on_at = hub_proj_basis%num_on_atom(par%distr_atom(hub%orig(hat)))

    allocate(hubbard_atoms_occupancy(n_on_at, n_on_at))

    call sparse_get_block(h_atoms_occupancy_tmp,hub%occupancy_matrix(spin),&
         par%distr_atom(hub%orig(hat)),par%distr_atom(hub%orig(hat)))

    call utils_assert(n_on_at<=hub_proj_basis%max_on_atom, &
         'hubbard_atoms_occupancy: inconsistent with number of orbitals '//&
         'hub_proj_basis%max_on_atom: '//&
         trim(utils_int_to_str(hub_proj_basis%max_on_atom))//&
         ', n_on_at: '//&
         trim(utils_int_to_str(n_on_at)))

    hubbard_atoms_occupancy(1:n_on_at,1:n_on_at) = h_atoms_occupancy_tmp(1:n_on_at,1:n_on_at)

  end function

  function hubbard_cmp_rot(rot_vec, l, ndim, hat)

    use utils, only: utils_assert, utils_abort, utils_unit

    implicit none

    ! Arguments
    real(kind=DP), intent(in)  :: rot_vec(3,3)
    integer,       intent(in) :: l
    integer,       intent(in) :: ndim
    integer,       intent(in) :: hat

    ! Local variables
    integer                   :: i, j, aa(7)
    complex(kind=DP)          :: hubbard_cmp_rot(ndim, ndim), cmp_tmp(ndim, ndim)
    real(kind=DP)             :: Rd(5, 5), Rf(7, 7), Rp(3, 3)
    integer                   :: ii, bb(7), iost
    logical                   :: schemerot, check_ngwfs
    integer                   :: unit

    !-----------------------------------------------------------------!
    ! transformed G/sigma is obtained by G' = ( cmp_rot^T G cmp_rot ) !
    !-----------------------------------------------------------------!

    schemerot=.false.

    inquire(file='mask_ngwfs',exist=check_ngwfs)

    if(check_ngwfs) then
       unit = utils_unit()
       open(unit=unit,file='mask_ngwfs')
       do ii=1,hat-1
          read(unit,*)
       end do
       read(unit,*,IOSTAT=iost) (bb(ii),ii=1,ndim)
       close(unit)

       call utils_assert(iost >= 0, 'hubbard_cmp_rot: '//&
            'NOT ENOUGH LINES IN MASK_NGWFS')
    else
       bb=(/( ii, ii=1,ndim )/)
    endif

    if (l==1) then
       !===========================!
       ! ONETEP IS WITH :
       ! -1: py
       !  0: pz
       ! +1: px
       !===========================!
       !==============================!
       ! Rp-language
       ! b1 : y
       ! b2 : z
       ! b3 : x
       !==============================!
       call rotate_p(Rp, rot_vec)
       aa(1)=1
       aa(2)=2
       aa(3)=3
       do i=1,3
          do j=1,3
             hubbard_cmp_rot(i,j)=Rp(aa(i),aa(j))
          enddo
       enddo
       call utils_assert(ndim == 3, 'hubbard_cmp_rot: ndim not matching, L=1')
    elseif (l==2) then
       !===========================!
       ! ONETEP IS WITH :
       ! 1 IS m=-2 : dxy
       ! 2 IS m=-1 : dzy
       ! 3 IS m= 0 : d3z2-r2
       ! 4 IS m= 1 : dxz
       ! 5 IS m= 2 : dx2-y2
       !===========================!
       ! THIS ROUTINE WITH NOTATIONS :
       !     yz
       !     zx
       !     xy
       !   x2-y2
       !   3z2-1
       !===========================!
       if(schemerot)then
          call rotate_d(Rd,rot_vec)
          aa(1)=3
          aa(2)=1
          aa(3)=5
          aa(4)=2
          aa(5)=4
          do i=1,5
             do j=1,5
                hubbard_cmp_rot(i,j)=Rd(aa(i),aa(j))
             enddo
          enddo
       else
          call rotate_our_notation(Rd,rot_vec)
          hubbard_cmp_rot=Rd
       endif
       call utils_assert(ndim == 5, 'hubbard_cmp_rot: ndim not matching, L=2')
    elseif (l==3) then
       call rotate_f(Rf, rot_vec)
       !onetep i is j in Rf-language
       aa(1)=1
       aa(2)=2
       aa(3)=3
       aa(4)=4
       aa(5)=5
       aa(6)=6
       aa(7)=7
       !==============================!
       ! Onetep
       ! b1 fy(3x^2-y^2)
       ! b2 fxyz
       ! b3 f y z^2
       ! b4 f   z^3
       ! b5 f x z^2
       ! b6 fz(x2-y2)
       ! b7 fx(x^2-3y^2)
       !==============================!
       ! Rf-language
       ! b1 fy(3x2-y2)
       ! b2 xyz
       ! b3 y(4z2 -x2   -y2)
       ! b4 z(2z2 -3x2 -3y2)
       ! b5 x(4z2 -x2 - y2)
       ! b6 z(x2 - y2)
       ! b7 x(x2 - 3y2)
       !==============================!
       do i=1,7
          do j=1,7
             hubbard_cmp_rot(i,j)=Rf(aa(i),aa(j))
          enddo
       enddo
       call utils_assert(ndim == 7, 'hubbard_cmp_rot: ndim not matching, L=3')
    else
       call utils_abort('hubbard_cmp_rot: '//&
            'rotation for angular momentum l/=2 not implemented')
    endif

    cmp_tmp=hubbard_cmp_rot
    do i=1,ndim
       do j=1,ndim
          hubbard_cmp_rot(i,j)=cmp_tmp(bb(i),bb(j))
       enddo
    enddo

  end function

  subroutine rotate_our_notation(Rd,R)

    implicit none

    real(kind=DP), intent(in)  :: R(3,3)
    real(kind=DP), intent(out) :: Rd(5,5)
    real(kind=DP)              :: s3,r1,r2,r3,r4,r5,r6,r7,r8,r9

    s3=sqrt(3.0_DP)

    r1=R(1,1)
    r2=R(1,2)
    r3=R(1,3)
    r4=R(2,1)
    r5=R(2,2)
    r6=R(2,3)
    r7=R(3,1)
    r8=R(3,2)
    r9=R(3,3)

    Rd(1,1)=r2*r4+r1*r5
    Rd(1,2)=r3*r5+r2*r6
    Rd(1,3)=s3*r3*r6
    Rd(1,4)=r3*r4+r1*r6
    Rd(1,5)=2.0_DP*r1*r4+r3*r6

    Rd(2,1)=r5*r7+r4*r8
    Rd(2,2)=r6*r8+r5*r9
    Rd(2,3)=s3*r6*r9
    Rd(2,4)=r6*r7+r4*r9
    Rd(2,5)= 2.0_DP*r4*r7+r6*r9

    Rd(3,1)= s3*r7*r8
    Rd(3,2)= s3*r8*r9
    Rd(3,3)= (3.* r9**2 -1.d0)/2.0_DP
    Rd(3,4)= s3*r7*r9
    Rd(3,5)= s3*(2.*r7**2 + r9**2 - 1.0_DP)/2.0_DP

    Rd(4,1)= r2*r7+r1*r8
    Rd(4,2)=  r3*r8+r2*r9
    Rd(4,3)=  s3*r3*r9
    Rd(4,4)= r3*r7+r1*r9
    Rd(4,5)= 2.*r1*r7 +r3*r9

    Rd(5,1)=r1*r2-r4*r5
    Rd(5,2)=r2*r3-r5*r6
    Rd(5,3)=s3*(r3**2-r6**2)/2.0_DP
    Rd(5,4)=r1*r3-r4*r6
    Rd(5,5)=(2.*r1**2 + r3**2 - 2.*r4**2 - r6**2)/2.0_DP

  end subroutine

  subroutine rotate_p(Rp, R)
    implicit none
    real(kind=DP) :: R(3,3),Rp(3,3)
    real(kind=DP) :: r1,r2,r3,r4,r5,r6,r7,r8,r9

    ! Rotation of the p electron cubic harmonic orbitals. Real space
    ! rotation matrix is R(3,3). The cubic wave functions are defined as follows:

    !  m=-1, y
    !  m= 0, z
    !  m= 1, x

    r1=R(1,1)
    r2=R(1,2)
    r3=R(1,3)
    r4=R(2,1)
    r5=R(2,2)
    r6=R(2,3)
    r7=R(3,1)
    r8=R(3,2)
    r9=R(3,3)

    Rp(1,1)= r5
    Rp(1,2)= r6
    Rp(1,3)= r4
    Rp(2,1)= r8
    Rp(2,2)= r9
    Rp(2,3)= r7
    Rp(3,1)= r2
    Rp(3,2)= r3
    Rp(3,3)= r1

  end subroutine

  subroutine rotate_d(Rd, R)

  ! Rotation of the d electron cubic harmonic orbitals. Real space
  ! rotation matrix is R(3,3). The cubic wave functions are defined as follows:

  ! b1 = 2yz/r^2 * a
  ! b2 = 2xz/r^2 * a
  ! b3 = 2xy/r^2 * a
  ! b4 = (x^2-y^2)/r^2 * a
  ! b5 = (3z^2-r^2)/r^2 * a

  ! Rotation changes them as follows

  ! R.b1 =   (R(2,2)*R(3,3) + R(2,3)*R(3,2))*b1 +
  !          (R(2,1)*R(3,3) + R(2,3)*R(3,1))*b2
  !        + (R(2,1)*R(3,2) + R(2,2)*R(3,1))*b3
  !        + (R(2,1)*R(3,1) - R(2,2)*R(3,2))*b4 +
  !          sqrt(3.)*R(2,3)*R(3,3)*b5

    implicit none
    real(kind=DP), intent(in)  :: R(3,3)
    real(kind=DP), intent(out) :: Rd(5,5)

    Rd(1,1) = R(2,2)*R(3,3) + R(2,3)*R(3,2);
    Rd(1,2) = R(2,1)*R(3,3) + R(2,3)*R(3,1);
    Rd(1,3) = R(2,1)*R(3,2) + R(2,2)*R(3,1);
    Rd(1,4) = R(2,1)*R(3,1) - R(2,2)*R(3,2);
    Rd(1,5) = sqrt(3.)*R(2,3)*R(3,3);

    Rd(2,1) = R(1,2)*R(3,3) + R(1,3)*R(3,2);
    Rd(2,2) = R(1,1)*R(3,3) + R(1,3)*R(3,1);
    Rd(2,3) = R(1,1)*R(3,2) + R(1,2)*R(3,1);
    Rd(2,4) = R(1,1)*R(3,1) - R(1,2)*R(3,2);
    Rd(2,5) = sqrt(3.)*R(1,3)*R(3,3);

    Rd(3,1) = R(1,2)*R(2,3) + R(1,3)*R(2,2);
    Rd(3,2) = R(1,1)*R(2,3) + R(1,3)*R(2,1);
    Rd(3,3) = R(1,1)*R(2,2) + R(1,2)*R(2,1);
    Rd(3,4) = R(1,1)*R(2,1) - R(1,2)*R(2,2);
    Rd(3,5) = sqrt(3.)*R(1,3)*R(2,3);

    Rd(4,1) = R(1,2)*R(1,3) - R(2,2)*R(2,3);
    Rd(4,2) = R(1,1)*R(1,3) - R(2,1)*R(2,3);
    Rd(4,3) = R(1,1)*R(1,2) - R(2,1)*R(2,2);
    Rd(4,4) = 0.5*(R(1,1)**2+R(2,2)**2) -0.5*(R(1,2)**2 + R(2,1)**2);
    Rd(4,5) = 0.5*sqrt(3.)*(R(1,3)**2-R(2,3)**2);

    Rd(5,1) = sqrt(3.)*R(3,2)*R(3,3);
    Rd(5,2) = sqrt(3.)*R(3,1)*R(3,3);
    Rd(5,3) = sqrt(3.)*R(3,1)*R(3,2);
    Rd(5,4) = 0.5*sqrt(3.)*(R(3,1)**2-R(3,2)**2);
    Rd(5,5) = -0.5 + 1.5*R(3,3)**2;

  end subroutine

  subroutine rotate_f(Rf, R)
    implicit none
    real(kind=DP) :: R(3,3),Rf(7,7)
    real(kind=DP) :: a(7),b(7),c(7),d(7),e(7),f(7),g(7)
    real(kind=DP) :: r1,r2,r3,r4,r5,r6,r7,r8,r9

    ! Rotation of the f electron cubic harmonic orbitals. Real space
    ! rotation matrix is R(3,3). The cubic wave functions are defined as follows:

    ! b1 m=-3, fy(3x^2-y^2)
    ! b2 m=-2, fxyz
    ! b3 m=-1, fyz^2
    ! b4 m= 0, fz^3
    ! b5 m= 1, fxz^2
    ! b6 m= 2, fz(x2-y2)
    ! b7 m= 3, fx(x^2-3y^2)

    r1=R(1,1)
    r2=R(1,2)
    r3=R(1,3)
    r4=R(2,1)
    r5=R(2,2)
    r6=R(2,3)
    r7=R(3,1)
    r8=R(3,2)
    r9=R(3,3)

    a(1) = ( 2.d0*r2*(4.d0*r1*r4+r3*r6)+4.d0*r5*( r1**2 - r4**2 )+ r5 * ( r3**2 - r6**2 ) ) / 4.d0
    a(2) = ( 4.d0*(r2*r4*r7 + r1*r5*r7 + r1*r4*r8) + r3*r6*r8 + r3*r5*r9 +r2*r6*r9 )/sqrt(6.d0)
    a(3) = sqrt(5.d0/3.d0) * (2.d0*r8*(4.d0*r4*r7 + r6*r9) + r5*(4.d0*r7**2 + r9**2 - 1.d0) ) /4.d0
    a(4) = sqrt(5.d0/2.d0) * r8 * (4.d0*r7**2 + r9**2 - 1.d0 ) / 2.d0
    a(5) = sqrt(5.d0/3.d0) * ( 2.d0*r8 * (4.d0*r1*r7+r3*r9 ) + r2*(4.d0*r7**2 +r9**2 - 1.d0) ) /4.d0
    a(6) = ( 8.d0*r7*( r1*r2 - r4*r5) + 4.d0*r8*(r1**2-r4**2) + r8*(r3**2-r6**2) + 2.d0*r9*(r2*r3 - r5*r6) ) / ( 2.d0*sqrt(6.d0) )
    a(7) = ( r2*(r3**2-4.d0*r4**2 -r6**2 +4.d0*r1**2) - 2.d0*r5*( 4.d0*r1*r4 + r3*r6 ) ) /4.d0

    b(1) = sqrt(3.d0/2.d0)* ( r3*(r2*r4 + r1*r5) + r6*(r1*r2 - r4*r5) )
    b(2) = (r3*r5 + r2*r6)*r7 + (r3*r4 + r1*r6)*r8 +(r2*r4 + r1*r5)*r9
    b(3) = sqrt(5.d0/2.d0) * ( r6*r7*r8 + r5*r7*r9 + r4*r8*r9 )
    b(4) = sqrt(15.d0)*r7*r8*r9
    b(5) = sqrt(5.d0/2.d0)*(r3*r7*r8 + r2*r7*r9 + r1*r8*r9 )
    b(6) = r7*(r2*r3 - r5*r6) + r8*(r1*r3 - r4*r6) + r9*(r1*r2 - r4*r5)
    b(7) = sqrt(3.d0/2.d0)*( r3*(r1*r2 - r4*r5) - r6*(r2*r4 + r1*r5) )

    c(1) = sqrt(15.d0)*( r3**2*r5 + 2.d0*r2*r3*r6 - r5*r6**2) / 4.d0
    c(2) = sqrt(5.d0/2.d0) * ( r3*r6*r8 + r3*r5*r9 + r2*r6*r9 )
    c(3) = (10.d0*r6*r8*r9 + r5*(5.d0*r9**2 - 1.d0) ) / 4.d0
    c(4) = sqrt(3.d0/2.d0)*r8*( 5.d0*r9**2 - 1.d0 ) / 2.d0
    c(5) = ( 10.d0*r3*r8*r9 + r2*(5.d0*r9**2 - 1.d0) ) / 4.d0
    c(6) = sqrt(5.d0/2.d0)* ( r3*(r3*r8+2.d0*r2*r9) - r6*(r6*r8 + 2.d0*r5*r9) ) / 2.d0
    c(7) = sqrt(15.d0)*( r2*(r3**2 - r6**2) - 2.d0*r3*r5*r6 ) / 4.d0

    d(1) = -sqrt(5.d0/2.d0)*r6*(r6**2 - 3.d0*r3**2) / 2.d0
    d(2) =  sqrt(15.d0) * r3*r6*r9
    d(3) =  sqrt(3.d0/2.d0)* r6 * (5.d0*r9**2 - 1.d0) / 2.d0
    d(4) =  r9 * (5.d0*r9**2 - 3.d0) / 2.d0
    d(5) =  sqrt(3.d0/2.d0) * r3 * (5.d0*r9**2 - 1.d0) / 2.d0
    d(6) =  sqrt(15.d0)*r9*(r3**2 - r6**2) / 2.d0
    d(7) =  sqrt(5.d0/2.d0) * r3 * ( r3**2 - 3.d0*r6**2 ) / 2.d0

    e(1) =  sqrt(15.d0)*(r3**2*r4 + 2.d0*r1*r3*r6 - r4*r6**2) / 4.d0
    e(2) =  sqrt(5.d0/2.d0) * (r3*r6*r7 + r3*r4*r9 + r1*r6*r9)
    e(3) =  (10.d0 * r6*r7*r9 + r4*(5.d0*r9**2 - 1.d0))/4.d0
    e(4) =  sqrt(3.d0/2.d0)*r7*( 5.d0*r9**2 - 1.d0 ) / 2.d0
    e(5) =  ( 10.d0*r3*r7*r9 + r1*(5.d0*r9**2 - 1.d0) ) / 4.d0
    e(6) =  sqrt(5.d0/2.d0) * ( r3**2*r7 + 2.d0*r1*r3*r9 - r6*(r6*r7 + 2.d0*r4*r9 ) ) / 2.d0
    e(7) =  sqrt(15.d0) * ( r1*(r3**2 - r6**2 ) - 2.d0*r3*r4*r6 ) / 4.d0

    f(1) =  sqrt(3.d0/2.d0) * ( 4.d0*r1*r3*r4 + r6*(2.d0*r1**2 - 2.d0*r4**2 - r6**2 + 3.d0*r3**2) ) / 2.d0
    f(2) =  2.d0*r7*(r3*r4 + r1*r6) + r9*(2.d0*r1*r4 + 3.d0*r3*r6)
    f(3) =  sqrt(5.d0/2.d0) * ( 4.d0*r4*r7*r9 + r6*(2.d0*r7**2 + 3.d0*r9**2 - 1.d0) ) / 2.d0
    f(4) =  sqrt(15.d0) * r9 * ( 2.d0*r7**2 + r9**2 - 1.d0 ) / 2.d0
    f(5) =  sqrt(5.d0/2.d0) * ( 4.d0*r1*r7*r9 + r3*(2.d0*r7**2 + 3.d0*r9**2 - 1.d0) ) / 2.d0
    f(6) =  2.d0*r7*( r1*r3 - r4*r6 ) + 3.d0*r9*( r3**2 - r6**2 )/2.d0 + r9*(r1**2 - r4**2 )
    f(7) =  sqrt(3.d0/2.d0) * ( r3*(2.d0*r1**2 + r3**2 - 2.d0*r4**2 - 3.d0*r6**2 ) - 4.d0*r1*r4*r6 ) / 2.d0

    g(1) =  ( r4*(4.d0*r5**2 + r6**2 - 4.d0*r2**2 - r3**2 ) - r1*( 8.d0*r2*r5 + 2.d0*r3*r6 ) ) / 4.d0
    g(2) = -(r3*(r6*r7 + r4*r9) + 4.d0*r2*( r5*r7 + r4*r8 ) + r1*( r6*r9 + 4.d0*r5*r8 ) ) / sqrt(6.d0)
    g(3) = -sqrt(5.d0/3.d0) * ( 2.d0*r7*( 4.d0*r5*r8 + r6*r9 ) + r4*(4.d0*r8**2 + r9**2 - 1.d0) ) / 4.d0
    g(4) = -sqrt(5.d0/2.d0) * r7 * ( 4.d0*r8**2 + r9**2 - 1.d0 ) / 2.d0
    g(5) = -sqrt(5.d0/3.d0) * (2.d0*r7*(4.d0*r2*r8 + r3*r9 ) + r1*(4.d0*r8**2 + r9**2 - 1.d0) ) /4.d0
    g(6) =  ( r7 * ( r6**2 - 4.d0*r2**2 - r3**2 + 4.d0*r5**2 ) - 8.d0*r8*( r1*r2 - r4*r5 ) - 2.d0*r9*( r1*r3 - r4*r6 ) ) &
       / sqrt(24.d0)
    g(7) =  ( 2.d0*r4 * (4.d0*r2*r5 + r3*r6 ) + r1*( r6**2 - 4.d0*r2**2 -r3**2 + 4.d0*r5**2 ) ) / 4.d0

    Rf(:,1)=a
    Rf(:,2)=b
    Rf(:,3)=c
    Rf(:,4)=d
    Rf(:,5)=e
    Rf(:,6)=f
    Rf(:,7)=g

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module hubbard_build
