! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   This file was created in August 2000 by Chris-Jriton Skylaris.
!
!   The subroutines in this file were written by Chris-Kriton Skylaris,
!   Arash A. Mostofi, Jacek Dziedzic, Nicholas Hine, Victor Milman,
!   Quintin Hill, and Jose M. Escartin.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module simulation_cell

  use geometry, only: point
  use constants, only: DP

  implicit none

  private

  type CELL_INFO
     ! cks: lattice vectors
     type(POINT)      :: a1, a2, a3
     ! cks: lattice vectors normalised to unit length
     type(POINT)      :: a1_unit, a2_unit, a3_unit
     ! cks: reciprocal lattice vectors
     type(POINT)      :: b1, b2, b3
     ! cks: regular grid spacing in atomic units
     real(kind=DP)   :: d1, d2, d3
     ! cks: weight for integrals on the standard grid
     real(kind=DP)   :: weight
     ! cks: plane-wave kinetic energy cutoff (in Eh)
     real(kind=DP)   :: ke_cutoff
     ! cks: total number of ppds in the simulation cell
     integer          :: n_ppds
     ! cks: total number of ppds in each lattice vector direction
     integer          :: n_ppds_a1,n_ppds_a2,n_ppds_a3
     ! cks: number of points per ppd on the standard grid
     integer          :: n_pts
     ! cks: number of points per ppd in each lattice vector dir on standard grid
     integer          :: n_pt1,n_pt2,n_pt3
     ! cks: total number of points in each lattice direction on the standard grid
     integer          :: total_pt1,total_pt2,total_pt3
  end type CELL_INFO

  ! Pointers in castep_cell and castep_model converted to allocatables
  ! by JM Escartin on 11/10/2016.

  type castep_cell                                         ! CASTEP cell definition
     real(kind=dp), dimension(3,3) :: real_lattice
     real(kind=dp), dimension(3,3) :: recip_lattice
     real(kind=dp) :: volume
     real(kind=dp), dimension(:,:,:), allocatable :: ionic_positions
     integer :: num_ions
     integer :: num_species
     integer :: max_ions_in_species
     integer, dimension(:), allocatable :: num_ions_in_species
     real(kind=dp), dimension(:), allocatable :: ionic_charge  ! pa: changed from integer to real
     real(kind=dp), dimension(:), allocatable :: species_mass
     character(len=2), dimension(:), allocatable:: species_symbol
  end type castep_cell


  type castep_model                                        ! CASTEP model definition
     type(castep_cell) :: cell
     type(castep_cell) :: orig_cell
     type(castep_cell) :: ref_cell
     type(castep_cell) :: reac_cell ! Reactant geometry for transition state runs
     type(castep_cell) :: prod_cell ! Product geometry for transition state runs
     type(castep_cell) :: intm_cell ! Intermediate geometry for transition state runs
     real(kind=dp) :: total_energy
     real(kind=dp), dimension(:,:,:), allocatable :: forces
     integer :: bfgs_iteration
     integer :: lbfgs_iteration
     integer :: tpsd_iteration
     real(kind=dp), dimension(:,:), allocatable :: bfgs_inv_hessian
     real(kind=dp), dimension(:,:), allocatable :: lbfgs_position_updates
     real(kind=dp), dimension(:,:), allocatable :: lbfgs_gradient_updates
     real(kind=dp) :: tpsd_alpha
     integer :: lbfgs_num_updates
     logical :: bfgs_optimisation
     logical :: lbfgs_optimisation
     logical :: tpsd_optimisation
     logical :: found_ground_state
     real(kind=dp), dimension(6) :: stress   ! stress tensor
     real(kind=dp), dimension(3,3) :: strain ! strain tensor
     ! aam: in presence of constraints we keep a copy of the uncontrained forces
     real(kind=dp), dimension(:,:,:), allocatable :: orig_forces
  end type castep_model

  public :: CELL_INFO
  public :: castep_cell
  public :: castep_model
  public :: simulation_cell_initialise
  public :: simulation_cell_add_padding
  public :: simulation_cell_num_periodic_images
  public :: simulation_cell_lattice_points_init
  public :: simulation_cell_lattice_points_exit
  public :: castep_cell_copy
  public :: castep_cell_alloc
  public :: castep_cell_dealloc
  public :: castep_cell_cart_to_frac
  public :: castep_cell_frac_to_cart
  public :: castep_model_alloc
  public :: castep_model_dealloc
  public :: castep_model_cell_changed
  public :: castep_cell_cart_lattice_to_abc
  public :: castep_cell_2_elements
  public :: copy_to_castep_cell
  public :: minimum_image_distance
  public :: minimum_image_distance_ortho
  public :: minimum_image_displacement
  public :: minimum_image_number_to_indices
  public :: minimum_image_number_to_displacement

  ! ndmh: Cutoff Coulomb extra padding required
  real(kind=DP), parameter, public :: cutoff_coulomb_tol = 10.0_DP

  !--------------------------------------------------------------------------!
  ! Overload subroutines
  !--------------------------------------------------------------------------!

  interface castep_cell_frac_to_cart
     module procedure castep_cell_frac_to_cart_cell
     module procedure castep_cell_frac_to_cart_vector
  end interface


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulation_cell_add_padding(d1,d2,d3,n_pt1,n_pt2,n_pt3, &
       a1,a2,a3,a1_pad,a2_pad,a3_pad)

    !================================================================!
    ! This subroutine checks that the padded cell is commensurate    !
    ! with the original cell, or sets up the padding if it is left   !
    ! to the default.                                                !
    !----------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                    !
    ! Modified by Gabriel Bramley in July 2019.                      !
    !================================================================!

    use constants, only: stdout
    use geometry, only: POINT, operator(+), operator(*), geometry_magnitude
    use rundat, only: pub_coulomb_radius, pub_coulomb_length, &
         pub_coulomb_cutoff_type, pub_debug_on_root
    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: d1,d2,d3
    integer, intent(in) :: n_pt1, n_pt2, n_pt3
    type(POINT), intent(in) :: a1, a2, a3
    type(POINT), intent(inout) :: a1_pad, a2_pad, a3_pad

    ! Local Variables
    type(POINT) :: da1, da2, da3
    real(kind=DP) :: len, len_pad
    real(kind=DP) :: cutoff

    cutoff = sqrt(pub_coulomb_radius**2 + pub_coulomb_length**2)

    ! ndmh: check if the default is still set, if so, override with a suitable
    ! ndmh: padded cell
    if ((a1_pad%x==-1.0_DP).and.(a1_pad%y== 0.0_DP).and.(a1_pad%z== 0.0_DP).and. &
        (a2_pad%x== 0.0_DP).and.(a2_pad%y==-1.0_DP).and.(a2_pad%z== 0.0_DP).and. &
        (a3_pad%x== 0.0_DP).and.(a3_pad%y== 0.0_DP).and.(a3_pad%z==-1.0_DP)) then

       ! Copy the original cell
       a1_pad = a1
       a2_pad = a2
       a3_pad = a3

       ! Increase a1 by one ppd at a time till we are large enough
       ! (only if 0D periodicity)
       if ((pub_coulomb_cutoff_type=='SPHERE') .or. &
            (pub_coulomb_cutoff_type=='CYLINDER')) then
          len = geometry_magnitude(a1)
          da1 = (real(n_pt1,kind=DP)/(len/d1)) * a1
          do
              len_pad = geometry_magnitude(a1_pad)
              if (len_pad >= len + cutoff + cutoff_coulomb_tol) exit
              a1_pad = a1_pad + da1
          end do
       end if

       ! Now increase a2 by one ppd at a time till we are large enough
       ! (only if 0D or 1D periodicity)
       if ((pub_coulomb_cutoff_type=='SPHERE') .or. &
            (pub_coulomb_cutoff_type=='CYLINDER').or. &
            (pub_coulomb_cutoff_type=='WIRE')) then
          len = geometry_magnitude(a2)
          da2 = (real(n_pt2,kind=DP)/(len/d2)) * a2
          do
              len_pad = geometry_magnitude(a2_pad)
              if (len_pad >= len + cutoff + cutoff_coulomb_tol) exit
              a2_pad = a2_pad + da2
          end do
       end if

       ! Now increase a3 by one ppd at a time till we are large enough
       ! (any periodicity)
       len = geometry_magnitude(a3)
       da3 = (real(n_pt3,kind=DP)/(len/d3)) * a3
       if ((pub_coulomb_cutoff_type=='SPHERE') .or. &
            (pub_coulomb_cutoff_type=='CYLINDER').or. &
            (pub_coulomb_cutoff_type=='WIRE')) then
          do
             len_pad = geometry_magnitude(a3_pad)
             if (len_pad >= len + cutoff + cutoff_coulomb_tol) exit
             a3_pad = a3_pad + da3
          end do

       else if (pub_coulomb_cutoff_type=='SLAB') then  ! gab: Slab case.
          do
             len_pad = geometry_magnitude(a3_pad)
             if ( (2*len) - len_pad <= 0.000001_DP ) exit ! pad a3 = cell a3.
             a3_pad = a3_pad + da3
          end do

       end if ! End slab case

    end if

    if (pub_debug_on_root) then
       write(stdout,*) 'DEBUG: Padded unit cell vectors:'
       write(stdout,*) a1_pad
       write(stdout,*) a2_pad
       write(stdout,*) a3_pad
       write(stdout,*)
    end if


    ! ndmh: Check for badly set-up padded cells
    if (geometry_magnitude(a1_pad)<geometry_magnitude(a1)) then
        call utils_abort('ERROR in simulation_cell_add_padding: Padded cell &
             &is smaller than original cell along a1')
    end if
    if (geometry_magnitude(a2_pad)<geometry_magnitude(a2)) then
        call utils_abort('ERROR in simulation_cell_add_padding: Padded cell &
             &is smaller than original cell along a2')
    end if
    if (geometry_magnitude(a3_pad)<geometry_magnitude(a3)) then
        call utils_abort('ERROR in simulation_cell_add_padding: Padded cell &
             &is smaller than original cell along a3')
    end if

  end subroutine simulation_cell_add_padding


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine simulation_cell_initialise(cell,a1,a2,a3,n_pt1,n_pt2,n_pt3,d1,d2,d3)

    !============================================================!
    ! This subroutine initialises the components of a CELL_INFO  !
    ! type which groups together information about the           !
    ! simulation cell.                                           !
    !------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 21/6/2001.             !
    ! Changed to accept cell as an argument to allow setting up  !
    ! of padded cell for cutoff coulomb by Nick Hine 01/04/2008  !
    !============================================================!

    use constants, only: DP, PI
    use geometry, only: point, unit_vector, operator(.cross.), operator(*),&
        operator(.dot.), geometry_magnitude

    implicit none

    type(CELL_INFO), intent(out) :: cell
    type(POINT), intent(in) :: a1, a2, a3
    real(kind=DP), intent(in) :: d1, d2, d3
    integer, intent(in) :: n_pt1, n_pt2, n_pt3

    ! cks: initilise grid spacing
    cell%d1 = d1
    cell%d2 = d2
    cell%d3 = d3

    ! cks: initialise primitive lattice vectors
    cell%a1 = a1
    cell%a2 = a2
    cell%a3 = a3

    ! cks: initialise number of points in ppd in each lattice vector direction
    cell%n_pt1 = n_pt1
    cell%n_pt2 = n_pt2
    cell%n_pt3 = n_pt3

    ! cks: initialise primitive lattice vectors with unit length
    cell%a1_unit = UNIT_VECTOR(cell%a1)
    cell%a2_unit = UNIT_VECTOR(cell%a2)
    cell%a3_unit = UNIT_VECTOR(cell%a3)

    ! cks: initialise reciprocal lattice vectors
    cell%b1=(2.0_DP*PI/( cell%a1.DOT.(cell%a2.CROSS.cell%a3) ) ) * &
         (cell%a2.CROSS.cell%a3)
    cell%b2=(2.0_DP*PI/( cell%a1.DOT.(cell%a2.CROSS.cell%a3) ) ) * &
         (cell%a3.CROSS.cell%a1)
    cell%b3=(2.0_DP*PI/( cell%a1.DOT.(cell%a2.CROSS.cell%a3) ) ) * &
         (cell%a1.CROSS.cell%a2)

    ! cks: initialise grid point weights
    cell%weight=abs(cell%d1*cell%d2*cell%d3 * &
         ( (cell%a1_unit.CROSS.cell%a2_unit).DOT.cell%a3_unit ))

    ! cks: initialise number of ppds in each lattice vector direction
    cell%n_ppds_a1=nint(geometry_MAGNITUDE(cell%a1) / &
         ( real(cell%n_pt1,DP)*cell%d1) )
    cell%n_ppds_a2=nint(geometry_MAGNITUDE(cell%a2) / &
         ( real(cell%n_pt2,DP)*cell%d2) )
    cell%n_ppds_a3=nint(geometry_MAGNITUDE(cell%a3) / &
         ( real(cell%n_pt3,DP)*cell%d3) )
    cell%n_ppds=cell%n_ppds_a1*cell%n_ppds_a2*cell%n_ppds_a3

    ! cks: initialise total number of points per ppd
    cell%n_pts = n_pt1 * n_pt2 * n_pt3

    ! cks: initialise total number of points in simulation cell in each
    !      lattice vector direction
    cell%total_pt1 = cell%n_ppds_a1 * cell%n_pt1
    cell%total_pt2 = cell%n_ppds_a2 * cell%n_pt2
    cell%total_pt3 = cell%n_ppds_a3 * cell%n_pt3

  end subroutine simulation_cell_initialise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_copy(in_cell,out_cell)
    !=========================================================================!
    ! Copy CASTEP cell data                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    ! Adapted to allocatable components by JM Escartin, 11/10/2016.           !
    !=========================================================================!

    use utils, only: utils_alloc_check, utils_assert

    implicit none

    type(castep_cell), intent(inout) :: out_cell
    type(castep_cell), intent(in) :: in_cell

    ! <<< local variables >>>
    integer :: i,j,ierr

    ! Check that 'in_cell' pointers are allocated - they should be! - abort if not.
    call utils_assert(&
         allocated(in_cell%ionic_positions) .and. &
         allocated(in_cell%num_ions_in_species) .and. &
         allocated(in_cell%ionic_charge) .and. &
         allocated(in_cell%species_mass) .and. &
         allocated(in_cell%species_symbol), &
         'Internal error in castep_cell_copy()')

    ! Real and reciprocal space lattices
    do i=1,3
       do j=1,3
          out_cell%real_lattice(i,j)  = in_cell%real_lattice(i,j)
          out_cell%recip_lattice(i,j) = in_cell%recip_lattice(i,j)
       enddo
    enddo

    ! Cell volume
    out_cell%volume = in_cell%volume

    ! Check that 'out_cell' pointers are allocated - they should be! - allocate if not.
    if (.not.allocated(out_cell%ionic_positions)) then
       allocate(out_cell%ionic_positions(1:3,1:1,1:in_cell%num_ions),stat=ierr)
       call utils_alloc_check('castep_cell_copy','out_cell%ionic_positions', ierr)
    endif
    if (.not.allocated(out_cell%num_ions_in_species)) then
       allocate(out_cell%num_ions_in_species(1:in_cell%num_ions),stat=ierr)
       call utils_alloc_check('castep_cell_copy','out_cell%num_ions_in_species', ierr)
    endif
    if (.not.allocated(out_cell%ionic_charge)) then
       allocate(out_cell%ionic_charge(1:in_cell%num_ions),stat=ierr)
       call utils_alloc_check('castep_cell_copy','out_cell%ionic_charge', ierr)
    endif
    if (.not.allocated(out_cell%species_mass)) then
       allocate(out_cell%species_mass(1:in_cell%num_ions),stat=ierr)
       call utils_alloc_check('castep_cell_copy','out_cell%species_mass', ierr)
    endif
    if (.not.allocated(out_cell%species_symbol)) then
       allocate(out_cell%species_symbol(1:in_cell%num_ions),stat=ierr)
       call utils_alloc_check('castep_cell_copy','out_cell%species_symbol', ierr)
    endif

    ! Ion data
    out_cell%num_ions            = in_cell%num_ions
    out_cell%num_species         = in_cell%num_species
    out_cell%max_ions_in_species = in_cell%max_ions_in_species
    do i=1,in_cell%num_ions
       out_cell%ionic_positions(1,1,i) = in_cell%ionic_positions(1,1,i)
       out_cell%ionic_positions(2,1,i) = in_cell%ionic_positions(2,1,i)
       out_cell%ionic_positions(3,1,i) = in_cell%ionic_positions(3,1,i)
       out_cell%num_ions_in_species(i) = in_cell%num_ions_in_species(i)
       out_cell%ionic_charge(i)        = in_cell%ionic_charge(i)
       out_cell%species_mass(i)        = in_cell%species_mass(i)
       out_cell%species_symbol(i)      = in_cell%species_symbol(i)
    end do

  end subroutine castep_cell_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_frac_to_cart_vector(current_cell,frac,cart)
    !=========================================================================!
    ! Obtain absolute Cartesian co-ordinates from fractional ones             !
    !-------------------------------------------------------------------------!
    ! Arash A Mostofi, 2004                                                   !
    !=========================================================================!

    implicit none

    type(castep_cell), intent(in)            :: current_cell
    real(kind=dp), dimension(3), intent(in)  :: frac
    real(kind=dp), dimension(3), intent(out) :: cart

    ! <<< local variables >>>
    integer :: ii

    do ii=1,3
       cart(ii) = current_cell%real_lattice(1,ii)*frac(1) &
            + current_cell%real_lattice(2,ii)*frac(2) &
            + current_cell%real_lattice(3,ii)*frac(3)
    end do

    return
  end subroutine castep_cell_frac_to_cart_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_frac_to_cart_cell(current_cell,cart)
    !==========================================================================!
    ! Obtain absolute Cartesian co-ordinates from fractional ones              !
    ! This one acts on the whole array of coordinates as stored in cell object !
    !--------------------------------------------------------------------------!
    ! Victor Milman 21 Apr 2006                                                !
    !==========================================================================!

    implicit none

    type(castep_cell), intent(in)                :: current_cell
    real(kind=dp), dimension(:,:,:), intent(out) :: cart

    ! <<< local variables >>>
    integer :: ispec,iatom

    do ispec = 1,current_cell%num_species
       do iatom = 1,current_cell%num_ions_in_species(ispec)
          call castep_cell_frac_to_cart(current_cell, &
                current_cell%ionic_positions(:,iatom,ispec),cart(:,iatom,ispec))
       enddo
    enddo

    return
  end subroutine castep_cell_frac_to_cart_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_cart_to_frac(current_cell,cart,frac)
    !=========================================================================!
    ! Obtain fractional coordinates from Cartesian ones                       !
    !-------------------------------------------------------------------------!
    ! Arash A Mostofi, 2004                                                   !
    !=========================================================================!

    use constants, only: two_pi
    implicit none

    type(castep_cell), intent(in)           :: current_cell
    real(kind=dp), dimension(3), intent(in) :: cart
    real(kind=dp), dimension(3), intent(out):: frac

    ! <<< local variables >>>
    integer :: ii

    do ii=1,3
       frac(ii)=current_cell%recip_lattice(ii,1)*cart(1) &
            + current_cell%recip_lattice(ii,2)*cart(2) &
            + current_cell%recip_lattice(ii,3)*cart(3)
    end do

    frac(:) = frac(:) / two_pi

    return
  end subroutine castep_cell_cart_to_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_model_cell_changed(mdl)
    !=========================================================================!
    ! Sets castep_model data when ionic positions change                      !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 2004                                        !
    !=========================================================================!

    implicit none

    type(castep_model), intent(inout) :: mdl

    mdl%found_ground_state = .false.
!!$    mdl%found_forces = .false.
!!$    mdl%bfgs_optimisation = .false.

    return
  end subroutine castep_model_cell_changed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_dealloc(cell)
    !=========================================================================!
    ! Dellocate CASTEP cell data                                              !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    !=========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    type(castep_cell), intent(inout) :: cell

    ! <<< local variables >>>
    integer :: ierr

    deallocate(cell%species_symbol,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%species_symbol', ierr)
    deallocate(cell%species_mass,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%species_mass', ierr)
    deallocate(cell%ionic_positions,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%ionic_positions', ierr)
    deallocate(cell%num_ions_in_species,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%num_ions_in_species', ierr)
    deallocate(cell%ionic_charge,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%ionic_charge', ierr)

  end subroutine castep_cell_dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_alloc(cell,nat)
    !=========================================================================!
    ! Allocate CASTEP cell data                                               !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    ! Adapted to allocatable components by JM Escartin, 11/10/2016.           !
    !=========================================================================!

    use utils, only: utils_alloc_check, utils_abort

    implicit none

    integer, intent(in) :: nat
    type(castep_cell), intent(inout) :: cell

    ! <<< local variables >>>
    integer :: ierr

    if  ((allocated(cell%ionic_positions)) .or. &
         (allocated(cell%num_ions_in_species)) .or. &
         (allocated(cell%ionic_charge)) .or. &
         (allocated(cell%species_mass)) .or. &
         (allocated(cell%species_symbol))) then
       call utils_abort('Internal error in castep_cell_alloc(): cell already &
            &allocated')
    end if

    ! For simplicity each atom has its own species
    allocate(cell%ionic_positions(1:3,1:1,1:nat),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%ionic_positions', ierr)
    allocate(cell%num_ions_in_species(1:nat),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%num_ions_in_species', ierr)
    allocate(cell%ionic_charge(1:nat),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%ionic_charge', ierr)
    allocate(cell%species_mass(1:nat),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%species_mass', ierr)
    allocate(cell%species_symbol(1:nat),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%species_symbol', ierr)

    return
  end subroutine castep_cell_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_model_alloc(mdl,nat,ndim,constrain_ions,do_lbfgs)
    !=========================================================================!
    ! Allocate CASTEP model data                                              !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    !=========================================================================!

    use comms, only: pub_on_root
    use utils, only: utils_alloc_check
    use rundat, only: pub_geom_lbfgs_block_length

    implicit none

    integer, intent(in)               :: nat
    type(castep_model), intent(inout) :: mdl
    integer, intent(in)               :: ndim
    logical, intent(in)               :: constrain_ions
    logical, intent(in)               :: do_lbfgs

    ! <<< local variables >>>
    integer :: length
    integer :: ierr

    ! Allocate mdl%cell
    call castep_cell_alloc(mdl%cell,nat)

    ! Allocate mdl%orig_cell
    call castep_cell_alloc(mdl%orig_cell,nat)

    ! Allocate mdl%ref_cell
    call castep_cell_alloc(mdl%ref_cell,nat)

    ! Allocate mdl%reac_cell
    call castep_cell_alloc(mdl%reac_cell,nat)

    ! Allocate mdl%prod_cell
    call castep_cell_alloc(mdl%prod_cell,nat)

    ! Allocate mdl%intm_cell
    call castep_cell_alloc(mdl%intm_cell,nat)

    ! Allocate mdl%forces
    ! For simplicity each atom has its own species
    allocate(mdl%forces(1:3,1:1,1:nat),stat=ierr)
    call utils_alloc_check('castep_model_alloc','mdl%forces', ierr)

    ! Allocate mdl%orig_forces if constrain_ions is TRUE
    ! For simplicity each atom has its own species
    if (constrain_ions) then
       allocate(mdl%orig_forces(1:3,1:1,1:nat),stat=ierr)
       call utils_alloc_check('castep_model_alloc','mdl%orig_forces', ierr)
    endif

    ! Allocate mdl%bfgs_inv_Hessian
    if (((.not.do_lbfgs).or.(mdl%bfgs_optimisation)).and.pub_on_root) then
       allocate(mdl%bfgs_inv_Hessian(1:ndim,1:ndim),stat=ierr)
       call utils_alloc_check('castep_model_alloc','mdl%bfgs_inv_Hessian', ierr)
    end if

!    mdl%lbfgs_optimisation = .false.
    if ((do_lbfgs.or.mdl%lbfgs_optimisation).and.pub_on_root) then
!       if(geom_lbfgs_max_updates.eq.0) then
!          length=geom_lbfgs_block_length+mdl%lbfgs_num_updates
!       else
!          length=geom_lbfgs_block_length
!       end if
       length=max(mdl%lbfgs_num_updates,pub_geom_lbfgs_block_length) ! might be too large, but this will be sorted out after...
       allocate(mdl%lbfgs_position_updates(1:ndim,1:length),stat=ierr)
       call utils_alloc_check('castep_model_alloc','mdl%lbfgs_position_updates', ierr)
       allocate(mdl%lbfgs_gradient_updates(1:ndim,1:length),stat=ierr)
       call utils_alloc_check('castep_model_alloc','mdl%lbfgs_gradient_updates', ierr)
    end if

    return
  end subroutine castep_model_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_model_dealloc(mdl,do_lbfgs)
    !=========================================================================!
    ! Deallocate CASTEP model data                                            !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    ! Adapted to allocatable components by JM Escartin, 11/10/2016.           !
    !=========================================================================!

    use comms, only : pub_on_root
    use utils, only : utils_dealloc_check

    implicit none

    type(castep_model), intent(inout) :: mdl
    logical :: do_lbfgs

    ! <<< local variables >>>
    integer :: ierr

    ! Deallocate mdl%bfgs_inv_Hessian
    if ((.not.do_lbfgs).and.pub_on_root) then
       deallocate(mdl%bfgs_inv_Hessian,stat=ierr)
       call utils_dealloc_check('castep_model_dealloc','mdl%bfgs_inv_Hessian', ierr)
    end if

    if (do_lbfgs.and.pub_on_root) then
       deallocate(mdl%lbfgs_position_updates,stat=ierr)
       call utils_dealloc_check('castep_model_dealloc','mdl%lbfgs_position_updates', ierr)
       deallocate(mdl%lbfgs_gradient_updates,stat=ierr)
       call utils_dealloc_check('castep_model_dealloc','mdl%lbfgs_gradient_updates', ierr)
    end if

    ! Dellocate mdl%orig_forces if constrain_ions is TRUE
    if (allocated(mdl%orig_forces)) then
       deallocate(mdl%orig_forces,stat=ierr)
       call utils_dealloc_check('castep_model_dealloc','mdl%orig_forces', ierr)
    endif

    ! Dellocate mdl%forces
    deallocate(mdl%forces,stat=ierr)
    call utils_dealloc_check('castep_model_dealloc','mdl%forces', ierr)

    ! Deallocate mdl%intm_cell
    call castep_cell_dealloc(mdl%intm_cell)

    ! Deallocate mdl%prod_cell
    call castep_cell_dealloc(mdl%prod_cell)

    ! Deallocate mdl%reac_cell
    call castep_cell_dealloc(mdl%reac_cell)

    ! Deallocate mdl%ref_cell
    call castep_cell_dealloc(mdl%ref_cell)

    ! Deallocate mdl%orig_cell
    call castep_cell_dealloc(mdl%orig_cell)

    ! Deallocate mdl%cell
    call castep_cell_dealloc(mdl%cell)

    return
  end subroutine castep_model_dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_cart_lattice_to_abc(cell,a,b,c,alpha,beta,gamma)

    !--------------------------------------------------------------------------!
    ! Purpose:                                                                 !
    ! This routine generates the lattice vectors of `cell' and                 !
    ! stores them in a,b,c,alpha,beta,gamma                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! cell(in): cell from which real_lattice will be used                      !
    ! a(out): a lattice parameter                                              !
    ! b(out): b lattice parameter                                              !
    ! c(out): c lattice parameter                                              !
    ! alpha(out): angle between b and c                                        !
    ! beta(out) : angle between c and a                                        !
    ! gamma(out): angle between a and b                                        !
    !--------------------------------------------------------------------------!
    ! Parent module variables used: unit_cell                                  !
    !--------------------------------------------------------------------------!
    ! Modules used: none                                                       !
    !--------------------------------------------------------------------------!
    ! Key internal variables: none                                             !
    !--------------------------------------------------------------------------!
    ! Necessary conditions: none                                               !
    !--------------------------------------------------------------------------!
    ! Transferred from CASTEP by Victor Milman 20/04/06                        !
    !--------------------------------------------------------------------------!

    use constants, only : pi

    implicit none

    type(castep_cell), intent(in)::cell    ! The cell
    real(kind=dp), intent(out)::a,b,c,alpha,beta,gamma  ! The lattice parameters

    ! Calculate a
    a=sqrt(cell%real_lattice(1,1)**2+ &
         cell%real_lattice(1,2)**2+ &
         cell%real_lattice(1,3)**2)

    ! Calculate b
    b=sqrt(cell%real_lattice(2,1)**2+ &
         cell%real_lattice(2,2)**2+ &
         cell%real_lattice(2,3)**2)

    ! Calculate c
    c=sqrt(cell%real_lattice(3,1)**2+ &
         cell%real_lattice(3,2)**2+ &
         cell%real_lattice(3,3)**2)

    ! Calculate cos(alpha), temp stored in alpha
    alpha=(cell%real_lattice(2,1)*cell%real_lattice(3,1)+ &
         cell%real_lattice(2,2)*cell%real_lattice(3,2)+ &
         cell%real_lattice(2,3)*cell%real_lattice(3,3))/(b*c)

    ! Take acos and tranfrom into degrees
    alpha=180.0_dp*acos(alpha)/pi

    ! Calculate cos(beta), temp stored in beta
    beta =(cell%real_lattice(3,1)*cell%real_lattice(1,1)+ &
         cell%real_lattice(3,2)*cell%real_lattice(1,2)+ &
         cell%real_lattice(3,3)*cell%real_lattice(1,3))/(c*a)

    ! Take acos and tranfrom into degrees
    beta =180.0_dp*acos(beta)/pi

    ! Calculate cos(gamma), temp stored in gamma
    gamma=(cell%real_lattice(1,1)*cell%real_lattice(2,1)+ &
         cell%real_lattice(1,2)*cell%real_lattice(2,2)+ &
         cell%real_lattice(1,3)*cell%real_lattice(2,3))/(a*b)

    ! Take acos and tranfrom into degrees
    gamma=180.0_dp*acos(gamma)/pi

    return

  end subroutine castep_cell_cart_lattice_to_abc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_2_elements(cell,elements,nat)
    !=========================================================================!
    ! Update elements with new ionic positions from castep cell               !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 2004                                        !
    ! AAM: Modified July 2005 to work with negative fractional co-ordinates   !
    !=========================================================================!

    use ion, only : element
    implicit none

    type(castep_cell), intent(in) :: cell
    integer, intent(in) :: nat
    type(element), intent(inout)  :: elements(1:nat)

    ! <<< local variable >>>
    integer :: atom,k
    real(kind=dp) :: frac(1:3), cart(1:3)

    do atom=1,nat
       frac(1:3) = cell%ionic_positions(1:3,1,atom)
       do k=1,3
          ! rationalise fractional positions to the interval [0,1]
          if ( frac(k).lt.0.0_dp ) then
             frac(k) = frac(k) + 1.0_dp
          else if ( frac(k).gt.1.0_dp ) then
             frac(k) = frac(k) - 1.0_dp
          else
             continue
          endif
       enddo

       call castep_cell_frac_to_cart(cell,frac,cart)

       elements(atom)%centre%x = cart(1)
       elements(atom)%centre%y = cart(2)
       elements(atom)%centre%z = cart(3)
    enddo

    return
  end subroutine castep_cell_2_elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine copy_to_castep_cell(current_cell,cell,nat,elements,nat_classical,classical_elements)
    !=========================================================================!
    ! Create and fill castep cell from ONETEP cell data structures            !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman 27/03/06                                       !
    ! Modified by Chris-Kriton Skylaris to inlcude classical info, 2010/08/30.!
    !=========================================================================!

    use constants, only: periodic_table_name, two_pi, &
         periodic_table_mass
    use geometry, only: operator(.dot.), operator(.cross.)
    use ion, only: element
    use utils, only: utils_assert

    implicit none

    ! Arguments
    type(CELL_INFO), intent(in)      :: cell
    type(castep_cell), intent(inout) :: current_cell
    integer, intent(in)              :: nat
    integer, intent(in)              :: nat_classical
    type(element), intent(in)        :: elements(1:nat)
    type(element), intent(in), optional :: &
         classical_elements(1:nat_classical)

    ! local variables
    integer :: ion_i, ion_k
    logical :: found_ion

    !***********************************************************************!
    !         Copy data from ONETEP structures into CASTEP structures       !
    !***********************************************************************!

    ! Allocate current_cell if necessary
    if (.not.allocated(current_cell%ionic_positions)) then
       call castep_cell_alloc(current_cell,nat+nat_classical)
    end if

    ! Real space lattice
    current_cell%real_lattice(1,1) = cell%a1%x
    current_cell%real_lattice(1,2) = cell%a1%y
    current_cell%real_lattice(1,3) = cell%a1%z
    current_cell%real_lattice(2,1) = cell%a2%x
    current_cell%real_lattice(2,2) = cell%a2%y
    current_cell%real_lattice(2,3) = cell%a2%z
    current_cell%real_lattice(3,1) = cell%a3%x
    current_cell%real_lattice(3,2) = cell%a3%y
    current_cell%real_lattice(3,3) = cell%a3%z

    ! Reciprocal lattice
    current_cell%recip_lattice(1,1) = cell%b1%x
    current_cell%recip_lattice(1,2) = cell%b1%y
    current_cell%recip_lattice(1,3) = cell%b1%z
    current_cell%recip_lattice(2,1) = cell%b2%x
    current_cell%recip_lattice(2,2) = cell%b2%y
    current_cell%recip_lattice(2,3) = cell%b2%z
    current_cell%recip_lattice(3,1) = cell%b3%x
    current_cell%recip_lattice(3,2) = cell%b3%y
    current_cell%recip_lattice(3,3) = cell%b3%z

    ! Cell volume
    current_cell%volume = abs( (cell%a1 .cross. cell%a2) .dot. cell%a3 )

    ! Ion data plus classical ion data
    current_cell%num_ions = nat + nat_classical
    current_cell%num_species = nat + nat_classical
    current_cell%max_ions_in_species = 1 ! Give each ion its own species for now

    ! ONETEP ionic positions are in absolute Cartesian coordinates (in bohr)
    ! CASTEP uses fractional co-ordinates
    do ion_i=1,nat
       current_cell%ionic_positions(1,1,ion_i) = &
            (cell%b1.dot.elements(ion_i)%centre) / two_pi
       current_cell%ionic_positions(2,1,ion_i) = &
            (cell%b2.dot.elements(ion_i)%centre) / two_pi
       current_cell%ionic_positions(3,1,ion_i) = &
            (cell%b3.dot.elements(ion_i)%centre) / two_pi
       current_cell%num_ions_in_species(ion_i) = 1
       current_cell%ionic_charge(ion_i)        = elements(ion_i)%ion_charge
       current_cell%species_symbol(ion_i)      = & ! Transform species_symbol
            species_symbol(elements(ion_i)%symbol) ! characters to standard form
    enddo

    ! cks: do the same for "classical" atoms
    ! ONETEP ionic positions are in absolute Cartesian coordinates (in bohr)
    ! CASTEP uses fractional co-ordinates
    if (present(classical_elements)) then
       do ion_i =nat +1, nat + nat_classical
          current_cell%ionic_positions(1,1,ion_i) = &
               (cell%b1.dot.classical_elements(ion_i-nat)%centre)/&
               two_pi
          current_cell%ionic_positions(2,1,ion_i) = &
               (cell%b2.dot.classical_elements(ion_i-nat)%centre)/&
               two_pi
          current_cell%ionic_positions(3,1,ion_i) = &
               (cell%b3.dot.classical_elements(ion_i-nat)%centre)/&
               two_pi
          current_cell%num_ions_in_species(ion_i) = 1
          current_cell%ionic_charge(ion_i)        = &
               classical_elements(ion_i -nat)%ion_charge
          current_cell%species_symbol(ion_i)      = & ! Transform species_symbol
               species_symbol(classical_elements(ion_i-nat)%symbol)
                                                   ! characters to standard form
       enddo
    endif

    ! Ionic masses
    do ion_i=1,nat
       found_ion =.false.
       periodic_table: do ion_k=1,size(periodic_table_mass)
          if ( current_cell%species_symbol(ion_i).eq.periodic_table_name(ion_k) ) then
             current_cell%species_mass(ion_i) = periodic_table_mass(ion_k)
             found_ion = .true.
             exit periodic_table
          endif
       enddo periodic_table
       call utils_assert(found_ion, 'Internal error in copy_to_castep_cell(): &
            &species_symbol '//current_cell%species_symbol(ion_i)//&
               ' in geometry_optimise not recognised')
    end do
    !*********************************************************************!
    !                       Cell data copy completed                      !
    !*********************************************************************!

  end subroutine copy_to_castep_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function species_symbol(string)
    !=========================================================================!
    ! Transform the species symbol character to a standard form (first letter !
    ! uppercase, second lowercase)                                            !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    character(*), intent(in) :: string
    character(2)             :: species_symbol

    ! <<< local variables >>>
    integer :: iA,iaa,iZ,izz,ishift,ic,ln

    iaa = ichar('a')
    izz = ichar('z')
    iA = ichar('A')
    iZ = ichar('Z')
    ishift = iA-iaa

    ln = len(string)
    if (ln.lt.1 .or. ln.gt.2) then
       call utils_abort('Error detected in species_symbol(): incorrect length &
            &of element symbol:'//trim(string)//'.')
    endif

    species_symbol(1:ln) = string(1:ln)
    if (ln.lt.2) species_symbol(2:2)=' '

    ! Make first character uppercase
    ic = ichar(species_symbol(1:1))
    if ((ic.ge.iaa).and.(ic.le.izz)) species_symbol(1:1) = char(ishift+ic)

    ! Make second character lowercase
    if (species_symbol(2:2).ne.' ') then
       ic = ichar(species_symbol(2:2))
       if ((ic.ge.iA).and.(ic.le.iZ)) species_symbol(2:2) = char(ic-ishift)
    endif

    return
  end function species_symbol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function minimum_image_distance(r1,r2,cell)
    !=========================================================================!
    ! Returns the distance between two points r1, r2 in the minimum image     !
    ! convention.                                                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   r1, r2 (input): the coordinates of the points in question.            !
    ! Returns:                                                                !
    !   |r1 - r2| in the minimum image convention.                            !
    ! Caveats, preconditions:                                                 !
    !   Both r1 and r2 are assumed to lie within the simulation cell.         !
    !   While this routine will sometimes work if the above is not satisified,!
    !   this should not be relied on.                                         !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use geometry, only: POINT, operator(+), operator(-), operator(*), &
         operator(.dot.), magnitude
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    type(POINT), intent(in) :: r1, r2
    type(CELL_INFO), intent(in) :: cell

    ! jd: Local variables
    type(POINT) :: rvec0, rvec, rvec_final
    integer :: i1, i2, i3
    real(kind=DP) :: r, rmin
    real(kind=DP) :: r_along_a1, r_along_a2, r_along_a3
    real(kind=DP) :: half_a1, half_a2, half_a3

    !------------------------------------------------------------------------

    ! jd: Calculate the original displacement between the two
    rvec0 = r2-r1

    ! jd: Check 27 candidate displacement vectors to see which one is the
    !     shortest one, including the original one. This is the simplest way
    !     in non-orthorhombic cells. Note that this fails if r1 or r2 lie
    !     sufficiently far away from the cell
    rmin = huge(1.0_DP)
    do i1 = -1, 1
       do i2 = -1, 1
          do i3 = -1, 1
             rvec = rvec0 + &
                  real(i1,kind=DP) * cell%a1 + &
                  real(i2,kind=DP) * cell%a2 + &
                  real(i3,kind=DP) * cell%a3
             r = magnitude(rvec)
             if (r < rmin) then
                rmin=r
                rvec_final = rvec
             end if
          end do
       end do
    end do

    r_along_a1 = rvec_final .dot. cell%a1_unit
    r_along_a2 = rvec_final .dot. cell%a2_unit
    r_along_a3 = rvec_final .dot. cell%a3_unit

    half_a1 = magnitude(cell%a1)*0.5_DP
    half_a2 = magnitude(cell%a2)*0.5_DP
    half_a3 = magnitude(cell%a3)*0.5_DP

    if(r_along_a1 > half_a1 .or. r_along_a1 < -half_a1) then
       call utils_assert(.false.,'minimum_image_distance: distance component &
            &along cell vector a1 longer than half of a1.', r_along_a1, half_a1)
    end if

    if(r_along_a2 > half_a2 .or. r_along_a2 < -half_a2) then
       call utils_assert(.false.,'minimum_image_distance: distance component &
            &along cell vector a2 longer than half of a2.', r_along_a2, half_a2)
    end if

    if(r_along_a3 > half_a3 .or. r_along_a3 < -half_a3) then
       call utils_assert(.false.,'minimum_image_distance: distance component &
            &along cell vector a3 longer than half of a3.', r_along_a3, half_a3)
    end if

    minimum_image_distance = rmin

  end function minimum_image_distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function minimum_image_distance_ortho(r1,r2,cell)
    !=========================================================================!
    ! Returns the distance between two points r1, r2 in the minimum image     !
    ! convention. Uses the fast algorithm applicable only for orthorhombic    !
    ! cells.                                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   r1, r2 (input): the coordinates of the points in question.            !
    ! Returns:                                                                !
    !   |r1 - r2| in the minimum image convention.                            !
    ! Caveats, preconditions:                                                 !
    !   Only guaranteed to work correctly if cell is orthorhombic.            !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use geometry, only: POINT, operator(+), operator(-), operator(*), magnitude
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    type(POINT), intent(in) :: r1, r2
    type(CELL_INFO), intent(in) :: cell

    ! jd: Local variables
    type(POINT) :: rvec0, rvec
    integer :: i1, i2, i3
    real(kind=DP) :: r, rmin
    real(kind=DP) :: cell_lx_half
    real(kind=DP) :: cell_ly_half
    real(kind=DP) :: cell_lz_half

    !------------------------------------------------------------------------

    cell_lx_half = 0.5_DP * cell%a1%X
    cell_ly_half = 0.5_DP * cell%a2%Y
    cell_lz_half = 0.5_DP * cell%a3%Z

    if(cell_lx_half == 0.0_DP .or. cell_ly_half == 0.0_DP .or. &
         cell_lz_half == 0.0_DP) then
       call utils_abort("minimum_image_distance_ortho: Your simulation cell &
            &vectors are unorthodox for an orthorhombic set-up. Please have &
            &a non-zero X component in the first vector, a non-zero Y compon&
            &ent in the second vector, and a non-zero Z component in the thi&
            &rd vector.")
    end if

    ! jd: Calculate the original displacement between the two
    rvec0 = r2-r1

    ! jd: Adjust all components to lie in (-L/2,L/2].
    do while(rvec0%X > cell_lx_half)
       rvec0%X = rvec0%X - cell%a1%X
    end do
    do while(rvec0%Y > cell_ly_half)
       rvec0%Y = rvec0%Y - cell%a2%Y
    end do
    do while(rvec0%Z > cell_lz_half)
       rvec0%Z = rvec0%Z - cell%a3%Z
    end do

    do while(rvec0%X <= -cell_lx_half)
       rvec0%X = rvec0%X + cell%a1%X
    end do
    do while(rvec0%Y <= -cell_ly_half)
       rvec0%Y = rvec0%Y + cell%a2%Y
    end do
    do while(rvec0%Z <= -cell_lz_half)
       rvec0%Z = rvec0%Z + cell%a3%Z
    end do

    minimum_image_distance_ortho = magnitude(rvec0)

  end function minimum_image_distance_ortho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine minimum_image_number_to_indices(which_image, &
       a1_neighbour, a2_neighbour, a3_neighbour)
    !=========================================================================!
    ! Converts a minimum image number (-13..13) to three indices (-1..1)      !
    ! corresponding to one-cell displacements along A1, A2, A3.               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   which_image (in): The minimum image number, e.g. one returned by      !
    !                     minimum_image_displacement.                         !
    !   a{1,2,3}_neighbour (out): Ret'd displacements (-1..1) along A{1,2,3}. !
    !-------------------------------------------------------------------------!
    ! Details:                                                                !
    ! which_image   a1  a2  a3                                                !
    !         -13   -1  -1  -1                                                !
    !         -12   -1  -1   0                                                !
    !         -11   -1  -1   1                                                !
    !         -10   -1   0  -1                                                !
    !          -9   -1   0   0                                                !
    !          -8   -1   0   1                                                !
    !          -7   -1   1  -1                                                !
    !          -6   -1   1   0                                                !
    !          -5   -1   1   1                                                !
    !          -4    0  -1  -1                                                !
    !          -3    0  -1   0                                                !
    !          -2    0  -1   1                                                !
    !          -1    0   0  -1                                                !
    !           0    0   0   0                                                !
    !           1    0   0   1                                                !
    !           2    0   1  -1                                                !
    !           3    0   1   0                                                !
    !           4    0   1   1                                                !
    !           5    1  -1  -1                                                !
    !           6    1  -1   0                                                !
    !           7    1  -1   1                                                !
    !           8    1   0  -1                                                !
    !           9    1   0   0                                                !
    !          10    1   0   1                                                !
    !          11    1   1  -1                                                !
    !          12    1   1   0                                                !
    !          13    1   1   1                                                !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 09/2022                                      !
    !=========================================================================!

    use rundat, only: pub_debug
    use utils, only: utils_abort, utils_int_to_str

    implicit none

    ! jd: Arguments
    integer, intent(in)  :: which_image
    integer, intent(out) :: a1_neighbour, a2_neighbour, a3_neighbour

    ! jd: Local variables
    integer :: i
    character(len=*), parameter :: myself = 'minimum_image_number_to_indices'

    ! -------------------------------------------------------------------------

    if(pub_debug) then
       if(.not. (which_image >= -13 .and. which_image <= 13)) then
          call utils_abort(myself//': Illegal periodic image index: '//&
               trim(utils_int_to_str(which_image))//'.')
       end if
    end if

    ! jd: We shift everything to non-negative first to avoid any mod weirdness
    !     with negative numbers.
    i = which_image + 13
    a1_neighbour = i/9-1
    a3_neighbour = mod(i,3)-1
    a2_neighbour = (which_image-9*a1_neighbour - a3_neighbour)/3

    if(pub_debug) then
       if(a1_neighbour /= -1 .and. a1_neighbour /= 0 .and. &
            a1_neighbour /= 1) then
          call utils_abort(myself//': Illegal a1_neighbour: '//&
               trim(utils_int_to_str(a1_neighbour))//', image index: '//&
               trim(utils_int_to_str(which_image))//'.')
       end if
       if(a2_neighbour /= -1 .and. a2_neighbour /= 0 .and. &
            a2_neighbour /= 1) then
          call utils_abort(myself//': Illegal a2_neighbour: '//&
               trim(utils_int_to_str(a2_neighbour))//', image index: '//&
               trim(utils_int_to_str(which_image))//'.')
       end if
       if(a3_neighbour /= -1 .and. a3_neighbour /= 0 .and. &
            a3_neighbour /= 1) then
          call utils_abort(myself//': Illegal a3_neighbour: '//&
               trim(utils_int_to_str(a3_neighbour))//', image index: '//&
               trim(utils_int_to_str(which_image))//'.')
       end if
    end if

  end subroutine minimum_image_number_to_indices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function minimum_image_number_to_displacement(which_image, cell)
    !=========================================================================!
    ! Converts a minimum image number (-13..13) to a displacement vector to   !
    ! the cell containing the image.                                          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   which_image (in): The minimum image number, e.g. one returned by      !
    !                     minimum_image_displacement.                         !
    ! Return value:                                                           !
    !   The displacement vector from the home cell to the image cell.         !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 09/2022                                      !
    !=========================================================================!

    use geometry, only: POINT, operator(+), operator(*)

    implicit none

    ! jd: Arguments
    integer, intent(in)         :: which_image
    type(CELL_INFO), intent(in) :: cell
    type(POINT) :: minimum_image_number_to_displacement

    ! jd: Local variables
    integer :: a1_neighbour, a2_neighbour, a3_neighbour

    ! -------------------------------------------------------------------------

    call minimum_image_number_to_indices(which_image, &
         a1_neighbour, a2_neighbour, a3_neighbour)

    minimum_image_number_to_displacement = &
         real(a1_neighbour,kind=DP)*cell%a1 + &
         real(a2_neighbour,kind=DP)*cell%a2 + &
         real(a3_neighbour,kind=DP)*cell%a3

  end function minimum_image_number_to_displacement

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine minimum_image_displacement(rad, disp, & ! out
       curpoint, centre_in_cell, cell, &             ! in
       which_image)                                  ! out, opt
    !==========================================================================!
    ! This subroutine finds the displacement vector using the minimum image    !
    ! convention.                                                              !
    ! @docme
    !--------------------------------------------------------------------------!
    ! Split from internal_swpot_point by Quintin Hill on 29/04/2010.           !
    ! Fixed by Jacek Dziedzic for the case of non-cubic boxes, 28.06.2012.     !
    ! Removed 'coulomb_cutoff' cruft, Jacek Dziedzic, 08.04.2014.              !
    ! Moved to simulation_cell_mod, Jacek Dziedzic 17.02.2015.                 !
    !==========================================================================!

    use constants, only: DP
    use geometry, only: POINT, operator(+), operator(*), operator(-), &
         magnitude, operator(.dot.)
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)  :: rad  ! Magnitude of displacement
    type(POINT), intent(out)    :: disp   ! Displacement vector of chosen point
    type(POINT), intent(in)     :: curpoint ! Curent point wrt cell
    type(POINT), intent(in)     :: centre_in_cell ! Centre of reference wrt cell
    type(CELL_INFO), intent(in) :: cell
    integer, intent(out), optional :: which_image

    ! jd: Local variables
    integer       :: a1_neighbour, a2_neighbour, a3_neighbour
    type(POINT)   :: periodic_centre
    type(POINT)   :: trial_disp
    real(kind=DP) :: trial_rad
    real(kind=DP) :: r_along_a1, r_along_a2, r_along_a3
    real(kind=DP) :: half_a1, half_a2, half_a3

    ! -------------------------------------------------------------------------

    rad = huge(1.0_DP)
    ! qoh: loop over periodic images to find appropriate one
    do a1_neighbour = -1,1
       do a2_neighbour = -1,1
          do a3_neighbour = -1,1

             periodic_centre = centre_in_cell &
                  + real(a1_neighbour,kind=DP)*cell%a1 &
                  + real(a2_neighbour,kind=DP)*cell%a2 &
                  + real(a3_neighbour,kind=DP)*cell%a3

             trial_disp = curpoint - periodic_centre
             trial_rad = magnitude(trial_disp)

             if(trial_rad < rad) then
                rad = trial_rad
                disp = trial_disp
                if(present(which_image)) then
                   which_image = &
                        a1_neighbour * 9 + a2_neighbour * 3 + a3_neighbour
                end if
             end if

          end do
       end do
    end do

    r_along_a1 = disp .dot. cell%a1_unit
    r_along_a2 = disp .dot. cell%a2_unit
    r_along_a3 = disp .dot. cell%a3_unit
    half_a1 = magnitude(cell%a1)*0.5_DP
    half_a2 = magnitude(cell%a2)*0.5_DP
    half_a3 = magnitude(cell%a3)*0.5_DP

    if(r_along_a1 > half_a1 .or. r_along_a1 < -half_a1) then
       call utils_assert(.false.,'minimum_image_displacement: distance component &
            &along cell vector a1 longer than half of a1.', r_along_a1, half_a1)
    end if

    if(r_along_a2 > half_a2 .or. r_along_a2 < -half_a2) then
       call utils_assert(.false.,'minimum_image_displacement: distance component &
            &along cell vector a2 longer than half of a2.', r_along_a2, half_a2)
    end if

    if(r_along_a3 > half_a3 .or. r_along_a3 < -half_a3) then
       call utils_assert(.false.,'minimum_image_displacement: distance component &
            &along cell vector a3 longer than half of a3.', r_along_a3, half_a3)
    end if

  end subroutine minimum_image_displacement

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulation_cell_num_periodic_images(lattice, cutoff, periodicity)
    !==========================================================================!
    ! This subroutine calculates the maximum number of periodic images of the  !
    ! cell with a given set of lattice vectors within a cutoff.                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   lattice(in) :  lattice vectors of the parent simulation cell.          !
    !   cutoff(in)  :  cutoff for calculating the number of periodic images.   !
    !   periodicity(out): number of periodic images along each lattice vector. !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Dec 2021                                  !
    !==========================================================================!

    use geometry, only: POINT, OPERATOR(.CROSS.), OPERATOR(.DOT.), unit_vector

    implicit none

    ! ab: Arguments
    type(POINT), intent(in) :: lattice(3) ! ab: array of lattice vectors
    real(kind=DP), intent(in) :: cutoff ! ab: cutoff for counting periodicity
    integer, intent(out) :: periodicity(3) ! ab: number of periodic images

    ! ab: local variables
    type(POINT) :: normal1, normal2, normal3 ! ab: normals to lattice vectors
    real(kind=DP) :: cos1, cos2, cos3 ! ab: angle between lattice vectors and
                                      ! normals
    ! ab: Get normals to lattice vectors
    normal1 = lattice(2) .CROSS. lattice(3)
    normal2 = lattice(3) .CROSS. lattice(1)
    normal3 = lattice(1) .CROSS. lattice(2)

    ! ab: Normalize to unit length
    normal1 = unit_vector(normal1)
    normal2 = unit_vector(normal2)
    normal3 = unit_vector(normal3)

    ! ab: angle between lattice vectors and normals
    cos1 = normal1 .DOT. lattice(1)
    cos2 = normal2 .DOT. lattice(2)
    cos3 = normal3 .DOT. lattice(3)

    ! ab: Determine the number of periodic images
    periodicity(1) = ceiling(abs(cutoff/cos1))
    periodicity(2) = ceiling(abs(cutoff/cos2))
    periodicity(3) = ceiling(abs(cutoff/cos3))

  end subroutine simulation_cell_num_periodic_images

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulation_cell_translations(lattice, cutoff, trans)
    !==========================================================================!
    ! This subroutine generates a distance-sorted list of valid translations of!
    ! cell within a given cutoff.                                              !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   lattice(in)     :  lattice vectors of the parent simulation cell.      !
    !   cutoff(in)      :  cutoff for calculating the valid translations.      !
    !   trans(inout)    :  distance-sorted list of valid translations.         !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Dec 2021                                  !
    !==========================================================================!

    use geometry, only: POINT, OPERATOR(.DOT.), OPERATOR(*), OPERATOR(+)
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort

    implicit none

    ! ab: Arguments
    type(POINT), intent(in) :: lattice(3)
    real(kind=DP) :: cutoff
    integer, intent(inout), allocatable :: trans(:,:)

    ! ab: local variables
    integer :: num_cells, ntr, xtr, ytr, ztr, ierr, itr
    character(len=*), parameter :: myself='simulation_cell_translations'
    integer :: periodicity(3) ! ab: max. number of periodic images
    real(kind=DP) :: sq_cutoff, sq_mag
    real(kind=DP), allocatable :: sq_dist(:)
    integer, allocatable :: id_dist(:)
    integer, allocatable :: translat(:,:)
    type(POINT) :: vec

    ! ab: determine the num of periodic images of lattice
    call simulation_cell_num_periodic_images(lattice, cutoff, periodicity)

    ! ab: max. number of cells
    num_cells = product(2*periodicity+1)

    ! ab: allocate space for possible translations
    allocate(translat(3,num_cells), stat=ierr)
    call utils_alloc_check(myself, 'translat', ierr)

    ! ab: allocate distance array
    allocate(sq_dist(num_cells), stat=ierr)
    call utils_alloc_check(myself, 'sq_dist', ierr)

    ! ab: square cutoff for valid translations
    sq_cutoff = cutoff**2

    ! ab: loop over all translations and store valid translations
    ntr = 0
    do xtr = -periodicity(1), periodicity(1)
       do ytr = -periodicity(2), periodicity(2)
          do ztr = -periodicity(3), periodicity(3)
             vec= real(xtr,kind=DP)*lattice(1) + &
                  real(ytr,kind=DP)*lattice(2) + &
                  real(ztr,kind=DP)*lattice(3)
             sq_mag = vec .DOT. vec
             if (sq_mag <= sq_cutoff) then
                ntr = ntr + 1
                sq_dist(ntr) = sq_mag
                translat(:,ntr) = [xtr, ytr, ztr]
             end if
          end do
       end do
    end do

    ! ab: sort the distances
    allocate(id_dist(ntr), stat=ierr)
    call utils_alloc_check(myself, 'id_dist', ierr)
    call utils_heapsort(ntr,sq_dist(1:ntr),id_dist)

    ! ab: sort the list of translations
    allocate(trans(3,ntr), stat=ierr)
    call utils_alloc_check(myself, 'trans', ierr)
    do itr = 1, ntr
       trans(:,itr) = translat(:,id_dist(itr))
    end do

    ! ab: deallocate arrays
    deallocate(id_dist,stat=ierr)
    call utils_dealloc_check(myself,'id_dist', ierr)
    deallocate(sq_dist,stat=ierr)
    call utils_dealloc_check(myself,'sq_dist', ierr)
    deallocate(translat,stat=ierr)
    call utils_dealloc_check(myself,'translat', ierr)

  end subroutine simulation_cell_translations

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulation_cell_lattice_points(lattice, trans, points)
    !==========================================================================!
    ! This subroutine generates lattice points given all translations of the   !
    ! cell.                                                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   lattice(in) :  lattice vectors of the parent simulation cell.          !
    !   trans(in)   :  list of valid translations.                             !
    !   points(out) :  list of lattice points.                                 !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Dec 2021                                  !
    !==========================================================================!

    use geometry, only: POINT, OPERATOR(*), OPERATOR(+)

    implicit none

    ! ab: Arguments
    type(POINT), intent(in)    :: lattice(3)
    integer, intent(in)        :: trans(:,:)
    type(POINT), intent(inout) :: points(:)

    ! ab: local variables
    integer :: ntr, itr, id

    ntr = size(trans, dim=2)
    do itr = 1, ntr
       points(itr) = &
            real(trans(1,itr),kind=DP)*lattice(1) + &
            real(trans(2,itr),kind=DP)*lattice(2) + &
            real(trans(3,itr),kind=DP)*lattice(3)
    end do

  end subroutine simulation_cell_lattice_points

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulation_cell_lattice_points_init(lattice, cutoff, points)
    !==========================================================================!
    ! This subroutine generates distance-sorted list of lattice points within  !
    ! a cutoff.                                                                !
    ! points(:) is allocated within this subroutine.                           !
    ! Must call simulation_cell_points_exit to deallocate the space allocated  !
    ! for lattice points.                                                      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   lattice(in)   :  lattice vectors of the parent simulation cell.        !
    !   cutoff(in)    :  cutoff for calculating the valid points.              !
    !   points(inout) :  list of lattice points.                               !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Dec 2021                                  !
    !==========================================================================!

    use geometry, only: POINT
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! ab: Arguments
    type(POINT), intent(in) :: lattice(3)
    real(kind=DP), intent(in) :: cutoff
    type(POINT), intent(inout), allocatable :: points(:)

    ! ab: local variables
    integer :: ierr, num_cells, ntr
    integer, allocatable :: trans(:,:)
    character(len=*), parameter :: myself='simulation_cell_lattice_points_init'

    ! ab: Generate distance-sorted list of valid translations
    call simulation_cell_translations(lattice, cutoff, trans)

    ! ab: number of valid translations
    ntr = size(trans, dim=2)

    ! ab: allocate space for lattice points
    allocate(points(ntr), stat=ierr)
    call utils_alloc_check(myself, 'points', ierr)

    ! ab: Generate lattice points
    call simulation_cell_lattice_points(lattice, trans, points)

    ! ab: deallocate translations
    deallocate(trans, stat=ierr)
    call utils_dealloc_check(myself,'trans', ierr)


  end subroutine simulation_cell_lattice_points_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulation_cell_lattice_points_exit(points)
    !==========================================================================!
    ! This subroutine deallocates space previously allocated for lattice points!
    ! in simulation_cell_lattice_points_init().                                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   points(inout) :  list of lattice points.                               !
    !--------------------------------------------------------------------------!
    ! Written by Arihant Bhandari in Dec 2021                                  !
    !==========================================================================!

    use geometry, only: POINT
    use utils, only: utils_dealloc_check

    implicit none

    ! ab: Arguments
    type(POINT), intent(inout), allocatable :: points(:)

    ! ab: local variables
    integer :: ierr
    character(len=*), parameter :: myself='simulation_cell_lattice_points_exit'

    ! ab: deallocate space for lattice points
    deallocate(points, stat=ierr)
    call utils_dealloc_check(myself, 'points', ierr)

  end subroutine simulation_cell_lattice_points_exit

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module simulation_cell
