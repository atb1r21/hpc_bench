! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Nicholas D.M. Hine
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module model_type

  use cell_grid, only: GRID_INFO
  use constants, only: DP
  use dftb_type, only: DFTB_STATE
  use dftb_parameters, only: DFTB_PARAMS
  use fft_box, only: FFTBOX_INFO
  use ion, only: ELEMENT, SPECIE, RADIAL_NGWF_TYPE, RADIAL_DENSITY_TYPE
  ! agreco: store information about unique k-points as well
  use k_points, only: KPOINT
  use parallel_strategy, only: PARAL_INFO
  use paw, only: PAW_SPECIES
  use pseudopotentials, only: PSEUDO_SPECIES
  use simulation_cell, only: CELL_INFO, castep_model

  implicit none

  private

  public :: model_allocate
  public :: model_deallocate
  public :: model_write_castep_cell_file

  ! rc2013: new type for  holding info related to multiple regions
!=============================================================
!================= REGION INFORMATION ========================
  type REGION

     ! rc2013: region counter within mdl (regions(isub)%reg_count = isub)
     integer                    :: reg_count
     type(ELEMENT), allocatable :: elements(:)
     ! TODO: fix classical_elements for embedding
     type(ELEMENT), allocatable :: classical_elements(:)
     type(SPECIE), allocatable :: species(:)
     type(RADIAL_NGWF_TYPE), allocatable :: radial_ngwfs(:,:)
     type(RADIAL_DENSITY_TYPE), allocatable :: radial_densities(:)

     ! rc2013: allow multiple pars in the model
     type(PARAL_INFO) :: par
     type(PSEUDO_SPECIES), pointer :: pseudo_sp(:)
     type(PAW_SPECIES), pointer :: paw_sp(:)

     ! rc2013: Ewald energy for subsystem (unused)
     real(kind=DP) :: ewald_energy, dispersion_energy
     ! rc2013: are there any nonlocal projectors in this region?
     logical :: any_nl_proj ! modified in pseudopotentials_mod
     ! jcap: core density on the fine grid for this region - only used
     ! extensively for the active region in an embedding calculation
     real(kind=DP), dimension(:,:,:), allocatable :: core_density_fine


  end type REGION
!================= REGION INFORMATION ========================
!=============================================================

  ! type for containing general model information that is not specific to
  ! the current electronic state
  type MODEL

     type(REGION), allocatable :: regions(:)

     ! ndmh 20_12_2014: species type is still not in use
     ! rc2013: full system parameters i.e. amalgamations of REGION

     type(ELEMENT), allocatable :: elements(:)
     type(ELEMENT), allocatable :: classical_elements(:)
     type(SPECIE), allocatable :: species(:)
     type(RADIAL_NGWF_TYPE), allocatable :: radial_ngwfs(:,:)
     type(RADIAL_DENSITY_TYPE), allocatable :: radial_densities(:) ! jd: <- Does not seem to be set

     type(PARAL_INFO), pointer :: par
     type(PSEUDO_SPECIES), pointer :: pseudo_sp(:)
     type(PAW_SPECIES), pointer :: paw_sp(:)
     real(kind=DP) :: ewald_energy, dispersion_energy
     ! rc2013: I think it's best to make localpseudo_fine system wide
     real(kind=DP), dimension(:,:,:), allocatable :: localpseudo_fine
     ! jcap: system-wide non-local projectors check
     ! rc2013: are there any nonlocal projectors?
     logical :: any_nl_proj
     real(kind=DP), dimension(:,:,:), allocatable :: core_density_fine
     integer :: num_pspecies, num_species

!================= UNIVERSAL INFORMATION ========================
     ! rc2013: the following are universal for all subsystems

     integer :: nat, nat_classical, nsp, nsub, nkpoints
     ! agreco: store information about unique k-points
     type(KPOINT), allocatable :: unique_kpoints(:)

     ! rc2013: all grid+cell info should be consistent across subsystems
     type(GRID_INFO) :: std_grid
     type(GRID_INFO) :: dbl_grid
     type(GRID_INFO) :: fine_grid
     type(GRID_INFO) :: padded_grid

     type(CELL_INFO) :: cell
     type(CELL_INFO) :: padded_cell

     ! rc2013: not sure about fftbox though...
     ! rc2013: it's set for the simulation cell so it should be fine
     type(FFTBOX_INFO) :: fftbox
     type(FFTBOX_INFO) :: uni_tightbox
     type(FFTBOX_INFO) :: aug_box

     type(castep_model) :: cmdl

     real(kind=DP) :: locps_gzero_term

     ! jd: Parameters of the DFTB method
     type(DFTB_PARAMS) :: dftb_par
     type(DFTB_STATE)  :: dftb_var

!================= UNIVERSAL INFORMATION ========================

  end type MODEL


  public :: MODEL, REGION

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine model_allocate(mdl,nat,class_nat,nsp)

    implicit none

    ! Arguments
    type(MODEL), intent(inout) :: mdl
    integer, intent(in) :: nat
    integer, intent(in) :: class_nat
    integer, intent(in) :: nsp

    mdl%nat_classical = class_nat
    mdl%nat = nat
    mdl%nsp = nsp

  end subroutine model_allocate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine model_deallocate(mdl)

    implicit none

    ! Arguments
    type(MODEL), intent(inout) :: mdl

  end subroutine model_deallocate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine model_write_castep_cell_file(filename,mdl)

    !==========================================================================!
    ! This routine writes a CASTEP cell file for use with OptaDOS.             !
    !--------------------------------------------------------------------------!
    ! Written by Jolyon Aarons. 23/07/2014                                     !
    ! Moved to model_type_mod so it can be used elsewhere by Nicholas Hine     !
    ! in February 2016.                                                        !
    ! Edited for multiple subsystems by Robert Charlton, 18/01/2017.           !
    !==========================================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Arguments
    character(len=*), intent(in) :: filename
    type(MODEL), intent(in) :: mdl

    ! Local Variables
    integer :: iat
    integer :: cellfile
    logical, allocatable :: species_written(:)
    integer :: num_species
    integer :: ierr


    if(pub_on_root) then

       write(stdout,*)''
       write(stdout,'(3a)', advance='no')&
            'Writing Castep output cell file to "',trim(adjustl(filename)),'" ...'

       cellfile = utils_unit()
       open(unit=cellfile, file=trim(adjustl(filename)), form='formatted', &
            action='write', iostat=ierr, status='replace')
       call utils_open_unit_check('model_write_castep_cell_file','cellfile',ierr)

       write(cellfile,*) "%block symmetry_ops"
       write(cellfile,*) "%endblock symmetry_ops"
       write(cellfile,*)

       ! ja531 -> Write cell vecs
       write(cellfile,*) "%block lattice_cart"
       write(cellfile,*) "Bohr"
       write(cellfile,*) mdl%cell%a1%x,mdl%cell%a1%y,mdl%cell%a1%z
       write(cellfile,*) mdl%cell%a2%x,mdl%cell%a2%y,mdl%cell%a2%z
       write(cellfile,*) mdl%cell%a3%x,mdl%cell%a3%y,mdl%cell%a3%z
       write(cellfile,*) "%endblock lattice_cart"

       write(cellfile,*)

       write(cellfile,*) "%block species_pot"
       num_species=maxval((/(mdl%regions(1)%elements(iat)%species_number,iat=1,mdl%regions(1)%par%nat)/))
       allocate(species_written(num_species),stat=ierr)
       call utils_alloc_check('model_write_castep_cell_file','species_written',ierr)
       species_written=.false.
       ! ja531 -> Write species pots
       do iat=1,mdl%regions(1)%par%nat
          if(.not.species_written(mdl%regions(1)%elements(iat)%species_number)) then
             write(cellfile,*) trim(mdl%regions(1)%elements(iat)%species_id), " ", &
                  trim(mdl%regions(1)%species(mdl%regions(1)%elements(iat)%species_number)%&
                     pseudo_name)
             species_written(mdl%regions(1)%elements(iat)%species_number)=.true.
          end if
       end do
       deallocate(species_written,stat=ierr)
       call utils_dealloc_check('model_write_castep_cell_file','species_written',ierr)
       write(cellfile,*) "%endblock species_pot"

       write(cellfile,*)

       write(cellfile,*) "%block positions_abs"
       write(cellfile,*) "Bohr"
       ! ja531 -> Write atomic pos
       do iat=1,mdl%regions(1)%par%nat
          write(cellfile,*) trim(mdl%regions(1)%elements(iat)%species_id), &
               mdl%regions(1)%elements(iat)%centre%x, &
               mdl%regions(1)%elements(iat)%centre%y, &
               mdl%regions(1)%elements(iat)%centre%z
       end do
       write(cellfile,*) "%endblock positions_abs"

       close(unit=cellfile,iostat=ierr)
       call utils_close_unit_check('model_write_castep_cell_file','cellfile',ierr)
       write(stdout,'(a)') ' done'
    end if

    call comms_barrier

  end subroutine model_write_castep_cell_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module model_type

