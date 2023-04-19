! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                       F R A G M E N T   D A T A                             !
!                            O B J E C T S                                    !
!=============================================================================!
!                                                                             !
! This module handles the pub_frag_data data type structure, the buffer       !
! storage of the EDA delta densities, and other fragment and energy           !
! decomposition analysis related data.                                        !
!                                                                             !
!-----------------------------------------------------------------------------!
! Originally written by Max Phipps August 2014                                !
!-----------------------------------------------------------------------------!

module fragment_data

  use model_type, only: MODEL
  use ngwf_representation, only: NGWF_REP, NGWF_HAM
  use function_basis, only: FUNC_BASIS
  use constants, only: dp
  use sparse_embed, only: SPAM3_EMBED
  use dense, only: DEM

  implicit none

  private


  type FRAG_DATA

     logical :: any_nl_proj

     integer :: num_spins

     integer :: nat, nat_classical

     type(MODEL) :: mdl

     type(NGWF_REP) :: rep

!     type(NGWF_HAM) :: ham  ! currently unused

     type(FUNC_BASIS), allocatable      :: ngwf_basis(:)

     integer, allocatable, dimension(:) :: basis_ref ! The fragment
                        !   ordering within the supermolecule order.
                        ! i.e. 1,4,3,2,6,5 for 2x 3atom fragments
                        !   would be indexes 1,4,3 of supermolecule
                        !   for fragment 1 of monomer ordering
                        !   and 2,6,5 of supermolecule for fragment 2
                        !   in monomer ordering.
                        ! This is used to
                        ! reference the fragment
                        ! ngwf_basis spheres in the supermolecule ngwf_basis
                        ! spheres array.

     type(DEM), allocatable, dimension(:) :: denskern_dens      ! DEM ScaLAPACK distributed
     type(DEM), allocatable, dimension(:) :: denskern_dens_loc  ! DEM on root  (ScaLAPACK) or
                                                                ! on all procs (Non-ScaLAPACK)

  end type FRAG_DATA

  !>
  !!  @brief Maps where array indices are NGWF indices.
  !!
  !!  @internal
  !!    @author hc5n17
  !!    @since 2019
  !!  @endinternal
  !!
  type, public :: NGWF_INDEX_MAP_TYPE
     !>
     !!  @brief Atom identifiers sorted by cumulative fragment NGWF index.
     !!
     !!  It is assumed that it is the same between calculations.
     !!
     integer, allocatable :: atom_ids_by_cum_frag(:)

     !>
     !!  @brief Atom identifiers sorted by supermolecule NGWF index.
     !!
     !!  It should not be assumed that it is the same between calculations.
     !!
     integer, allocatable :: atom_ids_by_super(:)

     !>
     !!  @brief Supermolecule NGWF indices sorted by cumulative fragment NGWF
     !!  index.
     !!
     !!  It should not be assumed that it is the same between calculations.
     !!
     integer, allocatable :: supers_by_cum_frag(:)
  end type NGWF_INDEX_MAP_TYPE


  ! public fragment and whole-system data
  type(FRAG_DATA), public, target, allocatable, dimension(:) :: pub_frag_data

  ! Array defining which supermolecule NGWFs belong on a particular
  ! fragment.
  ! (length=number of NGWFs on supermolecule, indexed as [NGWF_super, fragID])
  logical, public, dimension(:,:), allocatable :: pub_super2frag_ngwf_log

  ! The fragment ID of each NGWF on this proc (i.e. which fragment each NGWF
  ! belongs to).
  integer, public, dimension(:), allocatable :: pub_super2fragid_on_proc

  ! EDA delta density calculations:
  real(kind=DP), public, allocatable :: eda_frzIdem_density_fine(:,:,:,:)! frozen state density
  real(kind=DP), public, allocatable :: eda_pol_density_fine(:,:,:,:,:)! polarized state density (one density per fragment)
  real(kind=DP), public, allocatable :: eda_poliso_density_fine(:,:,:,:)! net polarized density for the individual fragments
  real(kind=DP), public, allocatable :: eda_full_density_fine(:,:,:,:)! fully optimized density

  ! EDA system electronic energies:
  real(kind=DP),public   :: complex_energy_frz   ! Frozen monomer system energy
  real(kind=DP),public   :: complex_energy_frzIdem ! MO occupancy-fixed frozen system energy
  real(kind=DP),public   :: complex_energy_pol_simul ! Simultaneously polarized system energies
  real(kind=DP),public   :: complex_energy_full  ! Full system energy
  real(kind=DP),allocatable, public   :: complex_energy_pol(:)    ! Fragment-wise polarized system energies
  real(kind=DP),allocatable, public   :: complex_energy_ct(:,:)   ! Fragment pair-wise charge transfer energies

  public :: fragment_data_destroy
  public :: fragment_data_add_fragment

  !>
  !!  @brief Gradient of exchange-correlation energy per unit volume with
  !!  respect to the kinetic-energy density.
  !!
  !!  It is initialized in @link eda_driver_supermol.eda_driver_supermol_run
  !!  @endlink and destroyed in @link fragment_data_destroy @endlink. It is
  !!  public, because it is needed in other modules, and procedures relying on
  !!  this value often have an argument intent of `inout` for this value.
  !!
  !!  According to @link hamiltonian.hamiltonian_lhxc_calculate @endlink, it is
  !!  considered not used if its size is `1` in each dimension. Otherwise, its
  !!  shape must conform to what a comment hidden inside @link
  !!  hamiltonian.hamiltonian_lhxc_calculate @endlink states.
  !!
  !!  @internal
  !!    @author hc5n17
  !!    @since 2019
  !!  @endinternal
  !!
  real(kind = DP), allocatable, public :: eda_dfdtau_fine(:, :, :, :)

  !>
  !!  @brief NGWF index map for second and third stages of EDA.
  !!
  !!  It is set by @link fragment_data_set_ngwf_index_map @endlink.
  !!
  !!  @internal
  !!    @author hc5n17
  !!    @since 2019
  !!  @endinternal
  !!
  type(NGWF_INDEX_MAP_TYPE), private :: ngwf_index_map

  public :: fragment_data_get_ngwf_index_map
  public :: fragment_data_set_ngwf_index_map
  public :: fragment_data_clear_ngwf_index_map
  public :: fragment_data_alloc_paw_core_densities
  public :: fragment_data_dealloc_paw_core_densities

contains

  subroutine fragment_data_destroy()

    !================================================================!
    ! This subroutine handles the deallocations of the fragment      !
    ! data container.                                                !
    !----------------------------------------------------------------!
    ! Written by Max Phipps, June 2015.                              !
    ! Modified for embedding by Joseph Prentice, September 2018      !
    !================================================================!

    use datatypes, only: data_functions_dealloc
    use utils, only: utils_dealloc_check
    use rundat, only: pub_frag_iatm, pub_eda_read_super
    use function_basis, only: function_basis_deallocate
    use parallel_strategy, only: parallel_strategy_exit

    implicit none

    integer :: it, ierr
    character(255) :: errmsg = ""

    ! iterate fragments
    do it=1,pub_frag_iatm(0)

       ! clean up parallel strategy
       if (associated(pub_frag_data(it)%mdl%par)) then
          call parallel_strategy_exit(pub_frag_data(it)%mdl%par)
       end if

       ! Deallocate n_occ:
       ! Fragments:
       deallocate(pub_frag_data(it)%rep%n_occ, stat=ierr)
       call utils_dealloc_check('energy_and_force_calculate','pfd(it)%rep%n_occ',&
                                  ierr)

       if (.not.(pub_eda_read_super)) then

          ! mjsp: deallocate storage for fragment NGWFs
          call data_functions_dealloc(pub_frag_data(it)%rep%ngwfs_on_grid(1))

          ! mjsp: deallocate storage for fragment density kernels
          deallocate(pub_frag_data(it)%denskern_dens_loc, stat=ierr)
          call utils_dealloc_check('fragment_destroy','pub_frag_data(it)%denskern_dens',ierr)
#ifdef SCALAPACK
          ! mjsp: ScaLAPACK distributed matrices:
          deallocate(pub_frag_data(it)%denskern_dens, stat=ierr)
          call utils_dealloc_check('fragment_destroy','pub_frag_data(it)%denskern_dens',ierr)
#endif

          ! mjsp: deallocate storage for fragment elements
          deallocate(pub_frag_data(it)%mdl%elements, stat=ierr)
          call utils_dealloc_check('fragment_destroy','pub_frag_data(it)%mdl%elements',ierr)

          ! mjsp: deallocate storage for supermolecule species
          deallocate(pub_frag_data(it)%mdl%species, stat=ierr)
          call utils_dealloc_check('fragment_destroy','pub_frag_data(it)%mdl%species',ierr)

          ! jcap: do this for the regions as well
          ! mjsp: deallocate storage for fragment elements
          deallocate(pub_frag_data(it)%mdl%regions(1)%elements, stat=ierr)
          call utils_dealloc_check('fragment_destroy',&
               'pub_frag_data(it)%mdl%regions%elements',ierr)

          ! mjsp: deallocate storage for supermolecule species
          deallocate(pub_frag_data(it)%mdl%regions(1)%species, stat=ierr)
          call utils_dealloc_check('fragment_destroy',&
               'pub_frag_data(it)%mdl%regions%species',ierr)

       end if

    end do

    ! iterate supermolecule + fragments
    do it=0,pub_frag_iatm(0)

       if (.not.(pub_eda_read_super)) then

          ! Deallocate the NGWFs
          call function_basis_deallocate(pub_frag_data(it)%ngwf_basis(1))

          if (allocated(pub_frag_data(it)%basis_ref)) then
             ! mjsp: deallocate basis referencing pivot tables
             deallocate(pub_frag_data(it)%basis_ref, stat = ierr)

             ! mjsp: NOTE two dealloc check calls: supermolecule and fragments
             call utils_dealloc_check( &
                  'fragment_destroy', &
                  'pub_frag_data(0)%basis_ref', &
                  ierr)

             call utils_dealloc_check( &
                  'fragment_destroy', &
                  'pub_frag_data(it)%basis_ref', &
                  ierr)
          end if

!          ! mjsp: hack to handle sp_overlap matrix destruction on/off in
!          ! ngwf_rep_destroy
!          logical_buffer = pub_any_nl_proj
!          pub_any_nl_proj = pub_frag_data(it)%any_nl_proj
!
!!          ! mjsp: deallocate NGWF representation
!!          call ngwf_rep_destroy(pub_frag_data(0)%rep)
!
!!          ! mjsp: deallocate supermolecule temporary Hamiltonian
!!          call ngwf_ham_destroy(pub_frag_data(it)%ham)
!
!          ! mjsp: restore pub_any_nl_proj
!          logical_buffer = pub_any_nl_proj

       end if

    end do

    ! mjsp: deallocate storage for supermolecule elements
    deallocate(pub_frag_data(0)%mdl%elements, stat=ierr)
    call utils_dealloc_check('fragment_destroy','pub_frag_data(it)%mdl%elements',ierr)

    ! mjsp: deallocate storage for supermolecule species
    deallocate(pub_frag_data(0)%mdl%species, stat=ierr)
    call utils_dealloc_check('fragment_destroy','pub_frag_data(0)%mdl%species',ierr)

    ! jcap: do this for the regions as well
    ! mjsp: deallocate storage for supermolecule elements
    deallocate(pub_frag_data(0)%mdl%regions(1)%elements, stat=ierr)
    call utils_dealloc_check('fragment_destroy',&
         'pub_frag_data(it)%mdl%regions%elements',ierr)

    ! mjsp: deallocate storage for supermolecule species
    deallocate(pub_frag_data(0)%mdl%regions(1)%species, stat=ierr)
    call utils_dealloc_check('fragment_destroy',&
         'pub_frag_data(0)%mdl%regions%species',ierr)

!    ! mjsp: deallocate NGWF representation
!    call ngwf_rep_destroy(pub_frag_data(0)%rep)

    ! Deallocate n_occ:
    ! Supermolecule:
    deallocate(pub_frag_data(0)%rep%n_occ, stat=ierr)
    call utils_dealloc_check('energy_and_force_calculate','pfd(0)%rep%n_occ',&
                               ierr)

    ! mjsp: deallocate fragment data container
    deallocate(pub_frag_data, stat=ierr)
    call utils_dealloc_check('fragment_destroy','pub_frag_data',ierr)


    ! mjsp: deallocate basis referencing pivot tables
    deallocate(ngwf_index_map%atom_ids_by_cum_frag, stat = ierr, errmsg = errmsg)
    call utils_dealloc_check('fragment_destroy', trim(errmsg), ierr)
    deallocate(ngwf_index_map%atom_ids_by_super, stat = ierr, errmsg = errmsg)
    call utils_dealloc_check('fragment_destroy', trim(errmsg), ierr)
    deallocate(ngwf_index_map%supers_by_cum_frag, stat = ierr, errmsg = errmsg)
    call utils_dealloc_check('fragment_destroy', trim(errmsg), ierr)
    deallocate(pub_super2fragid_on_proc, stat=ierr)
    call utils_dealloc_check('fragment_destroy','pub_super2fragid_on_proc',ierr)
    deallocate(pub_super2frag_ngwf_log, stat=ierr)
    call utils_dealloc_check('fragment_destroy','pub_super2frag_ngwf_log',ierr)

    deallocate(eda_dfdtau_fine, stat = ierr, errmsg = errmsg)
    call utils_dealloc_check('fragment_destroy', trim(errmsg), ierr)

  end subroutine fragment_data_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fragment_data_add_fragment(mdl, denskern, rep, ngwf_basis, &
       frag_id)

    !======================================================================!
    ! Stores a fragment's optimized and purified density kernel and NGWFs  !
    ! for later use in calculating the supermolecule energies.             !
    !----------------------------------------------------------------------!
    ! Written by Max Phipps January 2016 (from energy_and_force_calculate) !
    ! Modified for embedding by Joseph Prentice, September 2018            !
    !======================================================================!

    use datatypes, only: data_functions_alloc, data_functions_copy
    use dense, only: dense_create, dense_convert
    use function_basis, only: FUNC_BASIS
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP
    use rundat, only: PUB_1K, pub_num_spins
    use sparse_embed, only: SPAM3_EMBED_ARRAY
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(MODEL), intent(in)       :: mdl ! fragment model
    type(SPAM3_EMBED_ARRAY), intent(in) :: denskern ! fragment (optimised) density kernel
    type(FUNC_BASIS), intent(in)  :: ngwf_basis(1) ! fragment NGWFs
    type(NGWF_REP), intent(in)    :: rep ! fragment representation
    integer, intent(in)           :: frag_id  ! identifier for the fragment

    ! Local Variables
    integer :: ingwf, is, it, at_offset, ierr


    ! mjsp: Allocate fragment storage for density kernel
    allocate(pub_frag_data(frag_id)%denskern_dens_loc(pub_num_spins), stat=ierr)
    call utils_alloc_check('fragment_data_add_fragment','pub_frag_data(it)%denskern_dens_loc',ierr)
    do is=1,pub_num_spins

       ! mjsp: local_to_root for ScaLAPACK - redistribution later at
       ! supermolecule stage.
       call dense_create(pub_frag_data(frag_id)%denskern_dens_loc(is),&
            mdl%par%num_ngwfs,mdl%par%num_ngwfs,local_to_root=.true.)

       ! mjsp: Internal storage
       call dense_convert(pub_frag_data(frag_id)%denskern_dens_loc(is), &
            denskern%m(is,PUB_1K))

    end do

    ! jcap: allocate ngwfs_on_grid
    allocate(pub_frag_data(frag_id)%rep%ngwfs_on_grid(1),stat=ierr)
    call utils_alloc_check('fragment_data_add_fragment',&
         'pub_frag_data(it)%rep%ngwfs_on_grid',ierr)
    ! mjsp: Allocate fragment storage for NGWFs on grid
    call data_functions_alloc(pub_frag_data(frag_id)%rep%ngwfs_on_grid(1), &
         ngwf_basis(1)%n_ppds*mdl%cell%n_pts, &
         iscmplx=rep%ngwfs_on_grid(1)%iscmplx)

    ! mjsp: Store the NGWFs of this fragment for later frozen density calculation
    call data_functions_copy(pub_frag_data(frag_id)%rep%ngwfs_on_grid(1),&
         rep%ngwfs_on_grid(1))
    pub_frag_data(frag_id)%ngwf_basis(1) = ngwf_basis(1)


    ! mjsp: Allocate and fill the basis_ref array:
    ! This is a referencing scheme between the fragment NGWFs and the atoms in the input file
    ! and is used for constructing the frozen density state.
    allocate(pub_frag_data(frag_id)%basis_ref( &
         pub_frag_data(frag_id)%ngwf_basis(1)%num),stat=ierr)
    call utils_alloc_check('fragment_data_add_fragment','pub_frag_data(it)%basis_ref',ierr)

    ! mjsp: at_offset is used to offset the storage of the atoms for this fragment in
    ! the basis_ref array by the number of atoms in the previous fragment calculations.
    ! Count the number of atoms in all previous fragment calculations:
    at_offset = 0
    do it=1,frag_id-1
       at_offset = at_offset + pub_frag_data(it)%mdl%nat
    end do

    ! mjsp: Fill the basis referencing array:
    do ingwf=1,ngwf_basis(1)%num
       pub_frag_data(frag_id)%basis_ref(ingwf) = &
            mdl%regions(1)%par%orig_atom( ngwf_basis(1)%atom_of_func(ingwf) ) &
            + at_offset
    end do

  end subroutine fragment_data_add_fragment

  !>
  !!  @brief Gets the NGWF index map.
  !!
  !!  Program stops if @link fragment_data_set_ngwf_index_map @endlink has not
  !!  been called.
  !!
  !!  @return
  !!    NGWF index map.
  !!
  !!  @internal
  !!    @author hc5n17
  !!    @since 2019
  !!  @endinternal
  !!
  function fragment_data_get_ngwf_index_map()
    use utils, only: utils_abort

    implicit none

    character(len=*), parameter :: myself = &
         "fragment_data_get_ngwf_index_map()"

    type(NGWF_INDEX_MAP_TYPE) :: fragment_data_get_ngwf_index_map

    if (.not. allocated(ngwf_index_map%atom_ids_by_cum_frag)) then
       call utils_abort(myself // ": NGWF index map is not initialized.")
    end if

    fragment_data_get_ngwf_index_map = ngwf_index_map
  end function fragment_data_get_ngwf_index_map

  !>
  !!  @brief Sets the NGWF index map from an array of atom identifiers sorted
  !!  by cumulative fragment NGWF index.
  !!
  !!  Previous NGWF index map, if exists, is overwritten.
  !!
  !!  @param [in] atom_ids
  !!    One-based atom identifiers, where array indices are one-based
  !!    cumulative fragment NGWF indices.
  !!
  !!  @param [in] mdl
  !!    Model.
  !!
  !!  @param [in] super_ngwf_basis
  !!    Supermolecule NGWF basis.
  !!
  !!  @internal
  !!    @author hc5n17
  !!    @since 2019
  !!    @remark
  !!      Adapted from several places written by Max Phipps.
  !!  @endinternal
  !!
  subroutine fragment_data_set_ngwf_index_map( &
       atom_ids_by_cum_frag, &
       mdl, &
       super_ngwf_basis)

    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    integer, allocatable, intent(in) :: atom_ids_by_cum_frag(:)
    type(MODEL), intent(in) :: mdl
    type(FUNC_BASIS), intent(in) :: super_ngwf_basis

    character(len=*), parameter :: myself = &
         "fragment_data_set_ngwf_index_map()"

    integer, allocatable :: atom_ids_by_super_buf(:)
    integer :: cum_frag_ngwf_index
    integer :: super_ngwf_index
    integer :: ierr
    character(255) :: errmsg = ""

    if (lbound(atom_ids_by_cum_frag, 1) /= 1) then
       call utils_abort( &
            myself // &
            ": Lower bound of the array of atom identifiers is not 1.")
    end if

    call fragment_data_clear_ngwf_index_map()

    ngwf_index_map%atom_ids_by_cum_frag = atom_ids_by_cum_frag

    allocate( &
         ngwf_index_map%atom_ids_by_super(super_ngwf_basis%num), &
         stat = ierr, &
         errmsg = errmsg)

    call utils_alloc_check( &
         "fragment_data_set_ngwf_index_map", &
         trim(errmsg), &
         ierr)

    ! Populate the map of supermolecule NGWF index to atom identifier.
    do super_ngwf_index = 1, super_ngwf_basis%num
       ngwf_index_map%atom_ids_by_super(super_ngwf_index) = &
            mdl%par%orig_atom(super_ngwf_basis%atom_of_func(super_ngwf_index))
    end do

    !ja: modified to avoid error with gfortran 9.1 compiler
    allocate(atom_ids_by_super_buf(size(ngwf_index_map%atom_ids_by_super)),&
         stat = ierr, &
         errmsg = errmsg)

    call utils_alloc_check( &
         "fragment_data_set_atom_ids_by_super_buf", &
         trim(errmsg), &
         ierr)


    atom_ids_by_super_buf = ngwf_index_map%atom_ids_by_super

    allocate( &
         ngwf_index_map%supers_by_cum_frag(size(atom_ids_by_cum_frag)), &
         stat = ierr, &
         errmsg = errmsg)

    call utils_alloc_check( &
         "fragment_data_set_ngwf_index_map", &
         trim(errmsg), &
         ierr)

    ! Populate the map of cumulative fragment NGWF index to supermolecule NGWF
    ! index.
    do cum_frag_ngwf_index = 1, size(atom_ids_by_cum_frag)
       ! Find the first of the remaining supermolecular NGWF indices that has
       ! an atom identifier equal to that of the current cumulative fragment
       ! NGWF index.
       do super_ngwf_index = 1, super_ngwf_basis%num
          if (atom_ids_by_super_buf(super_ngwf_index) == &
               ngwf_index_map%atom_ids_by_cum_frag(cum_frag_ngwf_index)) then

             ngwf_index_map%supers_by_cum_frag(cum_frag_ngwf_index) = &
                  super_ngwf_index

             ! Nullify the atom identifier at the supermolecule NGWF index that
             ! is found.
             atom_ids_by_super_buf(super_ngwf_index) = 0

             exit
          end if
       end do
    end do

    deallocate( &
         atom_ids_by_super_buf, &
         stat = ierr, &
         errmsg = errmsg)

    call utils_dealloc_check( &
         "fragment_data_set_ngwf_index_map", &
         trim(errmsg), &
         ierr)
  end subroutine fragment_data_set_ngwf_index_map

  !>
  !!  @brief Clears and deallocates the NGWF index map.
  !!
  !!  @internal
  !!    @author hc5n17
  !!    @since 2019
  !!  @endinternal
  !!
  subroutine fragment_data_clear_ngwf_index_map()
    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr
    character(255) :: errmsg = ""

    if (allocated(ngwf_index_map%atom_ids_by_cum_frag)) then
       deallocate( &
            ngwf_index_map%atom_ids_by_cum_frag, &
            stat = ierr, &
            errmsg = errmsg)

       call utils_dealloc_check( &
            "fragment_data_clear_ngwf_index_map", &
            trim(errmsg), &
            ierr)

       deallocate( &
            ngwf_index_map%atom_ids_by_super, &
            stat = ierr, &
            errmsg = errmsg)

       call utils_dealloc_check( &
            "fragment_data_clear_ngwf_index_map", &
            trim(errmsg), &
            ierr)

       deallocate( &
            ngwf_index_map%supers_by_cum_frag, &
            stat = ierr, &
            errmsg = errmsg)

       call utils_dealloc_check( &
            "fragment_data_clear_ngwf_index_map", &
            trim(errmsg), &
            ierr)
    end if
  end subroutine fragment_data_clear_ngwf_index_map

  !>
  !!  @brief Allocates the core-density arrays of all regions for PAW.
  !!
  !!  If PAW is not being used, nothing is done. If the core-density array of a
  !!  region is already allocated, nothing is done for that array.
  !!
  !!  @param [inout] mdl
  !!    Model.
  !!
  !!  @internal
  !!    @author hc5n17
  !!    @since 2019
  !!  @endinternal
  !!
  subroutine fragment_data_alloc_paw_core_densities(mdl)
    use rundat, only: pub_paw
    use utils, only: utils_alloc_check

    implicit none

    type(MODEL), intent(inout) :: mdl

    integer :: region_index
    integer :: ierr
    character(255) :: errmsg = ""

    if (.not. pub_paw) then
       return
    end if

    do region_index = lbound(mdl%regions, 1), ubound(mdl%regions, 1)
       if (.not. allocated(mdl%regions(region_index)%core_density_fine)) then
          allocate( &
               mdl%regions(region_index)%core_density_fine( &
               mdl%fine_grid%ld1, &
               mdl%fine_grid%ld2, &
               mdl%fine_grid%max_slabs12), &
               stat = ierr, &
               errmsg = errmsg)

          call utils_alloc_check( &
               "fragment_data_alloc_paw_core_densities", &
               trim(errmsg), &
               ierr)
       end if
    end do
  end subroutine fragment_data_alloc_paw_core_densities

  !>
  !!  @brief Deallocates the core-density arrays of all regions for PAW.
  !!
  !!  If PAW is not being used, nothing is done. If the core-density array of a
  !!  region is already deallocated, nothing is done for that array.
  !!
  !!  @param [out] mdl
  !!    Model.
  !!
  !!  @internal
  !!    @author hc5n17
  !!    @since 2019
  !!  @endinternal
  !!
  subroutine fragment_data_dealloc_paw_core_densities(mdl)
    use rundat, only: pub_paw
    use utils, only: utils_dealloc_check

    implicit none

    type(MODEL), intent(out) :: mdl

    integer :: region_index
    integer :: ierr
    character(255) :: errmsg = ""

    if (.not. pub_paw) then
       return
    end if

    do region_index = lbound(mdl%regions, 1), ubound(mdl%regions, 1)
       if (allocated(mdl%regions(region_index)%core_density_fine)) then
          deallocate( &
               mdl%regions(region_index)%core_density_fine, &
               stat = ierr, &
               errmsg = errmsg)

          call utils_dealloc_check( &
               "fragment_data_dealloc_paw_core_densities", &
               trim(errmsg), &
               ierr)
       end if
    end do
  end subroutine fragment_data_dealloc_paw_core_densities

end module fragment_data
