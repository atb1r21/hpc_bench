! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                     Parallel strategy module                   !
!----------------------------------------------------------------!
! Written by Peter Haynes, 10/7/03                               !
! Modified and Maintained by Nicholas Hine, Oct 2007-present     !
! Modified to deal with multiple parallel strategies by Robert   !
! Charlton, 2017.                                                !
!================================================================!

module parallel_strategy

  use constants, only : DP, stdout, VERBOSE
  use geometry, only: POINT
  use ion, only : element
  use k_points, only: KPOINT

  implicit none

  private

  ! Public subroutines
  public :: parallel_strategy_distr_atoms
  public :: parallel_strategy_distr_funcs
  public :: parallel_strategy_check_atoms
  public :: parallel_strategy_list_overlaps
  public :: parallel_strategy_exit_overlaps
  public :: parallel_strategy_exit
  ! agreco: routine to distribute k-points across procs
  public :: parallel_strategy_distr_kpoints
  public :: parallel_strategy_list_cross_overlaps

  ! Public variables
  type OVERLAP_LIST

    ! Information about an overlap list created by the module
    integer :: max_overlaps                   ! Max # overlaps
    integer, allocatable :: num_overlaps(:)   ! Num overlaps
    integer, allocatable :: overlap_list(:,:) ! Overlap list

  end type OVERLAP_LIST

  public :: OVERLAP_LIST

  type PARAL_INFO

    ! cks: the total number of atoms in the cell
    integer          :: nat
    ! ndmh: the total number of Hubbard atoms in the cell
    integer          :: nat_hub
    ! gibo: the total number of cDFT charge-acceptor atoms in the cell
    integer          :: nat_hub_charge_acceptor
    ! gibo: the total number of cDFT charge-donor atoms in the cell
    integer          :: nat_hub_charge_donor
    ! gibo: the total number of cDFT spin-acceptor atoms in the cell
    integer          :: nat_hub_spin_acceptor
    ! gibo: the total number of cDFT spin-donor atoms in the cell
    integer          :: nat_hub_spin_donor
    ! gibo: the total number of cDFT-active atoms for cDFT+U simulations
    integer          :: nat_cdft_active
    ! ars: the total number of classical atoms in the cell
    integer          :: nat_classical
    ! max # atoms on any proc
    integer          :: max_atoms_on_proc

    ! cks: the total number of NGWFs in the cell
    integer          :: num_ngwfs
    ! ndmh: the total number of conduction NGWFs in the cell
    integer          :: num_ngwfs_cond
    ! ndmh: the total number of auxiliary NGWFs in the cell
    integer          :: num_ngwfs_aux
    ! cks: the total number of (non local pseudopotential) projectors in the cell
    integer          :: num_projectors
    ! ndmh: the total number of PAW partial waves in the cell
    integer          :: num_pawpws
    ! ndmh: the number of hubbard projectors in the cell
    integer          :: num_hub_proj
    ! lr408: the total number of core wavefunctions in the cell
    integer          :: num_corewfs
    ! qoh: the total number of distinct pseudopotential species in the cell
    integer          :: num_pspecies
    ! qoh: the total number of distinct atomic species in the cell
    integer          :: num_species
    ! ndmh: the total number of hubbard species in the cell
    integer          :: num_hub_species
    ! ndmh: the total number of cDFT-species in the cell
    integer          :: num_cdft_species
    ! qoh: the total number of SW basis functions in the cell
    integer          :: num_sw
    ! ja531: the max principal quantum number used in pDOS SW expansion.
    integer          :: max_pdos_n

    ! ndmh: public arrays giving information on the distribution of atoms
    integer, allocatable :: num_atoms_on_proc(:)  ! Num atoms on proc
    integer, allocatable :: first_atom_on_proc(:) ! First atom on proc
    integer, allocatable :: proc_of_atom(:)       ! Proc of each atom

    ! ddor: public arrays giving information on the distribution of Hubbard atoms
    integer, allocatable :: num_hub_atoms_on_proc(:)  ! Num Hubbard atoms on proc
    integer, allocatable :: proc_of_hub_atom(:)       ! Proc of each Hubbard atom
    ! ddor: List of Hubbard atoms on all procs
    integer, allocatable, dimension(:,:) :: hub_atoms_on_proc
    ! ddor: An array hub_hub_nat long, ordered as the input file,
    !       which contains the original atom number of the Hubbard atom.
    integer, allocatable, dimension(:) :: hub_atom_orig_atom

    ! Proc where atom is to be found:
    ! List of atoms on all procs
    integer, allocatable, dimension(:,:) :: atoms_on_proc
    ! cks: the new order of the atoms after distribution to procs
    ! cks: i.e first come all atoms of proc0, then all atoms of proc1, etc.
    integer, allocatable, dimension(:) :: orig_atom
    ! pdh: the opposite of orig_atom
    integer, allocatable, dimension(:) :: distr_atom
    ! List of atoms on this proc
    type(element), allocatable, dimension(:) :: elements_on_proc

    ! agreco: information about k-points
    ! number of unique k-points specified in input file
    integer :: num_kpoints
    ! max number of points on any proc
    integer :: max_kpoints_on_proc
    ! number of k-points on each proc
    integer, allocatable :: num_kpoints_on_proc(:)
    ! first k-point on each proc
    integer, allocatable :: first_kpoint_on_proc(:)
    ! proc of each k-point
    integer, allocatable :: proc_of_kpoint(:)
    ! list of k-points on all procs
    integer, allocatable, dimension(:,:) :: kp_indices_on_proc
    ! list of k-points on this proc
    type(KPOINT), allocatable, dimension(:) :: kpoints_on_proc
    ! rc2013: First block-column for this proc (shorthand)
    integer :: my_first_blk
    ! rc2013: Last block-column for this proc (shorthand)
    integer :: my_last_blk
    ! rc2013: info for sparse_mod -- removed to enable multiple par's
    ! rc2013: it might be best to replace these entirely with sparse
    ! functions. Easier to use + better consistency
    integer, allocatable :: first_elem_on_atom(:,:)
    integer, allocatable :: first_elem_on_proc(:,:)
    integer, allocatable :: num_elems_on_atom(:,:)
    integer, allocatable :: num_elems_on_proc(:,:)
    integer, allocatable :: atom_of_elem(:,:)
    integer, allocatable :: proc_of_elem(:,:)
    ! rc2013: identifier for this parallel strategy
    integer :: par_index
    ! rc2013: first proc allocated to this par
    integer :: first_proc
    integer :: num_procs

  end type PARAL_INFO

  public :: PARAL_INFO

  type(PARAL_INFO), public, pointer :: pub_par

  ! Private variables

  logical :: atoms_distributed = .false.               ! Flag whether atoms
                                                       ! have been distributed

  ! agreco: whether k-points have been distributed
  logical :: kpoints_distributed = .false.

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_distr_atoms(par,elements,cell)

    !=========================================================================!
    ! This subroutine sorts a list of elements according to a Hilbert space-  !
    ! filling curve so as to maximise the overlap between spheres on the same !
    ! proc, and then distributes the atoms across the procs.                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   elements (input) : The list of elements read from the input file      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   par%max_atoms_on_proc  : Maximum number of atoms on any proc          !
    !   par%num_atoms_on_proc  : Number of atoms on each proc                 !
    !   par%atoms_on_proc      : A list of atoms on each proc                 !
    !   par%elements_on_proc   : A list of elements on this proc              !
    !   par%first_atom_on_proc : First atom on each proc                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   constants       : For pi                                              !
    !   geometry        : For the point type and dot product                  !
    !   simulation_cell : For the CELL_INFO type                              !
    !   ion             : For the element type                                !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   fcoord   : list of fractional coordinates of the atoms wrt lattice    !
    !   gcoord   : list of atomic coordinates on a coarse grid                !
    !   graycode : the Gray codes for the atoms used for sorting according to !
    !              a space-filling curve                                      !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   par%nat > 0                                                           !
    !   size(elements) == par%nat                                             !
    !   pub_comms_initialised == .true.                                       !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/7/03                                        !
    ! Modified by Quintin Hill to use cell on 17/10/2008.                     !
    ! Modified by Nicholas Hine to take cell as argument on 12/08/2014.       !
    !=========================================================================!

    use comms, only: pub_comms_initialised, pub_on_root, &
        pub_my_proc_id, pub_total_num_procs
    use constants, only : PI, LONG
    use geometry, only: operator(.DOT.)
    use hubbard_init, only: h_species
    use ion, only: element
    use rundat, only: pub_use_space_filling_curve, pub_output_detail, &
        pub_hubbard, pub_extend_ngwf, pub_devel_code
    use simulation_cell, only : CELL_INFO
    use utils, only: utils_assert, utils_devel_code, utils_heapsort

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(inout) :: par
    type(element), intent(in) :: elements(:)  ! Unsorted atomic information
    type(CELL_INFO), intent(in) :: cell

    ! Local variables
    integer :: iat                            ! Atom loop counter
    integer :: hub_species
    integer :: orig_atom_count
    integer :: num
    integer :: ibit                           ! Bit loop counter
    integer :: proc                           ! Proc counter
    integer :: proc_iat
    integer :: ngrid(3),log2grid(3),abc       ! Grid sizes and their log_2s
    integer :: global_atom                    ! Counter for atom in system
    integer :: hub_atom, find_hub_atom        ! ddor: Counter used for distributing Hubbard atoms
    integer :: atom_count                     ! Counter for atoms
    integer, allocatable :: gcoord(:,:)       ! Coordinates in terms of grid
    integer, allocatable, dimension(:) :: sort_index ! Sorting index for atoms
    integer, allocatable, dimension(:) :: nfuncs     ! Number of functions/atom
    integer, allocatable, dimension(:) :: num_funcs_on_proc ! Num funcs/proc
    integer(kind=LONG), allocatable :: graycode(:) ! Gray code for atoms

    real(kind=DP) :: recip_twopi              ! 1 / (2 pi)
    real(kind=DP) :: av_funcs_on_proc         ! Average number of NGWFs per proc
    real(kind=DP) :: excess_funcs
    real(kind=DP), allocatable :: fcoord(:,:) ! Atomic coordinates in fractions
                                              ! of lattice vectors
    logical :: old_load_balance, debug_print
    real(kind=DP) :: av_funcs_on_proc_current
    integer :: num_funcs_old, num_funcs_current, num_current
    integer :: maxproc, skip
    integer :: max_procs_for_par

    ! Check arguments
    call utils_assert(par%nat > 0, &
         'Error in parallel_strategy_distr_atoms: no atoms in cell!')

    call utils_assert(size(elements) == par%nat, &
         'Error in parallel_strategy_distr_atoms: incompatible nat')

    ! Check comms module intialised
    call utils_assert(pub_comms_initialised, &
         'Error in parallel_strategy_distr_atoms: comms module not initialised')

    ! Check if user wants old logic for load balancing
    ! which can be activated with devel_code = PAR:OLD_LOAD_BALANCE=T:PAR
    old_load_balance = utils_devel_code(.false., 'PAR', 'OLD_LOAD_BALANCE', &
         pub_devel_code, no_bcast=.true., no_warn=.true.)

    ! Check if user wants debug printing
    ! which can be activated with devel_code = PAR:DEBUG_PRINT=T:PAR
    debug_print = utils_devel_code(.false., 'PAR', 'DEBUG_PRINT', &
         pub_devel_code, no_bcast=.true., no_warn=.true.)
    debug_print = debug_print.and.pub_on_root

    ! Reallocate module arrays if necessary
    call internal_allocate_mod_1

    ! cks: the "pub_use_space_filling_curve" is a logical that comes from the input
    if (pub_on_root .and. pub_output_detail >= VERBOSE) then
      write(stdout,'(/a)',advance='no') '... '
      if (.not. pub_use_space_filling_curve) write(stdout,'(a)',advance='no') 'not '
      write(stdout,'(a)') 'using space-filling curve'
      if (old_load_balance) write(stdout,'(a)') '... using old load balance algorithm'
    end if

    ! Allocate work arrays (and initialise sort index for atoms)
    call internal_allocate_work(pub_use_space_filling_curve)

    !---------------------------------------------!
    ! Sort atoms according to space-filling curve !
    !  -- only if required --                     !
    !---------------------------------------------!

    ! ebl: note for future DMFT development - mask_dimer code removed. May need to
    !      reintroduce this here.

    if (pub_use_space_filling_curve) then

       ! Convert Cartesian coordinates from cell into fractional coordinates in
       ! the interval [0,1)
       do iat=1,par%nat
          fcoord(1,iat) = elements(iat)%centre .dot. cell%b1
          fcoord(2,iat) = elements(iat)%centre .dot. cell%b2
          fcoord(3,iat) = elements(iat)%centre .dot. cell%b3
       end do
       recip_twopi = 0.5_DP / pi
       fcoord = fcoord * recip_twopi
       fcoord = modulo(fcoord,1.0_DP)

       ! Set up a grid which will be used to 'coarsen' the atomic coordinates
       ! which can then be used in the sorting procedure - for simplicity choose
       ! the nearest power of 2 above the simulation cell coarse grid
       ngrid(1) = cell%total_pt1
       ngrid(2) = cell%total_pt2
       ngrid(3) = cell%total_pt3
       do abc=1,3
          log2grid(abc) = int(log(real(ngrid(abc),kind=DP))/log(2.0_DP))+1
          ngrid(abc) = 2**log2grid(abc)
       end do

       ! 'Coarsen' the fractional coordinates onto this grid
       do iat=1,par%nat
          do abc=1,3
             gcoord(abc,iat) = int(fcoord(abc,iat)*ngrid(abc))
          end do
       end do

       ! Construct Gray code for each atom (interlace binary digits from the
       ! three dimensions) - use long integers to allow grid sizes > 1024
       graycode = 0_LONG

       do iat=1,par%nat   ! Loop over all atoms in cell

          ! agreco: same graycode for atoms aligned along an ext direction
          if (.not.pub_extend_ngwf(1)) then
             do ibit=0,log2grid(1)-1
                if (btest(gcoord(1,iat),ibit)) graycode(iat) = &
                     ior(graycode(iat),ishft(4_LONG,int(3*ibit,kind=LONG)))
             end do
          end if
          if (.not.pub_extend_ngwf(2)) then
             do ibit=0,log2grid(2)-1
                if (btest(gcoord(2,iat),ibit)) graycode(iat) = &
                     ior(graycode(iat),ishft(2_LONG,int(3*ibit,kind=LONG)))
             end do
          end if
          if (.not.pub_extend_ngwf(3)) then
             do ibit=0,log2grid(3)-1
                if (btest(gcoord(3,iat),ibit)) graycode(iat) = &
                     ior(graycode(iat),ishft(1_LONG,int(3*ibit,kind=LONG)))
             end do
          end if

       end do  ! End of loop over all atoms

       ! Sort atoms by their Gray codes
       call utils_heapsort(par%nat,graycode,sort_index)

    end if

    !-----------------------------------------------------------------------!
    ! Distribute atoms across procs, balancing number of funcs on each proc !
    !-----------------------------------------------------------------------!


    ! Fill nfuncs array with number of functions on each atom
    ! Can change this value to distribute by other types of function
    do iat=1,par%nat
       nfuncs(iat) = elements(iat)%nfunctions
    end do

    ! Find largest number of NGWFs on any atom
    num_funcs_current = maxval(nfuncs(:)) ! Current NGWFs / atom size
    num = 0                   ! Running count of NGWFs of current size or above
    par%num_atoms_on_proc(:) = 0 ! Counts number of atoms on each proc
    par%proc_of_atom(:) = -1
    num_funcs_on_proc(:) = 0     ! Counts number of functions on each proc
    ! rc2013: allow for a subset of procs to be allocated to par
    max_procs_for_par = par%first_proc + par%num_procs

    ! Descending loop over different numbers of NGWFs per atom
    ! This ensures larger atoms (by NGWF count) are distributed first, followed
    ! by smaller ones (essentially a "Modified Greedy Algorithm" for load
    ! balance)
    do
       ! Count total number of NGWFs on atoms with this many NGWFs
       num_current = 0
       do iat=1,par%nat
          if ((nfuncs(iat)==num_funcs_current).or.old_load_balance) then
             num = num + nfuncs(iat)
             num_current = num_current + nfuncs(iat)
          end if
       end do

       ! Calculate expected number of functions per proc for atoms of this size
       ! rc2013: no. of procs allocated to par may be less than total
       av_funcs_on_proc = num / real(par%num_procs,kind=DP)
       av_funcs_on_proc_current = num_current / real(par%num_procs,kind=DP)

       proc = par%first_proc     ! Stores current proc being filled up
       skip = 0

       ! Work out number of atoms on each proc - aim is to balance the number
       ! of NGWFs per proc to balance the work
       do iat=1,par%nat

          ! Skip atom if it is not the size we are currently working on
          if ((nfuncs(sort_index(iat))/=num_funcs_current) &
               .and.(.not.old_load_balance)) cycle

          ! Calculate what this proc's overall surplus/deficit of NGWFs is
          excess_funcs = num_funcs_on_proc(proc)-av_funcs_on_proc

          do
             ! If less than 1 atom per proc, check first if this is the best
             ! choice of atom within the next nfuncs / avfuncs procs
             ! If so, place atom then start a counter to skip the
             ! next (nfuncs / avfuncs) procs
             if ((av_funcs_on_proc_current <= nfuncs(sort_index(iat))) &
                  .and.(skip==0).and.(.not.old_load_balance)) then
                maxproc = proc+int(nfuncs(sort_index(iat))/av_funcs_on_proc_current)
                maxproc = min(maxproc,max_procs_for_par-1)
                if (.not.any(num_funcs_on_proc(proc:maxproc) < &
                     num_funcs_on_proc(proc))) then
                   skip = nfuncs(sort_index(iat)) / av_funcs_on_proc_current
                   exit
                end if
             end if

             ! Exit if placing this atom on this proc improves balance
             if ((abs(excess_funcs) > abs(excess_funcs+nfuncs(sort_index(iat)))) &
                  .and.(skip==0)) exit
             ! Exit if no more procs left so we must place atom here
             if (proc > max_procs_for_par-1) exit

             ! Otherwise move on to next proc
             proc = proc + 1
             if (skip>0) skip = skip - 1

             ! Update the expected average to take history into account:
             ! remaining number of NGWFs / remaining number of procs
             if (proc > par%first_proc) then
                av_funcs_on_proc = &
                  (num-sum(num_funcs_on_proc(0:proc-1))) / &
                  real(max_procs_for_par-proc,kind=DP)
                av_funcs_on_proc_current = num_current / &
                     real(max_procs_for_par-proc,kind=DP)
             end if
             excess_funcs = num_funcs_on_proc(proc)-av_funcs_on_proc
          end do

          ! Put the atom on the current proc
          !if (pub_on_root) print '(5(a,i5))','placing atom ',iat, &
          !     ' (',nfuncs(sort_index(iat)), &
          !     ') on proc ',proc, '(skip=',skip,')', av_funcs_on_proc_current
          num_funcs_on_proc(proc) = num_funcs_on_proc(proc) + &
               nfuncs(sort_index(iat))
          par%num_atoms_on_proc(proc) = par%num_atoms_on_proc(proc) + 1
          par%proc_of_atom(sort_index(iat)) = proc
          num_current = num_current - nfuncs(sort_index(iat))

          ! Update the expected average to take history into account:
          ! remaining number of NGWFs / remaining number of procs
          if (proc > 0) then
             av_funcs_on_proc = &
                  (num-sum(num_funcs_on_proc(0:proc-1))) / &
                  real(max_procs_for_par-proc,kind=DP)
             av_funcs_on_proc_current = num_current / &
                  real(max_procs_for_par-proc,kind=DP)
          end if

       end do  ! End of loop over atoms

       ! Now find the next NGWF / atom count to look at
       num_funcs_old = num_funcs_current
       num_funcs_current = -1
       do iat=1,par%nat
          if ((nfuncs(iat)>num_funcs_current).and.(nfuncs(iat)<num_funcs_old)) &
               num_funcs_current = nfuncs(iat)
       end do

       ! Optional load-balance debugging code
       if (debug_print) then
          print *,'nfuncs(:)'
          print *, nfuncs(sort_index(:))
          print *,'num_funcs_on_proc(:)'
          print *, num_funcs_on_proc(:)
          print *,'proc_of_atom(:)'
          print *, par%proc_of_atom(sort_index(:))
       end if

       ! Quit outer loop if no new NGWF/atom count is found
       if ((num_funcs_current==-1).or.(old_load_balance)) exit

    end do

    ! Find maximum number of atoms on any one proc
    par%max_atoms_on_proc = maxval(par%num_atoms_on_proc)

    ! Initialise array indicating first atom on each proc
    global_atom = 1
    par%first_atom_on_proc = 0
    ! rc2013: only loop over the procs this par has available
    do proc=par%first_proc,par%first_proc+par%num_procs-1
       par%first_atom_on_proc(proc) = global_atom
       global_atom = global_atom + par%num_atoms_on_proc(proc)
    end do
    par%first_atom_on_proc(par%first_proc+par%num_procs) = par%nat + 1

    ! Reallocate module arrays which depend upon par%max_atoms_on_proc
    ! if necessary
    call internal_allocate_mod_2

    ! Make up list of atoms on each proc
    do proc=par%first_proc,par%first_proc+par%num_procs-1
       proc_iat = 1
       do iat=1,par%nat
          if (par%proc_of_atom(sort_index(iat))==proc) then
             par%atoms_on_proc(proc_iat,proc) = sort_index(iat)
             proc_iat = proc_iat + 1
          end if
       end do
    end do

    ! cks: Global order of atoms as distributed on procs.
    ! cks: This is the order in which they will be arranged in matrices, etc.
    atom_count =0
    do proc=par%first_proc,par%first_proc+par%num_procs-1
       do iat=1,par%num_atoms_on_proc(proc)
          atom_count = atom_count+1

          orig_atom_count           = par%atoms_on_proc(iat,proc)
          par%orig_atom(atom_count) = orig_atom_count
          par%distr_atom(orig_atom_count) = atom_count

       end do
    end do

    ! Make up list of atoms on this proc
    par%elements_on_proc(:)%species_number = 0
    par%elements_on_proc(:)%pspecies_number = 0
    do iat=1,par%num_atoms_on_proc(pub_my_proc_id)
       par%elements_on_proc(iat) = &
           elements(par%atoms_on_proc(iat,pub_my_proc_id))
    end do

    ! ddor: Distribute Hubbard DFT+U atoms and projectors over procs if necessary,
    !       Hubbard atoms are distributed as per their host atoms,
    !       and it is necessary to have par%orig_atom at this stage.
    if (pub_hubbard) then

       hub_atom = 0
       do hub_species=1,par%num_hub_species
          do iat=1,par%nat
             if (elements(iat)%species_id == &
                  h_species(hub_species)%hub_species) then
                hub_atom = hub_atom + 1
                par%hub_atom_orig_atom(hub_atom) = iat
             end if
          end do
       end do

       par%num_hub_atoms_on_proc = 0

       do proc=par%first_proc,par%first_proc+par%num_procs-1

          ! ddor: Loop over Hubbard atoms
          do hub_atom = 1,par%nat_hub
             ! Loop over atoms on this proc
             do find_hub_atom=par%first_atom_on_proc(proc), &
                  par%first_atom_on_proc(proc+1)-1
                ! If it's the Hubbard atom we want
                if ( par%hub_atom_orig_atom(hub_atom) == par%orig_atom(find_hub_atom) ) then
                   par%num_hub_atoms_on_proc(proc) = par%num_hub_atoms_on_proc(proc) + 1

                   par%proc_of_hub_atom(hub_atom) = proc
                   par%hub_atoms_on_proc(par%num_hub_atoms_on_proc(proc),proc) = hub_atom
                endif
             enddo
          enddo

       end do

    endif

    ! Deallocate work arrays
    call internal_deallocate_work(pub_use_space_filling_curve)

    ! Set flag to show atoms have been distributed
    atoms_distributed = .true.

  contains

    !--------------------------------------------------------------------------

    subroutine internal_allocate_mod_1

      !=======================================================================!
      ! This subroutine reallocates module arrays as required by the parent   !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      ! Reorganised by Nicholas Hine, 18/05/09                                !
      !=======================================================================!

      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Local variables
      integer :: ierr     ! Error flag

      ! Deallocate module arrays if necessary
      if (allocated(par%num_atoms_on_proc)) then
         if (size(par%num_atoms_on_proc) /= pub_total_num_procs) then
            deallocate(par%num_atoms_on_proc,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','par%num_atoms_on_proc',ierr)
         end if
      end if
      if (allocated(par%first_atom_on_proc)) then
         if (size(par%first_atom_on_proc) /= pub_total_num_procs+1) then
            deallocate(par%first_atom_on_proc,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','par%first_atom_on_proc',ierr)
         end if
      end if
      if (allocated(par%proc_of_atom)) then
         if (size(par%proc_of_atom) /= par%nat) then
            deallocate(par%proc_of_atom,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','par%proc_of_atom',ierr)
         end if
      end if
      if (allocated(par%orig_atom)) then
         if (size(par%orig_atom) /= par%nat) then
            deallocate(par%orig_atom,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','par%orig_atom',ierr)
         end if
      end if
      if (allocated(par%distr_atom)) then
         if (size(par%distr_atom) /= par%nat) then
            deallocate(par%distr_atom,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','par%distr_atom',ierr)
         end if
      end if
      if (pub_hubbard) then ! ddor
         if (allocated(par%num_hub_atoms_on_proc)) then
            if (size(par%num_hub_atoms_on_proc) /= pub_total_num_procs) then
               deallocate(par%num_hub_atoms_on_proc,stat=ierr)
               call utils_dealloc_check('internal_allocate_mod_1 &
                    &(parallel_strategy_distr_atoms)','par%num_hub_atoms_on_proc',ierr)
            end if
         end if
         if (allocated(par%proc_of_hub_atom)) then
            if (size(par%proc_of_hub_atom) /= par%nat_hub) then
               deallocate(par%proc_of_hub_atom,stat=ierr)
               call utils_dealloc_check('internal_allocate_mod_1 &
                    &(parallel_strategy_distr_atoms)','par%proc_of_hub_atom',ierr)
            end if
         end if
         if (allocated(par%hub_atoms_on_proc)) then
            if (size(par%hub_atoms_on_proc,1) /= par%nat_hub .or. &
                 size(par%hub_atoms_on_proc,2) /= pub_total_num_procs) then
               deallocate(par%hub_atoms_on_proc,stat=ierr)
               call utils_dealloc_check('internal_allocate_mod_2 &
                    &(parallel_strategy_distr_atoms)','par%hub_atoms_on_proc',ierr)
            end if
         end if
         if (allocated(par%hub_atom_orig_atom)) then
            if (size(par%hub_atom_orig_atom) /= par%nat_hub) then
               deallocate(par%hub_atom_orig_atom,stat=ierr)
               call utils_dealloc_check('internal_allocate_mod_1 &
                   &(parallel_strategy_distr_atoms)','par%hub_atom_orig_atom',ierr)
            end if
         end if
      end if
      ! Allocate module arrays if necessary
      if (.not.allocated(par%num_atoms_on_proc)) then
         allocate(par%num_atoms_on_proc(0:pub_total_num_procs-1),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','par%num_atoms_on_proc',ierr)
      end if
      if (.not.allocated(par%first_atom_on_proc)) then
         allocate(par%first_atom_on_proc(0:pub_total_num_procs),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','par%first_atom_on_proc',ierr)
      end if
      if (.not.allocated(par%proc_of_atom)) then
         allocate(par%proc_of_atom(par%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','par%proc_of_atom',ierr)
      end if
      if (.not.allocated(par%orig_atom)) then
         allocate(par%orig_atom(par%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','par%orig_atom',ierr)
      end if
      if (.not.allocated(par%distr_atom)) then
         allocate(par%distr_atom(par%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','par%distr_atom',ierr)
      end if
      if (pub_hubbard) then ! ddor
         if (.not.allocated(par%num_hub_atoms_on_proc)) then
            allocate(par%num_hub_atoms_on_proc(0:pub_total_num_procs-1),stat=ierr)
            call utils_alloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','par%num_hub_atoms_on_proc',ierr)
         end if
         if (.not.allocated(par%proc_of_hub_atom)) then
            allocate(par%proc_of_hub_atom(par%nat_hub),stat=ierr)
            call utils_alloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','par%proc_of_hub_atom',ierr)
         end if
         if (.not. allocated(par%hub_atoms_on_proc)) then
            allocate(par%hub_atoms_on_proc(par%nat_hub, &
                 0:pub_total_num_procs-1),stat=ierr)
            call utils_alloc_check('internal_allocate_mod_2 &
                 &(parallel_strategy_distr_atoms)','par%hub_atoms_on_proc',ierr)
            ! ddor: initialise to something obvious
            par%hub_atoms_on_proc = -1
         end if
         if (.not.allocated(par%hub_atom_orig_atom)) then
            allocate(par%hub_atom_orig_atom(par%nat_hub),stat=ierr)
            call utils_alloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','par%hub_atom_orig_atom',ierr)
         end if
      end if

    end subroutine internal_allocate_mod_1

    !--------------------------------------------------------------------------

    subroutine internal_allocate_work(sfc)

      !=======================================================================!
      ! This subroutine allocates work arrays as required by the parent       !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! sfc (input) : flag whether to use space-filling curve or not          !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      use utils, only: utils_alloc_check

      implicit none

      ! Arguments
      logical, intent(in) :: sfc

      ! Local variables
      integer :: ierr     ! Error flag
      integer :: iat      ! Atom loop variable

      ! Allocate the index for the sorting of the atoms
      allocate(sort_index(par%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_work &
           &(parallel_strategy_distr_atoms)','sort_index',ierr)
      allocate(nfuncs(par%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_work &
           &(parallel_strategy_distr_atoms)','nfuncs',ierr)
      allocate(num_funcs_on_proc(0:pub_total_num_procs-1),stat=ierr)
      call utils_alloc_check('internal_allocate_work &
           &(parallel_strategy_distr_atoms)','num_funcs_on_proc',ierr)

      ! Allocate work arrays for space-filling curve if required
      if (sfc) then
         allocate(gcoord(3,par%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_work &
              &(parallel_strategy_distr_atoms)','gcoord',ierr)
         allocate(graycode(par%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_work &
              &(parallel_strategy_distr_atoms)','graycode',ierr)
         allocate(fcoord(3,par%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_work &
              &(parallel_strategy_distr_atoms)','fcoord',ierr)
      else
         ! cks: if no space-filling curve is used, the
         ! cks: sorting index needs to be initialised
         do iat=1, par%nat
            sort_index(iat)=iat
         end do
      endif


    end subroutine internal_allocate_work

    !--------------------------------------------------------------------------

    subroutine internal_allocate_mod_2

      !=======================================================================!
      ! This subroutine reallocates module arrays as required by the parent   !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Local variables
      integer :: ierr     ! Error flag

      ! Deallocate module arrays which depend upon par%max_atoms_on_proc
      ! if necessary
      if (allocated(par%atoms_on_proc)) then
         if (size(par%atoms_on_proc,1) /= par%max_atoms_on_proc .or. &
              size(par%atoms_on_proc,2) /= pub_total_num_procs) then
            deallocate(par%atoms_on_proc,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_2 &
                 &(parallel_strategy_distr_atoms)','par%atoms_on_proc',ierr)
         end if
      end if
      if (allocated(par%elements_on_proc)) then
         if (size(par%elements_on_proc) /= par%max_atoms_on_proc) then
            deallocate(par%elements_on_proc,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_2 &
                 &(parallel_strategy_distr_atoms)','par%elements_on_proc',ierr)
         end if
      end if

      ! Allocate module arrays which depend upon par%max_atoms_on_proc
      ! if necessary
      if (.not. allocated(par%atoms_on_proc)) then
         allocate(par%atoms_on_proc(par%max_atoms_on_proc, &
              0:pub_total_num_procs-1),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_2 &
              &(parallel_strategy_distr_atoms)','par%atoms_on_proc',ierr)
         ! cks: initialise to something obvious
         par%atoms_on_proc = -1
      end if
      if (.not. allocated(par%elements_on_proc)) then
         allocate(par%elements_on_proc(par%max_atoms_on_proc),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_2 &
              &(parallel_strategy_distr_atoms)','par%elements_on_proc',ierr)
      end if

    end subroutine internal_allocate_mod_2

    !--------------------------------------------------------------------------

    subroutine internal_deallocate_work(sfc)

      !=======================================================================!
      ! This subroutine deallocates work arrays as required by the parent     !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! sfc (input) : flag whether to use space-filling curve or not          !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      use utils, only: utils_dealloc_check

      implicit none

      ! Arguments
      logical, intent(in) :: sfc

      ! Local variables
      integer :: ierr     ! Error flag

      ! Deallocate work arrays required by space-filling-curve
      if (sfc) then
         deallocate(fcoord,stat=ierr)
         call utils_dealloc_check('internal_deallocate_work &
              &(parallel_strategy_distr_atoms)','fcoord',ierr)
         deallocate(graycode,stat=ierr)
         call utils_dealloc_check('internal_deallocate_work &
              &(parallel_strategy_distr_atoms)','graycode',ierr)
         deallocate(gcoord,stat=ierr)
         call utils_dealloc_check('internal_deallocate_work &
              &(parallel_strategy_distr_atoms)','gcoord',ierr)
      end if

      ! Deallocate sorting index
      deallocate(num_funcs_on_proc,stat=ierr)
      call utils_dealloc_check('internal_deallocate_work &
           &(parallel_strategy_distr_atoms)','num_funcs_on_proc',ierr)
      deallocate(nfuncs,stat=ierr)
      call utils_dealloc_check('internal_deallocate_work &
           &(parallel_strategy_distr_atoms)','nfuncs',ierr)
      deallocate(sort_index,stat=ierr)
      call utils_dealloc_check('internal_deallocate_work &
           &(parallel_strategy_distr_atoms)','sort_index',ierr)

    end subroutine internal_deallocate_work

  end subroutine parallel_strategy_distr_atoms

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_distr_funcs(num, nfuncs_orig, sort_index, par, &
       first_func_on_proc, num_funcs_on_proc, first_func_on_atom, &
       num_funcs_on_atom, proc_of_func, atom_of_func, max_funcs_on_proc, &
       max_funcs_on_atom)

    !=========================================================================!
    ! This subroutine takes a list of numbers of functions on each atom and a !
    ! sorting index and fills the arrays describing the distribution of       !
    ! functions over procs and atoms.                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   num (input)         : The total number of functions being distributed.!
    !   nfuncs_orig (input) : The number of functions on each atom in the     !
    !                         order given in the input file / elements array. !
    !   sort_index (input)  : The index of atoms after sorting over over the  !
    !                         procs of the calculation.                       !
    !   first_func_on_proc,                                                   !
    !   num_funcs_on_proc,                                                    !
    !   first_func_on_atom,                                                   !
    !   num_funcs_on_atom,                                                    !
    !   proc_of_func,                                                         !
    !   atom_of_func (output) : Arrays describing distribution of functions   !
    !                           over procs and atoms (meanings as above).     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   par%first_atom_on_proc : First atom on each proc                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   comms           : For pub_total_num_procs                             !
    !   simulation_cell : For the CELL_INFO type                              !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   par%nat > 0                                                      !
    !   par%first_atom_on_proc to be initialised and filled                   !
    !   pub_comms_initialised == .true.                                       !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine 9th July 2009.                            !
    !=========================================================================!

    use comms, only: pub_total_num_procs

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(in) :: par
    integer, intent(in) :: num
    integer, intent(in) :: nfuncs_orig(1:par%nat)
    integer, intent(in) :: sort_index(1:par%nat)
    integer, intent(out) :: first_func_on_proc(0:pub_total_num_procs)
    integer, intent(out) :: num_funcs_on_proc(0:pub_total_num_procs-1)
    integer, intent(out) :: first_func_on_atom(1:par%nat)
    integer, intent(out) :: num_funcs_on_atom(1:par%nat)
    integer, intent(out) :: proc_of_func(1:num)
    integer, intent(out) :: atom_of_func(1:num)
    integer, intent(out) :: max_funcs_on_proc
    integer, intent(out) :: max_funcs_on_atom

    ! Locals
    integer :: proc
    integer :: iat, ifunc, jfunc, nfuncs

    ! Work out number of functions on each proc
    proc = par%first_proc - 1 ! Stores current proc being filled up
    num_funcs_on_proc(:) = 0  ! Counts number of NGWFs on each proc
    ifunc = 1

    ! Loop over the atoms in their sorted order
    do iat=1,par%nat
       nfuncs = nfuncs_orig(sort_index(iat))

       ! Check if we are on a new proc
       if (iat == par%first_atom_on_proc(proc+1)) then
          proc = proc + 1
          first_func_on_proc(proc) = ifunc
       end if

       ! Add these functions to the current proc and atom
       num_funcs_on_proc(proc) = num_funcs_on_proc(proc) + nfuncs
       num_funcs_on_atom(iat) = nfuncs

       ! Record first function on this atom (if there are any functions)
       if (nfuncs > 0) then
          first_func_on_atom(iat) = ifunc
       else
          first_func_on_atom(iat) = 0
       end if

       ! Record the proc and atom of these functions
       do jfunc=ifunc,ifunc+nfuncs-1
          proc_of_func(jfunc) = proc
          atom_of_func(jfunc) = iat
       end do

       ! Increment current function index
       ifunc = ifunc + nfuncs

    end do  ! End of loop over atoms

    ! Record last function + 1 at end of array
    first_func_on_proc(par%first_proc+par%num_procs) = num + 1

    ! Find maximum number of functions on any proc
    max_funcs_on_proc = maxval(num_funcs_on_proc)

    ! Find maximum number of functions on any atom
    max_funcs_on_atom = maxval(num_funcs_on_atom)

  end subroutine parallel_strategy_distr_funcs

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_check_atoms(par,elements,cell)

    use comms, only: comms_reduce, pub_my_proc_id
    use constants, only: DP, TWO_PI, stdout
    use geometry, only: operator(.DOT.), operator(*), operator(+), &
        geometry_distance
    use ion, only: element
    use rundat, only: pub_check_atoms
    use simulation_cell, only : CELL_INFO
    use utils, only: utils_flush, utils_assert

    !====================================================================!
    ! Make sure the positions of all atoms are inside the simulation     !
    ! cell and there no nearly-overlapping atoms.                        !
    !--------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000 as                        !
    ! services_check_atom_positions and revised on 30/5/2001.            !
    ! Modified by Chris-Kriton Skylaris on 12/10/2004 to also check      !
    ! for nearly-overlapping atoms.                                      !
    ! Modified by Nick Hine on 06/09/2008 to only check the atoms on     !
    ! this proc, and moved to parallel_strategy_mod.F90                  !
    ! Modified by Nick Hine on 30/10/2008 to report which atoms are      !
    ! outside the unit cell in a more informative way.                   !
    ! Modified to prevent hanging when only some of the procs detect an  !
    ! atom outside of simulation cell or two atoms too close, by Nick    !
    ! Hine on 18/09/2009.                                                !
    !====================================================================!

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(inout) :: par
    type(ELEMENT), intent(in) :: elements(par%nat)
    type(CELL_INFO), intent(in) :: cell

    ! cks: internal variables
    integer :: atom, orig_atom
    integer :: row, local_row
    integer :: col, per_col
    integer :: per_col_start, per_col_end
    integer :: loc1
    integer :: loc2
    integer :: loc3
    integer :: abort_count
    real(kind=DP) :: f1, f2, f3 ! rab207: fractional co-ordinates
    logical :: in_a1, in_a2, in_a3
    real(kind=DP) :: row_col_dist
    real(kind=DP), parameter :: tolerance = 1.0e-10_DP


    abort_count = 0

    ! cks: make sure each atom is inside the simulation cell
    ! nmdh: only check atoms on this proc
    do atom=1,par%num_atoms_on_proc(pub_my_proc_id)

       ! rab207: find fractional position along a1, a2, a3
       f1 = (par%elements_on_proc(atom)%centre.DOT.cell%b1) / TWO_PI
       f2 = (par%elements_on_proc(atom)%centre.DOT.cell%b2) / TWO_PI
       f3 = (par%elements_on_proc(atom)%centre.DOT.cell%b3) / TWO_PI

       in_a1 = (f1 .ge. -tolerance) .and. (f1 .le. 1.0_DP + tolerance)
       in_a2 = (f2 .ge. -tolerance) .and. (f2 .le. 1.0_DP + tolerance)
       in_a3 = (f3 .ge. -tolerance) .and. (f3 .le. 1.0_DP + tolerance)

       if (.not. in_a1 .or. .not. in_a2 .or. .not. in_a3) then
          orig_atom = &
               par%orig_atom(atom + par%first_atom_on_proc(pub_my_proc_id) - 1)
          write(stdout,'(a,i6,a)') &
               'Error in parallel_strategy_check_atoms: atom number', &
               orig_atom,' is not in simulation cell.'
          write(stdout,'(2a,3f16.8)') par%elements_on_proc(atom)%symbol, &
               'abs position (bohr): ', &
               par%elements_on_proc(atom)%centre%x, &
               par%elements_on_proc(atom)%centre%y, &
               par%elements_on_proc(atom)%centre%z
          write(stdout,'(4x,a,3f16.8)') 'fractional position: ', f1, f2, f3
          call utils_flush
          abort_count = abort_count + 1
       end if

    enddo


    ! Abort if any pairs of atoms are too close
    call comms_reduce('SUM',abort_count)
    call utils_assert(abort_count == 0, &
         'Error in parallel_strategy_check_atoms: the following number of atoms&
         & are outside the simulation cell: ', abort_count)

    abort_count = 0

    if (pub_check_atoms) then
       abort_count = 0

       ! cks: check that there are no atoms that are almost overlapping
       ! ndmh: share checking equally over procs:
       ! ndmh:   only check rows on this column, and check same number
       ! ndmh:   of cols regardless of row, by modulo of periodic col
       do local_row=1,par%num_atoms_on_proc(pub_my_proc_id)
          row = local_row + par%first_atom_on_proc(pub_my_proc_id) - 1

          ! nmdh: find range of cols to check for this row (periodic)
          per_col_start = row + 1
          per_col_end = row + par%nat/2
          ! ndmh: avoid double counting for even-valued nat
          if (per_col_end > par%nat) then
             per_col_end = per_col_end - 1 + modulo(par%nat,2)
          end if
          do per_col=per_col_start,per_col_end

             ! ndmh: find real value of col
             col = modulo(per_col - 1,par%nat) + 1

             ! cks: initialise for current_pair
             row_col_dist =huge(1.0_DP)

             ! cks: 27 cases need to be considered for each distinct
             ! cks: pair of atoms because of periodic boundary conditions.
        xyz: do loc1 =-1, 1
                do loc2 =-1, 1
                   do loc3 =-1, 1

                      ! cks: accumulate minimum distance
                      row_col_dist=min(geometry_distance(elements(row)%centre, &
                           elements(col)%centre &
                           + real(loc1, kind=DP)*cell%a1 &
                           + real(loc2, kind=DP)*cell%a2 &
                           + real(loc3, kind=DP)*cell%a3 ), row_col_dist)

                      if (row_col_dist <= 0.5_DP) then

                         write(stdout,*)'ERROR! The distance between atoms', &
                              row,' and ',col
                         write(stdout,*)'is ', row_col_dist, ' bohr. '
                         abort_count = abort_count + 1
                         exit xyz ! rab207: break to prevent printing several times
                      endif

                   enddo
                enddo
             enddo xyz

          enddo
       enddo

       ! Abort if any pairs of atoms are too close
       call comms_reduce('SUM',abort_count)

       call utils_assert(abort_count == 0, &
            'Error in parallel_strategy_check_atoms: the following number of &
            &atoms are too close: ', abort_count)
    endif


  end subroutine parallel_strategy_check_atoms

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_exit(par)

    !=========================================================================!
    ! This subroutine frees up any memory used by the parallel strategy       !
    ! module.                                                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   par%num_atoms_on_proc : the number of atoms on each proc              !
    !   par%proc_of_atom      : the proc on which each atom is held           !
    !   par%atoms_on_proc     : list of atoms on each proc                    !
    !   par%elements_on_proc  : the elements on this proc                     !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   utils                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/07/03                                       !
    ! Reorganised by Nicholas Hine, 18/05/09                                  !
    ! Modified by Andrea Greco on 04/01/2016 to include k-points arrays       !
    !=========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(inout) :: par

    integer :: ierr ! Error flag

    ! Deallocate list of atoms on this proc
    if (allocated(par%atoms_on_proc)) then
       deallocate(par%atoms_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%atoms_on_proc',ierr)
    end if
    ! Deallocate list of elements on this proc
    if (allocated(par%elements_on_proc)) then
       deallocate(par%elements_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%elements_on_proc',ierr)
    end if

    !ddor: Deallocate information on DFT+U atom distribution
    if (allocated(par%hub_atom_orig_atom)) then
       deallocate(par%hub_atom_orig_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%hub_atom_orig_atom',ierr)
    end if
    if (allocated(par%hub_atoms_on_proc)) then
       deallocate(par%hub_atoms_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%hub_atoms_on_proc',ierr)
    end if
    if (allocated(par%proc_of_hub_atom)) then
       deallocate(par%proc_of_hub_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%proc_of_hub_atom',ierr)
    end if
    if (allocated(par%num_hub_atoms_on_proc)) then
       deallocate(par%num_hub_atoms_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%num_hub_atoms_on_proc',ierr)
    end if

    ! cks: Deallocate par%orig_atom and par%distr_atom on this proc
    if (allocated(par%distr_atom)) then
       deallocate(par%distr_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%distr_atom',ierr)
    end if
    if (allocated(par%orig_atom)) then
       deallocate(par%orig_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%orig_atom',ierr)
    end if
    ! Deallocate information on atom distribution
    if (allocated(par%proc_of_atom)) then
       deallocate(par%proc_of_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%proc_of_atom',ierr)
    end if
    if (allocated(par%first_atom_on_proc)) then
       deallocate(par%first_atom_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%first_atom_on_proc',ierr)
    end if
    if (allocated(par%num_atoms_on_proc)) then
       deallocate(par%num_atoms_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%num_atoms_on_proc',ierr)
    end if
    ! Deallocate list of k-points indices
    if (allocated(par%kp_indices_on_proc)) then
       deallocate(par%kp_indices_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%kp_indices_on_proc',ierr)
    end if
    ! Deallocate list of k-points on this proc
    if (allocated(par%kpoints_on_proc)) then
       deallocate(par%kpoints_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%kpoints_on_proc',ierr)
    end if
    ! Deallocate information on k-points distribution
    if (allocated(par%proc_of_kpoint)) then
       deallocate(par%proc_of_kpoint,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%kpoint_of_atom',ierr)
    end if
    if (allocated(par%first_kpoint_on_proc)) then
       deallocate(par%first_kpoint_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%first_kpoint_on_proc',ierr)
    end if
    if (allocated(par%num_kpoints_on_proc)) then
       deallocate(par%num_kpoints_on_proc,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'par%num_kpoints_on_proc',ierr)
    end if

  end subroutine parallel_strategy_exit

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_list_overlaps(overlaps,nat,elements,cell, &
       mode1,mode2,radius,tddft_region,multiplier)

    !=========================================================================!
    ! This subroutine calculates the list of overlaps between atoms, either   !
    ! according to the region/core radii or according to an optional fixed    !
    ! radius.                                                                 !
    ! The form of overlap list required is determined by the two mode         !
    ! arguments which can be any of the following options:                    !
    !   'C' - use nonlocal pseudopotential Core radii                         !
    !   'R' - use NGWF Region radii                                           !
    !   'A' - use Conduction NGWF Region radii                                !
    !   'F' - use Fixed radius                                                !
    ! mode1 determines the radius used for the primary atom (the index in the !
    ! module arrays) whereas mode2 determines the radius for the secondary    !
    ! atom (listed as overlapping with the primary atom in the list).         !
    ! A cell list algorithm is used to obtain linear scaling for large        !
    ! systems - these cells are referred to as subcells below to avoid        !
    ! confusion with the unit cell.                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   overlaps (inout) : An OVERLAP_LIST type with details of the overlaps  !
    !   nat (input)      : Number of atoms in elements list                   !
    !   elements (input) : A list of elements with positions and radii        !
    !   cell (input)     : The CELL_INFO type containing the cell parameters  !
    !   mode1 (input)    : Mode for primary atom (see above)                  !
    !   mode2 (input)    : Mode for secondary atom (see above)                !
    !   radius (input)   : The optional fixed radius to be used               !
    !   tddft_region(input) : Optional: Used if the TDDFT kernel is limited   !
    !                         To some region in real space                    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   constants       : For pi                                              !
    !   geometry        : For the point type and dot product                  !
    !   simulation_cell : For the CELL_INFO type                              !
    !   ion             : For the element type                                !
    !   utils           : For checking memory allocation/deallocation         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   fcoord   : list of fractional coordinates of the atoms wrt lattice    !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   par%nat > 0                                                           !
    !   size(elements) == par%nat                                             !
    !   mode1 and mode2 must be one of 'C','R' or 'F'                         !
    !   If either mode1 or mode2 == 'F', radius must be present               !
    !   Sphere radii should not be larger than unit cell                      !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/7/03                                        !
    ! Minor modification to reduce size of arrays by Nick Hine 13/03/08       !
    ! Modified by Quintin Hill to use cell on 17/10/2008.                     !
    ! Modified by Nicholas Hine to take cell as argument on 12/08/2014.       !
    !=========================================================================!

    use comms, only: comms_barrier, comms_bcast, comms_reduce, &
        pub_my_proc_id, pub_total_num_procs, pub_my_rank_in_group, &
        pub_comms_group_size, pub_group_comm, pub_rank_comm, &
        pub_my_rank_in_rank, pub_num_comms_groups
    use constants, only : PI, TWO_PI
    use geometry, only: point,operator(.DOT.),operator(-)
    use ion, only: element
    use rundat, only: pub_tddft_kernel_groups, pub_tddft_kernel_group_nsp,&
        pub_tddft_kernel_ngroups, pub_debug, pub_extend_ngwf,&
        pub_lr_tddft_ct_length,pub_tddft_ct_groups,&
        pub_tddft_ct_group_nsp,pub_tddft_ct_ngroups
    use simulation_cell, only : CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(OVERLAP_LIST), intent(inout) :: overlaps ! List of overlaps
    type(CELL_INFO), intent(in) :: cell
    integer, intent(in)       :: nat              ! Number of atoms
    type(element), intent(in) :: elements(:)      ! List of atoms
    character, intent(in) :: mode1                ! Mode for primary atom
    character, intent(in) :: mode2                ! Mode for secondary atom
    real(kind=DP), intent(in), optional :: radius ! Optional fixed radius
    logical, optional :: tddft_region             ! used if TDDFT kernel is only
                                                  ! defined within a region
    real(kind=dp), intent(in), optional :: multiplier ! ja531-> Used to multiply radii.
                                                      ! ja531-> for FOE sparsity.

    ! Local variables
    logical, parameter :: periodic(3) = (/.true.,.true.,.true./) ! Apply PBC
    integer :: ierr                           ! Error flag
    integer :: iat,jat                        ! Atom loop counters
    integer :: ii                             ! Additional counter
    integer :: idim                           ! Dimension loop counter
    integer :: proc                           ! Proc counter
    integer :: group                          ! Comms group counter
    integer :: isc,iscd(3),isc23              ! Subcell loop counters/labels
    integer :: jsc,jscd(3)                    ! Subcell loop counters/labels
    integer :: ksc,kscd(3),ksc23              ! Subcell loop counters/labels
    integer :: iatsc,jatsc                    ! Atom in subcell loop counters
    integer :: num_subcells(3)                ! Number of subcells in each dim
    integer :: total_num_subcells             ! Total number of subcells
    integer :: max_atoms_subcell              ! Max num atoms in any subcell
    integer :: max_my_overlaps                ! Max num overlaps on any proc
    integer :: novat                          ! Total num overlaps for atom
    integer :: nyrat                          ! Num overlaps for atom on proc
    integer :: iovlap                         ! Overlap loop counter
    real(kind=DP) :: recip_twopi              ! 1 / (2 pi)
    real(kind=DP) :: subcell_size             ! Minimum subcell size
    real(kind=DP) :: spacing(3)               ! Plane spacing of unit cell faces
    real(kind=DP) :: minspacing               ! minimum plane spacing
    real(kind=DP) :: a1(3),a2(3),a3(3)        ! Local copy of lattice vectors
    real(kind=DP) :: fextra(3)                ! Wrap-around for PBC
    real(kind=DP) :: fdiff(3)                 ! Fractional coordinate diffs
    real(kind=DP) :: adiff(3)                 ! Absolute Cartesian diffs
    real(kind=DP) :: dist_sq                  ! Squared distance between atoms
    real(kind=DP) :: cutoff                   ! Cutoff distance
    integer, allocatable :: atom_subcell(:)   ! Which subcell an atom belongs to
    integer, allocatable :: num_atoms_subcell(:) ! Num atoms in each subcell
    integer, allocatable :: subcell_list(:,:) ! Subcell list of atoms
    integer, allocatable :: num_my_overlaps(:)! Num overlaps on this proc
    integer, allocatable :: num_gp_overlaps(:)! Num overlaps in this group
    integer, allocatable :: num_yr_overlaps(:)! Num overlaps on other proc
    integer, allocatable :: my_overlap_list(:,:) ! Overlaps on this proc
    integer, allocatable :: gp_overlap_list(:,:) ! Overlaps in this group
    integer, allocatable :: yr_overlap_list(:,:) ! Overlaps on other proc
    logical, allocatable :: overlapped(:)     ! Flags to avoid multiple overlaps
    real(kind=DP), allocatable :: fcoord(:,:) ! Atomic coordinates in fractions
                                              ! of lattice vectors
    real(kind=DP), allocatable :: radii(:,:)  ! Radii to be used for primary and
                                              ! secondary atoms
    character(len=4)  :: dummy_id             ! tjz07: variables associated with
    logical  :: element_found                 ! TDDFT sparse region
    logical  :: element_found2
    integer  :: i_ngroup
    logical  :: loc_tddft_region
    integer  :: tddft_ngroup_index
    type(point) :: tddft_ct_vec
    logical :: ct_element_found
    logical :: ct_element_found2

    real(kind=DP), parameter :: tolerance = 1.0e-10_DP

    ! set default parameter
    loc_tddft_region=.false.
    if(present(tddft_region)) then
       loc_tddft_region=tddft_region
    endif

    ! Check arguments
    call internal_check_args(present(radius))

    ! Allocate module and work arrays
    call internal_allocate_1

    ! Make up lists of radii for primary and secondary atoms according to
    ! the desired modes
    select case (mode1)
    case('C','c')
       do iat=1,nat
          radii(1,iat) = elements(iat)%max_core_radius
       end do
    case('R','r')
       do iat=1,nat
          radii(1,iat) = elements(iat)%radius
       end do
    case('A','a')
       do iat=1,nat
          radii(1,iat) = elements(iat)%radius_cond
       end do
    case('F','f')
       radii(1,:) = radius
    case('L','l')
       do iat=1,nat
          radii(1,iat) = elements(iat)%max_core_wf_radius
       end do
    end select
    select case (mode2)
    case('C','c')
       do iat=1,nat
          radii(2,iat) = elements(iat)%max_core_radius
       end do
    case('R','r')
       do iat=1,nat
          radii(2,iat) = elements(iat)%radius
       end do
    case('A','a')
       do iat=1,nat
          radii(2,iat) = elements(iat)%radius_cond
       end do
    case('F','f')
       radii(2,:) = radius
    case('L','l')
       do iat=1,nat
          radii(2,iat) = elements(iat)%max_core_wf_radius
       end do
    end select

    !ja531-> if requested, then multiply the radii by multiplier.
    if(present(multiplier)) then
       radii=radii*multiplier
    end if

    ! tjz07: tddft kernel only has overlaps in a predefined region
    ! THIS is only true if only one group of atoms is defined as
    ! the region, ie. if pub_tddft_kernel_ngroups==1. If multiple
    ! atom groups are defined, then there are multiple regions in
    ! the system. In this case, tddft kernel is zero for any
    ! matrix elements corresponding to atoms between regions, but
    ! is not zero for matrix elements within a single region.
    if( loc_tddft_region .and. pub_tddft_kernel_ngroups==1) then
       ! loop over atoms:
       do iat=1, nat
          element_found=.false.
          do ii=1, pub_tddft_kernel_group_nsp(1)
             dummy_id =pub_tddft_kernel_groups(ii,1)
             if(elements(iat)%species_id == dummy_id) then
                element_found=.true.
             endif
          enddo
          if(.not. element_found) then
             radii(1,iat) =0.0_DP
             radii(2,iat) =0.0_DP
          endif
       enddo
    endif


    ! Convert Cartesian coordinates from cell into fractional coordinates in
    ! the interval [0,1)
    do iat=1,nat
       fcoord(1,iat) = elements(iat)%centre .dot. cell%b1
       fcoord(2,iat) = elements(iat)%centre .dot. cell%b2
       fcoord(3,iat) = elements(iat)%centre .dot. cell%b3
    end do

    ! jd: This ensures fcoord(:,:) are *really* within [0,1). Due to floating-precision
    !     rounding errors fcoord(:,:) can wind up very slightly negative on entry here.
    !     The subsequent modulo(,1_DP) operation then makes them exactly 1.0_DP, which
    !     leads to an off-by-one calc of subcell number.
    do iat=1,nat
       do jat = 1,3
           if (fcoord(jat,iat) < 0.0_DP .and. fcoord(jat,iat) > (-tolerance)) then
              fcoord(jat,iat) = 0.0_DP
           elseif (fcoord(jat,iat) > 1.0_DP .and. fcoord(jat,iat) < 1.0_DP+tolerance) then
              fcoord(jat,iat) = 1.0_DP
           end if
       end do
    end do

    recip_twopi = 0.5_DP / pi
    fcoord = fcoord * recip_twopi
    fcoord = modulo(fcoord,1.0_DP)

    ! Make local copies of lattice vectors
    a1(1) = cell%a1%x ; a1(2) = cell%a1%y ; a1(3) = cell%a1%z
    a2(1) = cell%a2%x ; a2(2) = cell%a2%y ; a2(3) = cell%a2%z
    a3(1) = cell%a3%x ; a3(2) = cell%a3%y ; a3(3) = cell%a3%z

    ! Figure out the maximum potential cutoff separation - this defines the
    ! minimum subcell size required
    subcell_size = maxval(radii(1,:)) + maxval(radii(2,:))

    ! Avoid subcell sizes which are unphysically small
    subcell_size = max(subcell_size,1.0_DP)

    ! Recall that the magnitude of reciprocal lattice vectors is related
    ! to plane spacings in real-space
    spacing(1) = TWO_PI / sqrt(cell%b1 .dot. cell%b1)
    spacing(2) = TWO_PI / sqrt(cell%b2 .dot. cell%b2)
    spacing(3) = TWO_PI / sqrt(cell%b3 .dot. cell%b3)

    ! agreco: check 'spheres' size do not exceed cell size only along the
    ! directions along which NGWFs are localised
    ! initialise variable minspacing
    minspacing = maxval(spacing)

    if (.not.pub_extend_ngwf(1)) then
       minspacing = spacing(1)
    end if
    if (.not.pub_extend_ngwf(2)) then
        minspacing = min(minspacing, spacing(2))
    end if
    if (.not.pub_extend_ngwf(3)) then
        minspacing = min(minspacing, spacing(3))
    end if

    if (.not. present(radius) .and. subcell_size > minspacing) then
      call utils_abort('Error in parallel_strategy_list_overlaps: &
            &spheres exceed cell size.')
    end if

    ! Work out how many of our cells will fit into the unit cell
    ! Do not allow more than 1024 subcells in any direction to avoid
    ! integer overflow
    ! agreco: modify subcells if extended NGWFs are used (28/10/2014)
    do idim=1,3
       num_subcells(idim) = min(max(int(spacing(idim) / subcell_size),1),1024)
       if ((num_subcells(idim) == 2).or.pub_extend_ngwf(idim)) num_subcells(idim) = 1
    end do
    total_num_subcells = num_subcells(1)*num_subcells(2)*num_subcells(3)

    ! Work out which subcell each atom belongs to - each subcell is labelled
    ! by a number from 0 to total_num_subcells-1
    do iat=1,nat
       do idim=1,3
          iscd(idim) = int(num_subcells(idim) * fcoord(idim,iat))
       end do
       atom_subcell(iat) = iscd(1)+num_subcells(1)*(iscd(2)+ &
            num_subcells(2)*iscd(3))
    end do

    ! Count the number of atoms in each subcell and the maximum
    allocate(num_atoms_subcell(0:total_num_subcells-1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell',ierr)
    num_atoms_subcell = 0
    do iat=1,nat
       isc = atom_subcell(iat)
       num_atoms_subcell(isc) = num_atoms_subcell(isc)+1
    end do
    max_atoms_subcell = maxval(num_atoms_subcell)

    ! Make up the cell-list i.e. the list of atoms in each subcell
    allocate(subcell_list(max_atoms_subcell,0:total_num_subcells-1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'subcell_list',ierr)
    num_atoms_subcell = 0
    do iat=1,nat
       isc = atom_subcell(iat)
       num_atoms_subcell(isc) = num_atoms_subcell(isc)+1
       subcell_list(num_atoms_subcell(isc),isc) = iat
    end do

    ! Count the number of overlaps for each atom
    num_my_overlaps = 0
    ! Loop over all subcells i, distributed across procs
    do isc=pub_my_proc_id,total_num_subcells-1,pub_total_num_procs
       iscd(1) = mod(isc,num_subcells(1))    ! Subcell i a1 location
       isc23 = isc/num_subcells(1)
       iscd(2) = mod(isc23,num_subcells(2))  ! Subcell i a2 location
       iscd(3) = isc23/num_subcells(2)       ! Subcell i a3 location
       ! Now loop over this and 26 neighbouring subcells j
       subcellj1: do ksc=0,26
          kscd(1) = mod(ksc,3)-1    ! Relative location of subcell j
          ksc23 = ksc/3             ! wrt subcell i (i.e. -1, 0, 1 in
          kscd(2) = mod(ksc23,3)-1  ! each dimension)
          kscd(3) = (ksc23/3)-1
          ! Calculate absolute location of subcell j
          jscd = iscd+kscd
          ! Apply appropriate boundary conditions
          fextra = 0.0_DP
          do idim=1,3 ! for each dimension
             if (jscd(idim) == -1) then
                ! agreco: no need to check for extra overlaps
                ! in periodic images when NGWFs are extended
                if ((periodic(idim)).and.(.not.pub_extend_ngwf(idim))) then
                   jscd(idim) = num_subcells(idim)-1
                   fextra(idim) = -1.0_DP
                else
                   cycle subcellj1
                end if
             end if
             if (jscd(idim) == num_subcells(idim)) then
                if ((periodic(idim)).and.(.not.pub_extend_ngwf(idim))) then
                   jscd(idim) = 0
                   fextra(idim) = 1.0_DP
                else
                   cycle subcellj1
                end if
             end if
          end do
          ! Calculate label for subcell j
          jsc = jscd(1)+num_subcells(1)*(jscd(2)+num_subcells(2)*jscd(3))
          ! Loop over atoms in subcell i
          do iatsc=1,num_atoms_subcell(isc)
             iat = subcell_list(iatsc,isc)
             if (radii(1,iat) == 0.0_DP) cycle ! radius 0 means nothing there
             ! check for tddft region
             if( loc_tddft_region .and. pub_tddft_kernel_ngroups>1) then
                element_found=.false.
                ! loop over TDDFT kernel groups
                do i_ngroup=1, pub_tddft_kernel_ngroups
                   ! loop over elements in ngroup
                   do ii=1, pub_tddft_kernel_group_nsp(i_ngroup)
                      ! check if atom in question is within that group
                      dummy_id =pub_tddft_kernel_groups(ii,i_ngroup)
                      if(elements(iat)%species_id == dummy_id) then
                         element_found=.true.
                         ! also set tddft_ngroup_index to the group
                         ! it was found in
                         tddft_ngroup_index=i_ngroup
                      endif
                   enddo
                enddo
             endif
             ! Check for tddft CT region:
             if( loc_tddft_region .and. pub_tddft_ct_ngroups>1) then
                ct_element_found=.false.
                ! loop over TDDFT kernel groups
                do i_ngroup=1, pub_tddft_ct_ngroups
                   ! loop over elements in ngroup
                   do ii=1, pub_tddft_ct_group_nsp(i_ngroup)
                      ! check if atom in question is within that group
                      dummy_id =pub_tddft_ct_groups(ii,i_ngroup)
                      if(elements(iat)%species_id == dummy_id) then
                         ct_element_found=.true.
                      endif
                   enddo
                enddo
             endif

             ! if current part is not part of any tddft region, cycle.
             if(loc_tddft_region .and.pub_tddft_kernel_ngroups>1 .and.&
                        .not. element_found) cycle

             ! Loop over atoms in subcell j
             do jatsc=1,num_atoms_subcell(jsc)
                jat = subcell_list(jatsc,jsc)
                if (radii(2,jat) == 0.0_DP) cycle ! radius 0 means nothing there
                ! check for the TDDFT gradient possiblility here:
                if(loc_tddft_region .and. pub_tddft_kernel_ngroups>1) then
                   ! check if jatom is part of the TDDFT group in which iatom
                   ! was found
                   element_found2=.false.
                   do ii=1, pub_tddft_kernel_group_nsp(tddft_ngroup_index)
                      dummy_id=pub_tddft_kernel_groups(ii,tddft_ngroup_index)
                      if(elements(jat)%species_id == dummy_id) then
                         element_found2=.true.
                      endif
                   enddo
                endif
                if(loc_tddft_region .and. pub_tddft_ct_ngroups>1) then
                   ! check if jatom is part of the TDDFT group in which iatom
                   ! was found
                   ct_element_found2=.false.
                   do i_ngroup=1, pub_tddft_ct_ngroups
                      ! loop over elements in ngroup
                      do ii=1, pub_tddft_ct_group_nsp(i_ngroup)
                         dummy_id=pub_tddft_ct_groups(ii,i_ngroup)
                         if(elements(jat)%species_id == dummy_id) then
                            ct_element_found2=.true.
                         endif
                      enddo
                   enddo
                endif

                ! cycle if jatom element is not part of iatom tddft region
                ! but do NOT cycle if they are both in CT groups and closer
                ! than the predefined CT radius.
                if(loc_tddft_region .and.pub_tddft_kernel_ngroups>1 .and.&
                        .not. element_found2) then
                    if(.not. (ct_element_found .and. ct_element_found2)) then
                       cycle
                    else
                       tddft_ct_vec=elements(iat)%centre-elements(jat)%centre
                       if( sqrt(tddft_ct_vec .dot. tddft_ct_vec)> &
                         pub_lr_tddft_ct_length) cycle
                    endif
                endif

                fdiff(:) = fcoord(:,jat)+fextra(:)-fcoord(:,iat)
                ! agreco: if extended NGWFs, only consider contributions
                ! from lattice directions along which NGWFs are localised
                if (pub_extend_ngwf(1)) then
                    fdiff(1) = 0.0_DP
                end if
                if (pub_extend_ngwf(2)) then
                    fdiff(2) = 0.0_DP
                end if
                if (pub_extend_ngwf(3)) then
                    fdiff(3) = 0.0_DP
                end if
                adiff(1) = fdiff(1)*a1(1)+fdiff(2)*a2(1)+fdiff(3)*a3(1)
                adiff(2) = fdiff(1)*a1(2)+fdiff(2)*a2(2)+fdiff(3)*a3(2)
                adiff(3) = fdiff(1)*a1(3)+fdiff(2)*a2(3)+fdiff(3)*a3(3)
                dist_sq = adiff(1)*adiff(1) + adiff(2)*adiff(2) + adiff(3)*adiff(3)
                cutoff = radii(1,iat) + radii(2,jat)
                if (dist_sq <= cutoff*cutoff) then
                   num_my_overlaps(iat) = num_my_overlaps(iat)+1
                end if
             end do ! Loop over atoms in subcell j
          end do ! Loop over atoms in subcell i
       end do subcellj1 ! Loop over subcell j
    end do ! Loop over subcell i

    ! Collect results in from all procs
    max_my_overlaps = maxval(num_my_overlaps)
    call comms_reduce('MAX',max_my_overlaps)

    ! ndmh: Number of overlaps of one atom will never be greater than nat
    max_my_overlaps = min(max_my_overlaps, nat)

    ! Allocate workspace my_overlap_list
    allocate(my_overlap_list(max_my_overlaps,nat),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'my_overlap_list',ierr)

    ! Make up the overlap list for this proc, this time counting multiple
    ! overlaps only once
    num_my_overlaps = 0
    overlapped = .false.
    ! Loop over all subcells i, distributed across procs
    do isc=pub_my_proc_id,total_num_subcells-1,pub_total_num_procs
       iscd(1) = mod(isc,num_subcells(1))    ! Subcell i a1 location
       isc23 = isc/num_subcells(1)
       iscd(2) = mod(isc23,num_subcells(2))  ! Subcell i a2 location
       iscd(3) = isc23/num_subcells(2)       ! Subcell i a3 location
       ! Loop over atoms in subcell i
       do iatsc=1,num_atoms_subcell(isc)
          iat = subcell_list(iatsc,isc)
          if (radii(1,iat) == 0.0_DP) cycle ! radius 0 means nothing there
          ! Check for tddft region in exactly the same way as before
          if( loc_tddft_region .and. pub_tddft_kernel_ngroups>1) then
             element_found=.false.
             ! loop over TDDFT kernel groups
             do i_ngroup=1, pub_tddft_kernel_ngroups
                ! loop over elements in ngroup
                do ii=1, pub_tddft_kernel_group_nsp(i_ngroup)
                   ! check if atom in question is within that group
                   dummy_id =pub_tddft_kernel_groups(ii,i_ngroup)
                   if(elements(iat)%species_id == dummy_id) then
                      element_found=.true.
                      ! also set tddft_ngroup_index to the group
                      ! it was found in
                      tddft_ngroup_index=i_ngroup
                   endif
                enddo
             enddo
          endif
          ! Check for ct TDDFT region
          if( loc_tddft_region .and. pub_tddft_ct_ngroups>1) then
             ct_element_found=.false.
             ! loop over TDDFT kernel groups
             do i_ngroup=1, pub_tddft_ct_ngroups
                ! loop over elements in ngroup
                do ii=1, pub_tddft_ct_group_nsp(i_ngroup)
                   ! check if atom in question is within that group
                   dummy_id =pub_tddft_ct_groups(ii,i_ngroup)
                   if(elements(iat)%species_id == dummy_id) then
                      ct_element_found=.true.
                   endif
                enddo
             enddo
          endif
          ! if current part is not part of any tddft region, cycle.
          if(loc_tddft_region .and.pub_tddft_kernel_ngroups>1 .and.&
                    .not. element_found) cycle

          ! Now loop over this subcell and 26 neighbouring subcells j
          subcellj2: do ksc=0,26
             kscd(1) = mod(ksc,3)-1    ! Relative location of subcell j
             ksc23 = ksc/3             ! wrt subcell i (i.e. -1, 0, 1 in
             kscd(2) = mod(ksc23,3)-1  ! each dimension)
             kscd(3) = (ksc23/3)-1
             ! Calculate absolute location of subcell j
             jscd = iscd+kscd
             ! Apply appropriate boundary conditions
             fextra = 0.0_DP
             do idim=1,3 ! for each dimension
                if (jscd(idim) == -1) then
                   ! agreco: no need to check for extra overlaps
                   ! in periodic images when NGWFs are extended
                   if ((periodic(idim)).and.(.not.pub_extend_ngwf(idim))) then
                      jscd(idim) = num_subcells(idim)-1
                      fextra(idim) = -1.0_DP
                   else
                      cycle subcellj2
                   end if
                end if
                if (jscd(idim) == num_subcells(idim)) then
                   if ((periodic(idim)).and.(.not.pub_extend_ngwf(idim))) then
                      jscd(idim) = 0
                      fextra(idim) = 1.0_DP
                   else
                      cycle subcellj2
                   end if
                end if
             end do
             ! Calculate label for subcell j
             jsc = jscd(1)+num_subcells(1)*(jscd(2)+num_subcells(2)*jscd(3))
             ! Loop over atoms in subcell j
             do jatsc=1,num_atoms_subcell(jsc)
                jat = subcell_list(jatsc,jsc)
                if (radii(2,jat) == 0.0_DP) cycle ! radius 0 means nothing there

                ! check for the TDDFT kernel possiblility here:
                if(loc_tddft_region .and. pub_tddft_kernel_ngroups>1) then
                   ! check if jatom is part of the TDDFT group in which iatom
                   ! was found
                   element_found2=.false.
                   do ii=1, pub_tddft_kernel_group_nsp(tddft_ngroup_index)
                      dummy_id=pub_tddft_kernel_groups(ii,tddft_ngroup_index)
                      if(elements(jat)%species_id == dummy_id) then
                         element_found2=.true.
                      endif
                   enddo
                endif
                ! check for charge_transfer TDDFT region
                if(loc_tddft_region .and. pub_tddft_ct_ngroups>1) then
                   ! check if jatom is part of the TDDFT group in which iatom
                   ! was found
                   ct_element_found2=.false.
                   do i_ngroup=1, pub_tddft_ct_ngroups
                      ! loop over elements in ngroup
                      do ii=1, pub_tddft_ct_group_nsp(i_ngroup)
                         dummy_id=pub_tddft_ct_groups(ii,i_ngroup)
                         if(elements(jat)%species_id == dummy_id) then
                            ct_element_found2=.true.
                         endif
                      enddo
                   enddo
                endif

                ! cycle if jatom element is not part of iatom tddft region
                ! but do NOT cycle if they are both in CT groups and closer
                ! than the predefined CT radius.
                if(loc_tddft_region .and.pub_tddft_kernel_ngroups>1 .and.&
                        .not. element_found2) then
                    if(.not. (ct_element_found .and. ct_element_found2)) then
                       cycle
                    else
                       tddft_ct_vec=elements(iat)%centre-elements(jat)%centre
                       if( sqrt(tddft_ct_vec .dot. tddft_ct_vec)> &
                         pub_lr_tddft_ct_length) cycle
                    endif
                endif

                fdiff(:) = fcoord(:,jat)+fextra(:)-fcoord(:,iat)
                ! agreco: if extended NGWFs, only consider contributions
                ! from lattice directions along which NGWFs are localised
                if (pub_extend_ngwf(1)) then
                    fdiff(1) = 0.0_DP
                end if
                if (pub_extend_ngwf(2)) then
                    fdiff(2) = 0.0_DP
                end if
                if (pub_extend_ngwf(3)) then
                    fdiff(3) = 0.0_DP
                end if
                adiff(1) = fdiff(1)*a1(1)+fdiff(2)*a2(1)+fdiff(3)*a3(1)
                adiff(2) = fdiff(1)*a1(2)+fdiff(2)*a2(2)+fdiff(3)*a3(2)
                adiff(3) = fdiff(1)*a1(3)+fdiff(2)*a2(3)+fdiff(3)*a3(3)
                dist_sq = adiff(1)*adiff(1) + adiff(2)*adiff(2) + adiff(3)*adiff(3)
                cutoff = radii(1,iat) + radii(2,jat)
                if (dist_sq <= cutoff*cutoff) then
                   if (.not. overlapped(jat)) then
                      num_my_overlaps(iat) = num_my_overlaps(iat)+1
                      my_overlap_list(num_my_overlaps(iat),iat) = jat
                      overlapped(jat) = .true.
                   end if
                end if
             end do ! Loop over atoms in subcell j
          end do subcellj2 ! Loop over subcell j
          ! Restore set of overlap flags
          do iovlap=1,num_my_overlaps(iat)
             overlapped(my_overlap_list(iovlap,iat)) = .false.
          end do
       end do ! Loop over atoms in subcell i
    end do ! Loop over subcell i

    ! Allocate workspaces for group overlap_list and remote overlap_list
    allocate(gp_overlap_list(max_my_overlaps,nat),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'gp_overlap_list',ierr)
    allocate(yr_overlap_list(max_my_overlaps,nat),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'yr_overlap_list',ierr)

    ! Collect results in from all procs
    max_my_overlaps = maxval(num_my_overlaps)
    call comms_reduce('MAX',max_my_overlaps)

    overlaps%num_overlaps = num_my_overlaps
    call comms_reduce('SUM',overlaps%num_overlaps)
    overlaps%max_overlaps = maxval(overlaps%num_overlaps)

    ! ebl: DMFT tweak, may need to reintroduce. But this module does not know about par%nat_hub!
    ! overlaps%max_overlaps=max(overlaps%max_overlaps,par%nat_hub)

    ! Allocate array overlaps%overlap_list
    ierr = 0
    if (allocated(overlaps%overlap_list)) then
       if (size(overlaps%overlap_list,1) /= overlaps%max_overlaps .or. &
            size(overlaps%overlap_list,2) /= nat) then
          deallocate(overlaps%overlap_list,stat=ierr)
          call utils_dealloc_check('parallel_strategy_list_overlaps', &
               'overlaps%overlap_list',ierr)
          allocate(overlaps%overlap_list(overlaps%max_overlaps,nat),stat=ierr)
          call utils_alloc_check('parallel_strategy_list_overlaps', &
               'overlaps%overlap_list',ierr)
       end if
    else
       allocate(overlaps%overlap_list(overlaps%max_overlaps,nat),stat=ierr)
       call utils_alloc_check('parallel_strategy_list_overlaps', &
            'overlaps%overlap_list',ierr)
    end if


    ! Make up global overlap list, looping over procs to gather up information
    num_gp_overlaps = 0
    gp_overlap_list = 0
    ! ndmh: first loop over procs in the comms group
    do proc=0,pub_comms_group_size-1
       if (pub_my_rank_in_group == proc) then
          num_yr_overlaps = num_my_overlaps
          yr_overlap_list = my_overlap_list
       end if
       call comms_bcast(proc,num_yr_overlaps,comm=pub_group_comm)
       call comms_barrier(comm=pub_group_comm)
       call comms_bcast(proc,yr_overlap_list,comm=pub_group_comm)
       do iat=1,nat
          nyrat = num_yr_overlaps(iat)
          novat = num_gp_overlaps(iat)
          gp_overlap_list(novat+1:novat+nyrat,iat) = &
               yr_overlap_list(1:nyrat,iat)
          num_gp_overlaps(iat) = novat+nyrat
       end do
    end do

    ! ndmh: now loop over equivalent-ranked procs in other comms groups
    overlaps%num_overlaps = 0
    overlaps%overlap_list = 0
    do group=0,pub_num_comms_groups-1
       if (pub_my_rank_in_rank == group) then
          num_yr_overlaps = num_gp_overlaps
          yr_overlap_list = gp_overlap_list
       end if
       call comms_bcast(group,num_yr_overlaps,comm=pub_rank_comm)
       call comms_barrier(comm=pub_rank_comm)
       call comms_bcast(group,yr_overlap_list,comm=pub_rank_comm)
       do iat=1,nat
          nyrat = num_yr_overlaps(iat)
          novat = overlaps%num_overlaps(iat)
          overlaps%overlap_list(novat+1:novat+nyrat,iat) = &
               yr_overlap_list(1:nyrat,iat)
          overlaps%num_overlaps(iat) = novat+nyrat
       end do
    end do

    ! Verify final result
    if(pub_debug) call internal_verify

    ! Deallocate work arrays
    deallocate(yr_overlap_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'yr_overlap_list',ierr)
    deallocate(my_overlap_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'my_overlap_list',ierr)
    deallocate(gp_overlap_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'gp_overlap_list',ierr)
    deallocate(subcell_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'subcell_list',ierr)
    deallocate(num_atoms_subcell,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell',ierr)

    call internal_deallocate_1

  contains

    !--------------------------------------------------------------------------
    subroutine internal_verify

      !=======================================================================!
      ! This subroutine verifies the overlap list generated by the parent     !
      ! subroutine for consistency                                            !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 16/03/04                                     !
      !=======================================================================!

      implicit none

      ! Local variables
      integer :: iat,jat   ! Atom counters
      integer :: novlaps   ! Number of overlaps
      integer :: iovlap    ! Overlap counter
      integer :: jovlap    ! Overlap counter
      logical :: found     ! Found flag

      ! Loop over all atoms in cell iat
      do iat=1,nat

         ! Check number of overlaps
         novlaps = overlaps%num_overlaps(iat)
         if (novlaps < 1 .and. radii(1,iat) > 0.0_DP .and. &
              radii(2,iat) > 0.0_DP) then
            call utils_abort(&
                 'Error verifying overlap list (parallel_strategy_mod.F90):&
                 & no overlaps detected for atom ',iat)
         end if
         if (novlaps > nat) then
            call utils_abort(&
                 'Error verifying overlap list (parallel_strategy_mod.F90):&
                 & too many overlaps detected for atom ',iat)
         end if

         ! Loop over overlapping atoms jat
         do iovlap=1,novlaps
            jat = overlaps%overlap_list(iovlap,iat)

            ! Check that this atom only appears once in the overlap list
            do jovlap=iovlap+1,novlaps
               if (overlaps%overlap_list(jovlap,iat) == jat) then
                  call utils_abort(&
                       'Error verifying overlap list &
                       &(parallel_strategy_mod.F90): &
                       &repeated overlap for atom pair: ',iat, jat)
               end if
            end do

            ! Check that cross-reference exists (if appropriate)
            if (mode1 == mode2) then

               ! Loop over atoms which overlap atom jat
               found = .false.
               do jovlap=1,overlaps%num_overlaps(jat)
                  if (overlaps%overlap_list(jovlap,jat) == iat) then
                     found = .true.
                     exit
                  end if
               end do
               if (.not. found) then
                  call utils_abort(&
                       'Error verifying overlap list &
                       &(parallel_strategy_mod.F90): &
                       &no cross-reference overlap for atom pair ',iat, jat)
               end if

            end if

         end do  ! jat

      end do  ! iat

    end subroutine internal_verify
    !--------------------------------------------------------------------------

    subroutine internal_check_args(have_radius)

      !=======================================================================!
      ! This subroutine checks the arguments to the parent subroutine         !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! have_radius (input) : flag indicating presence of optional argument   !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      use utils, only: utils_assert, utils_abort

      implicit none

      ! Arguments
      logical, intent(in) :: have_radius

      call utils_assert(nat > 0, &
           'Error in parallel_strategy_list_overlaps: no atoms in cell!')
      call utils_assert(size(elements) == nat, &
           'Error in parallel_strategy_list_overlaps: incompatible nat')

      if (mode1 /= 'C' .and. mode1 /= 'c' .and. mode1 /= 'R' .and. &
           mode1 /= 'r' .and. mode1 /= 'F' .and. mode1 /= 'f' .and. &
           mode1 /= 'A' .and. mode1 /= 'a' .and. mode1 /= 'L' .and. &
           mode1 /= 'l') then
         call utils_abort('Error in parallel_strategy_list_overlaps: &
              &unknown mode1 "'//trim(mode1)//'"')
      end if
      if (mode2 /= 'C' .and. mode2 /= 'c' .and. mode2 /= 'R' .and. &
           mode2 /= 'r' .and. mode2 /= 'F' .and. mode2 /= 'f' .and. &
           mode2 /= 'A' .and. mode2 /= 'a' .and. mode2 /= 'L' .and. &
           mode2 /= 'l') then
         call utils_abort('Error in parallel_strategy_list_overlaps: &
              &unknown mode2 "'//trim(mode2)//'"')
      end if
      if ((mode1 == 'F' .or. mode1 == 'f' .or. mode2 == 'F' .or. &
           mode2 == 'f') .and. (.not. have_radius)) then
         call utils_abort('Error in parallel_strategy_list_overlaps: &
              &radius must be present for fixed mode')
      end if

    end subroutine internal_check_args

    !--------------------------------------------------------------------------

    subroutine internal_allocate_1

      !=======================================================================!
      ! This subroutine (de)allocates arrays as required by the parent        !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      implicit none

      ! Allocate module array _num_overlaps
      ierr = 0
      if (allocated(overlaps%num_overlaps)) then
         if (size(overlaps%num_overlaps) /= nat) then
            deallocate(overlaps%num_overlaps,stat=ierr)
            call utils_dealloc_check('internal_allocate_1 &
                 &(parallel_strategy_list_overlaps)','overlaps%num_overlaps',ierr)
            allocate(overlaps%num_overlaps(nat),stat=ierr)
            call utils_alloc_check('internal_allocate_1 &
                 &(parallel_strategy_list_overlaps)','overlaps%num_overlaps',ierr)
         end if
      else
         allocate(overlaps%num_overlaps(nat),stat=ierr)
         call utils_alloc_check('internal_allocate_1 &
              &(parallel_strategy_list_overlaps)','overlaps%num_overlaps',ierr)
      end if


      ! Allocate work arrays
      allocate(radii(2,nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','radii',ierr)
      allocate(atom_subcell(nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','atom_subcell',ierr)
      allocate(fcoord(3,nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','fcoord',ierr)
      allocate(num_my_overlaps(nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','num_my_overlaps',ierr)
      allocate(num_gp_overlaps(nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','num_gp_overlaps',ierr)
      allocate(num_yr_overlaps(nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','num_yr_overlaps',ierr)
      allocate(overlapped(nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','overlapped',ierr)

    end subroutine internal_allocate_1

    subroutine internal_deallocate_1

      !=======================================================================!
      ! This subroutine deallocates arrays as required by the parent          !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Added by Nicholas Hine, 18/05/09 to ensure dealloc check marks        !
      ! the right arrays as having been deallocated when doing array checking !
      !=======================================================================!

      implicit none

      ! Allocate internal arrays
      deallocate(overlapped,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'overlapped',ierr)
      deallocate(num_yr_overlaps,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'num_yr_overlaps',ierr)
      deallocate(num_gp_overlaps,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'num_gp_overlaps',ierr)
      deallocate(num_my_overlaps,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'num_my_overlaps',ierr)
      deallocate(fcoord,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'fcoord',ierr)
      deallocate(atom_subcell,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'atom_subcell',ierr)
      deallocate(radii,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'radii',ierr)

    end subroutine internal_deallocate_1

  end subroutine parallel_strategy_list_overlaps

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_exit_overlaps(overlaps)

    !=========================================================================!
    ! This subroutine deallocates storage for the OVERLAP_LIST type.          !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 15/5/12                                       !
    !=========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(OVERLAP_LIST), intent(inout) :: overlaps ! List of overlaps

    ! Local Variables
    integer :: ierr

    if (allocated(overlaps%num_overlaps)) then
      deallocate(overlaps%num_overlaps,stat=ierr)
      call utils_dealloc_check('parallel_strategy_exit_overlaps', &
           'overlaps%num_overlaps',ierr)
    end if
    if (allocated(overlaps%overlap_list)) then
      deallocate(overlaps%overlap_list,stat=ierr)
      call utils_dealloc_check('parallel_strategy_exit_overlaps', &
           'overlaps%overlap_list',ierr)
    end if

   end subroutine parallel_strategy_exit_overlaps

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_distr_kpoints(par, kpoints_list)

    !=========================================================================!
    ! This subroutine distributes the k-points across the procs available.    !
    !                                                                         !
    ! Arguments:                                                              !
    !     kpoints_list: the list of kpoints read from the input file          !
    !                                                                         !
    ! Written by Andrea Greco on 03/01/2016.                                  !
    !=========================================================================!

    use comms, only: pub_comms_initialised, pub_on_root, &
        pub_my_proc_id, pub_total_num_procs
    use k_points, only: KPOINT
    use utils, only: utils_assert
    use rundat, only: pub_output_detail

    implicit none

    ! Arguments
    type(PARAL_INFO), intent(inout) :: par
    type(KPOINT), intent(in) :: kpoints_list(:)

    ! Local variables
    integer :: ik       ! k-point loop counter
    integer :: global_kpoint
    real(kind=DP) :: av_kpoints_on_proc
    real(kind=DP) :: excess_kpoints
    integer :: proc
    integer :: proc_ik
    integer :: kpoint_count


    ! agreco: check arguments
    call utils_assert(par%num_kpoints > 0, &
        'Error in parallel_strategy_distr_kpoints: no k-points specified')

    ! agreco: check size of kpoints_list is the same as num_kpoints
    call utils_assert(size(kpoints_list) == par%num_kpoints, &
        'Error in parallel_strategy_distr_kpoints: incompatible num_kpoints')

    ! agreco: check comms module intialised
    call utils_assert(pub_comms_initialised, &
        'Error in parallel_strategy_distr_kpoints: comms module not initialised')

    ! agreco: reallocate arrays if necessary
    call internal_allocate_mod_3

    ! agreco: calculate expected number of k-points per proc
    av_kpoints_on_proc = par%num_kpoints / real(pub_total_num_procs,kind=DP)

    ! agreco: work out number of k-points on each proc - aim is to balance the
    ! number of k-points per proc to balance the work
    proc = 0                    ! Stores current proc being filled up
    par%num_kpoints_on_proc = 0 ! Counts number of k-points on each proc

    ! Put first k-point on proc 0
    par%num_kpoints_on_proc(0) = 1

    ! Loop over remaining k-points
    do ik=2, par%num_kpoints

      ! Calculate what this proc's surplus/deficit of k-points is
      excess_kpoints = par%num_kpoints_on_proc(proc)-av_kpoints_on_proc

      ! If it's a surplus, move on to the next proc (unless we've run out)
      if (abs(excess_kpoints) < abs(excess_kpoints+1) &
         .and. proc < pub_total_num_procs-1) proc = proc+1

      ! Update number of k-points on current proc
      par%num_kpoints_on_proc(proc) = par%num_kpoints_on_proc(proc) + 1

      ! Update the expected average to take history into account
      ! remaining number of k-points / remaining number of procs
      if (proc > 0) av_kpoints_on_proc = &
         (par%num_kpoints-sum(par%num_kpoints_on_proc(0:proc-1))) / &
         real(pub_total_num_procs-proc,kind=DP)

    end do

    ! Find maximum number of k-points on any one proc
    par%max_kpoints_on_proc = maxval(par%num_kpoints_on_proc)

    global_kpoint = 1

    do proc=0,pub_total_num_procs-1
      par%first_kpoint_on_proc(proc) = global_kpoint
      global_kpoint = global_kpoint + par%num_kpoints_on_proc(proc)
    end do

    par%first_kpoint_on_proc(pub_total_num_procs) = par%num_kpoints + 1

    ! Reallocate module arrays which depend upon par%max_kpoints_on_proc
    ! if necessary
    call internal_allocate_mod_4

    proc = 0
    proc_ik = 0

    do ik=1,par%num_kpoints
      ! keep track of number of k-points found on current proc
      proc_ik = proc_ik + 1

      ! associate proc to current k-point
      par%proc_of_kpoint(ik) = proc
      ! list of k-points (indices in original k-points array) for each proc
      par%kp_indices_on_proc(proc_ik,proc) = ik

      ! move to next proc if required
      if (proc_ik==par%num_kpoints_on_proc(proc)) then
        proc = proc + 1
        ! reset counter of k-points on proc
        proc_ik = 0
      end if

    end do

    ! Make up list of k-points on this proc
    ! Loop over k-points on current proc
    do ik=1,par%num_kpoints_on_proc(pub_my_proc_id)
      par%kpoints_on_proc(ik) = &
        kpoints_list(par%kp_indices_on_proc(ik,pub_my_proc_id))
    end do

    ! Set flag to show k-points have been distributed
    kpoints_distributed = .true.

    ! agreco: warning if we have more procs than k-points
    if (par%num_kpoints / pub_total_num_procs < 1) then
      if (pub_on_root.and.pub_output_detail>=VERBOSE) then
        write(stdout,*) ' '
        write(stdout,'(a)') 'WARNING in parallel_strategy_distr_kpoints:'
        write(stdout,'(a)') 'Average of less than one k-point per proc, so &
                            &workload cannot be'
        write(stdout,'(a)') 'equally balanced between procs.'
      end if
    end if

  contains

    subroutine internal_allocate_mod_3

      !=======================================================================!
      ! This subroutine reallocates module arrays as required by the parent   !
      ! subroutine.                                                           !
      !=======================================================================!

      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Local variables
      integer :: ierr     ! Error flag

      if (allocated(par%num_kpoints_on_proc)) then
        if (size(par%num_kpoints_on_proc) /= pub_total_num_procs) then
          deallocate(par%num_kpoints_on_proc,stat=ierr)
          call utils_dealloc_check('internal_allocate_mod_3 &
            &(parallel_strategy_distr_kpoints)','par%num_kpoints_on_proc',ierr)
        end if
      end if
      if (allocated(par%first_kpoint_on_proc)) then
        if (size(par%first_kpoint_on_proc) /= pub_total_num_procs+1) then
          deallocate(par%first_kpoint_on_proc,stat=ierr)
          call utils_dealloc_check('internal_allocate_mod_3 &
            &(parallel_strategy_distr_kpoints)','par%first_kpoint_on_proc',ierr)
        end if
      end if
      if (allocated(par%proc_of_kpoint)) then
        if (size(par%proc_of_kpoint) /= par%num_kpoints) then
          deallocate(par%proc_of_kpoint,stat=ierr)
          call utils_dealloc_check('internal_allocate_mod_3 &
            &(parallel_strategy_distr_kpoints)','par%proc_of_kpoint',ierr)
        end if
      end if
      if (.not.allocated(par%num_kpoints_on_proc)) then
        allocate(par%num_kpoints_on_proc(0:pub_total_num_procs-1),stat=ierr)
        call utils_alloc_check('internal_allocate_mod_3 &
          &(parallel_strategy_distr_kpoints)','par%num_kpoints_on_proc',ierr)
      end if
      if (.not.allocated(par%first_kpoint_on_proc)) then
        allocate(par%first_kpoint_on_proc(0:pub_total_num_procs),stat=ierr)
        call utils_alloc_check('internal_allocate_mod_3 &
          &(parallel_strategy_distr_kpoints)','par%first_kpoint_on_proc',ierr)
      end if
      if (.not.allocated(par%proc_of_kpoint)) then
         allocate(par%proc_of_kpoint(par%num_kpoints),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_3 &
              &(parallel_strategy_distr_kpoints)','par%proc_of_kpoint',ierr)
      end if

    end subroutine internal_allocate_mod_3

    subroutine internal_allocate_mod_4

      !=======================================================================!
      ! This subroutine reallocates module arrays as required by the parent   !
      ! subroutine.                                                           !
      !=======================================================================!

      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Local variables
      integer :: ierr     ! Error flag

      ! Deallocate module arrays which depend upon par%max_kpoints_on_proc
      ! if necessary
      if (allocated(par%kp_indices_on_proc)) then
        if (size(par%kp_indices_on_proc,1) /= par%max_kpoints_on_proc .or. &
             size(par%kp_indices_on_proc,2) /= pub_total_num_procs) then
          deallocate(par%kp_indices_on_proc,stat=ierr)
          call utils_dealloc_check('internal_allocate_mod_4 &
               &(parallel_strategy_distr_kpoints)','par%kp_indices_on_proc',ierr)
        end if
      end if
      if (allocated(par%kpoints_on_proc)) then
        if (size(par%kpoints_on_proc) /= par%max_kpoints_on_proc) then
          deallocate(par%kpoints_on_proc,stat=ierr)
          call utils_dealloc_check('internal_allocate_mod_4 &
               &(parallel_strategy_distr_kpoints)','par%kpoints_on_proc',ierr)
        end if
      end if

      ! Allocate module arrays which depend upon par%max_kpoints_on_proc
      ! if necessary
      if (.not. allocated(par%kp_indices_on_proc)) then
        allocate(par%kp_indices_on_proc(par%max_kpoints_on_proc, &
                 0:pub_total_num_procs-1),stat=ierr)
        call utils_alloc_check('internal_allocate_mod_4 &
             &(parallel_strategy_distr_kpoints)','par%kp_indices_on_proc',ierr)
        ! agreco: initialise to something obvious
        par%kp_indices_on_proc = -1
      end if
      if (.not. allocated(par%kpoints_on_proc)) then
        allocate(par%kpoints_on_proc(par%max_kpoints_on_proc),stat=ierr)
        call utils_alloc_check('internal_allocate_mod_4 &
             &(parallel_strategy_distr_kpoints)','par%kpoints_on_proc',ierr)
      end if

    end subroutine internal_allocate_mod_4

  end subroutine parallel_strategy_distr_kpoints

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_list_cross_overlaps(overlaps,elements1, &
          elements2,cell,mode1,mode2,radius,tddft_region,multiplier)

    !=========================================================================!
    ! This subroutine calculates the list of overlaps between atoms, either   !
    ! according to the region/core radii or according to an optional fixed    !
    ! radius.                                                                 !
    ! The form of overlap list required is determined by the two mode         !
    ! arguments which can be any of the following options:                    !
    !   'C' - use nonlocal pseudopotential Core radii                         !
    !   'R' - use NGWF Region radii                                           !
    !   'A' - use Conduction NGWF Region radii                                !
    !   'F' - use Fixed radius                                                !
    ! mode1 determines the radius used for the primary atom (the index in the !
    ! module arrays) whereas mode2 determines the radius for the secondary    !
    ! atom (listed as overlapping with the primary atom in the list).         !
    ! A cell list algorithm is used to obtain linear scaling for large        !
    ! systems - these cells are referred to as subcells below to avoid        !
    ! confusion with the unit cell.                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   overlaps (inout)  : An OVERLAP_LIST type with details of the overlaps !
    !   elements1 (input) : A list of elements with positions and radii       !
    !   elements2 (input) : A 2nd list of elements with positions and radii   !
    !   cell (input)      : The CELL_INFO type containing the cell parameters !
    !   mode1 (input)     : Mode for primary atom (see above)                 !
    !   mode2 (input)     : Mode for secondary atom (see above)               !
    !   radius (input)    : The optional fixed radius to be used              !
    !   tddft_region(input) : Optional: Used if the TDDFT kernel is limited   !
    !                         To some region in real space                    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   constants       : For pi                                              !
    !   geometry        : For the point type and dot product                  !
    !   simulation_cell : For the CELL_INFO type                              !
    !   ion             : For the element type                                !
    !   utils           : For checking memory allocation/deallocation         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   fcoord   : list of fractional coordinates of the atoms wrt lattice    !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   par%nat > 0                                                           !
    !   size(elements) == par%nat                                             !
    !   mode1 and mode2 must be one of 'C','R' or 'F'                         !
    !   If either mode1 or mode2 == 'F', radius must be present               !
    !   Sphere radii should not be larger than unit cell                      !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/7/03                                        !
    ! Minor modification to reduce size of arrays by Nick Hine 13/03/08       !
    ! Modified by Quintin Hill to use cell on 17/10/2008.                     !
    ! Modified by Nicholas Hine to take cell as argument on 12/08/2014.       !
    ! Modified by Robert Charlton to handle separate element lists, May 2017. !
    !=========================================================================!

    use comms, only: comms_barrier, comms_bcast, comms_reduce, &
        pub_my_proc_id, pub_total_num_procs, pub_my_rank_in_group, &
        pub_comms_group_size, pub_group_comm, pub_rank_comm, &
        pub_my_rank_in_rank, pub_num_comms_groups
    use constants, only : PI, TWO_PI, stdout
    use geometry, only: point,operator(.DOT.),operator(-)
    use ion, only: element
    use rundat, only: pub_tddft_kernel_groups, pub_tddft_kernel_group_nsp,&
        pub_tddft_kernel_ngroups, pub_debug, pub_extend_ngwf,&
        pub_lr_tddft_ct_length,pub_tddft_ct_groups,&
        pub_tddft_ct_group_nsp,pub_tddft_ct_ngroups
    use simulation_cell, only : CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(OVERLAP_LIST), intent(inout) :: overlaps ! List of overlaps
    type(CELL_INFO), intent(in) :: cell
    type(element), intent(in) :: elements1(:)     ! List of atoms
    type(element), intent(in) :: elements2(:)     ! Another list of atoms
    character, intent(in) :: mode1                ! Mode for primary atom
    character, intent(in) :: mode2                ! Mode for secondary atom
    real(kind=DP), intent(in), optional :: radius ! Optional fixed radius
    logical, optional :: tddft_region             ! used if TDDFT kernel is only
                                                  ! defined within a region
    real(kind=dp), intent(in), optional :: multiplier ! ja531-> Used to multiply radii.
                                                      ! ja531-> for FOE sparsity.

    ! Local variables
    logical, parameter :: periodic(3) = (/.true.,.true.,.true./) ! Apply PBC
    integer :: ierr                           ! Error flag
    integer :: iat,jat                        ! Atom loop counters
    integer :: ii                             ! Additional counter
    integer :: idim                           ! Dimension loop counter
    integer :: proc                           ! Proc counter
    integer :: group                          ! Comms group counter
    integer :: isc,iscd(3),isc23              ! Subcell loop counters/labels
    integer :: jsc,jscd(3)                    ! Subcell loop counters/labels
    integer :: ksc,kscd(3),ksc23              ! Subcell loop counters/labels
    integer :: iatsc,jatsc                    ! Atom in subcell loop counters
    integer :: num_subcells(3)                ! Number of subcells in each dim
    integer :: total_num_subcells             ! Total number of subcells
    integer :: max_atoms_subcell              ! Max num atoms in any subcell
    integer :: max_atoms_subcell1             ! Max num atoms in any subcell
    integer :: max_atoms_subcell2             ! Max num atoms in any subcell
    integer :: max_my_overlaps                ! Max num overlaps on any proc
    integer :: novat                          ! Total num overlaps for atom
    integer :: nyrat                          ! Num overlaps for atom on proc
    integer :: iovlap                         ! Overlap loop counter
    real(kind=DP) :: recip_twopi              ! 1 / (2 pi)
    real(kind=DP) :: subcell_size             ! Minimum subcell size
    real(kind=DP) :: spacing(3)               ! Plane spacing of unit cell faces
    real(kind=DP) :: minspacing               ! minimum plane spacing
    real(kind=DP) :: a1(3),a2(3),a3(3)        ! Local copy of lattice vectors
    real(kind=DP) :: fextra(3)                ! Wrap-around for PBC
    real(kind=DP) :: fdiff(3)                 ! Fractional coordinate diffs
    real(kind=DP) :: adiff(3)                 ! Absolute Cartesian diffs
    real(kind=DP) :: dist_sq                  ! Squared distance between atoms
    real(kind=DP) :: cutoff                   ! Cutoff distance
    integer, allocatable :: atom_subcell1(:)  ! Which subcell an atom belongs to
    integer, allocatable :: atom_subcell2(:)  ! Which subcell an atom belongs to
    integer, allocatable :: num_atoms_subcell(:) ! Num atoms in each subcell
    integer, allocatable :: num_atoms_subcell1(:) ! Num atoms in each subcell
    integer, allocatable :: num_atoms_subcell2(:) ! Num atoms in each subcell
    integer, allocatable :: subcell_list(:,:) ! Subcell list of atoms
    integer, allocatable :: subcell_list1(:,:) ! Subcell list of atoms
    integer, allocatable :: subcell_list2(:,:) ! Subcell list of atoms
    integer, allocatable :: num_my_overlaps(:)! Num overlaps on this proc
    integer, allocatable :: num_gp_overlaps(:)! Num overlaps in this group
    integer, allocatable :: num_yr_overlaps(:)! Num overlaps on other proc
    integer, allocatable :: my_overlap_list(:,:) ! Overlaps on this proc
    integer, allocatable :: gp_overlap_list(:,:) ! Overlaps in this group
    integer, allocatable :: yr_overlap_list(:,:) ! Overlaps on other proc
    logical, allocatable :: overlapped(:)      ! Flags to avoid multiple overlaps
    real(kind=DP), allocatable :: fcoord1(:,:) ! Atomic coordinates in fractions
                                               ! of lattice vectors
    real(kind=DP), allocatable :: fcoord2(:,:) ! Atomic coordinates in fractions
                                               ! of lattice vectors
    real(kind=DP), allocatable :: radii(:,:)  ! Radii to be used for primary and
                                              ! secondary atoms
    character(len=4)  :: dummy_id             ! tjz07: variables associated with
    logical  :: element_found                 ! TDDFT sparse region
    logical  :: element_found2
    integer  :: i_ngroup
    logical  :: loc_tddft_region
    integer  :: tddft_ngroup_index
    type(point) :: tddft_ct_vec
    logical :: ct_element_found
    logical :: ct_element_found2
    integer :: nat1, nat2

    real(kind=DP), parameter :: tolerance = 1.0e-10_DP

    ! rc2013: to avoid bugs below, set nat here
    nat1 = size(elements1)
    nat2 = size(elements2)
    ! set default parameter
    loc_tddft_region=.false.
    if(present(tddft_region)) then
       loc_tddft_region=tddft_region
    endif

    ! Check arguments
    call internal_check_args(present(radius))

    ! Allocate module and work arrays
    call internal_allocate_1

    ! Make up lists of radii for primary and secondary atoms according to
    ! the desired modes
    select case (mode1)
    case('C','c')
        do iat=1,nat1
            radii(1,iat) = elements1(iat)%max_core_radius
        end do
    case('R','r')
        do iat=1,nat1
            radii(1,iat) = elements1(iat)%radius
        end do
    case('A','a')
        do iat=1,nat1
            radii(1,iat) = elements1(iat)%radius_cond
        end do
    case('F','f')
        radii(1,:) = radius
    case('L','l')
        do iat=1,nat1
            radii(1,iat) = elements1(iat)%max_core_wf_radius
        end do
    end select
    ! rc2013: now get the data from the 2nd list
    select case (mode2)
    case('C','c')
        do iat=1,nat2
            radii(2,iat) = elements2(iat)%max_core_radius
        end do
    case('R','r')
        do iat=1,nat2
            radii(2,iat) = elements2(iat)%radius
        end do
    case('A','a')
        do iat=1,nat2
            radii(2,iat) = elements2(iat)%radius_cond
        end do
    case('F','f')
        radii(2,:) = radius
    case('L','l')
        do iat=1,nat2
            radii(2,iat) = elements2(iat)%max_core_wf_radius
        end do
    end select

    ! rc2013: zero the extra elements of radii -- should probably split this up
    if (nat1 > nat2) then
        do iat=nat2+1,nat1
            radii(2,iat) = 0.0_DP
        end do
    else if (nat1 < nat2) then
        do iat=nat1+1,nat2
            radii(1,iat) = 0.0_DP
        end do
    end if

    !ja531-> if requested, then multiply the radii by multiplier.
    if(present(multiplier)) then
       radii=radii*multiplier
    end if

    ! tjz07: tddft kernel only has overlaps in a predefined region
    ! THIS is only true if only one group of atoms is defined as
    ! the region, ie. if pub_tddft_kernel_ngroups==1. If multiple
    ! atom groups are defined, then there are multiple regions in
    ! the system. In this case, tddft kernel is zero for any
    ! matrix elements corresponding to atoms between regions, but
    ! is not zero for matrix elements within a single region.
    if( loc_tddft_region .and. pub_tddft_kernel_ngroups==1) then
       ! loop over atoms:
       do iat=1, nat1
          element_found=.false.
          do ii=1, pub_tddft_kernel_group_nsp(1)
             dummy_id =pub_tddft_kernel_groups(ii,1)
             if(elements1(iat)%species_id == dummy_id) then
                element_found=.true.
             endif
          enddo
          if(.not. element_found) then
             radii(1,iat) =0.0_DP
          endif
       enddo
       ! rc2013: repeat for 2nd elements list
       do iat=1, nat2
          element_found=.false.
          do ii=1, pub_tddft_kernel_group_nsp(1)
             dummy_id =pub_tddft_kernel_groups(ii,1)
             if(elements2(iat)%species_id == dummy_id) then
                element_found=.true.
             endif
          enddo
          if(.not. element_found) then
             radii(2,iat) =0.0_DP
          endif
       enddo
    endif


    ! Convert Cartesian coordinates from cell into fractional coordinates in
    ! the interval [0,1)
    do iat=1,nat1
       fcoord1(1,iat) = elements1(iat)%centre .dot. cell%b1
       fcoord1(2,iat) = elements1(iat)%centre .dot. cell%b2
       fcoord1(3,iat) = elements1(iat)%centre .dot. cell%b3
    end do
    do iat=1,nat2
       fcoord2(1,iat) = elements2(iat)%centre .dot. cell%b1
       fcoord2(2,iat) = elements2(iat)%centre .dot. cell%b2
       fcoord2(3,iat) = elements2(iat)%centre .dot. cell%b3
    end do

    ! jd: This ensures fcoord[12](:,:) are *really* within [0,1). Due to floating-precision
    !     rounding errors fcoord[12](:,:) can wind up very slightly negative on entry here.
    !     The subsequent modulo(,1_DP) operation then makes them exactly 1.0_DP, which
    !     leads to an off-by-one calc of subcell number.
    do iat=1,nat1
       do jat = 1,3
           if (fcoord1(jat,iat) < 0.0_DP .and. fcoord1(jat,iat) > (-tolerance)) then
              fcoord1(jat,iat) = 0.0_DP
           elseif (fcoord1(jat,iat) > 1.0_DP .and. fcoord1(jat,iat) < 1.0_DP+tolerance) then
              fcoord1(jat,iat) = 1.0_DP
           end if
       end do
    end do
    do iat=1,nat2
       do jat = 1,3
           if (fcoord2(jat,iat) < 0.0_DP .and. fcoord2(jat,iat) > (-tolerance)) then
              fcoord2(jat,iat) = 0.0_DP
           elseif (fcoord2(jat,iat) > 1.0_DP .and. fcoord2(jat,iat) < 1.0_DP+tolerance) then
              fcoord2(jat,iat) = 1.0_DP
           end if
       end do
    end do

    recip_twopi = 0.5_DP / pi
    fcoord1 = fcoord1 * recip_twopi
    fcoord1 = modulo(fcoord1,1.0_DP)
    fcoord2 = fcoord2 * recip_twopi
    fcoord2 = modulo(fcoord2,1.0_DP)

    ! Make local copies of lattice vectors
    a1(1) = cell%a1%x ; a1(2) = cell%a1%y ; a1(3) = cell%a1%z
    a2(1) = cell%a2%x ; a2(2) = cell%a2%y ; a2(3) = cell%a2%z
    a3(1) = cell%a3%x ; a3(2) = cell%a3%y ; a3(3) = cell%a3%z

    ! Figure out the maximum potential cutoff separation - this defines the
    ! minimum subcell size required
    subcell_size = maxval(radii(1,:)) + maxval(radii(2,:))

    ! Avoid subcell sizes which are unphysically small
    subcell_size = max(subcell_size,1.0_DP)

    ! Recall that the magnitude of reciprocal lattice vectors is related
    ! to plane spacings in real-space
    spacing(1) = TWO_PI / sqrt(cell%b1 .dot. cell%b1)
    spacing(2) = TWO_PI / sqrt(cell%b2 .dot. cell%b2)
    spacing(3) = TWO_PI / sqrt(cell%b3 .dot. cell%b3)

    ! agreco: check 'spheres' size do not exceed cell size only along the
    ! directions along which NGWFs are localised
    ! initialise variable minspacing
    minspacing = maxval(spacing)

    if (.not.pub_extend_ngwf(1)) then
       minspacing = spacing(1)
    end if
    if (.not.pub_extend_ngwf(2)) then
        minspacing = min(minspacing, spacing(2))
    end if
    if (.not.pub_extend_ngwf(3)) then
        minspacing = min(minspacing, spacing(3))
    end if

    if (.not. present(radius) .and. subcell_size > minspacing) then
      call utils_abort('Error in parallel_strategy_list_cross_overlaps: &
            &spheres exceed cell size.')
    end if

    ! Work out how many of our cells will fit into the unit cell
    ! Do not allow more than 1024 subcells in any direction to avoid
    ! integer overflow
    ! agreco: modify subcells if extended NGWFs are used (28/10/2014)
    do idim=1,3
       num_subcells(idim) = min(max(int(spacing(idim) / subcell_size),1),1024)
       if ((num_subcells(idim) == 2).or.pub_extend_ngwf(idim)) num_subcells(idim) = 1
    end do
    total_num_subcells = num_subcells(1)*num_subcells(2)*num_subcells(3)

    ! Work out which subcell each atom belongs to - each subcell is labelled
    ! by a number from 0 to total_num_subcells-1
    ! rc2013: total will include atoms from both lists
    do iat=1,nat1
       do idim=1,3
          iscd(idim) = int(num_subcells(idim) * fcoord1(idim,iat))
       end do
       atom_subcell1(iat) = iscd(1)+num_subcells(1)*(iscd(2)+ &
            num_subcells(2)*iscd(3))
    end do
    ! rc2013: repeat for 2nd system
    do iat=1,nat2
       do idim=1,3
          iscd(idim) = int(num_subcells(idim) * fcoord2(idim,iat))
       end do
       atom_subcell2(iat) = iscd(1)+num_subcells(1)*(iscd(2)+ &
            num_subcells(2)*iscd(3))
    end do

    ! Count the number of atoms in each subcell and the maximum
    allocate(num_atoms_subcell(0:total_num_subcells-1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell',ierr)
    allocate(num_atoms_subcell1(0:total_num_subcells-1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell1',ierr)
    allocate(num_atoms_subcell2(0:total_num_subcells-1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell2',ierr)

    ! rc2013: initialise all these to 0
    num_atoms_subcell = 0
    num_atoms_subcell1 = 0
    num_atoms_subcell2 = 0
    do iat=1,nat1
       isc = atom_subcell1(iat)
       num_atoms_subcell1(isc) = num_atoms_subcell1(isc)+1
       num_atoms_subcell(isc) = num_atoms_subcell(isc)+1
    end do
    ! rc2013: now add the number from the 2nd system
    do iat=1,nat2
       isc = atom_subcell2(iat)
       ! rc2013: need to count both totals and each block for counting...
       num_atoms_subcell(isc) = num_atoms_subcell(isc)+1
       num_atoms_subcell2(isc) = num_atoms_subcell2(isc)+1
    end do
    max_atoms_subcell  = maxval(num_atoms_subcell)
    max_atoms_subcell1 = maxval(num_atoms_subcell1)
    max_atoms_subcell2 = maxval(num_atoms_subcell2)

    ! Make up the cell-lists i.e. the list of atoms in each subcell
    allocate(subcell_list1(max_atoms_subcell,0:total_num_subcells-1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'subcell_list1',ierr)
    allocate(subcell_list2(max_atoms_subcell,0:total_num_subcells-1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'subcell_list2',ierr)
    num_atoms_subcell = 0
    num_atoms_subcell1 = 0
    num_atoms_subcell2 = 0
    do iat=1,nat1
       isc = atom_subcell1(iat)
       num_atoms_subcell1(isc) = num_atoms_subcell1(isc)+1
       num_atoms_subcell(isc) = num_atoms_subcell(isc)+1
       subcell_list1(num_atoms_subcell1(isc),isc) = iat
    end do
    ! rc2013: now add the number from the 2nd system
    do iat=1,nat2
       isc = atom_subcell2(iat)
       ! rc2013: need to count both totals and each block for counting...
       num_atoms_subcell(isc) = num_atoms_subcell(isc)+1
       num_atoms_subcell2(isc) = num_atoms_subcell2(isc)+1
       subcell_list2(num_atoms_subcell2(isc),isc) = iat
    end do

    ! Count the number of overlaps for each atom
    num_my_overlaps = 0
    ! Loop over all subcells i, distributed across procs
    do isc=pub_my_proc_id,total_num_subcells-1,pub_total_num_procs
       iscd(1) = mod(isc,num_subcells(1))    ! Subcell i a1 location
       isc23 = isc/num_subcells(1)
       iscd(2) = mod(isc23,num_subcells(2))  ! Subcell i a2 location
       iscd(3) = isc23/num_subcells(2)       ! Subcell i a3 location
       ! Now loop over this and 26 neighbouring subcells j
       subcellj1: do ksc=0,26
          kscd(1) = mod(ksc,3)-1    ! Relative location of subcell j
          ksc23 = ksc/3             ! wrt subcell i (i.e. -1, 0, 1 in
          kscd(2) = mod(ksc23,3)-1  ! each dimension)
          kscd(3) = (ksc23/3)-1
          ! Calculate absolute location of subcell j
          jscd = iscd+kscd
          ! Apply appropriate boundary conditions
          fextra = 0.0_DP
          do idim=1,3 ! for each dimension
             if (jscd(idim) == -1) then
                ! agreco: no need to check for extra overlaps
                ! in periodic images when NGWFs are extended
                if ((periodic(idim)).and.(.not.pub_extend_ngwf(idim))) then
                   jscd(idim) = num_subcells(idim)-1
                   fextra(idim) = -1.0_DP
                else
                   cycle subcellj1
                end if
             end if
             if (jscd(idim) == num_subcells(idim)) then
                if ((periodic(idim)).and.(.not.pub_extend_ngwf(idim))) then
                   jscd(idim) = 0
                   fextra(idim) = 1.0_DP
                else
                   cycle subcellj1
                end if
             end if
          end do
          ! Calculate label for subcell j
          jsc = jscd(1)+num_subcells(1)*(jscd(2)+num_subcells(2)*jscd(3))
          ! Loop over atoms in subcell i
          ! rc2013: and in region 1!
          do iatsc=1,num_atoms_subcell1(isc)
             iat = subcell_list1(iatsc,isc)
             if (radii(1,iat) == 0.0_DP) cycle ! radius 0 means nothing there
             ! check for tddft region
             if( loc_tddft_region .and. pub_tddft_kernel_ngroups>1) then
                element_found=.false.
                ! loop over TDDFT kernel groups
                do i_ngroup=1, pub_tddft_kernel_ngroups
                   ! loop over elements in ngroup
                   do ii=1, pub_tddft_kernel_group_nsp(i_ngroup)
                      ! check if atom in question is within that group
                      dummy_id =pub_tddft_kernel_groups(ii,i_ngroup)
                      if(elements1(iat)%species_id == dummy_id) then
                         element_found=.true.
                         ! also set tddft_ngroup_index to the group
                         ! it was found in
                         tddft_ngroup_index=i_ngroup
                      endif
                   enddo
                enddo
             endif
             ! Check for tddft CT region:
             if( loc_tddft_region .and. pub_tddft_ct_ngroups>1) then
                ct_element_found=.false.
                ! loop over TDDFT kernel groups
                do i_ngroup=1, pub_tddft_ct_ngroups
                   ! loop over elements in ngroup
                   do ii=1, pub_tddft_ct_group_nsp(i_ngroup)
                      ! check if atom in question is within that group
                      dummy_id =pub_tddft_ct_groups(ii,i_ngroup)
                      if(elements1(iat)%species_id == dummy_id) then
                         ct_element_found=.true.
                      endif
                   enddo
                enddo
             endif

             ! if current part is not part of any tddft region, cycle.
             if(loc_tddft_region .and.pub_tddft_kernel_ngroups>1 .and.&
                        .not. element_found) cycle

             ! Loop over atoms in subcell j
             do jatsc=1,num_atoms_subcell2(jsc)
                jat = subcell_list2(jatsc,jsc)
                if (radii(2,jat) == 0.0_DP) cycle ! radius 0 means nothing there
                ! check for the TDDFT gradient possiblility here:
                if(loc_tddft_region .and. pub_tddft_kernel_ngroups>1) then
                   ! check if jatom is part of the TDDFT group in which iatom
                   ! was found
                   element_found2=.false.
                   do ii=1, pub_tddft_kernel_group_nsp(tddft_ngroup_index)
                      dummy_id=pub_tddft_kernel_groups(ii,tddft_ngroup_index)
                      if(elements2(jat)%species_id == dummy_id) then
                         element_found2=.true.
                      endif
                   enddo
                endif
                if(loc_tddft_region .and. pub_tddft_ct_ngroups>1) then
                   ! check if jatom is part of the TDDFT group in which iatom
                   ! was found
                   ct_element_found2=.false.
                   do i_ngroup=1, pub_tddft_ct_ngroups
                      ! loop over elements in ngroup
                      do ii=1, pub_tddft_ct_group_nsp(i_ngroup)
                         dummy_id=pub_tddft_ct_groups(ii,i_ngroup)
                         if(elements2(jat)%species_id == dummy_id) then
                            ct_element_found2=.true.
                         endif
                      enddo
                   enddo
                endif

                ! cycle if jatom element is not part of iatom tddft region
                ! but do NOT cycle if they are both in CT groups and closer
                ! than the predefined CT radius.
                if(loc_tddft_region .and.pub_tddft_kernel_ngroups>1 .and.&
                        .not. element_found2) then
                    if(.not. (ct_element_found .and. ct_element_found2)) then
                       cycle
                    else
                       tddft_ct_vec=elements1(iat)%centre-elements2(jat)%centre
                       if( sqrt(tddft_ct_vec .dot. tddft_ct_vec)> &
                         pub_lr_tddft_ct_length) cycle
                    endif
                endif

                fdiff(:) = fcoord2(:,jat)+fextra(:)-fcoord1(:,iat)
                ! agreco: if extended NGWFs, only consider contributions
                ! from lattice directions along which NGWFs are localised
                if (pub_extend_ngwf(1)) then
                    fdiff(1) = 0.0_DP
                end if
                if (pub_extend_ngwf(2)) then
                    fdiff(2) = 0.0_DP
                end if
                if (pub_extend_ngwf(3)) then
                    fdiff(3) = 0.0_DP
                end if
                adiff(1) = fdiff(1)*a1(1)+fdiff(2)*a2(1)+fdiff(3)*a3(1)
                adiff(2) = fdiff(1)*a1(2)+fdiff(2)*a2(2)+fdiff(3)*a3(2)
                adiff(3) = fdiff(1)*a1(3)+fdiff(2)*a2(3)+fdiff(3)*a3(3)
                dist_sq = adiff(1)*adiff(1) + adiff(2)*adiff(2) + adiff(3)*adiff(3)
                cutoff = radii(1,iat) + radii(2,jat)
                if (dist_sq <= cutoff*cutoff) then
                   num_my_overlaps(iat) = num_my_overlaps(iat)+1
                end if
             end do ! Loop over atoms in subcell j
          end do ! Loop over atoms in subcell i
       end do subcellj1 ! Loop over subcell j
    end do ! Loop over subcell i

    ! Collect results in from all procs
    max_my_overlaps = maxval(num_my_overlaps)
    call comms_reduce('MAX',max_my_overlaps)

    ! ndmh: Number of overlaps of one atom will never be greater than nat2
    max_my_overlaps = min(max_my_overlaps, nat2)

    ! Allocate workspace my_overlap_list
    allocate(my_overlap_list(max_my_overlaps,nat1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'my_overlap_list',ierr)

    ! Make up the overlap list for this proc, this time counting multiple
    ! overlaps only once
    num_my_overlaps = 0
    overlapped = .false.
    ! Loop over all subcells i, distributed across procs
    do isc=pub_my_proc_id,total_num_subcells-1,pub_total_num_procs
       iscd(1) = mod(isc,num_subcells(1))    ! Subcell i a1 location
       isc23 = isc/num_subcells(1)
       iscd(2) = mod(isc23,num_subcells(2))  ! Subcell i a2 location
       iscd(3) = isc23/num_subcells(2)       ! Subcell i a3 location
       ! Loop over atoms in subcell i
       do iatsc=1,num_atoms_subcell1(isc)
          iat = subcell_list1(iatsc,isc)
          if (radii(1,iat) == 0.0_DP) cycle ! radius 0 means nothing there
          ! Check for tddft region in exactly the same way as before
          if( loc_tddft_region .and. pub_tddft_kernel_ngroups>1) then
             element_found=.false.
             ! loop over TDDFT kernel groups
             do i_ngroup=1, pub_tddft_kernel_ngroups
                ! loop over elements in ngroup
                do ii=1, pub_tddft_kernel_group_nsp(i_ngroup)
                   ! check if atom in question is within that group
                   dummy_id =pub_tddft_kernel_groups(ii,i_ngroup)
                   if(elements1(iat)%species_id == dummy_id) then
                      element_found=.true.
                      ! also set tddft_ngroup_index to the group
                      ! it was found in
                      tddft_ngroup_index=i_ngroup
                   endif
                enddo
             enddo
          endif
          ! Check for ct TDDFT region
          if( loc_tddft_region .and. pub_tddft_ct_ngroups>1) then
             ct_element_found=.false.
             ! loop over TDDFT kernel groups
             do i_ngroup=1, pub_tddft_ct_ngroups
                ! loop over elements in ngroup
                do ii=1, pub_tddft_ct_group_nsp(i_ngroup)
                   ! check if atom in question is within that group
                   dummy_id =pub_tddft_ct_groups(ii,i_ngroup)
                   if(elements1(iat)%species_id == dummy_id) then
                      ct_element_found=.true.
                   endif
                enddo
             enddo
          endif
          ! if current part is not part of any tddft region, cycle.
          if(loc_tddft_region .and.pub_tddft_kernel_ngroups>1 .and.&
                    .not. element_found) cycle

          ! Now loop over this subcell and 26 neighbouring subcells j
          subcellj2: do ksc=0,26
             kscd(1) = mod(ksc,3)-1    ! Relative location of subcell j
             ksc23 = ksc/3             ! wrt subcell i (i.e. -1, 0, 1 in
             kscd(2) = mod(ksc23,3)-1  ! each dimension)
             kscd(3) = (ksc23/3)-1
             ! Calculate absolute location of subcell j
             jscd = iscd+kscd
             ! Apply appropriate boundary conditions
             fextra = 0.0_DP
             do idim=1,3 ! for each dimension
                if (jscd(idim) == -1) then
                   ! agreco: no need to check for extra overlaps
                   ! in periodic images when NGWFs are extended
                   if ((periodic(idim)).and.(.not.pub_extend_ngwf(idim))) then
                      jscd(idim) = num_subcells(idim)-1
                      fextra(idim) = -1.0_DP
                   else
                      cycle subcellj2
                   end if
                end if
                if (jscd(idim) == num_subcells(idim)) then
                   if ((periodic(idim)).and.(.not.pub_extend_ngwf(idim))) then
                      jscd(idim) = 0
                      fextra(idim) = 1.0_DP
                   else
                      cycle subcellj2
                   end if
                end if
             end do
             ! Calculate label for subcell j
             jsc = jscd(1)+num_subcells(1)*(jscd(2)+num_subcells(2)*jscd(3))
             ! Loop over atoms in subcell j
             do jatsc=1,num_atoms_subcell2(jsc)
                jat = subcell_list2(jatsc,jsc)
                if (radii(2,jat) == 0.0_DP) cycle ! radius 0 means nothing there

                ! check for the TDDFT kernel possiblility here:
                if(loc_tddft_region .and. pub_tddft_kernel_ngroups>1) then
                   ! check if jatom is part of the TDDFT group in which iatom
                   ! was found
                   element_found2=.false.
                   do ii=1, pub_tddft_kernel_group_nsp(tddft_ngroup_index)
                      dummy_id=pub_tddft_kernel_groups(ii,tddft_ngroup_index)
                      if(elements2(jat)%species_id == dummy_id) then
                         element_found2=.true.
                      endif
                   enddo
                endif
                ! check for charge_transfer TDDFT region
                if(loc_tddft_region .and. pub_tddft_ct_ngroups>1) then
                   ! check if jatom is part of the TDDFT group in which iatom
                   ! was found
                   ct_element_found2=.false.
                   do i_ngroup=1, pub_tddft_ct_ngroups
                      ! loop over elements in ngroup
                      do ii=1, pub_tddft_ct_group_nsp(i_ngroup)
                         dummy_id=pub_tddft_ct_groups(ii,i_ngroup)
                         if(elements2(jat)%species_id == dummy_id) then
                            ct_element_found2=.true.
                         endif
                      enddo
                   enddo
                endif

                ! cycle if jatom element is not part of iatom tddft region
                ! but do NOT cycle if they are both in CT groups and closer
                ! than the predefined CT radius.
                if(loc_tddft_region .and.pub_tddft_kernel_ngroups>1 .and.&
                        .not. element_found2) then
                    if(.not. (ct_element_found .and. ct_element_found2)) then
                       cycle
                    else
                       tddft_ct_vec=elements1(iat)%centre-elements2(jat)%centre
                       if( sqrt(tddft_ct_vec .dot. tddft_ct_vec)> &
                         pub_lr_tddft_ct_length) cycle
                    endif
                endif


                fdiff(:) = fcoord2(:,jat)+fextra(:)-fcoord1(:,iat)
                ! agreco: if extended NGWFs, only consider contributions
                ! from lattice directions along which NGWFs are localised
                if (pub_extend_ngwf(1)) then
                    fdiff(1) = 0.0_DP
                end if
                if (pub_extend_ngwf(2)) then
                    fdiff(2) = 0.0_DP
                end if
                if (pub_extend_ngwf(3)) then
                    fdiff(3) = 0.0_DP
                end if
                adiff(1) = fdiff(1)*a1(1)+fdiff(2)*a2(1)+fdiff(3)*a3(1)
                adiff(2) = fdiff(1)*a1(2)+fdiff(2)*a2(2)+fdiff(3)*a3(2)
                adiff(3) = fdiff(1)*a1(3)+fdiff(2)*a2(3)+fdiff(3)*a3(3)
                dist_sq = adiff(1)*adiff(1) + adiff(2)*adiff(2) + adiff(3)*adiff(3)
                cutoff = radii(1,iat) + radii(2,jat)
                if (dist_sq <= cutoff*cutoff) then
                   if (.not. overlapped(jat)) then
                      num_my_overlaps(iat) = num_my_overlaps(iat)+1
                      my_overlap_list(num_my_overlaps(iat),iat) = jat
                      overlapped(jat) = .true.
                   end if
                end if
             end do ! Loop over atoms in subcell j
          end do subcellj2 ! Loop over subcell j
          ! Restore set of overlap flags
          do iovlap=1,num_my_overlaps(iat)
             overlapped(my_overlap_list(iovlap,iat)) = .false.
          end do
       end do ! Loop over atoms in subcell i
    end do ! Loop over subcell i

    ! Allocate workspaces for group overlap_list and remote overlap_list
    allocate(gp_overlap_list(max_my_overlaps,nat1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'gp_overlap_list',ierr)
    allocate(yr_overlap_list(max_my_overlaps,nat1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'yr_overlap_list',ierr)

    ! Collect results in from all procs
    max_my_overlaps = maxval(num_my_overlaps)
    call comms_reduce('MAX',max_my_overlaps)

    overlaps%num_overlaps = num_my_overlaps
    call comms_reduce('SUM',overlaps%num_overlaps)
    overlaps%max_overlaps = maxval(overlaps%num_overlaps)

    ! ebl: DMFT tweak, may need to reintroduce. But this module does not know about par%nat_hub!
    ! overlaps%max_overlaps=max(overlaps%max_overlaps,par%nat_hub)

    ! Allocate array overlaps%overlap_list
    ierr = 0
    if (allocated(overlaps%overlap_list)) then
       if (size(overlaps%overlap_list,1) /= overlaps%max_overlaps .or. &
            size(overlaps%overlap_list,2) /= nat1) then
          deallocate(overlaps%overlap_list,stat=ierr)
          call utils_dealloc_check('parallel_strategy_list_overlaps', &
               'overlaps%overlap_list',ierr)
          allocate(overlaps%overlap_list(overlaps%max_overlaps,nat1),stat=ierr)
          call utils_alloc_check('parallel_strategy_list_overlaps', &
               'overlaps%overlap_list',ierr)
       end if
    else
       allocate(overlaps%overlap_list(overlaps%max_overlaps,nat1),stat=ierr)
       call utils_alloc_check('parallel_strategy_list_overlaps', &
            'overlaps%overlap_list',ierr)
    end if


    ! Make up global overlap list, looping over procs to gather up information
    num_gp_overlaps = 0
    gp_overlap_list = 0
    ! ndmh: first loop over procs in the comms group
    do proc=0,pub_comms_group_size-1
       if (pub_my_rank_in_group == proc) then
          num_yr_overlaps = num_my_overlaps
          yr_overlap_list = my_overlap_list
       end if
       call comms_bcast(proc,num_yr_overlaps,comm=pub_group_comm)
       call comms_barrier(comm=pub_group_comm)
       call comms_bcast(proc,yr_overlap_list,comm=pub_group_comm)
       do iat=1,nat1
          nyrat = num_yr_overlaps(iat)
          novat = num_gp_overlaps(iat)
          gp_overlap_list(novat+1:novat+nyrat,iat) = &
               yr_overlap_list(1:nyrat,iat)
          num_gp_overlaps(iat) = novat+nyrat
       end do
    end do

    ! ndmh: now loop over equivalent-ranked procs in other comms groups
    overlaps%num_overlaps = 0
    overlaps%overlap_list = 0
    do group=0,pub_num_comms_groups-1
       if (pub_my_rank_in_rank == group) then
          num_yr_overlaps = num_gp_overlaps
          yr_overlap_list = gp_overlap_list
       end if
       call comms_bcast(group,num_yr_overlaps,comm=pub_rank_comm)
       call comms_barrier(comm=pub_rank_comm)
       call comms_bcast(group,yr_overlap_list,comm=pub_rank_comm)
       do iat=1,nat1
          nyrat = num_yr_overlaps(iat)
          novat = overlaps%num_overlaps(iat)
          overlaps%overlap_list(novat+1:novat+nyrat,iat) = &
               yr_overlap_list(1:nyrat,iat)
          overlaps%num_overlaps(iat) = novat+nyrat
       end do
    end do

    ! Verify final result
    ! rc2013: internal_verify switched off for now due to us not having
    ! the structures needed to run it for every atom
    !if(pub_debug) call internal_verify

    ! Deallocate work arrays
    deallocate(yr_overlap_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'yr_overlap_list',ierr)
    deallocate(my_overlap_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'my_overlap_list',ierr)
    deallocate(gp_overlap_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'gp_overlap_list',ierr)
    deallocate(subcell_list1,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'subcell_list1',ierr)
    deallocate(subcell_list2,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'subcell_list2',ierr)
    deallocate(num_atoms_subcell,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell',ierr)
    deallocate(num_atoms_subcell1,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell1',ierr)
    deallocate(num_atoms_subcell2,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell2',ierr)

    call internal_deallocate_1

  contains

    !--------------------------------------------------------------------------
    subroutine internal_verify

      !=======================================================================!
      ! This subroutine verifies the overlap list generated by the parent     !
      ! subroutine for consistency                                            !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 16/03/04                                     !
      ! This routine is quite problematic, as it takes the original overlaps  !
      ! list.
      !=======================================================================!

      implicit none

      ! Local variables
      integer :: iat,jat   ! Atom counters
      integer :: novlaps   ! Number of overlaps
      integer :: iovlap    ! Overlap counter
      integer :: jovlap    ! Overlap counter
      logical :: found     ! Found flag

      ! Loop over all atoms in cell iat
      do iat=1,nat1

         ! Check number of overlaps
         novlaps = overlaps%num_overlaps(iat)
         do jat=1,nat2
            if (novlaps < 1 .and. radii(1,iat) > 0.0_DP .and. &
                 radii(2,jat) > 0.0_DP) then
               call utils_abort(&
                   'Error verifying overlap list (parallel_strategy_mod.F90):&
                   & no overlaps detected for atom ',iat)
            end if
         end do
         if (novlaps > (nat1+nat2)) then
            call utils_abort(&
                 'Error verifying overlap list (parallel_strategy_mod.F90):&
                 & too many overlaps detected for atom ',iat)
         end if

         ! Loop over overlapping atoms jat
         do iovlap=1,novlaps
            jat = overlaps%overlap_list(iovlap,iat)

            ! Check that this atom only appears once in the overlap list
            do jovlap=iovlap+1,novlaps
               if (overlaps%overlap_list(jovlap,iat) == jat) then
                  call utils_abort(&
                       'Error verifying overlap list &
                       &(parallel_strategy_mod.F90): &
                       &repeated overlap for atom pair: ',iat, jat)
               end if
            end do

            ! Check that cross-reference exists (if appropriate)
            if (mode1 == mode2) then

               ! Loop over atoms which overlap atom jat
               found = .false.
               do jovlap=1,overlaps%num_overlaps(jat)
                  if (overlaps%overlap_list(jovlap,jat) == iat) then
                     found = .true.
                     exit
                  end if
               end do
               if (.not. found) then
                  call utils_abort(&
                       'Error verifying overlap list &
                       &(parallel_strategy_mod.F90): &
                       &no cross-reference overlap for atom pair ',iat, jat)
               end if

            end if

         end do  ! jat

      end do  ! iat

    end subroutine internal_verify
    !--------------------------------------------------------------------------

    subroutine internal_check_args(have_radius)

      !=======================================================================!
      ! This subroutine checks the arguments to the parent subroutine         !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! have_radius (input) : flag indicating presence of optional argument   !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      use utils, only: utils_assert, utils_abort

      implicit none

      ! Arguments
      logical, intent(in) :: have_radius

      call utils_assert(nat1+nat2 > 0, &
           'Error in parallel_strategy_list_overlaps: no atoms in cell!')
      call utils_assert(size(elements1) == nat1, &
           'Error in parallel_strategy_list_cross_overlaps: incompatible nat1')
      call utils_assert(size(elements2) == nat2, &
           'Error in parallel_strategy_list_cross_overlaps: incompatible nat2')

      if (mode1 /= 'C' .and. mode1 /= 'c' .and. mode1 /= 'R' .and. &
           mode1 /= 'r' .and. mode1 /= 'F' .and. mode1 /= 'f' .and. &
           mode1 /= 'A' .and. mode1 /= 'a' .and. mode1 /= 'L' .and. &
           mode1 /= 'l') then
         call utils_abort('Error in parallel_strategy_list_overlaps: &
              &unknown mode1 "'//trim(mode1)//'"')
      end if
      if (mode2 /= 'C' .and. mode2 /= 'c' .and. mode2 /= 'R' .and. &
           mode2 /= 'r' .and. mode2 /= 'F' .and. mode2 /= 'f' .and. &
           mode2 /= 'A' .and. mode2 /= 'a' .and. mode2 /= 'L' .and. &
           mode2 /= 'l') then
         call utils_abort('Error in parallel_strategy_list_overlaps: &
              &unknown mode2 "'//trim(mode2)//'"')
      end if
      if ((mode1 == 'F' .or. mode1 == 'f' .or. mode2 == 'F' .or. &
           mode2 == 'f') .and. (.not. have_radius)) then
         call utils_abort('Error in parallel_strategy_list_overlaps: &
              &radius must be present for fixed mode')
      end if

    end subroutine internal_check_args

    !--------------------------------------------------------------------------

    subroutine internal_allocate_1

      !=======================================================================!
      ! This subroutine (de)allocates arrays as required by the parent        !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      implicit none

      ! Allocate module array _num_overlaps
      ierr = 0
      if (allocated(overlaps%num_overlaps)) then
         if (size(overlaps%num_overlaps) /= nat1) then
            deallocate(overlaps%num_overlaps,stat=ierr)
            call utils_dealloc_check('internal_allocate_1 &
                 &(parallel_strategy_list_overlaps)','overlaps%num_overlaps',ierr)
            allocate(overlaps%num_overlaps(nat1),stat=ierr)
            call utils_alloc_check('internal_allocate_1 &
                 &(parallel_strategy_list_overlaps)','overlaps%num_overlaps',ierr)
         end if
      else
         allocate(overlaps%num_overlaps(nat1),stat=ierr)
         call utils_alloc_check('internal_allocate_1 &
              &(parallel_strategy_list_overlaps)','overlaps%num_overlaps',ierr)
      end if


      ! Allocate work arrays
      ! rc2013: this could be messy...
      allocate(radii(2,max(nat1,nat2)),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','radii',ierr)
      allocate(atom_subcell1(nat1),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','atom_subcell1',ierr)
      allocate(atom_subcell2(nat2),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','atom_subcell2',ierr)
      allocate(fcoord1(3,nat1),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','fcoord1',ierr)
      allocate(fcoord2(3,nat2),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','fcoord2',ierr)
      allocate(num_my_overlaps(nat1),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','num_my_overlaps',ierr)
      allocate(num_gp_overlaps(nat1),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','num_gp_overlaps',ierr)
      allocate(num_yr_overlaps(nat1),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','num_yr_overlaps',ierr)
      allocate(overlapped(nat2),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','overlapped',ierr)

    end subroutine internal_allocate_1

    subroutine internal_deallocate_1

      !=======================================================================!
      ! This subroutine deallocates arrays as required by the parent          !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Added by Nicholas Hine, 18/05/09 to ensure dealloc check marks        !
      ! the right arrays as having been deallocated when doing array checking !
      !=======================================================================!

      implicit none

      ! Allocate internal arrays
      deallocate(overlapped,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'overlapped',ierr)
      deallocate(num_yr_overlaps,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'num_yr_overlaps',ierr)
      deallocate(num_gp_overlaps,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'num_gp_overlaps',ierr)
      deallocate(num_my_overlaps,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'num_my_overlaps',ierr)
      deallocate(fcoord1,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'fcoord1',ierr)
      deallocate(fcoord2,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'fcoord2',ierr)
      deallocate(atom_subcell1,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'atom_subcell1',ierr)
      deallocate(atom_subcell2,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'atom_subcell2',ierr)
      deallocate(radii,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'radii',ierr)

    end subroutine internal_deallocate_1

  end subroutine parallel_strategy_list_cross_overlaps


end module parallel_strategy
