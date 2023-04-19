!=============================================================================!
!                                                                             !
!  Green's function through sparse-selected-inverse                           !
!                                                                             !
!  The subroutines in this file were written by Robert Bell, May 2014         !
!                                                                             !
!  TCM Group, Cavendish laboratory, University of Cambridge                   !
!  Madingley Road, Cambridge CB3 0HE, UK                                      !
!                                                                             !
!=============================================================================!
!                                                                             !
!  This module implements the sparse-selected-inverse algorithm described in  !
!                                                                             !
!     Petersen et al., J. Comp. Phys. 227 (2008) 3174--3190                   !
!        (eqs. 4-7 and 17-19)                                                 !
!                                                                             !
!  for computing selected blocks of the Hamiltonian Green's function.         !
!                                                                             !
!                                                                             !
!  For one dimensional systems, the algorithm greatly reduces both the memory !
!  and operation count for producing any selected elements, scaling linearly  !
!  system size. Higher dimensional systems are still cubically scaling, but   !
!  the saving in memory/operation count may still be significant when compared!
!  to conventional matrix inversions.                                         !
!                                                                             !
!  The routines should not be used when the full inverse is required as the   !
!  computational cost is larger, and cannot make use of distributed memory.   !
!                                                                             !
!  After a tri-diagonal partitioning scheme is determined, the calculation is !
!  performed in serial on each proc. The intended scope of this module is,    !
!  therefore, to facilitate efficient calculation of the Green's function at  !
!  many points throughout the complex plane. Each complex energy should be    !
!  distributed over all the procs.                                            !
!                                                                             !
!  The selected inverse is used as following:                                 !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!                                                                             !
!     call greenf_init_blocking      ! (determine tri-diagonal partitioning)  !
!     call greenf_alloc_matrices     ! (allocate internal tri-diagonal arrays)!
!     call greenf_init_matrices      ! (fill internal matrices,with H/S data) !
!                                                                             !
!     energy_loop: do                                                         !
!        call greenf_prepare_inverse ! (prepare the inverse at this energy)   !
!                                                                             !
!        ! use selected results of the Green's function, e.g.                 !
!        call greenf_get_atom_block  ! (get a selected block of the matrix)   !
!        call greenf_dos             ! (compute the Density of states)        !
!        ...                                                                  !
!     enddo                                                                   !
!                                                                             !
!     call greenf_destroy            ! (destroy the Green's function storage) !
!                                                                             !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!                                                                             !
!=============================================================================!
!  Assumptions within this module:                                            !
!     1) the SPAM3 Hamiltonian sparsity pattern is the same for all spin      !
!        channels.                                                            !
!     2) Hamiltonian/Overlap matrices are assumed to be Hermitian. The system !
!        is assumed to be time-reversal symmetric.                            !
!                                                                             !
!=============================================================================!
!  Potential future extensions:                                               !
!     1) Add OpenMP level of parallelism. Currently possible through MKL      !
!        threading, but the implementation is not portable.                   !
!     2) Allow conversion of tri-diagonal matrix into a SPAM3 matrix, to allow!
!        the calculation, e.g., of the density kernel through an integral in  !
!        the complex plane.                                                   !
!     3) Remove requirement for requested sparse blocks to exist in the same  !
!        internal tri-diagonal block. The product of c^L, c^R matrices can be !
!        temporarily stored to allow multiple Green's function blocks to be   !
!        computed easily.                                                     !
!     4) Unconnected regions can be treated using a different tri-diagonal    !
!        matrix for each. Currently, a single tri-diagonal matrix is used,    !
!        with the upper/lower off-diagonal coupling elements set to zero.     !
!                                                                             !
!=============================================================================!


module greenf

   use constants, only: DP, LONG
   use rundat, only: pub_debug_on_root

   implicit none

   private

   ! public types
   public :: SPARGF

   ! public subroutines
   public :: greenf_init_blocking
   public :: greenf_init_matrices
   public :: greenf_alloc_blocks
   public :: greenf_destroy
   public :: greenf_prepare_inverse
   public :: greenf_dos
   public :: greenf_memory_estimate

   ! access subroutines
   public :: greenf_get_atom_block

   ! operation subroutines for tri-diagonal matrix
   public :: trimat_deposit_atom_block
   public :: trimat_scale
   ! public :: trimat_to_spam3 ! TODO

   ! debug routines
   public :: greenf_print_gf_obj
   public :: greenf_print_full_inverse
   public :: greenf_print_ham
   public :: tmat_print_matrix

   ! private simple container for arrays of arrays
   type ZMAT
      complex(kind=DP), allocatable :: zmtx(:,:)
   end type ZMAT

   ! private container for block tri-diagonal matrices
   ! note: as only the upper/lower off-diagonal is stored, this can only
   ! refer to Hermitian matrices (the Green's function is not explicitly stored)
   type TRIMAT
      type(ZMAT), allocatable :: d(:) ! the diagonal blocks
      type(ZMAT), allocatable :: c(:) ! the upper/lower off-diagonal blocks
      logical :: upper                ! true if upper off-diagonal stored
   end type TRIMAT

   ! public container for the sparse green function object
   type SPARGF
      !### information on Hamiltonian ###!
      integer              :: num_spins           ! number of spin channels
      character(len=20)    :: ham_type            ! valence/joint Hamiltonian

      !### description of the blocking scheme ###!
      integer              :: nat                 ! number of atoms contained
      integer              :: norb                ! total orbitals in matrix
      integer              :: num_blocks          ! number of diagonal blocks
      integer, allocatable :: block_atom_size(:)  ! atoms in diagonal blocks
      integer, allocatable :: block_norb_size(:)  ! # NGWFS in diagonal blocks
      integer, allocatable :: atom_ngwf_ndx(:)    ! 1st ngwf index for atom,
                                                  ! given as ndx in the block,
                                                  ! not the global index
      integer, allocatable :: force_group_ndx(:)  ! forced group index this
                                                  ! atom belongs in
      logical              :: groups_interact     ! false if force_groups must
                                                  ! have zero coupling elements

      !### locations of atoms in blocks ###!
      ! atom indices are in *internal* order, and array is (1:par%nat)
      ! (atoms that are ignored are given index IGNORED)
      integer, allocatable :: atom_order(:)       ! gf internal ordering of atoms
      integer, allocatable :: index_of_atom(:)    ! reverse of atom_order
      integer, allocatable :: block_of_atom(:)    ! block the atom is in
      integer, allocatable :: ngwfs_on_atom(:)    ! num ngwfs on atom


      !### main static data structures ###!
      type(TRIMAT), allocatable :: h(:)    ! hamiltonian
      type(TRIMAT)              :: s       ! overlap
      complex(kind=DP)          :: e_cmplx ! complex energy

      !### data storage required during the calculation ###!
      type(TRIMAT) :: L_pass ! left pass
      type(TRIMAT) :: R_pass ! right pass

      !### flags to indicate the state of the sparse inverse ###!
      logical :: is_blocked = .false.     ! true if blocking scheme calculated
      logical :: is_alloc = .false.       ! true if all arrays are allocated
      logical :: is_init = .false.        ! true if h/s matrices are filled
      logical :: pass_performed = .false. ! true if L/R passes performed
      integer :: ispin  ! spin channel being calculated in L_pass/R_pass
                        ! and stored in h(:), if only one present

   end type SPARGF

   ! internal labelling associated with the partitioning
   integer, parameter :: UNASSIGNED = -100
   integer, parameter :: IGNORED = -200
   integer, parameter :: NO_GROUP = -300


   contains

   !====================================================================!
   !====================================================================!

   subroutine greenf_init_blocking(gf_obj, ham, ham_type, force_groups, &
                           groups_interact, first_atom, last_atom, ispin)
      !================================================================!
      ! This routine generates the tri-diagonal blocking scheme for    !
      ! the sparse selected inverse using the sparsity pattern of the  !
      ! hamiltonian.                                                   !
      ! The partitioning is done by starting with a seed atom(s), and  !
      ! making a block group consisting of that atom(s) and all its    !
      ! direct neighbouring atoms.                                     !
      ! The next block group is defined as all the atoms that directly !
      ! neighbour any atom in the first group.                         !
      ! This is repeated recursively until all atoms are accounted for.!
      ! If there are subsets of atoms that are completely un-connected,!
      ! a new seed group is defined and the process repeated again.    !
      !                                                                !
      ! For quasi one-dimensional systems, this partitioning will      !
      ! produce blocks of similar size, and the sparse inverse will be !
      ! linear-scaling. For other system dimensionalities, the block   !
      ! sizes will increase each step making the procedure cubically   !
      ! scaling (but the prefactor should still be smaller than a      !
      ! conventional matrix inverse for most systems).                 !
      !----------------------------------------------------------------!
      ! Mandatory arguments:                                           !
      !   gf_obj       (output): a sparse gf object                    !
      !   ham          (input) : the SPAM3 hamiltonian(:)              !
      !   ham_type     (input) : valence or joint hamiltonian          !
      !----------------------------------------------------------------!
      ! Optional arguments:                                            !
      !   force_groups  (input) : list of atoms that are forced to be  !
      !       in the same block group. For efficiency reasons, any     !
      !       selected inverse block must contain atoms in the same    !
      !       block, which can be demanded here.                       !
      !       The format is: force_groups(2,ngroups), where            !
      !       force_groups(:,ii) = (/start, end/) atoms in orig. order !
      !       (default = <no groups>)                                  !
      !   groups_interact (inp.): if .false., the atoms defined in     !
      !       force_groups are defined not to interact (i.e. H, S = 0) !
      !       (default = .true.)                                       !
      !   first_atom    (input) : the first atom to include, in the    !
      !       original input order. (default = 1)                      !
      !   last_atom     (input) : the last  atom to include, in the    !
      !       original input order. (default = par%nat)           !
      !       Atoms outside these ranges are ignored.                  !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use comms, only: pub_on_root
      use constants, only: stdout, DP, NORMAL, VERBOSE, CRLF
      use neighbour_list, only: NL_NEIGHBOUR_LIST, neighbour_list_free, &
           neighbour_list_init_from_sparse
      use parallel_strategy, only: par=>pub_par
      use rundat, only: pub_output_detail, pub_num_spins
      use sparse, only: SPAM3, sparse_num_elems_on_atom
      use timer, only: timer_clock
      use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

      implicit none

      ! arguments
      type(SPARGF), intent(out)     :: gf_obj
      type(SPAM3), intent(in)       :: ham(pub_num_spins)
      character(len=*), intent(in)  :: ham_type
      integer, intent(in), optional :: force_groups(:,:)
      logical, intent(in), optional :: groups_interact
      integer, intent(in), optional :: first_atom
      integer, intent(in), optional :: last_atom
      integer, intent(in), optional :: ispin

      ! internal variables
      type(NL_NEIGHBOUR_LIST) :: nb_list

      integer :: first_atom_loc, last_atom_loc
      integer :: ispin_loc
      integer :: ierr
      integer :: iat, jat, kat ! atom counters
      integer :: this_atm, nb_atm
      integer :: idx
      integer :: iblock
      integer :: igroup
      integer :: orb_count
      integer :: start_ndx, end_ndx, end_ptr
      logical :: atom_found
      integer :: nat_to_assign
      real(kind=DP) :: sparsity_percent

      ! Temporary copy of gf_obj%block_atom_size: is same size as par%nat
      ! as number of blocks is initially unknown
      integer, allocatable :: loc_atoms_in_block(:)


      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering greenf_init_blocking'

      call timer_clock('greenf_init_blocking',1)

      gf_obj%groups_interact = .true.
      if (present(groups_interact)) gf_obj%groups_interact = groups_interact

      first_atom_loc = 1
      last_atom_loc  = 1
      if (present(first_atom)) first_atom_loc = first_atom
      if (present(last_atom))  last_atom_loc  = last_atom

      if (present(force_groups)) then
         do igroup = 1,size(force_groups,dim=2)
            call utils_assert(force_groups(1,igroup) .ge. first_atom_loc, &
                 'Error in greenf_init_blocking: force_groups(1,igroup) < &
                 &first_atom.',igroup)
            call utils_assert(force_groups(2,igroup) .le. last_atom_loc, &
                 'Error in greenf_init_blocking: force_groups(2,igroup) < &
                 &last_atom.',igroup)
         enddo
      endif

      ! pick if only selecting one spin channel
      gf_obj%num_spins = pub_num_spins
      ispin_loc = pub_num_spins
      if (present(ispin)) then
         gf_obj%num_spins = 1
         ispin_loc = 1
      endif


      ! fill in information about this greenf's function
      gf_obj%nat = last_atom - first_atom + 1
      gf_obj%ham_type = ham_type

      if(pub_on_root .and. pub_output_detail>=NORMAL) then
         write(stdout,'(/a)',advance='no') ' Computing blocking scheme for &
              &tridiagonal matrix inversion ...'
         if (pub_output_detail >= VERBOSE) write(stdout,*)
      endif


      ! from the hamiltonian, build the nearest neighbour
      ! assumes both spin channels have the same SPAM3 structure
      call neighbour_list_init_from_sparse(nb_list,'ham',ham(ispin_loc))

      ! allocate storage for the blocks the atoms are assigned to
      allocate(gf_obj%block_of_atom(par%nat),stat=ierr)
      call utils_alloc_check('greenf_init_blocking','gf_obj%block_of_atom',ierr)
      allocate(gf_obj%atom_order(par%nat),stat=ierr)
      call utils_alloc_check('greenf_init_blocking','gf_obj%atom_order',ierr)
      allocate(loc_atoms_in_block(par%nat),stat=ierr)
      call utils_alloc_check('greenf_init_blocking','loc_atoms_in_block',ierr)
      allocate(gf_obj%force_group_ndx(par%nat),stat=ierr)
      call utils_alloc_check('greenf_init_blocking',&
           'gf_obj%force_group_ndx',ierr)


      ! Assign each atom to a block based on the neighbour list
      ! initialise all atoms as unassigned
      gf_obj%block_of_atom(:) = UNASSIGNED
      gf_obj%atom_order(:) = IGNORED
      loc_atoms_in_block(:) = 0
      gf_obj%force_group_ndx(:) = NO_GROUP

      ! set atoms that are ignored
      if (present(first_atom) .and. present(last_atom)) then
         do iat=1,first_atom-1
            jat = par%distr_atom(iat)
            gf_obj%block_of_atom(jat) = IGNORED
         enddo
         do iat=last_atom+1,par%nat
            jat = par%distr_atom(iat)
            gf_obj%block_of_atom(jat) = IGNORED
         enddo
      endif

      ! set atoms that are in a lead group
      if (present(force_groups)) then
         do igroup = 1,size(force_groups,dim=2)
            do iat=force_groups(1,igroup),force_groups(2,igroup)
               jat = par%distr_atom(iat)
               gf_obj%force_group_ndx(jat) = igroup
            enddo
         enddo
      endif


      ! seed the first block...
      if (present(force_groups)) then
         ! ...using the seed block if present
         kat = 1
         do iat=force_groups(1,1),force_groups(2,1)
            jat = par%distr_atom(iat)
            gf_obj%block_of_atom(jat) = 1
            gf_obj%atom_order(kat) = jat
            kat = kat + 1
         enddo
         loc_atoms_in_block(1) = kat - 1
      else
         ! ...otherwise use the first atom and its neighbours
         kat = 1
         this_atm = par%distr_atom(first_atom_loc)
         gf_obj%block_of_atom(this_atm) = 1
         do idx = nb_list%first_idx(this_atm), nb_list%last_idx(this_atm)
            nb_atm = nb_list%neighbours(idx)
            gf_obj%block_of_atom(nb_atm) = 1
            gf_obj%atom_order(kat) = nb_atm
            kat = kat + 1
         enddo
         loc_atoms_in_block(1) = kat - 1
      endif



      ! recursively loop through atoms connected to the previous seed group
      iblock = 1
      nat_to_assign = count(gf_obj%block_of_atom(:) == UNASSIGNED)
      block_loop: &
      do

         atom_found = .false.

         ! note: in F95, sum of zero sized array = 0
         start_ndx = sum(loc_atoms_in_block(1:iblock-1)) + 1
         end_ndx = start_ndx + loc_atoms_in_block(iblock) - 1
         end_ptr = end_ndx + 1

         ! loop over all atoms in the current block
         do iat = start_ndx, end_ndx

            ! for each atom in this group
            this_atm = gf_obj%atom_order(iat)


            ! loop over all neighbouring atoms
            do idx=nb_list%first_idx(this_atm),nb_list%last_idx(this_atm)

               nb_atm = nb_list%neighbours(idx) ! a neighbouring atom

               ! if the atom is unassigned, assign it
               if (gf_obj%block_of_atom(nb_atm) == UNASSIGNED) then

                  ! if not in group, deal individually
                  if (gf_obj%force_group_ndx(nb_atm) == NO_GROUP) then
                     atom_found = .true.
                     gf_obj%block_of_atom(nb_atm) = iblock+1
                     gf_obj%atom_order(end_ptr) = nb_atm
                     end_ptr = end_ptr + 1
                     nat_to_assign = nat_to_assign - 1

                  ! otherwise, set all atoms in that group to this index
                  else

                     ! if the seed and neighbour are both lead atoms
                     ! then check if there should be no connection
                     ! note: if one NO_GROUP block neighbours two unassigned
                     ! forced groups, this will not remove the interaction
                     ! between the forced groups. This must be dealt with
                     ! in greenf_init_matrices
                     if (gf_obj%force_group_ndx(this_atm) /= NO_GROUP .and. &
                          .not.gf_obj%groups_interact) cycle

                     atom_found = .true.

                     ! find which group it is in
                     igroup = gf_obj%force_group_ndx(nb_atm)

                     ! set atoms
                     do kat = force_groups(1,igroup), force_groups(2,igroup)
                        jat = par%distr_atom(kat)
                        gf_obj%block_of_atom(jat) = iblock+1
                        gf_obj%atom_order(end_ptr) = jat
                        end_ptr = end_ptr + 1
                        nat_to_assign = nat_to_assign - 1
                     enddo
                  endif
               endif
            enddo

         enddo

         ! if no new atoms are found, then any remaining atoms are unconnected.
         ! The next block is chosen as the first unassigned atom and neighbours
         if (.not.atom_found) then
            if (pub_on_root) write(stdout,'(3x,a)') &
                 'Warning: unconnected atoms detected in system partitioning'

            ! find first unassigned atom
            do this_atm=1,par%nat
               if (gf_obj%block_of_atom(this_atm) == UNASSIGNED) exit
            enddo


            ! add this atom and all its neighbouring atoms
            do idx = nb_list%first_idx(this_atm), nb_list%last_idx(this_atm)

               nb_atm = nb_list%neighbours(idx) ! a neighbouring atom

               if (gf_obj%block_of_atom(nb_atm) == IGNORED) cycle

               gf_obj%block_of_atom(nb_atm) = iblock + 1
               gf_obj%atom_order(end_ptr) = nb_atm
               end_ptr = end_ptr + 1
               nat_to_assign = nat_to_assign - 1
            enddo

         endif


         ! increment block counter
         iblock = iblock + 1
         loc_atoms_in_block(iblock) = end_ptr - end_ndx - 1

         ! if all atoms are accounted for, exit iteration
         if (nat_to_assign == 0) exit block_loop
         call utils_assert(nat_to_assign>0,'Error in greenf_init_blocking: &
              &nat_to_assign<0. should not reach here.',nat_to_assign)


      enddo block_loop

      if(pub_on_root.and.pub_output_detail>=NORMAL) write(stdout,'(a/)') ' done'

      ! total number of blocks now given by iblock
      gf_obj%num_blocks = iblock
      call utils_assert(iblock > 1, 'Error in greenf_init_blocking: cannot &
           &construct a blocking scheme as all atoms connected.'//CRLF// &
           'Selected inverse is much less efficient than the direct inverse')


      ! fill in data about the storage
      allocate(gf_obj%block_atom_size(gf_obj%num_blocks),stat=ierr)
      call utils_alloc_check('greenf_init_blocking','gf_obj%block_atom_size',ierr)
      allocate(gf_obj%index_of_atom(par%nat),stat=ierr)
      call utils_alloc_check('greenf_init_blocking','gf_obj%index_of_atom',ierr)
      allocate(gf_obj%block_norb_size(gf_obj%num_blocks),stat=ierr)
      call utils_alloc_check('greenf_init_blocking','gf_obj%block_norb_size',ierr)
      allocate(gf_obj%atom_ngwf_ndx(par%nat),stat=ierr)
      call utils_alloc_check('greenf_init_blocking','gf_obj%atom_ngwf_ndx',ierr)
      allocate(gf_obj%ngwfs_on_atom(par%nat),stat=ierr)
      call utils_alloc_check('greenf_init_blocking','gf_obj%ngwfs_on_atom',ierr)

      gf_obj%index_of_atom(:) = IGNORED
      gf_obj%atom_ngwf_ndx(:) = IGNORED
      gf_obj%ngwfs_on_atom(:) = IGNORED
      gf_obj%block_norb_size(:) = UNASSIGNED

      do iat=1,gf_obj%nat
         gf_obj%index_of_atom(gf_obj%atom_order(iat)) = iat
      enddo
      gf_obj%block_atom_size(:) = loc_atoms_in_block(1:gf_obj%num_blocks)


      ! find number of orbitals per block, and first ngwf of atom
      iblock = 1
      orb_count = 0
      do iat = 1,gf_obj%nat
         jat = gf_obj%atom_order(iat)

         ! if this atom is in a new block, reset counters
         if (iblock /= gf_obj%block_of_atom(jat)) then
            orb_count = 0
            iblock = gf_obj%block_of_atom(jat)
         endif

         ! set the index of this atom
         gf_obj%atom_ngwf_ndx(jat) = orb_count + 1
         gf_obj%ngwfs_on_atom(jat) = &
              sparse_num_elems_on_atom(jat,ham(ispin_loc),'R')
         orb_count = orb_count + gf_obj%ngwfs_on_atom(jat)
         gf_obj%block_norb_size(iblock) = orb_count

      enddo
      gf_obj%norb = sum(gf_obj%block_norb_size(:))


      ! store sparsity percentage
      sparsity_percent = 0.0_DP
      do iblock = 1, gf_obj%num_blocks   ! diagonal blocks
         sparsity_percent = sparsity_percent + &
              real(gf_obj%block_norb_size(iblock),kind=DP)**2
      enddo
      do iblock = 1, gf_obj%num_blocks-1 ! off-diagonal blocks
         sparsity_percent = sparsity_percent + &
              2*0_DP*real(gf_obj%block_norb_size(iblock),kind=DP)*&
              real(gf_obj%block_norb_size(iblock+1),kind=DP)
      enddo
      sparsity_percent = sparsity_percent / real(gf_obj%norb,kind=DP)**2
      sparsity_percent = sparsity_percent * 100.0_DP



      ! print information about the blocking
      if (pub_on_root .and. pub_output_detail >= NORMAL) then
         write(stdout,'(a)')   ' Information on tri-diagonal blocking scheme:'
         write(stdout,'(a,i8)') ' Number of blocks found: ', gf_obj%num_blocks
         iblock = maxloc(gf_obj%block_norb_size,dim=1)
         write(stdout,'(a,2i8)')   ' Maximum block size (atoms, orbitals): ', &
              gf_obj%block_atom_size(iblock), gf_obj%block_norb_size(iblock)
         write(stdout,'(a,f8.3,a)') ' Tridiagonal sparsity ', &
              sparsity_percent ,' %'
         if (pub_output_detail >= VERBOSE) then
            write(stdout,'(a)') '   Block sizes: (atoms, orbitals)'
            do iblock=1,gf_obj%num_blocks
               write(stdout,'(2i8)')   gf_obj%block_atom_size(iblock), &
                                       gf_obj%block_norb_size(iblock)
            enddo
         endif
         write(stdout,*)
      endif

      if (pub_debug_on_root) then
         write(stdout,'(a)') '---DEBUG---'
         write(stdout,'(6a12)') 'orig_atm','distr_atm','atm_order', &
                                'atm_block','1st_ngwf','atm_ngwfs'
         do iat=1,par%nat
            jat = par%distr_atom(iat)
            kat = gf_obj%index_of_atom(jat)
            if (kat /= IGNORED) then
               write(stdout,'(6i12)',advance='no') iat, jat, &
                    gf_obj%index_of_atom(jat), gf_obj%block_of_atom(jat), &
                    gf_obj%atom_ngwf_ndx(jat), gf_obj%ngwfs_on_atom(jat)
               if (gf_obj%force_group_ndx(jat) /= NO_GROUP) then
                  write(stdout,'(a,i4)') ' <-- Forced group ', &
                       gf_obj%force_group_ndx(jat)
               else
                  write(stdout,*)
               endif
            else
               write(stdout,'(i12,a12)') iat, 'IGNORED'
            endif
         enddo
         write(stdout,'(a)') '---DEBUG---'
      endif

      ! print all info before checking that there are at least 2 blocks
      call utils_assert(gf_obj%num_blocks .gt. 1, 'Error in greenf_init_&
            &blocking: sparse selective inverse for less than two blocks&
            & is untested. It will also be very inefficient. Use with caution!')

      if (gf_obj%num_blocks .lt. 3 .and. pub_on_root) then
         write(stdout,'(a)') 'Warning: less than 3 blocks detected. Sparse &
              &selected inverse will be very inefficient, but will continue'
      endif

      ! free up neighbour list
      call neighbour_list_free(nb_list)

      ! set flag
      gf_obj%is_blocked = .true.

      deallocate(loc_atoms_in_block,stat=ierr)
      call utils_dealloc_check('greenf_init_blocking','loc_atoms_in_block',ierr)

      call timer_clock('greenf_init_blocking',2)

    if (pub_debug_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving greenf_init_blocking'

   end subroutine greenf_init_blocking

   !====================================================================!
   !====================================================================!

   subroutine greenf_destroy(gf_obj)
      !================================================================!
      ! Destroy the arrays containing the blocking structure.          !
      ! Also destroys the tri-diagonal matrices, if still allocated.   !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use utils, only: utils_dealloc_check, utils_assert

      implicit none

      ! arguments
      type(SPARGF), intent(inout) :: gf_obj

      ! internal
      integer :: ierr

      call utils_assert(gf_obj%is_blocked,"Error in greenf_destroy: &
           &Green's function hasn't been created.")

      ! deallocate storage of block tri-diagonal matrices
      if (gf_obj%is_alloc) then
         call greenf_dealloc_blocks(gf_obj)
      endif

      ! set state flags
      gf_obj%is_blocked = .false.
      gf_obj%is_init = .false.
      gf_obj%pass_performed = .false.

      ! deallocate storage for the blocks the atoms are assigned to
      deallocate(gf_obj%block_of_atom,stat=ierr)
      call utils_dealloc_check('greenf_destroy','gf_obj%block_of_atom',ierr)
      deallocate(gf_obj%atom_order,stat=ierr)
      call utils_dealloc_check('greenf_destroy','gf_obj%atom_order',ierr)
      deallocate(gf_obj%index_of_atom,stat=ierr)
      call utils_dealloc_check('greenf_destroy','gf_obj%index_of_atom',ierr)
      deallocate(gf_obj%atom_ngwf_ndx,stat=ierr)
      call utils_dealloc_check('greenf_destroy','gf_obj%atom_ngwf_ndx',ierr)
      deallocate(gf_obj%block_atom_size,stat=ierr)
      call utils_dealloc_check('greenf_destroy','gf_obj%block_atom_size',ierr)
      deallocate(gf_obj%block_norb_size,stat=ierr)
      call utils_dealloc_check('greenf_destroy','gf_obj%block_norb_size',ierr)
      deallocate(gf_obj%ngwfs_on_atom,stat=ierr)
      call utils_dealloc_check('greenf_destroy','gf_obj%ngwfs_on_atom',ierr)
      deallocate(gf_obj%force_group_ndx,stat=ierr)
      call utils_dealloc_check('greenf_destroy','gf_obj%force_group_ndx',ierr)

   end subroutine greenf_destroy

   !====================================================================!
   !====================================================================!

   subroutine greenf_init_matrices(gf_obj,ham,overlap,reuse_overlap,ispin)
      !================================================================!
      ! Routine to deposit the hamiltonian and overlap matrices into   !
      ! the tri-diagonal matrices generated according to the internal  !
      ! partitioning scheme.
      !----------------------------------------------------------------!
      ! Arguments:                                                     !
      !   gf_obj       (inout) : a sparse gf object                    !
      !   ham          (input) : the SPAM3 hamiltonian(:) to deposit   !
      !   overlap      (input) : the SPAM3 overlap, S, to deposit      !
      !   reuse_overlap (optional inp.): if .true., do not re-read S   !
      !   ispin (optional inp.): the spin channel of ham to use if     !
      !                          gf_obj does not contain both channels !
      !----------------------------------------------------------------!
      ! Necessary conditions:                                          !
      !    gf_obj%is_blocked = .true.                                  !
      !    gf_obj%is_alloc = .true.                                    !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use comms, only: pub_my_proc_id
      use constants, only: stdout, cmplx_0
      use parallel_strategy, only: par=>pub_par
      use rundat, only: pub_num_spins
      use sparse, only: SPAM3
      use timer, only: timer_clock
      use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

      implicit none

      ! arguments
      type(SPARGF), intent(inout)   :: gf_obj
      type(SPAM3), intent(in)       :: ham(pub_num_spins)
      type(SPAM3), intent(in)       :: overlap
      logical, intent(in), optional :: reuse_overlap ! do not gather overlap
      integer, intent(in), optional :: ispin     ! spin channel of ham to take

      ! internal
      integer :: iblock ,jblock
      integer :: iat, jat, ii, jj, iat_orig
      integer :: is
      integer :: igrp, jgrp
      real(kind=DP), allocatable :: dbuffer(:,:)
      integer :: ierr
      logical :: reuse_overlap_loc

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering greenf_init_matrices'

      call timer_clock('greenf_init_matrices',1)


      call utils_assert(gf_obj%is_blocked, 'Error in greenf_init_matrices: &
           & greenf_init_blocking has not been called.')

      call utils_assert(gf_obj%is_alloc, 'Error in greenf_init_matrices: &
           & greenf_alloc_blocks has not been called.')

      ! deal with optional arguments
      reuse_overlap_loc = .false.
      if (present(reuse_overlap)) reuse_overlap_loc = reuse_overlap

      if (present(ispin)) then
         gf_obj%ispin = ispin
      else
         gf_obj%ispin = -1 ! flag that multiple spins present
      endif

      call utils_assert( present(ispin).or.pub_num_spins == &
           gf_obj%num_spins,'Error in greenf_init_matrices: ispin must be &
           &defined if gf spins /= system spins')

      gf_obj%pass_performed = .false.


      ! clear current contents of h and s
      if (present(ispin)) then
         call trimat_scale(gf_obj%h(1),cmplx_0)
      else
         call trimat_scale(gf_obj%h(ispin),cmplx_0)
      endif
      if (.not.reuse_overlap_loc) call trimat_scale(gf_obj%s,cmplx_0)


      ! allocate workspace for dbuffer, will reallocate as necessary
      allocate(dbuffer(1,1),stat=ierr)
      call utils_alloc_check('greenf_init_matrices','dbuffer',ierr)


      ! loop through all atoms in each block
      ! this outer loop is a loop over COLUMNS
      do ii = 1, gf_obj%nat

         ! atom, its block index, and its forced group index
         iat = gf_obj%atom_order(ii)
         iblock = gf_obj%block_of_atom(iat)
         igrp = gf_obj%force_group_ndx(iat)

         ! if atom is not on my proc, cycle
         iat_orig = par%orig_atom(iat)
         if (par%proc_of_atom(iat_orig) .ne. pub_my_proc_id) cycle


         ! loop through all other atoms belonging to this/next block
         ! this inner loop is a loop over ROWS
         do jj = ii, gf_obj%nat
            jat = gf_obj%atom_order(jj)
            jblock = gf_obj%block_of_atom(jat)
            jgrp = gf_obj%force_group_ndx(jat)

            ! exit if not in next block
            ! (can exit as atoms are ordered by block index --> O(N) scaling)
            if (jblock .gt. iblock + 1) exit

            ! cycle if these atoms should not interact
            if (.not.gf_obj%groups_interact .and. igrp /= jgrp .and. &
                igrp /= NO_GROUP .and. jgrp /= NO_GROUP ) cycle

            ! otherwise, get block and deposit

            ! Hamiltonian
            if (present(ispin)) then ! single spin requested
               call internal_deposit(gf_obj%h(1),ham(ispin),iat,jat,iblock,jblock)
               gf_obj%ispin = ispin
            else ! get all spins
               do is=1,gf_obj%num_spins
                  call internal_deposit(gf_obj%h(is),ham(is),iat,jat,iblock,jblock)
               enddo
               gf_obj%ispin = -1 ! flag that multiple spins present
            endif

            ! overlap, if requested
            if (.not. reuse_overlap_loc) then
               call internal_deposit(gf_obj%s,overlap,iat,jat,iblock,jblock)
            endif

         enddo

      enddo

      ! sum up contributions
      do is = 1, gf_obj%num_spins
         call trimat_reduce(gf_obj%h(is))
      enddo
      if (.not. reuse_overlap_loc) then
         call trimat_reduce(gf_obj%s)
      endif

      ! tidy up
      if (allocated(dbuffer)) then
         deallocate(dbuffer,stat=ierr)
         call utils_dealloc_check('greenf_init_matrices','dbuffer',ierr)
      endif

      ! set flag
      gf_obj%is_init = .true.

      call timer_clock('greenf_init_matrices',2)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving greenf_init_matrices'

      contains

      subroutine internal_deposit(tmat,mat,col_at,row_at,col_block,row_block)
         !=============================================================!
         ! Deposits the (col_at,row_at) atom block of the SPAM3 mat in !
         ! the tri-diagonal matrix tmat.                               !
         !=============================================================!

         use sparse, only: SPAM3, sparse_get_block, sparse_num_elems_on_atom
         use utils, only: utils_abort

         implicit none

         ! arguments
         type(TRIMAT), intent(inout) :: tmat
         type(SPAM3), intent(in)     :: mat
         integer, intent(in)         :: col_at, row_at, col_block, row_block

         ! internal
         integer :: row_ngwfs, col_ngwfs
         integer :: istart, jstart, iend, jend


         ! get index of the subblock in the tridiagonal block
         row_ngwfs = sparse_num_elems_on_atom(row_at,mat,'R')
         col_ngwfs = sparse_num_elems_on_atom(col_at,mat,'R')
         istart = gf_obj%atom_ngwf_ndx(col_at)
         jstart = gf_obj%atom_ngwf_ndx(row_at)
         iend = istart + col_ngwfs - 1
         jend = jstart + row_ngwfs - 1

         ! check the index is found
         call utils_assert(istart /= UNASSIGNED,'Error in greenf_init_matrices:&
              & UNASSIGNED encountered for atom ',iat)
         call utils_assert(jstart /= UNASSIGNED,'Error in greenf_init_matrices:&
              &UNASSIGNED encountered for atom ',jat)

         ! check dbuffer is large enough
         ! (note: swapping of row/col is intentional)
         if(size(dbuffer,dim=1) .ne. row_ngwfs .or. &
            size(dbuffer,dim=2) .ne. col_ngwfs) then

            deallocate(dbuffer,stat=ierr)
            call utils_dealloc_check('greenf_init_matrices','dbuffer',ierr)
            allocate(dbuffer(row_ngwfs,col_ngwfs),stat=ierr)
            call utils_alloc_check('greenf_init_matrices','dbuffer',ierr)

         endif

         ! get block
         call sparse_get_block(dbuffer,mat,row_at,col_at)

         ! note: sparse_get_block returns in (row,col) order, not the
         ! usual Fortran (col,row) order.
         ! Wherever dbuffer is used, it is therefore transposed

         ! deposit
         if (row_block == col_block) then         ! diagonal block

            tmat%d(iblock)%zmtx(istart:iend,jstart:jend) = &
                 transpose(cmplx(dbuffer(1:row_ngwfs,1:col_ngwfs),kind=DP))

            ! deposit matrix on upper and lower diagonals in this block
            if (col_at /= row_at) then
               tmat%d(iblock)%zmtx(jstart:jend,istart:iend) = &
                    cmplx(dbuffer(1:row_ngwfs,1:col_ngwfs),kind=DP)
            endif
         elseif (row_block == col_block + 1) then ! off-diagonal block

            if (tmat%upper) then
               tmat%c(iblock)%zmtx(istart:iend,jstart:jend) = &
                    transpose(cmplx(dbuffer(1:row_ngwfs,1:col_ngwfs),kind=DP))
            else
               tmat%c(iblock)%zmtx(jstart:jend,istart:iend) = &
                    cmplx(dbuffer(1:row_ngwfs,1:col_ngwfs),kind=DP)
            endif
         else
            call utils_abort('Error in greenf_init_matrices: block indices do &
                 &not correspond to tri-diagonal structure', &
                 iat, jat, col_block, row_block)
         endif

      end subroutine internal_deposit

   end subroutine greenf_init_matrices

   !====================================================================!
   !====================================================================!

   subroutine greenf_alloc_blocks(gf_obj)
      !================================================================!
      ! Allocate the GF internal tri-diagonal matrix structures.       !
      !----------------------------------------------------------------!
      ! Arguments:                                                     !
      !   gf_obj       (input) : a sparse gf object                    !
      !----------------------------------------------------------------!
      ! Necessary conditions:                                          !
      !    gf_obj%is_blocked = .true.                                  !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: stdout
      use utils, only: utils_alloc_check, utils_assert

      implicit none

      ! arguments
      type(SPARGF), intent(inout) :: gf_obj

      ! internal
      integer :: ispin
      integer :: ierr

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering greenf_alloc_blocks'

      call utils_assert(gf_obj%is_blocked,'Error in greenf_alloc_blocks: &
           &greenf_init_blocking not called')

      ! hamiltonian
      allocate(gf_obj%h(gf_obj%num_spins),stat=ierr)
      call utils_alloc_check('greenf_alloc_blocks','gf_obj%h',ierr)

      do ispin=1,gf_obj%num_spins
         call internal_alloc_trimat(gf_obj%h(ispin),'h',.true.)
      enddo

      call internal_alloc_trimat(gf_obj%s,'s',.true.)
      call internal_alloc_trimat(gf_obj%L_pass,'L_pass',.false.)
      call internal_alloc_trimat(gf_obj%R_pass,'R_pass',.true.)

      ! set state flag
      gf_obj%is_alloc = .true.

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving greenf_alloc_blocks'

      contains

      subroutine internal_alloc_trimat(tmat,name,upper)
         !=============================================================!
         ! Allocates the tri-diagonal matrix tmat, with name 'name'    !
         ! If 'upper', off-diagonal blocks are the upper blocks        !
         !=============================================================!
         use constants, only: cmplx_0

         implicit none

         ! arguments
         type(TRIMAT), intent(inout)  :: tmat
         character(len=*), intent(in) :: name
         logical, intent(in)          :: upper

         ! internal
         integer :: iblock
         integer :: iorb, jorb
         integer :: ierr

         allocate(tmat%d(gf_obj%num_blocks),stat=ierr)
         call utils_alloc_check('greenf_alloc_blocks','tmat%d',ierr)
         allocate(tmat%c(gf_obj%num_blocks - 1),stat=ierr)
         call utils_alloc_check('greenf_alloc_blocks','tmat%c',ierr)

         ! store type
         tmat%upper = upper

         ! diagonal blocks
         do iblock = 1, gf_obj%num_blocks
            iorb = gf_obj%block_norb_size(iblock)

            allocate(tmat%d(iblock)%zmtx(iorb,iorb),stat=ierr)
            call utils_alloc_check('greenf_alloc_blocks','tmat%d(iblock)%zmtx',ierr)

            ! init to zero
            tmat%d(iblock)%zmtx = cmplx_0
         enddo

         ! off-diagonal blocks
         do iblock = 1, gf_obj%num_blocks - 1

            if (upper) then ! upper diagonal stored
               iorb = gf_obj%block_norb_size(iblock)
               jorb = gf_obj%block_norb_size(iblock + 1)
            else            ! lower diagonal stored
               iorb = gf_obj%block_norb_size(iblock + 1)
               jorb = gf_obj%block_norb_size(iblock)
            endif

            allocate(tmat%c(iblock)%zmtx(iorb,jorb),stat=ierr)
            call utils_alloc_check('greenf_alloc_blocks','tmat%c(iblock)%zmtx',ierr)

            ! init to zero
            tmat%c(iblock)%zmtx = cmplx_0
         enddo

      end subroutine internal_alloc_trimat

   end subroutine greenf_alloc_blocks

   !====================================================================!
   !====================================================================!

   subroutine greenf_dealloc_blocks(gf_obj)
      !================================================================!
      ! Deallocate the GF internal tri-diagonal matrix structures.     !
      !----------------------------------------------------------------!
      ! Arguments:                                                     !
      !   gf_obj       (input) : a sparse gf object                    !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: stdout
      use utils, only: utils_dealloc_check

      implicit none

      ! arguments
      type(SPARGF), intent(inout) :: gf_obj

      ! internal
      integer :: ispin
      integer :: ierr

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering greenf_dealloc_blocks'

      do ispin=1,gf_obj%num_spins
         call internal_dealloc_trimat(gf_obj%h(ispin),'h')
      enddo

      call internal_dealloc_trimat(gf_obj%s,'s')
      call internal_dealloc_trimat(gf_obj%L_pass,'L_pass')
      call internal_dealloc_trimat(gf_obj%R_pass,'R_pass')

      deallocate(gf_obj%h,stat=ierr)
      call utils_dealloc_check('greenf_dealloc_blocks','gf_obj%h',ierr)

      ! set state flag
      gf_obj%is_alloc = .false.

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving greenf_dealloc_blocks'


      contains

      subroutine internal_dealloc_trimat(tmat,name)
         !=============================================================!
         ! Deallocates the tri-diagonal matrix tmat, with name 'name'  !
         !=============================================================!

         implicit none

         ! arguments
         type(TRIMAT), intent(inout)  :: tmat
         character(len=*), intent(in) :: name

         ! internal
         integer :: iblock
         integer :: ierr

         do iblock = 1, gf_obj%num_blocks
            deallocate(tmat%d(iblock)%zmtx,stat=ierr)
            call utils_dealloc_check('greenf_dealloc_blocks','tmat%d(iblock)%zmtx',ierr)
         enddo
         do iblock = 1, gf_obj%num_blocks - 1
            deallocate(tmat%c(iblock)%zmtx,stat=ierr)
            call utils_dealloc_check('greenf_dealloc_blocks','tmat%c(iblock)%zmtx',ierr)
         enddo

         deallocate(tmat%d,stat=ierr)
         call utils_dealloc_check('greenf_dealloc_blocks','tmat%d',ierr)
         deallocate(tmat%c,stat=ierr)
         call utils_dealloc_check('greenf_dealloc_blocks','tmat%c',ierr)

      end subroutine internal_dealloc_trimat

   end subroutine greenf_dealloc_blocks

   !====================================================================!
   !====================================================================!

   subroutine greenf_prepare_inverse(gf_obj,cmplx_energy,info,ispin)
      !================================================================!
      ! Prepare the calculation of the Green's function using the left !
      ! and right passes defined in equations 4--7 of Petersen et al.  !
      ! J. Comp. Phys 227 (2008) 3174-3190.                            !
      !                                                                !
      ! Tri-diagonal blocking scheme must be initalised using          !
      ! greenf_init_blocking.                                          !
      ! Tri-diagonal matrices must be allocated and initialised using  !
      ! greenf_alloc_blocks, and greenf_init_matrices.                 !
      !----------------------------------------------------------------!
      ! Arguments:                                                     !
      !   gf_obj       (input) : a sparse gf object                    !
      !   cmplx_energy (input) : the GF complex energy                 !
      !   info         (output): zero if successful                    !
      !   ispin (optional inp.): the spin channel of gf_obj%h(:) to use!
      !                          If not present, will default to first !
      !                          channel.                              !
      !----------------------------------------------------------------!
      ! Necessary conditions:                                          !
      !   gf_obj%is_blocked = .true.                                   !
      !   gf_obj%is_alloc = .true.                                     !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: DP, cmplx_1, stdout
      use linalg, only: linalg_invert_serial, linalg_mat_mul_serial
      use timer, only: timer_clock
      use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check

      implicit none

      ! Arguments
      type(SPARGF), intent(inout)   :: gf_obj
      complex(kind=DP), intent(in)  :: cmplx_energy
      integer, intent(out)          :: info
      integer, intent(in), optional :: ispin

      ! internal
      integer :: ispin_loc
      integer :: iblock
      integer :: prev_norb     ! size of (n-1) diagonal block
      integer :: next_norb     ! size of n diagonal block
      integer :: ierr
      complex(kind=DP), allocatable :: dbuff(:,:) ! diagonal block buffer
      complex(kind=DP), allocatable :: cbuff(:,:) ! off-diagonal block buffer


      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering greenf_prepare_inverse'

      call timer_clock('greenf_prepare_inverse',1)

      ! reset the gf_obj: cannot get block unless this routine succeeds
      gf_obj%pass_performed = .false.

      ! ensure that blocks are allocated and populated
      call utils_assert(gf_obj%is_init.and.gf_obj%is_alloc, &
           'Error in greenf_prepare_inverse: tri-diagonal matrices are not &
           &initialised.')

      ! deal with optional spin parameter
      ispin_loc = 1
      if (present(ispin)) ispin_loc = ispin

      call utils_assert(ispin_loc .le. gf_obj%num_spins,'Error in greenf&
           &_prepare_inverse: incorrect ispin.',ispin_loc,gf_obj%num_spins)

      gf_obj%ispin = ispin_loc
      gf_obj%e_cmplx = cmplx_energy

      ! initialise L_pass/R_pass to the inverse Green's function = a = E*S-H
      call internal_build_inv_gf(gf_obj%L_pass,gf_obj%h(ispin_loc),gf_obj%s)
      call internal_build_inv_gf(gf_obj%R_pass,gf_obj%h(ispin_loc),gf_obj%s)

      ! allocate buffers, will be reallocated if necessary
      prev_norb = gf_obj%block_norb_size(1)
      next_norb = gf_obj%block_norb_size(2)
      allocate(dbuff(prev_norb,prev_norb),stat=ierr)
      call utils_alloc_check('greenf_prepare_inverse','dbuff',ierr)
      allocate(cbuff(next_norb,prev_norb),stat=ierr)
      call utils_alloc_check('greenf_prepare_inverse','cbuff',ierr)

      ! do left pass first
      ! will use the labelling assuming solving for iblock=3
      ! note that off-diagonal is lower diagonal block, i.e. a_32
      ! the upper off-diagonal is a_23 = conjg(transpose(a_32))
      do iblock = 2,gf_obj%num_blocks
         prev_norb = gf_obj%block_norb_size(iblock-1)
         next_norb = gf_obj%block_norb_size(iblock)

         ! reallocate buffers if necessary
         deallocate(dbuff,stat=ierr)
         call utils_dealloc_check('greenf_prepare_inverse','dbuff',ierr)
         allocate(dbuff(prev_norb,prev_norb),stat=ierr)
         call utils_alloc_check('greenf_prepare_inverse','dbuff',ierr)
         deallocate(cbuff,stat=ierr)
         call utils_dealloc_check('greenf_prepare_inverse','cbuff',ierr)
         allocate(cbuff(next_norb,prev_norb),stat=ierr)
         call utils_alloc_check('greenf_prepare_inverse','cbuff',ierr)


         ! get inverse of previous diagonal block
         ! dbuff = inv(d_22)
         dbuff(:,:) = gf_obj%L_pass%d(iblock-1)%zmtx(:,:)
         call linalg_invert_serial(dbuff,info)

         ! if inverse fails, return without setting pass_performed
         if (info /= 0) goto 999

         ! compute previous off-diagonal block
         ! cbuff = -a_32 * inv(d_22)
         call linalg_mat_mul_serial(cbuff, &                            ! out
              gf_obj%L_pass%c(iblock-1)%zmtx, dbuff, alpha=-cmplx_1)

         ! compute new diagonal block
         ! d_33 = a_33 - a_32 * inv(d_22) * a_23
         !      = a_33 + cbuff * a_23
         ! must transpose a_32 to get a_23
         call linalg_mat_mul_serial(gf_obj%L_pass%d(iblock)%zmtx, &     ! out
              cbuff, gf_obj%L_pass%c(iblock-1)%zmtx, beta=cmplx_1, opB='T')

         ! deposit cbuff = new a_32 (i.e. lower off-diagonal)
         gf_obj%L_pass%c(iblock-1)%zmtx(:,:) = cbuff(:,:)

      enddo

      ! now do right pass
      ! will use the labelling assuming solving for iblock=3
      do iblock = gf_obj%num_blocks-1,1,-1

         prev_norb = gf_obj%block_norb_size(iblock+1)
         next_norb = gf_obj%block_norb_size(iblock)

         ! reallocate buffers
         deallocate(dbuff,stat=ierr)
         call utils_dealloc_check('greenf_prepare_inverse','dbuff',ierr)
         allocate(dbuff(prev_norb,prev_norb),stat=ierr)
         call utils_alloc_check('greenf_prepare_inverse','dbuff',ierr)
         deallocate(cbuff,stat=ierr)
         call utils_dealloc_check('greenf_prepare_inverse','cbuff',ierr)
         allocate(cbuff(next_norb,prev_norb),stat=ierr)
         call utils_alloc_check('greenf_prepare_inverse','cbuff',ierr)

         ! get inverse of previous diagonal block
         ! dbuff = inv(d_44)
         dbuff(:,:) = gf_obj%R_pass%d(iblock+1)%zmtx(:,:)
         call linalg_invert_serial(dbuff,info)

         ! if inverse fails, return without setting pass_performed
         if (info /= 0) goto 999

         ! compute previous off-diagonal block
         ! cbuff = -a_34 * inv(d_4)
         call linalg_mat_mul_serial(cbuff, &                            ! out
              gf_obj%R_pass%c(iblock)%zmtx, dbuff, alpha=-cmplx_1)

         ! compute new diagonal block
         ! d_33 = a_33 - a_34 * inv(d_44) * a_43
         !      = a_33 + cbuff * transpose(a_34)
         ! must transpose upper diagonal
         call linalg_mat_mul_serial(gf_obj%R_pass%d(iblock)%zmtx, &     ! out
              cbuff, gf_obj%R_pass%c(iblock)%zmtx, beta=cmplx_1,opB='T')

         ! deposit cbuff = new a_34 (i.e. upper off-diagonal)
         gf_obj%R_pass%c(iblock)%zmtx(:,:) = cbuff(:,:)

      enddo

      ! set state flag
      gf_obj%pass_performed = .true.

  999 continue ! if error found

      ! destroy workspace
      deallocate(dbuff,stat=ierr)
      call utils_dealloc_check('greenf_prepare_inverse','dbuff',ierr)
      deallocate(cbuff,stat=ierr)
      call utils_dealloc_check('greenf_prepare_inverse','cbuff',ierr)

      call timer_clock('greenf_prepare_inverse',2)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving greenf_prepare_inverse'

      contains

      subroutine internal_build_inv_gf(inv_gf,h,s)
         !=============================================================!
         ! Constructs the inverse Green's function = E*S - H and       !
         ! deposits it in the tri-diagonal matrix inv_gf.              !
         !=============================================================!
         implicit none

         ! Arguments
         type(TRIMAT), intent(inout)  :: inv_gf
         type(TRIMAT), intent(in)     :: h
         type(TRIMAT), intent(in)     :: s

         ! internal
         integer :: iblock, nblocks

         nblocks = size(h%d)

         ! diagonal part
         do iblock = 1,nblocks
            inv_gf%d(iblock)%zmtx = cmplx_energy * s%d(iblock)%zmtx - &
                 h%d(iblock)%zmtx
         enddo

         ! upper off-diagonal
         do iblock = 1,nblocks-1
            if (inv_gf%upper) then
               inv_gf%c(iblock)%zmtx = &
                    cmplx_energy*s%c(iblock)%zmtx - h%c(iblock)%zmtx
            else
               inv_gf%c(iblock)%zmtx = transpose( &
                    cmplx_energy*s%c(iblock)%zmtx - h%c(iblock)%zmtx )
            endif
         enddo


      end subroutine internal_build_inv_gf

   end subroutine greenf_prepare_inverse

   !====================================================================!
   !====================================================================!

   subroutine greenf_dos(dos,gf_obj,info)
      !================================================================!
      ! Calculates the density of states of this Green's function      !
      ! using: DOS = - 1/pi * Im[ tr(G*S) ]                            !
      !                                                                !
      ! greenf_prepare_inverse must have been called successfully      !
      ! before this routine can be called.                             !
      !----------------------------------------------------------------!
      ! Arguments:                                                     !
      !   dos         (output): computed DOS                           !
      !   gf_obj      (input) : a sparse gf object                     !
      !   info        (output): zero if successful                     !
      !----------------------------------------------------------------!
      ! Necessary conditions:                                          !
      !   gf_obj%pass_performed = .true.                               !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: DP, cmplx_1, PI, stdout
      use linalg, only: linalg_invert_serial, linalg_mat_mul_serial
      use timer, only: timer_clock
      use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

      implicit none

      ! Arguments
      type(SPARGF), intent(in)   :: gf_obj
      real(kind=DP), intent(out) :: dos
      integer, intent(out)       :: info

      ! internal
      integer :: iblock
      integer :: iorb, jorb
      integer :: cols
      complex(kind=DP), allocatable :: cs(:,:)
      complex(kind=DP), allocatable :: g_ii(:,:)
      integer :: ierr

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering greenf_dos'

      call timer_clock('greenf_dos',1)

      call utils_assert(gf_obj%pass_performed,'Error in greenf_dos: &
           &greenf_prepare_inverse has not been called.')

      dos = 0.0_DP

      ! allocate temporary workspace
      ! will be reallocated as needed
      allocate(cs(1,1),stat=ierr)
      call utils_alloc_check('greenf_dos','cs',ierr)
      allocate(g_ii(1,1),stat=ierr)
      call utils_alloc_check('greenf_dos','g_ii',ierr)


      !---------------------------------------------
      ! compute contribution from each group
      ! will annotate as if for the 3rd block
      ! dos = tr[ g_33 * (s_33 + c^L_2*s_23 + c^R_4*s_43) ]

      do iblock = 1, gf_obj%num_blocks

         ! allocate required workspace
         cols = gf_obj%block_norb_size(iblock)
         call internal_reallocate(cs,cols,cols,'cs')
         call internal_reallocate(g_ii,cols,cols,'g_ii')


         ! compute inverse of diagonal block
         g_ii(:,:) = - gf_obj%e_cmplx * gf_obj%s%d(iblock)%zmtx(:,:) & ! - E*s
                     + gf_obj%h(gf_obj%ispin)%d(iblock)%zmtx(:,:)    & ! + h
                     + gf_obj%L_pass%d(iblock)%zmtx(:,:)             & ! + d^L
                     + gf_obj%R_pass%d(iblock)%zmtx(:,:)               ! + d^R
         call linalg_invert_serial(g_ii,info)

         ! exit early if failure
         if (info /= 0) then
            dos = 0.0_DP
            goto 999
         endif


         ! compute s_33 + c^L_2*s_23 + c^R_4*s_43
         ! additions of the products are done using beta=cmplx_1 in mat muls
         cs(:,:) = gf_obj%s%d(iblock)%zmtx(:,:)

         ! if first block, add c^L_2 * s_23
         if (iblock > 1) then
            call linalg_mat_mul_serial(cs,gf_obj%L_pass%c(iblock-1)%zmtx, &
                 gf_obj%s%c(iblock-1)%zmtx,beta=cmplx_1)
         endif

         ! if not final block, add c^R_4 * s_43
         ! s_32 = transposed as only upper diagonal of s stored
         if (iblock < gf_obj%num_blocks) then
            call linalg_mat_mul_serial(cs,gf_obj%R_pass%c(iblock)%zmtx, &
                 gf_obj%s%c(iblock)%zmtx,beta=cmplx_1,opB='T')
         endif


         ! get contribution to trace im[ tr( g_33 * cs ) ]
         do iorb=1,cols
         do jorb=1,cols
            dos = dos + aimag(g_ii(iorb,jorb)*cs(jorb,iorb))
         enddo
         enddo

      enddo

      ! multiply DOS by - 1.0 / pi
      dos = -dos / PI


 999  continue ! branch for inverse failure

      ! tidy up
      deallocate(cs,stat=ierr)
      call utils_dealloc_check('greenf_dos','cs',ierr)
      deallocate(g_ii,stat=ierr)
      call utils_dealloc_check('greenf_dos','g_ii',ierr)


      call timer_clock('greenf_dos',2)

      if (pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving greenf_dos'

      contains

      subroutine internal_reallocate(buff,cols,rows,name)
         !=============================================================!
         ! Checks that buff is of the correct size, and reallocates if !
         ! necessary.                                                  !
         !=============================================================!
         implicit none

         ! arguments
         complex(kind=DP), allocatable, intent(inout) :: buff(:,:)
         integer, intent(in) :: cols, rows
         character(len=*), intent(in) :: name

         ! internal
         integer :: ierr

         if( size(buff,dim=1) .ne. cols .or. size(buff,dim=2) .ne. rows) then
            deallocate(buff,stat=ierr)
            call utils_dealloc_check('greenf_dos','buff',ierr)
            allocate(buff(cols,rows),stat=ierr)
            call utils_alloc_check('greenf_dos','buff',ierr)
         endif

      end subroutine internal_reallocate

   end subroutine greenf_dos

   !====================================================================!
   !====================================================================!

   subroutine greenf_get_block(gf_block, gf_obj, col_block, row_block, info)
      !================================================================!
      ! Routine to compute a selected block of the Green's function,   !
      ! using equations 17--19 of Petersen et al. J. Comp. Phys 227    !
      ! (2008) 3174-3190.                                              !
      !                                                                !
      ! greenf_prepare_inverse must have been called successfully      !
      ! before this routine can be called.                             !
      !----------------------------------------------------------------!
      ! Arguments:                                                     !
      !   gf_block    (output): Requested Green's function block       !
      !   gf_obj      (input) : a sparse gf object                     !
      !   col_block   (input) : the block column of the required block !
      !   row_block   (input) : the block row of the required block    !
      !   info        (output): zero if successful                     !
      !----------------------------------------------------------------!
      ! Key Internal Variables:                                        !
      ! gc_prod, gc_prev : pointers to workspace for computing product !
      !    of many matrices m_1 * m_2 * m_3 * m_4 ...                  !
      !    Works by pointing gc_prev => m_1, and calculating           !
      !    gc_prod = gc_prev * m_2, then points gc_prev at gc_prod and !
      !    repeating the process recursively until all products done.  !
      !----------------------------------------------------------------!
      ! Necessary conditions:                                          !
      !   gf_obj%pass_performed = .true.                               !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: DP, cmplx_0, stdout
      use timer, only: timer_clock
      use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check
      use linalg, only: linalg_invert_serial, linalg_mat_mul_serial

      implicit none

      ! Arguments
      complex(kind=DP), intent(out) :: gf_block(:,:)
      type(SPARGF), intent(in)      :: gf_obj
      integer, intent(in)           :: col_block, row_block
      integer, intent(out)          :: info

      ! internal
      integer :: iblock
      integer :: col_norb, row_norb
      integer :: inorb
      integer :: ierr

      ! workspace for calculating the products, referenced with pointers
      complex(kind=DP), pointer :: buff1(:,:) => null() ! jd: Setting these to null() or else gfortran
      complex(kind=DP), pointer :: buff2(:,:) => null() !     complains at [*] 'buff1' is used uninitialized in this function
      complex(kind=DP), pointer :: gc_prev(:,:)  ! previous products
      complex(kind=DP), pointer :: gc_prod(:,:)  ! current products
      complex(kind=DP), pointer :: gc_swap(:,:)  ! temp for swapping

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Entering greenf_get_block'

      call timer_clock('greenf_get_block',1)

      ! init pointers
      gc_prev => null()
      gc_prod => null()
      gc_swap => null()

      ! check that the pass has been done
      call utils_assert(gf_obj%pass_performed,'Error in greenf_get_block: &
           &greenf_prepare_inverse has not been called.')

      ! check blocks are valid and assigned
      call utils_assert(col_block.le.gf_obj%num_blocks .and. col_block.ge.1, &
           'Error in greenf_get_block: col_block out of bounds', &
           col_block,gf_obj%num_blocks)
      call utils_assert(row_block.le.gf_obj%num_blocks .and. row_block.ge.1, &
           'Error in greenf_get_block: row_block out of bounds', &
           row_block,gf_obj%num_blocks)

      ! init result to zero
      gf_block(:,:) = cmplx_0


      ! size of diagonal block
      inorb = gf_obj%block_norb_size(col_block)

      ! point to initial buffers
      gc_prev => buff1 ! [*]
      gc_prod => buff2 ! [*]

      ! allocate workspace
      ! buffers will be reallocated as needed
      allocate(gc_prod(inorb,inorb),stat=ierr)
      call utils_alloc_check('greenf_get_block','=>gc_prod',ierr)
      allocate(gc_prev(1,1),stat=ierr)
      call utils_alloc_check('greenf_get_block','=>gc_prev',ierr)


      gc_prod(:,:) = cmplx_0
      gc_prev(:,:) = cmplx_0


      ! get diagonal block
      ! gf_ii = inv( -a_ii + d_ii^L + d_ii^R )
      ! stored in gc_prod
      gc_prod = -gf_obj%e_cmplx*gf_obj%s%d(col_block)%zmtx + &   ! - E*s_ii
                 gf_obj%h(gf_obj%ispin)%d(col_block)%zmtx + &    ! + h_ii
                 gf_obj%L_pass%d(col_block)%zmtx + &             ! + d_ii^L
                 gf_obj%R_pass%d(col_block)%zmtx                 ! + d_ii^R

      ! invert
      call linalg_invert_serial(gc_prod,info)

      ! exit early if inverse fails, returning info flag
      if (info/=0) goto 999


      ! compute gc_prod = g_ii * product of all off-diagonal matrices
      ! will comment with assuming col/row = 6,1
      if( col_block > row_block) then ! use left pass


         ! compute g_16 = { [ (gf_66 * c^L_5) * c^L_4] * ... }

         ! recursively do matrix products, with decreasing index
         do iblock=col_block-1,row_block,-1

            ! swap pointers
            gc_swap => gc_prev
            gc_prev => gc_prod
            gc_prod => gc_swap
            gc_swap => null()

            ! reallocate if necessary
            col_norb = gf_obj%block_norb_size(col_block)
            row_norb = gf_obj%block_norb_size(iblock)
            call internal_reallocate(gc_prod,col_norb,row_norb,'gc_prod')


            ! do product: gc_prod := gc_prev * c_iblock
            call linalg_mat_mul_serial(gc_prod, &
                 gc_prev, gf_obj%L_pass%c(iblock)%zmtx)

         enddo

      elseif( col_block < row_block ) then ! use right pass


         ! compute g_61 = { [ (gf_11 * c^R_2) * c^R_3] * ... }

         ! recursively do matrix products
         do iblock=col_block,row_block-1

            ! swap pointers
            gc_swap => gc_prev
            gc_prev => gc_prod
            gc_prod => gc_swap
            gc_swap => null()

            ! reallocate if necessary
            col_norb = gf_obj%block_norb_size(col_block)
            row_norb = gf_obj%block_norb_size(iblock+1)
            call internal_reallocate(gc_prod,col_norb,row_norb,'gc_prod')


            ! do product: gc_prod = gc_prev * c_iblock
            call linalg_mat_mul_serial(gc_prod, &
                 gc_prev, gf_obj%R_pass%c(iblock)%zmtx)

         enddo


      else ! col_block == row_block

         continue ! nothing required

      endif

      ! the required block is now in gc_prod, copy into gf_block
      call utils_assert(all(shape(gf_block) == shape(gc_prod)), &
           'Error in greenf_get_block: shape(gf_block) /= shape(gc_prod)', &
           size(gf_block,dim=1),size(gf_block,dim=2),&
           size(gc_prod,dim=1),size(gc_prod,dim=1))
      gf_block(:,:) = gc_prod(:,:)

 999  continue ! if info /= 0, make sure deallocation is done

      ! tidy up
      deallocate(gc_prod,stat=ierr)
      call utils_dealloc_check('greenf_get_block','=>gc_prod',ierr)
      deallocate(gc_prev,stat=ierr)
      call utils_dealloc_check('greenf_get_block','=>gc_prev',ierr)
      nullify(gc_swap)
      nullify(gc_prev)
      nullify(gc_prod)


      call timer_clock('greenf_get_block',2)

      if (pub_debug_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving greenf_get_block'

      contains

      subroutine internal_reallocate(buff,cols,rows,name)
         !=============================================================!
         ! Checks that buff is of the correct size, and reallocates if !
         ! necessary.                                                  !
         !=============================================================!

         implicit none

         ! arguments
         complex(kind=DP), pointer, intent(inout) :: buff(:,:)
         integer, intent(in) :: cols, rows
         character(len=*), intent(in) :: name

         ! internal
         integer :: ierr

         if( size(buff,dim=1) .ne. cols .or. size(buff,dim=2) .ne. rows) then
            deallocate(buff,stat=ierr)
            call utils_dealloc_check('greenf_get_block','=>buff',ierr)
            allocate(buff(cols,rows),stat=ierr)
            call utils_alloc_check('greenf_get_block','=>buff',ierr)
         endif

      end subroutine internal_reallocate

   end subroutine greenf_get_block

   !====================================================================!
   !====================================================================!

   subroutine greenf_get_atom_block(gf_block, gf_obj, col_start, col_end, &
                                    row_start, row_end, info)
      !================================================================!
      ! Returns the Green's function block for column atoms between    !
      ! col_start, col_end, and row atoms between row_start, row_end.  !
      !                                                                !
      ! The atoms are indexed by their original ordering in the input  !
      ! file, and the range gives all atoms between those two values,  !
      ! as found in the original order.                                !
      !                                                                !
      ! For efficiency reasons, all column atoms must be contained     !
      ! within a tridiagonal block (similarly for row atoms), and all  !
      ! atoms must be sequential in this internal odering. This can be !
      ! enforced in greenf_init_blocking using the optional            !
      ! force_groups argument.                                         !
      !                                                                !
      ! greenf_prepare_inverse must have been called successfully      !
      ! before this routine can be called.                             !
      !----------------------------------------------------------------!
      ! Arguments:                                                     !
      !   gf_block    (output): Requested Green's function block       !
      !   gf_obj      (input) : a sparse gf object                     !
      !   col_start   (input) : the first column atom (original order) !
      !   col_end     (input) : the last column atom  (original order) !
      !   row_start   (input) : the first row atom    (original order) !
      !   row_end     (input) : the last row atom     (original order) !
      !   info        (output): zero if successful                     !
      !----------------------------------------------------------------!
      ! Necessary conditions:                                          !
      !   gf_obj%pass_performed = .true.                               !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: DP, CRLF, cmplx_0
      use parallel_strategy, only: par=>pub_par
      use rundat, only: pub_debug
      use utils, only: utils_assert, utils_alloc_check, utils_dealloc_check
      use linalg, only: linalg_invert_serial, linalg_mat_mul_serial

      implicit none

      ! Arguments
      complex(kind=DP), intent(out) :: gf_block(:,:)
      type(SPARGF), intent(in)      :: gf_obj
      integer, intent(in)           :: col_start, col_end, row_start, row_end
      integer, intent(out)          :: info

      ! internal
      integer :: iblock, col_block, row_block
      integer :: col_norb, row_norb
      integer :: istart, iend, jstart, jend, inorb, jnorb
      integer :: iat, iat_distr, ndx, idx
      integer :: ierr
      complex(kind=DP), allocatable :: full_gf_block(:,:)


      if(pub_debug) then
         ! check that all row/col atoms are in the same block
         do iat=col_start,col_end
            col_block = gf_obj%block_of_atom(par%distr_atom(col_start))
            iat_distr = par%distr_atom(iat)
            iblock = gf_obj%block_of_atom(iat_distr)
            call utils_assert( col_block == iblock, 'Error in greenf_get_atom_&
              &block: column atoms are not in same block.'//CRLF//'A block must &
              &be forced using force_groups when calling greenf_init_blocking.', &
              col_start, col_end, col_block, iblock)
         enddo
         do iat=row_start,row_end
            row_block = gf_obj%block_of_atom(par%distr_atom(row_start))
            iat_distr = par%distr_atom(iat)
            iblock = gf_obj%block_of_atom(iat_distr)
            call utils_assert( row_block == iblock, 'Error in greenf_get_atom_&
              &block: row atoms are not in same block.'//CRLF//'A block must &
              &be forced using force_groups when calling greenf_init_blocking.', &
              row_start, row_end, row_block, iblock)
         enddo

         ! check that row/col atoms are sequential in the block
         ndx = gf_obj%index_of_atom(par%distr_atom(col_start))
         do idx=0,col_end - col_start
            iat_distr = par%distr_atom(col_start+idx)
            call utils_assert(iat_distr == gf_obj%atom_order(ndx+idx), &
                'Error in greenf_get_atom_block: column atoms must be &
                &sequential in internal ordering.')
         enddo
         ndx = gf_obj%index_of_atom(par%distr_atom(row_start))
         do idx=0,row_end - row_start
            iat_distr = par%distr_atom(row_start+idx)
            call utils_assert(iat_distr == gf_obj%atom_order(ndx+idx), &
                'Error in greenf_get_atom_block: row atoms must be &
                &sequential in internal ordering.')
         enddo
      end if


      ! find the row/column block to calculate
      col_block = gf_obj%block_of_atom(par%distr_atom(col_start))
      row_block = gf_obj%block_of_atom(par%distr_atom(row_start))

      ! check blocks are valid and assigned
      call utils_assert(col_block.le.gf_obj%num_blocks .and. col_block.ge.1, &
           'Error in greenf_get_atom_block: col_block out of bounds', &
           col_block,gf_obj%num_blocks)
      call utils_assert(row_block.le.gf_obj%num_blocks .and. row_block.ge.1, &
           'Error in greenf_get_atom_block: row_block out of bounds', &
           row_block,gf_obj%num_blocks)


      ! allocate workspace
      col_norb = gf_obj%block_norb_size(col_block)
      row_norb = gf_obj%block_norb_size(row_block)
      allocate(full_gf_block(col_norb,row_norb),stat=ierr)
      call utils_alloc_check('greenf_get_atom_block','full_gf_block',ierr)

      ! get entire block
      call greenf_get_block(full_gf_block,gf_obj,col_block,row_block,info)

      gf_block = cmplx_0

      if (info == 0) then

         ! get required sub-block of this gf block
         ! (uses the fact that requested atoms are sequential in the block)
         istart = gf_obj%atom_ngwf_ndx(par%distr_atom(col_start))
         jstart = gf_obj%atom_ngwf_ndx(par%distr_atom(row_start))
         iend   = gf_obj%atom_ngwf_ndx(par%distr_atom(col_end)) + &
                  gf_obj%ngwfs_on_atom(par%distr_atom(col_end)) - 1
         jend   = gf_obj%atom_ngwf_ndx(par%distr_atom(row_end)) + &
                  gf_obj%ngwfs_on_atom(par%distr_atom(row_end)) - 1
         inorb  = iend - istart + 1
         jnorb  = jend - jstart + 1
         col_norb = size(gf_block,dim=1)
         row_norb = size(gf_block,dim=2)
         call utils_assert(inorb==col_norb .and. jnorb==row_norb, 'Error in &
              &greenf_get_atom_block: argument gf_block incorrect size', &
              col_norb, row_norb, inorb, jnorb)

         gf_block(:,:) = full_gf_block(istart:iend,jstart:jend)
      endif

      ! tidy up
      deallocate(full_gf_block,stat=ierr)
      call utils_dealloc_check('greenf_get_atom_block','full_gf_block',ierr)

   end subroutine greenf_get_atom_block

   !====================================================================!
   !====================================================================!

   integer(kind=LONG) function greenf_memory_estimate(gf_obj)
      !================================================================!
      ! Returns an estimate of the memory required to store the gf_obj !
      ! data structures on each proc.                                  !
      ! Ignores the information describing the blocking scheme which   !
      ! contributes only a small part of the memory used.              !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: LONG, cmplx_size

      implicit none

      ! Arguments
      type(SPARGF), intent(in) :: gf_obj

      ! internal
      integer :: iblock
      integer(kind=LONG) :: norb1, norb2

      integer, parameter :: num_work_mats = 2
      integer :: num_trimats
      integer(kind=LONG), parameter :: cmplx_size_long=int(cmplx_size,kind=LONG)


      greenf_memory_estimate = 0_LONG

      ! return 0 if no blocking performed
      if (.not.gf_obj%is_blocked) return

      ! compute memory of tri-diagonal matrices
      ! diagonal storage
      do iblock = 1, gf_obj%num_blocks
         norb1 = int(gf_obj%block_norb_size(iblock),kind=LONG)
         greenf_memory_estimate = greenf_memory_estimate + norb1*norb1
      enddo

      ! off-diagonal storage
      do iblock = 1, gf_obj%num_blocks-1
         norb1 = int(gf_obj%block_norb_size(iblock),kind=LONG)
         norb2 = int(gf_obj%block_norb_size(iblock+1),kind=LONG)
         greenf_memory_estimate = greenf_memory_estimate + norb1*norb2
      enddo

      ! multiply by the number of trimats
      num_trimats = 3 + gf_obj%num_spins ! s, L_pass, R_pass, h(spins)
      greenf_memory_estimate = greenf_memory_estimate * num_trimats

      ! include memory from largest temporary array = largest diagonal block
      norb1 = int(maxval(gf_obj%block_norb_size),kind=LONG)
      greenf_memory_estimate = greenf_memory_estimate + &
           norb1*norb1*num_work_mats

      ! all are complex double precision
      greenf_memory_estimate = greenf_memory_estimate*cmplx_size_long

      return

   end function greenf_memory_estimate

   !====================================================================!
   !====================================================================!

   subroutine trimat_reduce(tmat)
      !================================================================!
      ! Routine that sums a tri-diagonal matrix across all procs.      !
      ! Afterwards, the tri-diagonal matrix is the same on all procs.  !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use comms, only: comms_reduce

      implicit none

      ! Arguments
      type(TRIMAT), intent(inout)  :: tmat

      ! internal
      integer :: iblock

      ! perform comms_reduce on all blocks
      do iblock=1,size(tmat%d)
         call comms_reduce('SUM',tmat%d(iblock)%zmtx)
      enddo
      do iblock=1,size(tmat%c)
         call comms_reduce('SUM',tmat%c(iblock)%zmtx)
      enddo

   end subroutine trimat_reduce

   !====================================================================!
   !====================================================================!

   subroutine trimat_deposit_atom_block(block,tmat,gf_obj,col_start, &
                                        col_end,row_start,row_end)
      !================================================================!
      ! Subroutine to deposit matrix elements into the tri-diagonal    !
      ! matrix tmat. The block is between start/end column and row     !
      ! atoms (in original order), and must be the correct shape/size. !
      ! The block must be contained within the tri-diagonal sparsity   !
      ! pattern, and the column atoms must all be in the same tri-     !
      ! diagonal block, similarly for row atoms.                       !
      !----------------------------------------------------------------!
      ! Arguments:                                                     !
      !   block       (output): Block to deposit into tmat             !
      !   tmat        (input) : the tri-digonal matrix to deposit into !
      !   col_start   (input) : the first column atom (original order) !
      !   col_end     (input) : the last column atom  (original order) !
      !   row_start   (input) : the first row atom    (original order) !
      !   row_end     (input) : the last row atom     (original order) !
      !----------------------------------------------------------------!
      ! Necessary conditions:                                          !
      !   the tmat must be allocated                                   !
      !   column atoms must be in the same tri-diagonal block          !
      !   row atoms must be in the same tri-diagonal block             !
      !   the column/row blocks must be adjacent, maintaining the tri- !
      !   diagonal sparsity pattern.                                   !
      !                                                                !
      !   ***This routine assumes that tmat is Hermitian!***           !
      !                                                                !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: DP, CRLF
      use parallel_strategy, only: par=>pub_par
      use rundat, only: pub_debug
      use utils, only: utils_assert

      implicit none

      ! Arguments
      complex(kind=DP), intent(in) :: block(:,:)
      type(SPARGF), intent(in)     :: gf_obj
      type(TRIMAT), intent(inout)  :: tmat
      integer, intent(in)          :: col_start,col_end,row_start,row_end

      ! internal
      integer :: col_block, row_block, iblock
      integer :: istart, iend, jstart, jend
      integer :: iat, iat_distr, idx, ndx
      integer :: inorb, jnorb, col_norb, row_norb
      logical :: do_ctranspose

      if(pub_debug) then
         ! check that all row/col atoms are in the same block
         do iat=col_start,col_end
            col_block = gf_obj%block_of_atom(par%distr_atom(col_start))
            iat_distr = par%distr_atom(iat)
            iblock = gf_obj%block_of_atom(iat_distr)
            call utils_assert( col_block == iblock, 'Error in trimat_deposit_atom_&
              &block: column atoms are not in same block.'//CRLF//'A block must &
              &be forced using force_groups when calling greenf_init_blocking.', &
              col_start, col_end, col_block, iblock)
         enddo
         do iat=row_start,row_end
            row_block = gf_obj%block_of_atom(par%distr_atom(row_start))
            iat_distr = par%distr_atom(iat)
            iblock = gf_obj%block_of_atom(iat_distr)
            call utils_assert( row_block == iblock, 'Error in trimat_deposit_atom_&
              &block: row atoms are not in same block.'//CRLF//'A block must &
              &be forced using force_groups when calling greenf_init_blocking.', &
              row_start, row_end, row_block, iblock)
         enddo

         ! check that row/col atoms are sequential in the block
         ndx = gf_obj%index_of_atom(par%distr_atom(col_start))
         do idx=0,col_end - col_start
            iat_distr = par%distr_atom(col_start+idx)
            call utils_assert(iat_distr == gf_obj%atom_order(ndx+idx), &
                'Error in greenf_get_atom_block: column atoms must be &
                &sequential in internal ordering.')
         enddo
         ndx = gf_obj%index_of_atom(par%distr_atom(row_start))
         do idx=0,row_end - row_start
            iat_distr = par%distr_atom(row_start+idx)
            call utils_assert(iat_distr == gf_obj%atom_order(ndx+idx), &
                'Error in greenf_get_atom_block: row atoms must be &
                &sequential in internal ordering.')
         enddo
      end if

      ! make sure that the block is within the tridiagonal scheme
      col_block = gf_obj%block_of_atom(par%distr_atom(col_start))
      row_block = gf_obj%block_of_atom(par%distr_atom(row_start))
      call utils_assert(abs(col_block-row_block) .le. 1, 'Error in trimat_&
           &deposit_atom_block: requested block is outside the tridiagonal&
           & scheme.', col_block, row_block)

      ! find the indices of where to deposit
      istart = gf_obj%atom_ngwf_ndx(par%distr_atom(col_start))
      jstart = gf_obj%atom_ngwf_ndx(par%distr_atom(row_start))
      iend   = gf_obj%atom_ngwf_ndx(par%distr_atom(col_end)) + &
               gf_obj%ngwfs_on_atom(par%distr_atom(col_end)) - 1
      jend   = gf_obj%atom_ngwf_ndx(par%distr_atom(row_end)) + &
               gf_obj%ngwfs_on_atom(par%distr_atom(row_end)) - 1
      inorb  = iend - istart + 1
      jnorb  = jend - jstart + 1
      col_norb = size(block,dim=1)
      row_norb = size(block,dim=2)
      call utils_assert(inorb==col_norb .and. jnorb==row_norb, 'Error in &
           &trimat_deposit_atom_block: argument block incorrect size', &
           col_norb, row_norb, inorb, jnorb)


      ! deposit
      if (col_block == row_block) then
         ! deposit in diagonal block
         tmat%d(col_block)%zmtx(istart:iend,jstart:jend) = block(:,:)
      else
         ! must transpose if this is an lower diagonal block, and tmat%upper
         ! or this is upper diagonal block, and .not. tmat%upper
         do_ctranspose = ((row_block < col_block) .and. tmat%upper) .or. &
                         ((row_block > col_block) .and. .not.tmat%upper)

         ! tmat is assumed to be Hermitian here
         if (do_ctranspose) then
            tmat%c(row_block)%zmtx(jstart:jend,istart:iend) = &
                 conjg(transpose(block(:,:)))
         else
            tmat%c(col_block)%zmtx(istart:iend,jstart:jend) = block(:,:)
         endif

      endif

   end subroutine trimat_deposit_atom_block

   !====================================================================!
   !====================================================================!

   subroutine trimat_scale(tmat,alpha)
      !================================================================!
      ! Scales a tri-diagonal matrix by the complex constant alpha     !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, Nov 2014                               !
      !================================================================!

      use constants, only: DP

      implicit none

      ! Arguments
      type(TRIMAT), intent(inout)  :: tmat
      complex(kind=DP), intent(in) :: alpha

      ! internal
      integer :: iblock, nblocks

      nblocks = size(tmat%d)

      ! diagonal part
      do iblock = 1,nblocks
         tmat%d(iblock)%zmtx = alpha * tmat%d(iblock)%zmtx
      enddo

      ! upper off-diagonal
      do iblock = 1,nblocks-1
         tmat%c(iblock)%zmtx = alpha*tmat%c(iblock)%zmtx
      enddo


   end subroutine trimat_scale

   !====================================================================!
   !====================================================================!

   subroutine greenf_print_full_inverse(gf_obj)
      !================================================================!
      ! This routine computes the full Green's function by calculating !
      ! every block. This is provided for debug purposes only and      !
      ! should not be used in production code as it is slower, and uses!
      ! more memory, than the conventional matrix inverse.             !
      !                                                                !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: DP, cmplx_0, stdout
      use utils, only: utils_unit, utils_alloc_check, utils_dealloc_check

      implicit none

      ! Arguments
      type(SPARGF), intent(in) :: gf_obj

      ! internal
      complex(kind=DP) :: mat(gf_obj%norb,gf_obj%norb)
      complex(kind=DP), allocatable :: zbuff(:,:)
      integer :: istart, iend, jstart, jend, iorb, jorb
      integer :: iblock, jblock
      integer :: info
      integer :: iunit
      integer :: ierr

      ! init
      mat(:,:) = cmplx_0

      ! recursively get each block, and deposit in matrix
      istart = 1
      do iblock = 1,gf_obj%num_blocks

         iorb = gf_obj%block_norb_size(iblock)
         iend = istart + iorb - 1

         jstart = 1
         do jblock = 1,gf_obj%num_blocks

            jorb = gf_obj%block_norb_size(jblock)
            jend = jstart + jorb - 1

            write(stdout,'(a,2i3,a,4i5)') 'Full inverse: ', iblock, jblock, &
                 '|', istart, iend, jstart, jend

            allocate(zbuff(iorb,jorb),stat=ierr)
            call utils_alloc_check('greenf_print_full_inverse','zbuff',ierr)

            call greenf_get_block(zbuff,gf_obj,iblock,jblock,info)

            mat(istart:iend,jstart:jend) = zbuff(:,:)
            deallocate(zbuff,stat=ierr)
            call utils_dealloc_check('greenf_print_full_inverse','zbuff',ierr)

            jstart = jstart + jorb
         enddo
         istart = istart + iorb
      enddo

      ! write out to file
      iunit = utils_unit()
      open(unit=iunit,file='full_inv.fbin',form='UNFORMATTED')
      write(iunit) gf_obj%e_cmplx
      write(iunit) gf_obj%norb
      write(iunit) mat(:,:)
      close(iunit)


   end subroutine greenf_print_full_inverse

   !====================================================================!
   !====================================================================!

   subroutine tmat_print_matrix(tmat,filename)
      !================================================================!
      ! This routine prints a tri-diagonal matrix of TRIMAT type for   !
      ! debugging purposes. The output file is in:                     !
      !     row_ndx, col_ndx, real(el), imag(el)                       !
      ! format, where el = tmat(col_ndx,row_ndx).                      !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!
      use constants, only: DP
      use utils, only: utils_unit

      implicit none

      ! Arguments
      type(TRIMAT), intent(in)     :: tmat
      character(len=*), intent(in) :: filename

      ! internal
      integer :: iunit
      integer :: iorb, jorb, offset, iblock, borb
      complex(kind=DP) :: el

      iunit = utils_unit()

      open(unit=iunit,file=trim(filename))

      ! write diagonal block
      offset = 0
      do iblock = 1,size(tmat%d)
         borb = size(tmat%d(iblock)%zmtx,dim=1)
         do iorb = 1, size(tmat%d(iblock)%zmtx,dim=1)
         do jorb = 1, size(tmat%d(iblock)%zmtx,dim=2)
            el = tmat%d(iblock)%zmtx(iorb,jorb)
            write(iunit,'(2i6,2f20.12)') iorb+offset, jorb+offset, &
                 real(el,kind=DP), aimag(el)
         enddo
         enddo
         offset = offset + borb
      enddo

      ! write upper off-diagonal block
      offset = 0
      do iblock = 1,size(tmat%c)
         borb = size(tmat%c(iblock)%zmtx,dim=1)
         do iorb = 1, size(tmat%c(iblock)%zmtx,dim=1)
         do jorb = 1, size(tmat%c(iblock)%zmtx,dim=2)
            el = tmat%c(iblock)%zmtx(iorb,jorb)
            write(iunit,'(2i6,2f20.12)') iorb+offset, jorb+offset+borb, &
                 real(el,kind=DP), aimag(el)
            write(iunit,'(2i6,2f20.12)') jorb+offset+borb, iorb+offset, &
                 real(el,kind=DP), aimag(el)
         enddo
         enddo
         offset = offset + borb
      enddo

      close(iunit)

   end subroutine tmat_print_matrix

   !====================================================================!
   !====================================================================!

   subroutine greenf_print_ham(gf_obj,filename)
      !================================================================!
      ! This routine prints a tri-diagonal matrix of TRIMAT type for   !
      ! debugging purposes. The output is unformatted, in the original !
      ! atom order.                                                    !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: DP, cmplx_0
      use parallel_strategy, only: par=>pub_par
      use utils, only: utils_unit

      implicit none

      ! Arguments
      type(SPARGF), intent(in)     :: gf_obj
      character(len=*), intent(in) :: filename

      ! internal
      integer :: iunit
      integer :: iat, jat
      integer :: idx, jdx, iorb, jorb
      integer :: iat_distr, jat_distr
      complex(kind=DP), target  :: mat(gf_obj%norb,gf_obj%norb)
      complex(kind=DP), pointer :: ptr(:,:)

      mat = cmplx_0
      ptr => null()

      idx = 1

      ! loop over all atoms, in internal order
      do iat = 1,par%nat
         iat_distr = par%distr_atom(iat)
         if (gf_obj%index_of_atom(iat_distr) == IGNORED) cycle

         iorb = gf_obj%ngwfs_on_atom(iat_distr)

         jdx = 1
         do jat = 1,par%nat
            jat_distr = par%distr_atom(jat)
            if (gf_obj%index_of_atom(jat_distr) == IGNORED) cycle

            jorb = gf_obj%ngwfs_on_atom(jat_distr)

            ptr => mat(idx:idx+iorb-1,jdx:jdx+jorb-1)
            call internal_get_atom_block(ptr,gf_obj%h(gf_obj%ispin),iat_distr,jat_distr)

            jdx = jdx + jorb
         enddo
         idx = idx + iorb
      enddo


      nullify(ptr)

      iunit = utils_unit()
      open(unit=iunit,file=trim(filename),form='unformatted')
      write(iunit) mat(:,:)
      close(iunit)

      return

      contains
         subroutine internal_get_atom_block(block,tmat,iat,jat)

            use constants, only: DP, cmplx_0
            use utils, only: utils_assert

            implicit none

            ! arguments
            type(TRIMAT), intent(in)        :: tmat
            complex(kind=DP), intent(inout) :: block(:,:)
            integer, intent(in)             :: iat, jat

            ! internal
            integer :: iblock, jblock
            integer :: iostart, iostop, jostart, jostop
            integer :: iorb, jorb

            ! set block to 0
            block = cmplx_0

            ! find blocks of atoms
            iblock = gf_obj%block_of_atom(iat)
            jblock = gf_obj%block_of_atom(jat)

            ! find orbital indices of atoms
            iorb    = gf_obj%ngwfs_on_atom(iat)
            iostart = gf_obj%atom_ngwf_ndx(iat)
            iostop  = iostart + iorb - 1

            jostart = gf_obj%atom_ngwf_ndx(jat)
            jorb    = gf_obj%ngwfs_on_atom(jat)
            jostop  = jostart + jorb - 1


            call utils_assert(all(shape(block) == (/iorb,jorb/)),'Error in &
                 &internal_get_atom_block: block incorrect shape',iat,jat)

            if (abs(iblock - jblock) .gt. 1) then
               ! exit if outside tri-diagonal partitioning
               return

            else if (iblock == jblock) then
               ! diagonal block
               block(:,:) = tmat%d(iblock)%zmtx(iostart:iostop,jostart:jostop)

            else if (iblock == jblock - 1) then
               ! upper diagonal block
               block(:,:) = tmat%c(iblock)%zmtx(iostart:iostop,jostart:jostop)

            else
               ! lower diagonal block
               block(:,:) = conjg(transpose( &
                    tmat%c(jblock)%zmtx(jostart:jostop,iostart:iostop) ))
            endif

         end subroutine internal_get_atom_block


   end subroutine greenf_print_ham

   !====================================================================!
   !====================================================================!

   subroutine greenf_print_gf_obj(gf_obj,filename)
      !================================================================!
      ! This routine prints the tri-diagonal blocking structure of the !
      ! Green's function for debug purposes.                           !
      ! If 'filename' absent, information will be written to stdout.   !
      !----------------------------------------------------------------!
      ! Written by Robert Bell, May 2014                               !
      !================================================================!

      use constants, only: stdout
      use utils, only: utils_unit

      implicit none

      ! arguments
      type(SPARGF), intent(in)               :: gf_obj
      character(len=*), intent(in), optional :: filename

      ! internal
      integer :: iunit
      integer :: ii

      iunit = stdout
      if (present(filename)) then
         iunit = utils_unit()
         open(unit=iunit,file=trim(filename))
      endif

      write(iunit,'(a,i8)') 'num_spins       ', gf_obj%num_spins
      write(iunit,'(a,a)')  'ham_type        ', gf_obj%ham_type
      write(iunit,'(a,i8)') 'nat             ', gf_obj%nat
      write(iunit,'(a,i8)') 'norb            ', gf_obj%norb
      write(iunit,'(a,i8)') 'ispin           ', gf_obj%ispin
      write(iunit,'(a,l8)') 'is_blocked      ', gf_obj%is_blocked
      write(iunit,'(a,l8)') 'is_alloc        ', gf_obj%is_alloc
      write(iunit,'(a,l8)') 'is_init         ', gf_obj%is_init
      write(iunit,'(a,l8)') 'pass_performed  ', gf_obj%pass_performed

      write(iunit,'(a,i8)') 'num_blocks      ', gf_obj%num_blocks
      write(iunit,'(a8,2a18)') 'iblock', 'block_atom_size', 'block_norb_size'
      do ii = 1,gf_obj%num_blocks
         write(iunit,'(i8,2i18)') ii, gf_obj%block_atom_size(ii), &
              gf_obj%block_norb_size(ii)
      enddo

      write(iunit,'(a8,5a18)') 'iat', 'atom_order', 'index_of_atom', &
           'block_of_atom', 'atom_ngwf_ndx', 'ngwfs_on_atom'
      do ii = 1, size(gf_obj%atom_order)
         write(iunit,'(i8,5i18)') ii, gf_obj%atom_order(ii), &
                                      gf_obj%index_of_atom(ii), &
                                      gf_obj%block_of_atom(ii), &
                                      gf_obj%atom_ngwf_ndx(ii), &
                                      gf_obj%ngwfs_on_atom(ii)
      enddo

      if (present(filename)) close(iunit)

   end subroutine greenf_print_gf_obj

   !====================================================================!
   !====================================================================!

end module greenf

