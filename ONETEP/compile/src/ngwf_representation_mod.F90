! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-   !
!==================================================================!
!                                                                  !
!                  NGWF Representation module                      !
!                                                                  !
! This module contains routines associated with the use of the     !
! NGWF_REP and NGWF_HAM types. NGWF_REP is a container to simplify !
! the use of a given NGWF representation, its associated matrices  !
! (overlap, inverse overlap, projector overlap, and the parts of   !
! the Hamiltonian that do not depend on the density). NGWF_HAM is  !
! a container to store the matrices associated with the            !
! Hamiltonian in a given NGWF representation.                      !
!------------------------------------------------------------------!
! This module was created by Nicholas Hine in May 2010.            !
! Modified for embedding structures by Robert Charlton, Oct 2018.  !
!==================================================================!

module ngwf_representation

   use datatypes, only: FUNCTIONS
   use constants, only: DP, PI, SQRT_PI, NORMAL, VERBOSE, stdout, REP_N_SWEXES
   use hash_table, only: HT_HASH_TABLE
   use sparse_embed, only: SPAM3_EMBED
   ! agreco: add dependence on sparse_array
   !use sparse_array, only: SPAM3_ARRAY
   use sw_expansion_type, only: SW_EX

   implicit none

   private

   type NGWF_REP

      ! The NGWFs themselves, in PPDs on the standard grid
      !real(kind=DP), allocatable     :: ngwfs_on_grid(:)
      ! agrecocmplx: use new FUNCTIONS type
      ! rc2013: each region has its own ngwfs_on_grid
      type(FUNCTIONS), allocatable :: ngwfs_on_grid(:)

      ! Overlap matrix (augmented, in PAW or USP formalisms)
      type(SPAM3_EMBED) :: overlap
      ! Original overlap <phi_a|phi_b> before augmentation
      type(SPAM3_EMBED) :: ngwf_overlap
      ! Overlap between conduction NGWFs and valence NGWFs
      type(SPAM3_EMBED) :: cross_overlap, cross_overlap_tr
      ! Original versions before augmentation
      type(SPAM3_EMBED) :: ngwf_cross_overlap, ngwf_cross_overlap_tr
      ! Inverse overlap matrix (and whether it is initialised)
      type(SPAM3_EMBED) :: inv_overlap
      logical :: inv_overlap_init
      ! NGWF-projector overlap matrix
      type(SPAM3_EMBED) :: sp_overlap

      ! mjsp: Fragment matrices: for EDA and/or SCF-MI calculations,
      ! mjsp: 'overlap' holds the fragmented overlap matrix, and
      ! mjsp: 'overlap_scfmi_full' holds the full supermolecule overlap matrix
      ! Whether the fragment overlap matrix (%overlap) is initialised
      logical           :: inv_overlap_b_init
      ! Supermolecule (full) overlap matrix
      type(SPAM3_EMBED) :: overlap_scfmi_full
      ! Supermolecule (full) inverse overlap matrix
      type(SPAM3_EMBED) :: inv_overlap_scfmi_full

      ! NGWF-Hubbard projector overlap matrix
      type(SPAM3_EMBED) :: hub_overlap, hub_overlap_t, hub_overlap_t_bare

      ! ddor: Hubbard projector - PAW projector overlaps for PAW+U
      type(SPAM3_EMBED) :: hub_proj_paw_overlap, hub_ngwf_paw_overlap

      ! Occupation number of each spin channel and k-point.
      integer, allocatable :: n_occ(:,:)

      ! Parts of Hamiltonian indep of density kernel
      type(SPAM3_EMBED) :: kinet
      type(SPAM3_EMBED) :: nonlocpot ! only for nonlocal NCPPs

      ! rc2013: block-orthogonalised projector matrices
      type(SPAM3_EMBED) :: bo_projector, bo_projector_t

      ! jd: Spherical wave expansion of this NGWF rep, for DMA and HFx
      type(SW_EX) :: swexes(REP_N_SWEXES)

      ! jd: Handle to cache of remote NGWFs (currently only used with SWx)
      integer     :: ngwf_cache_handle

      ! jd: Flag indicating if the NGWF_REP should be invalidated when NGWFs
      !     change. This is normally .true., but will be .false. for a rep
      !     corresponding to a fixed state (e.g. in-vacuum NGWFs).
      logical :: invalidatable = .true.

      ! Postfix string for the structure code of matrices of this rep
      character(len=9) :: postfix

      ! rc2013: number of subsystems in this rep
      integer          :: nsub

      ! jd: Reflect any changes to this type in ngwf_rep_register_change()
      ! ******************************************************************

   end type NGWF_REP

   type NGWF_HAM

      ! Components of the Hamiltonian dependent on the density kernel
      type(SPAM3_EMBED), allocatable :: nonlocpot(:) ! only in PAW/USP
      !type(SPAM3_ARRAY) :: nonlocpot
      type(SPAM3_EMBED), allocatable :: hubbard_ham(:)
      !type(SPAM3_ARRAY) :: hubbard_ham
      type(SPAM3_EMBED), allocatable :: hfexchange(:)
      !type(SPAM3_ARRAY) :: hfexchange
      type(SPAM3_EMBED), allocatable :: polemb(:) ! jd: only in polarisable embedding
                                            ! 1:induced, 2:permanent. No spins.
      type(SPAM3_EMBED), allocatable :: vacuum_polemb(:) ! jd: only in polarisable embedding with vacuum QM*
                                            ! 1:induced, 2:permanent. No spins.
      type(SPAM3_EMBED), allocatable :: vacuum_aux_polemb(:) ! jd: only in polarisable embedding with vacuum QM*
                                            ! 1:induced, 2:permanent. No spins.
      type(SPAM3_EMBED), allocatable :: lhxc(:)
      type(SPAM3_EMBED), allocatable :: scissor(:) ! ny: species dependent scissor shift
      type(SPAM3_EMBED), allocatable :: dfdtau(:) ! JCW: matrix element over dfdtau
                                            ! JCW: for use in FDO XC potential
      type(SPAM3_EMBED), allocatable :: confinement(:) ! gcc32: for use in
                                               ! confining the NGWFS
      !type(SPAM3_ARRAY) :: confinement
      type(SPAM3_EMBED), allocatable :: ham(:)
      !type(SPAM3_ARRAY) :: ham
      type(SPAM3_EMBED), allocatable :: cond_non_proj_ham(:) ! only in conduction
      !type(SPAM3_ARRAY) :: cond_non_proj_ham
      type(SPAM3_EMBED), allocatable :: dijhat(:)
      !type(SPAM3_ARRAY) :: dijhat
      type(SPAM3_EMBED), allocatable :: proj_embed_ham(:) ! rc2013: for embedding
      type(SPAM3_EMBED), allocatable :: unproj_ham(:) ! rc2013: for embedding
      real(kind=DP) :: cond_shift ! only used in conduction

      ! ******************************************************************
      ! jd: Reflect any changes to this type in ngwf_ham_register_change()
      ! ******************************************************************

   end type NGWF_HAM

   public :: NGWF_REP
   public :: NGWF_HAM

   public :: ngwf_rep_create
   public :: ngwf_rep_destroy
   public :: ngwf_rep_register_change
   public :: ngwf_rep_build_subsystem

   public :: ngwf_ham_create
   public :: ngwf_ham_destroy
   public :: ngwf_ham_copy
   public :: ngwf_ham_register_change
   public :: ngwf_ham_build_subsystem

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_rep_create(rep, postfix, mdl, is_cmplx)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the NGWF_REP type,   !
    ! apart from the density-dependent matrices of the Hamiltonian (which are !
    ! done separately to save memory).                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rep (inout) : Container for NGWFs and matrices to be allocated.       !
    !   postfix (input)   : Label for this rep.                               !
    !   mdl (input)       : Holds info relating to entire system.             !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 07/10/2010                                  !
    ! Modified by Andrea Greco on 04/05/2015 to allow use of complex matrices !
    ! Modified by Jacek Dziedzic in August 2018 to set ngwf_cache_handle.     !
    ! Modified for embedding structures by Robert Charlton, October 2018.     !
    !=========================================================================!

    use comms, only: pub_on_root
    use model_type, only: MODEL
    use rundat, only: pub_any_nl_proj, pub_paw, pub_aug, pub_output_detail, &
         pub_hubbard, pub_num_kpoints, pub_num_spins, pub_eda_scfmi_any, &
         pub_block_orthogonalise
    use sparse_embed, only: sparse_embed_create
    use sparse_initialise, only: sparse_init_rep
    use utils, only: utils_alloc_check, utils_abort

    implicit none

    ! Arguments
    ! rc2013: EMBED_FIX!
    type(NGWF_REP), intent(inout), target :: rep
    character(*), intent(in) :: postfix
    type(MODEL), intent(in)  :: mdl    ! rc2013: needed for sparse_initialise
    ! agrecocmplx: optional argument for complex matrices and NGWFs
    logical, intent(in), optional :: is_cmplx

    ! Local variables
    ! agrecocmplx
    logical :: loc_cmplx

    ! Local variables
    integer :: ierr

    ! -------------------------------------------------------------------------

    ! Set postfix, if present
    rep%postfix = postfix

    ! agrecocmplx
    if (present(is_cmplx)) then
        loc_cmplx = is_cmplx
    else
        loc_cmplx = .false.
    end if

    ! Initialise matrices for this representation
    call sparse_init_rep(mdl, postfix)

    ! Banner before allocating all these matrices so it is clear where
    ! it goes wrong if it runs out of memory
    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) then
       if (postfix=='c') then
          if (pub_on_root) write(stdout,'(a)',advance='no') &
               'Conduction NGWF sparse matrix initialisation ...'
       else if (postfix=='j') then
          if (pub_on_root) write(stdout,'(a)',advance='no') &
               'Joint Valence + Conduction NGWF sparse matrix initialisation ...'
       else if (postfix=='a') then
          if (pub_on_root) write(stdout,'(a)',advance='no') &
               'Auxiliary NGWF sparse matrix initialisation ...'
       else if (postfix=='z') then
          if (pub_on_root) write(stdout,'(a)',advance='no') &
               'JointPH NGWF sparse matrix initialisation ...'
       else
          if (pub_on_root) write(stdout,'(a)',advance='no') &
               'Sparse matrix initialisation ...'
       end if
    end if

    ! rc2013: store the number of subsystems for this rep
    rep%nsub = mdl%nsub
    if( .not. allocated(rep%ngwfs_on_grid)) then
       allocate(rep%ngwfs_on_grid(rep%nsub), stat=ierr)
       call utils_alloc_check('ngwf_rep_create', 'rep%ngwfs_on_grid', ierr)
    end if

    ! Create the matrix structure for overlap
    rep%overlap%structure = 'S'//rep%postfix
    call sparse_embed_create(rep%overlap, iscmplx=loc_cmplx, &
         mrows=rep%nsub, ncols=rep%nsub)

    ! mjsp: Create the matrix structure for SCF-MI full overlap
    if (pub_eda_scfmi_any) then
       rep%overlap_scfmi_full%structure = 'S'//trim(rep%postfix)
       call sparse_embed_create(rep%overlap_scfmi_full,iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
    end if

    ! Create the matrix structure for direct NGWF-only overlap
    ! We do not actually need the matrix to be allocated though, as
    ! it will only be used for its structure code and library entry
    if (pub_aug) then
       rep%ngwf_overlap%structure = 'bS'//rep%postfix
    else
       rep%ngwf_overlap%structure = 'S'//rep%postfix
    end if
    call sparse_embed_create(rep%ngwf_overlap, iscmplx=loc_cmplx, &
         mrows=rep%nsub, ncols=rep%nsub)

    ! Create the matrix structure for inverse overlap
    rep%inv_overlap%structure = 'K'//rep%postfix
    call sparse_embed_create(rep%inv_overlap, iscmplx=loc_cmplx, &
         mrows=rep%nsub, ncols=rep%nsub)

    ! mjsp: Create the matrix structure for SCF-MI full inverse overlap
    if (pub_eda_scfmi_any) then
       rep%inv_overlap_scfmi_full%structure = 'K'//rep%postfix
       call sparse_embed_create(rep%inv_overlap_scfmi_full,iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
    end if

    ! Create the NGWF-projector overlap matrix if required
    if (pub_any_nl_proj.or.pub_paw) then
       rep%sp_overlap%structure = 'Q'//rep%postfix
       call sparse_embed_create(rep%sp_overlap,iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
    else
       ! rc2013: allocate the arrays so that we can pass the SPAM3 to interfaces
       ! even if it's never used
       allocate(rep%sp_overlap%m(1,1), stat=ierr)
       call utils_alloc_check('ngwf_rep_create', 'sp_overlap%m', ierr)
       rep%sp_overlap%p => rep%sp_overlap%m(1,1)
    end if

    ! Create the NGWF-projector overlap matrix if required
    if (pub_hubbard) then
       rep%hub_overlap%structure = 'V'//rep%postfix
       call sparse_embed_create(rep%hub_overlap,iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
       rep%hub_overlap_t%structure = 'W'//rep%postfix
       call sparse_embed_create(rep%hub_overlap_t,iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
       if (pub_aug) then ! ddor: PAW+U
          rep%hub_proj_paw_overlap%structure = 'Y'//rep%postfix
          call sparse_embed_create(rep%hub_proj_paw_overlap,iscmplx=loc_cmplx, &
               mrows=rep%nsub, ncols=rep%nsub)
          rep%hub_ngwf_paw_overlap%structure = 'Q'//rep%postfix
          call sparse_embed_create(rep%hub_ngwf_paw_overlap,iscmplx=loc_cmplx, &
               mrows=rep%nsub, ncols=rep%nsub)
       else
          allocate(rep%hub_proj_paw_overlap%m(1,1), stat=ierr)
          call utils_alloc_check('ngwf_rep_create', 'hub_proj_paw_overlap%m', ierr)
          allocate(rep%hub_ngwf_paw_overlap%m(1,1), stat=ierr)
          call utils_alloc_check('ngwf_rep_create', 'hub_ngwf_proj_overlap%m', ierr)
          rep%hub_proj_paw_overlap%p => rep%hub_proj_paw_overlap%m(1,1)
          rep%hub_ngwf_paw_overlap%p => rep%hub_ngwf_paw_overlap%m(1,1)
       endif
       ! Create empty version to store unmodified sparsity pattern for
       ! extended-metric cDFT
       rep%hub_overlap_t_bare%structure = 'W'//rep%postfix
       call sparse_embed_create(rep%hub_overlap_t_bare,iscmplx=loc_cmplx,&
            mrows=rep%nsub, ncols=rep%nsub)
    else
       ! rc2013: hubbard is not yet fully compatible with embedding;
       ! allocate the arrays so that we can pass the SPAM3 to relevant subroutines
       allocate(rep%hub_overlap%m(1,1), stat=ierr)
       call utils_alloc_check('ngwf_rep_create', 'hub_overlap%m', ierr)
       allocate(rep%hub_overlap_t%m(1,1), stat=ierr)
       call utils_alloc_check('ngwf_rep_create', 'hub_overlap_t%m', ierr)
       allocate(rep%hub_proj_paw_overlap%m(1,1), stat=ierr)
       call utils_alloc_check('ngwf_rep_create', 'hub_proj_paw_overlap%m', ierr)
       allocate(rep%hub_ngwf_paw_overlap%m(1,1), stat=ierr)
       call utils_alloc_check('ngwf_rep_create', 'hub_ngwf_proj_overlap%m', ierr)
       allocate(rep%hub_overlap_t_bare%m(1,1), stat=ierr)
       call utils_alloc_check('ngwf_rep_create', 'hub_overlap_t_bare%m', ierr)
       rep%hub_overlap%p => rep%hub_overlap%m(1,1)
       rep%hub_overlap_t%p => rep%hub_overlap_t%m(1,1)
       rep%hub_proj_paw_overlap%p => rep%hub_proj_paw_overlap%m(1,1)
       rep%hub_ngwf_paw_overlap%p => rep%hub_ngwf_paw_overlap%m(1,1)
       rep%hub_overlap_t_bare%p => rep%hub_overlap_t_bare%m(1,1)
    end if

    ! If we have norm-conserving PPs and nonlocal projectors, create
    ! nonlocal potential matrix
    if (pub_any_nl_proj.and.(.not.pub_aug)) then
       rep%nonlocpot%structure = 'H'//rep%postfix
       call sparse_embed_create(rep%nonlocpot, iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
    end if

    ! Create the kinetic part of the hamiltonian
    call sparse_embed_create(rep%kinet, rep%ngwf_overlap)

    ! ndmh: ensure inverse overlap is (re-)initialised for this rep
    rep%inv_overlap_init = .false.
    if (pub_eda_scfmi_any) rep%inv_overlap_b_init = .false.

    ! lr408: If this is a conduction representation, create the cross overlap
    ! lr408: matrix
    if ((rep%postfix=='c').or.(rep%postfix=='j').or.(rep%postfix=='a').or. &
         (rep%postfix=='z')) then
       if (pub_aug) then
          rep%ngwf_cross_overlap%structure = 'bT'//rep%postfix
          rep%ngwf_cross_overlap_tr%structure = 'bU'//rep%postfix
       else
          rep%ngwf_cross_overlap%structure = 'T'//rep%postfix
          rep%ngwf_cross_overlap_tr%structure = 'U'//rep%postfix
       end if
       rep%cross_overlap%structure = 'T'//rep%postfix
       rep%cross_overlap_tr%structure = 'U'//rep%postfix
       call sparse_embed_create(rep%cross_overlap,iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
       call sparse_embed_create(rep%cross_overlap_tr,iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
       call sparse_embed_create(rep%ngwf_cross_overlap,iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
       call sparse_embed_create(rep%ngwf_cross_overlap_tr,iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
    end if

    ! rc2013: create block-orthogonalised projectors if needed
    if(pub_block_orthogonalise) then
       call sparse_embed_create(rep%bo_projector, rep%inv_overlap, rep%overlap)
       call sparse_embed_create(rep%bo_projector_t, rep%overlap, rep%inv_overlap)
    end if


    ! Allocate n_occ according to public numbers of spins and k-points.
    ! agreco: currently pub_num_kpoints=1 as a parameter, remember to
    ! merge pub_num_kpoints_temp and pub_num_kpoints when everything
    ! is working
    ! allocate(rep%n_occ(pub_num_spins, pub_num_kpoints_temp), stat=ierr)
    allocate(rep%n_occ(pub_num_spins, pub_num_kpoints), stat=ierr)
    call utils_alloc_check('ngwf_rep_create', 'rep%n_occ', ierr)

    ! jd: Set NGWF cache handle
    if(postfix == '') then
       rep%ngwf_cache_handle = 1  ! valence
    else if(postfix == 'c') then
       rep%ngwf_cache_handle = 2  ! conduction
    else if(postfix == 'j') then
       rep%ngwf_cache_handle = 3  ! joint
    else if(postfix == 'a') then
       rep%ngwf_cache_handle = 4  ! auxiliary
    else if(postfix == 'z') then
       rep%ngwf_cache_handle = 5  ! jointPH
    ! jd: 6 is reserved for vacuum state NGWFs in polarisable embedding
    !       (the relevant %ngwf_cache_handle is set there).
    ! jd: If this list is ever extended, be sure
    !     - not to use 6
    !     - to increase remote::MAX_N_NGWF_CACHES
    else
       call utils_abort('ngwf_rep_create: Unrecognized postfix "'//&
            trim(postfix)//'".')
    end if

    ! jd: By default every NGWF_REP is invalidatable
    rep%invalidatable = .true.

  end subroutine ngwf_rep_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_rep_destroy(rep)

    !=========================================================================!
    ! This subroutine deallocates memory for the arrays in the NGWF_REP type, !
    ! apart from the density-dependent matrices of the Hamiltonian (which are !
    ! done separately to save memory).                                        !
    ! Also frees SW_EX containers holding NGWF_REP's SW expansions and the    !
    ! NGWF cache.                                                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rep (inout) : Container for NGWFs and matrices to be deallocated.     !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 07/10/2010                                  !
    ! Extended by Jacek Dziedzic in October 2016.                             !
    ! Modified for embedding structures by Robert Charlton, May 2017.         !
    !=========================================================================!

    use constants, only: REP_N_SWEXES
    use hash_table, only: hash_table_purge
    use remote, only: remote_ngwf_cache_hts, MAX_N_NGWF_CACHES
    use rundat, only: pub_any_nl_proj, pub_aug, pub_paw, pub_hubbard, &
         pub_eda_scfmi_any, pub_use_remote_ngwfs, pub_block_orthogonalise
    use sparse_embed, only: sparse_embed_destroy
    use sw_expansion_type, only: swx_destroy_container
    use utils, only: utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout) :: rep

    ! Local variables
    integer :: swex_h
    integer :: ierr

    deallocate(rep%ngwfs_on_grid, stat=ierr)
    call utils_dealloc_check('ngwf_rep_create', 'rep%ngwfs_on_grid', ierr)

    if ((rep%postfix=='c').or.(rep%postfix=='j').or.(rep%postfix=='a').or. &
         (rep%postfix=='z')) then
       call sparse_embed_destroy(rep%cross_overlap_tr)
       call sparse_embed_destroy(rep%cross_overlap)
    end if

    ! cks: memory deallocation for sparse matrices
    if (pub_hubbard) then
       if (pub_aug) then ! ddor: PAW+U
          call sparse_embed_destroy(rep%hub_ngwf_paw_overlap)
          call sparse_embed_destroy(rep%hub_proj_paw_overlap)
       endif
       call sparse_embed_destroy(rep%hub_overlap_t)
       call sparse_embed_destroy(rep%hub_overlap)
    ! jcap: need to deallocate these matrices even if pub_hubbard is
    ! false, to undo what is done in ngwf_rep_create
    else
       deallocate(rep%hub_overlap%m,stat=ierr)
       call utils_dealloc_check('ngwf_rep_destroy','hub_overlap%m',ierr)
       deallocate(rep%hub_overlap_t%m,stat=ierr)
       call utils_dealloc_check('ngwf_rep_destroy','hub_overlap_t%m',ierr)
       deallocate(rep%hub_proj_paw_overlap%m,stat=ierr)
       call utils_dealloc_check('ngwf_rep_destroy','hub_proj_paw_overlap%m',ierr)
       deallocate(rep%hub_ngwf_paw_overlap%m,stat=ierr)
       call utils_dealloc_check('ngwf_rep_destroy','hub_ngwf_proj_overlap%m',ierr)
       deallocate(rep%hub_overlap_t_bare%m,stat=ierr)
       call utils_dealloc_check('ngwf_rep_destroy','hub_overlap_t_bare%m',ierr)
    end if
    if (pub_any_nl_proj.or.pub_paw) then
       call sparse_embed_destroy(rep%sp_overlap)
    end if
    if (pub_any_nl_proj.and.(.not.pub_aug)) then
       call sparse_embed_destroy(rep%nonlocpot)
    end if

    rep%inv_overlap_init = .false.
    call sparse_embed_destroy(rep%inv_overlap)

    call sparse_embed_destroy(rep%kinet)
    call sparse_embed_destroy(rep%overlap)
    if (pub_eda_scfmi_any) then
       rep%inv_overlap_b_init = .false.
       call sparse_embed_destroy(rep%inv_overlap_scfmi_full)
       call sparse_embed_destroy(rep%overlap_scfmi_full)
    end if

    if(pub_block_orthogonalise) then
       call sparse_embed_destroy(rep%bo_projector)
       call sparse_embed_destroy(rep%bo_projector_t)
    end if

    ! Deallocate n_occ.
    deallocate(rep%n_occ, stat=ierr)
    call utils_dealloc_check('ngwf_rep_destroy', 'rep%n_occ', ierr)

    ! jd: Destroy all SW_EX containers in use
    do swex_h = 1, REP_N_SWEXES
       if(rep%swexes(swex_h)%initialised) then
          call swx_destroy_container(rep%swexes(swex_h), rep%postfix)
       end if
    end do

    ! jd: Purge the NGWF cache corresponding to this rep. The structure stays,
    !     only the contents are purged.
    if(pub_use_remote_ngwfs) then
       call utils_assert(rep%ngwf_cache_handle <= MAX_N_NGWF_CACHES .and. &
            rep%ngwf_cache_handle > 0, 'ngwf_rep_destroy: &
            &ngwf_cache_handle out of bounds', rep%ngwf_cache_handle)
       call hash_table_purge(remote_ngwf_cache_hts(rep%ngwf_cache_handle))
    end if

  end subroutine ngwf_rep_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_rep_register_change(rep, where, updated_by_nnho)

    !=========================================================================!
    ! This subroutine should be called any time NGWFs change. Actions that    !
    ! should be taken each time NGWFs change can be collected here.           !
    !                                                                         !
    ! Currently performed actions:                                            !
    !  - Rep's SPAM3 matrices are set to garbage.                             !
    !    This will help detect erroneous uses of invalidated NGWFs.           !
    !  - Current SW expansions of NGWFs are invalidated, if in use.           !
    !  - The NGWF cache is invalidated.                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rep (inout): Container to operate on.                                 !
    !   where (in): Human-readable description of where the change happened.  !
    !   updated_by_nnho(in, opt): Set to .true. if the NGWF change is due to  !
    !                             an nnho transformation. In these case only  !
    !                             a subset of matrices is overwritten with    !
    !                             garbage, because nnho takes care to update  !
    !                             the others.                                 !
    !-------------------------------------------------------------------------!
    ! Points of note:                                                         !
    !   - inv_overlap is not smashed, as that results in a performance hit    !
    !     from re-calculating it afresh (rather than starting from previous)  !
    !     every time NGWFs change.                                            !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2016.10.14.                                !
    ! Adapted for SPAM3_EMBED by Robert Charlton, 10th August 2018.           !
    !=========================================================================!

    use constants, only: stdout, REP_N_SWEXES
    use hash_table, only: hash_table_purge
    use remote, only: remote_ngwf_cache_hts
    use rundat, only: pub_eda_scfmi_any, pub_any_nl_proj, pub_paw, pub_aug, &
         pub_hubbard, pub_debug_on_root, pub_use_remote_ngwfs
    use sparse_embed, only: sparse_embed_set_to_garbage
    use sw_expansion_type, only: swx_invalidate_expansion
    use utils, only: utils_postfix_to_ngwf_set_name

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout) :: rep
    character(len=*), intent(in)  :: where
    logical, intent(in), optional :: updated_by_nnho

    ! Local variables
    integer :: swex_h
    logical :: local_updated_by_nnho

    ! -------------------------------------------------------------------------

    if(.not. rep%invalidatable) return

    if(pub_debug_on_root) then
       write(stdout,'(a,a,a)') 'DEBUG: NGWFs in NGWF_REP '''//&
       trim(utils_postfix_to_ngwf_set_name(rep%postfix))//''' changed: ', &
       where, '.'
    end if

    local_updated_by_nnho = .false.
    if(present(updated_by_nnho)) then
       local_updated_by_nnho = updated_by_nnho
    end if

    ! jd: These are matrices that an nnho transformation updates, so smash them
    !     only if this NGWF change did not result from an nnho transformation
    if(.not. local_updated_by_nnho) then
       call sparse_embed_set_to_garbage(rep%overlap)
       call sparse_embed_set_to_garbage(rep%kinet)
       if (pub_any_nl_proj .or. pub_paw) then
          call sparse_embed_set_to_garbage(rep%sp_overlap)
       end if
       if (pub_hubbard) then
          call sparse_embed_set_to_garbage(rep%hub_overlap)
       end if
       if (pub_any_nl_proj.and.(.not.pub_aug)) then
          call sparse_embed_set_to_garbage(rep%nonlocpot)
       end if
    end if

    ! jd: These matrices are not updated by nnho transformation, so need to be
    !     smashed in any case
    if (pub_eda_scfmi_any) then
       call sparse_embed_set_to_garbage(rep%overlap_scfmi_full)
       call sparse_embed_set_to_garbage(rep%inv_overlap_scfmi_full)
       rep%inv_overlap_b_init = .false.
    end if
    if (pub_hubbard) then
       call sparse_embed_set_to_garbage(rep%hub_overlap_t)
       if (pub_aug) then
          call sparse_embed_set_to_garbage(rep%hub_proj_paw_overlap)
          call sparse_embed_set_to_garbage(rep%hub_ngwf_paw_overlap)
       endif
    end if
    if ((rep%postfix=='c').or.(rep%postfix=='j').or.(rep%postfix=='a')) then
       call sparse_embed_set_to_garbage(rep%cross_overlap)
       call sparse_embed_set_to_garbage(rep%cross_overlap_tr)
    end if

    ! jd: Invalidate the expansions in SW_EX containers
    do swex_h = 1, REP_N_SWEXES
       if(rep%swexes(swex_h)%initialised) then
          call swx_invalidate_expansion(rep%swexes(swex_h), rep%postfix)
       end if
    end do

    ! jd: Invalidate the NGWF cache
    if(pub_use_remote_ngwfs) then
       call hash_table_purge(remote_ngwf_cache_hts(rep%ngwf_cache_handle))
    end if

  end subroutine ngwf_rep_register_change

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_ham_create(ham,rep)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the NGWF_HAM type    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ham (inout) : Container for Hamiltonian matrices to be allocated.     !
    !   rep (in): The NGWF_REP associated with this NGWF_HAM.                 !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 07/10/2010                                  !
    ! Modified by Andrea Greco on 04/05/2015 to allow complex matrices.       !
    ! Got rid of 'lhxc_nl_only' (not used anywhere) -- Jacek Dziedzic 2016.10.!
    ! Modified for embedding structures by Robert Charlton, August 2017.      !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: VERBOSE
    use rundat, only: pub_aug, pub_hubbard, cond_init_shift, pub_num_spins, &
         pub_use_hfx, pub_confined_ngwfs, pub_pol_emb_qmstar, &
         pub_pol_emb_vacuum_qmstar, pub_xc_ke_density_required, &
         pub_output_detail, pub_project_embed, pub_use_activehfx, &
         pub_scissor_ngroups
    use sparse_embed, only: sparse_embed_create
    use utils, only: utils_alloc_check, utils_postfix_to_ngwf_set_name
    ! agreco: dependence on sparse_array module
    !use sparse_array, only: sparse_array_create

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout)     :: ham
    type(NGWF_REP), intent(in)        :: rep

    ! Local Variables
    integer :: is
    integer :: ierr
    logical :: loc_cond
    ! agrecocmplx: local variable to use complex matrices
    logical :: loc_cmplx

    ! -------------------------------------------------------------------------

    if (rep%postfix=='c') then
       loc_cond = .true.
    else
       loc_cond = .false.
    end if

    ! agrecocmplx: loc_cmplx true if using rep with complex NGWFs
    loc_cmplx = rep%ngwf_overlap%iscmplx

    ! Allocate matrices for local potential and final Hamiltonian
    ! agreco: include k-point dependence, use SPAM3_ARRAY type
    ! lhxc should be k-point independent in kp method, but k-point
    ! dependent in tight-binding method
    allocate(ham%ham(pub_num_spins),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','ham%ham',ierr)
    allocate(ham%lhxc(pub_num_spins),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','ham%lhxc',ierr)

    ! ny: For species dependent scissor shifted hamiltonian
    if (pub_scissor_ngroups > 0) then
       allocate(ham%scissor(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%scissor',ierr)
    end if

    ! JCW: For tau-dependent XC functional, allocate matrix for dfdtau
    if (pub_xc_ke_density_required) then
      allocate(ham%dfdtau(pub_num_spins),stat=ierr)
      call utils_alloc_check('ngwf_ham_create','ham%dfdtau',ierr)
    end if

    ! gcc32: create confinement hamiltonian if confined ngwfs is used
    ! gcc32: make it block diagonal, as only Hconf_{\alpha\alpha} != 0
    if(pub_confined_ngwfs) then
       allocate(ham%confinement(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%confinement',ierr)
       do is=1,pub_num_spins
          ham%confinement(is)%structure = 'D'//rep%postfix
          call sparse_embed_create(ham%confinement(is),iscmplx=loc_cmplx)
       end do
    end if

    ! ndmh: in PAW, nonlocpot is density-dependent
    if (pub_aug) then
       allocate(ham%nonlocpot(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%nonlocpot',ierr)
       allocate(ham%dijhat(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%dijhat',ierr)
    else
       allocate(ham%nonlocpot(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%nonlocpot',ierr)
       allocate(ham%dijhat(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%dijhat',ierr)
    end if

    ! qoh: Allocate Hartree-Fock exchange sparse matrix if necessary
    if (pub_use_hfx .or. pub_use_activehfx) then
       allocate(ham%hfexchange(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%hfexchange',ierr)
    else
       allocate(ham%hfexchange(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%hfexchange',ierr)
    end if

    ! jd: Allocate polarisable embedding sparse matrix.
    !     NB: In the absence of polarisable embedding this still needs to be
    !         allocated (but is not sparse_created), because we're passing
    !         this around as an inout arg ('attempt to use pointer...').
    allocate(ham%polemb(2),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','ham%polemb',ierr)
    allocate(ham%vacuum_polemb(2),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','ham%vacuum_polemb',ierr)
    allocate(ham%vacuum_aux_polemb(2),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','ham%vacuum_aux_polemb',ierr)

    ! ddor: Allocate DFT+U Hamiltonian if necessary
    if (pub_hubbard) then
       allocate(ham%hubbard_ham(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%hubbard_ham',ierr)
    else
       allocate(ham%hubbard_ham(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%hubbard_ham',ierr)
    endif

    ! lr408: Allocate conduction projected Hamiltonian if necessary -
    ! lr408: only needed if this is for the cond rep
    if (loc_cond) then
       allocate(ham%cond_non_proj_ham(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%cond_non_proj_ham',ierr)
       ham%cond_shift = cond_init_shift
    else
       allocate(ham%cond_non_proj_ham(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%cond_non_proj_ham',ierr)
       ham%cond_shift = 0.0_DP
    end if

    ! rc2013: Allocate projected embedding Hamiltonian if needed
    !if (pub_project_embed) then
    !   allocate(ham%proj_embed_ham(pub_num_spins), stat=ierr)
    !   call utils_alloc_check('ngwf_ham_create', 'ham%proj_embed_ham', ierr)
    !   allocate(ham%unproj_ham(pub_num_spins), stat=ierr)
    !   call utils_alloc_check('ngwf_ham_create', 'ham%proj_ham', ierr)
    !endif

    do is=1,pub_num_spins
       if (loc_cond) then
          ham%ham(is)%structure = 'L'//rep%postfix
       else
          ham%ham(is)%structure = 'H'//rep%postfix
       end if
       call sparse_embed_create(ham%ham(is), iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
       call sparse_embed_create(ham%lhxc(is), rep%ngwf_overlap, &
            iscmplx=loc_cmplx)
       if (pub_scissor_ngroups > 0) then
          call sparse_embed_create(ham%scissor(is), rep%overlap, &
               iscmplx=loc_cmplx)
       end if
       if (pub_xc_ke_density_required) then
          ! JCW: Create dfdtau matrix with sparsity pattern of overlap matrix
          call sparse_embed_create(ham%dfdtau(is),rep%ngwf_overlap, &
               iscmplx=loc_cmplx)
       end if
       if (pub_aug) then
          ham%nonlocpot%structure = 'H'//rep%postfix
          call sparse_embed_create(ham%nonlocpot(is), iscmplx=loc_cmplx, &
               mrows=rep%nsub, ncols=rep%nsub)
          ham%dijhat%structure = 'E'
          call sparse_embed_create(ham%dijhat(is),mrows=rep%nsub,ncols=rep%nsub)
       end if
       ! rc2013: create the full system exchange matrix even if
       !         we only populate part of it
       if (pub_use_hfx .or. pub_use_activehfx) then
          if(pub_output_detail >= VERBOSE .and. pub_on_root) then
             write(stdout,'(a,a,a,i0,a)') 'HFx: Creating X matrix structure for ', &
                  trim(utils_postfix_to_ngwf_set_name(rep%postfix)), ' rep (spin ', is, ').'
          end if
          ham%hfexchange(is)%structure = 'HFx'//rep%postfix
          call sparse_embed_create(ham%hfexchange(is), iscmplx=loc_cmplx, &
               mrows=rep%nsub, ncols=rep%nsub)
       end if
       if (pub_hubbard) then
          ham%hubbard_ham(is)%structure = 'H'//rep%postfix
          call sparse_embed_create(ham%hubbard_ham(is), iscmplx=loc_cmplx, &
              mrows=rep%nsub, ncols=rep%nsub)
       end if
       if (loc_cond) then
          ham%cond_non_proj_ham(is)%structure = 'H'//rep%postfix
          call sparse_embed_create(ham%cond_non_proj_ham(is),iscmplx=loc_cmplx,&
              mrows=rep%nsub, ncols=rep%nsub)
       end if
       !if (pub_project_embed) then
       !   ! rc2013: new sparsity pattern?
       !   ham%proj_embed_ham(is)%structure = 'H'//rep%postfix
       !   call sparse_embed_create(ham%proj_embed_ham(is), iscmplx=loc_cmplx, &
       !       mrows=rep%nsub, ncols=rep%nsub)
       !   ham%unproj_ham(is)%structure = 'H'//rep%postfix
       !   call sparse_embed_create(ham%unproj_ham(is), iscmplx=loc_cmplx, &
       !       mrows=rep%nsub, ncols=rep%nsub)
       !end if
    end do

    if(pub_pol_emb_qmstar) then
       ham%polemb(1)%structure = 'S'
       ham%polemb(2)%structure = 'S'
       call sparse_embed_create(ham%polemb(1),iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
       call sparse_embed_create(ham%polemb(2),iscmplx=loc_cmplx, &
            mrows=rep%nsub, ncols=rep%nsub)
       if(pub_pol_emb_vacuum_qmstar) then
          call sparse_embed_create(ham%vacuum_polemb(1), ham%polemb(1))
          call sparse_embed_create(ham%vacuum_polemb(2), ham%polemb(2))
          call sparse_embed_create(ham%vacuum_aux_polemb(1), ham%polemb(1))
          call sparse_embed_create(ham%vacuum_aux_polemb(2), ham%polemb(2))
       end if
    end if

    ! agreco: equivalent versions using sparse_array routines
    ! jd: now out of date wrt code above 2016.10.18
    !if (loc_cond) then
    !   call sparse_array_create(ham%ham, n_spins=pub_num_spins, &
    !      n_kpoints=pub_num_kpoints, structure='L'//rep%postfix, iscmplx=loc_cmplx)
    !else
    !   call sparse_array_create(ham%ham, n_spins=pub_num_spins, &
    !      n_kpoints=pub_num_kpoints, structure='H'//rep%postfix, iscmplx=loc_cmplx)
    !end if
    !
    !call sparse_array_create(ham%lhxc, n_spins=pub_num_spins, &
    !     n_kpoints=pub_num_kpoints, spammata=rep%ngwf_overlap, iscmplx=loc_cmplx)
    !
    !if (pub_aug) then
    !   call sparse_array_create(ham%nonlocpot, n_spins=pub_num_spins, &
    !      n_kpoints=pub_num_kpoints, structure='H'//rep%postfix, iscmplx=loc_cmplx)
    !   call sparse_array_create(ham%dijhat, n_spins=pub_num_spins, &
    !      n_kpoints=pub_num_kpoints, structure='E')
    ! agreco: dummy allocation
    !else
    !   call sparse_array_create(ham%nonlocpot, n_spins=0, &
    !      n_kpoints=0, structure='H'//rep%postfix, iscmplx=loc_cmplx)
    !   call sparse_array_create(ham%dijhat, n_spins=0, &
    !      n_kpoints=0, structure='E')
    !end if
    !
    !if (loc_use_hfx) then
    !   call sparse_array_create(ham%hfexchange, n_spins=pub_num_spins, &
    !      n_kpoints=pub_num_kpoints, structure='HFx', iscmplx=loc_cmplx)
    !else
    !   call sparse_array_create(ham%hfexchange, n_spins=0, &
    !      n_kpoints=0, structure='HFx', iscmplx=loc_cmplx)
    !end if
    !
    !if (loc_hubbard) then
    !   call sparse_array_create(ham%hubbard_ham, n_spins=pub_num_spins, &
    !      n_kpoints=pub_num_kpoints, structure='H'//rep%postfix, iscmplx=loc_cmplx)
    !else
    !   call sparse_array_create(ham%hubbard_ham, n_spins=0, &
    !      n_kpoints=0, structure='H'//rep%postfix, iscmplx=loc_cmplx)
    !end if
    !
    !if (loc_cond) then
    !   call sparse_array_create(ham%cond_non_proj_ham, n_spins=pub_num_spins, &
    !      n_kpoints=pub_num_kpoints, structure='H'//rep%postfix, iscmplx=loc_cmplx)
    !   ham%cond_shift = cond_init_shift
    !else
    !   call sparse_array_create(ham%cond_non_proj_ham, n_spins=0, &
    !      n_kpoints=0, structure='H'//rep%postfix, iscmplx=loc_cmplx)
    !end if
    !
    !if(pub_confined_ngwfs) then
    !   call sparse_array_create(ham%confinement, n_spins=pub_num_spins, &
    !      n_kpoints=pub_num_kpoints, structure='D'//rep%postfix, iscmplx=loc_cmplx)
    !end if
    !
    ! agreco: is this k-point dependent?
    !if (pub_pol_emb_pot) then
    !   call sparse_array_create(ham%polemb, n_spins=1, &
    !      n_kpoints=pub_num_kpoints, structure='S', iscmplx=loc_cmplx)
    !end if

  end subroutine ngwf_ham_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_ham_destroy(ham)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the NGWF_REP type,   !
    ! apart from the density-dependent matrices of the Hamiltonian (which are !
    ! done separately to save memory).                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ham (inout) : Container for Hamiltonian matrices to be deallocated.   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 07/10/2010                                  !
    ! Modified for embedding by Robert Charlton, July 2017.                   !
    ! Modified to destroy everything created in ngwf_ham_create by      !
    ! Joseph Prentice, June 2018                                              !
    !=========================================================================!

    use rundat, only: pub_aug, pub_hubbard, pub_num_spins, pub_use_hfx, &
         pub_pol_emb_qmstar, pub_pol_emb_vacuum_qmstar, &
         pub_xc_ke_density_required, pub_confined_ngwfs, pub_project_embed, &
         pub_scissor_ngroups
    use sparse_embed, only: sparse_embed_destroy
    use utils, only: utils_dealloc_check
    ! agreco: add dependence on sparse_array
    !use sparse_array, only: sparse_array_destroy

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout) :: ham

    ! Local Variables
    integer :: is
    integer :: ierr

    ! Deallocate storage in matrices
    if(pub_pol_emb_qmstar) then
       call sparse_embed_destroy(ham%polemb(1))
       call sparse_embed_destroy(ham%polemb(2))
       if(pub_pol_emb_vacuum_qmstar) then
          call sparse_embed_destroy(ham%vacuum_polemb(1))
          call sparse_embed_destroy(ham%vacuum_polemb(2))
          call sparse_embed_destroy(ham%vacuum_aux_polemb(1))
          call sparse_embed_destroy(ham%vacuum_aux_polemb(2))
       end if
    end if
    do is=pub_num_spins,1,-1
       if (size(ham%cond_non_proj_ham)>0) &
            call sparse_embed_destroy(ham%cond_non_proj_ham(is))
       if (pub_hubbard.and.size(ham%hubbard_ham)>0) &
            call sparse_embed_destroy(ham%hubbard_ham(is))
       if (pub_use_hfx.and.size(ham%hfexchange)>0) &
            call sparse_embed_destroy(ham%hfexchange(is))
       if (pub_aug) then
          call sparse_embed_destroy(ham%dijhat(is))
          call sparse_embed_destroy(ham%nonlocpot(is))
       end if
       if (pub_xc_ke_density_required) then
          ! JCW: Destroy dfdtau matrix, if KE-density-dependent functional is
          ! JCW: used
          call sparse_embed_destroy(ham%dfdtau(is))
       end if
       !if (pub_project_embed) then
       !   call sparse_embed_destroy(ham%proj_embed_ham(is))
       !   call sparse_embed_destroy(ham%unproj_ham(is))
       !end if
       if (pub_scissor_ngroups > 0) call sparse_embed_destroy(ham%scissor(is))
       call sparse_embed_destroy(ham%lhxc(is))
       call sparse_embed_destroy(ham%ham(is))
    end do

    ! Deallocate all matrix arrays
    if (allocated(ham%cond_non_proj_ham)) then
       deallocate(ham%cond_non_proj_ham,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%cond_non_proj_ham',ierr)
    end if
    if (allocated(ham%hubbard_ham)) then
       deallocate(ham%hubbard_ham,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%hubbard_ham',ierr)
    end if
    if (allocated(ham%polemb)) then
       deallocate(ham%polemb,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%polemb',ierr)
    end if
    if (allocated(ham%vacuum_polemb)) then
       deallocate(ham%vacuum_polemb,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%vacuum_polemb',ierr)
    end if
    if (allocated(ham%vacuum_aux_polemb)) then
       deallocate(ham%vacuum_aux_polemb,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%vacuum_aux_polemb',ierr)
    end if
    if (allocated(ham%hfexchange)) then
       deallocate(ham%hfexchange,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%hfexchange',ierr)
    end if
    if (allocated(ham%nonlocpot)) then
       deallocate(ham%nonlocpot,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%nonlocpot',ierr)
    end if
    if (allocated(ham%dijhat)) then
       deallocate(ham%dijhat,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%dijhat',ierr)
    end if
    ! JCW: For tau-dependent XC functional, allocate matrix for dfdtau
    if (allocated(ham%dfdtau)) then
       deallocate(ham%dfdtau,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%dfdtau',ierr)
    end if
    if (allocated(ham%lhxc)) then
       deallocate(ham%lhxc,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%lhxc',ierr)
    end if
    if (allocated(ham%scissor)) then
       deallocate(ham%scissor,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%scissor',ierr)
    end if
    if (allocated(ham%ham)) then
       deallocate(ham%ham,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%ham',ierr)
    end if

    if(allocated(ham%confinement)) then
       deallocate(ham%confinement,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%confinement',ierr)
    end if
    if (allocated(ham%proj_embed_ham)) then
       deallocate(ham%proj_embed_ham,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%proj_embed_ham',ierr)
    end if
    if (allocated(ham%unproj_ham)) then
       deallocate(ham%unproj_ham,stat=ierr)
       call utils_dealloc_check('ngwf_ham_destroy','ham%unproj_ham',ierr)
    end if

    ! agreco: equivalent version using sparse_array routines
    !if (pub_pol_emb_pot .and. ham%polemb%num_spins > 0 .and. &
    !   ham%polemb%num_kpoints > 0) then
    !   call sparse_array_destroy(ham%polemb)
    !end if
    !
    !if
    !
    !if (ham%cond_non_proj_ham%num_spins > 0 .and. &
    !   ham%cond_non_proj_ham%num_kpoints > 0) then
    !   call sparse_array_destroy(ham%cond_non_proj_ham)
    !end if
    !
    !if (pub_hubbard .and. ham%hubbard_ham%num_spins > 0 .and. &
    !   ham%hubbard_ham%num_kpoints > 0) then
    !   call sparse_array_destroy(ham%hubbard_ham)
    !end if
    !
    !if (pub_use_hfx .and. ham%hfexchange%num_spins > 0 .and. &
    !   ham%hfexchange%num_kpoints > 0) then
    !   call sparse_array_destroy(ham%hfexchange)
    !end if
    !
    !if (pub_confined_ngwfs .and. ham%confinement%num_spins > 0 .and. &
    !   ham%confinement%num_kpoints > 0) then
    !   call sparse_array_destroy(ham%confinement)
    !end if
    !
    !if (pub_aug) then
    !   call sparse_array_destroy(ham%dijhat)
    !   call sparse_array_destroy(ham%nonlocpot)
    !end if
    !
    !call sparse_array_destroy(ham%ham)
    !call sparse_array_destroy(ham%lhxc)


  end subroutine ngwf_ham_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_ham_copy(ham_dst,ham_src)

    !=========================================================================!
    ! This subroutine copies an NGWF_HAM type object (matrices only).         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  ham_src (inout) : Container for Hamiltonian matrices to be copied into.!
    !  ham_dst (inout) : Container for Hamiltonian matrices to be copied from.!
    !-------------------------------------------------------------------------!
    ! Written by Max Phipps April 2015.                                       !
    ! Modified from ngwf_ham_create.
    !=========================================================================!

    use rundat, only: pub_aug, pub_hubbard, pub_num_spins, pub_use_hfx
    use sparse_embed, only: sparse_embed_copy

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout)     :: ham_dst
    type(NGWF_HAM), intent(in)        :: ham_src

    ! Local Variables
    integer :: is

    do is=1,pub_num_spins
       ham_dst%ham(is)%structure = ham_src%ham(is)%structure
       call sparse_embed_copy(ham_dst%ham(is),ham_src%ham(is))
       call sparse_embed_copy(ham_dst%lhxc(is),ham_src%lhxc(is))
       if (pub_aug) then
          ham_dst%nonlocpot%structure = ham_src%nonlocpot%structure
          call sparse_embed_copy(ham_dst%nonlocpot(is),ham_src%nonlocpot(is))
          ham_dst%dijhat%structure = ham_src%dijhat%structure
          call sparse_embed_copy(ham_dst%dijhat(is),ham_src%dijhat(is))
       end if
       if (pub_use_hfx) then
          ham_dst%hfexchange(is)%structure = 'HFx'
          call sparse_embed_copy(ham_dst%hfexchange(is),ham_src%hfexchange(is))
       end if
       if (pub_hubbard) then
          ham_dst%hubbard_ham(is)%structure = ham_src%hubbard_ham(is)%structure
          call sparse_embed_copy(ham_dst%hubbard_ham(is),ham_src%hubbard_ham(is))
       end if
       if (size(ham_src%cond_non_proj_ham)>0) then
          ham_dst%cond_non_proj_ham(is)%structure = &
            ham_src%cond_non_proj_ham(is)%structure
          call sparse_embed_copy(ham_dst%cond_non_proj_ham(is), &
            ham_src%cond_non_proj_ham(is))
       end if
    end do

  end subroutine ngwf_ham_copy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_ham_register_change(ham, rep, updated_by_nnho)

    !=========================================================================!
    ! This subroutine should be called any time NGWFs change. Actions that    !
    ! should be taken each time NGWFs change can be collected here.           !
    !                                                                         !
    ! Currently performed actions:                                            !
    !  - Ham's SPAM3 matrices are set to garbage.                             !
    !    This will help detect erroneous uses of invalidated NGWFs.           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ham (inout): Container to operate on.
    !   rep (in): NGWF_REP corresponding to ham. Needed only for its postfix, !
    !             so that it we know if ham is a conduction Hamiltonian.      !
    !   updated_by_nnho(in, opt): Set to .true. if the NGWF change is due to  !
    !                             an nnho transformation. In these case only  !
    !                             a subset of matrices is overwritten with    !
    !                             garbage, because nnho takes care to update  !
    !                             the others.                                 !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic on 2016.10.18.                                !
    ! Adapted for SPAM3_EMBED by Robert Charlton, 28th August 2018.           !
    !=========================================================================!

    use constants, only: stdout
    use rundat, only: pub_debug_on_root, pub_xc_ke_density_required, &
         pub_aug, pub_use_hfx, pub_pol_emb_qmstar, pub_pol_emb_vacuum_qmstar, &
         pub_confined_ngwfs, pub_hubbard, pub_use_activehfx
    use sparse_embed, only: sparse_embed_set_to_garbage
    use utils, only: utils_postfix_to_ngwf_set_name

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout) :: ham
    type(NGWF_REP), intent(in)    :: rep
    logical, intent(in), optional :: updated_by_nnho

    ! Local variables
    logical :: local_updated_by_nnho
    logical :: loc_cond
    integer :: i

    ! -------------------------------------------------------------------------

    if(pub_debug_on_root) then
       write(stdout,'(a)') 'DEBUG: Discarding NGWF_HAM of NGWF_REP '''//&
            trim(utils_postfix_to_ngwf_set_name(rep%postfix))//'''.'
    end if

    local_updated_by_nnho = .false.
    if(present(updated_by_nnho)) then
       local_updated_by_nnho = updated_by_nnho
    end if

    if (rep%postfix=='c') then
       loc_cond = .true.
    else
       loc_cond = .false.
    end if

    ! jd: These are matrices that an nnho transformation updates, so smash them
    !     only if this NGWF change did not result from an nnho transformation
    if(.not. local_updated_by_nnho) then
       if(allocated(ham%nonlocpot) .and. pub_aug) then
          do i = lbound(ham%nonlocpot,1), ubound(ham%nonlocpot,1)
             call sparse_embed_set_to_garbage(ham%nonlocpot(i))
          end do
       end if
       if(allocated(ham%lhxc)) then
          do i = lbound(ham%lhxc,1), ubound(ham%lhxc,1)
             call sparse_embed_set_to_garbage(ham%lhxc(i))
          end do
       end if
    end if

    ! jd: These matrices are not updated by nnho transformation, so need to be
    !     smashed in any case
    if(allocated(ham%dijhat) .and. pub_aug) then
       do i = lbound(ham%dijhat,1), ubound(ham%dijhat,1)
          ! jd: @@@ Breaks QC test #48, where dijhat is not updated after
          !     cond NGWFs are changed. Uncomment when this is fixed
          ! call sparse_embed_set_to_garbage(ham%dijhat(i))
       end do
    end if
    if(allocated(ham%hubbard_ham) .and. pub_hubbard) then
       do i = lbound(ham%hubbard_ham,1), ubound(ham%hubbard_ham,1)
          call sparse_embed_set_to_garbage(ham%hubbard_ham(i))
       end do
    end if
    if(allocated(ham%hfexchange) .and. (pub_use_hfx .or. pub_use_activehfx)) then
       do i = lbound(ham%hfexchange,1), ubound(ham%hfexchange,1)
          call sparse_embed_set_to_garbage(ham%hfexchange(i))
       end do
    end if
    if(allocated(ham%polemb) .and. pub_pol_emb_qmstar) then
       do i = lbound(ham%polemb,1), ubound(ham%polemb,1)
          call sparse_embed_set_to_garbage(ham%polemb(i))
       end do
    end if
    if(allocated(ham%vacuum_polemb) .and. pub_pol_emb_vacuum_qmstar) then
       do i = lbound(ham%vacuum_polemb,1), ubound(ham%vacuum_polemb,1)
          call sparse_embed_set_to_garbage(ham%vacuum_polemb(i))
       end do
    end if
    if(allocated(ham%vacuum_aux_polemb) .and. pub_pol_emb_vacuum_qmstar) then
       do i = lbound(ham%vacuum_aux_polemb,1), ubound(ham%vacuum_aux_polemb,1)
          call sparse_embed_set_to_garbage(ham%vacuum_aux_polemb(i))
       end do
    end if
    if(allocated(ham%dfdtau) .and. pub_xc_ke_density_required) then
       do i = lbound(ham%dfdtau,1), ubound(ham%dfdtau,1)
          call sparse_embed_set_to_garbage(ham%dfdtau(i))
       end do
    end if
    if(allocated(ham%confinement) .and. pub_confined_ngwfs) then
       do i = lbound(ham%confinement,1), ubound(ham%confinement,1)
          call sparse_embed_set_to_garbage(ham%confinement(i))
       end do
    end if
    if(allocated(ham%scissor)) then
       do i = lbound(ham%scissor,1), ubound(ham%scissor,1)
          call sparse_embed_set_to_garbage(ham%scissor(i))
       end do
    end if
    if(allocated(ham%ham)) then
       do i = lbound(ham%ham,1), ubound(ham%ham,1)
          call sparse_embed_set_to_garbage(ham%ham(i))
       end do
    end if
    if(allocated(ham%cond_non_proj_ham) .and. loc_cond) then
       do i = lbound(ham%cond_non_proj_ham,1), ubound(ham%cond_non_proj_ham,1)
          call sparse_embed_set_to_garbage(ham%cond_non_proj_ham(i))
       end do
    end if

  end subroutine ngwf_ham_register_change


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_ham_allocate(ham, use_hfx, use_hubbard, use_cond)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the NGWF_HAM type    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ham (inout) : Container for Hamiltonian matrices to be allocated.     !
    !-------------------------------------------------------------------------!
    ! Written by Robert Charlton, 26/04/2018.                                 !
    !=========================================================================!

    use model_type, only: MODEL
    use rundat, only: pub_aug, cond_init_shift, pub_num_spins, &
         pub_confined_ngwfs, &
         pub_xc_ke_density_required, pub_project_embed
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout) :: ham
    logical, intent(in)                 :: use_hfx
    logical, intent(in)                 :: use_hubbard
    logical, intent(in)                 :: use_cond

    ! Local Variables
    integer :: ierr
    ! agrecocmplx: local variable to use complex matrices


    ! Allocate matrices for local potential and final Hamiltonian
    ! agreco: include k-point dependence, use SPAM3_ARRAY type
    ! lhxc should be k-point independent in kp method, but k-point
    ! dependent in tight-binding method
    allocate(ham%ham(pub_num_spins),stat=ierr)
    call utils_alloc_check('ngwf_ham_allocate','ham%ham',ierr)
    allocate(ham%lhxc(pub_num_spins),stat=ierr)
    call utils_alloc_check('ngwf_ham_allocate','ham%lhxc',ierr)

    ! JCW: For tau-dependent XC functional, allocate matrix for dfdtau
    if (pub_xc_ke_density_required) then
      allocate(ham%dfdtau(pub_num_spins),stat=ierr)
      call utils_alloc_check('ngwf_ham_allocate','ham%dfdtau',ierr)
    end if

    ! gcc32: create confinement hamiltonian if confined ngwfs is used
    ! gcc32: make it block diagonal, as only Hconf_{\alpha\alpha} != 0
    if(pub_confined_ngwfs) then
       allocate(ham%confinement(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_allocate','ham%confinement',ierr)
    end if

    ! ndmh: in PAW, nonlocpot is density-dependent
    if (pub_aug) then
       allocate(ham%nonlocpot(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_allocate','ham%nonlocpot',ierr)
       allocate(ham%dijhat(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_allocate','ham%dijhat',ierr)
    else
       allocate(ham%nonlocpot(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_allocate','ham%nonlocpot',ierr)
       allocate(ham%dijhat(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%dijhat',ierr)
    end if

    ! qoh: Allocate Hartree-Fock exchange sparse matrix if necessary
    if (use_hfx) then
       allocate(ham%hfexchange(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%hfexchange',ierr)
    else
       allocate(ham%hfexchange(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%hfexchange',ierr)
    end if

    ! jd: Allocate polarisable embedding sparse matrix.
    !     NB: In the absence of polarisable embedding this still needs to be
    !         allocated (but is not sparse_created), because we're passing
    !         this around as an inout arg ('attempt to use pointer...').
    allocate(ham%polemb(2),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','ham%polemb',ierr)
    allocate(ham%vacuum_polemb(2),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','ham%vacuum_polemb',ierr)
    allocate(ham%vacuum_aux_polemb(2),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','ham%vacuum_aux_polemb',ierr)

    ! ddor: Allocate DFT+U Hamiltonian if necessary
    if (use_hubbard) then
       allocate(ham%hubbard_ham(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%hubbard_ham',ierr)
    else
       allocate(ham%hubbard_ham(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%hubbard_ham',ierr)
    endif

    ! lr408: Allocate conduction projected Hamiltonian if necessary -
    ! lr408: only needed if this is for the cond rep
    if (use_cond) then
       allocate(ham%cond_non_proj_ham(pub_num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%cond_non_proj_ham',ierr)
       ham%cond_shift = cond_init_shift
    else
       allocate(ham%cond_non_proj_ham(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','ham%cond_non_proj_ham',ierr)
       ham%cond_shift = 0.0_DP
    end if

    ! rc2013: Allocate projected embedding Hamiltonain if needed
    !if (pub_project_embed) then allocate(ham%proj_embed_ham(pub_num_spins), stat=ierr) call utils_alloc_check('ngwf_ham_create',
    !'ham%proj_embed_ham', ierr) allocate(ham%unproj_ham(pub_num_spins), stat=ierr) call utils_alloc_check('ngwf_ham_create',
    !'ham%proj_ham', ierr) endif

  end subroutine ngwf_ham_allocate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_ham_build_subsystem(ham_dst, ham_src, lhxc_nl_only)

    !=========================================================================!
    ! This subroutine copies the components of a NGWF_HAM type object to a    !
    ! subsystem matrix.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  ham_src (inout) : Container for Hamiltonian matrices to be copied into.!
    !  ham_dst (inout) : Container for Hamiltonian matrices to be copied from.!
    !-------------------------------------------------------------------------!
    ! Written by Robert Charlton, 26/04/2018.                                 !
    !=========================================================================!

    use rundat, only: pub_aug, pub_hubbard, pub_num_spins, pub_use_hfx, &
        pub_xc_ke_density_required, pub_project_embed, pub_pol_emb_pot, &
        pub_debug_on_root
    use sparse_embed, only: sparse_embed_create, sparse_embed_extract_sub

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout)     :: ham_dst(:)
    type(NGWF_HAM), intent(in)        :: ham_src
    logical, intent(in), optional           :: lhxc_nl_only

    ! Local Variables
    integer :: is, isub
    logical :: loc_use_hfx, loc_hubbard
    logical :: loc_cond

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering ngwf_ham_build_subsystem'

    ! Deal with optional argument (used to suppress hfx and hubbard parts)
    if (present(lhxc_nl_only)) then
       if (lhxc_nl_only) then
          loc_use_hfx = .false.
          loc_hubbard = .false.
       else
          loc_use_hfx = pub_use_hfx
          loc_hubbard = pub_hubbard
       end if
    else
       loc_use_hfx = pub_use_hfx
       loc_hubbard = pub_hubbard
    end if

    ! rc2013: EMBED_FIX
    !if (rep%postfix=='c') then
    !   loc_cond = .true.
    !else
       loc_cond = .false.
    !end if

    do isub=1,ham_src%ham(1)%mrows

       ! rc2013: allocate all matrices in the CROSS_HAM
       call ngwf_ham_allocate(ham_dst(isub),loc_use_hfx,loc_hubbard,loc_cond)

       ! rc2013: now let's create the matrices
       do is=1,pub_num_spins
          ! rc2013: use existing ham structures to form new matrices
          ! Cond structures are already dealt with
          call sparse_embed_create(ham_dst(isub)%ham(is), &
               ham_src%ham(is), arow=isub, bcol=isub)
          call sparse_embed_extract_sub(ham_dst(isub)%ham(is), &
               ham_src%ham(is), row=isub, col=isub)
          call sparse_embed_create(ham_dst(isub)%lhxc(is), &
               ham_src%lhxc(is), arow=isub, bcol=isub)
          call sparse_embed_extract_sub(ham_dst(isub)%lhxc(is), &
               ham_src%lhxc(is), row=isub, col=isub)
          if (pub_xc_ke_density_required) then
             call sparse_embed_create(ham_dst(isub)%dfdtau(is), &
                  ham_src%dfdtau(is), arow=isub, bcol=isub)
             call sparse_embed_extract_sub(ham_dst(isub)%dfdtau(is), &
                  ham_src%dfdtau(is), row=isub, col=isub)
          end if
          if (pub_aug) then
             call sparse_embed_create(ham_dst(isub)%nonlocpot(is), &
                  ham_src%nonlocpot(is), arow=isub, bcol=isub)
             call sparse_embed_extract_sub(ham_dst(isub)%nonlocpot(is), &
                  ham_src%nonlocpot(is), row=isub, col=isub)
             call sparse_embed_create(ham_dst(isub)%dijhat(is), &
                  ham_src%dijhat(is), arow=isub, bcol=isub)
             call sparse_embed_extract_sub(ham_dst(isub)%dijhat(is), &
                  ham_src%dijhat(is), row=isub, col=isub)
          end if
          if (loc_use_hfx) then
             call sparse_embed_create(ham_dst(isub)%hfexchange(is), &
                  ham_src%hfexchange(is), arow=isub, bcol=isub)
             call sparse_embed_extract_sub(ham_dst(isub)%hfexchange(is), &
                  ham_src%hfexchange(is), row=isub, col=isub)
          end if
          if (loc_hubbard) then
             call sparse_embed_create(ham_dst(isub)%hubbard_ham(is), &
                  ham_src%hubbard_ham(is), arow=isub, bcol=isub)
             call sparse_embed_extract_sub(ham_dst(isub)%hubbard_ham(is), &
                  ham_src%hubbard_ham(is), row=isub, col=isub)
          end if
          if (loc_cond) then
             call sparse_embed_create(ham_dst(isub)%cond_non_proj_ham(is), &
                  ham_src%hubbard_ham(is), arow=isub, bcol=isub)
             call sparse_embed_extract_sub(ham_dst(isub)%cond_non_proj_ham(is), &
                  ham_src%hubbard_ham(is), row=isub, col=isub)
          end if
          !if (pub_project_embed) then
          !   call sparse_embed_create(ham_dst(isub)%cond_non_proj_ham(is), &
          !        ham_src%hubbard_ham(is), arow=isub, bcol=isub)
          !   call sparse_embed_extract_sub(ham_dst(isub)%cond_non_proj_ham(is), &
          !        ham_src%hubbard_ham(is), row=isub, col=isub)
          !end if
       end do

       if (pub_pol_emb_pot) then
          call sparse_embed_create(ham_dst(isub)%polemb(1), &
               ham_src%polemb(1), arow=isub, bcol=isub)
          call sparse_embed_extract_sub(ham_dst(isub)%polemb(1), &
               ham_src%polemb(1), row=isub, col=isub)
          call sparse_embed_create(ham_dst(isub)%polemb(2), &
               ham_src%polemb(2), arow=isub, bcol=isub)
          call sparse_embed_extract_sub(ham_dst(isub)%polemb(2), &
               ham_src%polemb(2), row=isub, col=isub)
       end if

    end do
    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving ngwf_ham_build_subsystem'

  end subroutine ngwf_ham_build_subsystem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ngwf_rep_build_subsystem(rep_dst,rep_src)

    !=========================================================================!
    ! This subroutine copies the components of a NGWF_REP type object   !
    ! to a subsystem matrix.                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  rep_src (inout) : Container for Hamiltonian matrices to be copied into.!
    !  rep_dst (inout) : Container for Hamiltonian matrices to be copied from.!
    !-------------------------------------------------------------------------!
    ! Written by Robert Charlton, 03/05/2018.                                 !
    !=========================================================================!

    use datatypes, only: data_functions_alloc, data_functions_copy, &
         data_size
    use model_type, only: MODEL
    use rundat, only: pub_any_nl_proj, pub_paw, pub_aug, &
         pub_hubbard, pub_num_kpoints, pub_num_spins, pub_eda_scfmi_any, &
         pub_debug_on_root
    use sparse_embed, only: sparse_embed_create, sparse_embed_extract_sub
    use sparse_initialise, only: sparse_init_rep
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout) :: rep_dst(:)
    type(NGWF_REP), intent(in)    :: rep_src

    ! Local variables
    ! agrecocmplx
    integer :: ierr, isub

    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Entering ngwf_rep_build_subsystem'

    ! rc2013: build the various matrices for the subsystem rep
    do isub=1,rep_src%overlap%mrows
       ! set postfix
       rep_dst(isub)%postfix = rep_src%postfix
       ! rc2013: only 1 subsystem in this subsystem rep
       rep_dst(isub)%nsub = 1

       ! Overlap
       call sparse_embed_create(rep_dst(isub)%overlap, rep_src%overlap, &
            arow=isub, bcol=isub)
       call sparse_embed_extract_sub(rep_dst(isub)%overlap, rep_src%overlap, &
            row=isub, col=isub)
       if (pub_eda_scfmi_any) then
          call sparse_embed_create(rep_dst(isub)%overlap_scfmi_full, &
               rep_src%overlap_scfmi_full, arow=isub, bcol=isub)
          call sparse_embed_extract_sub(rep_dst(isub)%overlap_scfmi_full, &
               rep_src%overlap_scfmi_full, row=isub, col=isub)
       endif
       ! rc2013: structure already specified within rep_src, so
       ! no need to distinguish with pub_aug
       call sparse_embed_create(rep_dst(isub)%ngwf_overlap, &
            rep_src%ngwf_overlap, arow=isub, bcol=isub)
       call sparse_embed_extract_sub(rep_dst(isub)%ngwf_overlap, &
            rep_src%ngwf_overlap, row=isub, col=isub)
       ! Inverse overlap
       call sparse_embed_create(rep_dst(isub)%inv_overlap, &
            rep_src%inv_overlap, arow=isub, bcol=isub)
       call sparse_embed_extract_sub(rep_dst(isub)%inv_overlap, &
            rep_src%inv_overlap, row=isub, col=isub)
       ! SCF-MI full inverse overlap
       if (pub_eda_scfmi_any) then
          call sparse_embed_create(rep_dst(isub)%inv_overlap_scfmi_full, &
               rep_src%inv_overlap_scfmi_full, arow=isub, bcol=isub)
          call sparse_embed_extract_sub(rep_dst(isub)%inv_overlap_scfmi_full, &
               rep_src%inv_overlap_scfmi_full, row=isub, col=isub)
       endif
       ! NGWF-Projector overlap
       if (pub_any_nl_proj.or.pub_paw) then
          call sparse_embed_create(rep_dst(isub)%sp_overlap, &
               rep_src%sp_overlap, arow=isub, bcol=isub)
          call sparse_embed_extract_sub(rep_dst(isub)%sp_overlap, &
               rep_src%sp_overlap, row=isub, col=isub)
       endif
       ! Create the NGWF-projector overlap matrix if required
       if (pub_hubbard) then
          call sparse_embed_create(rep_dst(isub)%hub_overlap, &
               rep_src%hub_overlap, arow=isub, bcol=isub)
          call sparse_embed_extract_sub(rep_dst(isub)%hub_overlap, &
               rep_src%hub_overlap, row=isub, col=isub)
          call sparse_embed_create(rep_dst(isub)%hub_overlap_t, &
               rep_src%hub_overlap_t, arow=isub, bcol=isub)
          call sparse_embed_extract_sub(rep_dst(isub)%hub_overlap_t, &
               rep_src%hub_overlap_t, row=isub, col=isub)
          if(pub_aug) then
             call sparse_embed_create(rep_dst(isub)%hub_proj_paw_overlap, &
                  rep_src%hub_proj_paw_overlap, arow=isub, bcol=isub)
             call sparse_embed_extract_sub(rep_dst(isub)%hub_proj_paw_overlap, &
                  rep_src%hub_proj_paw_overlap, row=isub, col=isub)
             call sparse_embed_create(rep_dst(isub)%hub_ngwf_paw_overlap, &
                  rep_src%hub_ngwf_paw_overlap, arow=isub, bcol=isub)
             call sparse_embed_extract_sub(rep_dst(isub)%hub_ngwf_paw_overlap, &
                  rep_src%hub_ngwf_paw_overlap, row=isub, col=isub)
          endif
       endif
       ! If we have norm-conserving PPs and nonlocal projectors, create
       ! nonlocal potential matrix
       if (pub_any_nl_proj.and.(.not.pub_aug)) then
          call sparse_embed_create(rep_dst(isub)%nonlocpot, &
               rep_src%nonlocpot, arow=isub, bcol=isub)
          call sparse_embed_extract_sub(rep_dst(isub)%nonlocpot, &
               rep_src%nonlocpot, row=isub, col=isub)
       endif

       ! Kinetic matrix
       call sparse_embed_create(rep_dst(isub)%kinet, &
            rep_src%kinet, arow=isub, bcol=isub)
       call sparse_embed_extract_sub(rep_dst(isub)%kinet, &
            rep_src%kinet, row=isub, col=isub)
       ! lr408: If this is a conduction representation, create the cross overlap
       ! lr408: matrix
       ! rc2013: EMBED FIX
       if ((rep_src%postfix=='c').or.(rep_src%postfix=='j') &
            .or.(rep_src%postfix=='a')) then
          call sparse_embed_create(rep_dst(isub)%cross_overlap, &
               rep_src%cross_overlap, arow=isub, bcol=isub)
          call sparse_embed_extract_sub(rep_dst(isub)%cross_overlap, &
               rep_src%cross_overlap, row=isub, col=isub)
          call sparse_embed_create(rep_dst(isub)%cross_overlap_tr, &
               rep_src%cross_overlap_tr, arow=isub, bcol=isub)
          call sparse_embed_extract_sub(rep_dst(isub)%cross_overlap_tr, &
               rep_src%cross_overlap_tr, row=isub, col=isub)
       endif

       rep_dst(isub)%inv_overlap_init = .false.
       ! Allocate n_occ according to public numbers of spins and k-points.
       ! agreco: currently pub_num_kpoints=1 as a parameter, remember to
       ! merge pub_num_kpoints_temp and pub_num_kpoints when everything
       ! is working
       allocate(rep_dst(isub)%n_occ(pub_num_spins, pub_num_kpoints), stat=ierr)
       call utils_alloc_check('ngwf_rep_build_subsystem', 'rep_dst%n_occ', ierr)

       ! rc2013: now let's copy the NGWFs across too
       if( .not. allocated(rep_dst(isub)%ngwfs_on_grid)) then
          allocate(rep_dst(isub)%ngwfs_on_grid(1), stat=ierr)
          call utils_alloc_check('ngwf_rep_build_subsystem', 'rep_dst%ngwfs_on_grid', ierr)
       end if
       call data_functions_alloc(rep_dst(isub)%ngwfs_on_grid(1), &
            data_size(rep_src%ngwfs_on_grid(isub)), iscmplx=.false.)
       call data_functions_copy(rep_dst(isub)%ngwfs_on_grid(1), &
            rep_src%ngwfs_on_grid(isub))


    enddo
    if(pub_debug_on_root) write(stdout,'(a)') 'DEBUG: Leaving ngwf_rep_build_subsystem'

  end subroutine ngwf_rep_build_subsystem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ngwf_representation
