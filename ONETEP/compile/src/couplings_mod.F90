! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                    Couplings module for                        !
!                Configuration Interaction (CI)                  !
!                                                                !
!----------------------------------------------------------------!
! This module was created by David Turban in 2014.               !
! Updated in July 2015                                           !
! Modified by Robert Charlton for embedding structures, Sep 2018.!
!================================================================!

module couplings

  use constants, only: DP
  use rundat, only: pub_debug_on_root

  implicit none

  private

  ! public variables, used for reading in couplings block
  ! in rundat module
  character(len=4), allocatable  :: couplings_state_ids(:)
  character(len=5), allocatable  :: couplings_state_types(:)
  character(len=64), allocatable :: couplings_state_names(:)
  real(kind=DP), allocatable     :: couplings_state_energies(:)
  integer, allocatable           :: couplings_state_exc_idxs(:)
  integer :: couplings_nstates

  ! Subroutines
  public :: couplings_calculate
  public :: couplings_exit

  ! Global Variables
  public :: couplings_state_ids
  public :: couplings_state_names
  public :: couplings_state_types
  public :: couplings_state_energies
  public :: couplings_state_exc_idxs
  public :: couplings_nstates

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Modified by Andrea Greco on 28/06/2015 to allow use of complex NGWFs. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine couplings_calculate(ngwf_basis, rep, hub, mdl, hub_proj_basis, &
       cond_ngwf_basis, is_cmplx)

    use datatypes, only: FUNCTIONS
    use cdft_intermediate_cg, only: cdft_intermediate_restart_read, &
         cdft_intermediate_update_matrices
    use comms, only: pub_on_root, comms_barrier
    use constants, only: stdout, NORMAL, VERBOSE, UP, DN
    use dense, only: DEM, dense_create, dense_destroy, dense_eigensolve, &
         dense_product, dense_get_col, dense_convert, dense_put_element, &
         dense_write, dense_determinant, dense_put_col, dense_transpose, &
         dense_invert, dense_get_element, dense_scale, dense_axpy, &
         dense_copy
    use function_basis, only: FUNC_BASIS
    use hubbard_build, only: HUBBARD_MODEL, hubbard_model_init, &
         hubbard_projection_mtx, hubbard_ham_matrix, &
         hubbard_projector_update
    use function_ops, only: function_ops_brappd_ketppd
    use model_type, only: MODEL
    use ngwf_representation, only: NGWF_REP, ngwf_rep_create, &
         ngwf_rep_destroy
    use projectors, only: projectors_func_ovlp_box
    use restart, only: restart_ngwfs_tightbox_input, restart_kernel_read
    use rundat, only: pub_rootname, pub_eigensolver_orfac, &
         pub_eigensolver_abstol, pub_num_spins, pub_lr_tddft_calculate, &
         cond_num_states, pub_spin, pub_output_detail, pub_maxit_hotelling, &
         pub_max_resid_hotelling, pub_hub_proj_mixing, pub_hubbard_atomsolve, &
         pub_print_qc, pub_num_kpoints, PUB_1K
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_copy, &
         sparse_product, sparse_read, &
         sparse_hotelling_init, sparse_hotelling_invert
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort, &
         utils_devel_code, utils_assert, utils_qc_print

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: ngwf_basis
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis
    type(NGWF_REP), intent(inout) :: rep
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(MODEL), intent(inout) :: mdl
    ! agrecocmplx: optional argument for complex NGWFs
    logical, intent(in), optional :: is_cmplx

    ! Local Variables
    integer :: ierr                    ! memory allocation error flag
    integer :: is
    integer :: ii, jj, kk
    integer :: orbital_index
    type(NGWF_REP) :: cond_rep ! need to create conduction rep here
    type(FUNCTIONS), allocatable :: ngwfs_on_grid(:), cond_ngwfs_on_grid(:)
    ! sparse matrices
    type(SPAM3), allocatable   :: hubbard_ham(:)
    type(SPAM3)                :: hub_overlap_temp
    type(SPAM3)                :: dummy_sparse
    ! dense matrices
    type(DEM), allocatable     :: hubbard_ham_dense(:)
    ! overlaps of KS determinants + overlaps with CDFT potentials
    type(DEM)                  :: states_overlap, states_V_overlap
    ! eigenstates of CI Hamiltonian
    type(DEM)                  :: states_eigenvecs
    ! CI Hamitlonian + symmetrized version
    type(DEM)                  :: coupling_ham, coupling_ham_sym
    ! dense NGWF overlap, NGWF-COND cross overlap
    type(DEM)                  :: ngwf_olp_dense, cross_olp_dense
    ! auxillary matrices
    type(DEM)                  :: ks_orbs_occ_t, dummy_dense
    !type(DEM)                  :: up_matrix_dense, down_matrix_dense
    ! dense matrices for storing expansion coefficients of
    ! KS orbitals in NGWF basis, response kernels of
    ! TDDFT excitations in KS basis
    type(DEM), allocatable :: eigenvecs(:), ks_orbitals_occ(:,:)
    type(DEM), allocatable :: cond_ks_orbitals_occ(:,:)
    type(DEM), allocatable :: resp_dkns_ks(:)
    real(kind=DP), allocatable :: eigenvals(:,:)
    ! stores eigenenergies of CI Hamiltonian
    real(kind=DP), allocatable :: states_eigenvals(:)
    ! buffer rootname
    character(len=80)            :: rootname_buff
    ! dhpt: for test purposes
    !type(DEM) :: o_matrix_dense, hub_overlap_dense, &
    !             hub_overlap_t_dense
    !type(DEM), allocatable     :: occ_mat(:)
    !type(SPAM3), allocatable   :: occ_mat_sparse(:)
    ! agrecocmplx
    logical :: loc_cmplx

    ! cks: start timer
    call timer_clock("couplings_calculate",1)

    ! jme: KPOINTS_DANGER
    ! Some parts of this subroutine enforce a single k-point.
    call utils_assert(pub_num_kpoints == PUB_1K, &
         'Subroutine couplings_calculate not ready yet for more&
         & than one k-point.')

    ! agrecocmplx
    if (present(is_cmplx)) then
        loc_cmplx = is_cmplx
    else
        loc_cmplx = .false.
    end if

    call internal_initialise
    !write(stdout,'(a)') pub_rootname

    call internal_read_files_calc_eigen

    call internal_calc_overlaps

    call internal_build_ham

    call internal_write_matrices

    call internal_destroy

    ! cks: stop timer
    call timer_clock("couplings_calculate",2)

    contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      subroutine internal_initialise
        use datatypes, only: data_functions_alloc
        implicit none

        allocate(hubbard_ham(pub_num_spins),stat=ierr)
        call utils_alloc_check('couplings_calculate','hubbard_ham',ierr)
        allocate(hubbard_ham_dense(pub_num_spins),stat=ierr)
        call utils_alloc_check('couplings_calculate','hubbard_ham_dense',ierr)

        !allocate(ngwfs_on_grid(couplings_nstates,ngwf_basis%size_on_grid),stat=ierr)
        ! agrecocmplx: allocate using appropriate routines
        allocate(ngwfs_on_grid(couplings_nstates), stat=ierr)
        call utils_alloc_check('couplings_calculate','ngwfs_on_grid',ierr)
        do ii=1, couplings_nstates
            call data_functions_alloc(ngwfs_on_grid(ii), &
                 ngwf_basis%size_on_grid, iscmplx=loc_cmplx)
        end do
        allocate(resp_dkns_ks(couplings_nstates),stat=ierr)
        call utils_alloc_check('couplings_calculate','resp_dkns_ks',ierr)

        !do is=1,pub_num_spins
        !   hubbard_ham(is)%structure = 'H'
        !   call sparse_create(hubbard_ham(is), iscmplx=loc_cmplx)
        !   call dense_create(hubbard_ham_dense(is),&
        !        ngwf_basis%num,ngwf_basis%num, iscmplx_hubbard_ham(is)%iscmplx)
        !end do

        allocate(ks_orbitals_occ(couplings_nstates,pub_num_spins),stat=ierr)
        call utils_alloc_check('couplings_calculate','ks_orbitals_occ',ierr)
        do ii=1,couplings_nstates
          do is=1,pub_num_spins
            if (rep%n_occ(is,PUB_1K)>0) then
              call dense_create(ks_orbitals_occ(ii,is), ngwf_basis%num, &
                   rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
            end if
          end do
        end do

        if (pub_lr_tddft_calculate) then
           ! initialise conduction representation
           call ngwf_rep_create(cond_rep, 'c', mdl)
           ! set conduction occupation numbers
           cond_rep%n_occ(UP,PUB_1K) = (cond_num_states - pub_spin)/2
           cond_rep%n_occ(DN,PUB_1K) = (cond_num_states + pub_spin)/2
           !allocate(cond_ngwfs_on_grid(couplings_nstates,cond_ngwf_basis%size_on_grid),stat=ierr)
           ! agrecocmplx: allocate using appropriate routines
           allocate(cond_ngwfs_on_grid(couplings_nstates), stat=ierr)
           call utils_alloc_check('couplings_calculate','cond_ngwfs_on_grid',ierr)
           do ii=1, couplings_nstates
               call data_functions_alloc(cond_ngwfs_on_grid(ii), &
                    cond_ngwf_basis%size_on_grid, iscmplx=loc_cmplx)
           end do
           allocate(cond_ks_orbitals_occ(couplings_nstates,pub_num_spins),stat=ierr)
           call utils_alloc_check('couplings_calculate','cond_ks_orbitals_occ',ierr)
           do ii=1,couplings_nstates
             do is=1,pub_num_spins
               call dense_create(cond_ks_orbitals_occ(ii,is),cond_ngwf_basis%num, &
                    cond_rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
             end do
           end do
           ! create dense NGWF-COND cross overlap
           call dense_create(cross_olp_dense,ngwf_basis%num,cond_ngwf_basis%num, &
                iscmplx=loc_cmplx)
        end if

        ! create dense NGWF overlap
        call dense_create(ngwf_olp_dense,ngwf_basis%num,ngwf_basis%num, &
             iscmplx=loc_cmplx)

        call dense_create(states_overlap,couplings_nstates,couplings_nstates, &
             iscmplx=loc_cmplx)
        call dense_create(states_V_overlap,couplings_nstates,couplings_nstates, &
             iscmplx=loc_cmplx)
        call dense_create(coupling_ham,couplings_nstates,couplings_nstates, &
             iscmplx=loc_cmplx)
        call dense_create(coupling_ham_sym,couplings_nstates,couplings_nstates, &
             iscmplx=loc_cmplx)
        call dense_create(states_eigenvecs,couplings_nstates,couplings_nstates, &
             iscmplx=loc_cmplx)
        allocate(states_eigenvals(couplings_nstates),stat=ierr)
        call utils_alloc_check('couplings_calculate','states_eigenvals',ierr)

        write(rootname_buff,'(a)') trim(pub_rootname)

        !!!! CONSISTENCY CHECK !!!!
        !allocate(occ_mat(pub_num_spins),stat=ierr)
        !call utils_alloc_check('couplings_calculate','occ_mat',ierr)
        !allocate(occ_mat_sparse(pub_num_spins),stat=ierr)
        !call utils_alloc_check('couplings_calculate','occ_mat_sparse',ierr)
        !!!! CONSISTENCY CHECK !!!!
      end subroutine internal_initialise

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      subroutine internal_destroy
        use datatypes, only: data_functions_dealloc
        implicit none

        !!!! DEALLOCATE COND/TDDFT STRUCTURES !!!!
        ! cond_rep, cross_overlap, cond_ks_orbitals_occ, ...

        ! write original rootname back
        write(pub_rootname,'(a)') trim(rootname_buff)

        !do is=1,pub_num_spins
        !   call sparse_destroy(hubbard_ham(is))
        !   call dense_destroy(hubbard_ham_dense(is))
        !end do

        do ii=1,couplings_nstates
          do is=1,pub_num_spins
            if (rep%n_occ(is,PUB_1K)>0) then
              call dense_destroy(ks_orbitals_occ(ii,is))
            end if
          end do
          if (couplings_state_types(ii) == 'tddft') then
             call dense_destroy(resp_dkns_ks(ii))
          end if
        end do

        deallocate(ks_orbitals_occ,stat=ierr)
        call utils_dealloc_check('couplings_calculate','ks_orbitals_occ',ierr)
        deallocate(hubbard_ham,stat=ierr)
        call utils_dealloc_check('couplings_calculate','hubbard_ham',ierr)
        deallocate(hubbard_ham_dense,stat=ierr)
        call utils_dealloc_check('couplings_calculate','hubbard_ham_dense',ierr)
        ! agrecocmplx: deallocate using appropriate routine
        do ii=1,couplings_nstates
            call data_functions_dealloc(ngwfs_on_grid(ii))
        end do
        deallocate(ngwfs_on_grid,stat=ierr)
        call utils_dealloc_check('couplings_calculate','ngwfs_on_grid',ierr)
        deallocate(resp_dkns_ks,stat=ierr)
        call utils_dealloc_check('couplings_calculate','resp_dkns_ks',ierr)

        call dense_destroy(ngwf_olp_dense)

        call dense_destroy(states_overlap)
        call dense_destroy(states_V_overlap)
        call dense_destroy(coupling_ham)
        call dense_destroy(coupling_ham_sym)
        call dense_destroy(states_eigenvecs)
        deallocate(states_eigenvals,stat=ierr)
        call utils_dealloc_check('couplings_calculate','states_eigenvals',ierr)

        if (pub_lr_tddft_calculate) then
           call ngwf_rep_destroy(cond_rep)
           do ii=1,couplings_nstates
             do is=1,pub_num_spins
               call dense_destroy(cond_ks_orbitals_occ(ii,is))
             end do
           end do
           ! create dense NGWF-COND cross overlap
           call dense_destroy(cross_olp_dense)
           deallocate(cond_ngwfs_on_grid,stat=ierr)
           call utils_dealloc_check('couplings_calculate','cond_ngwfs_on_grid',ierr)
           deallocate(cond_ks_orbitals_occ,stat=ierr)
           call utils_dealloc_check('couplings_calculate','cond_ks_orbitals_occ',ierr)
        end if

      end subroutine internal_destroy
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      subroutine internal_read_files_calc_eigen

        use sparse_embed, only: SPAM3_EMBED_ARRAY, sparse_embed_array_create, &
             sparse_embed_array_destroy
        implicit none

        character(len=80) :: filename
        logical :: fileexists
        type(SPAM3_EMBED_ARRAY)    :: denskern
        type(SPAM3)                :: resp_denskern
        type(DEM), allocatable     :: dkn_dense(:)
        type(DEM)                  :: resp_dkn_dense
        type(DEM)                  :: unity ! dense unit matrix
        real(kind=DP), allocatable, dimension(:) :: col_dummy
        character(len=10) :: structure

        ! initialise
        allocate(eigenvecs(pub_num_spins),stat=ierr)
        call utils_alloc_check('couplings_calculate','eigenvecs',ierr)
        allocate(dkn_dense(pub_num_spins),stat=ierr)
        call utils_alloc_check('couplings_calculate','dkn_dense',ierr)
        allocate(eigenvals(ngwf_basis%num,pub_num_spins),stat=ierr)
        call utils_alloc_check('couplings_calculate','eigenvals',ierr)


        ! read valence states and construct KS orbitals
        allocate(col_dummy(ngwf_basis%num),stat=ierr)
        call utils_alloc_check('couplings_calculate','col_dummy',ierr)
        structure = 'K'
        call sparse_embed_array_create(denskern, n_spins=pub_num_spins, &
             n_kpoints=PUB_1K, structure=structure, &
             iscmplx=loc_cmplx)
        do is=1,pub_num_spins
           call dense_create(eigenvecs(is),ngwf_basis%num,ngwf_basis%num, &
                iscmplx=loc_cmplx)
           ! create dense version of denskern
           call dense_create(dkn_dense(is),ngwf_basis%num,ngwf_basis%num, &
                iscmplx=loc_cmplx)
        end do

        ! loop over states to couple
        do ii=1,couplings_nstates
           ! overwrite rootname with name of current state
           write(pub_rootname,'(a)') trim(couplings_state_names(ii))
           ! read ngwfs for this state in tightbox format
           call restart_ngwfs_tightbox_input(ngwfs_on_grid(ii), ngwf_basis, &
                mdl%cell,mdl%fftbox,mdl%elements, 'tightbox_ngwfs',mdl%regions(1))
           ! read density kernel
           call restart_kernel_read(denskern)
           ! necessary?
           call comms_barrier
           !call services_flush
           ! calculate NGWF overlap
           call function_ops_brappd_ketppd(rep%overlap%p, &
                ngwfs_on_grid(ii), ngwf_basis, ngwfs_on_grid(ii), &
                ngwf_basis, mdl%cell)
           ! convert into dense format
           call dense_convert(ngwf_olp_dense, rep%overlap%p)

           do is=1,pub_num_spins
              if (rep%n_occ(is,PUB_1K)==0) cycle
              ! calculate eigenstates of kernel for each spin -> KS orbitals
              call dense_convert(dkn_dense(is),denskern%m(is,PUB_1K))
              ! initialise unity
              call dense_create(unity,ngwf_basis%num,ngwf_basis%num, &
                   iscmplx=loc_cmplx)
              do kk=1,ngwf_basis%num
                 call dense_put_element(1.0_DP,unity,kk,kk)
              end do
              ! WARNING: assuming that eigenvecs are normalised to 1
              ! is this the case?? -> YES
              call dense_eigensolve(ngwf_basis%num,eigenvals(:,is),dkn_dense(is), &
                   unity,1,eigenvecs(is),pub_eigensolver_orfac, &
                   pub_eigensolver_abstol)
              ! destory unity
              call dense_destroy(unity)
              !call dense_eigensolve(ngwf_basis%num,eigenvals(:,is),dkn_dense(is), &
              !     ngwf_olp_dense, 2,eigenvecs(is))
              ! select occupied KS orbitals and store in ks_orbitals_occ
              do orbital_index=ngwf_basis%num-rep%n_occ(is,PUB_1K)+1,ngwf_basis%num
                 call dense_get_col(col_dummy,eigenvecs(is),orbital_index)
                 ! scale eigenvectors with sqrt eigenvalue
                 col_dummy = col_dummy*sqrt(eigenvals(orbital_index,is))
                 call dense_put_col(col_dummy,ks_orbitals_occ(ii,is),orbital_index-&
                                    &ngwf_basis%num+rep%n_occ(is,PUB_1K))
              end do
           end do
        end do

        ! read conduction states and construct KS orbitals
        if (pub_lr_tddft_calculate) then

           call sparse_embed_array_destroy(denskern)
           do is=1,pub_num_spins
              call dense_destroy(eigenvecs(is))
              call dense_destroy(dkn_dense(is))
           end do
           deallocate(eigenvals,stat=ierr)
           call utils_dealloc_check('couplings_calculate','eigenvals',ierr)
           deallocate(col_dummy,stat=ierr)
           call utils_dealloc_check('couplings_calculate','col_dummy',ierr)
           allocate(col_dummy(cond_ngwf_basis%num),stat=ierr)
           call utils_alloc_check('couplings_calculate','col_dummy',ierr)
           allocate(eigenvals(cond_ngwf_basis%num,pub_num_spins),stat=ierr)
           call utils_alloc_check('couplings_calculate','eigenvals',ierr)
           call sparse_embed_array_create(denskern, n_spins=pub_num_spins, &
                n_kpoints=PUB_1K, structure='K'//cond_rep%postfix, &
                iscmplx=loc_cmplx)
           do is=1,pub_num_spins
              call dense_create(eigenvecs(is),cond_ngwf_basis%num,cond_ngwf_basis%num, &
                   iscmplx=loc_cmplx)
              ! create dense version of denskern
              call dense_create(dkn_dense(is),cond_ngwf_basis%num,cond_ngwf_basis%num, &
                   iscmplx=loc_cmplx)
           end do

           ! assuming not using joint set (->check!!)
           ! row index: COND, column index: NGWFs
           resp_denskern%structure = 'TDR'//cond_rep%postfix
           call sparse_create(resp_denskern, iscmplx=loc_cmplx)
           call dense_create(resp_dkn_dense,cond_ngwf_basis%num,ngwf_basis%num, &
                iscmplx=loc_cmplx)

           ! loop over states to couple
           do ii=1,couplings_nstates
              ! overwrite rootname with name of current state
              write(pub_rootname,'(a)') trim(couplings_state_names(ii))

              ! read conduction NGWFs, conduction kernel and response kernel in case
              ! of tddft state
              ! and convert into KS orbital representation
              if (couplings_state_types(ii) == 'tddft') then
                 ! read conduction NGWFs
                 call restart_ngwfs_tightbox_input(cond_ngwfs_on_grid(ii), &
                       cond_ngwf_basis, mdl%cell,mdl%fftbox, mdl%elements, &
                       'tightbox_ngwfs_cond',mdl%regions(1))
                 ! read conduction density kernel
                 call restart_kernel_read(denskern, read_cond=.true.)
                 ! convert to dense matrices
                 do is=1,pub_num_spins
                    call dense_convert(dkn_dense(is),denskern%m(is,PUB_1K))
                 end do
                 ! calculate eigenstates -> obtain conduction KS orbitals
                 do is=1,pub_num_spins
                    ! calculate eigenstates of kernel for each spin -> KS orbitals
                    ! re-init unity
                    call dense_create(unity,cond_ngwf_basis%num,cond_ngwf_basis%num, &
                         iscmplx=loc_cmplx)
                    do kk=1,cond_ngwf_basis%num
                       call dense_put_element(1.0_DP,unity,kk,kk)
                    end do
                    call dense_eigensolve(cond_ngwf_basis%num,eigenvals(:,is),dkn_dense(is), &
                         unity,1,eigenvecs(is),pub_eigensolver_orfac, &
                         pub_eigensolver_abstol)
                    ! destroy unit matrix
                    call dense_destroy(unity)
                    ! select 'occupied' KS orbitals and store in cond_ks_orbitals_occ
                    do orbital_index = &
                         cond_ngwf_basis%num-cond_rep%n_occ(is,PUB_1K)+1, &
                         cond_ngwf_basis%num
                       call dense_get_col(col_dummy,eigenvecs(is),orbital_index)
                       !do jj=1,ngwf_basis%num
                       !   col_dummy(jj) = col_dummy(jj)*sqrt(eigenvals(orbital_index,is))
                       !end do
                       col_dummy = col_dummy*sqrt(eigenvals(orbital_index,is))
                       call dense_put_col(col_dummy,cond_ks_orbitals_occ(ii,is),orbital_index-&
                                          &cond_ngwf_basis%num+cond_rep%n_occ(is,PUB_1K))
                    end do
                 end do

                 ! read response densisty kernel
                 write(filename,'(2a,i0,a)') trim(pub_rootname), &
                      '_response_denskern_', couplings_state_exc_idxs(ii), &
                      '.dkn'
#ifdef HDF5
                 filename=trim(filename)//'.h5'
#endif

                 if (pub_on_root) write(stdout,'(/3a)',advance='no') &
                      'Reading density kernel from file "', trim(filename),'" ...'

                 if (pub_on_root) then
                    inquire(file=filename,exist=fileexists)
                 else
                    fileexists = .true.
                 end if

                 if (fileexists) then

                    ! Read density kernel from this file
                    call sparse_read(resp_denskern,trim(filename))

                 else

                    call utils_abort('File "'//trim(filename)//'" not found in &
                         &restart_response_kernel_batch_read().' // &
                         'If the file that ONETEP was trying to read had been &
                         &written by a version of ONETEP < v5.0, this error &
                         &may be due to a discrepancy in the name of the file, &
                         &the actual file name having too many spaces for the &
                         &current version of ONETEP; if this is the case, &
                         &please manually rename the file appropriately and &
                         &run ONETEP again.')

                 endif

                 ! convert into KS representation and store in resp_dkns_ks
                 ! WATCH OUT: row index: COND, column index: NGWFs
                 call dense_convert(resp_dkn_dense,resp_denskern)
                 ! NOTE: this only works for closed shell systems
                 call dense_create(ks_orbs_occ_t, &
                      cond_rep%n_occ(1,PUB_1K),cond_ngwf_basis%num, &
                      iscmplx=loc_cmplx)
                 call dense_create(dummy_dense, &
                      cond_ngwf_basis%num, rep%n_occ(1,PUB_1K), &
                      iscmplx=loc_cmplx)
                 call dense_create(resp_dkns_ks(ii), &
                      cond_rep%n_occ(1,PUB_1K), &
                      rep%n_occ(1,PUB_1K), iscmplx=loc_cmplx)
                 call dense_transpose(ks_orbs_occ_t,cond_ks_orbitals_occ(ii,1))
                 call dense_product(dummy_dense,resp_dkn_dense, &
                                    ks_orbitals_occ(ii,1))
                 call dense_product(resp_dkns_ks(ii),ks_orbs_occ_t,dummy_dense)
                 call dense_destroy(ks_orbs_occ_t)
                 call dense_destroy(dummy_dense)
              end if
           end do

        end if

        ! deallocate
        do is=1,pub_num_spins
           call dense_destroy(eigenvecs(is))
           call dense_destroy(dkn_dense(is))
        end do
        deallocate(eigenvecs,stat=ierr)
        call utils_dealloc_check('couplings_calculate','eigenvecs',ierr)
        deallocate(eigenvals,stat=ierr)
        call utils_dealloc_check('couplings_calculate','eigenvals',ierr)
        call sparse_embed_array_destroy(denskern)
        deallocate(dkn_dense,stat=ierr)
        call utils_dealloc_check('couplings_calculate','dkn_dense',ierr)
        deallocate(col_dummy,stat=ierr)
        call utils_dealloc_check('couplings_calculate','col_dummy',ierr)

        if (pub_lr_tddft_calculate) then
           call sparse_destroy(resp_denskern)
           call dense_destroy(resp_dkn_dense)
        end if

      end subroutine internal_read_files_calc_eigen


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      subroutine internal_calc_overlaps

      use rundat, only: pub_num_kpoints, PUB_1K
      use projectors, only: projectors_func_ovlp_box, PROJECTOR_SET
      implicit none

      ! auxillary variables
      integer        :: A, B
      logical        :: V_flag

      real(kind=DP)  :: overlap, V_overlap, det, dummy
      real(kind=DP)  :: det_S_vv, tr_VSinv
      real(kind=DP)  :: tr_TRtSinv, tr_VSinvTRtSinv, tr_WRtSinv
      type(DEM) :: S_vv, S_vv_inv, T_vc, V_vv, W_vc, Rt
      type(DEM) :: VSinv, TRtSinv, WRtSinv

      ! needed fo hubbard_projector_update call
      type(FUNC_BASIS)    :: projector_basis
      type(PROJECTOR_SET) :: nl_projectors
      !type(SPAM3)         :: sp_overlap
      !type(SPAM3)         :: hub_proj_paw_overlap
      !type(SPAM3)         :: hub_ngwf_paw_overlap
      ! required for calculating inverse NGWF overlap
      type(DEM) :: inv_overlap_dens

      ! jme: KPOINTS_DANGER
      ! Some parts of this subroutine enforce a single k-point.
      call utils_assert(pub_num_kpoints == PUB_1K, &
           'Subroutine internal_calc_overlaps not ready yet for more&
           & than one k-point.')


      do ii=1,couplings_nstates

         if (couplings_state_types(ii) == 'cdft') then
            ! read cdft potentials
            write(pub_rootname,'(a)') trim(couplings_state_names(ii))
            call cdft_intermediate_restart_read(hub,mdl)
            call cdft_intermediate_update_matrices(hub, hub_proj_basis)
            ! prepare reading cdft projectors from file
            pub_hub_proj_mixing = -1.0_DP
            hub%consistency_iteration = 1
            pub_hubbard_atomsolve = .false.
            ! calculate NGWF overlap
            call function_ops_brappd_ketppd(rep%overlap%p, &
                 ngwfs_on_grid(ii), ngwf_basis, ngwfs_on_grid(ii), &
                 ngwf_basis, mdl%cell)

            ! copied for hamiltonian_dens_indep_matrices

            ! ======================== INVERSE OVERLAP MATRIX ================

            ! cks: initialise inverse overlap guess if Hotelling recursion
            ! cks: will be used to approximate it.
            ! ndmh: call sparse_mod version of hotelling's algorithm initialisation
            if (pub_maxit_hotelling > 0) call sparse_hotelling_init(rep%inv_overlap%p, &
                 rep%overlap%p)

            ! cks : approximate inverse overlap by Hotelling's recursion
            if (pub_maxit_hotelling > 0) then

               if (pub_output_detail>=VERBOSE .and. pub_on_root) then
                  write(stdout,'(a)')'============ Calculation of NGWF S^-1 using &
                       &Hotelling algorithm ================ '
                  write(stdout,'(a)')'         Iteration  |  Frobenius norm and Maximum &
                       &abs value    (of I-S*S_n^-1)   '
               end if

               ! ndmh: call sparse_mod version of hotelling's algorithm
               call sparse_hotelling_invert(rep%inv_overlap%p, & !inout
                    rep%overlap%p,show_output=(pub_output_detail>=VERBOSE), &
                    max_resid_converged=pub_max_resid_hotelling, &
                    num_iter=pub_maxit_hotelling)!,final_max_resid=final_max_resid)

               if (pub_output_detail>=VERBOSE.and.pub_on_root) then
                  write(stdout,'(a)')'===================================&
                       &============================================='
               end if

            else if (pub_maxit_hotelling == 0) then

               ! ndmh: allocate storage for dense inverse
               call dense_create(inv_overlap_dens, ngwf_basis%num, &
                    ngwf_basis%num, iscmplx=loc_cmplx)

               ! ndmh: copy current overlap matrix into dense inverse overlap
               call dense_convert(inv_overlap_dens, rep%overlap%p)

               ! ndmh: invert dense matrix
               call dense_invert(inv_overlap_dens)

               ! ndmh: copy dense inverse overlap back to sparse matrix
               call dense_convert(rep%inv_overlap%p,inv_overlap_dens)

               ! ndmh: deallocate storage for dense inverse
               call dense_destroy(inv_overlap_dens)

            else
               call utils_abort('Error in couplings_calculate: &
                    &negative maxit_hotelling ',pub_maxit_hotelling)
            end if

            ! ==================== END INVERSE OVERLAP MATRIX =================
            call hubbard_projector_update( &
                                    ngwfs_on_grid(ii), &        !inout
                                    ngwf_basis, &               !in
                                    nl_projectors, &            !paw
                                    projector_basis, &          !paw
                                    hub_proj_basis, &           !inout
                                    hub, &                      !inout
                                    rep%inv_overlap%p, &          !in
                                    rep%hub_overlap%p, &          !inout
                                    rep%hub_overlap_t%p, &        !inout
                                    rep%sp_overlap%p, &           !paw
                                    rep%hub_proj_paw_overlap%p, & !paw
                                    rep%hub_ngwf_paw_overlap%p, & !paw
                                    mdl=mdl)
            ! this is now taken care of by hubbard_projector_update
            ! projector overlap wit ii ngwfs
            !call function_ops_brappd_ketppd(rep%hub_overlap, &         !input-output
            !     ngwfs_on_grid(ii), ngwf_basis, hub%consistency_projs, & !input
            !     hub_proj_basis,mdl%cell)                                  !input
            V_flag = .true.
         else
            V_flag = .false.
         end if

         do jj=1,couplings_nstates
            if (pub_on_root) then
               write(stdout,'(a)',advance='no') 'Calculating matrix element '
               write(stdout,'(i3)',advance='no') (ii-1)*couplings_nstates+jj
               write(stdout,'(a)',advance='no') ' of '
               write(stdout,'(i3)') couplings_nstates*couplings_nstates
            end if

            if (V_flag) then
               A = ii
               B = jj
            else
               A = jj
               B = ii
            end if

            ! calculate NGWF overlap
            call function_ops_brappd_ketppd(rep%overlap%p, &
                 ngwfs_on_grid(A), ngwf_basis, ngwfs_on_grid(B), &
                 ngwf_basis,mdl%cell)
            ! convert into dense format
            call dense_convert(ngwf_olp_dense, rep%overlap%p)


            ! Case-sensitive part
            if((couplings_state_types(A)=='tddft').and. &
               (couplings_state_types(B)=='tddft')) then

               if(A==B) then
                  overlap = 1.0_DP
               else
                  overlap = 0.0_DP
                  if (pub_on_root) write(stdout,'(a)') 'WARNING: Cannot couple &
                     &TDDFT states; setting matrix element to zero.'
               end if
               V_overlap = 0.0_DP

            else

               det_S_vv = 1.0_DP

               tr_VSinv        = 0.0_DP
               tr_TRtSinv      = 0.0_DP
               tr_WRtSinv      = 0.0_DP
               tr_VSinvTRtSinv = 0.0_DP

               if (V_flag) then
                  ! in this case necessarily A==ii, B==jj
                  ! projector overlap with jj ngwfs
                  call sparse_create(hub_overlap_temp,rep%hub_overlap%p)
                  call function_ops_brappd_ketppd(hub_overlap_temp, &      !input-output
                       ngwfs_on_grid(B), ngwf_basis, hub%consistency_projs, & !input
                       hub_proj_basis, mdl%cell)                                     !input
                  call hubbard_projection_mtx(.true.,.false., hub_overlap_temp, &
                        rep%hub_overlap_t%p, hub%o_matrix)
                  call sparse_destroy(hub_overlap_temp)

                  do is=1,pub_num_spins
                     if (is == 2) then
                        call sparse_copy(hub%projector_ham(is),hub%down_matrix)
                     else
                        call sparse_copy(hub%projector_ham(is),hub%up_matrix)
                     end if
                  end do

                  ! build hubbard_ham in NGWF ii, jj basis
                  !call hubbard_ham_matrix(hub,hubbard_ham,rep%hub_overlap,rep%hub_overlap_t)
                  do is = 1, pub_num_spins
                     ! create buffer
                     call sparse_create(dummy_sparse,rep%hub_overlap%p,&
                          hub%projector_ham(is))
                     call sparse_product(dummy_sparse,&
                          rep%hub_overlap%p,hub%projector_ham(is))
                     call sparse_create(hubbard_ham(is),dummy_sparse,&
                                        rep%hub_overlap_t%p)
                     call sparse_product(hubbard_ham(is),&
                          dummy_sparse,rep%hub_overlap_t%p)
                     call sparse_destroy(dummy_sparse)
                  end do
               end if

               do is=1,pub_num_spins
                  if (rep%n_occ(is,PUB_1K)==0) then
                     if (V_flag) call sparse_destroy(hubbard_ham(is))
                     cycle
                  endif
                  ! calculate overlap matrix of KS orbitals
                  ! corresponding to states A, B
                  call dense_create(ks_orbs_occ_t, &
                       rep%n_occ(is,PUB_1K),ngwf_basis%num, &
                       iscmplx=loc_cmplx)
                  call dense_create(dummy_dense, ngwf_basis%num, &
                       rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
                  call dense_create(S_vv, rep%n_occ(is,PUB_1K), &
                       rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
                  call dense_transpose(ks_orbs_occ_t,ks_orbitals_occ(A,is))
                  call dense_product(dummy_dense,ngwf_olp_dense, &
                                     ks_orbitals_occ(B,is))
                  call dense_product(S_vv,ks_orbs_occ_t,dummy_dense)
                  call dense_destroy(dummy_dense)
                  ! calculate determinants
                  call dense_determinant(S_vv,det)
                  ! update overlap
                  det_S_vv = det_S_vv*det

                  call dense_create(V_vv,rep%n_occ(is,PUB_1K), &
                       rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
                  ! if ii is cdft caculate matrix elements V_vv of potential
                  ! else leave it all 0
                  if (V_flag) then
                     ! calculate KS A, B matrix elements of CDFT potential A==ii
                     call dense_create(hubbard_ham_dense(is),&
                          ngwf_basis%num,ngwf_basis%num, iscmplx=loc_cmplx)
                     call dense_convert(hubbard_ham_dense(is),hubbard_ham(is))
                     call dense_create(dummy_dense, &
                          ngwf_basis%num, rep%n_occ(is,PUB_1K), &
                          iscmplx=loc_cmplx)
                     call dense_product(dummy_dense,hubbard_ham_dense(is), &
                                        ks_orbitals_occ(B,is))
                     call dense_product(V_vv,ks_orbs_occ_t,dummy_dense)
                     call dense_destroy(dummy_dense)
                     call dense_destroy(hubbard_ham_dense(is))
                     call sparse_destroy(hubbard_ham(is))
                  end if

                  ! calculate inverse overlap
                  call dense_create(S_vv_inv,S_vv, iscmplx=loc_cmplx)
                  call dense_copy(S_vv_inv,S_vv)
                  call dense_invert(S_vv_inv)
                  call dense_create(VSinv, rep%n_occ(is,PUB_1K), &
                       rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
                  ! calculate VSinv
                  call dense_product(VSinv,V_vv,S_vv_inv)
                  ! calculate trace
                  do kk=1,rep%n_occ(is,PUB_1K)
                     call dense_get_element(dummy,VSinv,kk,kk)
                     tr_VSinv = tr_VSinv + dummy
                  end do

                  ! coupling to tddft state
                  if (couplings_state_types(B)=='tddft') then

                     ! cross overlap (ngwf-cond) T_vc
                     call function_ops_brappd_ketppd(cond_rep%cross_overlap%p, &
                          ngwfs_on_grid(A),ngwf_basis,cond_ngwfs_on_grid(B),&
                          cond_ngwf_basis, mdl%cell)
                     call dense_convert(cross_olp_dense,cond_rep%cross_overlap%p)

                     call dense_create(dummy_dense,ngwf_basis%num, &
                          cond_rep%n_occ(is,PUB_NUM_SPINS), iscmplx=loc_cmplx)
                     call dense_product(dummy_dense,cross_olp_dense, &
                                        cond_ks_orbitals_occ(B,is))
                     call dense_create(T_vc,rep%n_occ(is,PUB_1K), &
                          cond_rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
                     call dense_product(T_vc,ks_orbs_occ_t,dummy_dense)
                     call dense_destroy(dummy_dense)

                     ! calculate W_vc
                     call dense_create( W_vc,rep%n_occ(is,PUB_1K), &
                          cond_rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
                     ! Calculate if potential overlap included
                     ! leave W_vc = 0 otherwise
                     if (V_flag) then
                        !! CONSTRUCT HUBBARD_HAM IN NGWF-COND BASIS !!
                        call sparse_create(hub_overlap_temp,cond_rep%hub_overlap%p)
                        call function_ops_brappd_ketppd(hub_overlap_temp, &      !input-output
                             cond_ngwfs_on_grid(B), cond_ngwf_basis, &
                             hub%consistency_projs, hub_proj_basis, mdl%cell)    !input
                        call hubbard_projection_mtx(.true.,.false., hub_overlap_temp, &
                              cond_rep%hub_overlap_t%p, hub%o_matrix)
                        call sparse_destroy(hub_overlap_temp)

                        call sparse_create(dummy_sparse,rep%hub_overlap%p,hub%projector_ham(is))
                        call sparse_product(dummy_sparse,&
                             rep%hub_overlap%p,hub%projector_ham(is))
                        call sparse_create(hubbard_ham(is),dummy_sparse,&
                                           cond_rep%hub_overlap_t%p)
                        call sparse_product(hubbard_ham(is),&
                             dummy_sparse,cond_rep%hub_overlap_t%p)
                        call sparse_destroy(dummy_sparse)

                        call dense_create(hubbard_ham_dense(is), &
                             ngwf_basis%num, cond_ngwf_basis%num, &
                             iscmplx=loc_cmplx)
                        call dense_convert(hubbard_ham_dense(is),hubbard_ham(is))
                        call dense_create(dummy_dense, ngwf_basis%num, &
                             cond_rep%n_occ(is,PUB_NUM_SPINS), iscmplx=loc_cmplx)
                        call dense_product(dummy_dense,hubbard_ham_dense(is), &
                                           cond_ks_orbitals_occ(B,is))
                        call dense_product(W_vc,ks_orbs_occ_t,dummy_dense)
                        call dense_destroy(dummy_dense)
                        call dense_destroy(hubbard_ham_dense(is))
                        call sparse_destroy(hubbard_ham(is))
                     end if

                     ! transpose of response kernel
                     ! NOTE: already cond-ngwf(!)
                     ! NOTE: assuming closed shell
                     call dense_create(Rt, cond_rep%n_occ(is,PUB_1K), &
                          rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
                     call dense_copy(Rt,resp_dkns_ks(B))

                     call dense_create( dummy_dense,cond_rep%n_occ(is,PUB_1K), &
                          rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
                     call dense_product(dummy_dense,Rt,S_vv_inv)
                     call dense_create( TRtSinv,rep%n_occ(is,PUB_1K), &
                          rep%n_occ(is,PUB_1K), iscmplx=loc_cmplx)
                     call dense_product(TRtSinv,T_vc,dummy_dense)
                     call dense_create( WRtSinv,TRtSinv, iscmplx=loc_cmplx)
                     call dense_product(WRtSinv,T_vc,dummy_dense)
                     call dense_destroy(dummy_dense)

                     call dense_create( dummy_dense,VSinv, iscmplx=loc_cmplx)
                     call dense_product(dummy_dense,VSinv,TRtSinv)

                     ! calculate traces
                     do kk=1,rep%n_occ(is,PUB_1K)
                        call dense_get_element(dummy,TRtSinv,kk,kk)
                        tr_TRtSinv = tr_TRtSinv + dummy
                        call dense_get_element(dummy,WRtSinv,kk,kk)
                        tr_WRtSinv = tr_WRtSinv + dummy
                        call dense_get_element(dummy,dummy_dense,kk,kk)
                        tr_VSinvTRtSinv = tr_VSinvTRtSinv + dummy
                     end do

                     call dense_destroy(dummy_dense)
                     call dense_destroy(T_vc)
                     call dense_destroy(W_vc)
                     call dense_destroy(Rt)
                     call dense_destroy(TRtSinv)
                     call dense_destroy(WRtSinv)
                  end if

                  ! destroy matrices
                  call dense_destroy(ks_orbs_occ_t)
                  call dense_destroy(S_vv)
                  call dense_destroy(S_vv_inv)
                  call dense_destroy(V_vv)
                  call dense_destroy(VSinv)

               end do

               ! calculate overlaps
               if (couplings_state_types(B)=='tddft') then
                  overlap   = det_S_vv*tr_TRtSinv
                  V_overlap = det_S_vv*(   tr_VSinv*tr_TRtSinv &
                                       & - tr_VSinvTRtSinv     &
                                       & + tr_WRtSinv)
               else
                  overlap   = det_S_vv
                  V_overlap = det_S_vv*tr_VSinv
               end if

            endif
            ! write states overlap matrices
            call dense_put_element(overlap,states_overlap,ii,jj)
            call dense_put_element(V_overlap,states_V_overlap,ii,jj)
         end do
      end do


      end subroutine internal_calc_overlaps

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      subroutine internal_build_ham

        real(kind=DP) :: dummy_elm, dummy_ovlp, dummy_V_ovlp, dummy_VN

        do ii=1,couplings_nstates
          do jj=1,couplings_nstates
            call dense_get_element(dummy_ovlp,states_overlap,ii,jj)
            call dense_get_element(dummy_V_ovlp,states_V_overlap,ii,jj)
            call dense_get_element(dummy_VN,states_V_overlap,ii,ii)
            dummy_elm = (couplings_state_energies(ii)+dummy_VN)*dummy_ovlp
            dummy_elm = dummy_elm - dummy_V_ovlp
            call dense_put_element(dummy_elm,coupling_ham,ii,jj)
          end do
        end do

        ! symmetrize
        call dense_create(dummy_dense,couplings_nstates,couplings_nstates, &
             iscmplx=loc_cmplx)
        call dense_transpose(dummy_dense,coupling_ham)
        call dense_copy(coupling_ham_sym,coupling_ham)
        call dense_axpy(coupling_ham_sym,dummy_dense,1.0_DP)
        call dense_scale(coupling_ham_sym,0.5_DP)
        call dense_destroy(dummy_dense)

      end subroutine internal_build_ham

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      subroutine internal_write_matrices

        logical :: tostdout
        character(len=12) :: ii_string

        if (pub_output_detail > NORMAL) then
           tostdout = .true.
        else
           tostdout = .false.
        end if

        if (pub_on_root) then
           if(tostdout) write(stdout,'(a)') 'CI Hamiltonian:'
           call dense_write(coupling_ham,trim(rootname_buff)//'_ci_ham', &
                to_stdout=tostdout)
           if(tostdout) write(stdout,'(a)') 'Symmetrised CI Hamiltonian:'
           call dense_write(coupling_ham_sym,trim(rootname_buff)//'_ci_ham_sym', &
                to_stdout=tostdout)
        end if

        ! calculate eigenstates of coupling Hamiltonian
        call dense_eigensolve(couplings_nstates,states_eigenvals,coupling_ham_sym, &
             states_overlap,1,states_eigenvecs,pub_eigensolver_orfac, &
             pub_eigensolver_abstol)

        call dense_create(dummy_dense,couplings_nstates,1, iscmplx=loc_cmplx)
        call dense_put_col(states_eigenvals,dummy_dense,1)
        if (pub_on_root) then
           if(tostdout) write(stdout,'(a)') 'Eigenvalues of symmetrised CI Hamiltonian:'
           call dense_write(dummy_dense,trim(rootname_buff)//'_ci_eigvals', &
                to_stdout=tostdout)
           if(tostdout) write(stdout,'(a)') 'Eigenvectors of symmetrised CI Hamiltonian (columns):'
           call dense_write(states_eigenvecs,trim(rootname_buff)//'_ci_eigstates', &
                to_stdout=tostdout)
           ! print eigenvalues as a QC test
           if (pub_print_qc) then
              do ii=1,couplings_nstates
                 write(ii_string,'(i10)') ii
                 ii_string = trim(adjustl(ii_string))
                 call utils_qc_print('eigenvalue_'//trim(ii_string), &
                       states_eigenvals(ii))
              end do
           end if
        end if
        call dense_destroy(dummy_dense)

      end subroutine internal_write_matrices

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  end subroutine couplings_calculate


  subroutine couplings_exit

    !=======================================================================!
    ! This subroutine deallocates the allocatable module variables in       !
    ! the couplings module (if required).                                   !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine on 18/12/2014.                               !
    !=======================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr

    ! ndmh: Couplings state parameters
    if (allocated(couplings_state_ids)) then
        deallocate(couplings_state_ids,stat=ierr)
        call utils_dealloc_check('couplings_exit','couplings_state_ids',ierr)
        deallocate(couplings_state_types,stat=ierr)
        call utils_dealloc_check('couplings_exit','couplings_state_types',ierr)
        deallocate(couplings_state_names,stat=ierr)
        call utils_dealloc_check('couplings_exit','couplings_state_names',ierr)
        deallocate(couplings_state_energies,stat=ierr)
        call utils_dealloc_check('couplings_exit','couplings_state_energies', &
             ierr)
        deallocate(couplings_state_exc_idxs,stat=ierr)
        call utils_dealloc_check('couplings_exit','couplings_state_exc_idxs', &
             ierr)
    end if
    ! ndmh: Couplings state parameters

  end subroutine couplings_exit

end module couplings
